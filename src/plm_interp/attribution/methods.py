"""Per-residue attributions for a single masked-marginal variant score.

Target: f(x) = log p(mt | x with position i masked) - log p(wt | same).

integrated_gradients: path integral from an "average residue" sequence (every context
    residue replaced by the mean of the 20 amino-acid embeddings) to the real sequence.
    A zero baseline (how ESM-2 embeds <mask>) is degenerate here: the pre-LayerNorm blocks
    are scale invariant, so f barely changes along the path until alpha ~ 0 and the
    integral does not converge. Completeness: sum of attributions = f(input) - f(baseline).
occlusion: change in f when one extra context residue is masked. Model-agnostic and
    needs no gradients; used as a second, independent method.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

from plm_interp.baseline.scoring import AA_INDEX
from plm_interp.models.hooks import replace_word_embeddings
from plm_interp.models.loader import ModelBundle


@dataclass
class IGResult:
    attribution: np.ndarray  # (L,) signed, per residue
    f_input: float
    f_baseline: float
    steps: int = 0

    @property
    def completeness_gap(self) -> float:
        return float(self.attribution.sum() - (self.f_input - self.f_baseline))


def _target(bundle: ModelBundle, logits: torch.Tensor, tok_pos: int, wt: str, mt: str) -> torch.Tensor:
    lp = logits[:, tok_pos].log_softmax(-1)
    ids = bundle.aa_ids
    return lp[:, ids[AA_INDEX[mt]]] - lp[:, ids[AA_INDEX[wt]]]


def _ig_pass(bundle, ids, E, E0, tok_pos, wt, mt, alphas: torch.Tensor, batch_size: int) -> torch.Tensor:
    """Sum over alphas of grad f(E0 + a (E - E0)), shape (T, D)."""
    grad_sum = torch.zeros_like(E[0])
    for b in range(0, len(alphas), batch_size):
        a = alphas[b : b + batch_size].view(-1, 1, 1)
        X = (E0 + a * (E - E0)).detach().requires_grad_(True)
        with replace_word_embeddings(bundle, X):
            logits = bundle.model(input_ids=ids.expand(len(a), -1)).logits
        f = _target(bundle, logits, tok_pos, wt, mt)
        (g,) = torch.autograd.grad(f.sum(), X)
        grad_sum += g.sum(0)
    return grad_sum


def integrated_gradients(
    bundle: ModelBundle,
    seq: str,
    pos: int,
    wt: str,
    mt: str,
    steps: int = 32,
    max_steps: int = 128,
    rel_tol: float = 0.1,
    batch_size: int = 16,
) -> IGResult:
    """Midpoint-rule IG, doubling the step count until completeness holds within rel_tol.

    f can change abruptly along the path (an attention pattern flipping), so a fixed small
    step count can miss most of the change. Worst case costs steps * (1 + 2 + 4 ...) passes.
    """
    ids = bundle.encode(seq)
    tok_pos = pos + 1
    ids[0, tok_pos] = bundle.mask_id
    with torch.no_grad():
        E = bundle.model.esm.embeddings.word_embeddings(ids)  # (1, T, D)
        mean_aa = bundle.model.esm.embeddings.word_embeddings.weight[bundle.aa_ids].mean(0)
    E0 = E.clone()
    E0[:, 1:-1] = mean_aa  # every residue -> "average residue"; <cls>/<eos> unchanged
    E0[:, tok_pos] = E[:, tok_pos]  # the masked site is zeroed by token dropout either way

    with torch.no_grad():
        f_in = _target(bundle, bundle.model(input_ids=ids).logits, tok_pos, wt, mt).item()
        with replace_word_embeddings(bundle, E0):
            f_base = _target(bundle, bundle.model(input_ids=ids).logits, tok_pos, wt, mt).item()
    delta = f_in - f_base

    # Midpoints of n uniform intervals; doubling n adds the midpoints of the new halves
    n = steps
    grad_sum = _ig_pass(bundle, ids, E, E0, tok_pos, wt, mt, (torch.arange(n) + 0.5) / n, batch_size)
    while True:
        attr = ((E - E0)[0] * grad_sum / n).sum(-1)
        gap = attr[1:-1].sum().item() - delta
        if abs(gap) <= rel_tol * max(abs(delta), 0.5) or n >= max_steps:
            break
        # Midpoints for 2n intervals = old midpoints shifted by +-1/(4n); old ones are not reused
        new_alphas = torch.cat([(torch.arange(n) + 0.25) / n, (torch.arange(n) + 0.75) / n])
        grad_sum = _ig_pass(bundle, ids, E, E0, tok_pos, wt, mt, new_alphas, batch_size)
        n *= 2

    return IGResult(attribution=attr[1:-1].detach().numpy(), f_input=f_in, f_baseline=f_base, steps=n)


@torch.no_grad()
def occlusion_logprobs(bundle: ModelBundle, seq: str, pos: int, batch_size: int = 32) -> np.ndarray:
    """Log-probs (L, 20) at `pos` when each residue j is additionally masked (row pos = only pos masked)."""
    L = len(seq)
    base = bundle.encode(seq)[0]
    base[pos + 1] = bundle.mask_id
    out = np.zeros((L, 20), dtype=np.float32)
    for b in range(0, L, batch_size):
        js = list(range(b, min(L, b + batch_size)))
        ids = base.repeat(len(js), 1)
        ids[torch.arange(len(js)), torch.tensor(js) + 1] = bundle.mask_id
        lp = bundle.model(input_ids=ids).logits[:, pos + 1].log_softmax(-1)[:, bundle.aa_ids]
        out[js] = lp.numpy()
    return out


def occlusion_attribution(occ_lp: np.ndarray, pos: int, wt: str, mt: str) -> np.ndarray:
    """f(full) - f(residue j masked): positive = residue j raised the score. Zero at pos."""
    llr = occ_lp[:, AA_INDEX[mt]] - occ_lp[:, AA_INDEX[wt]]
    attr = llr[pos] - llr
    attr[pos] = 0.0
    return attr


@torch.no_grad()
def masked_llr(bundle: ModelBundle, seq: str, pos: int, wt: str, mt: str, extra_masked: list[list[int]]) -> np.ndarray:
    """LLR at pos with extra sets of residues masked (one forward per set). Used for deletion tests."""
    base = bundle.encode(seq)[0]
    base[pos + 1] = bundle.mask_id
    ids = base.repeat(len(extra_masked), 1)
    for r, js in enumerate(extra_masked):
        if js:
            ids[r, torch.tensor(js) + 1] = bundle.mask_id
    return _target(bundle, bundle.model(input_ids=ids).logits, pos + 1, wt, mt).numpy()
