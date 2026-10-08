"""Masked-marginal variant scoring (Meier et al. 2021).

For a substitution wt -> mt at position i, mask position i and score
    LLR = log p(mt | x_{-i}) - log p(wt | x_{-i}).
Multi-mutants are scored additively, as in ProteinGym.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch

from plm_interp.data.proteingym import parse_mutant
from plm_interp.models.loader import AMINO_ACIDS, MAX_RESIDUES, ModelBundle

AA_INDEX = {a: i for i, a in enumerate(AMINO_ACIDS)}


def window_for(pos: int, seq_len: int, max_len: int = MAX_RESIDUES) -> tuple[int, int]:
    """[start, end) window of at most max_len residues centred on pos (ProteinGym convention)."""
    if seq_len <= max_len:
        return 0, seq_len
    half = max_len // 2
    start = max(0, min(pos - half, seq_len - max_len))
    return start, start + max_len


@torch.no_grad()
def masked_logprobs(
    bundle: ModelBundle,
    seq: str,
    positions: list[int] | None = None,
    batch_size: int = 16,
) -> pd.DataFrame:
    """Log-probabilities of the 20 amino acids at each position with that position masked.

    Returns a DataFrame indexed by 0-based position with one column per amino acid.
    Sequences longer than the model limit are scored in a window centred on each position.
    """
    positions = list(range(len(seq))) if positions is None else sorted(set(positions))
    aa_ids = bundle.aa_ids
    rows: dict[int, np.ndarray] = {}

    # Every window has the same length (the whole sequence, or max_len for long ones),
    # so positions batch together even when their windows differ.
    windows = {p: window_for(p, len(seq)) for p in positions}
    encoded: dict[tuple[int, int], torch.Tensor] = {}
    for b in range(0, len(positions), batch_size):
        chunk = positions[b : b + batch_size]
        rows_ids = []
        for p in chunk:
            w = windows[p]
            if w not in encoded:
                encoded[w] = bundle.encode(seq[w[0] : w[1]])[0]
            rows_ids.append(encoded[w])
        ids = torch.stack(rows_ids).clone()
        tok_pos = torch.tensor([p - windows[p][0] + 1 for p in chunk])
        ar = torch.arange(len(chunk))
        ids[ar, tok_pos] = bundle.mask_id
        logits = bundle.model(input_ids=ids).logits
        lp = logits[ar, tok_pos].log_softmax(-1)[:, aa_ids]
        for p, v in zip(chunk, lp.numpy()):
            rows[p] = v
        if len(encoded) > 4 * batch_size:
            encoded.clear()

    return pd.DataFrame.from_dict(rows, orient="index", columns=list(AMINO_ACIDS)).sort_index()


def score_mutants(mutants: pd.Series | list[str], logprobs: pd.DataFrame, seq: str) -> np.ndarray:
    """Additive masked-marginal LLR for each mutant string, using a precomputed logprob table."""
    lp = logprobs.to_numpy()
    row_of = {p: r for r, p in enumerate(logprobs.index)}
    scores = np.empty(len(mutants))
    for k, mutant in enumerate(mutants):
        total = 0.0
        for wt, pos, mt in parse_mutant(mutant):
            if seq[pos] != wt:
                raise ValueError(f"{mutant}: sequence has {seq[pos]} at {pos + 1}, not {wt}")
            r = row_of[pos]
            total += lp[r, AA_INDEX[mt]] - lp[r, AA_INDEX[wt]]
        scores[k] = total
    return scores


def mutated_positions(mutants) -> list[int]:
    return sorted({pos for m in mutants for _, pos, _ in parse_mutant(m)})


def llr_matrix(logprobs: pd.DataFrame, seq: str) -> pd.DataFrame:
    """Position x amino-acid table of single-substitution LLRs (wild type column = 0)."""
    wt_lp = np.array([logprobs.loc[p, seq[p]] for p in logprobs.index])
    return logprobs.sub(wt_lp, axis=0)
