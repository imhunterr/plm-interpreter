"""Logit lens: decode the masked position from every layer's residual stream.

Each layer's residual is passed through ESM-2's own final LayerNorm and LM head, giving a
per-layer masked-marginal score. The final entry equals the model's real output.
"""

from __future__ import annotations

import numpy as np
import torch

from plm_interp.baseline.scoring import window_for
from plm_interp.models.hooks import capture_residual
from plm_interp.models.loader import ModelBundle


@torch.no_grad()
def layerwise_logprobs(bundle: ModelBundle, seq: str, positions: list[int], batch_size: int = 16) -> np.ndarray:
    """Log-probs over the 20 amino acids at each masked position, per layer: (n_layers+1, P, 20)."""
    positions = list(positions)
    out = np.zeros((bundle.n_layers + 1, len(positions), 20), dtype=np.float32)
    aa = bundle.aa_ids
    for b in range(0, len(positions), batch_size):
        chunk = positions[b : b + batch_size]
        rows, tok = [], []
        for p in chunk:
            s, e = window_for(p, len(seq))
            ids = bundle.encode(seq[s:e])[0]
            ids[p - s + 1] = bundle.mask_id
            rows.append(ids)
            tok.append(p - s + 1)
        ids = torch.stack(rows)
        ar = torch.arange(len(chunk))
        with capture_residual(bundle) as resid:
            bundle.model(input_ids=ids)
        for layer, r in enumerate(resid):
            lp = bundle.unembed(r[ar, torch.tensor(tok)]).log_softmax(-1)[:, aa]
            out[layer, b : b + len(chunk)] = lp.numpy()
    return out
