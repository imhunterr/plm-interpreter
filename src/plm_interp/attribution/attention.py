"""Attention heads as contact detectors (after Vig et al. 2021 and Rao et al. 2021)."""

from __future__ import annotations

import numpy as np
import torch


def apc(x: torch.Tensor) -> torch.Tensor:
    """Average product correction over the last two dims."""
    a1 = x.sum(-1, keepdim=True)
    a2 = x.sum(-2, keepdim=True)
    a12 = x.sum((-1, -2), keepdim=True)
    return x - a1 * a2 / a12


def head_contact_scores(att: torch.Tensor, contacts: np.ndarray, valid: np.ndarray, min_sep: int = 6) -> dict:
    """Score every head against a contact map.

    att: (layers, heads, L, L) attention probabilities over residues.
    contacts: (L, L) bool. valid: (L, L) bool pairs to evaluate (confident residues, |i-j| >= min_sep).
    Returns arrays (layers, heads):
        precision_at_L: fraction of the head's top-L long-range pairs (symmetrised + APC) that are contacts
        contact_attention: share of the head's attention mass (over valid pairs) that lands on contacts
    """
    nl, nh, L, _ = att.shape
    sym = apc(att + att.transpose(-1, -2))
    iu = np.triu_indices(L, k=min_sep)
    pair_valid = valid[iu]
    pair_contact = contacts[iu][pair_valid]
    k = min(L, int(pair_valid.sum()))

    flat_sym = sym[..., iu[0], iu[1]][..., torch.from_numpy(pair_valid)]  # (nl, nh, P)
    top = flat_sym.topk(k, dim=-1).indices
    prec = torch.from_numpy(pair_contact)[top].float().mean(-1).numpy()

    v = torch.from_numpy(valid)
    c = torch.from_numpy(contacts & valid)
    mass_valid = (att * v).sum((-1, -2))
    mass_contact = (att * c).sum((-1, -2))
    share = (mass_contact / mass_valid.clamp_min(1e-12)).numpy()
    return {"precision_at_L": prec, "contact_attention": share, "base_rate": float(pair_contact.mean())}
