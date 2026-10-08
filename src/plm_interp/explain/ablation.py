"""Component ablation: which heads / MLPs does variant-effect prediction depend on?"""

from __future__ import annotations

from contextlib import contextmanager

import torch

from plm_interp.models.hooks import capture_head_outputs
from plm_interp.models.loader import ModelBundle


@torch.no_grad()
def head_means(bundle: ModelBundle, seq: str) -> dict[int, torch.Tensor]:
    """Mean per-head output over residue tokens of the unmasked sequence: {layer: (H, d_head)}."""
    with capture_head_outputs(bundle) as store:
        bundle.model(input_ids=bundle.encode(seq))
    return {layer: h[0, 1:-1].mean(0) for layer, h in enumerate(store)}


@torch.no_grad()
def mlp_means(bundle: ModelBundle, seq: str) -> dict[int, torch.Tensor]:
    """Mean MLP output (before the residual add) over residue tokens: {layer: (D,)}."""
    store = {}
    handles = []
    for i, layer in enumerate(bundle.layers):
        def save(_m, _inp, out, i=i):
            store[i] = out[0, 1:-1].mean(0)
        handles.append(layer.output.dense.register_forward_hook(save))
    try:
        bundle.model(input_ids=bundle.encode(seq))
    finally:
        for h in handles:
            h.remove()
    return store


@contextmanager
def ablate_mlps(bundle: ModelBundle, layers: list[int], means: dict[int, torch.Tensor] | None = None):
    """Replace the MLP output of the given layers by its mean (or zero)."""
    handles = []
    for i in layers:
        fill = None if means is None else means[i]

        def edit(_m, _inp, out, fill=fill):
            return torch.zeros_like(out) if fill is None else fill.expand_as(out).clone()

        handles.append(bundle.layers[i].output.dense.register_forward_hook(edit))
    try:
        yield
    finally:
        for h in handles:
            h.remove()
