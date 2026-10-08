"""PyTorch hooks for reading and editing ESM-2 internals.

Residual stream index convention used throughout the project:
    resid[0]   = embeddings (input to layer 0)
    resid[l+1] = output of layer l   (before the final emb_layer_norm_after)

Per-head outputs are the input to each layer's attention output projection
(layer.attention.output.dense), shaped (batch, tokens, n_heads * d_head).
"""

from __future__ import annotations

from contextlib import contextmanager

import torch

from plm_interp.models.loader import ModelBundle


@contextmanager
def capture_residual(bundle: ModelBundle, detach: bool = True):
    """Collect the residual stream after embeddings and after every layer.

    Yields a list that is filled during the forward pass: n_layers + 1 tensors (B, T, D).
    """
    store: list[torch.Tensor] = []

    def save(_module, _inputs, output):
        out = output[0] if isinstance(output, tuple) else output
        store.append(out.detach() if detach else out)

    handles = [bundle.model.esm.embeddings.register_forward_hook(save)]
    handles += [layer.register_forward_hook(save) for layer in bundle.layers]
    try:
        yield store
    finally:
        for h in handles:
            h.remove()


@contextmanager
def capture_head_outputs(bundle: ModelBundle):
    """Collect per-head attention outputs for every layer: list of (B, T, H, d_head)."""
    store: list[torch.Tensor] = []
    H, dh = bundle.n_heads, bundle.d_head

    def save(_module, args):
        x = args[0]
        store.append(x.detach().view(*x.shape[:-1], H, dh))

    handles = [layer.attention.output.dense.register_forward_pre_hook(save) for layer in bundle.layers]
    try:
        yield store
    finally:
        for h in handles:
            h.remove()


@contextmanager
def ablate_heads(
    bundle: ModelBundle,
    heads: dict[int, list[int]],
    means: dict[int, torch.Tensor] | None = None,
):
    """Replace the outputs of the given heads during the forward pass.

    heads: {layer: [head, ...]}. means: optional {layer: (H, d_head)} replacement values
    (mean ablation); without it the heads are zeroed.
    """
    H, dh = bundle.n_heads, bundle.d_head
    handles = []

    for layer_idx, head_list in heads.items():
        if not head_list:
            continue
        idx = torch.tensor(sorted(head_list), dtype=torch.long)
        fill = None if means is None else means[layer_idx][idx]

        def edit(_module, args, idx=idx, fill=fill):
            x = args[0]
            x = x.view(*x.shape[:-1], H, dh).clone()
            x[..., idx, :] = 0.0 if fill is None else fill.to(x.dtype)
            return (x.view(*x.shape[:-2], H * dh),) + tuple(args[1:])

        dense = bundle.layers[layer_idx].attention.output.dense
        handles.append(dense.register_forward_pre_hook(edit))
    try:
        yield
    finally:
        for h in handles:
            h.remove()


@contextmanager
def replace_word_embeddings(bundle: ModelBundle, embeds: torch.Tensor):
    """Substitute the token-embedding lookup output (used for integrated gradients).

    The ESM token-dropout rescaling and masking that follow the lookup still apply, so
    a zero embedding behaves like a masked token.
    """

    def swap(_module, _inputs, _output):
        return embeds

    h = bundle.model.esm.embeddings.word_embeddings.register_forward_hook(swap)
    try:
        yield
    finally:
        h.remove()


@torch.no_grad()
def attention_maps(bundle: ModelBundle, seq: str) -> torch.Tensor:
    """Attention probabilities (n_layers, n_heads, L, L) over residues only (special tokens dropped).

    The bundle must be loaded with attn="eager".
    """
    ids = bundle.encode(seq)
    out = bundle.model(input_ids=ids, output_attentions=True)
    att = torch.stack(out.attentions, dim=0)[:, 0]  # (layers, heads, T, T)
    return att[:, :, 1:-1, 1:-1]
