import numpy as np
import pytest
import torch

from plm_interp.attribution.attention import apc, head_contact_scores
from plm_interp.attribution.methods import integrated_gradients, occlusion_attribution, occlusion_logprobs
from plm_interp.baseline.scoring import masked_logprobs
from plm_interp.explain.ablation import ablate_mlps, head_means, mlp_means
from plm_interp.explain.logit_lens import layerwise_logprobs
from plm_interp.models.hooks import ablate_heads, capture_head_outputs, capture_residual


def test_capture_residual_shapes(bundle, seq):
    with capture_residual(bundle) as resid:
        bundle.model(input_ids=bundle.encode(seq))
    assert len(resid) == bundle.n_layers + 1
    assert resid[0].shape == (1, len(seq) + 2, bundle.d_model)


def test_head_outputs_shape(bundle, seq):
    with capture_head_outputs(bundle) as heads:
        bundle.model(input_ids=bundle.encode(seq))
    assert len(heads) == bundle.n_layers
    assert heads[0].shape == (1, len(seq) + 2, bundle.n_heads, bundle.d_head)


def test_logit_lens_last_layer_matches_model(bundle, seq):
    positions = [3, 20, 50]
    lens = layerwise_logprobs(bundle, seq, positions)
    ref = masked_logprobs(bundle, seq, positions).to_numpy()
    np.testing.assert_allclose(lens[-1], ref, atol=1e-4)


def test_empty_ablation_is_noop_and_real_ablation_changes_output(bundle, seq):
    ref = masked_logprobs(bundle, seq, [10]).to_numpy()
    means = head_means(bundle, seq)
    with ablate_heads(bundle, {1: []}, means):
        same = masked_logprobs(bundle, seq, [10]).to_numpy()
    with ablate_heads(bundle, {1: list(range(bundle.n_heads))}, means):
        diff = masked_logprobs(bundle, seq, [10]).to_numpy()
    with ablate_mlps(bundle, [0], mlp_means(bundle, seq)):
        diff_mlp = masked_logprobs(bundle, seq, [10]).to_numpy()
    np.testing.assert_array_equal(ref, same)
    assert np.abs(ref - diff).max() > 1e-3
    assert np.abs(ref - diff_mlp).max() > 1e-3
    # hooks are removed after the context exits
    np.testing.assert_array_equal(ref, masked_logprobs(bundle, seq, [10]).to_numpy())


def test_integrated_gradients_completeness(bundle, seq):
    pos = 10
    r = integrated_gradients(bundle, seq, pos, seq[pos], "W", steps=32, max_steps=256)
    assert r.attribution.shape == (len(seq),)
    assert abs(r.completeness_gap) <= 0.1 * max(abs(r.f_input - r.f_baseline), 0.5)


def test_occlusion_is_zero_at_site(bundle, seq):
    pos = 5
    occ = occlusion_attribution(occlusion_logprobs(bundle, seq, pos), pos, seq[pos], "A")
    assert occ[pos] == 0.0 and np.abs(occ).sum() > 0


def test_apc_removes_row_column_means():
    x = torch.rand(6, 6) + 1
    rank1 = torch.outer(x.sum(1), x.sum(0)) / x.sum()
    assert torch.allclose(apc(rank1), torch.zeros(6, 6), atol=1e-6)


def test_head_contact_scores_perfect_head():
    L = 30
    contacts = np.zeros((L, L), dtype=bool)
    for i in range(0, L - 10, 3):
        contacts[i, i + 10] = contacts[i + 10, i] = True
    valid = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :]) >= 6
    att = torch.full((1, 2, L, L), 1e-4)
    att[0, 0][torch.from_numpy(contacts)] = 1.0  # head 0 attends exactly to contacts
    sc = head_contact_scores(att, contacts, valid)
    assert sc["contact_attention"][0, 0] > 0.9
    assert sc["contact_attention"][0, 0] > sc["contact_attention"][0, 1]
