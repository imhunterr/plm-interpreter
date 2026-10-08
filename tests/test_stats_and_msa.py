import numpy as np

from plm_interp.data.msa import conservation
from plm_interp.verification.stats import benjamini_hochberg, bootstrap_ci, paired_permutation_test


def test_conservation_extremes():
    # column 0 fully conserved, column 1 maximally variable over 20 residues
    aas = "ACDEFGHIKLMNPQRSTVWY"
    seqs = [f"M{a}" for a in aas]
    cons, neff = conservation(seqs)
    assert cons[0] > 0.99
    assert cons[1] < cons[0]
    assert neff == len(seqs)  # all sequences are < 80% identical (50%), so weight 1 each


def test_identical_sequences_share_weight():
    cons, neff = conservation(["MKV"] * 10 + ["WWW"])
    assert neff == 2.0


def test_bh_monotone_and_bounded():
    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.5])
    assert np.all(q <= 1) and q[0] <= q[1]
    assert np.isclose(q[0], 0.04)


def test_bootstrap_contains_mean():
    est, lo, hi = bootstrap_ci(np.arange(100.0))
    assert lo < est < hi


def test_permutation_detects_shift():
    rng = np.random.default_rng(0)
    a = rng.normal(1, 1, 50)
    assert paired_permutation_test(a, np.zeros(50), n_perm=2000) < 0.01
    assert paired_permutation_test(a, a, n_perm=2000) == 1.0
