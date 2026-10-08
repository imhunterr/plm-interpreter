"""Statistics used across the trust checks: bootstrap CIs, paired permutation tests, BH-FDR."""

from __future__ import annotations

import numpy as np


def bootstrap_ci(values, n_boot: int = 1000, alpha: float = 0.05, seed: int = 0, stat=np.nanmean) -> tuple[float, float, float]:
    """(estimate, low, high) percentile bootstrap CI of a statistic over independent units."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if len(v) == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.array([stat(v[rng.integers(0, len(v), len(v))]) for _ in range(n_boot)])
    return float(stat(v)), float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2))


def paired_permutation_test(a, b, n_perm: int = 10000, seed: int = 0) -> float:
    """Two-sided sign-flip test of mean(a - b) = 0 for paired samples."""
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    d = d[~np.isnan(d)]
    if len(d) == 0:
        return np.nan
    rng = np.random.default_rng(seed)
    obs = abs(d.mean())
    signs = rng.choice([-1.0, 1.0], size=(n_perm, len(d)))
    null = np.abs((signs * d).mean(1))
    return float((1 + (null >= obs).sum()) / (n_perm + 1))


def one_sample_permutation_p(values, null_mean: float, n_perm: int = 10000, seed: int = 0) -> float:
    """Sign-flip test of mean(values) = null_mean."""
    v = np.asarray(values, dtype=float)
    return paired_permutation_test(v, np.full_like(v, null_mean), n_perm, seed)


def benjamini_hochberg(p) -> np.ndarray:
    """BH-adjusted q-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out
