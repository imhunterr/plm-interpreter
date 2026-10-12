"""Does the model score function, or evolutionary plausibility? Two tests across all assays.

Within each assay, variants are ranked by their measured effect (DMS percentile) and by the
model's score (LLR percentile).

1. Dip test. If the score tracked function, the model's percentile would keep rising into the
   best-performing tail. Measurement noise alone can flatten the top of the curve but cannot make
   it dip (with independent noise the expected true effect still rises with the measured one), so
       top_dip = mean LLR pct of the top 5% of DMS - mean LLR pct of the 80-95% band
   top_dip < 0 means the best-measured variants look less plausible than merely good ones.

2. Asymmetry test. An evolutionary prior says which changes are forbidden, not which are better.
   Then the score should rank variants well among the damaging ones and poorly among the best ones:
       asymmetry = Spearman within the bottom 25% of DMS - Spearman within the top 25%
   asymmetry > 0 means the model separates damaging variants far better than beneficial ones.

Outputs (outputs/runs/evolution_vs_function/): curve.csv, per_assay.csv, summary.csv

Usage:
    python scripts/14_evolution_vs_function.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from scipy.stats import spearmanr

from plm_interp.data.proteingym import dms_reference
from plm_interp.paths import RUNS_DIR, SCORES_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config
from plm_interp.verification.stats import bootstrap_ci, one_sample_permutation_p

log = get_logger("evo_vs_function", "evolution_vs_function.log")

N_BINS = 20  # 5% bins of DMS percentile
MIN_SINGLES = 200


def main() -> None:
    cfg = load_config()
    out = ensure_dir(RUNS_DIR / "evolution_vs_function")
    ref = dms_reference().set_index("DMS_id")
    curve_rows, assay_rows = [], []

    for size in cfg["models"]:
        for path in sorted((SCORES_DIR / size / "dms").glob("*.parquet")):
            dms_id = path.stem
            df = pd.read_parquet(path)
            df = df[~df.mutant.str.contains(":")]
            if len(df) < MIN_SINGLES:
                continue
            dms_pct = df.DMS_score.rank(pct=True).to_numpy()
            llr_pct = df.llr.rank(pct=True).to_numpy() * 100
            bins = np.minimum((dms_pct * N_BINS).astype(int), N_BINS - 1)
            means = pd.Series(llr_pct).groupby(bins).mean()
            for b, v in means.items():
                curve_rows.append(dict(model=size, dms_id=dms_id, bin=int(b), dms_pct_mid=(b + 0.5) * 100 / N_BINS, llr_pct=v))
            top = llr_pct[dms_pct > 0.95]
            band = llr_pct[(dms_pct > 0.80) & (dms_pct <= 0.95)]
            bottom = llr_pct[dms_pct <= 0.05]
            lo_q, hi_q = dms_pct <= 0.25, dms_pct > 0.75
            rho_bottom = spearmanr(df.llr.to_numpy()[lo_q], df.DMS_score.to_numpy()[lo_q])[0]
            rho_top = spearmanr(df.llr.to_numpy()[hi_q], df.DMS_score.to_numpy()[hi_q])[0]
            assay_rows.append(
                dict(
                    model=size,
                    dms_id=dms_id,
                    selection_type=ref.loc[dms_id, "coarse_selection_type"],
                    n=len(df),
                    llr_pct_top5=top.mean(),
                    llr_pct_80_95=band.mean(),
                    llr_pct_bottom5=bottom.mean(),
                    top_dip=top.mean() - band.mean(),
                    rho_bottom_quarter=rho_bottom,
                    rho_top_quarter=rho_top,
                    asymmetry=rho_bottom - rho_top,
                    # share of the best 5% that the model puts in its own bottom quarter
                    top5_flagged=float((top < 25).mean()),
                    band_flagged=float((band < 25).mean()),
                )
            )

    curve = pd.DataFrame(curve_rows)
    per_assay = pd.DataFrame(assay_rows)
    curve.to_csv(out / "curve.csv", index=False)
    per_assay.to_csv(out / "per_assay.csv", index=False)

    rows = []
    for (size, group), g in [((m, "all"), per_assay[per_assay.model == m]) for m in cfg["models"]] + [
        ((m, t), g) for (m, t), g in per_assay.groupby(["model", "selection_type"])
    ]:
        est, lo, hi = bootstrap_ci(g.top_dip, cfg["stats"]["n_bootstrap"], seed=cfg["seed"])
        a_est, a_lo, a_hi = bootstrap_ci(g.asymmetry, cfg["stats"]["n_bootstrap"], seed=cfg["seed"])
        rows.append(
            dict(
                model=size,
                group=group,
                n_assays=len(g),
                top_dip=est,
                ci_low=lo,
                ci_high=hi,
                p_value=one_sample_permutation_p(g.top_dip, 0.0, n_perm=10000, seed=cfg["seed"]),
                rho_bottom_quarter=g.rho_bottom_quarter.mean(),
                rho_top_quarter=g.rho_top_quarter.mean(),
                asymmetry=a_est,
                asymmetry_ci_low=a_lo,
                asymmetry_ci_high=a_hi,
                asymmetry_p=one_sample_permutation_p(g.asymmetry, 0.0, n_perm=10000, seed=cfg["seed"]),
                share_assays_asymmetric=float((g.asymmetry.dropna() > 0).mean()),
                share_assays_with_dip=float((g.top_dip < 0).mean()),
                top5_flagged=g.top5_flagged.mean(),
                band_flagged=g.band_flagged.mean(),
            )
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "summary.csv", index=False)
    log.info("\n" + summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
