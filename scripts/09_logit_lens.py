"""Stage 7: at which layer does the variant-effect signal appear?

Decodes the masked position from every layer (logit lens) and correlates the per-layer score
with the DMS measurements. Run on all panel assays (<= 400 residues, one assay per protein).

Output: outputs/runs/logit_lens/layers.csv (model x assay x layer: Spearman, top-1 accuracy)

Usage:
    python scripts/09_logit_lens.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from plm_interp.baseline.scoring import AA_INDEX, mutated_positions
from plm_interp.data.proteingym import dms_info, load_dms, parse_mutant
from plm_interp.explain.logit_lens import layerwise_logprobs
from plm_interp.models.loader import AMINO_ACIDS, load_model
from plm_interp.paths import ANNOTATIONS_DIR, RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("logit_lens", "logit_lens.log")


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=cfg["models"])
    ap.add_argument("--threads", type=int, default=cfg["threads"])
    args = ap.parse_args()
    set_threads(args.threads)
    out = ensure_dir(RUNS_DIR / "logit_lens")
    proteins = pd.read_csv(ANNOTATIONS_DIR / "proteins.csv")

    rows = []
    for size in args.models:
        bundle = load_model(size)
        for dms_id in proteins.dms_id:
            seq = dms_info(dms_id).target_seq
            dms = load_dms(dms_id)
            dms = dms[~dms.mutant.str.contains(":")]
            positions = mutated_positions(dms.mutant)
            lp = layerwise_logprobs(bundle, seq, positions)  # (layers+1, P, 20)
            row_of = {p: i for i, p in enumerate(positions)}
            muts = [parse_mutant(m)[0] for m in dms.mutant]
            r_idx = np.array([row_of[p] for _, p, _ in muts])
            wt_idx = np.array([AA_INDEX[w] for w, _, _ in muts])
            mt_idx = np.array([AA_INDEX[m] for _, _, m in muts])
            true_aa = np.array([AMINO_ACIDS.index(seq[p]) if seq[p] in AMINO_ACIDS else -1 for p in positions])
            for layer in range(lp.shape[0]):
                llr = lp[layer, r_idx, mt_idx] - lp[layer, r_idx, wt_idx]
                rows.append(
                    dict(
                        model=size,
                        dms_id=dms_id,
                        layer=layer,
                        rel_depth=layer / (lp.shape[0] - 1),
                        spearman=spearmanr(llr, dms.DMS_score)[0],
                        wt_top1=float((lp[layer].argmax(-1) == true_aa).mean()),
                        wt_logprob=float(lp[layer, np.arange(len(positions)), true_aa].mean()),
                        n=len(dms),
                    )
                )
            log.info(f"{size} {dms_id}: final rho {rows[-1]['spearman']:.3f}")
            pd.DataFrame(rows).to_csv(out / "layers.csv", index=False)

    df = pd.DataFrame(rows)
    log.info("\n" + df.groupby(["model", "layer"]).spearman.mean().unstack(0).round(3).to_string())


if __name__ == "__main__":
    main()
