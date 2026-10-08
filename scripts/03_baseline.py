"""Stage 1 baseline: zero-shot masked-marginal scores for every ProteinGym assay and clinical set.

Writes per-assay scores to data/processed/scores/<model>/{dms,clinical}/<id>.parquet and
summaries to outputs/runs/baseline/. Re-running skips anything already scored.

Usage:
    python scripts/03_baseline.py                 # both models, DMS + clinical
    python scripts/03_baseline.py --models 35M --only dms
    python scripts/03_baseline.py --assays RASH_HUMAN_Bandaru_2017
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from plm_interp.baseline.scoring import masked_logprobs, mutated_positions, score_mutants
from plm_interp.data.proteingym import (
    clinical_reference,
    dms_reference,
    load_clinical,
    load_dms,
    published_benchmark,
)
from plm_interp.models.loader import load_model
from plm_interp.paths import RUNS_DIR, SCORES_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("baseline", "baseline.log")


def score_dms(bundle, size: str, assays: list[str]) -> None:
    out_dir = ensure_dir(SCORES_DIR / size / "dms")
    ref = dms_reference().set_index("DMS_id")
    assays = sorted(assays, key=lambda a: ref.loc[a, "seq_len"])
    for k, dms_id in enumerate(assays):
        path = out_dir / f"{dms_id}.parquet"
        if path.exists():
            continue
        t = time.time()
        seq = ref.loc[dms_id, "target_seq"]
        df = load_dms(dms_id)
        lp = masked_logprobs(bundle, seq, mutated_positions(df.mutant))
        df["llr"] = score_mutants(df.mutant, lp, seq)
        df.drop(columns="mutated_sequence").to_parquet(path)
        log.info(f"[{size}] dms {k + 1}/{len(assays)} {dms_id} L={len(seq)} n={len(df)} {time.time() - t:.0f}s")


def score_clinical(bundle, size: str) -> None:
    out_dir = ensure_dir(SCORES_DIR / size / "clinical")
    ref = clinical_reference()
    ref = ref.assign(L=ref.target_seq.str.len()).sort_values("L")
    for k, row in enumerate(ref.itertuples()):
        path = out_dir / f"{row.DMS_id}.parquet"
        if path.exists():
            continue
        df = load_clinical(row.DMS_id)
        seq = row.target_seq
        lp = masked_logprobs(bundle, seq, mutated_positions(df.mutant))
        df["llr"] = score_mutants(df.mutant, lp, seq)
        df[["protein", "mutant", "label", "llr"]].to_parquet(path)
        if k % 100 == 0:
            log.info(f"[{size}] clinical {k + 1}/{len(ref)} {row.DMS_id} L={len(seq)}")


def summarise(sizes: list[str]) -> None:
    out_dir = ensure_dir(RUNS_DIR / "baseline")
    ref = dms_reference().set_index("DMS_id")
    pub_sp = published_benchmark("Spearman").set_index("DMS_id")
    pub_auc = published_benchmark("AUC").set_index("DMS_id")
    rows, clin_rows = [], []
    for size in sizes:
        for path in sorted((SCORES_DIR / size / "dms").glob("*.parquet")):
            df = pd.read_parquet(path)
            dms_id = path.stem
            rows.append(
                dict(
                    model=size,
                    DMS_id=dms_id,
                    seq_len=ref.loc[dms_id, "seq_len"],
                    selection_type=ref.loc[dms_id, "coarse_selection_type"],
                    n=len(df),
                    spearman=spearmanr(df.llr, df.DMS_score)[0],
                    auc=roc_auc_score(df.DMS_score_bin, df.llr) if df.DMS_score_bin.nunique() == 2 else np.nan,
                    published_spearman=pub_sp.loc[dms_id, f"ESM2 ({size})"],
                    published_auc=pub_auc.loc[dms_id, f"ESM2 ({size})"],
                )
            )
        clin_files = sorted((SCORES_DIR / size / "clinical").glob("*.parquet"))
        if clin_files:
            clin = pd.concat([pd.read_parquet(p) for p in clin_files])
            per_gene = clin.groupby("protein").filter(lambda g: g.label.nunique() == 2)
            gene_auc = per_gene.groupby("protein").apply(lambda g: roc_auc_score(g.label, -g.llr))
            clin_rows.append(
                dict(
                    model=size,
                    n_proteins=clin.protein.nunique(),
                    n_variants=len(clin),
                    n_pathogenic=int(clin.label.sum()),
                    pooled_auc=roc_auc_score(clin.label, -clin.llr),
                    n_genes_both_classes=len(gene_auc),
                    mean_per_gene_auc=gene_auc.mean(),
                    median_per_gene_auc=gene_auc.median(),
                )
            )
            gene_auc.rename("auc").to_frame().assign(model=size).to_csv(out_dir / f"clinical_per_gene_auc_{size}.csv")

    dms = pd.DataFrame(rows)
    dms["abs_diff_vs_published"] = (dms.spearman - dms.published_spearman).abs()
    dms.to_csv(out_dir / "dms_per_assay.csv", index=False)
    summary = dms.groupby("model").agg(
        n_assays=("DMS_id", "size"),
        mean_spearman=("spearman", "mean"),
        mean_published_spearman=("published_spearman", "mean"),
        max_abs_diff=("abs_diff_vs_published", "max"),
        mean_auc=("auc", "mean"),
        mean_published_auc=("published_auc", "mean"),
    )
    summary.to_csv(out_dir / "dms_summary.csv")
    log.info("\n" + summary.to_string())
    by_type = dms.groupby(["model", "selection_type"]).spearman.mean().unstack(0)
    by_type.to_csv(out_dir / "dms_by_selection_type.csv")
    if clin_rows:
        clin_df = pd.DataFrame(clin_rows)
        clin_df.to_csv(out_dir / "clinical_summary.csv", index=False)
        log.info("\n" + clin_df.to_string())


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=cfg["models"])
    ap.add_argument("--only", choices=["dms", "clinical", "summary"])
    ap.add_argument("--assays", nargs="+")
    ap.add_argument("--threads", type=int, default=cfg["threads"])
    args = ap.parse_args()
    set_threads(args.threads)

    if args.only != "summary":
        for size in args.models:
            bundle = load_model(size)
            if args.only in (None, "dms"):
                score_dms(bundle, size, args.assays or dms_reference().DMS_id.tolist())
            if args.only in (None, "clinical") and not args.assays:
                score_clinical(bundle, size)
    summarise(args.models)


if __name__ == "__main__":
    main()
