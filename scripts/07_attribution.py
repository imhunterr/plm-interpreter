"""Stage 5: which residues drove each variant score, and are they the variant's 3D neighbours?

For each showcase protein, a fixed sample of single mutants (one per position, stratified over
the DMS score range, confident structure only) is explained with integrated gradients and,
for a subset, occlusion. Each attribution is scored by how well it picks out 3D contacts of
the mutated site, and by a deletion test (masking top-attributed residues vs random ones).
A random-init model gives the chance level for the method.

Outputs (outputs/runs/attribution/):
    variants.csv        one row per (model, variant): scores, completeness, contact AUROCs, deletion results
    attributions.npz    raw per-residue attribution vectors, keyed "<model>|<dms_id>|<mutant>|<method>"

Usage:
    python scripts/07_attribution.py [--models 35M] [--proteins RASH_HUMAN_Bandaru_2017]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from plm_interp.attribution.methods import (
    integrated_gradients,
    masked_llr,
    occlusion_attribution,
    occlusion_logprobs,
)
from plm_interp.data.proteingym import dms_info, load_dms, parse_mutant
from plm_interp.models.loader import load_model
from plm_interp.paths import ANNOTATIONS_DIR, RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("attribution", "attribution.log")


def select_variants(dms: pd.DataFrame, confident: np.ndarray, n: int, seed: int) -> pd.DataFrame:
    """One single mutant per position, spread evenly over the DMS score quantiles."""
    df = dms[~dms.mutant.str.contains(":")].copy()
    df["pos"] = [parse_mutant(m)[0][1] for m in df.mutant]
    df = df[confident[df.pos.to_numpy()]]
    df["bin"] = pd.qcut(df.DMS_score.rank(method="first"), q=min(n, len(df)), labels=False)
    rng = np.random.default_rng(seed)
    picked, used = [], set()
    for _, g in df.sample(frac=1.0, random_state=seed).groupby("bin"):
        g = g[~g.pos.isin(used)]
        if len(g):
            row = g.iloc[rng.integers(len(g))]
            picked.append(row)
            used.add(row.pos)
    return pd.DataFrame(picked).reset_index(drop=True)


def contact_aurocs(attr: np.ndarray, dist: np.ndarray, pos: int, confident: np.ndarray, cutoff: float, min_sep: int):
    L = len(attr)
    j = np.arange(L)
    ok = confident & (j != pos)
    contact = dist[pos] < cutoff
    imp = np.abs(attr)
    out = {}
    for name, mask in (("all", ok), ("long_range", ok & (np.abs(j - pos) >= min_sep))):
        y = contact[mask]
        out[f"auroc_{name}"] = roc_auc_score(y, imp[mask]) if 0 < y.sum() < len(y) else np.nan
    # Trivial baseline: closeness in sequence
    y = contact[ok]
    out["auroc_seqdist_baseline"] = roc_auc_score(y, -np.abs(j - pos)[ok]) if 0 < y.sum() < len(y) else np.nan
    # Share of the top-10 attributed residues that are 3D contacts
    top = j[ok][np.argsort(-imp[ok])[:10]]
    out["top10_contact_frac"] = contact[top].mean()
    return out


def deletion_test(bundle, seq, pos, wt, mt, attr, ks, n_random, rng) -> dict:
    """|change in LLR| after masking the top-k attributed residues vs k random / k nearest-in-sequence."""
    L = len(seq)
    others = np.array([j for j in range(L) if j != pos])
    ranked = others[np.argsort(-np.abs(attr[others]))]
    nearest = others[np.argsort(np.abs(others - pos))]
    sets, labels = [[]], [("none", 0, 0)]
    for k in ks:
        sets.append(ranked[:k].tolist())
        labels.append(("top", k, 0))
        sets.append(nearest[:k].tolist())
        labels.append(("seq_nearest", k, 0))
        for r in range(n_random):
            sets.append(rng.choice(others, size=k, replace=False).tolist())
            labels.append(("random", k, r))
    llr = masked_llr(bundle, seq, pos, wt, mt, sets)
    base = llr[0]
    out = {}
    for k in ks:
        top = [abs(llr[i] - base) for i, lab in enumerate(labels) if lab[0] == "top" and lab[1] == k][0]
        near = [abs(llr[i] - base) for i, lab in enumerate(labels) if lab[0] == "seq_nearest" and lab[1] == k][0]
        rand = np.mean([abs(llr[i] - base) for i, lab in enumerate(labels) if lab[0] == "random" and lab[1] == k])
        out[f"del_top_{k}"] = top
        out[f"del_near_{k}"] = near
        out[f"del_rand_{k}"] = rand
    return out


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=cfg["models"])
    ap.add_argument("--proteins", nargs="+", default=cfg["showcase"])
    ap.add_argument("--threads", type=int, default=cfg["threads"])
    args = ap.parse_args()
    set_threads(args.threads)
    a_cfg, s_cfg = cfg["attribution"], cfg["structure"]
    out = ensure_dir(RUNS_DIR / "attribution")
    residues = pd.read_parquet(ANNOTATIONS_DIR / "residues.parquet")

    rows, vectors = [], {}
    runs = [(m, False) for m in args.models] + [(m, True) for m in args.models]
    for size, random_init in runs:
        bundle = load_model(size, random_init=random_init, seed=cfg["seed"])
        for dms_id in args.proteins:
            info = dms_info(dms_id)
            seq = info.target_seq
            res = residues[residues.uniprot == info.UniProt_ID]
            confident = res.confident.to_numpy()
            dist = np.load(ANNOTATIONS_DIR / "dist" / f"{info.UniProt_ID}.npy").astype(np.float32)
            variants = select_variants(load_dms(dms_id), confident, a_cfg["n_variants_per_protein"], cfg["seed"])
            n_ig = a_cfg["n_random_init_per_protein"] if random_init else len(variants)
            rng = np.random.default_rng(cfg["seed"])
            t0 = time.time()
            for k, v in enumerate(variants.head(n_ig).itertuples()):
                (wt, pos, mt) = parse_mutant(v.mutant)[0]
                ig = integrated_gradients(
                    bundle, seq, pos, wt, mt, a_cfg["ig_steps"], a_cfg["ig_max_steps"], a_cfg["ig_rel_tol"]
                )
                delta = ig.f_input - ig.f_baseline
                row = dict(
                    model=bundle.name,
                    random_init=random_init,
                    dms_id=dms_id,
                    mutant=v.mutant,
                    pos=pos,
                    DMS_score=v.DMS_score,
                    DMS_score_bin=v.DMS_score_bin,
                    llr=ig.f_input,
                    ig_baseline_llr=ig.f_baseline,
                    ig_steps=ig.steps,
                    ig_gap=ig.completeness_gap,
                    ig_converged=abs(ig.completeness_gap) <= a_cfg["ig_rel_tol"] * max(abs(delta), 0.5),
                    rsa=res.rsa.iloc[pos],
                    ss=res.ss.iloc[pos],
                    conservation=res.conservation.iloc[pos],
                )
                row.update({f"ig_{key}": val for key, val in contact_aurocs(
                    ig.attribution, dist, pos, confident, s_cfg["contact_cutoff"], s_cfg["min_seq_sep"]).items()})
                row.update(deletion_test(bundle, seq, pos, wt, mt, ig.attribution,
                                         a_cfg["deletion_k"], a_cfg["n_random_deletions"], rng))
                vectors[f"{bundle.name}|{dms_id}|{v.mutant}|ig"] = ig.attribution.astype(np.float32)

                if not random_init and k < a_cfg["n_occlusion_per_protein"]:
                    occ = occlusion_attribution(occlusion_logprobs(bundle, seq, pos), pos, wt, mt)
                    row.update({f"occ_{key}": val for key, val in contact_aurocs(
                        occ, dist, pos, confident, s_cfg["contact_cutoff"], s_cfg["min_seq_sep"]).items()})
                    row["ig_occ_spearman"] = spearmanr(np.abs(ig.attribution), np.abs(occ))[0]
                    vectors[f"{bundle.name}|{dms_id}|{v.mutant}|occlusion"] = occ.astype(np.float32)
                rows.append(row)
            log.info(f"{bundle.name} {dms_id}: {n_ig} variants in {time.time() - t0:.0f}s")
            pd.DataFrame(rows).to_csv(out / "variants.csv", index=False)
            np.savez_compressed(out / "attributions.npz", **vectors)

    df = pd.DataFrame(rows)
    summary = df.groupby(["model"]).agg(
        n=("mutant", "size"),
        converged=("ig_converged", "mean"),
        ig_auroc_all=("ig_auroc_all", "mean"),
        ig_auroc_long_range=("ig_auroc_long_range", "mean"),
        seqdist_baseline=("ig_auroc_seqdist_baseline", "mean"),
        ig_top10_contact_frac=("ig_top10_contact_frac", "mean"),
    )
    log.info("\n" + summary.to_string())


if __name__ == "__main__":
    main()
