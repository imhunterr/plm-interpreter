"""Stage 8: trust checks. Turns the raw outputs of Stages 1-7 into tested, paper-ready claims.

Every claim gets an effect size, a 95% bootstrap CI over independent units (variants or
proteins), and a permutation-test p-value against its control. Claims:

  C1 baseline   our ESM-2 scores reproduce ProteinGym's published numbers
  C2 probing    trained layers encode structure/conservation beyond residue identity and a random-init model
  C3 attention  some heads attend to 3D contacts far above the base rate, unlike random-init heads
  C4 attribution IG importance picks out 3D contacts of the mutated site (incl. long-range), above chance
                 and above the random-init model; IG agrees with occlusion
  C5 deletion   masking the top-attributed residues changes the score more than random / sequence-nearest ones
  C6 ablation   the top-k heads (ranked on held-out positions) matter more than random k heads
  C7 lens       the variant signal emerges in late layers

Output: outputs/runs/verification/claims.csv and claims.md

Usage:
    python scripts/10_verification.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from plm_interp.paths import RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config
from plm_interp.verification.stats import (
    benjamini_hochberg,
    bootstrap_ci,
    one_sample_permutation_p,
    paired_permutation_test,
)

log = get_logger("verification", "verification.log")


def claim(rows, cid, model, text, values, null=None, paired=None, unit="", cfg=None):
    est, lo, hi = bootstrap_ci(values, cfg["stats"]["n_bootstrap"], seed=cfg["seed"])
    if paired is not None:
        p = paired_permutation_test(values, paired, n_perm=cfg["stats"]["n_permutations"] * 10, seed=cfg["seed"])
        ctrl = float(np.nanmean(paired))
    elif null is not None:
        p = one_sample_permutation_p(values, null, n_perm=cfg["stats"]["n_permutations"] * 10, seed=cfg["seed"])
        ctrl = null
    else:
        p, ctrl = np.nan, np.nan
    rows.append(
        dict(claim=cid, model=model, statement=text, estimate=est, ci_low=lo, ci_high=hi, control=ctrl,
             p_value=p, n=int(np.sum(~np.isnan(np.asarray(values, dtype=float)))), unit=unit)
    )


def main() -> None:
    cfg = load_config()
    out = ensure_dir(RUNS_DIR / "verification")
    rows: list[dict] = []

    # C1 baseline reproduction
    f = RUNS_DIR / "baseline" / "dms_per_assay.csv"
    if f.exists():
        dms = pd.read_csv(f)
        for m, g in dms.groupby("model"):
            claim(rows, "C1", m, "mean |our Spearman - published Spearman| per assay",
                  (g.spearman - g.published_spearman).abs(), unit="assays", cfg=cfg)
            claim(rows, "C1", m, "mean Spearman with DMS (zero-shot)", g.spearman, null=0.0, unit="assays", cfg=cfg)

    # C2 probing: best trained layer vs aa-identity and vs best random-init layer (point estimates per task)
    f = RUNS_DIR / "probing" / "probe_results.csv"
    if f.exists():
        pr = pd.read_csv(f)
        aa = pr[pr.model == "aa_identity"].set_index("task").value
        for m in cfg["models"]:
            trained = pr[(pr.model == m) & ~pr.random_init]
            rnd = pr[(pr.model == m) & pr.random_init]
            for task, g in trained.groupby("task"):
                best = g.loc[g.value.idxmax()]
                rows.append(dict(
                    claim="C2", model=m,
                    statement=f"probe {task}: best layer {int(best.layer)} ({best.metric}) vs aa identity / random-init best",
                    estimate=best.value, ci_low=np.nan, ci_high=np.nan,
                    control=max(aa.get(task, np.nan), rnd[rnd.task == task].value.max()),
                    p_value=np.nan, n=int(best.layer), unit="held-out residues",
                ))

    # C3 attention heads vs contacts (per-protein precision of the best head vs random-init best head)
    f = RUNS_DIR / "attention" / "head_contacts.csv"
    if f.exists():
        hc = pd.read_csv(f)
        summ = pd.read_csv(RUNS_DIR / "attention" / "head_summary.csv")
        for m in cfg["models"]:
            s = summ[(summ.model == m) & ~summ.random_init]
            best = s.loc[s.precision_at_L.idxmax()]
            sel = hc[(hc.model == m) & ~hc.random_init & (hc.layer == best.layer) & (hc["head"] == best["head"])]
            s_r = summ[(summ.model == m) & summ.random_init]
            best_r = s_r.loc[s_r.precision_at_L.idxmax()]
            sel_r = hc[(hc.model == m) & hc.random_init & (hc.layer == best_r.layer) & (hc["head"] == best_r["head"])]
            merged = sel.merge(sel_r, on="uniprot", suffixes=("", "_r"))
            claim(rows, "C3", m, f"precision@L of best contact head L{int(best.layer)}H{int(best['head'])} vs best random-init head",
                  merged.precision_at_L, paired=merged.precision_at_L_r, unit="proteins", cfg=cfg)
            claim(rows, "C3", m, "same head: precision@L vs contact base rate",
                  sel.precision_at_L, paired=sel.base_rate, unit="proteins", cfg=cfg)

    # C4 / C5 attribution and deletion
    f = RUNS_DIR / "attribution" / "variants.csv"
    if f.exists():
        va = pd.read_csv(f)
        for m in cfg["models"]:
            t = va[va.model == m]
            r = va[va.model == f"{m}-random"]
            claim(rows, "C4", m, "IG |attribution| AUROC for 3D contacts of the mutated site", t.ig_auroc_all, null=0.5, unit="variants", cfg=cfg)
            claim(rows, "C4", m, "IG AUROC for long-range contacts (|i-j|>=6)", t.ig_auroc_long_range, null=0.5, unit="variants", cfg=cfg)
            if len(r):
                pair = t.merge(r, on=["dms_id", "mutant"], suffixes=("", "_rand"))
                claim(rows, "C4", m, "IG contact AUROC, trained vs random-init model (same variants)",
                      pair.ig_auroc_all, paired=pair.ig_auroc_all_rand, unit="variants", cfg=cfg)
            claim(rows, "C4", m, "share of top-10 IG residues that are 3D contacts", t.ig_top10_contact_frac, unit="variants", cfg=cfg)
            claim(rows, "C4", m, "fraction of IG runs meeting completeness (<=10% gap)", t.ig_converged.astype(float), unit="variants", cfg=cfg)
            if "ig_occ_spearman" in t:
                claim(rows, "C4", m, "Spearman(|IG|, |occlusion|) per variant", t.ig_occ_spearman, null=0.0, unit="variants", cfg=cfg)
                claim(rows, "C4", m, "occlusion AUROC for 3D contacts", t.occ_auroc_all, null=0.5, unit="variants", cfg=cfg)
            for k in cfg["attribution"]["deletion_k"]:
                claim(rows, "C5", m, f"|dLLR| masking top-{k} IG residues vs {k} random", t[f"del_top_{k}"], paired=t[f"del_rand_{k}"], unit="variants", cfg=cfg)
                claim(rows, "C5", m, f"|dLLR| masking top-{k} IG residues vs {k} sequence-nearest", t[f"del_top_{k}"], paired=t[f"del_near_{k}"], unit="variants", cfg=cfg)

    # C6 ablation
    f = RUNS_DIR / "ablation" / "topk.csv"
    if f.exists():
        tk = pd.read_csv(f)
        for m, g in tk.groupby("model"):
            for k, gk in g.groupby("k"):
                drop_top = gk.clean_rho_eval - gk.top_rho_eval
                drop_rand = gk.clean_rho_eval - gk.random_rho_eval_mean
                claim(rows, "C6", m, f"held-out Spearman drop: ablate top-{k} heads vs {k} random heads", drop_top, paired=drop_rand, unit="proteins", cfg=cfg)
                if "contact_rho_eval" in gk:
                    claim(rows, "C6", m, f"held-out Spearman drop: ablate top-{k} contact heads", gk.clean_rho_eval - gk.contact_rho_eval, paired=drop_rand, unit="proteins", cfg=cfg)

    # C7 logit lens: fraction of final Spearman reached at each relative depth
    f = RUNS_DIR / "logit_lens" / "layers.csv"
    if f.exists():
        ll = pd.read_csv(f)
        for m, g in ll.groupby("model"):
            final = g[g.layer == g.layer.max()].set_index("dms_id").spearman
            for layer, gl in g.groupby("layer"):
                frac = gl.set_index("dms_id").spearman / final
                frac = frac[final > 0.1]
                claim(rows, "C7", m, f"layer {layer}: Spearman as a fraction of the final layer's", frac, unit="assays", cfg=cfg)

    df = pd.DataFrame(rows)
    has_p = df.p_value.notna()
    df.loc[has_p, "q_value"] = benjamini_hochberg(df.loc[has_p, "p_value"])
    df.to_csv(out / "claims.csv", index=False)

    lines = ["# Verified claims", "", "| claim | model | statement | estimate [95% CI] | control | q | n |", "|---|---|---|---|---|---|---|"]
    for r in df.itertuples():
        ci = f"[{r.ci_low:.3f}, {r.ci_high:.3f}]" if not np.isnan(r.ci_low) else ""
        q = f"{r.q_value:.2g}" if not pd.isna(getattr(r, "q_value", np.nan)) else ""
        ctrl = f"{r.control:.3f}" if not pd.isna(r.control) else ""
        lines.append(f"| {r.claim} | {r.model} | {r.statement} | {r.estimate:.3f} {ci} | {ctrl} | {q} | {r.n} {r.unit} |")
    (out / "claims.md").write_text("\n".join(lines) + "\n")
    log.info(f"{len(df)} claims written to {out}")


if __name__ == "__main__":
    main()
