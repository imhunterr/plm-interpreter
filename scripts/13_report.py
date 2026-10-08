"""Collect headline numbers from every stage into outputs/runs/report.md and the README results section.

Usage:
    python scripts/13_report.py
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from plm_interp.paths import EXPLANATIONS_DIR, ROOT, RUNS_DIR
from plm_interp.utils.common import load_config


def _fmt(x, nd=3):
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


def section_baseline() -> list[str]:
    f = RUNS_DIR / "baseline" / "dms_summary.csv"
    if not f.exists():
        return []
    s = pd.read_csv(f).set_index("model")
    out = ["### 1. Baseline: the scorer is correct", "",
           "| model | assays | mean Spearman (ours) | mean Spearman (published) | max per-assay difference | mean AUROC |",
           "|---|---|---|---|---|---|"]
    for m, r in s.iterrows():
        out.append(f"| ESM-2 {m} | {int(r.n_assays)} | {_fmt(r.mean_spearman)} | {_fmt(r.mean_published_spearman)} | {_fmt(r.max_abs_diff, 4)} | {_fmt(r.mean_auc)} |")
    per = RUNS_DIR / "baseline" / "dms_per_assay.csv"
    if per.exists():
        d = pd.read_csv(per)
        off = d[d.abs_diff_vs_published > 0.01]
        n_ok = (d.abs_diff_vs_published <= 0.01).groupby(d.model).sum()
        out += ["", "Per-assay agreement within 0.01 Spearman: " + ", ".join(f"{int(v)}/217 ({m})" for m, v in n_ok.items())
                + ". The exceptions are " + ", ".join(sorted(set(off.DMS_id)))
                + ", proteins longer than ESM-2's 1,022-residue limit, which are scored in windows; our window placement"
                " differs from ProteinGym's for these."]
    c = RUNS_DIR / "baseline" / "clinical_summary.csv"
    if c.exists():
        cl = pd.read_csv(c)
        out += ["", "Clinical variants (pathogenic vs benign, score = -LLR):", "",
                "| model | genes | variants | pooled AUROC | mean per-gene AUROC |", "|---|---|---|---|---|"]
        for r in cl.itertuples():
            out.append(f"| ESM-2 {r.model} | {r.n_proteins} | {r.n_variants} | {_fmt(r.pooled_auc)} | {_fmt(r.mean_per_gene_auc)} (n={r.n_genes_both_classes}) |")
    return out + [""]


def section_probing() -> list[str]:
    f = RUNS_DIR / "probing" / "probe_results.csv"
    if not f.exists():
        return []
    df = pd.read_csv(f)
    aa = df[df.model == "aa_identity"].set_index("task").value
    out = ["### 2. Probing: what each layer encodes (held-out proteins)", "",
           "| task | metric | residue identity only | 8M best layer | 35M best layer | 35M random-init best |", "|---|---|---|---|---|---|"]
    for task, g in df[df.model != "aa_identity"].groupby("task"):
        cells = []
        for m, rnd in (("8M", False), ("35M", False), ("35M", True)):
            gg = g[(g.model == m) & (g.random_init == rnd)]
            if gg.empty:
                cells.append("n/a")
                continue
            b = gg.loc[gg.value.idxmax()]
            cells.append(f"{b.value:.3f} (L{int(b.layer)})")
        metric = g.metric.iloc[0]
        out.append(f"| {task} | {metric} | {_fmt(aa.get(task, np.nan))} | " + " | ".join(cells) + " |")
    return out + [""]


def section_attention() -> list[str]:
    f = RUNS_DIR / "attention" / "head_summary.csv"
    if not f.exists():
        return []
    df = pd.read_csv(f)
    out = ["### 3. Attention heads that track 3D contacts", ""]
    for m in ("8M", "35M"):
        g = df[(df.model == m) & ~df.random_init].nlargest(3, "precision_at_L")
        r = df[(df.model == m) & df.random_init].precision_at_L.max()
        heads = ", ".join(f"L{int(x.layer)}H{int(x['head'])} ({x.precision_at_L:.2f})" for _, x in g.iterrows())
        out.append(f"- ESM-2 {m}: best heads {heads}; contact base rate {g.base_rate.mean():.3f}; best random-init head {r:.3f}.")
    return out + [""]


def section_attribution() -> list[str]:
    f = RUNS_DIR / "attribution" / "variants.csv"
    if not f.exists():
        return []
    df = pd.read_csv(f)
    out = ["### 4. Residue attributions point at 3D neighbours", "",
           "| model | variants | IG converged | contact AUROC | long-range AUROC | seq-distance AUROC | top-10 in contact | IG vs occlusion (Spearman) | random-init AUROC |",
           "|---|---|---|---|---|---|---|---|---|"]
    for m in ("8M", "35M"):
        t = df[df.model == m]
        r = df[df.model == f"{m}-random"]
        if t.empty:
            continue
        occ = t.ig_occ_spearman.mean() if "ig_occ_spearman" in t else np.nan
        out.append(
            f"| ESM-2 {m} | {len(t)} | {t.ig_converged.mean():.0%} | {_fmt(t.ig_auroc_all.mean())} | {_fmt(t.ig_auroc_long_range.mean())} | "
            f"{_fmt(t.ig_auroc_seqdist_baseline.mean())} | {t.ig_top10_contact_frac.mean():.0%} | {_fmt(occ)} | {_fmt(r.ig_auroc_all.mean())} |"
        )
    conv = df[df.ig_converged.astype(bool) & ~df.random_init].groupby("model").ig_auroc_all.mean()
    out += ["", "Robustness: restricted to IG runs whose attributions sum to the score change within 10%, the contact AUROC is "
            + ", ".join(f"{_fmt(v)} ({m})" for m, v in conv.items()) + "."]
    t = df[df.model == "35M"]
    if not t.empty:
        out += ["", "Deletion test (ESM-2 35M): mean |change in score| after masking k residues", "",
                "| k | top attributed | sequence-nearest | random |", "|---|---|---|---|"]
        for k in load_config()["attribution"]["deletion_k"]:
            out.append(f"| {k} | {t[f'del_top_{k}'].mean():.3f} | {t[f'del_near_{k}'].mean():.3f} | {t[f'del_rand_{k}'].mean():.3f} |")
    return out + [""]


def section_ablation() -> list[str]:
    f = RUNS_DIR / "ablation" / "topk.csv"
    if not f.exists():
        return []
    tk = pd.read_csv(f)
    out = ["### 5. Causal tests: ablating heads", "",
           "Drop in held-out Spearman (heads ranked on other positions), mean over showcase proteins:", "",
           "| model | k | top-k heads | random k heads | top-k contact heads |", "|---|---|---|---|---|"]
    for (m, k), g in tk.groupby(["model", "k"]):
        contact = (g.clean_rho_eval - g.contact_rho_eval).mean() if "contact_rho_eval" in g else np.nan
        out.append(f"| ESM-2 {m} | {k} | {(g.clean_rho_eval - g.top_rho_eval).mean():.3f} | "
                   f"{(g.clean_rho_eval - g.random_rho_eval_mean).mean():.3f} | {_fmt(contact)} |")
    comp = RUNS_DIR / "ablation" / "components.csv"
    if comp.exists():
        c = pd.read_csv(comp)
        lay = c[c.kind.isin(["attn_layer", "mlp"])].groupby(["model", "kind"]).d_rho_all.min()
        out += ["", "Largest drop in Spearman from ablating one whole layer, in any single protein: " + "; ".join(f"{m} {k} {v:.3f}" for (m, k), v in lay.items())]
    return out + [""]


def section_lens() -> list[str]:
    f = RUNS_DIR / "logit_lens" / "layers.csv"
    if not f.exists():
        return []
    df = pd.read_csv(f)
    out = ["### 6. Logit lens: when the variant signal appears", "", "Mean Spearman with DMS when decoding from each layer:", ""]
    for m, g in df.groupby("model"):
        means = g.groupby("layer").spearman.mean()
        out.append(f"- ESM-2 {m}: " + " → ".join(f"L{l} {v:.2f}" for l, v in means.items()))
    return out + [""]


def section_cases() -> list[str]:
    f = EXPLANATIONS_DIR / load_config()["main_model"] / "cases.csv"
    if not f.exists():
        return []
    c = pd.read_csv(f)
    out = ["### 7. Case studies", "", "Open `outputs/explanations/35M/index.html` for the interactive cards.", "",
           "| protein | variant | category | model LLR (pct) | assay score (pct) | top-10 residues in 3D contact |", "|---|---|---|---|---|---|"]
    for r in c.itertuples():
        out.append(f"| {r.protein} | {r.mutant} | {r.category.replace('_', ' ')} | {r.llr:+.2f} ({r.llr_percentile:.0f}) | "
                   f"{r.dms_score:+.3f} ({r.dms_percentile:.0f}) | {r.top10_contact_fraction:.0%} |")
    return out + [""]


def main() -> None:
    lines = []
    for fn in (section_baseline, section_probing, section_attention, section_attribution, section_ablation, section_lens, section_cases):
        lines += fn()
    claims = RUNS_DIR / "verification" / "claims.md"
    if claims.exists():
        lines += ["Every claim above, with 95% bootstrap CIs, controls and FDR-corrected permutation tests: "
                  "[`outputs/runs/verification/claims.md`](outputs/runs/verification/claims.md).", ""]
    body = "\n".join(lines)
    (RUNS_DIR / "report.md").write_text("# PLM Interpreter results\n\n" + body)

    readme = ROOT / "README.md"
    text = readme.read_text()
    text = re.sub(r"## Results\n.*?\n## Pipeline", lambda _m: "## Results\n\n<!-- RESULTS -->\n" + body + "\n## Pipeline", text, flags=re.S)
    readme.write_text(text)


if __name__ == "__main__":
    main()
