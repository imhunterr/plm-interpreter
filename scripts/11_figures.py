"""Paper figures from the Stage 1-8 outputs. Each figure is skipped if its inputs are missing.

Output: outputs/figures/*.png and *.pdf

Usage:
    python scripts/11_figures.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plm_interp.paths import FIGURES_DIR, RUNS_DIR
from plm_interp.utils.common import get_logger, load_config
from plm_interp.viz import style
from plm_interp.viz.style import AQUA, DIVERGING, GRAY, LIGHT_GRAY, MODEL_COLORS, SEQUENTIAL, TEXT_2, VIOLET

log = get_logger("figures", "figures.log")

TASK_TITLES = {
    "secondary_structure": "Secondary structure (balanced acc.)",
    "buried": "Buried vs exposed (AUROC)",
    "rsa": "Solvent accessibility (Spearman)",
    "conservation": "Conservation (Spearman)",
    "n_contacts": "Contact count (Spearman)",
}


def fig_baseline() -> None:
    f = RUNS_DIR / "baseline" / "dms_per_assay.csv"
    if not f.exists():
        return
    df = pd.read_csv(f)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ax = axes[0]
    lim = [-0.2, 0.8]
    ax.plot(lim, lim, color=LIGHT_GRAY, lw=1, zorder=0)
    for m in ["8M", "35M"]:
        g = df[df.model == m]
        ax.scatter(g.published_spearman, g.spearman, s=12, color=MODEL_COLORS[m], label=f"ESM-2 {m}", alpha=0.8, lw=0)
    ax.set(xlim=lim, ylim=lim, xlabel="Published Spearman (ProteinGym)", ylabel="Our Spearman",
           title="Our scorer reproduces ProteinGym")
    ax.legend(loc="upper left")
    ax = axes[1]
    order = df.groupby("selection_type").spearman.mean().sort_values().index
    y = np.arange(len(order))
    for i, m in enumerate(["8M", "35M"]):
        g = df[df.model == m].groupby("selection_type").spearman.mean().reindex(order)
        ax.barh(y + (i - 0.5) * 0.38, g.values, height=0.36, color=MODEL_COLORS[m], label=f"ESM-2 {m}")
    ax.set_yticks(y, order)
    ax.set(xlabel="Mean Spearman with DMS", title="Zero-shot accuracy by assay type")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right")
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig1_baseline")


def fig_probing() -> None:
    f = RUNS_DIR / "probing" / "probe_results.csv"
    if not f.exists():
        return
    df = pd.read_csv(f)
    aa = df[df.model == "aa_identity"].set_index("task").value
    tasks = list(TASK_TITLES)
    fig, axes = plt.subplots(1, len(tasks), figsize=(12, 2.8), sharex=True)
    for ax, task in zip(axes, tasks):
        for m in ["8M", "35M"]:
            for rnd, ls, alpha in ((False, "-", 1.0), (True, "--", 0.6)):
                g = df[(df.model == m) & (df.random_init == rnd) & (df.task == task)].sort_values("layer")
                if g.empty:
                    continue
                depth = g.layer / g.layer.max()
                ax.plot(depth, g.value, ls=ls, color=MODEL_COLORS[m] if not rnd else GRAY, alpha=alpha,
                        lw=2 if not rnd else 1.2, marker="o" if not rnd else None,
                        label=f"{m}" + (" random init" if rnd else ""))
        if task in aa:
            ax.axhline(aa[task], color=TEXT_2, lw=1, ls=":", label="residue identity only")
        ax.set_title(TASK_TITLES[task], fontsize=8.5)
        ax.set_xlabel("Relative depth (0 = embeddings)")
    axes[0].set_ylabel("Held-out probe score")
    handles, labels = axes[0].get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    fig.legend(uniq.values(), uniq.keys(), loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.08))
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig2_probing")


def _layer_head_heatmap(ax, mat, cmap, vmin, vmax, title, cbar_label, fig):
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto", interpolation="nearest")
    ax.set(xlabel="Head", ylabel="Layer", title=title)
    ax.set_xticks(range(0, mat.shape[1], 2))
    ax.set_yticks(range(mat.shape[0]))
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(cbar_label, color=TEXT_2)
    cb.outline.set_visible(False)


def fig_attention() -> None:
    f = RUNS_DIR / "attention" / "head_summary.csv"
    if not f.exists():
        return
    df = pd.read_csv(f)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), gridspec_kw={"width_ratios": [1, 1.6, 1]})
    vmax = df.precision_at_L.max()
    for ax, m in zip(axes[:2], ["8M", "35M"]):
        g = df[(df.model == m) & ~df.random_init]
        mat = g.pivot(index="layer", columns="head", values="precision_at_L").to_numpy()
        _layer_head_heatmap(ax, mat, SEQUENTIAL, 0, vmax, f"ESM-2 {m}: contact precision@L per head", "precision@L", fig)
    ax = axes[2]
    for m in ["8M", "35M"]:
        for rnd in (False, True):
            g = df[(df.model == m) & (df.random_init == rnd)]
            vals = np.sort(g.precision_at_L.to_numpy())[::-1]
            ax.plot(np.arange(1, len(vals) + 1), vals, color=GRAY if rnd else MODEL_COLORS[m],
                    ls="--" if rnd else "-", lw=1.2 if rnd else 2, label=f"{m}" + (" random init" if rnd else ""))
    ax.axhline(df.base_rate.mean(), color=TEXT_2, ls=":", lw=1, label="contact base rate")
    ax.set(xscale="log", xlabel="Head rank", ylabel="precision@L", title="Few heads are contact detectors")
    ax.legend()
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig3_attention_contacts")


def fig_attribution(cfg) -> None:
    f = RUNS_DIR / "attribution" / "variants.csv"
    if not f.exists():
        return
    df = pd.read_csv(f)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.3))
    ax = axes[0]
    groups, labels, colors = [], [], []
    for m in ["8M", "35M"]:
        t = df[df.model == m]
        r = df[df.model == f"{m}-random"]
        groups += [t.ig_auroc_seqdist_baseline.dropna(), t.ig_auroc_all.dropna(), t.ig_auroc_long_range.dropna(), r.ig_auroc_all.dropna()]
        labels += [f"{m}\nseq. distance", f"{m}\nIG", f"{m}\nIG long-range", f"{m}\nIG random init"]
        colors += [LIGHT_GRAY, MODEL_COLORS[m], MODEL_COLORS[m], GRAY]
    bp = ax.boxplot(groups, patch_artist=True, widths=0.6, showfliers=False, medianprops=dict(color="white", lw=1.5))
    for patch, c in zip(bp["boxes"], colors):
        patch.set(facecolor=c, edgecolor=c)
    for el in ("whiskers", "caps"):
        for line, c in zip(bp[el], np.repeat(colors, 2)):
            line.set(color=c)
    ax.axhline(0.5, color=TEXT_2, ls=":", lw=1)
    ax.set_xticks(range(1, len(labels) + 1), labels, fontsize=6.5)
    ax.set(ylabel="AUROC: attribution vs 3D contact", title="Attributions find 3D neighbours")
    ax.grid(axis="x", visible=False)

    ks = cfg["attribution"]["deletion_k"]
    for ax, m in zip(axes[1:], ["8M", "35M"]):
        t = df[df.model == m]
        for prefix, label, color in (("del_top", "top attributed", MODEL_COLORS[m]), ("del_near", "sequence-nearest", AQUA), ("del_rand", "random", GRAY)):
            means = [t[f"{prefix}_{k}"].mean() for k in ks]
            sems = [t[f"{prefix}_{k}"].std() / np.sqrt(len(t)) for k in ks]
            ax.errorbar(ks, means, yerr=sems, color=color, marker="o", lw=2, capsize=0, label=label)
        ax.set(xscale="log", xlabel="Residues masked (k)", ylabel="|change in variant score|",
               title=f"Deletion test, ESM-2 {m}")
        ax.set_xticks(ks, [str(k) for k in ks])
        ax.legend()
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig4_attribution")


def fig_ablation() -> None:
    f = RUNS_DIR / "ablation" / "components.csv"
    if not f.exists():
        return
    df = pd.read_csv(f)
    heads = df[df.kind == "head"].groupby(["model", "layer", "head"]).d_rho_all.mean().reset_index()
    lim = max(0.02, heads.d_rho_all.abs().max())
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), gridspec_kw={"width_ratios": [1, 1.6, 1]})
    for ax, m in zip(axes[:2], ["8M", "35M"]):
        g = heads[heads.model == m]
        if g.empty:
            ax.set_visible(False)
            continue
        mat = g.pivot(index="layer", columns="head", values="d_rho_all").to_numpy()
        _layer_head_heatmap(ax, mat, DIVERGING.reversed(), -lim, lim, f"ESM-2 {m}: Spearman change when head ablated",
                            "change in Spearman", fig)
    ax = axes[2]
    tk_f = RUNS_DIR / "ablation" / "topk.csv"
    if tk_f.exists():
        tk = pd.read_csv(tk_f)
        for m in ["8M", "35M"]:
            g = tk[tk.model == m].groupby("k").mean(numeric_only=True)
            if g.empty:
                continue
            ax.plot(g.index, g.clean_rho_eval - g.top_rho_eval, color=MODEL_COLORS[m], marker="o", label=f"{m} top-k heads")
            ax.plot(g.index, g.clean_rho_eval - g.random_rho_eval_mean, color=MODEL_COLORS[m], ls="--", lw=1.2, label=f"{m} random k heads")
            if "contact_rho_eval" in g:
                ax.plot(g.index, g.clean_rho_eval - g.contact_rho_eval, color=AQUA if m == "35M" else VIOLET, ls=":", lw=1.5, label=f"{m} top-k contact heads")
        ax.set(xscale="log", xlabel="Heads ablated (k)", ylabel="Drop in held-out Spearman", title="Few heads carry the signal")
        ax.set_xticks(sorted(tk.k.unique()), [str(k) for k in sorted(tk.k.unique())])
        ax.legend(fontsize=7)
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig5_ablation")

    # Layer-level: attention layer vs MLP
    lay = df[df.kind.isin(["attn_layer", "mlp"])].groupby(["model", "kind", "layer"]).d_rho_all.mean().reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(8, 2.8), sharey=True)
    for ax, m in zip(axes, ["8M", "35M"]):
        for kind, color, label in (("attn_layer", MODEL_COLORS[m], "whole attention layer"), ("mlp", AQUA, "MLP")):
            g = lay[(lay.model == m) & (lay.kind == kind)]
            ax.plot(g.layer, g.d_rho_all, color=color, marker="o", label=label)
        ax.axhline(0, color=TEXT_2, lw=0.8)
        ax.set(xlabel="Layer", title=f"ESM-2 {m}: layer ablations")
        ax.legend()
    axes[0].set_ylabel("Change in Spearman")
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig5b_layer_ablation")



def fig_logit_lens() -> None:
    f = RUNS_DIR / "logit_lens" / "layers.csv"
    if not f.exists():
        return
    df = pd.read_csv(f)
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    for m in ["8M", "35M"]:
        g = df[df.model == m]
        for _, a in g.groupby("dms_id"):
            axes[0].plot(a.rel_depth, a.spearman, color=MODEL_COLORS[m], alpha=0.12, lw=0.8)
        mean = g.groupby("rel_depth").spearman.mean()
        axes[0].plot(mean.index, mean.values, color=MODEL_COLORS[m], marker="o", label=f"ESM-2 {m} (mean)")
        acc = g.groupby("rel_depth").wt_top1.mean()
        axes[1].plot(acc.index, acc.values, color=MODEL_COLORS[m], marker="o", label=f"ESM-2 {m}")
    axes[0].set(xlabel="Relative depth", ylabel="Spearman with DMS", title="Logit lens: when the variant signal appears")
    axes[1].set(xlabel="Relative depth", ylabel="Wild-type residue top-1 accuracy", title="Logit lens: recovering the masked residue")
    for ax in axes:
        ax.legend()
    fig.tight_layout()
    style.save(fig, FIGURES_DIR / "fig6_logit_lens")


def main() -> None:
    cfg = load_config()
    style.apply()
    for fn in (fig_baseline, fig_probing, fig_attention, lambda: fig_attribution(cfg), fig_ablation, fig_logit_lens):
        try:
            fn()
        except Exception as e:  # keep going so one missing input does not block the rest
            log.warning(f"{getattr(fn, '__name__', 'figure')} failed: {e!r}")
    log.info(f"figures in {FIGURES_DIR}")


if __name__ == "__main__":
    main()
