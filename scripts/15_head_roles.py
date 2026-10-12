"""What do the heads that matter for variant scores look at?

Variant scores are read at a masked site, so for every head we measure where the masked token's
attention goes, averaged over sampled sites in the ablation proteins:
    self         the masked token itself
    special      <cls> / <eos>
    local        residues 1-4 positions away in sequence
    contact      3D contacts of the site (C-beta < 8 A, >= 6 apart in sequence, confident structure)
    conserved    the most conserved quarter of the protein's positions
    seq_distance attention-weighted mean |i - j| over residues
    entropy      spread of attention over residues (normalised, 0 = one residue, 1 = uniform)

Heads are then grouped by their ablation importance (Stage 6, mean drop in Spearman across proteins
where the model works) and compared with the contact-detector heads from Stage 4.

Outputs (outputs/runs/head_roles/): head_properties.csv (one row per model x head), groups.csv,
correlations.csv

Usage:
    python scripts/15_head_roles.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
from scipy.stats import mannwhitneyu, spearmanr

from plm_interp.data.proteingym import dms_info
from plm_interp.models.loader import load_model
from plm_interp.paths import ANNOTATIONS_DIR, RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("head_roles", "head_roles.log")

N_SITES = 16
TOP_N = 10
MIN_CLEAN_RHO = 0.2
PROPS = ["self", "special", "local", "contact", "conserved", "seq_distance", "entropy"]


@torch.no_grad()
def masked_attention_rows(bundle, seq: str, positions: list[int]) -> np.ndarray:
    """Attention from each masked site to every token: (P, layers, heads, T)."""
    base = bundle.encode(seq)[0]
    ids = base.repeat(len(positions), 1)
    tok = torch.tensor(positions) + 1
    ids[torch.arange(len(positions)), tok] = bundle.mask_id
    att = bundle.model(input_ids=ids, output_attentions=True).attentions  # layers x (P, H, T, T)
    rows = torch.stack([a[torch.arange(len(positions)), :, tok, :] for a in att], dim=1)
    return rows.numpy()


def site_properties(rows: np.ndarray, pos: int, dist: np.ndarray, conf: np.ndarray, conserved: np.ndarray, cfg) -> np.ndarray:
    """rows: (layers, heads, T) for one site -> (layers, heads, len(PROPS))."""
    L = dist.shape[0]
    j = np.arange(L)
    d = np.abs(j - pos)
    res = rows[..., 1:-1]  # residues only
    special = rows[..., 0] + rows[..., -1]
    local = (d >= 1) & (d <= 4)
    contact = (dist[pos] < cfg["contact_cutoff"]) & (d >= cfg["min_seq_sep"]) & conf & conf[pos]
    cons = conserved & (j != pos)
    res_mass = res.sum(-1).clip(1e-9)
    p = res / res_mass[..., None]
    entropy = -(p * np.log(p.clip(1e-12))).sum(-1) / np.log(L)
    return np.stack(
        [
            res[..., pos],
            special,
            res[..., local].sum(-1),
            res[..., contact].sum(-1),
            res[..., cons].sum(-1),
            (p * d).sum(-1),
            entropy,
        ],
        axis=-1,
    )


def main() -> None:
    cfg = load_config()
    set_threads(cfg["threads"])
    out = ensure_dir(RUNS_DIR / "head_roles")
    residues = pd.read_parquet(ANNOTATIONS_DIR / "residues.parquet")
    comp = pd.read_csv(RUNS_DIR / "ablation" / "components.csv")
    contact_summary = pd.read_csv(RUNS_DIR / "attention" / "head_summary.csv")
    proteins = sorted(comp.dms_id.unique())
    rng = np.random.default_rng(cfg["seed"])

    prop_rows = []
    for size in cfg["models"]:
        bundle = load_model(size, attn="eager")
        acc = np.zeros((bundle.n_layers, bundle.n_heads, len(PROPS)))
        n = 0
        for dms_id in proteins:
            info = dms_info(dms_id)
            res = residues[residues.uniprot == info.UniProt_ID].reset_index(drop=True)
            if res.empty:
                continue
            seq = "".join(res.aa)
            dist = np.load(ANNOTATIONS_DIR / "dist" / f"{info.UniProt_ID}.npy").astype(np.float32)
            conf = res.confident.to_numpy()
            cons = res.conservation.to_numpy()
            conserved = cons >= np.nanquantile(cons, 0.75) if np.isfinite(cons).any() else np.zeros(len(seq), bool)
            candidates = np.flatnonzero(conf)
            sites = sorted(rng.choice(candidates, size=min(N_SITES, len(candidates)), replace=False).tolist())
            rows = masked_attention_rows(bundle, seq, sites)
            for k, pos in enumerate(sites):
                acc += site_properties(rows[k], pos, dist, conf, conserved, cfg["structure"])
                n += 1
        acc /= n
        for layer in range(bundle.n_layers):
            for head in range(bundle.n_heads):
                prop_rows.append(dict(model=size, layer=layer, head=head, **dict(zip(PROPS, acc[layer, head]))))
        log.info(f"{size}: {n} masked sites over {len(proteins)} proteins")

    props = pd.DataFrame(prop_rows)

    # Importance = mean drop in Spearman when the head is ablated, over proteins where the model works
    heads = comp[(comp.kind == "head") & (comp.clean_rho_all >= MIN_CLEAN_RHO)]
    importance = (-heads.groupby(["model", "layer", "head"]).d_rho_all.mean()).rename("importance")
    n_prot = heads.groupby("model").dms_id.nunique()
    contact = contact_summary[~contact_summary.random_init].set_index(["model", "layer", "head"]).precision_at_L.rename("contact_precision")
    props = props.join(importance, on=["model", "layer", "head"]).join(contact, on=["model", "layer", "head"])
    props.to_csv(out / "head_properties.csv", index=False)

    group_rows, corr_rows = [], []
    for size, g in props.groupby("model"):
        g = g.copy()
        mut = g.nlargest(TOP_N, "importance").index
        con = g.nlargest(TOP_N, "contact_precision").index
        g["group"] = "other"
        g.loc[con, "group"] = "contact heads"
        g.loc[mut, "group"] = "mutation heads"
        overlap = len(set(mut) & set(con))
        for grp, gg in g.groupby("group"):
            group_rows.append(dict(model=size, group=grp, n_heads=len(gg), mean_layer=gg.layer.mean(),
                                   **{p: gg[p].mean() for p in PROPS + ["importance", "contact_precision"]}))
        others = g[g.group == "other"]
        for p in PROPS + ["contact_precision"]:
            rho, pval = spearmanr(g.importance, g[p])
            mw = mannwhitneyu(g.loc[mut, p], others[p]).pvalue
            corr_rows.append(dict(model=size, property=p, spearman_with_importance=rho, p_value=pval,
                                  mutation_heads_mean=g.loc[mut, p].mean(), other_heads_mean=others[p].mean(),
                                  mannwhitney_p_vs_other=mw))
        log.info(f"{size}: {overlap}/{TOP_N} mutation heads are also top contact heads; "
                 f"importance from {n_prot.get(size, 0)} proteins with clean Spearman >= {MIN_CLEAN_RHO}")
    groups = pd.DataFrame(group_rows)
    corrs = pd.DataFrame(corr_rows)
    groups.to_csv(out / "groups.csv", index=False)
    corrs.to_csv(out / "correlations.csv", index=False)
    log.info("\n" + groups.round(3).to_string(index=False))
    log.info("\n" + corrs.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
