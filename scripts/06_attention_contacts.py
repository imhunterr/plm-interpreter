"""Stage 4: which attention heads look at residues that touch in 3D?

For every protein in the panel, compares each head's attention map with the AF2 contact map
(C-beta < 8 A, |i-j| >= 6, both residues pLDDT >= 70). Repeated on a random-init model.

Output: outputs/runs/attention/head_contacts.csv (one row per model x protein x layer x head)
        outputs/runs/attention/head_summary.csv (averaged over proteins)

Usage:
    python scripts/06_attention_contacts.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from plm_interp.attribution.attention import head_contact_scores
from plm_interp.models.hooks import attention_maps
from plm_interp.models.loader import load_model
from plm_interp.paths import ANNOTATIONS_DIR, RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("attention", "attention.log")


def main() -> None:
    cfg = load_config()
    set_threads(cfg["threads"])
    out = ensure_dir(RUNS_DIR / "attention")
    s_cfg = cfg["structure"]
    residues = pd.read_parquet(ANNOTATIONS_DIR / "residues.parquet")

    rows = []
    for size in cfg["models"]:
        for random_init in (False, True):
            bundle = load_model(size, attn="eager", random_init=random_init, seed=cfg["seed"])
            for uniprot, res in residues.groupby("uniprot", sort=False):
                seq = "".join(res.aa)
                L = len(seq)
                dist = np.load(ANNOTATIONS_DIR / "dist" / f"{uniprot}.npy").astype(np.float32)
                conf = res.confident.to_numpy()
                sep = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :])
                valid = conf[:, None] & conf[None, :] & (sep >= s_cfg["min_seq_sep"])
                contacts = dist < s_cfg["contact_cutoff"]
                if (contacts & valid).sum() < 10:
                    continue
                att = attention_maps(bundle, seq)
                sc = head_contact_scores(att, contacts, valid, s_cfg["min_seq_sep"])
                nl, nh = sc["precision_at_L"].shape
                for layer in range(nl):
                    for head in range(nh):
                        rows.append(
                            dict(
                                model=size,
                                random_init=random_init,
                                uniprot=uniprot,
                                layer=layer,
                                head=head,
                                precision_at_L=sc["precision_at_L"][layer, head],
                                contact_attention=sc["contact_attention"][layer, head],
                                base_rate=sc["base_rate"],
                            )
                        )
            log.info(f"{bundle.name}: done")

    df = pd.DataFrame(rows)
    df["precision_lift"] = df.precision_at_L / df.base_rate
    df["attention_lift"] = df.contact_attention / df.base_rate
    df.to_csv(out / "head_contacts.csv", index=False)
    summary = (
        df.groupby(["model", "random_init", "layer", "head"])
        .agg(
            precision_at_L=("precision_at_L", "mean"),
            precision_lift=("precision_lift", "mean"),
            contact_attention=("contact_attention", "mean"),
            attention_lift=("attention_lift", "mean"),
            base_rate=("base_rate", "mean"),
            n_proteins=("uniprot", "nunique"),
        )
        .reset_index()
    )
    summary.to_csv(out / "head_summary.csv", index=False)
    for (size, rnd), g in summary.groupby(["model", "random_init"]):
        top = g.nlargest(5, "precision_at_L")
        log.info(
            f"{size} random={rnd}: base rate {g.base_rate.mean():.3f}; top heads "
            + ", ".join(f"L{r.layer}H{r.head}={r.precision_at_L:.2f}" for r in top.itertuples())
        )


if __name__ == "__main__":
    main()
