"""Stage 3: linear probes per layer for secondary structure, burial, RSA, conservation, contact count.

Trained on the panel's train proteins, evaluated on held-out proteins (protein-level split).
Controls: the same architecture with random weights, and a probe on amino-acid identity alone.

Output: outputs/runs/probing/probe_results.csv

Usage:
    python scripts/05_probing.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from plm_interp.explain.probing import TASKS, aa_onehot, residual_features, run_probe, task_mask
from plm_interp.models.loader import load_model
from plm_interp.paths import ANNOTATIONS_DIR, RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, seed_everything, set_threads

log = get_logger("probing", "probing.log")


def main() -> None:
    cfg = load_config()
    set_threads(cfg["threads"])
    seed_everything(cfg["seed"])
    out = ensure_dir(RUNS_DIR / "probing")

    residues = pd.read_parquet(ANNOTATIONS_DIR / "residues.parquet")
    proteins = residues.groupby("uniprot", sort=False).aa.apply("".join)
    is_train = (residues.split == "train").to_numpy()

    variants = [(size, False) for size in cfg["models"]] + [(size, True) for size in cfg["models"]]
    rows = []

    # Amino-acid identity baseline: how much of each property is explained by residue type alone
    X_aa = np.concatenate([aa_onehot(s) for s in proteins])
    for task, (col, kind, _) in TASKS.items():
        m = task_mask(residues, task)
        y = residues[col].to_numpy()
        tr, te = m & is_train, m & ~is_train
        r = run_probe(X_aa[tr], y[tr], X_aa[te], y[te], kind, cfg["seed"])
        rows.append(dict(model="aa_identity", random_init=False, layer=-1, task=task, **r))
        log.info(f"aa_identity {task}: {r['value']:.3f}")

    for size, random_init in variants:
        bundle = load_model(size, random_init=random_init, seed=cfg["seed"])
        feats = [residual_features(bundle, s) for s in proteins]  # list of (layers+1, L, D)
        n_layers = feats[0].shape[0]
        for layer in range(n_layers):
            X = np.concatenate([f[layer] for f in feats])
            for task, (col, kind, _) in TASKS.items():
                m = task_mask(residues, task)
                y = residues[col].to_numpy()
                tr, te = m & is_train, m & ~is_train
                r = run_probe(X[tr], y[tr], X[te], y[te], kind, cfg["seed"])
                rows.append(dict(model=size, random_init=random_init, layer=layer, task=task, **r))
            log.info(f"{bundle.name} layer {layer}: " + ", ".join(f"{t}={x['value']:.3f}" for t, x in zip(TASKS, rows[-len(TASKS):])))
        del feats

    df = pd.DataFrame(rows)
    df.to_csv(out / "probe_results.csv", index=False)
    best = df[~df.random_init & (df.layer >= 0)].loc[lambda d: d.groupby(["model", "task"]).value.idxmax()]
    log.info("\nBest layer per task:\n" + best[["model", "task", "layer", "metric", "value"]].to_string(index=False))


if __name__ == "__main__":
    main()
