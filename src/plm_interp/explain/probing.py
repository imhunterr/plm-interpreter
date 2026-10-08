"""Linear probes: which layer's residual stream linearly encodes which residue property."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from plm_interp.models.hooks import capture_residual
from plm_interp.models.loader import AMINO_ACIDS, ModelBundle

# task name -> (label column, kind, uses structure-confident residues only)
TASKS = {
    "secondary_structure": ("ss", "multiclass", True),
    "buried": ("buried", "binary", True),
    "rsa": ("rsa", "regression", True),
    "conservation": ("conservation", "regression", False),
    "n_contacts": ("n_contacts", "regression", True),
}


@torch.no_grad()
def residual_features(bundle: ModelBundle, seq: str) -> np.ndarray:
    """Residual stream for an unmasked sequence: (n_layers + 1, L, D), special tokens dropped."""
    with capture_residual(bundle) as resid:
        bundle.model(input_ids=bundle.encode(seq))
    return torch.stack([r[0, 1:-1] for r in resid]).numpy().astype(np.float32)


def aa_onehot(seq: str) -> np.ndarray:
    idx = {a: i for i, a in enumerate(AMINO_ACIDS)}
    X = np.zeros((len(seq), 20), dtype=np.float32)
    for i, a in enumerate(seq):
        if a in idx:
            X[i, idx[a]] = 1.0
    return X


def _model(kind: str, seed: int):
    if kind == "regression":
        return make_pipeline(StandardScaler(), Ridge(alpha=10.0))
    return make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=2000, random_state=seed))


def _score(kind: str, model, X: np.ndarray, y: np.ndarray) -> dict:
    if kind == "regression":
        pred = model.predict(X)
        return {"metric": "spearman", "value": spearmanr(pred, y)[0]}
    if kind == "binary":
        return {"metric": "auroc", "value": roc_auc_score(y, model.predict_proba(X)[:, 1])}
    return {"metric": "balanced_accuracy", "value": balanced_accuracy_score(y, model.predict(X))}


def run_probe(
    X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray, y_test: np.ndarray, kind: str, seed: int = 0
) -> dict:
    model = _model(kind, seed)
    model.fit(X_train, y_train)
    return _score(kind, model, X_test, y_test)


def task_mask(residues: pd.DataFrame, task: str) -> np.ndarray:
    col, _, needs_conf = TASKS[task]
    m = residues[col].notna().to_numpy().copy()
    if needs_conf:
        m &= residues.confident.to_numpy()
    return m
