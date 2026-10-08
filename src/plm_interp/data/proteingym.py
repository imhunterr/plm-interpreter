"""Loaders for the ProteinGym v1.2 files under data/raw/proteingym."""

from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from plm_interp.paths import BENCHMARK_DIR, CLINICAL_DIR, DMS_DIR, REFERENCE_DIR

MUTANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


@lru_cache(maxsize=1)
def dms_reference() -> pd.DataFrame:
    return pd.read_csv(REFERENCE_DIR / "DMS_substitutions.csv")


@lru_cache(maxsize=1)
def clinical_reference() -> pd.DataFrame:
    return pd.read_csv(REFERENCE_DIR / "clinical_substitutions.csv")


def dms_info(dms_id: str) -> pd.Series:
    ref = dms_reference()
    row = ref[ref.DMS_id == dms_id]
    if row.empty:
        raise KeyError(f"unknown DMS_id {dms_id!r}")
    return row.iloc[0]


def load_dms(dms_id: str) -> pd.DataFrame:
    """One DMS assay: columns mutant, mutated_sequence, DMS_score, DMS_score_bin."""
    return pd.read_csv(DMS_DIR / dms_info(dms_id).DMS_filename)


def load_clinical(protein_id: str) -> pd.DataFrame:
    """Clinical variants for one RefSeq protein; adds label (1 = pathogenic, 0 = benign)."""
    df = pd.read_csv(CLINICAL_DIR / f"{protein_id}.csv", index_col=0)
    df["label"] = (df.DMS_bin_score == "Pathogenic").astype(int)
    return df


def parse_mutant(mutant: str) -> list[tuple[str, int, str]]:
    """'A12G:L30P' -> [('A', 11, 'G'), ('L', 29, 'P')] with 0-based positions."""
    out = []
    for m in mutant.split(":"):
        match = MUTANT_RE.match(m)
        if match is None:
            raise ValueError(f"cannot parse mutation {m!r}")
        wt, pos, mt = match.groups()
        out.append((wt, int(pos) - 1, mt))
    return out


@lru_cache(maxsize=2)
def published_benchmark(metric: str = "Spearman") -> pd.DataFrame:
    """ProteinGym's published per-assay zero-shot scores (downloaded from their GitHub)."""
    df = pd.read_csv(BENCHMARK_DIR / f"DMS_substitutions_{metric}_DMS_level.csv")
    return df.rename(columns={"DMS ID": "DMS_id"})
