"""Per-position conservation from ProteinGym's a2m alignments, read straight from the zip.

Conservation = 1 - H / log(21), where H is the sequence-weighted Shannon entropy of a
focus column (20 amino acids + gap). Sequences are weighted by 1 / (number of neighbours
at >= 80% identity), the standard EVcouplings/EVE scheme.
"""

from __future__ import annotations

import io
import zipfile

import numpy as np

from plm_interp.paths import MSA_ZIP

ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
_LUT = np.full(256, 20, dtype=np.int8)  # unknown letters count as gaps
for _i, _c in enumerate(ALPHABET):
    _LUT[ord(_c)] = _i


def read_a2m(msa_filename: str, max_seqs: int = 5000, seed: int = 0) -> tuple[list[str], list[int]]:
    """Focus-column sequences (insertion columns dropped), query first, plus focus positions.

    Keeps the query plus a uniform random sample (reservoir) of up to max_seqs - 1 others.
    Focus positions are 0-based offsets into the query's residues: in a2m, the query's
    uppercase letters are the match (focus) columns and lowercase letters are not.
    """
    rng = np.random.default_rng(seed)
    with zipfile.ZipFile(MSA_ZIP) as z, z.open(f"DMS_msa_files/{msa_filename}") as raw:
        handle = io.TextIOWrapper(raw)
        query, reservoir, seen, current = None, [], 0, []
        focus: list[int] = []

        def flush():
            nonlocal query, seen
            if not current:
                return
            raw_seq = "".join(current)
            s = "".join(c for c in raw_seq if not (c.islower() or c == "."))
            if query is None:
                query = s
                residues = [c for c in raw_seq if c not in "-."]
                focus.extend(i for i, c in enumerate(residues) if c.isupper())
                return
            seen += 1
            if len(reservoir) < max_seqs - 1:
                reservoir.append(s)
            else:
                j = rng.integers(seen)
                if j < max_seqs - 1:
                    reservoir[j] = s

        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                flush()
                current = []
            elif line:
                current.append(line)
        flush()
    return [query] + reservoir, focus


def conservation(seqs: list[str], identity: float = 0.8, max_gap_fraction: float = 0.5) -> tuple[np.ndarray, float]:
    """Return (conservation per focus column, effective number of sequences)."""
    X = np.stack([_LUT[np.frombuffer(s.encode(), dtype=np.uint8)] for s in seqs])  # (N, C)
    # Drop gappy sequences (fragments) except the query
    keep = (X == 20).mean(1) <= max_gap_fraction
    keep[0] = True
    X = X[keep]
    N, C = X.shape

    # One-hot over the 20 amino acids only, so dot products count matching non-gap residues
    oh_nogap = np.eye(21, dtype=np.float32)[X][..., :20].reshape(N, C * 20)
    nongap = (X != 20).sum(1).astype(np.float32)

    # Identity relative to each row sequence's non-gap length, computed in chunks
    weights = np.zeros(N, dtype=np.float64)
    for start in range(0, N, 1024):
        block = oh_nogap[start : start + 1024] @ oh_nogap.T  # matching non-gap residues
        ident = block / np.maximum(nongap[start : start + 1024, None], 1.0)
        weights[start : start + 1024] = 1.0 / (ident >= identity).sum(1)

    counts = np.zeros((C, 21))
    for a in range(21):
        counts[:, a] = ((X == a) * weights[:, None]).sum(0)
    p = counts / counts.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
    return 1.0 - H / np.log(21), float(weights.sum())
