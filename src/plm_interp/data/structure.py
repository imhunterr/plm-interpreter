"""Per-residue structural labels from ProteinGym's AlphaFold2 models.

Labels: pLDDT, DSSP 3-state secondary structure (pydssp), relative solvent accessibility
(Shrake-Rupley on heavy atoms / Tien et al. 2013 theoretical maxima), C-beta contact map,
and disulfide bonds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pydssp
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.SASA import ShrakeRupley

from plm_interp.paths import STRUCTURES_DIR

# Theoretical max ASA (A^2), Tien et al. 2013
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0, "Q": 225.0,
    "G": 104.0, "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0,
    "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


@dataclass
class Structure:
    seq: str
    plddt: np.ndarray  # (L,)
    backbone: np.ndarray  # (L, 4, 3) N, CA, C, O
    cb: np.ndarray  # (L, 3) C-beta, C-alpha for glycine
    rsa: np.ndarray  # (L,)
    ss: np.ndarray  # (L,) of 'H', 'E', '-'
    disulfide: np.ndarray  # (L,) bool

    def distances(self) -> np.ndarray:
        d = self.cb[:, None, :] - self.cb[None, :, :]
        return np.sqrt((d**2).sum(-1))

    def contacts(self, cutoff: float = 8.0, min_sep: int = 6) -> np.ndarray:
        L = len(self.seq)
        sep = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :])
        return (self.distances() < cutoff) & (sep >= min_sep)


def load_structure(pdb_file: str) -> Structure:
    model = PDBParser(QUIET=True).get_structure("s", STRUCTURES_DIR / pdb_file)[0]
    chain = next(iter(model))
    residues = [r for r in chain if r.id[0] == " "]

    # Heavy atoms only, so SASA matches the theoretical maxima
    for res in residues:
        for atom in [a for a in res if a.element == "H"]:
            res.detach_child(atom.id)
    ShrakeRupley().compute(model, level="R")

    seq, plddt, backbone, cb, rsa = [], [], [], [], []
    sg = []
    for res in residues:
        aa = index_to_one(three_to_index(res.get_resname()))
        seq.append(aa)
        plddt.append(res["CA"].get_bfactor())
        backbone.append([res[n].coord for n in ("N", "CA", "C", "O")])
        cb.append(res["CB"].coord if "CB" in res else res["CA"].coord)
        rsa.append(min(1.0, res.sasa / MAX_ASA[aa]))
        sg.append(res["SG"].coord if aa == "C" and "SG" in res else None)

    backbone = np.asarray(backbone, dtype=np.float32)
    ss = pydssp.assign(backbone, out_type="c3")

    disulfide = np.zeros(len(seq), dtype=bool)
    cys = [i for i, c in enumerate(sg) if c is not None]
    for a in cys:
        for b in cys:
            if a < b and np.linalg.norm(sg[a] - sg[b]) < 2.5:
                disulfide[[a, b]] = True

    return Structure(
        seq="".join(seq),
        plddt=np.asarray(plddt),
        backbone=backbone,
        cb=np.asarray(cb, dtype=np.float32),
        rsa=np.asarray(rsa),
        ss=np.asarray(ss),
        disulfide=disulfide,
    )
