"""Build per-residue biology labels (the "answer key") for the showcase proteins and a wider panel.

Outputs (data/processed/annotations/):
    proteins.csv        one row per protein: UniProt id, DMS id used, length, split, Neff
    residues.parquet    one row per residue: ss, rsa, buried, pLDDT, conservation, contacts, ...
    dist/<uniprot>.npy  C-beta distance matrix (float16)

Usage:
    python scripts/04_annotate.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from plm_interp.data.msa import conservation, read_a2m
from plm_interp.data.proteingym import dms_reference
from plm_interp.data.structure import load_structure
from plm_interp.paths import ANNOTATIONS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config

log = get_logger("annotate", "annotate.log")


def choose_panel(cfg: dict) -> pd.DataFrame:
    """Showcase proteins plus a deterministic sample of other single-structure proteins."""
    ref = dms_reference()
    full_range = ref.pdb_range == "1-" + ref.seq_len.astype(str)
    ok = ref[full_range & ~ref.pdb_file.str.contains(r"\|") & (ref.seq_len <= cfg["panel"]["max_len"])]
    showcase = ref[ref.DMS_id.isin(cfg["showcase"])]
    # One assay per protein; prefer the one with the most single mutants
    others = (
        ok[~ok.UniProt_ID.isin(showcase.UniProt_ID)]
        .sort_values("DMS_number_single_mutants", ascending=False)
        .drop_duplicates("UniProt_ID")
        .sort_values("UniProt_ID")
    )
    n_other = cfg["panel"]["n_proteins"] - len(showcase)
    others = others.sample(n=min(n_other, len(others)), random_state=cfg["seed"])
    return pd.concat([showcase, others], ignore_index=True)


def main() -> None:
    cfg = load_config()
    out = ensure_dir(ANNOTATIONS_DIR)
    dist_dir = ensure_dir(out / "dist")
    panel = choose_panel(cfg)
    s_cfg = cfg["structure"]

    # Protein-level split for probing: showcase proteins are always held out (test)
    rng = np.random.default_rng(cfg["seed"])
    non_show = panel.UniProt_ID[~panel.DMS_id.isin(cfg["showcase"])].tolist()
    n_test = int(round(cfg["panel"]["test_fraction"] * len(panel))) - len(cfg["showcase"])
    test_ids = set(rng.choice(non_show, size=max(n_test, 0), replace=False)) | set(
        panel.UniProt_ID[panel.DMS_id.isin(cfg["showcase"])]
    )

    residue_rows, protein_rows = [], []
    for row in panel.itertuples():
        st = load_structure(row.pdb_file)
        if st.seq != row.target_seq:
            log.warning(f"{row.UniProt_ID}: structure sequence differs from target, skipped")
            continue
        L = len(st.seq)
        dist = st.distances()
        np.save(dist_dir / f"{row.UniProt_ID}.npy", dist.astype(np.float16))
        contacts = st.contacts(s_cfg["contact_cutoff"], s_cfg["min_seq_sep"])
        confident = st.plddt >= s_cfg["plddt_min"]
        contacts_conf = contacts & confident[:, None] & confident[None, :]

        cons = np.full(L, np.nan)
        neff = np.nan
        try:
            seqs, focus = read_a2m(row.MSA_filename, max_seqs=cfg["msa"]["max_seqs"], seed=cfg["seed"])
            c, neff = conservation(seqs, max_gap_fraction=cfg["msa"]["max_gap_fraction"])
            idx = np.asarray(focus) + int(row.MSA_start) - 1
            valid = (idx >= 0) & (idx < L)
            cons[idx[valid]] = c[valid]
        except KeyError:
            log.warning(f"{row.UniProt_ID}: MSA {row.MSA_filename} not in zip")

        split = "test" if row.UniProt_ID in test_ids else "train"
        protein_rows.append(
            dict(
                uniprot=row.UniProt_ID,
                dms_id=row.DMS_id,
                length=L,
                split=split,
                showcase=row.DMS_id in cfg["showcase"],
                mean_plddt=float(st.plddt.mean()),
                neff=neff,
                n_contacts=int(contacts_conf.sum() // 2),
            )
        )
        residue_rows.append(
            pd.DataFrame(
                dict(
                    uniprot=row.UniProt_ID,
                    pos=np.arange(L),
                    aa=list(st.seq),
                    plddt=st.plddt,
                    confident=confident,
                    ss=st.ss,
                    rsa=st.rsa,
                    buried=st.rsa < s_cfg["buried_rsa"],
                    n_contacts=contacts_conf.sum(1),
                    conservation=cons,
                    disulfide=st.disulfide,
                    split=split,
                )
            )
        )
        log.info(f"{row.UniProt_ID} L={L} split={split} Neff={neff:.0f} pLDDT={st.plddt.mean():.0f}")

    proteins = pd.DataFrame(protein_rows)
    residues = pd.concat(residue_rows, ignore_index=True)
    proteins.to_csv(out / "proteins.csv", index=False)
    residues.to_parquet(out / "residues.parquet")
    log.info(
        f"{len(proteins)} proteins ({(proteins.split == 'test').sum()} test), {len(residues)} residues; "
        f"SS fractions: {residues.ss.value_counts(normalize=True).round(2).to_dict()}"
    )


if __name__ == "__main__":
    main()
