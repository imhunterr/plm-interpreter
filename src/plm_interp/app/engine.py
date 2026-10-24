"""Model-side logic for the demo. Every function reuses the tested pipeline code."""

from __future__ import annotations

import hashlib
import threading
from functools import lru_cache

import numpy as np
import pandas as pd
import torch

from plm_interp.attribution.methods import integrated_gradients
from plm_interp.baseline.scoring import AA_INDEX, masked_logprobs
from plm_interp.data.proteingym import dms_info, load_dms
from plm_interp.explain.ablation import head_means
from plm_interp.explain.logit_lens import layerwise_logprobs
from plm_interp.models.loader import AMINO_ACIDS, MAX_RESIDUES, ModelBundle, load_model
from plm_interp.paths import ANNOTATIONS_DIR, INTERIM_DIR, TABLES_DIR, STRUCTURES_DIR

# One inference at a time: torch already uses every core, parallel requests would only thrash
_LOCK = threading.Lock()

DEMO_PROTEINS = [
    dict(key="hras", dms_id="RASH_HUMAN_Bandaru_2017", name="HRas", gene="HRAS",
         blurb="A molecular switch that controls cell growth. Mutations at G12, G13 and Q61 are among the most common in cancer.",
         suggest=[(12, "V"), (10, "I"), (61, "L")]),
    dict(key="p53", dms_id="P53_HUMAN_Kotler_2018", name="p53", gene="TP53",
         blurb="The 'guardian of the genome', a tumour suppressor. R175H, R248W and R273H are classic cancer hotspots.",
         suggest=[(175, "H"), (248, "W"), (273, "H")]),
    dict(key="tpk1", dms_id="TPK1_HUMAN_Weile_2017", name="TPK1", gene="TPK1",
         blurb="Thiamin pyrophosphokinase, which activates vitamin B1. Defects cause a rare neurological disease.",
         suggest=[(73, "R"), (23, "R")]),
    dict(key="otc", dms_id="OTC_HUMAN_Lo_2023", name="OTC", gene="OTC",
         blurb="Ornithine transcarbamylase, a urea-cycle enzyme. Its deficiency is the most common urea-cycle disorder.",
         suggest=[(303, "W"), (267, "G")]),
    dict(key="hem3", dms_id="HEM3_HUMAN_Loggerenberg_2023", name="HMBS", gene="HMBS",
         blurb="Porphobilinogen deaminase, a step in making heme. Mutations cause acute intermittent porphyria.",
         suggest=[(42, "D"), (22, "W")]),
    dict(key="phot", dms_id="PHOT_CHLRE_Chen_2023", name="Phototropin LOV", gene="PHOT",
         blurb="A light-sensing domain from green algae. ESM-2 agrees with this lab experiment better than almost any other (Spearman 0.82).",
         suggest=[]),
    dict(key="dnaj", dms_id="DNJA1_HUMAN_Tsuboyama_2023_2LO1", name="DNAJA1 domain", gene="DNAJA1",
         blurb="A small 65-residue domain of a human chaperone. Its stability was measured for every possible mutation.",
         suggest=[]),
]
_BY_KEY = {p["key"]: p for p in DEMO_PROTEINS}


class DemoError(ValueError):
    """A problem with the request that the user can fix."""


@lru_cache(maxsize=2)
def bundle(size: str) -> ModelBundle:
    if size not in ("8M", "35M"):
        raise DemoError("Model must be 8M or 35M.")
    return load_model(size)


def clean_sequence(seq: str) -> str:
    seq = "".join(seq.split()).upper()
    bad = sorted(set(seq) - set(AMINO_ACIDS))
    if bad:
        raise DemoError(f"The sequence contains letters that are not standard amino acids: {', '.join(bad)}.")
    if not 10 <= len(seq) <= MAX_RESIDUES:
        raise DemoError(f"Use a sequence of 10 to {MAX_RESIDUES} amino acids (this one has {len(seq)}).")
    return seq


def _check_position(seq: str, position: int) -> None:
    if not 0 <= position < len(seq):
        raise DemoError(f"Position must be between 1 and {len(seq)}.")


# ---------------------------------------------------------------- proteins

def protein_list() -> list[dict]:
    out = []
    for p in DEMO_PROTEINS:
        info = dms_info(p["dms_id"])
        out.append(dict(key=p["key"], name=p["name"], gene=p["gene"], blurb=p["blurb"], length=int(info.seq_len)))
    return out


@lru_cache(maxsize=16)
def protein_detail(key: str) -> dict:
    p = _BY_KEY.get(key)
    if p is None:
        raise DemoError(f"Unknown protein {key!r}.")
    info = dms_info(p["dms_id"])
    seq = info.target_seq

    residues = None
    ann = ANNOTATIONS_DIR / "residues.parquet"
    if ann.exists():
        r = pd.read_parquet(ann)
        r = r[r.uniprot == info.UniProt_ID]
        if len(r) == len(seq):
            residues = [
                dict(ss=str(x.ss), rsa=round(float(x.rsa), 3), plddt=round(float(x.plddt), 1),
                     conservation=None if pd.isna(x.conservation) else round(float(x.conservation), 3))
                for x in r.itertuples()
            ]

    pdb = None
    pdb_path = STRUCTURES_DIR / str(info.pdb_file)
    if "|" not in str(info.pdb_file) and pdb_path.exists():
        pdb = pdb_path.read_text()

    dms = load_dms(p["dms_id"])
    dms = dms[~dms.mutant.str.contains(":")]
    pct = dms.DMS_score.rank(pct=True)
    lab = {m: [round(float(s), 4), round(float(q) * 100, 1)] for m, s, q in zip(dms.mutant, dms.DMS_score, pct)}

    return dict(
        key=key, name=p["name"], gene=p["gene"], blurb=p["blurb"], uniprot=info.UniProt_ID,
        dms_id=p["dms_id"], assay=f"{info.first_author} {int(info.year)}", assay_type=info.coarse_selection_type,
        sequence=seq, residues=residues, pdb=pdb, lab=lab,
        suggest=[dict(position=pos, mutant=mt) for pos, mt in p["suggest"]],
    )


def tokens(seq: str, masked: int | None = None) -> list[dict]:
    """The model's view of the sequence: token strings and vocabulary ids."""
    b = bundle("8M")  # same vocabulary for every ESM-2 size
    ids = b.encode(seq)[0]
    if masked is not None:
        ids[masked + 1] = b.mask_id
    toks = b.tokenizer.convert_ids_to_tokens(ids.tolist())
    return [dict(token=t, id=int(i)) for t, i in zip(toks, ids.tolist())]


def _distances(key: str | None, seq: str) -> np.ndarray | None:
    if not key or key not in _BY_KEY:
        return None
    info = dms_info(_BY_KEY[key]["dms_id"])
    if info.target_seq != seq:
        return None
    path = ANNOTATIONS_DIR / "dist" / f"{info.UniProt_ID}.npy"
    return np.load(path).astype(np.float32) if path.exists() else None


# ---------------------------------------------------------------- predictions

def predict(seq: str, position: int, size: str) -> dict:
    """Mask one position and return the model's probabilities for all 20 amino acids."""
    seq = clean_sequence(seq)
    _check_position(seq, position)
    with _LOCK:
        lp = masked_logprobs(bundle(size), seq, [position]).loc[position]
    probs = np.exp(lp.to_numpy())
    order = np.argsort(-probs)
    wt = seq[position]
    entropy = float(-(probs * np.log(probs.clip(1e-12))).sum() / np.log(20))
    return dict(
        position=position, wt=wt, model=size,
        probs=[dict(aa=AMINO_ACIDS[i], p=round(float(probs[i]), 5), logp=round(float(lp.iloc[i]), 4)) for i in order],
        top=AMINO_ACIDS[order[0]], wt_rank=int(list(order).index(AA_INDEX[wt])) + 1,
        wt_prob=round(float(probs[AA_INDEX[wt]]), 5), uncertainty=round(entropy, 3),
        tokens=tokens(seq, position),
    )


def landscape(seq: str, size: str) -> dict:
    """LLR of every possible single substitution (positions x 20), cached on disk."""
    seq = clean_sequence(seq)
    cache_dir = INTERIM_DIR / "app_cache" / size
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{hashlib.sha1(seq.encode()).hexdigest()}.npy"
    if path.exists():
        lp = np.load(path)
    else:
        with _LOCK:
            lp = masked_logprobs(bundle(size), seq).to_numpy().astype(np.float32)
        np.save(path, lp)
    wt_idx = np.array([AA_INDEX[a] for a in seq])
    llr = lp - lp[np.arange(len(seq)), wt_idx][:, None]
    return dict(model=size, amino_acids=AMINO_ACIDS, llr=np.round(llr, 3).tolist())


def score(seq: str, position: int, mutant: str, size: str) -> dict:
    seq = clean_sequence(seq)
    _check_position(seq, position)
    if mutant not in AMINO_ACIDS:
        raise DemoError("Pick one of the 20 standard amino acids.")
    land = np.array(landscape(seq, size)["llr"])
    llr = float(land[position, AA_INDEX[mutant]])
    # Percentile among all possible substitutions in this protein (wild-type cells excluded)
    is_wt = np.zeros_like(land, dtype=bool)
    is_wt[np.arange(len(seq)), [AA_INDEX[a] for a in seq]] = True
    pct = float((land[~is_wt] < llr).mean() * 100) if mutant != seq[position] else None
    return dict(position=position, wt=seq[position], mutant=mutant, llr=round(llr, 3), percentile=pct, model=size)


# ---------------------------------------------------------------- explanations

def explain(seq: str, position: int, mutant: str, size: str, key: str | None = None) -> dict:
    """Integrated gradients: how much each residue pushed this variant's score up or down."""
    seq = clean_sequence(seq)
    _check_position(seq, position)
    wt = seq[position]
    if mutant == wt or mutant not in AMINO_ACIDS:
        raise DemoError("Pick an amino acid different from the original one.")
    with _LOCK:
        ig = integrated_gradients(bundle(size), seq, position, wt, mutant, steps=16, max_steps=64)
    attr = ig.attribution
    dist = _distances(key, seq)
    order = [int(j) for j in np.argsort(-np.abs(attr)) if j != position][:10]
    top = []
    for j in order:
        row = dict(position=j, aa=seq[j], value=round(float(attr[j]), 4), seq_separation=abs(j - position))
        if dist is not None:
            row["distance"] = round(float(dist[position, j]), 1)
            row["contact"] = bool(dist[position, j] < 8.0)
        top.append(row)
    delta = ig.f_input - ig.f_baseline
    return dict(
        position=position, wt=wt, mutant=mutant, model=size, llr=round(ig.f_input, 3),
        attribution=np.round(attr, 4).tolist(), top=top, steps=ig.steps,
        completeness_gap=round(ig.completeness_gap, 3),
        converged=abs(ig.completeness_gap) <= 0.1 * max(abs(delta), 0.5),
        has_structure=dist is not None,
    )


@torch.no_grad()
def _head_effects(b: ModelBundle, seq: str, position: int, wt: str, mt: str, chunk: int = 40) -> np.ndarray:
    """Change in this variant's LLR when each head is mean-ablated, all heads batched."""
    means = head_means(b, seq)
    H, dh = b.n_heads, b.d_head
    heads = [(l, h) for l in range(b.n_layers) for h in range(H)]
    base = b.encode(seq)[0]
    base[position + 1] = b.mask_id
    tok = position + 1
    aa = b.aa_ids

    def llr_of(logits):
        lp = logits[:, tok].log_softmax(-1)
        return lp[:, aa[AA_INDEX[mt]]] - lp[:, aa[AA_INDEX[wt]]]

    clean = llr_of(b.model(input_ids=base[None]).logits)[0].item()
    out = np.zeros(len(heads))
    for start in range(0, len(heads), chunk):
        batch = heads[start : start + chunk]
        layer_of = torch.tensor([l for l, _ in batch])
        head_of = torch.tensor([h for _, h in batch])
        handles = []
        for layer in set(layer_of.tolist()):
            rows = (layer_of == layer).nonzero().flatten()
            fill = means[layer][head_of[rows]]  # (n_rows, dh)

            def edit(_m, args, rows=rows, hidx=head_of[rows], fill=fill):
                x = args[0].view(*args[0].shape[:-1], H, dh).clone()
                x[rows, :, hidx, :] = fill[:, None, :]
                return (x.view(*x.shape[:-2], H * dh),) + tuple(args[1:])

            handles.append(b.layers[layer].attention.output.dense.register_forward_pre_hook(edit))
        try:
            vals = llr_of(b.model(input_ids=base.repeat(len(batch), 1)).logits)
        finally:
            for h in handles:
                h.remove()
        out[start : start + len(batch)] = vals.numpy() - clean
    return out.reshape(b.n_layers, H)


@lru_cache(maxsize=4)
def _known_heads(size: str) -> dict:
    """Top contact-detector heads (06_attention_contacts.py) and top variant-score heads (08/11), if computed."""
    out = {"contact": [], "mutation": []}
    f = TABLES_DIR / "attention" / "head_summary.csv"
    if f.exists():
        d = pd.read_csv(f)
        d = d[(d.model == size) & ~d.random_init].nlargest(10, "precision_at_L")
        out["contact"] = [[int(r.layer), int(r["head"])] for _, r in d.iterrows()]
    f = TABLES_DIR / "head_roles" / "head_properties.csv"
    if f.exists():
        d = pd.read_csv(f)
        d = d[d.model == size].nlargest(10, "importance")
        out["mutation"] = [[int(r.layer), int(r["head"])] for _, r in d.iterrows()]
    return out


def inside(seq: str, position: int, mutant: str, size: str) -> dict:
    """Layer-by-layer score (logit lens) and the effect of switching off each head."""
    seq = clean_sequence(seq)
    _check_position(seq, position)
    wt = seq[position]
    if mutant == wt or mutant not in AMINO_ACIDS:
        raise DemoError("Pick an amino acid different from the original one.")
    b = bundle(size)
    with _LOCK:
        lens = layerwise_logprobs(b, seq, [position])[:, 0]
        heads = _head_effects(b, seq, position, wt, mutant)
    layer_llr = (lens[:, AA_INDEX[mutant]] - lens[:, AA_INDEX[wt]]).tolist()
    layer_top = [AMINO_ACIDS[int(i)] for i in lens.argmax(-1)]
    return dict(
        position=position, wt=wt, mutant=mutant, model=size,
        layers=[round(v, 3) for v in layer_llr], layer_top_guess=layer_top,
        heads=np.round(heads, 3).tolist(), n_layers=b.n_layers, n_heads=b.n_heads,
        known_heads=_known_heads(size),
    )
