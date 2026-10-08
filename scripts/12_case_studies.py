"""Per-variant explanation cards: residue, head and layer-level explanations of single predictions.

For each showcase protein, picks:
  agree_damaging  - both the model and the assay rate it among the most damaging
  model_too_harsh - the model calls it very unlikely, the assay says it works fine
  model_missed    - the assay says damaging, the model barely objects
plus named hotspots from the config (e.g. HRas G12V, Q61L).

Each card combines:
  residues: integrated-gradients and occlusion importance, with 3D distance to the site
  heads:    change in this variant's score when each attention head is mean-ablated
  layers:   logit-lens score after each layer
  context:  what the model thinks fits at this position vs what the assay measured

Outputs (outputs/explanations/<model>/<dms_id>/): <mutant>.json, <mutant>.html (3D view + card),
and outputs/explanations/<model>/index.html.

Usage:
    python scripts/12_case_studies.py [--model 35M]
"""

from __future__ import annotations

import argparse
import html
import json

import numpy as np
import pandas as pd
import py3Dmol
import torch

from plm_interp.attribution.methods import integrated_gradients, occlusion_attribution, occlusion_logprobs
from plm_interp.baseline.scoring import AA_INDEX, llr_matrix, masked_logprobs, score_mutants
from plm_interp.data.proteingym import dms_info, load_dms, parse_mutant
from plm_interp.explain.ablation import head_means
from plm_interp.explain.logit_lens import layerwise_logprobs
from plm_interp.models.hooks import ablate_heads
from plm_interp.models.loader import AMINO_ACIDS, load_model
from plm_interp.paths import ANNOTATIONS_DIR, EXPLANATIONS_DIR, STRUCTURES_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("case_studies", "case_studies.log")

HOTSPOTS = {"RASH_HUMAN_Bandaru_2017": ["G12V", "G13D", "Q61L"], "P53_HUMAN_Kotler_2018": ["R175H", "R248W", "R273H"]}
SS_NAMES = {"H": "helix", "E": "strand", "-": "loop"}


def pick_variants(df: pd.DataFrame) -> dict[str, str]:
    """Rank-based selection; df has DMS_score and llr (singles only)."""
    r_dms = df.DMS_score.rank(pct=True)
    r_llr = df.llr.rank(pct=True)
    picks = {}
    picks["agree_damaging"] = df.mutant[(r_dms + r_llr).idxmin()]
    picks["model_too_harsh"] = df.mutant[((1 - r_llr) + r_dms).idxmax()]
    picks["model_missed"] = df.mutant[(r_llr + (1 - r_dms)).idxmax()]
    return picks


@torch.no_grad()
def head_effects(bundle, seq, pos, wt, mt, means) -> np.ndarray:
    """(layers, heads) change in this variant's LLR when each head is mean-ablated."""
    base = bundle.encode(seq)
    base[0, pos + 1] = bundle.mask_id
    tok = pos + 1
    aa = bundle.aa_ids

    def llr():
        lp = bundle.model(input_ids=base).logits[0, tok].log_softmax(-1)
        return (lp[aa[AA_INDEX[mt]]] - lp[aa[AA_INDEX[wt]]]).item()

    clean = llr()
    out = np.zeros((bundle.n_layers, bundle.n_heads))
    for layer in range(bundle.n_layers):
        for head in range(bundle.n_heads):
            with ablate_heads(bundle, {layer: [head]}, means):
                out[layer, head] = llr() - clean
    return out


def structure_html(pdb_file: str, importance: np.ndarray, pos: int, top: list[int]) -> str:
    """3Dmol view: cartoon coloured by |importance| (white -> blue), mutated site in orange sticks."""
    lines = (STRUCTURES_DIR / pdb_file).read_text().splitlines()
    imp = np.abs(importance)
    imp = imp / imp.max() if imp.max() > 0 else imp
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            resi = int(line[22:26]) - 1
            b = imp[resi] if 0 <= resi < len(imp) else 0.0
            line = f"{line[:60]}{b * 100:6.2f}{line[66:]}"
        out.append(line)
    view = py3Dmol.view(width=560, height=440)
    view.addModel("\n".join(out), "pdb")
    view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "linear", "colors": ["#f0efec", "#86b6ef", "#184f95"], "min": 0, "max": 100}}})
    view.addStyle({"resi": [j + 1 for j in top]}, {"stick": {"colorscheme": "grayCarbon", "radius": 0.18}})
    view.addStyle({"resi": [pos + 1]}, {"stick": {"color": "#eb6834", "radius": 0.3}, "sphere": {"color": "#eb6834", "radius": 0.9, "opacity": 0.5}})
    view.zoomTo({"resi": [pos + 1] + [j + 1 for j in top[:5]]})
    return view._make_html()


def card_html(card: dict, viewer: str) -> str:
    rows = "".join(
        f"<tr><td>{r['residue']}</td><td>{r['ig']:+.3f}</td><td>{r['occlusion']:+.3f}</td><td>{r['distance_A']:.1f}</td>"
        f"<td>{r['seq_separation']}</td><td>{'yes' if r['contact'] else ''}</td></tr>"
        for r in card["top_residues"]
    )
    heads = ", ".join(f"L{h['layer']}H{h['head']} ({h['delta_llr']:+.2f})" for h in card["top_heads"])
    layers = " → ".join(f"{v:+.1f}" for v in card["logit_lens_llr"])
    fits = "".join(
        f"<tr><td>{a}</td><td>{p:.3f}</td><td>{'' if d is None else f'{d:+.2f}'}</td></tr>"
        for a, p, d in card["position_table"]
    )
    s = card["site"]
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(card['protein'])} {card['mutant']}</title>
<style>
body{{font-family:system-ui,sans-serif;margin:24px;color:#1b2230;background:#fcfcfb;max-width:1200px}}
.grid{{display:grid;grid-template-columns:minmax(0,560px) minmax(0,1fr);gap:24px}}
@media(max-width:900px){{.grid{{grid-template-columns:1fr}}}}
table{{border-collapse:collapse;font-size:13px;font-variant-numeric:tabular-nums}}td,th{{padding:3px 10px;border-bottom:1px solid #e3e2dd;text-align:left}}
.pill{{display:inline-block;padding:2px 10px;border-radius:99px;background:#e6edfb;color:#184f95;font-size:12px}}
.muted{{color:#5b6577;font-size:13px}} h1{{margin:0 0 4px}} h2{{font-size:15px;margin:18px 0 6px}}
</style></head><body>
<h1>{html.escape(card['protein'])} · {card['mutant']}</h1>
<p><span class="pill">{card['category'].replace('_', ' ')}</span> &nbsp; model score (LLR) <b>{card['llr']:+.2f}</b>
 (percentile {card['llr_percentile']:.0f}) · assay score <b>{card['dms_score']:+.3f}</b> (percentile {card['dms_percentile']:.0f}) · model {card['model']}</p>
<p class="muted">Site: {SS_NAMES.get(s['ss'], s['ss'])}, relative solvent accessibility {s['rsa']:.2f}, conservation {s['conservation']:.2f}, AF2 pLDDT {s['plddt']:.0f}.
IG completeness gap {card['ig_completeness_gap']:+.3f} ({card['ig_steps']} steps).</p>
<div class="grid"><div>{viewer}<p class="muted">Colour = |IG importance| (light → dark blue); orange = mutated site; grey sticks = top-10 residues.</p></div>
<div><h2>Which residues drove the score</h2><table><tr><th>residue</th><th>IG</th><th>occlusion</th><th>3D dist (Å)</th><th>seq sep</th><th>contact</th></tr>{rows}</table>
<p class="muted">Positive = pushed the score up (made the mutation look more acceptable); negative = made it look worse.</p>
<h2>Which heads matter for this variant</h2><p>{heads}</p>
<h2>When the model decides (score after each layer)</h2><p>{layers}</p>
<h2>What fits at position {s['position']}</h2><table><tr><th>aa</th><th>model prob</th><th>assay score</th></tr>{fits}</table></div></div>
</body></html>"""


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=cfg["main_model"])
    ap.add_argument("--threads", type=int, default=cfg["threads"])
    args = ap.parse_args()
    set_threads(args.threads)
    s_cfg = cfg["structure"]
    bundle = load_model(args.model)
    residues = pd.read_parquet(ANNOTATIONS_DIR / "residues.parquet")
    root = ensure_dir(EXPLANATIONS_DIR / args.model)
    index_rows = []

    for dms_id in cfg["showcase"]:
        info = dms_info(dms_id)
        seq = info.target_seq
        res = residues[residues.uniprot == info.UniProt_ID].reset_index(drop=True)
        dist = np.load(ANNOTATIONS_DIR / "dist" / f"{info.UniProt_ID}.npy").astype(np.float32)
        dms = load_dms(dms_id)
        dms = dms[~dms.mutant.str.contains(":")].reset_index(drop=True)
        lp = masked_logprobs(bundle, seq)
        dms["llr"] = score_mutants(dms.mutant, lp, seq)
        llr_tab = llr_matrix(lp, seq)
        probs = np.exp(lp)
        means = head_means(bundle, seq)
        out_dir = ensure_dir(root / dms_id)

        picks = pick_variants(dms)
        for hs in HOTSPOTS.get(dms_id, []):
            if hs in set(dms.mutant):
                picks[f"hotspot_{hs}"] = hs

        for category, mutant in picks.items():
            wt, pos, mt = parse_mutant(mutant)[0]
            row = dms[dms.mutant == mutant].iloc[0]
            ig = integrated_gradients(bundle, seq, pos, wt, mt, cfg["attribution"]["ig_steps"],
                                      cfg["attribution"]["ig_max_steps"], cfg["attribution"]["ig_rel_tol"])
            occ = occlusion_attribution(occlusion_logprobs(bundle, seq, pos), pos, wt, mt)
            heads = head_effects(bundle, seq, pos, wt, mt, means)
            lens = layerwise_logprobs(bundle, seq, [pos])[:, 0]
            lens_llr = (lens[:, AA_INDEX[mt]] - lens[:, AA_INDEX[wt]]).tolist()

            order = [j for j in np.argsort(-np.abs(ig.attribution)) if j != pos][:10]
            top_res = [
                dict(residue=f"{seq[j]}{j + 1}", ig=float(ig.attribution[j]), occlusion=float(occ[j]),
                     distance_A=float(dist[pos, j]), seq_separation=int(abs(j - pos)),
                     contact=bool(dist[pos, j] < s_cfg["contact_cutoff"]))
                for j in order
            ]
            flat = [(l, h, heads[l, h]) for l in range(heads.shape[0]) for h in range(heads.shape[1])]
            top_heads = [dict(layer=l, head=h, delta_llr=float(d)) for l, h, d in sorted(flat, key=lambda x: -abs(x[2]))[:6]]
            pos_dms = dms[dms.mutant.str.match(rf"^{wt}{pos + 1}[A-Z]$")].set_index(dms.mutant.str[-1])
            table = [
                (a, float(probs.loc[pos, a]), float(pos_dms.DMS_score[a]) if a in pos_dms.index else None)
                for a in sorted(AMINO_ACIDS, key=lambda a: -probs.loc[pos, a])
            ]
            card = dict(
                protein=info.UniProt_ID, dms_id=dms_id, model=args.model, mutant=mutant, category=category,
                llr=float(row.llr), dms_score=float(row.DMS_score),
                llr_percentile=float((dms.llr < row.llr).mean() * 100),
                dms_percentile=float((dms.DMS_score < row.DMS_score).mean() * 100),
                site=dict(position=pos + 1, wt=wt, ss=str(res.ss[pos]), rsa=float(res.rsa[pos]),
                          conservation=float(res.conservation[pos]) if not pd.isna(res.conservation[pos]) else float("nan"),
                          plddt=float(res.plddt[pos])),
                ig_completeness_gap=ig.completeness_gap, ig_steps=ig.steps,
                top_residues=top_res, top10_contact_fraction=float(np.mean([r["contact"] for r in top_res])),
                top_heads=top_heads, logit_lens_llr=lens_llr, position_table=table,
                most_likely_at_site=table[0][0], model_rank_of_wt=int([t[0] for t in table].index(wt)) + 1,
                llr_row=llr_tab.loc[pos].round(3).to_dict(),
            )
            (out_dir / f"{mutant}.json").write_text(json.dumps(card, indent=2))
            viewer = structure_html(info.pdb_file, ig.attribution, pos, order)
            (out_dir / f"{mutant}.html").write_text(card_html(card, viewer))
            index_rows.append(card)
            log.info(f"{dms_id} {category} {mutant}: llr {row.llr:+.2f}, dms {row.DMS_score:+.3f}, "
                     f"top-10 contacts {card['top10_contact_fraction']:.0%}")

    items = "".join(
        f"<tr><td>{c['protein']}</td><td><a href='{c['dms_id']}/{c['mutant']}.html'>{c['mutant']}</a></td>"
        f"<td>{c['category'].replace('_', ' ')}</td><td>{c['llr']:+.2f}</td><td>{c['dms_score']:+.3f}</td>"
        f"<td>{c['top10_contact_fraction']:.0%}</td></tr>"
        for c in index_rows
    )
    (root / "index.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>PLM Interpreter case studies</title>"
        "<style>body{font-family:system-ui,sans-serif;margin:24px}td,th{padding:4px 12px;border-bottom:1px solid #ddd;text-align:left}</style>"
        f"<h1>Case studies (ESM-2 {args.model})</h1><table><tr><th>protein</th><th>variant</th><th>category</th>"
        f"<th>model LLR</th><th>assay score</th><th>top-10 residues in 3D contact</th></tr>{items}</table>"
    )
    pd.DataFrame([{k: v for k, v in c.items() if k not in ("top_residues", "top_heads", "position_table", "llr_row", "site", "logit_lens_llr")}
                  for c in index_rows]).to_csv(root / "cases.csv", index=False)


if __name__ == "__main__":
    main()
