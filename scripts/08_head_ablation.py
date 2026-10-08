"""Stage 6: which attention heads and MLPs does variant-effect prediction causally depend on?

For each showcase protein, a fixed sample of mutated positions is split into a discovery half
and an evaluation half. Every head (and every whole attention layer / MLP) is mean-ablated in
turn and the masked-marginal scores are recomputed. Heads are ranked on the discovery half; the
top-k heads are then ablated together and evaluated on the held-out half against random sets of
k heads and against the top "contact heads" from Stage 4.

Outputs (outputs/runs/ablation/):
    components.csv   one row per (model, protein, component): delta Spearman on discovery/eval/all
    topk.csv         top-k vs random-k vs contact-k ablations on the evaluation half

Usage:
    python scripts/08_head_ablation.py [--models 35M] [--proteins RASH_HUMAN_Bandaru_2017]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from plm_interp.baseline.scoring import masked_logprobs, score_mutants
from plm_interp.data.proteingym import dms_info, load_dms, parse_mutant
from plm_interp.explain.ablation import ablate_mlps, head_means, mlp_means
from plm_interp.models.hooks import ablate_heads
from plm_interp.models.loader import load_model
from plm_interp.paths import RUNS_DIR
from plm_interp.utils.common import ensure_dir, get_logger, load_config, set_threads

log = get_logger("ablation", "ablation.log")

N_POSITIONS = 48
TOP_K = [1, 2, 4, 8, 16]


def main() -> None:
    cfg = load_config()
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=cfg["models"])
    ap.add_argument("--proteins", nargs="+", default=cfg["showcase"])
    ap.add_argument("--threads", type=int, default=cfg["threads"])
    args = ap.parse_args()
    set_threads(args.threads)
    out = ensure_dir(RUNS_DIR / "ablation")
    n_random = cfg["ablation"]["n_random_sets"]
    use_mean = cfg["ablation"]["mode"] == "mean"

    contact_heads = None
    head_summary = RUNS_DIR / "attention" / "head_summary.csv"
    if head_summary.exists():
        contact_heads = pd.read_csv(head_summary)
        contact_heads = contact_heads[~contact_heads.random_init]

    comp_rows, topk_rows = [], []
    for size in args.models:
        bundle = load_model(size)
        all_heads = [(l, h) for l in range(bundle.n_layers) for h in range(bundle.n_heads)]
        for dms_id in args.proteins:
            t0 = time.time()
            seq = dms_info(dms_id).target_seq
            dms = load_dms(dms_id)
            dms = dms[~dms.mutant.str.contains(":")].copy()
            dms["pos"] = [parse_mutant(m)[0][1] for m in dms.mutant]
            rng = np.random.default_rng(cfg["seed"])
            positions = sorted(rng.choice(dms.pos.unique(), size=min(N_POSITIONS, dms.pos.nunique()), replace=False))
            disc_pos = set(positions[0::2])
            dms = dms[dms.pos.isin(positions)].reset_index(drop=True)
            is_disc = dms.pos.isin(disc_pos).to_numpy()
            y = dms.DMS_score.to_numpy()

            def evaluate(llr: np.ndarray) -> dict:
                return dict(
                    rho_disc=spearmanr(llr[is_disc], y[is_disc])[0],
                    rho_eval=spearmanr(llr[~is_disc], y[~is_disc])[0],
                    rho_all=spearmanr(llr, y)[0],
                )

            def run(ctx=None) -> np.ndarray:
                if ctx is None:
                    lp = masked_logprobs(bundle, seq, positions)
                else:
                    with ctx:
                        lp = masked_logprobs(bundle, seq, positions)
                return score_mutants(dms.mutant, lp, seq)

            clean_llr = run()
            clean = evaluate(clean_llr)
            hm = head_means(bundle, seq) if use_mean else None
            mm = mlp_means(bundle, seq) if use_mean else None

            def record(kind, layer, head, llr):
                r = evaluate(llr)
                comp_rows.append(
                    dict(
                        model=size, dms_id=dms_id, kind=kind, layer=layer, head=head,
                        **{f"d_{k}": r[k] - clean[k] for k in r},
                        mean_abs_dllr=float(np.mean(np.abs(llr - clean_llr))),
                        **{f"clean_{k}": v for k, v in clean.items()},
                    )
                )

            for layer, head in all_heads:
                record("head", layer, head, run(ablate_heads(bundle, {layer: [head]}, hm)))
            for layer in range(bundle.n_layers):
                record("attn_layer", layer, -1, run(ablate_heads(bundle, {layer: list(range(bundle.n_heads))}, hm)))
                record("mlp", layer, -1, run(ablate_mlps(bundle, [layer], mm)))

            heads_df = pd.DataFrame([r for r in comp_rows if r["model"] == size and r["dms_id"] == dms_id and r["kind"] == "head"])
            ranked = heads_df.sort_values("d_rho_disc")[["layer", "head"]].to_numpy().tolist()

            def ablate_set(heads) -> float:
                spec: dict[int, list[int]] = {}
                for l, h in heads:
                    spec.setdefault(int(l), []).append(int(h))
                return evaluate(run(ablate_heads(bundle, spec, hm)))["rho_eval"]

            for k in TOP_K:
                row = dict(model=size, dms_id=dms_id, k=k, clean_rho_eval=clean["rho_eval"])
                row["top_rho_eval"] = ablate_set(ranked[:k])
                rand = [ablate_set([all_heads[i] for i in rng.choice(len(all_heads), k, replace=False)]) for _ in range(n_random)]
                row["random_rho_eval_mean"] = float(np.mean(rand))
                row["random_rho_eval_std"] = float(np.std(rand))
                row["random_rho_eval_min"] = float(np.min(rand))
                if contact_heads is not None:
                    ch = contact_heads[contact_heads.model == size].nlargest(k, "precision_at_L")
                    row["contact_rho_eval"] = ablate_set(ch[["layer", "head"]].to_numpy().tolist())
                topk_rows.append(row)

            pd.DataFrame(comp_rows).to_csv(out / "components.csv", index=False)
            pd.DataFrame(topk_rows).to_csv(out / "topk.csv", index=False)
            log.info(
                f"{size} {dms_id}: clean rho {clean['rho_all']:.3f}; strongest head "
                f"L{ranked[0][0]}H{ranked[0][1]}; {time.time() - t0:.0f}s"
            )


if __name__ == "__main__":
    main()
