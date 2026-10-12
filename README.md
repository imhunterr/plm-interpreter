# PLM Interpreter

Mechanistic interpretability of ESM-2 protein language models: **why** does the model think a
mutation is harmful, and can we trust that explanation?

ESM-2 scores a missense variant zero-shot by masking the site and comparing the model's
probability of the mutant and wild-type amino acids (the masked-marginal log-likelihood ratio,
LLR). This project explains that score at three levels and checks every explanation against
experimental data and protein structure:

| Level | Question | Method | Checked against |
|---|---|---|---|
| Residues | Which other residues drove this variant's score? | Integrated gradients, occlusion | AlphaFold2 3D contacts; deletion test |
| Components | Which attention heads / MLPs does variant prediction depend on? | Mean-ablation of every head, layer and MLP | Held-out positions; random-head controls; contact heads |
| Layers | What does each layer know, and when does the variant signal appear? | Linear probes; logit lens | DSSP, solvent accessibility, MSA conservation; DMS data |

Models: ESM-2 8M (6 layers × 20 heads) and 35M (12 × 20), CPU only.
Data: ProteinGym v1.2 (217 deep mutational scanning assays, 2,525 clinical variant sets,
AlphaFold2 structures, MSAs).

## Live demo

A local web app walks through the whole idea on real proteins: pick a protein, hide one residue,
see ESM-2's guess for all 20 amino acids, score a mutation against the lab measurement, then ask
why (residue importance painted on the AlphaFold2 structure) and look inside the model (score after
every layer, effect of switching off each attention head). It also draws the model's and the lab's
mutation maps side by side, and includes the findings and figures.

```bash
conda activate plm-interp
python -m plm_interp.app --warm   # first time: precomputes the showcase examples (cached on disk)
python -m plm_interp.app          # afterwards: starts instantly and opens http://127.0.0.1:8000
```

Everything runs offline on the laptop CPU (the 3D viewer and fonts are bundled). Code:
`src/plm_interp/app/` (`engine.py` reuses the pipeline functions, `server.py` is a small FastAPI app,
`static/` is the plain HTML/JS front end).

## Key findings

1. **The scores are reproduced exactly.** Our masked-marginal scores match ProteinGym's published
   Spearman within 0.01 for all 217 assays and both models (mean 0.319 for 35M, 0.206 for 8M). On 62,727
   curated clinical variants across 2,525 human genes, the 35M score separates pathogenic from benign with
   AUROC 0.744 (8M: 0.649).
2. **Structure is learned, not built in.** Linear probes recover secondary structure (balanced accuracy
   0.80 vs 0.45 from residue identity alone), burial, solvent accessibility, contact number and
   conservation from held-out proteins. The signal grows with depth and is absent in a randomly initialised model.
3. **A handful of late heads are contact detectors.** In ESM-2 35M, head L11H13 ranks true 3D contacts with
   precision@L 0.39, 11x the 3.5% base rate; no random-init head exceeds 0.06.
4. **Single-variant explanations point at 3D neighbours.** Integrated-gradients importance identifies
   residues in 3D contact with the mutated site (AUROC 0.74-0.75; 0.60 for long-range contacts only; 0.46-0.52
   for a random-init model), and masking the top-attributed residue moves the score about 10x more than
   masking a random one. Simple sequence proximity is a strong baseline (AUROC 0.86): much of what the
   model uses is local sequence context.
5. **Variant prediction depends on a few heads, and they are not the contact heads.** Across 24 proteins,
   ablating the 16 heads that matter most on one half of the positions lowers Spearman on the other half by
   0.28 (ESM-2 35M), against 0.05 for 16 random heads; the top heads beat random heads in 23 of 24 proteins
   (sign test p = 1.5e-6). Ablating the 16 best contact-detector heads costs only 0.04, about the same as
   random heads.
6. **The heads that matter look locally.** From the masked site, the 10 most important heads of ESM-2 35M put
   60% of their attention on residues 1-4 positions away (other heads: 23%), look about 5 residues away on
   average (others: 25) and attend to 3D contacts less than other heads. None of them is among the top 10
   contact heads, in either model. This explains finding 4: the variant verdict is built mainly from local
   sequence context, while the fold-sensing machinery sits elsewhere.
7. **The variant signal forms late.** Decoded with the logit lens, the variant score only correlates with
   experiments in the last third of the network (35M: Spearman about 0 at layers 4-8, 0.34 at layer 12).
8. **The model detects damage, not improvement.** Across 209 assays, the 35M score ranks variants well within
   the most damaging quarter (Spearman 0.18) but barely within the best quarter (0.04); this holds in 74% of
   assays and 88% of stability assays (p < 0.001). The score reflects evolutionary plausibility rather than
   function: in HRas, the oncogenic hotspots G12V, G13D and Q61L work normally or better in the lab, yet the
   model rates them unusual. (Measurement noise among near-neutral variants may account for part of the
   asymmetry.)

## Results

<!-- RESULTS -->
### 1. Baseline: the scorer is correct

| model | assays | mean Spearman (ours) | mean Spearman (published) | max per-assay difference | mean AUROC |
|---|---|---|---|---|---|
| ESM-2 35M | 217 | 0.319 | 0.319 | 0.0096 | 0.675 |
| ESM-2 8M | 217 | 0.206 | 0.206 | 0.0016 | 0.612 |

Per-assay agreement within 0.01 Spearman: 217/217 (35M), 217/217 (8M). Proteins longer than ESM-2's 1,022-residue limit are scored as ProteinGym does: on the MSA region when it covers only part of the protein, otherwise in windows centred on each site.

Clinical variants (pathogenic vs benign, score = -LLR):

| model | genes | variants | pooled AUROC | mean per-gene AUROC |
|---|---|---|---|---|
| ESM-2 8M | 2525 | 62727 | 0.649 | 0.641 (n=2525) |
| ESM-2 35M | 2525 | 62727 | 0.744 | 0.726 (n=2525) |

### 2. Probing: what each layer encodes (held-out proteins)

| task | metric | residue identity only | 8M best layer | 35M best layer | 35M random-init best |
|---|---|---|---|---|---|
| buried | auroc | 0.734 | 0.814 (L6) | 0.847 (L12) | 0.734 (L0) |
| conservation | spearman | 0.176 | 0.338 (L3) | 0.521 (L11) | 0.259 (L12) |
| n_contacts | spearman | 0.393 | 0.657 (L6) | 0.720 (L12) | 0.411 (L7) |
| rsa | spearman | 0.450 | 0.627 (L6) | 0.691 (L12) | 0.450 (L0) |
| secondary_structure | balanced_accuracy | 0.448 | 0.737 (L5) | 0.797 (L11) | 0.455 (L9) |

### 3. Attention heads that track 3D contacts

- ESM-2 8M: best heads L5H4 (0.23), L5H10 (0.19), L5H3 (0.12); contact base rate 0.035; best random-init head 0.050.
- ESM-2 35M: best heads L11H13 (0.39), L11H14 (0.31), L10H15 (0.24); contact base rate 0.035; best random-init head 0.055.

### 4. Residue attributions point at 3D neighbours

| model | variants | IG converged | contact AUROC | long-range AUROC | seq-distance AUROC | top-10 in contact | IG vs occlusion (Spearman) | random-init AUROC |
|---|---|---|---|---|---|---|---|---|
| ESM-2 8M | 200 | 96% | 0.748 | 0.596 | 0.859 | 31% | 0.300 | 0.460 |
| ESM-2 35M | 200 | 50% | 0.741 | 0.597 | 0.859 | 31% | 0.259 | 0.517 |

Robustness: restricted to IG runs whose attributions sum to the score change within 10%, the contact AUROC is 0.751 (35M), 0.748 (8M).

Deletion test (ESM-2 35M): mean |change in score| after masking k residues

| k | top attributed | sequence-nearest | random |
|---|---|---|---|
| 1 | 0.370 | 0.275 | 0.038 |
| 2 | 0.520 | 0.456 | 0.062 |
| 4 | 0.661 | 0.655 | 0.098 |
| 8 | 0.947 | 1.076 | 0.148 |
| 16 | 1.365 | 1.654 | 0.228 |

### 5. Causal tests: ablating heads

Drop in held-out Spearman (heads ranked on other positions), mean over showcase proteins:

| model | k | top-k heads | random k heads | top-k contact heads |
|---|---|---|---|---|
| ESM-2 35M | 1 | 0.055 | 0.002 | 0.005 |
| ESM-2 35M | 2 | 0.091 | 0.004 | 0.006 |
| ESM-2 35M | 4 | 0.112 | 0.008 | 0.008 |
| ESM-2 35M | 8 | 0.177 | 0.018 | 0.020 |
| ESM-2 35M | 16 | 0.278 | 0.048 | 0.038 |
| ESM-2 8M | 1 | 0.043 | 0.003 | 0.000 |
| ESM-2 8M | 2 | 0.067 | 0.008 | 0.005 |
| ESM-2 8M | 4 | 0.108 | 0.017 | 0.010 |
| ESM-2 8M | 8 | 0.146 | 0.034 | 0.015 |
| ESM-2 8M | 16 | 0.193 | 0.072 | 0.046 |

Largest drop in Spearman from ablating one whole layer, in any single protein: 35M attn_layer -0.706; 35M mlp -0.672; 8M attn_layer -0.548; 8M mlp -0.795

### 6. Logit lens: when the variant signal appears

Mean Spearman with DMS when decoding from each layer:

- ESM-2 35M: L0 0.08 → L1 0.12 → L2 0.13 → L3 0.03 → L4 0.04 → L5 -0.02 → L6 -0.05 → L7 -0.02 → L8 0.05 → L9 0.11 → L10 0.17 → L11 0.26 → L12 0.34
- ESM-2 8M: L0 0.02 → L1 -0.16 → L2 0.01 → L3 0.03 → L4 0.01 → L5 0.10 → L6 0.23

### 7. Damage, not improvement: what the score measures

Within each assay, Spearman between model score and lab result, computed separately inside the worst and the best quarter of variants (209 assays with >= 200 single mutants):

| model | assays | within worst 25% | within best 25% | difference [95% CI] | assays where worst > best | top-5% dip |
|---|---|---|---|---|---|---|
| ESM-2 8M | 209 | 0.119 | 0.032 | 0.088 [0.062, 0.114] | 65% | +1.17 |
| ESM-2 35M | 209 | 0.179 | 0.042 | 0.138 [0.109, 0.169] | 74% | +1.45 |
| ESM-2 35M (stability) | 66 | 0.287 | 0.049 | 0.237 [0.190, 0.283] | 88% | +1.52 |
| ESM-2 8M (stability) | 66 | 0.178 | 0.036 | 0.142 [0.093, 0.198] | 76% | +1.08 |

The model ranks damaging variants much better than beneficial ones. On average it does not rate the best variants as worse than good ones (the top-5% dip is positive), but individual assays such as HRas do show that dip. Part of the asymmetry may come from measurement noise among near-neutral variants, which dominate the best quarter.

### 8. What the important heads look at

Attention from the masked site, averaged over sampled sites in the ablation proteins. Mutation heads = the 10 heads whose ablation hurts variant scores most; contact heads = the 10 best contact detectors from section 3.

| model | group | mean layer | to neighbours (1-4 away) | to 3D contacts | to conserved positions | mean distance | entropy |
|---|---|---|---|---|---|---|---|
| ESM-2 35M | contact heads | 10.2 | 0.17 | 0.087 | 0.28 | 27.7 | 0.78 |
| ESM-2 35M | mutation heads | 9.0 | 0.60 | 0.018 | 0.29 | 5.1 | 0.35 |
| ESM-2 35M | other | 5.1 | 0.23 | 0.033 | 0.26 | 25.2 | 0.70 |
| ESM-2 8M | contact heads | 4.4 | 0.21 | 0.046 | 0.27 | 23.5 | 0.81 |
| ESM-2 8M | mutation heads | 3.2 | 0.55 | 0.017 | 0.27 | 6.6 | 0.41 |
| ESM-2 8M | other | 2.2 | 0.25 | 0.026 | 0.24 | 25.2 | 0.70 |

### 9. Case studies

Open `outputs/explanations/35M/index.html` for the interactive cards.

| protein | variant | category | model LLR (pct) | assay score (pct) | top-10 residues in 3D contact |
|---|---|---|---|---|---|
| RASH_HUMAN | G10I | agree damaging | -13.02 (0) | -0.899 (0) | 50% |
| RASH_HUMAN | K117G | model too harsh | -9.52 (5) | +0.424 (99) | 50% |
| RASH_HUMAN | H27R | model missed | +2.78 (100) | -0.439 (19) | 10% |
| RASH_HUMAN | G12V | hotspot G12V | -4.22 (41) | +0.022 (66) | 40% |
| RASH_HUMAN | G13D | hotspot G13D | -4.81 (35) | +0.254 (93) | 40% |
| RASH_HUMAN | Q61L | hotspot Q61L | -3.14 (54) | +0.262 (94) | 20% |
| TPK1_HUMAN | D73R | agree damaging | -11.16 (0) | -0.551 (0) | 10% |
| TPK1_HUMAN | I23R | model too harsh | -11.62 (0) | +1.410 (95) | 30% |
| TPK1_HUMAN | W237I | model missed | +3.15 (100) | -0.353 (2) | 30% |
| OTC_HUMAN | C303W | agree damaging | -10.19 (0) | -0.037 (0) | 30% |
| OTC_HUMAN | S267G | model too harsh | -5.42 (11) | +1.153 (98) | 60% |
| OTC_HUMAN | W298L | model missed | +3.70 (100) | -0.018 (2) | 20% |
| HEM3_HUMAN | L42D | agree damaging | -15.94 (0) | +0.000 (0) | 50% |
| HEM3_HUMAN | R22W | model too harsh | -10.13 (6) | +1.210 (94) | 60% |
| HEM3_HUMAN | G203T | model missed | +2.73 (100) | +0.000 (0) | 20% |
| P53_HUMAN | D259Y | agree damaging | -6.01 (0) | -0.503 (3) | 40% |
| P53_HUMAN | P153H | model too harsh | -2.40 (3) | +3.105 (97) | 20% |
| P53_HUMAN | E271P | model missed | +1.51 (99) | -0.569 (2) | 10% |
| P53_HUMAN | R175H | hotspot R175H | -0.52 (33) | -0.240 (13) | 30% |
| P53_HUMAN | R248W | hotspot R248W | -1.90 (6) | -0.051 (22) | 10% |
| P53_HUMAN | R273H | hotspot R273H | -1.14 (15) | -0.393 (6) | 30% |

Every claim above, with 95% bootstrap CIs, controls and FDR-corrected permutation tests: [`outputs/runs/verification/claims.md`](outputs/runs/verification/claims.md).

## Pipeline

Each script reads the previous stages' outputs, can be re-run safely, and writes to
`outputs/runs/<stage>/` (tables) or `outputs/figures/` (figures). Settings live in
`configs/default.yaml`.

| Script | Stage | Output |
|---|---|---|
| `01_download_models.py` | Download ESM-2 8M/35M to `models/` | |
| `02_download_data.py` | Download ProteinGym v1.2 to `data/raw/` (md5-verified) | `data/raw/MANIFEST.json` |
| `03_baseline.py` | Zero-shot scores for every assay and clinical set; compare with ProteinGym's published numbers | `runs/baseline/` |
| `04_annotate.py` | Per-residue labels from AF2 structures and MSAs for 60 proteins | `data/processed/annotations/` |
| `05_probing.py` | Linear probes per layer (held-out proteins) | `runs/probing/` |
| `06_attention_contacts.py` | Attention heads vs 3D contact maps | `runs/attention/` |
| `07_attribution.py` | Integrated gradients + occlusion for sampled variants; contact AUROC; deletion test | `runs/attribution/` |
| `08_head_ablation.py` | Ablate every head / layer / MLP; top-k vs random-k on held-out positions | `runs/ablation/` |
| `09_logit_lens.py` | Variant signal decoded from every layer | `runs/logit_lens/` |
| `10_verification.py` | Effect sizes, bootstrap CIs, permutation tests, FDR for every claim | `runs/verification/claims.md` |
| `11_figures.py` | Paper figures | `outputs/figures/` |
| `12_case_studies.py` | Per-variant explanation cards with 3D views | `outputs/explanations/` |
| `14_evolution_vs_function.py` | Does the score rank damage and improvement equally well? (all assays) | `runs/evolution_vs_function/` |
| `15_head_roles.py` | What the heads that matter for variant scores attend to | `runs/head_roles/` |
| `13_report.py` | Collect headline numbers into this README and `runs/report.md` | `runs/report.md` |

```bash
conda activate plm-interp
python scripts/01_download_models.py
python scripts/02_download_data.py
for s in 03_baseline 04_annotate 05_probing 06_attention_contacts 07_attribution \
         08_head_ablation 09_logit_lens 14_evolution_vs_function 15_head_roles \
         10_verification 11_figures 12_case_studies 13_report; do
  python scripts/$s.py
done
pytest
```

The full baseline (both models, all assays and clinical sets) takes several hours on an
8-core CPU; every other stage takes minutes to about an hour.

## Design decisions

- **What is explained.** Always the masked-marginal LLR at the mutated site, the same quantity
  ProteinGym benchmarks. Because the site is masked, wild-type and mutant inputs are identical,
  so explanations perturb the *context* (other residues, heads, layers), never the site itself.
- **Correctness gate.** Our scores reproduce ProteinGym's published per-assay Spearman for both
  models before any interpretation is run (see `tests/test_scoring.py`).
- **Integrated gradients baseline.** ESM-2 embeds `<mask>` as a zero vector, but its pre-LayerNorm
  blocks are scale invariant, so a zero baseline makes the path integral degenerate. We use the
  mean amino-acid embedding as the baseline, adapt the number of steps until the completeness
  axiom holds within 10%, and report the remaining gap for every variant. The model output can
  jump abruptly along the path, which is why a fixed small step count is not enough.
- **Controls everywhere.** Every analysis is repeated on a randomly initialised model with the
  same architecture; probes are compared with residue identity alone; heads are ranked on one half
  of the positions and evaluated on the other; attributions are compared with random and
  sequence-nearest residues.
- **Structural labels only where the structure is trustworthy.** Residues with AF2 pLDDT < 70
  are excluded from structure-based labels.
- **Plain PyTorch hooks** on HuggingFace `EsmForMaskedLM` (`src/plm_interp/models/hooks.py`);
  TransformerLens does not support the ESM architecture.

## Layout

```
configs/default.yaml          shared settings
src/plm_interp/
  models/    loader.py (ESM-2 + helpers), hooks.py (residual/head capture, head ablation)
  baseline/  scoring.py (masked-marginal scores, long-sequence windows)
  data/      proteingym.py, structure.py (DSSP, RSA, contacts), msa.py (weighted conservation)
  attribution/ methods.py (IG, occlusion, deletion), attention.py (heads vs contacts)
  explain/   probing.py, ablation.py, logit_lens.py
  verification/ stats.py (bootstrap, permutation tests, BH-FDR)
  viz/       style.py
scripts/                      numbered pipeline stages
tests/                        pytest suite (uses the 8M model)
outputs/                      figures, result tables, explanation cards
```

## Environment

```bash
conda env create -f environment.yml
conda activate plm-interp
pip install -r requirements.lock.txt --extra-index-url https://download.pytorch.org/whl/cpu
pip install -e .
```

## References

- Lin et al. 2023, *Evolutionary-scale prediction of atomic-level protein structure with a language model*, Science.
- Meier et al. 2021, *Language models enable zero-shot prediction of the effects of mutations on protein function*, NeurIPS.
- Notin et al. 2023, *ProteinGym: large-scale benchmarks for protein fitness prediction and design*, NeurIPS.
- Vig et al. 2021, *BERTology meets biology: interpreting attention in protein language models*, ICLR.
- Rao et al. 2021, *Transformer protein language models are unsupervised structure learners*, ICLR.
- Zhang et al. 2024, *Protein language models learn evolutionary statistics of interacting sequence motifs*, PNAS.
- Simon & Zou 2025, *InterPLM: discovering interpretable features in protein language models via sparse autoencoders*, Nature Methods.
- Sundararajan et al. 2017, *Axiomatic attribution for deep networks*, ICML.
- Adebayo et al. 2018, *Sanity checks for saliency maps*, NeurIPS.
- nostalgebraist 2020, *interpreting GPT: the logit lens*; Belrose et al. 2023, *Eliciting latent predictions from transformers with the tuned lens*.
