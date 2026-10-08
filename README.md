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

## Results

<!-- RESULTS -->

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
| `13_report.py` | Collect headline numbers into this README and `runs/report.md` | `runs/report.md` |

```bash
conda activate plm-interp
python scripts/01_download_models.py
python scripts/02_download_data.py
for s in 03_baseline 04_annotate 05_probing 06_attention_contacts 07_attribution \
         08_head_ablation 09_logit_lens 10_verification 11_figures 12_case_studies 13_report; do
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
