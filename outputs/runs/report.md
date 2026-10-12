# PLM Interpreter results

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
