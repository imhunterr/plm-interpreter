# PLM Interpreter results

### 1. Baseline: the scorer is correct

| model | assays | mean Spearman (ours) | mean Spearman (published) | max per-assay difference | mean AUROC |
|---|---|---|---|---|---|
| ESM-2 35M | 217 | 0.319 | 0.319 | 0.1112 | 0.675 |
| ESM-2 8M | 217 | 0.203 | 0.206 | 0.3963 | 0.611 |

Per-assay agreement within 0.01 Spearman: 215/217 (35M), 215/217 (8M). The exceptions are KCNH2_HUMAN_Kozek_2020, SCN5A_HUMAN_Glazer_2019, proteins longer than ESM-2's 1,022-residue limit, which are scored in windows; our window placement differs from ProteinGym's for these.

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
| ESM-2 35M | 1 | 0.007 | 0.002 | -0.004 |
| ESM-2 35M | 2 | 0.045 | 0.004 | -0.005 |
| ESM-2 35M | 4 | 0.067 | 0.010 | 0.000 |
| ESM-2 35M | 8 | 0.114 | 0.021 | 0.004 |
| ESM-2 35M | 16 | 0.164 | 0.033 | 0.027 |
| ESM-2 8M | 1 | 0.026 | 0.000 | -0.004 |
| ESM-2 8M | 2 | 0.027 | -0.000 | 0.004 |
| ESM-2 8M | 4 | 0.069 | 0.004 | 0.021 |
| ESM-2 8M | 8 | 0.079 | 0.004 | 0.037 |
| ESM-2 8M | 16 | 0.102 | 0.031 | 0.015 |

Largest drop in Spearman from ablating one whole layer, in any single protein: 35M attn_layer -0.487; 35M mlp -0.351; 8M attn_layer -0.472; 8M mlp -0.597

### 6. Logit lens: when the variant signal appears

Mean Spearman with DMS when decoding from each layer:

- ESM-2 35M: L0 0.08 → L1 0.12 → L2 0.13 → L3 0.03 → L4 0.04 → L5 -0.02 → L6 -0.05 → L7 -0.02 → L8 0.05 → L9 0.11 → L10 0.17 → L11 0.26 → L12 0.34
- ESM-2 8M: L0 0.02 → L1 -0.16 → L2 0.01 → L3 0.03 → L4 0.01 → L5 0.10 → L6 0.23

### 7. Case studies

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
