# Verified claims

| claim | model | statement | estimate [95% CI] | control | q | n |
|---|---|---|---|---|---|---|
| C1 | 35M | mean |our Spearman - published Spearman| per assay | 0.001 [0.000, 0.002] |  |  | 217 assays |
| C1 | 35M | mean Spearman with DMS (zero-shot) | 0.319 [0.289, 0.347] | 0.000 | 0.00023 | 217 assays |
| C1 | 8M | mean |our Spearman - published Spearman| per assay | 0.003 [0.000, 0.007] |  |  | 217 assays |
| C1 | 8M | mean Spearman with DMS (zero-shot) | 0.203 [0.174, 0.231] | 0.000 | 0.00023 | 217 assays |
| C2 | 8M | probe buried: best layer 6 (auroc) vs aa identity / random-init best | 0.814  | 0.734 |  | 2592 held-out residues |
| C2 | 8M | probe conservation: best layer 3 (spearman) vs aa identity / random-init best | 0.338  | 0.266 |  | 2733 held-out residues |
| C2 | 8M | probe n_contacts: best layer 6 (spearman) vs aa identity / random-init best | 0.657  | 0.410 |  | 2592 held-out residues |
| C2 | 8M | probe rsa: best layer 6 (spearman) vs aa identity / random-init best | 0.627  | 0.450 |  | 2592 held-out residues |
| C2 | 8M | probe secondary_structure: best layer 5 (balanced_accuracy) vs aa identity / random-init best | 0.737  | 0.451 |  | 2592 held-out residues |
| C2 | 35M | probe buried: best layer 12 (auroc) vs aa identity / random-init best | 0.847  | 0.734 |  | 2592 held-out residues |
| C2 | 35M | probe conservation: best layer 11 (spearman) vs aa identity / random-init best | 0.521  | 0.259 |  | 2733 held-out residues |
| C2 | 35M | probe n_contacts: best layer 12 (spearman) vs aa identity / random-init best | 0.720  | 0.411 |  | 2592 held-out residues |
| C2 | 35M | probe rsa: best layer 12 (spearman) vs aa identity / random-init best | 0.691  | 0.450 |  | 2592 held-out residues |
| C2 | 35M | probe secondary_structure: best layer 11 (balanced_accuracy) vs aa identity / random-init best | 0.797  | 0.455 |  | 2592 held-out residues |
| C3 | 8M | precision@L of best contact head L5H4 vs best random-init head | 0.232 [0.204, 0.259] | 0.050 | 0.00023 | 58 proteins |
| C3 | 8M | same head: precision@L vs contact base rate | 0.232 [0.204, 0.259] | 0.035 | 0.00023 | 58 proteins |
| C3 | 35M | precision@L of best contact head L11H13 vs best random-init head | 0.389 [0.349, 0.426] | 0.055 | 0.00023 | 58 proteins |
| C3 | 35M | same head: precision@L vs contact base rate | 0.389 [0.349, 0.426] | 0.035 | 0.00023 | 58 proteins |
| C4 | 8M | IG |attribution| AUROC for 3D contacts of the mutated site | 0.748 [0.730, 0.764] | 0.500 | 0.00023 | 200 variants |
| C4 | 8M | IG AUROC for long-range contacts (|i-j|>=6) | 0.596 [0.568, 0.621] | 0.500 | 0.00023 | 182 variants |
| C4 | 8M | robustness: IG contact AUROC on runs meeting completeness only | 0.748 [0.730, 0.763] | 0.500 | 0.00023 | 191 variants |
| C4 | 8M | IG contact AUROC, trained vs random-init model (same variants) | 0.717 [0.689, 0.745] | 0.460 | 0.00023 | 50 variants |
| C4 | 8M | share of top-10 IG residues that are 3D contacts | 0.309 [0.290, 0.328] |  |  | 200 variants |
| C4 | 8M | fraction of IG runs meeting completeness (<=10% gap) | 0.955 [0.925, 0.980] |  |  | 200 variants |
| C4 | 8M | Spearman(|IG|, |occlusion|) per variant | 0.300 [0.285, 0.316] | 0.000 | 0.00023 | 100 variants |
| C4 | 8M | occlusion AUROC for 3D contacts | 0.769 [0.749, 0.789] | 0.500 | 0.00023 | 100 variants |
| C5 | 8M | |dLLR| masking top-1 IG residues vs 1 random | 0.330 [0.278, 0.383] | 0.028 | 0.00023 | 200 variants |
| C5 | 8M | |dLLR| masking top-1 IG residues vs 1 sequence-nearest | 0.330 [0.278, 0.383] | 0.213 | 0.00023 | 200 variants |
| C5 | 8M | |dLLR| masking top-2 IG residues vs 2 random | 0.410 [0.342, 0.484] | 0.049 | 0.00023 | 200 variants |
| C5 | 8M | |dLLR| masking top-2 IG residues vs 2 sequence-nearest | 0.410 [0.342, 0.484] | 0.317 | 0.0024 | 200 variants |
| C5 | 8M | |dLLR| masking top-4 IG residues vs 4 random | 0.511 [0.426, 0.599] | 0.081 | 0.00023 | 200 variants |
| C5 | 8M | |dLLR| masking top-4 IG residues vs 4 sequence-nearest | 0.511 [0.426, 0.599] | 0.438 | 0.07 | 200 variants |
| C5 | 8M | |dLLR| masking top-8 IG residues vs 8 random | 0.602 [0.506, 0.703] | 0.127 | 0.00023 | 200 variants |
| C5 | 8M | |dLLR| masking top-8 IG residues vs 8 sequence-nearest | 0.602 [0.506, 0.703] | 0.614 | 0.81 | 200 variants |
| C5 | 8M | |dLLR| masking top-16 IG residues vs 16 random | 0.754 [0.629, 0.886] | 0.190 | 0.00023 | 200 variants |
| C5 | 8M | |dLLR| masking top-16 IG residues vs 16 sequence-nearest | 0.754 [0.629, 0.886] | 0.838 | 0.11 | 200 variants |
| C4 | 35M | IG |attribution| AUROC for 3D contacts of the mutated site | 0.741 [0.724, 0.758] | 0.500 | 0.00023 | 200 variants |
| C4 | 35M | IG AUROC for long-range contacts (|i-j|>=6) | 0.597 [0.569, 0.626] | 0.500 | 0.00023 | 182 variants |
| C4 | 35M | robustness: IG contact AUROC on runs meeting completeness only | 0.751 [0.728, 0.774] | 0.500 | 0.00023 | 100 variants |
| C4 | 35M | IG contact AUROC, trained vs random-init model (same variants) | 0.732 [0.705, 0.758] | 0.517 | 0.00023 | 50 variants |
| C4 | 35M | share of top-10 IG residues that are 3D contacts | 0.308 [0.288, 0.330] |  |  | 200 variants |
| C4 | 35M | fraction of IG runs meeting completeness (<=10% gap) | 0.500 [0.435, 0.570] |  |  | 200 variants |
| C4 | 35M | Spearman(|IG|, |occlusion|) per variant | 0.259 [0.239, 0.278] | 0.000 | 0.00023 | 100 variants |
| C4 | 35M | occlusion AUROC for 3D contacts | 0.778 [0.759, 0.798] | 0.500 | 0.00023 | 100 variants |
| C5 | 35M | |dLLR| masking top-1 IG residues vs 1 random | 0.370 [0.310, 0.442] | 0.038 | 0.00023 | 200 variants |
| C5 | 35M | |dLLR| masking top-1 IG residues vs 1 sequence-nearest | 0.370 [0.310, 0.442] | 0.275 | 0.0042 | 200 variants |
| C5 | 35M | |dLLR| masking top-2 IG residues vs 2 random | 0.520 [0.434, 0.615] | 0.062 | 0.00023 | 200 variants |
| C5 | 35M | |dLLR| masking top-2 IG residues vs 2 sequence-nearest | 0.520 [0.434, 0.615] | 0.456 | 0.13 | 200 variants |
| C5 | 35M | |dLLR| masking top-4 IG residues vs 4 random | 0.661 [0.563, 0.764] | 0.098 | 0.00023 | 200 variants |
| C5 | 35M | |dLLR| masking top-4 IG residues vs 4 sequence-nearest | 0.661 [0.563, 0.764] | 0.655 | 0.92 | 200 variants |
| C5 | 35M | |dLLR| masking top-8 IG residues vs 8 random | 0.947 [0.809, 1.073] | 0.148 | 0.00023 | 200 variants |
| C5 | 35M | |dLLR| masking top-8 IG residues vs 8 sequence-nearest | 0.947 [0.809, 1.073] | 1.076 | 0.1 | 200 variants |
| C5 | 35M | |dLLR| masking top-16 IG residues vs 16 random | 1.365 [1.155, 1.566] | 0.228 | 0.00023 | 200 variants |
| C5 | 35M | |dLLR| masking top-16 IG residues vs 16 sequence-nearest | 1.365 [1.155, 1.566] | 1.654 | 0.0023 | 200 variants |
| C6 | 35M | held-out Spearman drop: ablate top-1 heads vs 1 random heads | 0.007 [-0.010, 0.033] | 0.002 | 0.94 | 5 proteins |
| C6 | 35M | proteins where top-1 heads beat random heads (sign test p, one-sided) | 2.000  | 2.500 | 0.84 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-1 contact heads | -0.004 [-0.013, 0.004] | 0.002 | 0.57 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-2 heads vs 2 random heads | 0.045 [0.013, 0.094] | 0.004 | 0.11 | 5 proteins |
| C6 | 35M | proteins where top-2 heads beat random heads (sign test p, one-sided) | 5.000  | 2.500 | 0.061 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-2 contact heads | -0.005 [-0.012, 0.001] | 0.004 | 0.18 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-4 heads vs 4 random heads | 0.067 [0.033, 0.133] | 0.010 | 0.11 | 5 proteins |
| C6 | 35M | proteins where top-4 heads beat random heads (sign test p, one-sided) | 5.000  | 2.500 | 0.061 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-4 contact heads | 0.000 [-0.004, 0.006] | 0.010 | 0.11 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-8 heads vs 8 random heads | 0.114 [0.041, 0.211] | 0.021 | 0.11 | 5 proteins |
| C6 | 35M | proteins where top-8 heads beat random heads (sign test p, one-sided) | 5.000  | 2.500 | 0.061 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-8 contact heads | 0.004 [-0.006, 0.014] | 0.021 | 0.18 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-16 heads vs 16 random heads | 0.164 [0.061, 0.274] | 0.033 | 0.18 | 5 proteins |
| C6 | 35M | proteins where top-16 heads beat random heads (sign test p, one-sided) | 4.000  | 2.500 | 0.25 | 5 proteins |
| C6 | 35M | held-out Spearman drop: ablate top-16 contact heads | 0.027 [-0.002, 0.056] | 0.033 | 0.84 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-1 heads vs 1 random heads | 0.026 [-0.002, 0.056] | 0.000 | 0.31 | 5 proteins |
| C6 | 8M | proteins where top-1 heads beat random heads (sign test p, one-sided) | 3.000  | 2.500 | 0.57 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-1 contact heads | -0.004 [-0.014, 0.005] | 0.000 | 0.57 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-2 heads vs 2 random heads | 0.027 [-0.003, 0.064] | -0.000 | 0.31 | 5 proteins |
| C6 | 8M | proteins where top-2 heads beat random heads (sign test p, one-sided) | 4.000  | 2.500 | 0.25 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-2 contact heads | 0.004 [-0.015, 0.018] | -0.000 | 0.84 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-4 heads vs 4 random heads | 0.069 [0.022, 0.115] | 0.004 | 0.18 | 5 proteins |
| C6 | 8M | proteins where top-4 heads beat random heads (sign test p, one-sided) | 4.000  | 2.500 | 0.25 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-4 contact heads | 0.021 [-0.013, 0.060] | 0.004 | 0.63 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-8 heads vs 8 random heads | 0.079 [-0.010, 0.202] | 0.004 | 0.38 | 5 proteins |
| C6 | 8M | proteins where top-8 heads beat random heads (sign test p, one-sided) | 4.000  | 2.500 | 0.25 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-8 contact heads | 0.037 [0.003, 0.074] | 0.004 | 0.46 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-16 heads vs 16 random heads | 0.102 [-0.003, 0.223] | 0.031 | 0.31 | 5 proteins |
| C6 | 8M | proteins where top-16 heads beat random heads (sign test p, one-sided) | 3.000  | 2.500 | 0.57 | 5 proteins |
| C6 | 8M | held-out Spearman drop: ablate top-16 contact heads | 0.015 [-0.031, 0.075] | 0.031 | 0.75 | 5 proteins |
| C7 | 35M | layer 0: Spearman as a fraction of the final layer's | 0.322 [0.235, 0.414] |  |  | 52 assays |
| C7 | 35M | layer 1: Spearman as a fraction of the final layer's | 0.397 [0.270, 0.529] |  |  | 52 assays |
| C7 | 35M | layer 2: Spearman as a fraction of the final layer's | 0.420 [0.265, 0.572] |  |  | 52 assays |
| C7 | 35M | layer 3: Spearman as a fraction of the final layer's | 0.080 [-0.038, 0.182] |  |  | 52 assays |
| C7 | 35M | layer 4: Spearman as a fraction of the final layer's | 0.142 [0.060, 0.239] |  |  | 52 assays |
| C7 | 35M | layer 5: Spearman as a fraction of the final layer's | -0.024 [-0.102, 0.062] |  |  | 52 assays |
| C7 | 35M | layer 6: Spearman as a fraction of the final layer's | -0.113 [-0.223, -0.014] |  |  | 52 assays |
| C7 | 35M | layer 7: Spearman as a fraction of the final layer's | -0.020 [-0.132, 0.088] |  |  | 52 assays |
| C7 | 35M | layer 8: Spearman as a fraction of the final layer's | 0.202 [0.112, 0.294] |  |  | 52 assays |
| C7 | 35M | layer 9: Spearman as a fraction of the final layer's | 0.354 [0.278, 0.433] |  |  | 52 assays |
| C7 | 35M | layer 10: Spearman as a fraction of the final layer's | 0.524 [0.427, 0.616] |  |  | 52 assays |
| C7 | 35M | layer 11: Spearman as a fraction of the final layer's | 0.780 [0.714, 0.852] |  |  | 52 assays |
| C7 | 35M | layer 12: Spearman as a fraction of the final layer's | 1.000 [1.000, 1.000] |  |  | 52 assays |
| C7 | 8M | layer 0: Spearman as a fraction of the final layer's | 0.155 [0.047, 0.272] |  |  | 46 assays |
| C7 | 8M | layer 1: Spearman as a fraction of the final layer's | -0.609 [-0.800, -0.429] |  |  | 46 assays |
| C7 | 8M | layer 2: Spearman as a fraction of the final layer's | 0.045 [-0.061, 0.146] |  |  | 46 assays |
| C7 | 8M | layer 3: Spearman as a fraction of the final layer's | 0.252 [0.120, 0.432] |  |  | 46 assays |
| C7 | 8M | layer 4: Spearman as a fraction of the final layer's | 0.128 [-0.005, 0.286] |  |  | 46 assays |
| C7 | 8M | layer 5: Spearman as a fraction of the final layer's | 0.498 [0.400, 0.603] |  |  | 46 assays |
| C7 | 8M | layer 6: Spearman as a fraction of the final layer's | 1.000 [1.000, 1.000] |  |  | 46 assays |
