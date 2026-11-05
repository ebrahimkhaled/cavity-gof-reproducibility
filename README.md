# Reproduction materials: a leave-one-out goodness-of-fit test for penalised logistic regression in high dimensions

Ebrahim Khaled Ebrahim (Alexandria University) and Mohammad Arashi (Ferdowsi University of Mashhad).

This repository holds the code, the simulation outputs and the figures behind every number, table and
figure of the article and its online supplement. Every simulation writes each replicate to disk, so the
`results/` folder contains the per-replicate p-values as well as the summaries.

## Layout

| Folder | Contents |
|---|---|
| `R/cavity_engine.R` | The cavity test: the one-step leave-one-out (cavity) index, the two- and three-parameter recalibration, the grouped statistic and its reference law (`cav_gof_fast` for ridge, `cav_gof_glmnet` for the lasso and the elastic net, `cav_split_gof` for the sample-split and cross-fitted tests, `cav_from_index` for any externally computed index). |
| `R/cluster_correction.R` | The twin-share diagnostic (`cav_clusters`, `cav_twin_share`) and the cluster correction tried in Supplementary Section S9.4 (`cav_from_index_cluster`). |
| `R/closedform_engine.R` | The comparator called the *in-sample statistic* in the article: the shrinkage-corrected decile statistic of the companion paper, with its closed-form reference (function `calm.gof`). It also holds the ridge fitter `cf_ridge` used by every script. |
| `sims/` | One script per simulation study (`E18` to `E44`), the companion-paper data sets reused for power (`E4b`, `E4b2`), and short checks. |
| `results/` | Outputs: `*_perrep.csv` (one row per replicate), `*_summary.csv` and the console `*.log` of each run. |
| `figures/` | `make_figures.R`, the figures of the article as vector PDF, and the supplementary QQ plot (`figS_qq.pdf`). |
| `data_raw/internet_ads/` | The Internet Advertisements data of the UCI Machine Learning Repository (Kushmerick, 1999), as distributed by UCI under CC BY 4.0. |

## Where each result comes from

| Article item | Script(s) | Output(s) |
|---|---|---|
| Figure 1 | `figures/make_figures.R` (panel A simulated inside the script; panel B from `E23_cavity_verify_prop1.R`) | `E23_verify_prop1_perobs_seed20261103.csv` |
| Figure 2 | `E24_selection_formula.R` | `E24_selection_formula_pergroup.csv` |
| Figure 3 (two-parameter version; the default is in Supplementary Section S4) | `E18_cavity_size.R`, `E32_textbook_on_E18.R` (panel A); `E19_cavity_stress.R` (B); `E22_cavity_lasso.R` (C) | `E32_*`, `E19_*`, `E22_*` |
| Figure 4 (two-parameter version; the default is in Supplementary Section S9) | `E30_power_strong.R` | `E30_power_strong_*` |
| Figure 5 and Section 7 (glaucoma, leave-one-out fits with standardisation inside each fit) | `E36_glaucoma_noleak.R` (P-values, deciles and figure), `E36b_glaucoma_cv_seeds.R` (cross-validated penalty over 50 fold seeds), `E36c_glaucoma_penalty_curve.R` (penalty grid) | `E36_*`, `E36b_*`, `E36c_*`; `figures/fig5_glaucoma_v2.pdf` |
| Table 1 (the algorithm) | `R/cavity_engine.R` | — |
| Section 4.1 and Supplementary Section S12 (check of the calibration law) | `E23_cavity_verify_prop1.R`, `E23_decompose.R` | `E23_*` |
| Supplementary Section S1 (size in designs A and B and seven further designs) | `E18_cavity_size.R`, `E32_textbook_on_E18.R` | `E18_*`, `E32_*` |
| Supplementary Section S2 (designs chosen to break the test) | `E19_cavity_stress.R`, `E28_fixedp_large_n.R` | `E19_*`, `E28_*` |
| Supplementary Section S3 (lasso and elastic net) | `E22_cavity_lasso.R`, `E26_confirm_borderline.R` | `E22_*`, `E26_*` |
| Supplementary Section S4 and Sections 5.3 and 6.2 (the three-parameter default: size at 1, 5 and 10 per cent in all 32 designs, QQ plot) | `E38_q3_size_all.R`, `E38_summarise.R`; `E42_q3_repair.R` (refits the few logistic-normal fits whose information matrix was singular) | `E38_*`, `E42_*`; `figures/figS_qq.pdf` |
| Supplementary Section S5 (covariance of the grouped residuals) | `E18_covariance_check.R` | reads `E18_cavity_size_perrep.csv` |
| Supplementary Section S6 and Section 5.3 (two or three parameters) | `E25_q2_vs_q3.R`, `E29_q3_excess_diagnosis.R`, `q3_reference_check.R` | `E25_*`, `E29_*` |
| Supplementary Section S7 (sample-split and cross-fitted tests) | `E31_split_crossfit.R` (two parameters); `E43_split_q3.R`, `E43_summarise.R` (three parameters, same data and splits) | `E31_*`, `E43_*` |
| Supplementary Section S8 (penalty chosen by cross-validation) | `E27_cv_penalty.R` | `E27_*` |
| Supplementary Section S9 (power) | `E30_power_strong.R`; `E20_cavity_power.R` on the data sets of `E4b_final_design_grid.R` and `E4b2_calm_final_on_E4b.R` | `E30_*`, `E20_*`, `E4b_*`, `E4b2_*` |
| Supplementary Section S9 (power of the default and the full parametric bootstrap on the strong departures) | `E39_power_q3_boot.R`, `E39_summarise.R` | `E39_*` |
| Supplementary Section S9 (identification of the noise parameter tau) | `E41_tau_identification.R` | `E41_*` |
| Supplementary Section S10.1 (K-fold and residual-prediction comparators) | `E33_comparators.R` (needs the CRAN package `GRPtests`) | `E33_*` |
| Supplementary Section S10.2 and Section 6.4 (the residual-prediction test in sparse designs) | `E40_grp_sparse.R`, `E40_summarise.R` (needs the CRAN package `GRPtests`) | `E40_*` |
| Supplementary Section S10.3 (sparse binary and skewed predictors) | `E34_sparse_binary_stress.R`, `E34b_confirm_rare_binary.R` | `E34_*`, `E34b_*` |
| Supplementary Section S10.3 (the three-parameter default on the sparse binary designs) and S10.5 (Internet Advertisements placebo) | `E44_q3_binary_ads.R`, `E44_summarise.R` | `E44_*` |
| Supplementary Section S10.4 (placebo checks on real predictors) | `E35_real_data_placebo.R` | `E35_*` |
| Supplementary Section S10.4 and Section 7 (glaucoma placebo with standardisation inside each fit) | `E37_glaucoma_placebo_noleak.R` | `E37_*` |
| Supplementary Section S10.5 and Section 8 (Internet Advertisements: the placebo failure, its cause and the twin-share diagnostic) | `E35_real_data_placebo.R` (placebo); `E35b_ads_placebo_diagnosis.R` (exact leave-one-out and ten-fold check); `E35c_ads_split_placebo.R` (split test); `E35d_cluster_correction.R` (cluster correction); `E35e_twin_share.R` (diagnostic and simulation-calibrated P); `E35f_ads_misfit_decomposition.R` (misfit and noise parts, kernel correction); `E35g_ads_conditional_calibration.R` (size of the reference simulated on the observed design) | `E35_*`, `E35b_*` to `E35g_*`, `E35e_twin_share.txt` |
| Supplementary Section S11 (the self-selection formula) | `E24_selection_formula.R` | `E24_*` |
| Supplementary Section S13 (the logistic-normal gap, the eigenvalues of the covariance and the conservative reference) | `logistic_normal_gap.R`, `q3_eigen_check.R`, `chi2_dominance_check.R`, `chi2_dominance_three_weights.R`, `q3_reference_check.R` | console |
| Section 5.2 and Supplementary Section S14 (dependence between leave-one-out indices, before and after the projection) | `dependence_check.py` (Python 3 with NumPy; arguments: seed, penalty, replications, output file, n, p), summarised by `dependence_analyse.py` and `dependence_excess_se.py` | console |

## How to run

Install R (version 4.4.1 was used) and the CRAN packages `CompQuadForm`, `glmnet`, `GRPtests`, `TH.data`
and `ebrahim.gof` (the last is needed only for the parametric-bootstrap comparator in `E4b`). Run every
script from the repository root with `Rscript`, for example

```
Rscript sims/E18_cavity_size.R 500 12
Rscript figures/make_figures.R
```

The first argument is the number of replications and the second the number of parallel workers; the
number used for each published result is recorded in its per-replicate file and log. Seeds are fixed
inside each script. A run that is interrupted resumes from the replicates or cells already written to
`results/`. The scripts locate the repository from their own path, so no directory needs to be set.
`sessionInfo.txt` lists the package versions.

## Licence

Code is released under the MIT licence (`LICENSE`). Simulation outputs and figures are released under
CC BY 4.0 (`LICENSE-DATA.txt`). The glaucoma data are those of the CRAN package `TH.data` and are not
redistributed here. The Internet Advertisements data are redistributed unchanged from the UCI Machine
Learning Repository under CC BY 4.0 (Kushmerick, N. (1999). Learning to remove Internet advertisements.
Proceedings of the Third Annual Conference on Autonomous Agents, 175–181.
https://doi.org/10.1145/301136.301186; data DOI https://doi.org/10.24432/C5V011).

## Changes

- 1.3.0: the three-parameter recalibration is the default of the article. New studies E36 to E44: the glaucoma analysis redone with standardisation inside each leave-one-out fit (E36, E36b, E36c) and its placebo (E37); the size of the default at three levels in all 32 designs (E38, with E42 for singular fits); power of the default and the full parametric bootstrap (E39); the residual-prediction test in sparse designs (E40); identification of tau (E41); the sample-split tests with three parameters (E43); the default on the sparse binary designs and the Internet Advertisements placebo (E44). Two new figures; supplementary sections renumbered.
- 1.2.0: the misfit decomposition and kernel correction (E35f) and the size of the reference simulated
  on the observed design (E35g; stopped after 24 of 60 outer data sets, summary computed from the
  per-replicate file); the twin-share and cluster-correction functions in `R/cluster_correction.R`
  gain the kernel version.
- 1.1.0: new title of the article; further comparators (E33), sparse binary and skewed designs (E34,
  E34b), placebo checks on real predictors (E35), the Internet Advertisements boundary case with its diagnosis
  and the twin-share diagnostic (E35b to E35e, `R/cluster_correction.R`; the E35d run was stopped after
  24 of its 60 power replications once the size arm had shown that the correction fails, so its log ends
  there and `E35d_cluster_summary.csv` was computed from the per-replicate file), the eigenvalue
  and three-weight checks of the reference law; supplementary sections renumbered.
- 1.0.0: first release.
