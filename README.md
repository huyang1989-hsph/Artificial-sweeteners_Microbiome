# Artificial sweetener intake, human gut microbiome, plasma metabolome, and cardiometabolic disease risk in free-living individuals

This repository contains the analysis code accompanying the manuscript **“Artificial sweetener intake, human gut microbiome, plasma metabolome, and cardiometabolic disease risk in free-living individuals,”**.

The analyses integrate repeated dietary assessments, gut microbiome profiles, plasma metabolomics, participant characteristics, and prospective cardiometabolic outcomes. They include external replication in the MicroCardio Consortium, sensitivity analyses using centered log-ratio (CLR) transformation, comparisons of zero-replacement fractions, and assessments of statistical stability in feature selection and sparse canonical correlation analysis (sCCA).

## Sweeteners evaluated

The reported analyses include nine individual artificial sweeteners:

- **Synthetic sweeteners:** acesulfame K, aspartame, saccharin, and sucralose.
- **Sugar alcohols:** erythritol, maltitol, mannitol, sorbitol, and xylitol.

The corresponding prefixes used in the code are `acesk_`, `aspart_`, `sach_`, `sucral_`, `eryth_`, `maltitol_`, `mani_`, `sorb_`, and `xyli_`. 

## Code files

### Primary analysis, data preparation, and replication

| File | Purpose | Main output |
| --- | --- | --- |
| `LVS_7DDR.sas` | Reads and cleans the 7-day dietary-record data, aligns daily records across two recording weeks, merges nutrient and food variables, and prepares sweetener intake variables. | Processed dietary-record data. |
| `ddr.R` | Integrates dietary variables, biomarkers, covariates, fecal-sample timing, and microbiome profiles. | `TAXON_ONLY.RData`. |
| `Taxon_AVG.R` | Averages repeated measurements by participant, recodes selected indicators, and defines analysis phenotypes. | `Taxon_use_AVG.RData`. |
| `ELN.R` | Original elastic-net regression code for species feature selection across ten mixing-parameter values (`alpha = 0.1, 0.2, ..., 1.0`), with leave-one-out cross-validation for lambda tuning. | Elastic-net model results used in the primary analysis. |
| `batch_0.01%.R` | Runs the elastic-net feature-selection pipeline for the original arcsine square-root (AST)–transformed microbiome analysis across the ten alpha settings. | Ten `INT_alpha_*_0.3.RData` files containing AST model results. |
| `figures.R` | Combines elastic-net features across alpha settings, builds microbial scores, performs community-level, sCCA, metabolite, pathway, and enzyme analyses, and generates primary visualization panels. | Figures 2, 3, and 5; Extended Data Figures 1 and 3; AST reference results for sensitivity comparisons. |
| `surv.R` | Conducts Cox proportional-hazards and inverse-probability-weighted survival analyses of metabolite composite scores in relation to incident type 2 diabetes (T2D) and coronary heart disease (CHD). | Figure 6. |
| `replicate.R` | Summarizes and meta-analyzes study-specific T2D odds ratios from the MicroCardio Consortium and generates forest plots and sensitivity/heterogeneity summaries. | Figure 4; Extended Data Figure 2; Extended Data Table 1. |

### CLR sensitivity analyses and statistical stability

| File | Purpose | Main output or figure contribution |
| --- | --- | --- |
| `batch_CLR.R` | Runs elastic-net feature selection for the original CLR-based analysis across the ten alpha settings. | Ten CLR-specific `INT_alpha_*_0.3.RData` files. |
| `CLR_ELN.R` | Original CLR sensitivity-analysis script for reconstructing microbial scores, fitting sCCA, and assessing metabolomic associations. | Extended Data Figures 4–6; Extended Data Table 2. |
| `CLR_ELN_update.R` | Updates the CLR sensitivity analysis with species structure loadings, metabolite comparisons, and statistical-stability assessments. Sources `ELN_LOOCV_stability.R`, `verify_SCCA_LOO_stability.R`, and `SCCA_LOO_stability.R`. | CLR analysis outputs used with AST reference outputs to visualize Supplemental Figures 1–3 and 6–11; updated CLR sensitivity results. |
| `ELN_LOOCV_stability.R` | Defines helpers for leave-one-participant-out elastic-net refitting at the original alpha and previously tuned lambda settings. Checks full-data refits against saved model selections and coefficients before starting the omission analyses. | Selection frequencies, species retention summaries, and per-refit results for Supplemental Figures 8 and 9. |
| `SCCA_LOO_stability.R` | Defines helpers for leave-one-participant-out sCCA refitting with fixed candidate species and penalties. Aligns the first two components to their reference counterparts and summarizes canonical correlations, species retention, loadings, profiles, and score agreement. | Component and species stability outputs for Supplemental Figures 10 and 11. |
| `verify_SCCA_LOO_stability.R` | Verifies the sCCA stability helper using simulated data, including component order/sign alignment, participant alignment, direct-refit agreement, and detection of mismatched reference inputs. | Verification checks; no participant-level study data are required for this script. |
| `batch_CLR_0.40.R` | Runs CLR elastic-net feature selection using a zero-replacement fraction of 0.40. | Model results in `zero_0.40/`; inputs to `CLR_ELN_0.4.R`. |
| `CLR_ELN_0.4.R` | Loads the 0.40 elastic-net results and performs the subsequent score-building and sCCA pipeline using the same replacement fraction. | `ELN_CLR_0.40.csv` and `CLR_sCCA_coef_0.40.csv`; inputs for Supplemental Figures 4 and 5. |
| `batch_CLR_0.90.R` | Runs CLR elastic-net feature selection using a zero-replacement fraction of 0.90. | Model results in `zero_0.90/`; inputs to `CLR_ELN_0.9.R`. |
| `CLR_ELN_0.9.R` | Loads the 0.90 elastic-net results and performs the subsequent score-building and sCCA pipeline using the same replacement fraction. | `ELN_CLR_0.90.csv` and `CLR_sCCA_coef_0.90.csv`; inputs for Supplemental Figures 4 and 5. |

The stability helper files define functions; sourcing them alone does not run the study analyses. `verify_SCCA_LOO_stability.R` is an exception: sourcing it runs the simulated-data verification. Calls using the corresponding AST or CLR data and saved reference fits are required to generate study stability results. The supplied `CLR_ELN_update.R` includes the CLR calls; the same helpers are used with AST inputs for the primary-analysis stability assessment.

## Supplemental figure mapping

| Supplemental figure(s) | Analysis summarized | Analysis inputs |
| --- | --- | --- |
| 1 | AST–CLR comparisons of canonical correlations and correlations between individual sweetener intakes and microbial component scores. | Original AST results and updated CLR results from `CLR_ELN_update.R`. |
| 2–3 | Agreement in species structure loadings, loading directions, and loading magnitudes across AST and CLR analyses. | AST and CLR structure loadings and sCCA species-selection results; CLR exports are generated by `CLR_ELN_update.R`. |
| 4–5 | Elastic-net species selection and sCCA canonical-weight comparisons across zero-replacement fractions 0.40, 0.65 (reference), and 0.90. | The paired 0.40 and 0.90 batch/downstream scripts and the reference CLR analysis. |
| 6–7 | Direct AST–CLR comparisons of metabolite association coefficients, directions, and false discovery rate (FDR)–adjusted P values. | Original AST metabolite results and CLR metabolite outputs from `CLR_ELN_update.R`. |
| 8–9 | Leave-one-participant-out elastic-net selection stability for AST and CLR, respectively. | Outputs from `ELN_LOOCV_stability.R` using the respective data and saved models. |
| 10–11 | Leave-one-participant-out sCCA component stability for AST and CLR, respectively. | Outputs from `SCCA_LOO_stability.R` using the respective data and reference fits. |

These scripts generate the analysis results used for the supplemental visualizations. Final plotting and assembly are separate visualization steps. Working image and CSV filenames within the scripts do not necessarily match the final manuscript figure numbers.

## Analysis workflow

### 1. Prepare the shared analysis data

1. Run `LVS_7DDR.sas` to prepare the dietary-record data.
2. Run `ddr.R` to integrate dietary, biomarker, covariate, fecal-timing, and microbiome data and create `TAXON_ONLY.RData`.
3. Run `Taxon_AVG.R` to average repeated measurements by participant and create `Taxon_use_AVG.RData`.

### 2. Run the primary AST analysis

1. Run `batch_0.01%.R` to perform AST-based elastic-net feature selection across the ten alpha settings. `ELN.R` contains the original elastic-net modeling code.
2. Run `figures.R` with the AST model-result files to construct scores and generate the primary analysis results and figure panels.
3. Run `surv.R` after the required metabolite-score and outcome inputs are available.
4. Retain the AST model objects, sCCA reference fit, structure loadings, and metabolite association results for the AST–CLR comparisons and stability assessments.

### 3. Run the reference CLR analysis and updated comparisons

1. Run `batch_CLR.R` to generate the original CLR elastic-net model-result files.
2. Run `CLR_ELN_update.R` with the CLR model results, shared analysis data, metabolomics objects, and required AST reference outputs available. This is the updated analysis entry point; `CLR_ELN.R` is retained as the original CLR sensitivity-analysis script.
3. Ensure that all three stability helper scripts are available at the paths used by the `source()` calls. The updated script sources the verification and fitting helpers and invokes the CLR stability analyses.
4. Apply the stability helpers to the corresponding AST inputs and saved reference fits to obtain the AST stability outputs.
5. Use the exported results and AST references for Supplemental Figures 1–3 and 6–11.

For elastic-net stability, each participant is omitted in turn, models are refitted at all ten alpha settings with previously tuned lambda values held fixed, and the original signature criterion is reapplied: a species must have a nonzero coefficient at every alpha setting. For sCCA stability, candidate species and penalty parameters are held fixed, and the two retained components are aligned in order and sign to the full-data reference. These assessments quantify stability conditional on the supplied data transformations, candidate-species sets, and tuning choices.

The helpers compare a full-data refit with the saved reference before proceeding. Both `eln_loo_from_saved()` and `scca_loo_stability()` support `check_only = TRUE` to inspect reference checks before running the leave-one-out refits. The verification script can also be run independently using:

```r
source("verify_SCCA_LOO_stability.R")
```

### 4. Run the zero-replacement sensitivity analyses

| Fraction | Feature selection: run first | Score building and sCCA: run next | Batch result directory |
| --- | --- | --- | --- |
| 0.40 | `batch_CLR_0.40.R` | `CLR_ELN_0.4.R` | `zero_0.40/` |
| 0.90 | `batch_CLR_0.90.R` | `CLR_ELN_0.9.R` | `zero_0.90/` |

Create the fraction-specific result directories before running the batch scripts. Each downstream script loads the ten model-result files from its matching directory. Both stages of each alternative branch explicitly use the matching replacement fraction. Compare their species-selection and canonical-weight exports with the reference CLR results to prepare Supplemental Figures 4 and 5.

### 5. Run the replication analysis

Run `replicate.R` with the study-specific MicroCardio odds-ratio input files. This branch uses summary estimates and is independent of the primary individual-level analysis.

## Key exported results

- **Elastic-net batches:** ten `INT_alpha_*_0.3.RData` files for each analysis branch.
- **Updated CLR analysis:** `species_structure_loadings_CLR.csv`, `CLR_met_M1.csv`, and `CLR_met_M2.csv`, along with the fitted scores, correlations, and species-weight outputs used in comparative visualizations.
- **Elastic-net stability:** `full_data_reference_check.csv`, `LOO_selection_summary.csv`, outcome-specific species-frequency and alpha-specific frequency CSV files, per-refit summaries, and RDS results. The CLR call uses `LOO_CLR_all_sweeteners_v3/`.
- **sCCA stability:** `SCCA_LOO_summary.csv`, `SCCA_LOO_component_overview.csv`, `SCCA_LOO_fold_results.csv`, species and sweetener-profile summaries, and `SCCA_LOO_results.rds`. The CLR call uses `SCCA_LOO_CLR/`.
- **Alternative zero-replacement analyses:** `ELN_CLR_0.40.csv`, `CLR_sCCA_coef_0.40.csv`, `ELN_CLR_0.90.csv`, and `CLR_sCCA_coef_0.90.csv`.

## Required data and intermediate objects

Individual-level study data are not included in this code collection. The scripts expect locally available source data and intermediate objects, including:

- Participant characteristics and biospecimen information (`Basic.csv`).
- Biomarker data (`Biomarkers.csv`).
- Food-frequency questionnaire data (`FFQD.csv`).
- Raw and derived 7-day dietary-record datasets used by `LVS_7DDR.sas`.
- The processed SAS dietary-record file (`mlvs_ddr.sas7bdat`).
- Taxonomic abundance and annotation objects (`ReadIn_taxon.RData` and `TAXON_ONLY.RData`).
- Participant-level averaged data (`Taxon_use_AVG.RData`).
- Plasma metabolomics objects (`MET_pm.RData` and `MET_fea.RData`).
- Saved AST and CLR elastic-net models and the corresponding sCCA reference fits for stability assessment.
- AST reference coefficients, species structure loadings, and metabolite association results for comparison with CLR.
- Pooled outcome data for survival analysis (`AS_outcome.csv`).
- MicroCardio study-specific odds-ratio, sex-stratified, sample-size, and sensitivity-analysis CSV files used by `replicate.R`.

## Software requirements

- **SAS** is required for `LVS_7DDR.sas`.
- **R** is required for all `.R` scripts.

Major R packages include `tidyverse`, `data.table`, `caret`, `glmnet`, `vegan`, `CCA`, `CCP`, `PMA`, `ComplexHeatmap`, `zCompositions`, `survival`, `survminer`, `mice`, `meta`, `ggplot2`, `patchwork`, and `cowplot`. Additional packages are loaded within individual scripts; consult their `library()` calls for the complete list. Some original scripts also reference the project-specific `chanmetab` package.

A complete pinned software environment is not supplied. The stability helpers retain selected version information, including `caret`/`glmnet` versions for elastic-net refitting and R/`PMA` versions for sCCA refitting. Record `sessionInfo()` for the complete R environment and retain the SAS version and configuration used for dietary-data preparation.

## Running and adapting the scripts

The scripts preserve the directory structure and interactive workflow used for the analyses. They contain project-specific paths, calls to `setwd()`, and references to saved or in-memory objects. Before running them:

1. Place the authorized input data and intermediate objects in the expected locations, then update file paths and working directories for the local environment, including paths in `source()` calls.
2. Keep the AST and CLR batch results in separate analysis directories: both original batch scripts use the same `INT_alpha_*_0.3.RData` filenames. Keep the alternative fractions in their respective `zero_0.40/` and `zero_0.90/` directories.
3. Ensure that the required project-specific helpers, packages, and reference model objects are available and that participant IDs and feature names match the saved fits.
4. In the elastic-net stability section of the supplied `CLR_ELN_update.R`, the model list is created as `eln_results` but passed to the helper as `eln_results_CLR`. Use the same name in the assignment and call, or assign `eln_results_CLR <- eln_results` immediately after constructing the list.
5. Preserve the zero-replacement implementation used for each reference fit. 
6. Follow the relevant workflow above and inspect the full-data reference checks before interpreting stability summaries.

The scripts are provided as an analytic record rather than a fully self-contained, one-command workflow. Access to the underlying participant-level data is governed by the applicable cohort consent, institutional approvals, and data-use agreements.

## Citation

Please cite the associated manuscript when using this code. Full citation information will be added when the article is published.
