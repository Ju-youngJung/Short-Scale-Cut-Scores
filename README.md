# Enumeration-Based Evaluation of Cut Score Decisions in Short Tests — analysis package

One-shot reproduction of every number, figure, and supplementary table in the manuscript.

## Run (R >= 4.3; Windows, macOS, or Linux)
1. Unzip; the folder contains `00_setup.R`, `R/`, `scripts/`, `docs/`.
2. Restart R, then:
   ```r
   setwd("<this folder>")
   source("scripts/run_everything.R")
   ```
3. Wait. Exact study ~15 min (9 cores), ICAR ~2 min, ICAR calibration bootstrap ~30–60 min
   (set `B_BOOT <- 0` in `scripts/run_everything.R` to skip it, or `20` for a quicker run), supplement ~2 min.

## Outputs
- `results/final_v2/`: `ALL_summary.csv`, `ALL_marginal_summary.csv`, `ALL_screen_size.csv`, `ALL_cell_info.csv`,
  `REPORT_final.md`, `run_log.txt` (implementation checks), `protocol_snapshot.txt`, `figures/G1–G4_*.png`.
- `results/empirical_icar/`: `icar_bank_2pl.csv`, `ALL_summary.csv`, `ALL_marginal_summary.csv`,
  `respondent_scores.csv`, `bootstrap_calibration.csv`, `run_log.txt`, `figures/G1–G4_ICAR.png`.
- `results/Supplementary_Materials.xlsx` (README, A_Article_Values, V_Verification, S1–S22).

## Manuscript figures
Copy to the Overleaf folder `CUTSCORE/figures/`:
`results/final_v2/figures/G1_wrong_two_outcome_{wide,peaked}.png`, `G2_three_way_{wide,peaked}.png`,
`G3_coverage_up0.2_{wide,peaked}.png`, and `results/empirical_icar/figures/G{1,2,3}_*_ICAR.png`.

## Individual steps (optional)
`scripts/run_final_study.R` (exact study, resumable), `scripts/run_empirical_icar.R`,
`scripts/make_final_figures.R` (figures from `ALL_summary.csv`), `scripts/make_supplement_xlsx.R`.

## Checks
The exact study stops if any identity check fails (exponential-family identity, monotone likelihood ratio,
exact-interval coverage >= 0.95, Lord–Wingersky agreement, Huber moments, roots, mirt agreement).
Sheet `A_Article_Values` in the supplement recomputes every number quoted in the article from the output files.
