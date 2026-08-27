# Redistricting Ensemble Algorithm Comparison: SMC vs. Flip vs. Merge-Split

Compares three redistricting ensemble-generation algorithms: Sequential Monte Carlo (SMC), Flip MCMC, and Merge-Split MCMC, on New Mexico's 2020 congressional map, using the `redist` package. Each algorithm's ensemble is
checked for convergence, then compared on compactness, partisan symmetry metrics, and sampling efficiency.

## Files

- `NM Ensemble Comparison.R`: the ensemble comparison itself (SMC vs. Flip vs. Merge-split on New Mexico). See below.
- `Ensemble Generation Algorithm Figures.R`: builds the dual graph representation (`shp`, `g`, `coords`) from a small subset of Florida's precincts, draws the state-to-graph figure, and then draws step-by-step diagrams of how each of the three algorithms above actually operates. See "Figure Generation" below.

## Requirements

R packages: `redist`, `coda`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `alarmdata`.

`Ensemble Generation Algorithm Figures.R` has a slightly separate set of requirements. See "Figure Generation" below.

## Running

The script has two entry points, controlled by which block you run:

- **Fresh run**: run the whole script top to bottom. This calibrates each sampler (targeting a combined effective sample size of `target_ess`), runs the full ensembles, and saves all outputs (`.rds` files) to `out_dir`. Flip's full run in particular is computationally expensive (multi-hour on a laptop).
- **Resume from saved output**: if `.rds` files from a previous run are present, skip over "Run SMC Ensemble" through until "District-Level Metrics", resuming from the latter. Everything needed is reconstructed from the saved files in the loading block near the top of the script.

## Outputs

- Convergence diagnostics (R-hat) per algorithm via `summary()`
- District-level compactness (Polsby-Popper, Reock) vs. the enacted plan
- Partisan symmetry: efficiency gap, mean-median difference, declination
- Sampler efficiency: effective sample size, integrated autocorrelation time, and runtime per algorithm
- Comparison plots: compactness violins, partisan-symmetry ECDFs, and a three-panel district-rank boxplot, all referenced against the enacted plan

## Figure Generation

`Ensemble Generation Algorithm Figures.R` is an illustrative companion to the ensemble comparison above, used in the dissertation write-up to show how SMC, Flip, and Merge-split actually move through the state space. It
runs on a small 25-county subset of Florida (`fl25`, shipped with `redist`) rather than the full New Mexico map, so the underlying graph is small enough to read at a glance.

### Requirements

R packages: `redist`, `sf`, `spdep`, `igraph`.

### Running

Run the whole script top to bottom. It loads `fl25`, builds the dual graph (`shp`, `g`, `coords`), and saves the state-to-graph figure, then reuses those same objects to draw the three algorithm step figures.

### Outputs

- `state_graph_figure.png`: the state map next to its underlying adjacency graph (2 panels).
- `recom_steps.png`: Merge-split/ReCom. District graphs, spanning trees, merging two districts, sampling a new tree, cutting an edge, and splitting back into two trees (6 panels).
- `flip_steps.png`: Flip MCMC. Turning on boundary edges, gathering connected components, proposing a swap, and accepting or rejecting it (4 panels).
- `smc_steps.png`: Sequential Monte Carlo. Repeatedly cutting a spanning tree to peel off one new district at a time from the remaining map (4 panels).

## Known limitations

- **Flip MCMC does not converge** at the sample budget used here (R-hat up to ~2.2 on several statistics). Due to computational budget, `target_ess` was set relatively low.
- **Flip MCMC does not match the target distribution of SMC/Merge-split** in the above. `redist_flip` does not expose the exact `compactness` metric that `redist_smc` or `redist_mergesplit` do. Matching the compactness of the Flip algorithm can be achieved through setting `flip_constraints = redist_constr(state_map) %>% add_constr_log_st(strength = 1)` and passing this to the `constraints` argument of `redist_flip`, however, this appeared to increase runtime further.
- **Efficiency gap is coarsely discrete at n = 3 districts** (New Mexico's district count), a known property of the metric at low seat counts, not an artifact of any one sampler. Declination and mean-median difference are included as less-discretised companions. See `egap_discreteness` in the script output.
- **SMC's effective sample size** is a single particle-weight (Kish) ESS, identical across metrics by construction. It measures resampling degeneracy, not per-metric autocorrelation, so it is not directly comparable to the `coda`-based ESS reported for Flip/Merge-split. See the note printed alongside `ess_table`.
- The enacted plan's own **declination is undefined** (its vote-share pattern is itself a seat sweep), so `pct_declination_more_extreme` has no reference value and is omitted from `outlier_summary_bias`.
