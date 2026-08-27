#========================================================

# Redistricting Ensemble Generation Algorithm Comparison
# SMC vs Flip vs Merge-Split

#========================================================

#Required packages

library(redist)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(alarmdata)

set.seed(2)

#==============
# Configuration
#==============

state_abbr = "NM" #Abbreviation of US state
year = 2020

pop_tol = 0.05 # Population tolerance, typically seen as delta

nsims_smc = 5000
nruns_smc = 8

target_ess = 5000 # Target effective sample size
warmup_fraction = 0.2 # Burn-in fraction of total sims

adapt_k_thresh_ms = 0.99 # Adaptive k threshold for Merge-split algorithm
out_dir = "." # Points saveRDS calls to working directory

n_cores = 8 # n_cores set to physical cores of my laptop
state_name = ifelse(state_abbr %in% state.abb,
                    state.name[match(state_abbr, state.abb)],
                    state_abbr)

# Load state data

state_map = alarm_50state_map(state_abbr, year = year) # Retrieve state map from alarmdata package
state_map = set_pop_tol(state_map, pop_tol) # Set population tolerance

ndists = attr(state_map, "ndists") # Retrieve number of districts
target_pop = get_target(state_map) # Retrive target population per district
existing_col = attr(state_map, "existing_col") # Retrieve enacted plan column

#===============================================================
# Load previously saved run outputs + reconstruct run parameters
#===============================================================

# Skip straight to "Run SMC ensemble" below to regenerate
# everything from scratch instead.

#Load saved ensembles
plans_smc = readRDS("plans_smc_NM_2020.rds")
plans_flip = readRDS("plans_flip_NM_2020.rds")
plans_ms = readRDS("plans_ms_NM_2020.rds")

#Load wall-clock time
time_smc = readRDS("time_smc.rds")
time_flip = readRDS("time_flip.rds")
time_ms = readRDS("time_ms.rds")

#Load calibration metadata
flip_calib_meta = readRDS("flip_calib_meta.rds")
ms_calib_meta = readRDS("ms_calib_meta.rds")
smc_calib_meta = readRDS("smc_calib_meta.rds")

# Chain counts
nchains_flip = 4
nchains_ms = 8

# Calibration diagnostics pulled from metadata save files
flip_iat = flip_calib_meta$iat
flip_thin = flip_calib_meta$thin
warmup_flip = flip_calib_meta$warmup
flip_calib_time = flip_calib_meta$calib_time
flip_calib_acceptance = flip_calib_meta$acceptance

ms_iat = ms_calib_meta$iat
ms_thin = ms_calib_meta$thin
warmup_ms = ms_calib_meta$warmup
ms_calib_time = ms_calib_meta$calib_time

smc_calib_kish_ess = smc_calib_meta$kish_ess
nsims_smc_calib = smc_calib_meta$nsims_calib
smc_calib_time = smc_calib_meta$calib_time

#Retained draw counts
n_draws_flip = plans_flip %>% filter(draw != existing_col) %>% distinct(chain, draw) %>% nrow()
n_draws_ms = plans_ms %>% filter(draw != existing_col) %>% distinct(chain, draw) %>% nrow()
n_draws_smc = plans_smc %>% filter(draw != existing_col) %>% distinct(chain, draw) %>% nrow()

nsims_flip = n_draws_flip / nchains_flip
nsims_ms = n_draws_ms / nchains_ms
nsims_smc = n_draws_smc / nruns_smc

# If working from saved outputs, skip down to District-Level Metrics section to resume

#===========================
# SMC: Calibrate, Then Run
#===========================

# Pilot run at a small particle count, used to estimate how Kish ESS scales
# with nsims and allocate nsims for full run
nsims_smc_calib = 2000

time_smc_calib = system.time(
  plans_smc_calib <- redist_smc(
    state_map,
    nsims = nsims_smc_calib,
    runs = nruns_smc,
    ncores = n_cores
  )
)

smc_calib_draw_chain = plans_smc_calib %>% distinct(draw, chain) %>% pull(chain)
smc_calib_weights = get_plans_weights(plans_smc_calib)
smc_calib_kish_ess = sum(tapply(smc_calib_weights, smc_calib_draw_chain, function(w) sum(w)^2 / sum(w^2)))

smc_calib_sec_per_particle = time_smc_calib[["elapsed"]] / nsims_smc_calib

nsims_smc = ceiling(target_ess / smc_calib_kish_ess * nsims_smc_calib)
smc_hours_needed = (nsims_smc * smc_calib_sec_per_particle) / 3600

saveRDS(list(kish_ess = smc_calib_kish_ess, nsims_calib = nsims_smc_calib, nsims = nsims_smc,
             calib_time = time_smc_calib),
        file.path(out_dir, "smc_calib_meta.rds"))

cat("SMC calibration Kish ESS:", round(smc_calib_kish_ess, 1),
    "(at nsims =", nsims_smc_calib, ")\n")
cat("SMC: nsims =", nsims_smc, "| est. hours =", round(smc_hours_needed, 1), "\n")

time_smc = system.time(
  plans_smc <- redist_smc(
    state_map,
    nsims = nsims_smc,
    runs = nruns_smc,
    ncores = n_cores
  )
)

saveRDS(plans_smc, file.path(out_dir, paste0("plans_smc_", state_abbr, "_", year, ".rds")))
saveRDS(time_smc, file.path(out_dir, "time_smc.rds"))

# Integrated autocorrelation time (IAT) per chain for compactness and dem_share, from an
# unthinned calibration run (works for both samplers, since with thin=1, warmup=0, nsims is
# the raw step count in both)

estimate_iat = function(plans_calib, nsims_calib, map) {
  metrics = plans_calib %>%
    mutate(compactness = comp_polsby(pl(), map), dem_share = group_frac(map, ndv, ndv + nrv)) %>%
    group_by(chain, draw) %>%
    summarise(compactness = mean(compactness), dem_share = mean(dem_share), .groups = "drop")
  
  chains = unique(metrics$chain)
  ess = sapply(c("compactness", "dem_share"), function(m) {
    chain_list = lapply(chains, function(ch) coda::as.mcmc(metrics[[m]][metrics$chain == ch]))
    unname(coda::effectiveSize(coda::mcmc.list(chain_list)))
  })
  
  nsims_calib / ess # IAT per metric
}

#===========================
# Flip: Calibration and Run
#===========================

nchains_flip = 4
nsims_flip_calib = 2000

time_flip_calib = system.time(
  plans_flip_calib <- redist_flip(
    state_map,
    nsims = nsims_flip_calib,
    warmup = 0,
    chains = nchains_flip,
    ncores = nchains_flip,
    init_plan = flip_init_plan,
    thin = 1,
    eprob = 0.01,
    lambda = 3,
    verbose = FALSE
  )
)

get_mh_acceptance_rate(plans_flip_calib) #quick check for 20-40% acceptance across chains

flip_iat = estimate_iat(plans_flip_calib, nsims_flip_calib, state_map)
flip_thin = ceiling(max(flip_iat)) # thin by the slowest-mixing metric

calib_sec_per_step = time_flip_calib[["elapsed"]] / nsims_flip_calib

# Raw steps per chain needed so the combined, post-thin ESS reaches target_ess,
# given the slowest-mixing metric's IAT from calibration

flip_steps_needed = (target_ess * max(flip_iat)) / nchains_flip
flip_hours_needed = (flip_steps_needed * calib_sec_per_step) / 3600 # Rough lower bound
warmup_flip = round(warmup_fraction * flip_steps_needed)
nsims_flip = ceiling((flip_steps_needed - warmup_flip) / flip_thin)

saveRDS(list(iat = flip_iat, thin = flip_thin, warmup = warmup_flip, nsims = nsims_flip,
             calib_time = time_flip_calib, acceptance = get_mh_acceptance_rate(plans_flip_calib)),
        file.path(out_dir, "flip_calib_meta.rds"))

cat("Flip calibration IAT:", paste(names(flip_iat), round(flip_iat, 1), collapse = ", "), "\n")
cat("Flip: thin =", flip_thin, "| warmup =", warmup_flip, "| nsims =", nsims_flip,
    "| raw steps =", warmup_flip + nsims_flip * flip_thin,
    "| est. hours =", round(flip_hours_needed, 1), "\n")

time_flip = system.time(
  plans_flip <- redist_flip(
    state_map,
    nsims = nsims_flip,
    warmup = warmup_flip,
    chains = nchains_flip,
    ncores = nchains_flip,
    init_plan = flip_init_plan,
    constraints = flip_constraints,
    thin = flip_thin,
    eprob = 0.01,
    lambda = 3,
    verbose = TRUE
  )
)

saveRDS(plans_flip, file.path(out_dir, paste0("plans_flip_", state_abbr, "_", year, ".rds")))
saveRDS(time_flip, file.path(out_dir, "time_flip.rds"))

#=================================
# Merge-Split: Calibration and Run
#=================================

nchains_ms = n_cores
nsims_ms_calib = 2000 # unthinned pilot, same length as flip's for a fair comparison

time_ms_calib = system.time(
  plans_ms_calib <- redist_mergesplit_parallel(
    state_map,
    nsims = nsims_ms_calib,
    chains = nchains_ms,
    init_plan = flip_init_plan,
    init_name = FALSE,
    compactness = 1,
    adapt_k_thresh = adapt_k_thresh_ms,
    warmup = 0,
    thin = 1,
    verbose = FALSE
  )
)

ms_iat = estimate_iat(plans_ms_calib, nsims_ms_calib, state_map)
ms_thin = ceiling(max(ms_iat)) # thin by the slowest-mixing metric

ms_calib_sec_per_step = time_ms_calib[["elapsed"]] / nsims_ms_calib

ms_steps_needed = (target_ess * max(ms_iat)) / nchains_ms
ms_hours_needed = (ms_steps_needed * ms_calib_sec_per_step) / 3600
warmup_ms = round(warmup_fraction * ms_steps_needed)
nsims_ms = ceiling(ms_steps_needed / ms_thin) * ms_thin

saveRDS(list(iat = ms_iat, thin = ms_thin, warmup = warmup_ms, nsims = ms_steps_needed, 
             calib_time = time_ms_calib), file.path(out_dir, "ms_calib_meta.rds"))

cat("Merge-split calibration IAT:", paste(names(ms_iat), round(ms_iat, 1), collapse = ", "), "\n")
cat("Merge-split: thin =", ms_thin, "| warmup =", warmup_ms, "| nsims =", nsims_ms,
    "| est. hours =", round(ms_hours_needed, 1), "\n")

time_ms = system.time(
  plans_ms <- redist_mergesplit_parallel(
    state_map,
    nsims = nsims_ms,
    chains = nchains_ms,
    init_plan = flip_init_plan,
    init_name = FALSE,
    compactness = 1,
    adapt_k_thresh = adapt_k_thresh_ms,
    warmup = warmup_ms,
    thin = ms_thin,
    verbose = TRUE
  )
)

saveRDS(plans_ms, file.path(out_dir, paste0("plans_ms_", state_abbr, "_", year, ".rds")))
saveRDS(time_ms, file.path(out_dir, "time_ms.rds"))

#========================
# District-Level Metrics
#========================

# Extract the enacted plan's district assignment

enacted_vec = as.integer(get_plans_matrix(plans_smc)[, existing_col])

# Flip's saved output already contains four auto-registered reference draws
# from the original run. These aren't meaningful reference plans, so drop them here.

plans_flip = plans_flip %>% filter(!grepl("^<init", as.character(draw)))

plans_flip = add_reference(plans_flip, enacted_vec, name = existing_col)
plans_ms = add_reference(plans_ms, enacted_vec, name = existing_col)

# Polsby-Popper and Reock compactness, population deviation, Democratic vote share
add_district_metrics = function(plans, map) {
  plans %>%
    mutate(
      compactness = comp_polsby(pl(), map),
      compactness_reock = comp_reock(pl(), map),
      pop_dev = plan_parity(map),
      dem_share = group_frac(map, ndv, ndv + nrv)
    )
}

plans_smc = add_district_metrics(plans_smc, state_map)
plans_flip = add_district_metrics(plans_flip, state_map)
plans_ms = add_district_metrics(plans_ms, state_map)

# Aligns district numbering across ensembles to the enacted plan's, since raw district
# labels are otherwise arbitrary and not comparable across samplers
plans_ms = match_numbers(plans_ms, enacted_vec) %>% mutate(district = as.integer(as.character(district)))
plans_flip = match_numbers(plans_flip, enacted_vec) %>% mutate(district = as.integer(as.character(district)))

#===================
# Convergence Check
#===================

# redist's built-in summary() reports Rhat per sampler

summary(plans_smc)
summary(plans_flip)
summary(plans_ms)

#=======================
# Partisan Bias Metrics
#=======================

# Efficiency gap, mean-median, declination: the three most consistently
# literature-supported symmetry measures

compute_bias_metrics = function(plans, map, algorithm_name) {
  plans %>%
    mutate(
      efficiency_gap = part_egap(plans = plans, shp = map, dvote = ndv, rvote = nrv),
      mean_median = part_mean_median(plans = plans, shp = map, dvote = ndv, rvote = nrv),
      declination = part_decl(plans = plans, shp = map, dvote = ndv, rvote = nrv)
    ) %>%
    group_by(chain, draw) %>%
    mutate(dem_seats = sum(dem_share > 0.5)) %>%
    ungroup() %>%
    mutate(sweep = dem_seats %in% c(0, ndists)) %>% 
    distinct(chain, draw, efficiency_gap, mean_median, declination, sweep) %>%
    `class<-`(c("tbl_df", "tbl", "data.frame")) %>%
    mutate(algorithm = algorithm_name)
}

bias_smc = compute_bias_metrics(plans_smc, state_map, "SMC")
bias_flip = compute_bias_metrics(plans_flip, state_map, "Flip")
bias_ms = compute_bias_metrics(plans_ms, state_map, "Merge-split")

bias_all = bind_rows(bias_smc, bias_flip, bias_ms)
enacted_bias = bias_smc %>% filter(draw == existing_col)

#===================================
# Compare Enacted Plan to Ensembles
#===================================

enacted = plans_smc %>% filter(draw == existing_col)

sim_smc = plans_smc %>%
  filter(draw != existing_col) %>%
  `class<-`(c("tbl_df", "tbl", "data.frame")) %>%
  mutate(algorithm = "SMC")

sim_flip = plans_flip %>%
  filter(draw != existing_col) %>%
  `class<-`(c("tbl_df", "tbl", "data.frame")) %>%
  mutate(algorithm = "Flip")

sim_ms = plans_ms %>%
  filter(draw != existing_col) %>%
  `class<-`(c("tbl_df", "tbl", "data.frame")) %>%
  mutate(algorithm = "Merge-split")

sim_all = bind_rows(sim_smc, sim_flip, sim_ms)

#===========================================
# Sampler Efficiency: ESS, IAT, and Runtime
#===========================================

# Per-metric effective sample size, all three samplers

flip_draw_metrics = sim_flip %>%
  group_by(chain, draw) %>%
  summarise(compactness = mean(compactness), compactness_reock = mean(compactness_reock),
            dem_share = mean(dem_share), .groups = "drop") %>%
  inner_join(bias_flip, by = c("chain", "draw"))

ms_draw_metrics = sim_ms %>%
  group_by(chain, draw) %>%
  summarise(compactness = mean(compactness), compactness_reock = mean(compactness_reock),
            dem_share = mean(dem_share), .groups = "drop") %>%
  inner_join(bias_ms, by = c("chain", "draw"))

# Chain-wise ESS via coda, pooled across chains as an mcmc.list
compute_ess_by_metric = function(df, metric_cols) {
  chains = unique(df$chain)
  sapply(metric_cols, function(m) {
    chain_list = lapply(chains, function(ch) coda::as.mcmc(df[[m]][df$chain == ch]))
    unname(coda::effectiveSize(coda::mcmc.list(chain_list)))
  })
}

# declination excluded here as it is undefined (NA) whenever a chain sweeps all seats,
# which coda::effectiveSize() can't handle. Particularly applicable for low district number states

metric_cols_ess = c("compactness", "compactness_reock", "dem_share", "efficiency_gap", "mean_median")

flip_ess = compute_ess_by_metric(flip_draw_metrics, metric_cols_ess)
ms_ess = compute_ess_by_metric(ms_draw_metrics, metric_cols_ess)

# SMC design effective sample size (Kish ESS from importance weights). This measures
# particle weight degeneracy, not per-metric autocorrelation, so it's the same value
# for every metric by construction and is not directly comparable to the coda-based
# ESS values above

draw_chain = plans_smc %>% distinct(draw, chain) %>% pull(chain)
smc_weights = get_plans_weights(plans_smc)
smc_ess_kish = sum(tapply(smc_weights, draw_chain, function(w) sum(w)^2 / sum(w^2)))

smc_ess = setNames(rep(smc_ess_kish, length(metric_cols_ess)), metric_cols_ess)

ess_table = bind_rows(
  as_tibble(t(smc_ess)) %>% mutate(algorithm = "SMC"),
  as_tibble(t(flip_ess)) %>% mutate(algorithm = "Flip"),
  as_tibble(t(ms_ess)) %>% mutate(algorithm = "Merge-split")
) %>%
  relocate(algorithm)

ess_table

n_draws_table = data.frame(
  algorithm = c("SMC", "Flip", "Merge-split"),
  n_draws = c(NA, n_draws_flip, n_draws_ms) # NA for SMC as draws are independent by construction, so IAT is undefined
)

iat_table = ess_table %>%
  left_join(n_draws_table, by = "algorithm") %>%
  mutate(across(all_of(metric_cols_ess), ~ n_draws / .x)) %>%
  select(-n_draws)

iat_table

runtime_table = data.frame(
  algorithm = c("SMC", "Flip", "Merge-split"),
  runtime_sec = c(time_smc[["elapsed"]], time_flip[["elapsed"]], time_ms[["elapsed"]])
)

runtime_table

ess_per_sec_table = ess_table %>%
  left_join(runtime_table, by = "algorithm") %>%
  mutate(across(all_of(metric_cols_ess), ~ .x / runtime_sec)) %>%
  select(-runtime_sec)

ess_per_sec_table

# Display-only relabel, applied last so it doesn't break the joins above

ess_table_display = ess_table %>%
  mutate(algorithm = ifelse(algorithm == "SMC", "SMC (particle/Kish ESS)", algorithm))
ess_table_display

#================================
# Outlier / Extremity Summaries
#================================

outlier_summary_district = sim_all %>%
  group_by(algorithm, draw) %>%
  summarise(mean_compactness = mean(compactness), mean_compactness_reock = mean(compactness_reock), .groups = "drop") %>%
  group_by(algorithm) %>%
  summarise(
    pct_more_compact_than_enacted = mean(mean_compactness > mean(enacted$compactness)),
    pct_more_compact_reock_than_enacted = mean(mean_compactness_reock > mean(enacted$compactness_reock))
  )

outlier_summary_district

# Share of ensemble draws more extreme than the enacted plan, by symmetry metric

outlier_summary_bias = bias_all %>%
  filter(draw != existing_col) %>%
  group_by(algorithm) %>%
  summarise(
    pct_sweep = mean(sweep),
    pct_egap_more_extreme = mean(abs(efficiency_gap) > abs(enacted_bias$efficiency_gap)),
    pct_mean_median_more_extreme = mean(abs(mean_median) > abs(enacted_bias$mean_median)),
    pct_declination_more_extreme = mean(abs(declination) > abs(enacted_bias$declination), na.rm = TRUE),
    pct_declination_na = mean(is.na(declination)) # driven by the sweep rate, since declination is undefined under a sweep
  )

outlier_summary_bias

# Efficiency gap is a coarsely discrete statistic at low district counts (n = 3 here),
# so this counts how many distinct values each ensemble actually realizes

egap_discreteness = bias_all %>%
  filter(draw != existing_col) %>%
  group_by(algorithm) %>%
  summarise(n_distinct_egap = n_distinct(round(efficiency_gap, 4)))

egap_discreteness

#===================
# Comparison Plots
#===================

# Violin and boxplot with enacted plan as reference line, compactness, then the three bias metrics

ggplot(sim_all, aes(x = algorithm, y = compactness)) +
  geom_violin(aes(fill = algorithm), alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width = 0.15) +
  geom_hline(yintercept = mean(enacted$compactness), color = "red", linetype = "dashed") +
  labs(
    title = paste("Compactness (Polsby-Popper) by algorithm"),
    x = NULL, y = "Polsby-Popper compactness"
  ) +
  theme_minimal()

ggplot(sim_all, aes(x = algorithm, y = compactness_reock)) +
  geom_violin(aes(fill = algorithm), alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width = 0.15) +
  geom_hline(yintercept = mean(enacted$compactness_reock), color = "red", linetype = "dashed") +
  labs(
    title = paste("Compactness (Reock) by algorithm"),
    x = NULL, y = "Reock compactness"
  ) +
  theme_minimal()

plot_bias_ecdf = function(metric_col, metric_label) {
  ggplot(bias_all %>% filter(draw != existing_col),
         aes(x = .data[[metric_col]], color = algorithm)) +
    stat_ecdf(geom = "step", linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = enacted_bias[[metric_col]], color = "red", linetype = "dashed") +
    labs(
      title = paste0(metric_label, " ECDF by algorithm"),
      subtitle = "Dashed line = enacted plan's value, y-axis = percentile within each ensemble",
      x = metric_label, y = "Cumulative probability", color = NULL
    ) +
    theme_minimal()
}

plot_bias_ecdf("efficiency_gap", "Efficiency gap")
plot_bias_ecdf("mean_median", "Mean-median difference")
plot_bias_ecdf("declination", "Declination")

# Native redist diagnostic: boxplot of a quantity across districts, sorted,
# with the enacted plan overlaid automatically.

# Shared y-axis range so panel heights are visually comparable across algorithms
dem_share_range = range(c(sim_all$dem_share, enacted$dem_share), na.rm = TRUE)

# Display-only copies with the enacted plan's draw label renamed for a plain-English
# legend, so nothing downstream that filters or joins on existing_col's literal
# value (sim_smc, enacted, bias_all, etc.) is affected
rename_enacted_level = function(plans, new_label) {
  levels(plans$draw)[levels(plans$draw) == existing_col] = new_label
  plans
}

plans_smc_display = rename_enacted_level(plans_smc, "Enacted Plan")
plans_flip_display = rename_enacted_level(plans_flip, "Enacted Plan")
plans_ms_display = rename_enacted_level(plans_ms, "Enacted Plan")

p_smc = redist.plot.distr_qtys(plans_smc_display, dem_share, sort = "asc",
                               geom = ggplot2::geom_boxplot) +
  labs(title = "SMC")
p_flip = redist.plot.distr_qtys(plans_flip_display, dem_share, sort = "asc",
                                geom = ggplot2::geom_boxplot) +
  labs(title = "Flip")
p_ms = redist.plot.distr_qtys(plans_ms_display, dem_share, sort = "asc",
                              geom = ggplot2::geom_boxplot) +
  labs(title = "Merge-split")

((p_smc | p_flip | p_ms) & coord_cartesian(ylim = dem_share_range) & labs(y = "Democratic vote share")) +
  plot_annotation(
    title = paste("Democratic vote share by district")
  )

#Flip trace plot of all 4 chains, confirming lack of convergence

redist.plot.trace(plans_flip, dem_share) +
  labs(title = "Flip MCMC trace: Democratic vote share by chain")