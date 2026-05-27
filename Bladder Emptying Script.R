# =============================================================================
# Amoxicillin population simulation — dynamic urinary "bladder" volume + voids
# =============================================================================
#
# Re-processes the urinary (Kidney|Urine) compartment of a PK-Sim population
# simulation so that:
#   1. Urine volume V starts at 50 mL.
#   2. V grows at 1 mL / kg / h, scaled to each individual's body weight.
#   3. At user-supplied bladder-emptying times, V is reset to 50 mL and the
#      drug AMOUNT is reduced in proportion to the volume removed
#      (i.e. only the residual 50 mL worth of drug, at the prevailing
#      concentration, remains).
#   4. Final bladder drug concentration vs. time for every individual is
#      written to a long-format CSV.
#
# Important physiological assumption:
#   The original PK-Sim model has a fixed-volume Urine compartment that simply
#   accumulates renally excreted drug. We treat the time derivative of that
#   accumulated amount, dA_orig/dt, as the *input rate* of drug into the
#   "real" bladder. This holds because in PK-Sim the renal-excretion flux into
#   the Urine compartment depends on kidney processes upstream and not on what
#   is sitting in the Urine compartment itself.
#
# Requirements:
#   - R >= 4.1
#   - ospsuite (https://www.open-systems-pharmacology.org/OSPSuite-R/)
#   - data.table
#   - (optional) ggplot2 for the diagnostic plot at the end
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER CONFIGURATION — edit these paths/values to match your environment
# -----------------------------------------------------------------------------
setwd("C:/OSP R")
# Simulation file exported from PK-Sim (Population simulation).
# The simulation .pkml carries the model structure + dosing/output schemas;
# the population .csv carries the per-individual parameter values.
sim_pkml_path       <- "XXXX"
population_csv_path <- "XXXX"   # exported from PK-Sim:
# right-click population -> "Export for simulation..."

# Where to write the per-individual concentration profile
output_csv_path <- "XXXX"

# ---- Bladder dynamics --------------------------------------------------------
V_initial_mL              <- 50    # Urine volume at t = 0
V_residual_mL             <- 50    # Urine volume retained immediately after a void
urine_rate_mL_per_kg_per_h <- 1    # Urine production rate per kg body weight

# User-supplied emptying schedule (hours after dose).
# Same schedule applied to every individual. Edit freely.
# To use a per-individual schedule instead, see the "PER_INDIVIDUAL_SCHEDULE"
# block further down.
emptying_times_h <- c(3, 6, 9, 12, 15, 20)

# ---- Drug / unit conversion --------------------------------------------------
# OSP returns the Urine amount in µmol. We convert to mg using the molecular
# weight so the CSV carries both µmol/mg amounts and µmol/mL & mg/mL concentrations.
molecule_name        <- "Amoxicillin"
molecular_weight_g_per_mol <- 365.4         # amoxicillin MW = 365.4 g/mol
# 1 µmol -> MW (g/mol) µg  ->  MW/1000 mg
umol_to_mg <- molecular_weight_g_per_mol / 1000   # = 0.3654 mg per µmol

# Quantity paths in the simulation
amount_path        <- "Organism|Kidney|Urine|Amoxicillin"   # confirmed from your PKML
bodyweight_param   <- "Organism|Weight"                     # parameter path in the population CSV (kg)


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------
library(ospsuite)
library(data.table)


# -----------------------------------------------------------------------------
# 2. Load simulation, declare outputs, load population
# -----------------------------------------------------------------------------
message("Loading simulation: ", sim_pkml_path)
sim <- loadSimulation(sim_pkml_path)

# Make sure the Urine amount is in the output schema.
# (Harmless if it is already an output.)
clearOutputs(sim)
addOutputs(quantitiesOrPaths = amount_path, simulation = sim)

message("Loading population: ", population_csv_path)
population <- loadPopulation(population_csv_path)

n_individuals <- population$count
message("Population size: ", n_individuals)


# -----------------------------------------------------------------------------
# 3. Run the population simulation
# -----------------------------------------------------------------------------
message("Running population simulation (", n_individuals, " individuals)...")
sim_results <- runSimulations(simulations = sim, population = population)[[1]]


# -----------------------------------------------------------------------------
# 4. Extract results into a tidy data.table
# -----------------------------------------------------------------------------
# getOutputValues returns:
#   $data     -> data.frame with columns IndividualId, Time (seconds),
#                and one column per requested quantity path
#   $metaData -> per-quantity unit / dimension info
out <- getOutputValues(
  simulationResults = sim_results,
  quantitiesOrPaths = amount_path
)

raw <- as.data.table(out$data)
amount_unit <- out$metaData[[amount_path]]$unit       # OSP base unit for amount = µmol
time_unit   <- "h"                                    # we convert below

# Sanity check the amount unit before applying the µmol -> mg conversion
if (!is.null(amount_unit) && !(tolower(amount_unit) %in% c("µmol", "umol", "μmol"))) {
  warning(sprintf(
    "Expected amount unit µmol but OSP reports '%s'. Adjust 'umol_to_mg' if needed.",
    amount_unit))
}

# OSPSuite-R reports the Time column from getOutputValues() in MINUTES
# (the OSP simulation engine's base time unit), so convert to hours by /60.
# (Previously divided by 3600 — that's seconds, which gave a 60x truncation:
#  a 24 h simulation appeared to end at 0.4 h.)
raw[, Time_h := Time / 60]
setnames(raw, amount_path, "AmountOrig_umol")
raw <- raw[, .(IndividualId, Time_h, AmountOrig_umol)]
setkey(raw, IndividualId, Time_h)

# Defensive sanity check — flag if the resulting end time looks wrong.
sim_end_h <- max(raw$Time_h)
message(sprintf("Simulation time range: 0 to %.2f h (%d output points per individual)",
                sim_end_h, length(unique(raw$Time_h))))
if (sim_end_h < 1) {
  warning("Simulation end time is < 1 h. Check that getOutputValues() returns ",
          "Time in minutes (the OSP default). If your build returns seconds, ",
          "change `Time / 60` above to `Time / 3600`.")
}


# -----------------------------------------------------------------------------
# 5. Per-individual body weight (kg)
# -----------------------------------------------------------------------------
# In a PK-Sim population CSV, body weight is stored as the *parameter*
# "Organism|Weight" (kg) — NOT as a covariate. We pull it via the Population
# API. Pass the path *without* the "[kg]" unit suffix — the OSP API requires
# the bare path. Values are returned in SI base units (kg for weight).
#
# IMPORTANT: the returned vector is aligned to `population$allIndividualIds`,
# not to 0..(N-1). We build a lookup table keyed on IndividualId so the
# body weight is always matched to the correct simulated individual.
if (!(bodyweight_param %in% population$allParameterPaths)) {
  stop("Body weight parameter '", bodyweight_param,
       "' not found in population. Available paths include: ",
       paste(head(population$allParameterPaths, 5), collapse = ", "), " ...")
}
bw_values  <- population$getParameterValues(bodyweight_param)
bw_ind_ids <- population$allIndividualIds
stopifnot(length(bw_values) == n_individuals,
          length(bw_ind_ids) == n_individuals)
bw_lookup <- setNames(bw_values, as.character(bw_ind_ids))   # name = IndividualId
message(sprintf("Body weight range: %.1f - %.1f kg (median %.1f)",
                min(bw_values), max(bw_values), median(bw_values)))


# -----------------------------------------------------------------------------
# 6. Optional: per-individual emptying schedule
# -----------------------------------------------------------------------------
# Default: every individual uses `emptying_times_h` defined at the top.
# Stored as a named list keyed by IndividualId so each schedule pairs
# correctly with its simulated trajectory.
emptying_schedule <- setNames(
  replicate(n_individuals, emptying_times_h, simplify = FALSE),
  as.character(bw_ind_ids)
)
# ---- PER_INDIVIDUAL_SCHEDULE -------------------------------------------------
# Example override (uncomment & edit). Each list element is a numeric vector
# of voiding times (in hours) for that individual. Names MUST be IndividualIds.
# set.seed(1)
# emptying_schedule <- setNames(
#   lapply(seq_len(n_individuals), function(i) {
#     sort(cumsum(rexp(8, rate = 1/3)))   # mean inter-void interval ~3 h
#   }),
#   as.character(bw_ind_ids)
# )
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# 7. Bladder re-simulation core
# -----------------------------------------------------------------------------
# For one individual:
#   - march through native output times
#   - between two consecutive output times t_{k-1} and t_k:
#         dV = (1 mL/kg/h * BW) * (t_k - t_{k-1})
#         dA = AmountOrig(t_k) - AmountOrig(t_{k-1})    (drug input from kidney)
#     V_k = V_{k-1} + dV
#     A_k = A_{k-1} + dA
#   - if any emptying time falls within (t_{k-1}, t_k]:
#         A_k <- A_k * V_residual / V_k
#         V_k <- V_residual
#     (i.e. only the residual 50 mL worth of drug, at the prevailing
#      concentration, is retained.  Multiple events in one step are collapsed
#      to one — output spacing is fine enough that this is not a concern.)
# -----------------------------------------------------------------------------
process_individual <- function(ind_data, bw, void_times,
                               V0    = V_initial_mL,
                               V_res = V_residual_mL,
                               rate_per_kg = urine_rate_mL_per_kg_per_h,
                               mw_factor   = umol_to_mg) {
  t       <- ind_data$Time_h
  A_orig  <- ind_data$AmountOrig_umol
  n       <- length(t)
  growth  <- rate_per_kg * bw                  # mL / h
  
  V <- numeric(n); A_umol <- numeric(n)
  V[1]      <- V0
  A_umol[1] <- A_orig[1]                       # usually 0 at t=0
  
  for (k in seq.int(2L, n)) {
    dt   <- t[k] - t[k - 1L]
    V_k  <- V[k - 1L] + growth * dt
    A_k  <- A_umol[k - 1L] + (A_orig[k] - A_orig[k - 1L])
    
    # Apply any voiding event(s) that fall inside (t_{k-1}, t_k]
    if (any(void_times > t[k - 1L] & void_times <= t[k])) {
      A_k <- A_k * V_res / V_k
      V_k <- V_res
    }
    
    V[k]      <- V_k
    A_umol[k] <- A_k
  }
  
  A_mg          <- A_umol * mw_factor          # µmol  -> mg          (mw_factor = MW/1000)
  conc_umol_mL  <- A_umol / V                  # µmol  / mL urine
  conc_mg_mL    <- A_mg   / V                  # mg    / mL urine  ==  g/L  ==  mg/mL
  
  data.table(
    Time_h                 = t,
    BodyWeight_kg          = bw,
    UrineVolume_mL         = V,
    AmountInBladder_umol   = A_umol,
    AmountInBladder_mg     = A_mg,
    ConcInBladder_umol_per_mL = conc_umol_mL,
    ConcInBladder_mg_per_mL   = conc_mg_mL,
    ConcInBladder_mg_per_L    = conc_mg_mL * 1000   # = µg/mL, common clinical unit
  )
}


# -----------------------------------------------------------------------------
# 8. Apply to all individuals
# -----------------------------------------------------------------------------
message("Re-processing bladder dynamics for ", n_individuals, " individuals...")

# IDs in the simulation results — should match the population's IndividualIds.
ind_ids <- sort(unique(raw$IndividualId))
missing_ids <- setdiff(as.character(ind_ids), names(bw_lookup))
if (length(missing_ids) > 0) {
  stop("Simulation contains IndividualIds with no body-weight in the ",
       "population (first few): ",
       paste(head(missing_ids, 5), collapse = ", "))
}

result_list <- vector("list", length(ind_ids))

for (i in seq_along(ind_ids)) {
  id_k     <- ind_ids[i]
  key      <- as.character(id_k)
  ind_data <- raw[IndividualId == id_k]
  bw_k     <- bw_lookup[[key]]
  voids_k  <- emptying_schedule[[key]]
  
  res_k <- process_individual(ind_data, bw = bw_k, void_times = voids_k)
  res_k[, IndividualId := id_k]
  result_list[[i]] <- res_k
}

result_dt <- rbindlist(result_list)
setcolorder(result_dt, c("IndividualId", "Time_h", "BodyWeight_kg",
                         "UrineVolume_mL",
                         "AmountInBladder_umol", "AmountInBladder_mg",
                         "ConcInBladder_umol_per_mL",
                         "ConcInBladder_mg_per_mL",
                         "ConcInBladder_mg_per_L"))


# -----------------------------------------------------------------------------
# 9. Write CSV
# -----------------------------------------------------------------------------
fwrite(result_dt, output_csv_path)
message("Wrote: ", normalizePath(output_csv_path, mustWork = FALSE),
        "  (", format(nrow(result_dt), big.mark = ","), " rows)")


# -----------------------------------------------------------------------------
# 10. (Optional) Quick diagnostic plot — first 25 individuals
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  plot_ids <- head(ind_ids, 25)
  
  p_conc <- ggplot(result_dt[IndividualId %in% plot_ids],
                   aes(Time_h, ConcInBladder_mg_per_L, group = IndividualId)) +
    geom_line(alpha = 0.6) +
    labs(title = "Amoxicillin bladder concentration (first 25 individuals)",
         x = "Time (h)", y = "Concentration (mg/L = µg/mL)") +
    theme_bw()
  p_vol  <- ggplot(result_dt[IndividualId %in% plot_ids],
                   aes(Time_h, UrineVolume_mL, group = IndividualId)) +
    geom_line(alpha = 0.6) +
    labs(title = "Urine volume (first 25 individuals)",
         x = "Time (h)", y = "Urine volume (mL)") +
    theme_bw()
  
  ggsave("bladder_concentration_preview.png", p_conc, width = 7, height = 4, dpi = 150)
  ggsave("urine_volume_preview.png",          p_vol,  width = 7, height = 4, dpi = 150)
  message("Saved diagnostic plots: bladder_concentration_preview.png, urine_volume_preview.png")
}

# =============================================================================
# End of script
# =============================================================================

# =============================================================================
# Amoxicillin urinary "bladder" — population summary statistics
# =============================================================================
#
# Companion script to `amoxicillin_bladder_postprocess.R`.
#
# Reads the per-individual long-format CSV produced by the main script and
# computes, at each native output time point:
#     - median
#     - 2.5th percentile
#     - 97.5th percentile
# for every concentration column (µmol/mL, mg/mL, mg/L) and the urine volume.
#
# Output: a single wide-format CSV with one row per time point, e.g.
#   Time_h, n,
#   UrineVolume_mL_median, UrineVolume_mL_p2.5, UrineVolume_mL_p97.5,
#   ConcInBladder_mg_per_L_median, ConcInBladder_mg_per_L_p2.5, ConcInBladder_mg_per_L_p97.5,
#   ...
#
# An optional diagnostic plot (median + 95% population interval ribbon) is
# saved if ggplot2 is available.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER CONFIGURATION
# -----------------------------------------------------------------------------
input_csv_path   <- "amoxicillin_bladder_concentration.csv"   # from main script
output_csv_path  <- "amoxicillin_bladder_summary.csv"

# Which percentiles to compute (low, high). Median is always added.
percentile_low   <- 0.025
percentile_high  <- 0.975

# Columns to summarise. Defaults to all numeric value columns produced by the
# main script — edit if you only want a subset.
summary_columns  <- c(
  "UrineVolume_mL",
  "AmountInBladder_umol",
  "AmountInBladder_mg",
  "ConcInBladder_umol_per_mL",
  "ConcInBladder_mg_per_mL",
  "ConcInBladder_mg_per_L"
)


# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------
library(data.table)


# -----------------------------------------------------------------------------
# 2. Load the per-individual CSV
# -----------------------------------------------------------------------------
if (!file.exists(input_csv_path)) {
  stop("Input file not found: ", normalizePath(input_csv_path, mustWork = FALSE),
       "\nRun amoxicillin_bladder_postprocess.R first.")
}

dt <- fread(input_csv_path)
message(sprintf("Loaded %s: %s rows, %d individuals, %d time points",
                input_csv_path,
                format(nrow(dt), big.mark = ","),
                length(unique(dt$IndividualId)),
                length(unique(dt$Time_h))))

# Sanity-check that requested summary columns exist
missing_cols <- setdiff(summary_columns, names(dt))
if (length(missing_cols) > 0) {
  stop("Summary columns not found in input CSV: ",
       paste(missing_cols, collapse = ", "))
}


# -----------------------------------------------------------------------------
# 3. Per-time-point summary
# -----------------------------------------------------------------------------
# Strategy: melt to long form, aggregate once, then cast back to wide.
# This is far easier to reason about than nested data.table list expressions.
#
#   step 1: keep only IndividualId, Time_h, and the value columns
#   step 2: melt -> columns (IndividualId, Time_h, variable, value)
#   step 3: group by (Time_h, variable), compute median / p_low / p_high
#   step 4: melt the three stats to long, then dcast to wide:
#           one row per Time_h, columns "<variable>_<stat>"

# Friendly suffixes for column names (e.g. "p2.5", "p97.5")
suffix_low  <- sprintf("p%g", percentile_low  * 100)
suffix_high <- sprintf("p%g", percentile_high * 100)

long_dt <- melt(
  dt[, c("IndividualId", "Time_h", summary_columns), with = FALSE],
  id.vars       = c("IndividualId", "Time_h"),
  measure.vars  = summary_columns,
  variable.name = "variable",
  value.name    = "value"
)

# n is the same for every variable at a given Time_h, so grab it once
n_dt <- dt[, .(n = .N), by = Time_h]

agg_dt <- long_dt[, .(
  median = stats::median(value, na.rm = TRUE),
  p_low  = stats::quantile(value, probs = percentile_low,  na.rm = TRUE, names = FALSE),
  p_high = stats::quantile(value, probs = percentile_high, na.rm = TRUE, names = FALSE)
),
by = .(Time_h, variable)]

# Long -> wide: one column per (variable, stat)
agg_long <- melt(agg_dt,
                 id.vars       = c("Time_h", "variable"),
                 measure.vars  = c("median", "p_low", "p_high"),
                 variable.name = "stat",
                 value.name    = "value")

# Map internal stat name -> output suffix
stat_suffix_map <- c(median = "median", p_low = suffix_low, p_high = suffix_high)
agg_long[, colname := paste0(variable, "_", stat_suffix_map[as.character(stat)])]

summary_dt <- dcast(agg_long, Time_h ~ colname, value.var = "value")

# Add per-time-point sample size and order columns/rows
summary_dt <- merge(n_dt, summary_dt, by = "Time_h", all = TRUE)

# Reorder value columns so each variable's three stats stay together,
# in the requested column order (median, p_low, p_high)
ordered_cols <- c("Time_h", "n",
                  unlist(lapply(summary_columns, function(col) {
                    paste0(col, "_", c("median", suffix_low, suffix_high))
                  })))
ordered_cols <- intersect(ordered_cols, names(summary_dt))
setcolorder(summary_dt, ordered_cols)
setorder(summary_dt, Time_h)


# -----------------------------------------------------------------------------
# 4. Write the summary CSV
# -----------------------------------------------------------------------------
fwrite(summary_dt, output_csv_path)
message("Wrote: ", normalizePath(output_csv_path, mustWork = FALSE),
        "  (", nrow(summary_dt), " time points, ",
        ncol(summary_dt) - 2L, " summary columns)")


# -----------------------------------------------------------------------------
# 5. Optional diagnostic plot — median + 95% population interval ribbon
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  median_col <- paste0("ConcInBladder_mg_per_L_", "median")
  low_col    <- paste0("ConcInBladder_mg_per_L_", suffix_low)
  high_col   <- paste0("ConcInBladder_mg_per_L_", suffix_high)
  
  if (all(c(median_col, low_col, high_col) %in% names(summary_dt))) {
    p <- ggplot(summary_dt, aes(x = Time_h)) +
      geom_ribbon(aes(ymin = .data[[low_col]], ymax = .data[[high_col]]),
                  fill = "steelblue", alpha = 0.25) +
      geom_line(aes(y = .data[[median_col]]),
                colour = "steelblue", linewidth = 0.8) +
      labs(title = "Amoxicillin bladder concentration — population summary",
           subtitle = sprintf("Median + %g–%g%% interval (n = %d individuals)",
                              percentile_low * 100, percentile_high * 100,
                              length(unique(dt$IndividualId))),
           x = "Time (h)",
           y = "Concentration (mg/L = µg/mL)") +
      theme_bw()
    
    ggsave("bladder_concentration_summary.png", p,
           width = 7, height = 4, dpi = 150)
    message("Saved diagnostic plot: bladder_concentration_summary.png")
  }
}

# =============================================================================
# End of script
# =============================================================================