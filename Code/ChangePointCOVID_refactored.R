################################################################################
## Wastewater COVID-19 Surveillance: Change Point and VAR Network Analysis
##
## This script implements the piecewise VAR Granger causality framework applied
## to Oregon wastewater COVID-19 concentration data collected since Sept 2020.
##
## Pipeline:
##   1. Load and preprocess weekly wastewater data
##   2. Impute missing values (non-participation gaps) via Kalman filter
##   3. Detect structural change points per time series (Bai & Perron via strucchange)
##   4. Align segment boundaries across locations and impute boundary padding
##   5. Assess stationarity (ADF test); apply first differencing if needed
##   6. Fit lasso-penalized VAR model per segment (NGC package)
##   7. Construct Granger causality networks per segment
##   8. Compute node importance measures and visualise on geographic layout
##
## NOTE on differencing: first-differencing is applied uniformly prior to VAR
## fitting. This means Granger causality results reflect predictive relationships
## between *changes* in wastewater concentrations, not levels. To model level
## relationships, a VECM approach with prior cointegration testing would be needed.
################################################################################


# ------------------------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------------------------
library(mcp)
library(changepoint)
library(RColorBrewer)
library(glmnet)
#library(NGC)
library(simts)
library(tidyverse)
library(strucchange)   # Bai & Perron breakpoint detection
library(readxl)
library(imputeTS)      # na_kalman for Kalman filter imputation
library(igraph)
library(influential)   # IVI and other node importance measures
library(tseries)       # adf.test for stationarity
#library(MTS)           # mq() multivariate Portmanteau test for residual autocorrelation
library(fs)
options(warn = -1)
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
args_1 <- 1#args[1]

# NGC package is called via NGC:: namespace below; ensure it is installed.
# install.packages("NGC")  # uncomment if needed

## Setting seed for reproducibility
#set.seed(42)

# ------------------------------------------------------------------------------
# 1. Global parameters
# ------------------------------------------------------------------------------

N_BRK    <- 7      # Number of change points P (selected based on RSS criterion)
LAG_D    <- 2      # VAR lag order p (COVID spread ~2-3 weeks; Jiang et al. 2023)
PCT_THRS <- 40     # Missingness threshold (%): locations above this are excluded
TOP_N    <- 3      # Number of top/bottom nodes to label in network plots

# File paths — update these to match your local setup
DATA_PATH <- "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/Combined_aggregated_OHA_normalized_2026-04-24.xlsx"
GRAPH_RDS <- "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/GraphData.rds"
SAVE_DIR  <- "C:/Users/gauph/Box/FinalTestingParamComboFiles/FinalPlots_TS/"


# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

#' Fit lasso-penalized VAR to each aligned segment and return igraph networks
#'
#' @param df_ip   List of length (N_BRK+1), each element an N x l_j matrix
#' @param d       VAR lag order
#' @param n_brk   Number of change points (P)
#' @param ret_fit Logical; if TRUE also return raw NGC fit objects
#' @return List of (graph_list, fit_list)

fit_lasso_var <- function(df_ip, d, n_brk, ret_fit = FALSE) {
  grphs <- list()
  fits  <- list()

  for (m in 1:(n_brk + 1)) {
    #message(sprintf("Fitting VAR for segment %d of %d", m, n_brk + 1))

    seg        <- df_ip[[m]]
    city_names <- colnames(seg)
    seg_t      <- t(seg)  # NGC expects series as rows

    if (any(is.na(seg_t))) {
      warning(sprintf("Segment %d contains missing values after imputation.", m))
    }

    fit <- NGC::ngc(seg_t, d)
    # How many non-zero coefficients?
    print(sum(abs(fit$estMat) > 1e-8))
    
    # Range of values
    print(range(fit$estMat))

    # Attach location names to graph vertices
    igraph::V(fit$ring)$names <- city_names
    grphs[[m]] <- fit$ring

    if (ret_fit) {
      fits[[m]] <- fit
    }
  }

  return(list(graphs = grphs, fits = fits))
}


#' Compute residuals from a fitted NGC lasso VAR model
#'
#' The NGC fit object contains a coefficient matrix B of dimension N x (N*d),
#' where N is the number of series and d is the lag order. This function
#' reconstructs the fitted values and returns the residual matrix.
#'
#' @param fit  A single NGC fit object (one element of var_fits)
#' @param seg  The data matrix used to fit this segment (T x N, rows = time)
#' @param d    VAR lag order
#' @return     Matrix of residuals, dimension (T - d) x N

compute_ngc_residuals <- function(fit, seg, d) {
  T_obs <- nrow(seg)
  N     <- ncol(seg)

  # NGC stores coefficients in fit$B: an N x (N*d) matrix
  # Row i of B gives the coefficients for the equation of series i
  B <- fit$estMat

  # Build the lagged regressor matrix Z: each row is (y_{t-1}, ..., y_{t-d})
  # Dimensions: (T_obs - d) x (N * d)
  Z <- do.call(cbind, lapply(seq_len(d), function(lag) {
    seg[(d - lag + 1):(T_obs - lag), , drop = FALSE]
  }))

  # Fitted values: Y_hat = Z %*% t(B), dimension (T_obs - d) x N
  Y_actual <- seg[(d + 1):T_obs, , drop = FALSE]
  Y_hat    <- Z %*% t(B)

  residuals <- Y_actual - Y_hat
  colnames(residuals) <- colnames(seg)
  return(residuals)
}


#' Check residual autocorrelation for all VAR segments
#'
#' Applies the multivariate Portmanteau (Ljung-Box) test to the residuals
#' of each fitted segment using MTS::mq(). A non-significant result
#' (p > 0.05) indicates no remaining autocorrelation and supports model
#' adequacy. Per-equation univariate Box.test results are also returned
#' for detailed inspection.
#'
#' @param fits     List of NGC fit objects (one per segment)
#' @param stat_dat List of stationary data matrices (one per segment)
#' @param d        VAR lag order
#' @param n_brk    Number of change points (P)
#' @param lags     Number of lags for the Portmanteau test (default 10)
#' @return         Data frame of per-equation Ljung-Box p-values per segment

check_residual_autocorrelation <- function(fits, stat_dat, d, n_brk, lags = 10) {
  results <- data.frame(
    Segment   = integer(),
    Location  = character(),
    LB_Stat   = numeric(),
    LB_Pval   = numeric(),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(n_brk + 1)) {
    message(sprintf("Residual autocorrelation check: segment %d", j))

    resid_mat <- compute_ngc_residuals(
      fit = fits[[j]],
      seg = stat_dat[[j]],
      d   = d
    )

    # Multivariate Portmanteau test on full residual matrix
    # Prints test statistic and p-value for the segment
    cat(sprintf("\n--- Segment %d: Multivariate Portmanteau test ---\n", j))
    MTS::mq(resid_mat, lag = lags)

    # Univariate Ljung-Box test per equation for detailed inspection
    for (loc in colnames(resid_mat)) {
      lb <- Box.test(resid_mat[, loc], lag = lags, type = "Ljung-Box")
      results <- rbind(results, data.frame(
        Segment  = j,
        Location = loc,
        LB_Stat  = round(lb$statistic, 4),
        LB_Pval  = round(lb$p.value,   4)
      ))
    }
  }

  # Flag any equations where residual autocorrelation is significant
  n_flagged <- sum(results$LB_Pval < 0.05, na.rm = TRUE)
  if (n_flagged > 0) {
    warning(sprintf(
      paste0("%d series-segment combinations show significant residual ",
             "autocorrelation (p < 0.05). Consider increasing lag order."),
      n_flagged
    ))
  } else {
    message("All equations pass residual autocorrelation check (p > 0.05).")
  }

  return(results)
}


#' Check stability of fitted VAR models across all segments
#'
#' Constructs the companion matrix from the NGC lasso VAR coefficient
#' matrices and checks that all eigenvalues lie strictly inside the unit
#' circle (modulus < 1). A stable VAR implies the process is stationary
#' and that impulse responses decay to zero over time.
#'
#' The companion matrix for a VAR(d) with N series is (N*d) x (N*d):
#'   [ B1  B2  ...  Bd ]
#'   [ I   0   ...  0  ]
#'   [ 0   I   ...  0  ]
#'   [         ...     ]
#' where B1,...,Bd are the N x N coefficient matrices at each lag,
#' stacked from the NGC coefficient matrix B (N x N*d).
#'
#' @param fits   List of NGC fit objects (one per segment)
#' @param N      Number of time series
#' @param d      VAR lag order
#' @param n_brk  Number of change points (P)
#' @return       Data frame with max eigenvalue modulus per segment

check_var_stability <- function(fits, N, d, n_brk) {
  results <- data.frame(
    Segment      = integer(),
    Max_Modulus  = numeric(),
    Stable       = logical(),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(n_brk + 1)) {
    B <- fits[[j]]$B   # N x (N*d) coefficient matrix from NGC

    # Build top block of companion matrix: [B1 | B2 | ... | Bd]
    # B is already arranged as [B1 | B2 | ... | Bd] (columns ordered by lag)
    top_block <- B   # N x (N*d)

    # Build identity and zero blocks for lower rows of companion matrix
    if (d > 1) {
      lower_block <- cbind(
        diag(N * (d - 1)),                    # identity block
        matrix(0, nrow = N * (d - 1), ncol = N)  # zero block
      )
      companion <- rbind(top_block, lower_block)
    } else {
      companion <- top_block
    }

    # Eigenvalues of the companion matrix
    eig_vals   <- eigen(companion, only.values = TRUE)$values
    max_modulus <- max(Mod(eig_vals))
    is_stable   <- max_modulus < 1

    if (!is_stable) {
      warning(sprintf(
        "Segment %d is UNSTABLE: max eigenvalue modulus = %.4f (>= 1).",
        j, max_modulus
      ))
    } else {
      message(sprintf(
        "Segment %d stable: max eigenvalue modulus = %.4f.",
        j, max_modulus
      ))
    }

    results <- rbind(results, data.frame(
      Segment     = j,
      Max_Modulus = round(max_modulus, 4),
      Stable      = is_stable
    ))
  }

  return(results)
}


#' Compute node importance measures for each segment network
#'
#' Computes indegree, outdegree, betweenness centrality, in-strength,
#' out-strength, and several measures from the influential package.
#' In- and out-strength use edge weights from the lasso VAR coefficients.
#'
#' @param nw      List of igraph objects (one per segment)
#' @param n_brk   Number of change points (P)
#' @return Data frame of node importance measures across all segments

compute_node_importance <- function(nw, n_brk) {
  col_names <- c(
    "Section", "City",
    "Indegree", "Outdegree", "Betweenness",
    "Strength_in", "Strength_out",
    "NeighConnect_in", "NeighConnect_out",
    "H_Index_in", "H_Index_out",
    "Coll_Inf_in", "Coll_Inf_out",
    "IVI_in", "IVI_out", "IVI_all",
    "Closeness", "Eigenvector"
  )

  node_imp <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(node_imp) <- col_names

  for (i in 1:(n_brk + 1)) {
    g_i <- igraph::simplify(nw[[i]], remove.multiple = FALSE)

    rows <- cbind(
      Section            = i,
      City               = igraph::V(g_i)$names,
      Indegree           = igraph::degree(g_i, mode = "in"),
      Outdegree          = igraph::degree(g_i, mode = "out"),
      Betweenness        = round(igraph::betweenness(g_i, normalized = TRUE), 4),
      # In- and out-strength: sum of absolute lasso VAR edge weights
      Strength_in        = igraph::strength(g_i, mode = "in",
                                             weights = igraph::E(g_i)$weight),
      Strength_out       = igraph::strength(g_i, mode = "out",
                                             weights = igraph::E(g_i)$weight),
      NeighConnect_in    = round(influential::neighborhood.connectivity(g_i, mode = "in"), 4),
      NeighConnect_out   = round(influential::neighborhood.connectivity(g_i, mode = "out"), 4),
      H_Index_in         = influential::h_index(g_i, mode = "in"),
      H_Index_out        = influential::h_index(g_i, mode = "out"),
      Coll_Inf_in        = influential::collective.influence(g_i, mode = "in"),
      Coll_Inf_out       = influential::collective.influence(g_i, mode = "out"),
      IVI_in             = round(influential::ivi(g_i, directed = TRUE, mode = "in"), 4),
      IVI_out            = round(influential::ivi(g_i, directed = TRUE, mode = "out"), 4),
      IVI_all            = round(influential::ivi(g_i, directed = TRUE, mode = "all"), 4),
      Closeness          = igraph::closeness(g_i, mode = "out"),
      Eigenvector        = igraph::eigen_centrality(g_i, directed = TRUE)$vector
    )

    node_imp <- rbind(node_imp, rows)
  }

  # Ensure numeric columns are stored as numeric
  numeric_cols <- setdiff(col_names, c("Section", "City"))
  node_imp[numeric_cols] <- lapply(node_imp[numeric_cols], as.numeric)
  node_imp$Section <- as.integer(node_imp$Section)

  return(node_imp)
}


##OHa regions
Regions <- tibble::tribble(
  ~City,           ~Regions,
  # Region 1 — Northwest / Portland Metro
  "Astoria", 1, "St. Helens", 1, "Tillamook", 1,
  "Forest.Grove", 1, "Hillsboro", 1, "Rock.Creek", 1,
  "Durham", 1, "Portland", 1, "Sandy", 1,
  "Wilsonville", 1, "Canby", 1,
  # Region 2 — Mid-Willamette Valley / Central Coast
  "McMinnville", 2, "Newberg", 2, "Sheridan", 2,
  "Dallas", 2, "Salem", 2, "Woodburn", 2,
  "Silverton", 2, "Stayton", 2, "Albany", 2,
  "Corvallis", 2, "Lincoln.City", 2, "Newport", 2,
  "Siletz", 2,
  # Region 3 — South Willamette / South Coast
  "Eugene", 3, "Cottage.Grove", 3, "Florence", 3,
  "Roseburg", 3, "North.Bend", 3, "Port.Orford", 3,
  "Gold.Beach", 3,
  # Region 5 — Southern Oregon
  "Grants.Pass", 5, "Medford", 5, "Ashland", 5,
  # Region 6 — Columbia Gorge / North Central
  "Hood.River", 6, "The.Dalles", 6,
  # Region 7 — Central / South Central
  "Warm.Springs", 7, "Redmond", 7, "Bend", 7,
  "Sunriver", 7, "Klamath.Falls", 7,
  # Region 9 — Eastern Oregon
  "Boardman", 4, "Umatilla", 4, "Hermiston", 4,
  "Pendleton", 4, "La.Grande", 4, "Baker.City",4,
  "Ontario", 4
)

# ------------------------------------------------------------------------------
# 3. Load and preprocess data
# ------------------------------------------------------------------------------

dat <- read_excel(DATA_PATH, sheet = "COVID", guess_max = 10000)

# Retain relevant columns
df_n <- dat  %>% 
  subset(select = c(Sample_Date,LogCopiesPerL,LogCopiesPerDayPerPerson_NormtoFlowPopRec, Location, County, Site))

df <- df_n %>%
  # Anchor each row to the start of its ISO week (Monday)
  mutate(week_date = floor_date(Sample_Date, unit = "week", week_start = 1)) %>%
  # Average across multiple samples for the same location in the same week
  group_by(Location, week_date) %>%
  summarise(mean_log_copies = mean(LogCopiesPerL, na.rm = TRUE), .groups = "drop") %>% ungroup() %>%
  # Sort so the wide output is in date order
  arrange(week_date) %>%
  # Pivot: one column per location
  pivot_wider(names_from = Location, values_from = mean_log_copies)

# df <- df_n %>%
#   mutate(week_date = floor_date(Sample_Date, unit = "week", week_start = 1)) %>%
#   group_by(Location, week_date) %>%
#   summarise(mean_log_copies = mean(LogCopiesPerL, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from = Location, values_from = mean_log_copies) %>%
#   arrange(week_date)

# Restrict to active surveillance window (weeks 11–290)
# Week 0 = first week of September 2020; week 229 = last week January 2025
#df <- df[11:290, ]

# Lookup table: location -> county
df_county <- na.omit(
  distinct(data.frame(loc = make.names(dat$Location), county = dat$County))
)


# ------------------------------------------------------------------------------
# 4. Missingness assessment and location filtering
# ------------------------------------------------------------------------------

# Compute % missing per location
nas <- data.frame(
  names = colnames(df),
  val   = as.numeric(colSums(is.na(df)) / nrow(df)) * 100
)

# Plot missingness histogram to inform threshold choice
hist(nas$val[-(1)], breaks = 70,
     main = "% Missing values per location",
     xlab = "% Missing", col = "steelblue")

# Exclude locations exceeding the missingness threshold
# (PCT_THRS = 40%: majority of locations fall below this threshold)
df <-  df %>% subset(select = which(nas[,2] < PCT_THRS))
colnames(df) <- make.names(colnames(df))

# Additional manual exclusions: locations with data quality issues
# NOTE: these exclusions should be documented in the paper supplementary material
df <- df %>% subset(select = -c(St..Helens,
                                Silverton,
                                Siletz))

# Recompute missingness after filtering
nas <- data.frame(
  names = colnames(df),
  val   = as.numeric(colSums(is.na(df)) / nrow(df)) * 100
)
nas <- nas[-1, ]  # drop week_start row

N_LOC    <- ncol(df) - 1  # number of retained locations (N = 27)
message(sprintf("Retained %d locations after missingness filtering.", N_LOC))

# Extract numeric matrix of wastewater measurements (locations as columns)
ww_samp        <- as.matrix(df[, -1])
colnames(ww_samp) <- colnames(df)[-1]

# State average (used for plotting change points against context)
avg_val <- rowMeans(ww_samp, na.rm = TRUE)


# ------------------------------------------------------------------------------
# 5. Imputation step 1: missing values due to non-participation
#
#    Missing values arising from voluntary non-participation are assumed
#    missing at random (MAR). They are imputed via Kalman filter using a
#    StructTS structural time series model (level + trend + seasonal components)
#    implemented via na_kalman(..., model = "StructTS") in the imputeTS package.
#    This step occurs BEFORE change point detection.
# ------------------------------------------------------------------------------

ww_ip <- matrix(NA, nrow = nrow(ww_samp), ncol = ncol(ww_samp))
colnames(ww_ip) <- colnames(ww_samp)

for (i in seq_len(ncol(ww_samp))) {
  ww_ip[, i] <- na_kalman(ww_samp[, i], model = "StructTS")
}
ww_ip <- as.data.frame(ww_ip)


# ------------------------------------------------------------------------------
# 6. Change point detection (Bai & Perron via strucchange)
#
#    Structural change points are detected for each univariate time series
#    separately using the breakpoints() function in the strucchange package.
#    The number of change points P = N_BRK = 7 was selected based on the
#    RSS criterion: the value at which RSS showed diminishing returns across
#    the majority of locations.
# ------------------------------------------------------------------------------
seg_colors <- c(
  "brown4", "blue4", "chartreuse4", "darkorchid", "cornflowerblue",
  "coral2", "salmon4", "cyan4", "brown4","maroon","pink4","violet"
)

op_bp  <- list()   # raw breakpoints objects per location
bp_nat <- data.frame(matrix(ncol = 4, nrow = 0))  # tidy breakpoint table
col    <- data.frame(loc = colnames(ww_ip))
col    <- left_join(col, df_county, by = c("loc" = "loc"))

for (i in seq_len(N_LOC)) {
  # Prepare regression data for breakpoints: y ~ t (simple linear trend model)
  bp_data <- data.frame(
    y = ww_ip[, i],
    t = seq_len(nrow(ww_ip))
   # y_1 = lag(ww_ip[, i], 1)
  )

  # First call: determine optimal breakpoint structure via RSS
  bp_full <- breakpoints(y ~ t, h = N_LOC, data = bp_data)

  # Second call: extract exactly N_BRK breakpoints
  op_bp[[i]] <- breakpoints(bp_full, breaks = N_BRK)

  # Store breakpoints including series start (1) and end (nrow)
  bp_nat <- rbind(
    bp_nat,
    cbind(
      col[i, 1], col[i, 2],
      c(1, op_bp[[i]]$breakpoints, nrow(ww_ip)),
      as.character(df$week_date[as.numeric(c(1, op_bp[[i]]$breakpoints, nrow(ww_ip)))])
    )
  )
}

colnames(bp_nat) <- c("Location", "County", "Breakpoint","Date")
bp_nat$Breakpoint <- as.numeric(bp_nat$Breakpoint)
bp_nat$brkpt      <- rep(0:(N_BRK + 1), times = N_LOC)
bp_nat$Location   <- factor(bp_nat$Location,
                             levels = bp_nat$Location[bp_nat$brkpt == 2])

bp_dt <- as.data.frame(bp_nat) %>% select(c(Date, brkpt, Breakpoint)) %>% mutate(Date = as.Date(Date))
Date_result <- bp_dt %>%
  group_by(brkpt) %>%
  summarise(
    first_date  = min(Date),
    last_date   = max(Date),
    mean_date   = mean(Date),
    median_date = median(Date),
    median_point  = Breakpoint[which(Date  == median(Date))]
  )

# Identify the earliest and latest change points per breakpoint index
# (used as vertical reference lines in the change point plot)
highest_rows <- bp_nat %>%
  filter(Breakpoint != 1 & Breakpoint != nrow(df)) %>%
  group_by(brkpt) %>%
  slice_max(Breakpoint)

lowest_rows <- bp_nat %>%
  filter(Breakpoint != 1 & Breakpoint != nrow(df)) %>%
  group_by(brkpt) %>%
  slice_min(Breakpoint)

# Collect BIC values for each location and each candidate number of breakpoints
max_breaks <- 10
bic_mat <- matrix(NA, nrow = N_LOC, ncol = max_breaks)
RSS_mat <- matrix(NA, nrow = N_LOC, ncol = max_breaks)
rownames(bic_mat) <- col$loc
rownames(RSS_mat) <- col$loc

for (i in seq_len(N_LOC)) {
  bp_data <- data.frame(y = ww_ip[, i], t = seq_len(nrow(ww_ip)))
  bp_full <- breakpoints(y ~ t, h = N_LOC, data = bp_data)
  bic_vals <- scales::rescale(AIC(bp_full, k = log(nrow(ww_ip))), to = c(-1,1))  # BIC)
  RSS_vals <- scales::rescale(summary(bp_full)$RSS[1,], to=c(-1,1)) #RSS
  # pad with NA if fewer breaks were estimable
  bic_mat[i, seq_along(bic_vals)] <- bic_vals
  RSS_mat[i, seq_along(RSS_vals)] <- RSS_vals
  
}

# Plot as heatmap
bic_df <- as.data.frame(bic_mat) %>%
  rownames_to_column("Location") %>%
  pivot_longer(-Location, names_to = "Breaks", values_to = "BIC") %>%
  mutate(Breaks = as.integer(gsub("V", "", Breaks)))

ggplot(bic_df, aes(x = as.factor(Breaks), y = Location, fill = BIC)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", na.value = "grey80") +
  geom_vline(xintercept = N_BRK, linetype = "dashed", color = "white") +
  labs(title = "BIC by number of breakpoints per location",
       x = "Number of breakpoints", y = "Location") + theme(legend.position = "none")+
  theme_minimal()
ggsave(paste0(SAVE_DIR, "BICBreakPoints.png"))

#Plot RSS as heatmap
RSS_df <- as.data.frame(RSS_mat) %>%
  rownames_to_column("Location") %>%
  pivot_longer(-Location, names_to = "Breaks", values_to = "RSS") %>%
  mutate(Breaks = as.integer(gsub("V", "", Breaks)))

ggplot(RSS_df, aes(x = Breaks, y = Location, fill = RSS)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", na.value = "grey80") +
  geom_vline(xintercept = N_BRK, linetype = "dashed", color = "white") +
  labs(title = "RSS by number of breakpoints per location",
       x = "Number of breakpoints", y = "Location") +
  theme_minimal()
ggsave(paste0(SAVE_DIR, "RSSBreakPoints.png"))


# ------------------------------------------------------------------------------
# 7. Change point visualisation (Figure: brkPts.png)
# ------------------------------------------------------------------------------
bp_nat$Location <- as.character(bp_nat$Location)

bp_nat %>%
  ggplot() +
  geom_point(aes(y = Location, x = Breakpoint, color = as.factor(brkpt))) +
  scale_color_manual(values = seg_colors)+
  geom_line(aes(y = Location,  x = Breakpoint)) +
  #geom_vline(data = highest_rows, aes(xintercept = Breakpoint),
  #           color = "darkblue", linetype = 2) +
  #geom_vline(data = lowest_rows,  aes(xintercept = Breakpoint),
  #           color = "darkred",  linetype = 2) +
  scale_x_continuous(guide = guide_axis(angle = 45),
                     limits = c(0, nrow(ww_samp)),
                     breaks = scales::breaks_pretty(n = 20)) +
  labs(title = "Change point distribution across locations",
       y = "Location", x = "Week index") + guides(color = "none")+
  theme_minimal()

ggsave(paste0(SAVE_DIR, "brkPts.png"))

# --- week <-> date mapping (week 1 = 2020-08-31, weekly spacing) ---
anchor <- as.Date("2020-08-31")
week_to_date <- function(w) anchor + (w - 1) * 7          # numeric week -> Date
date_to_week <- function(d) as.numeric(d - anchor) / 7 + 1 # Date -> numeric week

wk_max <- nrow(ww_samp)   # 295

# first-of-month dates that fall inside the plotted window, for the bottom axis
month_dates  <- seq(from = as.Date(cut(week_to_date(1),  "month")),
                    to   = week_to_date(wk_max),
                    by   = "3 month")
month_breaks <- date_to_week(month_dates)   # positions on the (week) x-scale

bp_nat %>%
  ggplot() +
  geom_line(aes(y = Location, x = Breakpoint)) +
  geom_point(aes(y = Location, x = Breakpoint, color = as.factor(brkpt))) +
  scale_color_manual(values = seg_colors) +
  scale_x_continuous(
    # PRIMARY (bottom) axis = dates, one tick per first-of-month
    name   = "Date",
    limits = c(0, wk_max),
    breaks = month_breaks,
    labels = format(month_dates, "%b %Y"),
    guide  = guide_axis(angle = 45),
    # SECONDARY (top) axis = week index, as before
    sec.axis = sec_axis(
      transform = ~ .,                       # identity: same underlying scale
      name   = "Week index",
      breaks = scales::breaks_pretty(n = 20)
    )
  ) + 
  labs(title = "Change point distribution across locations", y = "Location", x = " ") +
  guides(color = "none") +
  theme_minimal() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

ggsave(paste0(SAVE_DIR, "brkPtsNew.png"))

## Plot the time series for each location and the breakpont for that time series
#Color based on imputed vs not 
cts <- c("Hermiston", "Corvallis", "Portland", "Salem" )

for(loc in colnames(ww_samp)){
  temp_df <- as.data.frame(cbind(t= 1:295, 
                                 val = ww_ip %>% select(c(loc)), 
                                 imp = as.numeric(is.na(as.data.frame(ww_samp) %>% select(c(loc))))
                                 ))
  colnames(temp_df) <- c("t", "val", "imp")
  bp_tmp <- bp_nat %>% filter(Location == loc)
  
  print(temp_df %>% ggplot()+
    geom_line(aes(x = t , y= val), alpha =0.5)+
    geom_point(aes(x = t , y= val, color = as.factor(imp)))+
    scale_color_manual(values = seg_colors)+
    geom_vline(xintercept = bp_tmp$Breakpoint, color = "darkgreen")+
    ggtitle(loc)+ 
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA), # Panel area
          plot.background = element_rect(fill = "white", color = NA),  # Entire image area
          panel.grid.major = element_blank(),                          # Remove major grids
          panel.grid.minor = element_blank()))
 # ggsave(paste0(SAVE_DIR,loc,"BP.png"))
  
}
rowSums(!is.na(ww_samp))
sort(colSums(!is.na(ww_samp)))

##Combined plots
library(patchwork)

cts <- c("Hermiston", "Corvallis", "Portland", "Salem","Eugene") 
# "Woodburn", "Stayton","McMinnville", "Astoria","Bend" 

plot_list <- lapply(cts, function(loc) {
  temp_df <- as.data.frame(cbind(
    t   = 1:295,
    val = ww_ip %>% select(all_of(loc)),
    imp = as.numeric(is.na(as.data.frame(ww_samp) %>% select(all_of(loc))))
  ))
  colnames(temp_df) <- c("t", "val", "imp")
  bp_tmp <- bp_nat %>% filter(Location == loc)
  
  ggplot(temp_df) +
    geom_line(aes(x = t, y = val), alpha = 0.5) +
    geom_point(aes(x = t, y = val, color = as.factor(imp))) +
    scale_color_manual(values = seg_colors) +
    geom_vline(xintercept = bp_tmp$Breakpoint, color = "darkgreen") +
    scale_x_continuous(breaks = seq(0, 295, by = 30)) +
    labs(y = "Mean log copies per L", x = "Week") +
    ggtitle(loc) +
    theme(
      legend.position  = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
})

# remove x-axis title from all but the last plot
plot_list[-length(plot_list)] <- lapply(
  plot_list[-length(plot_list)],
  function(p) p + theme(axis.title.x = element_blank())
)

combined <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(axis_titles = "collect_y")

print(combined)
ggsave(paste0(SAVE_DIR, "combined_BP1.png"), combined,width = 8, height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# 8. Segment boundary table
#
#    For each segment group j, the common window is defined as
#    [min_i(q_{i,j-1}), max_i(q_{ij})], i.e. the union across all series.
# ------------------------------------------------------------------------------

# Build a tidy segment boundary table (start/end per series per segment)
bp <- data.frame(matrix(nrow = 0, ncol = 3))
for (i in seq_len(N_LOC)) {
  tmp <- c(1, op_bp[[i]]$breakpoints, nrow(ww_ip))
  bp  <- rbind(bp, cbind(rep(i, length(tmp)), tmp, rep(col$loc[i], length(tmp))))
}
colnames(bp) <- c("TS", "BkPtS", "Location")
bp$BkPtS     <- as.numeric(bp$BkPtS)
bp$Date      <- df$week_date[bp$BkPtS]

bp <- bp %>%
  group_by(TS) %>%
  mutate(BkPtE = as.numeric(lead(BkPtS, 1, default = NA))) %>%
  filter(!is.na(BkPtE)) %>%
  mutate(
    range  = paste0(BkPtS, "-", BkPtE),
    id     = row_number()
  ) %>%
  ungroup()

# For each segment group j, compute the common window boundaries:
# minV = min_i(q_{i,j-1}), maxV = max_i(q_{ij})
rng <- bp %>%
  group_by(id) %>%
  summarise(
    minV    = min(BkPtS),
    maxV    = max(BkPtE),
    .groups = "drop"
  ) %>%
  mutate(
    minDate = df$week_start[minV],
    maxDate = df$week_start[maxV]
  )


# ------------------------------------------------------------------------------
# 9. Segment alignment visualisation (Figure: CPperCity.png)
# ------------------------------------------------------------------------------

#bp$Location <- sort(bp$Location)
bp_plot <- bp %>%
  mutate(
    xRangeL = round(
      rep(seq(min(rowMeans(ww_ip, na.rm = TRUE)),
              max(rowMeans(ww_ip, na.rm = TRUE)),
              length.out = N_LOC),
          each = N_BRK + 1), 2),
    xRangeH = rep(c(unique(xRangeL)[-1], 5.5), each = N_BRK + 1)
  )

# ggplot() +
#   geom_rect(aes(
#     xmin = as.numeric(bp_plot$BkPtS), xmax = as.numeric(bp_plot$BkPtE),
#     ymin = bp_plot$xRangeL,           ymax = bp_plot$xRangeH,
#     fill = as.factor(bp_plot$id)
#   ), alpha = 1) +
#   geom_line(aes(y = rowMeans(ww_ip, na.rm = TRUE),
#                 x = seq_len(nrow(ww_ip)))) +
#   geom_text(aes(x = -40, y = as.numeric(unique(bp_plot$xRangeL)),
#             label = sort(unique(bp_plot$Location)) ),
#             size = 3, vjust = 0, hjust = 0, color = "blue4") +
#   scale_fill_manual(values = seg_colors) +
#   scale_x_continuous(guide = guide_axis(angle = 45),
#                      #limits = c(0, nrow(ww_samp)),
#                      breaks = scales::breaks_pretty(n = 10)
#                      ) +
#   labs(x = "Sample week", y = "Mean log copies per L",
#        title = "Change points per location against state average") +
#   #scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
#   theme_minimal() +
#   theme(legend.position = "none")
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# one row per location, positioned at the middle of its band
# vertical lanes currently in the plot, ordered top -> bottom
lanes <- unique(bp_plot[c("xRangeL", "xRangeH")])
lanes <- lanes[order(lanes$xRangeL, decreasing = FALSE), ]

# locations A -> Z  (as.character guards against factor-level ordering)
locs <- sort(unique(as.character(bp_plot$Location)))

# pair alphabetical rank with lane position: first name -> top lane
lane_map <- data.frame(Location = locs,
                       xRangeL  = lanes$xRangeL,
                       xRangeH  = lanes$xRangeH,
                       stringsAsFactors = FALSE)

# swap old lanes for the alphabetical ones
bp_plot$xRangeL <- bp_plot$xRangeH <- NULL
bp_plot <- merge(bp_plot, lane_map, by = "Location")

loc_axis <- unique(bp_plot[c("Location", "xRangeL", "xRangeH")])
loc_axis$ymid <- (loc_axis$xRangeL + loc_axis$xRangeH) / 2

x_max <- max(as.numeric(bp_plot$BkPtE), nrow(ww_ip), na.rm = TRUE)

ggplot() +
  geom_rect(aes(
    xmin = as.numeric(bp_plot$BkPtS), xmax = as.numeric(bp_plot$BkPtE),
    ymin = bp_plot$xRangeL,           ymax = bp_plot$xRangeH,
    fill = as.factor(bp_plot$id)
  ), alpha = 1) +
  geom_line(aes(y = rowMeans(ww_ip, na.rm = TRUE),
                x = seq_len(nrow(ww_ip)))) +
  scale_fill_manual(values = seg_colors) +
  scale_x_continuous(breaks = seq(0, x_max, by = 20),
                     guide  = guide_axis(angle = 45)) +
  scale_y_continuous(
    breaks   = loc_axis$ymid,
    labels   = loc_axis$Location,
    sec.axis = sec_axis(~ ., name   = "Mean log copies per L",
                        breaks = scales::breaks_pretty(n = 6))
  ) +
  labs(x = "Sample week", y = NULL,
       title = "Change points per location against state average") +
  theme_minimal() +
  theme(legend.position    = "none",
        axis.text.y.left   = element_text(color = "blue4", size = 7))

ggsave(paste0(SAVE_DIR, "CPperCity.png"))

bp_plot <- left_join(bp_plot, Regions, join_by(Location == City)) 
bp_plot <- left_join(bp_plot, col, join_by( Location == loc))
bp_plot$county <- as.factor(bp_plot$county)

# Reusable: reorder y-bands by an arbitrary grouping column
build_layout <- function(df, group_col) {
  band_w <- median(df$xRangeH - df$xRangeL)
  y0     <- min(df$xRangeL)
  
  loc_order <- df %>%
    distinct(Location, !!sym(group_col)) %>%
    arrange(!!sym(group_col), Location) %>%
    mutate(rank  = row_number(),
           new_L = y0 + (rank - 1) * band_w,
           new_H = y0 +  rank      * band_w)
  
  df_out <- df %>%
    select(-xRangeL, -xRangeH) %>%
    left_join(loc_order %>% select(Location, !!sym(group_col),
                                   xRangeL = new_L, xRangeH = new_H),
              by = c("Location", group_col))
  
  list(data = df_out, labels = loc_order %>% arrange(new_L))
}

# Region palette (7 regions used in your data)
region_colors <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3",
                   "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02",
                   "7" = "#a6761d")

# County palette: build once from the 18 county levels
county_levels <- levels(bp_plot$county)
county_colors <- setNames(
  scales::hue_pal()(length(county_levels)),  # or viridis::viridis(length(county_levels))
  county_levels
)

make_plot <- function(layout, group_col, group_colors, title) {
  bp   <- layout$data
  labs <- layout$labels
  label_cols <- group_colors[as.character(labs[[group_col]])]
  
  ggplot(bp) +
    geom_rect(aes(xmin = as.numeric(BkPtS), xmax = as.numeric(BkPtE),
                  ymin = xRangeL,           ymax = xRangeH,
                  fill = as.factor(id)), alpha = 1) +
    geom_line(data = data.frame(x = seq_len(nrow(ww_ip)),
                                y = rowMeans(ww_ip, na.rm = TRUE)),
              aes(x = x, y = y)) +
    geom_text(data = labs,
              aes(x = -45, y = new_L, label = Location),
              color = label_cols,
              size = 3, vjust = 0, hjust = 0) +
    scale_fill_manual(values = seg_colors) +
    labs(x = "Sample week", y = "Mean log copies per L", title = title) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Build both plots
p_region <- make_plot(build_layout(bp_plot, "Regions"),
                      "Regions", region_colors,
                      "Change points by location, grouped by OHA region")
p_region
ggsave(paste0(SAVE_DIR, "CPperCity_OHAregion.png"))

p_county <- make_plot(build_layout(bp_plot, "county"),
                      "county", county_colors,
                      "Change points by location, grouped by county")
p_county
ggsave(paste0(SAVE_DIR, "CPperCity_county.png"))

# ------------------------------------------------------------------------------
# 10. Imputation step 2: boundary padding for segment alignment
#
#     Each series is extended to the common window [minV, maxV] for its
#     segment group. Missing boundary values (where a series does not
#     naturally span the full window) are imputed via Kalman filter using
#     the StructTS model, consistent with the framework description.
#     NOTE: imputation occurs at segment boundaries near detected change
#     points; the structural model fitted to the segment interior may not
#     perfectly represent the boundary dynamics.
# ------------------------------------------------------------------------------

TSsec <- list()

for (i in seq_len(N_LOC)) {
  seg_list <- list()

  for (j in seq_len(N_BRK + 1)) {
    bp_ss <- bp[bp$TS == i, ]
    strt  <- bp_ss[bp_ss$id == j, ]$BkPtS
    end   <- bp_ss[bp_ss$id == j, ]$BkPtE

    seg_vals <- ww_ip[strt:end, i]

    # Pad front (left boundary): series starts after common window start
    front_pad <- strt - rng[j, ]$minV
    if (front_pad > 0) seg_vals <- c(rep(NA, front_pad), seg_vals)

    # Pad end (right boundary): series ends before common window end
    end_pad <- rng[j, ]$maxV - end
    if (end_pad > 0) seg_vals <- c(seg_vals, rep(NA, end_pad))

    # Impute boundary NAs using Kalman filter with StructTS model
    seg_list[[j]] <- na_kalman(seg_vals, model = "StructTS")
  }

  TSsec[[i]] <- seg_list
}

# Reshape into a list of N_LOC x l_j matrices, one per segment group j
df_split <- list()
for (j in seq_len(N_BRK + 1)) {
  mat <- as.matrix(sapply(TSsec, "[[", j))
  colnames(mat) <- col$loc
  df_split[[j]] <- mat
}


# ------------------------------------------------------------------------------
# 11. Stationarity assessment and first differencing
#
#     ADF test is applied to each series in each segment. If any series has
#     a unit root (non-stationary), first differencing is applied to that
#     segment. Here we apply first differencing uniformly since the VAR model
#     requires stationarity. Results therefore reflect Granger causality
#     between *changes* in wastewater concentration, not levels.
#     For level-based inference, a VECM with cointegration testing would
#     be appropriate (see Johansen test via urca::ca.jo).
# ------------------------------------------------------------------------------

adf_results <- data.frame(
  matrix(ncol = 6, nrow = 0)
)
colnames(adf_results) <- c(
  "Subsection", "Location",
  "Orig_Stat", "Orig_Pval",
  "Diff_Stat", "Diff_Pval"
)

stat_dat <- list()

for (i in seq_len(N_BRK + 1)) {
  message(sprintf("Stationarity check: segment %d", i))
  seg    <- df_split[[i]]
  cn     <- colnames(seg)
  diff_mat <- matrix(ncol = ncol(seg), nrow = nrow(seg) - 1)

  for (j in seq_len(ncol(seg))) {
    diff_vec <- diff(seg[, j], lag = 1)
    diff_mat[, j] <- diff_vec

    orig_test <- adf.test(seg[, j])
    diff_test <- adf.test(diff_vec)

    adf_results <- rbind(
      adf_results,
      data.frame(
        Subsection = i,
        Location   = cn[j],
        Orig_Stat  = round(orig_test$statistic[[1]], 4),
        Orig_Pval  = round(orig_test$p.value, 4),
        Diff_Stat  = round(diff_test$statistic[[1]], 4),
        Diff_Pval  = round(diff_test$p.value, 4)
      )
    )
  }
  colnames(diff_mat) <- cn
  stat_dat[[i]] <- diff_mat
}

# Summary of stationarity results
message(sprintf(
  "Non-stationary segments (original, p > 0.05): %d",
  sum(adf_results$Orig_Pval > 0.05, na.rm = TRUE)
))
message(sprintf(
  "Non-stationary segments after differencing (p > 0.05): %d",
  sum(adf_results$Diff_Pval > 0.05, na.rm = TRUE)
))


# ------------------------------------------------------------------------------
# 12. Lasso-penalized VAR fitting and Granger causality network construction
#
#     A lasso-penalized VAR(p) model is fitted to each aligned and stationary
#     segment using the NGC package (Shojaie & Michailidis 2010). The lasso
#     penalty parameter is selected via cross-validation internally by NGC
#     (via cv.glmnet from the glmnet package). The resulting sparse coefficient
#     matrices define the directed Granger causality network for each segment.
# ------------------------------------------------------------------------------
outdegree_imp <- matrix(0, nrow = 0, ncol =4) 
indegree_imp <- matrix(0, nrow = 0, ncol =4)
outstrength_imp <- matrix(0, nrow = 0, ncol =4) 
betweenness_imp <- matrix(0, nrow = 0, ncol =4)

for(i in 1:1000){
print(paste0("The ", i, "th iteration"))
var_output <- fit_lasso_var(
  df_ip   = stat_dat,
  d       = LAG_D,
  n_brk   = N_BRK,
  ret_fit = TRUE
)

G        <- var_output$graphs   # list of igraph objects, one per segment
var_fits <- var_output$fits     # list of raw NGC fit objects

{
# ------------------------------------------------------------------------------
# 13. VAR model diagnostics
#
#     Two checks are applied to each fitted segment model:
#
#     (a) Residual autocorrelation — Portmanteau test (multivariate) and
#         Ljung-Box test (per equation). A non-significant result (p > 0.05)
#         indicates the model has adequately captured the temporal dependence
#         structure. Significant autocorrelation suggests the lag order p
#         may need to be increased.
#
#     (b) Stability — all eigenvalues of the companion matrix must have
#         modulus strictly less than 1. An unstable model implies the fitted
#         VAR process is explosive and results should be interpreted with
#         caution. Note that lasso penalisation shrinks coefficients toward
#         zero, which tends to promote stability; instability would be
#         unusual but should still be verified.
# ------------------------------------------------------------------------------

# #(a) Residual autocorrelation check
# autocorr_results <- check_residual_autocorrelation(
#   fits     = var_fits,
#   stat_dat = stat_dat,
#   d        = LAG_D,
#   n_brk    = N_BRK,
#   lags     = 10
# )
# 
# # # Inspect flagged equations (significant autocorrelation at p < 0.05)
# flagged_autocorr <- autocorr_results %>% filter(LB_Pval < 0.05)
# if (nrow(flagged_autocorr) > 0) {
#   message("Series-segment pairs with significant residual autocorrelation:")
#   print(flagged_autocorr)
# } else {
#   message("No significant residual autocorrelation detected in any segment.")
# }
# 
# # (b) VAR stability check
# stability_results <- check_var_stability(
#   fits  = var_fits,
#   N     = N_LOC,
#   d     = LAG_D,
#   n_brk = N_BRK
# )
# 
# # # Summary table of stability results
# print(stability_results)
# 
# unstable_segs <- stability_results %>% filter(!Stable)
# if (nrow(unstable_segs) > 0) {
#   warning(sprintf(
#     "Unstable VAR found in segment(s): %s. Interpret network results with caution.",
#     paste(unstable_segs$Segment, collapse = ", ")
#   ))
# } else {
#   message("All VAR segment models are stable (max eigenvalue modulus < 1).")
# }

}
# ------------------------------------------------------------------------------
# 14. Node importance measures
# ------------------------------------------------------------------------------

node_imp <- compute_node_importance(G, N_BRK)

## finding the highest outdegree fr each segment
outdegree_imp <- rbind(outdegree_imp, as.matrix(node_imp %>% group_by( Section) %>% 
  filter(Outdegree == max(Outdegree)) %>%
  select(Section, City, Outdegree) %>% mutate( i = i)))

##finding the highest indegree for each segment
indegree_imp <- rbind(indegree_imp, as.matrix(node_imp %>% group_by(Section) %>% 
  filter(Indegree == max(Indegree)  )%>%
  select(Section, City, Indegree)%>% mutate( i = i) ))

##finding the out strength for each segment
outstrength_imp <- rbind(outstrength_imp, as.matrix(node_imp %>% group_by(Section) %>% 
  filter(Strength_out ==max(Strength_out))%>%
  select(Section, City, Strength_out ) %>% mutate( i = i)))

## finding the betweenness value for all segments
betweenness_imp <- rbind(betweenness_imp, as.matrix(node_imp %>% group_by(Section) %>% 
  filter(Betweenness == max(Betweenness) & Betweenness > 0 )%>%
  select(Section, City, Betweenness)%>% mutate( i = i)))

}

################## reading the data from ######################
fp <- "C:/Users/gauph/Box/ChangepointOP/"
nodeImpOP <- readRDS(paste0(fp, "nodeImpOP_imp1.rds"))
indegree_imp <- nodeImpOP[[1]]
outdegree_imp <- nodeImpOP[[2]]
outstrength_imp <- nodeImpOP[[3]]
betweenness_imp <- nodeImpOP[[4]]

for(i in 2:49){
  nodeImpOP <- readRDS(paste0(fp, "nodeImpOP_imp",i,".rds"))
  indegree_imp <- rbind(indegree_imp, nodeImpOP[[1]])
  outdegree_imp <- rbind(outdegree_imp, nodeImpOP[[2]])
  outstrength_imp <- rbind(outstrength_imp, nodeImpOP[[3]])
  betweenness_imp <- rbind(betweenness_imp, nodeImpOP[[4]])
  
}

indegree_sum <- as.data.frame(indegree_imp) %>%
  group_by(Section, City) %>%
  summarise(
    count = n(),
    avg_indegree = floor(mean(as.numeric(Indegree), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  group_by(Section) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  filter(count >= quantile(count, 0.90, na.rm = TRUE))

outdegree_sum <- as.data.frame(outdegree_imp) %>%
  group_by(Section, City) %>%
  summarise(
    count = n(),
    avg_outdegree = floor(mean(as.numeric(Outdegree), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  group_by(Section) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  filter(count >= quantile(count, 0.90, na.rm = TRUE))

outstrength_sum <- as.data.frame(outstrength_imp) %>%
  group_by(Section, City) %>%
  summarise(
    count = n(),
    avg_strengthout = mean(as.numeric(Strength_out), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Section) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  filter(count >= quantile(count, 0.90, na.rm = TRUE))

betweenness_sum <- as.data.frame(betweenness_imp) %>%
  group_by(Section, City) %>%
  summarise(
    count = n(),
    avg_betweenness = mean(as.numeric(Betweenness), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Section) %>%
  arrange(desc(count), .by_group = TRUE) %>%
  filter(count >= quantile(count, 0.90, na.rm = TRUE))

#nodeImpOP <- list(indegree_sum, outdegree_sum, outstrength_sum,betweenness_sum)
#saveRDS(nodeImpOP, paste0(SAVE_DIR,"nodeImpOP",args_1,".rds"))
#gbg<- readRDS( paste0(SAVE_DIR,"nodeImpOP.rds"))

################## reading the data from ######################
fp <- "C:/Users/gauph/Box/ChangepointOP2/"
nodeImpOP <- readRDS(paste0(fp, "nodeImpOP_imp1.rds"))
indegree_sum <- cbind(nodeImpOP[[1]],1)
outdegree_sum <- cbind(nodeImpOP[[2]],1)
outstrength_sum <- cbind(nodeImpOP[[3]],1)
betweenness_sum <- cbind(nodeImpOP[[4]],1)

G_edges <- matrix(0,nrow=0, ncol =6)#replicate(8, matrix(0,nrow=0, ncol =6))

for(i in 2:130){
  if(file.exists(paste0(fp, "nodeImpOP_imp",i,".rds"))){
  nodeImpOP <- readRDS(paste0(fp, "nodeImpOP_imp",i,".rds"))
  indegree_sum <- rbind(indegree_sum, cbind(nodeImpOP[[1]], i ))
  outdegree_sum <- rbind(outdegree_sum, cbind(nodeImpOP[[2]], i ))
  outstrength_sum <- rbind(outstrength_sum, cbind(nodeImpOP[[3]], i ))
  betweenness_sum <- rbind(betweenness_sum, cbind(nodeImpOP[[4]], i ))
  }
}

for(i in 1:130){
  ## read the graph rds
  if(file.exists(paste0(fp,"GraphList",i,".rds"))){
    G_tmp <- readRDS(paste0(fp,"GraphList",i,".rds")) 
    for(j in 1:20){
      for(k in 1:8){
        if(is.igraph(G_tmp[[j]][[k]])){
          G_edges <- rbind(G_edges, 
                           cbind(as_data_frame(G_tmp[[j]][[k]], 
                                               what = "edges"), 
                                 k, j, i)) 
        }
      }
    }
  }
}

nms <- cbind(as.data.frame(vertex_attr(G_tmp[[j]][[k]])) ,id= 1:27)
G_edges$from <- nms$names[G_edges$from]
G_edges$to <- nms$names[G_edges$to]
G_edges_total <- G_edges %>% 
  group_by(from, to, k) %>% 
  summarise( total  = n(),
             weight  = mean(weight)) %>%
  ungroup() %>% group_by(k) %>%
  filter(total >= quantile(total, 0.75, na.rm = TRUE))

for(i in 1:8){
  print(G_edges_total %>% filter( k == i) %>%
    ggplot()+ geom_histogram(aes(x = total)))
}

G <- list()

for(i in 1:8){
  EL <- G_edges_total %>% filter(k == i)
  G[[i]] <- simplify(igraph::graph_from_data_frame(EL))
  G[[i]] <- set_vertex_attr(G[[i]], "names",index = V(G[[i]]),V(G[[i]])$name )
  plot(G[[i]],edge.arrow.size = 0.5)
}
node_imp <- compute_node_importance(G, 7)

# Summary: maximum out-strength per segment (key driver identification)
max_imp <- node_imp %>%
  group_by(Section) %>%
  summarise(
    Max_Outdegree   = max(Outdegree),
    Max_Indegree    = max(Indegree),
    Max_Betweenness = max(Betweenness),
    Top_Betweenness = City[which.max(Betweenness)],
    Max_Strength_out = max(Strength_out),
    Top_Strength_out = City[which.max(Strength_out)],
    .groups = "drop"
  )

# print(max_imp)

# ------------------------------------------------------------------------------
# 15. Geographic network visualisation
#
#     Networks are plotted using geographic lat/lon coordinates as node layout.
#     Node labels and colours highlight the top-ranked nodes by out-strength,
#     enabling visual comparison of influential locations across segments.
# ------------------------------------------------------------------------------

# Load geographic layout data
grph_dat <- readRDS(GRAPH_RDS)
g_geo    <- grph_dat$WWgrph

layout_base <- data.frame(
  names = make.names(igraph::V(g_geo)$name),
  lat   = as.numeric(igraph::V(g_geo)$lat),
  lon   = as.numeric(igraph::V(g_geo)$lon)
)

# Add locations not present in the base graph
extra_locs <- data.frame(
  names = c("Rock.Creek", "Dallas",   "Sunriver"),
  lat   = c(-122.8808,    -123.3170,  -121.4334),
  lon   = c(45.5554,       44.9193,    43.8694)
)
layout_base <- rbind(layout_base, extra_locs)


# Plot networks for out-strength (primary node importance measure)
# Additional measures can be added to the 'plot_vars' vector
plot_vars <- c("Closeness") #c("Outdegree","Indegree","Betweenness",
              #  "Strength_out", "Strength_in","NeighConnect_in",
               # "NeighConnect_out","H_Index_in", "H_Index_out",
              #  "Coll_Inf_in", "Coll_Inf_out" ,"IVI_in", 
               # "IVI_out", "IVI_all","Closeness","Eigenvector")#
#c("Strength_out","Outdegree","Indegree", "Betweenness")
Ranked_top <- matrix(0, nrow= 0 , ncol = 10)

for (v in plot_vars) {
  # Rank nodes within each segment by the chosen measure
  ranked <- node_imp %>%
    group_by(Section) %>%
    subset(select = c(Section, get(v), City)) %>%
    mutate(
      max_rank = get(v) >= quantile(get(v), 0.95, na.rm = TRUE),#dense_rank(-as.numeric(.data[[v]])),
      min_rank = get(v) < quantile(get(v), 0.95, na.rm = TRUE),#dense_rank( as.numeric(.data[[v]]))
    ) %>%
    mutate(
      label_top    = ifelse(max_rank == T, City, ""),
      label_bottom = ifelse(min_rank == T, City, ""),
      color_top    = ifelse(max_rank == T, "green4", NA),
      color_bottom = ifelse(min_rank == T, "red",    NA)
    )
  
  Ranked_top <- rbind(Ranked_top, 
                      cbind(as.matrix(ranked %>% 
                                        filter(max_rank ==T)),v))
  

  for (i in seq_len(N_BRK + 1)) {
    # Match node order to network vertex names
    name_val <- data.frame(names = make.names(igraph::V(G[[i]])$name))
    layout   <- left_join(name_val, layout_base, by = "names")
    
    seg_ranked <- ranked %>% filter(Section == i)

    igraph::V(G[[i]])$label_top   <- seg_ranked$label_top
    igraph::V(G[[i]])$color_top   <- seg_ranked$color_top

    g_seg  <- igraph::simplify(G[[i]])
    edg_wt <- igraph::E(g_seg)$weight * 10

    out_path <- paste0(SAVE_DIR, v,"_segment_", i, ".png")
    png(out_path, width = 800, height = 800)
    par(mar = c(0.8, 0, 0.8, 0))

    # NOTE: igraph layout expects (x, y) = (longitude, latitude)
    plot.igraph(
      g_seg,
      directed            = TRUE,
      layout              = cbind(as.numeric(layout$lat),
                                  as.numeric(layout$lon)),
      edge.arrow.size     = 0.3,
      edge.width          = edg_wt,
      edge.color          = "gray60",
      vertex.label        = igraph::V(g_seg)$label_top,
      vertex.size         = 5,
      vertex.color        = igraph::V(G[[i]])$color_top,
      vertex.label.dist   = 1,
      vertex.frame.color  = NULL,
      vertex.label.font   = 2
      #main                = sprintf("Segment %d — %s", i, v)
    )

    dev.off()
    message(sprintf("Saved: %s", out_path))
  }
}

Date_Range <- c("2020-08-31 : 2021-07-12", "2021-07-12 : 2022-04-18", 
                "2022-04-18 : 2022-12-19", "2022-12-19 : 2023-09-18", 
                "2023-09-18 : 2024-06-17", "2024-06-17 : 2025-02-17",
                "2025-02-17 : 2025-09-29", "2025-09-29 : 2026-04-20")

colnames(Ranked_top) <- c("Section","Value","City","max_rank","min_rank",
                          "label_top","label_bottom","color_top","color_bottom",
                          "metric" ) 

Ranked_top<- as.data.frame(Ranked_top)

Ranked_top <- Ranked_top %>% 
  select(c(Section, Value, City, metric))

for (i in seq_len(N_BRK + 1)) {
  # Match node order to network vertex names
  name_val <- data.frame(names = make.names(igraph::V(G[[i]])$names))
  layout   <- left_join(name_val, layout_base, by = "names")
  
  g_seg  <- igraph::simplify(G[[i]])
  edg_wt <- igraph::E(g_seg)$weight * 10
  
  #out_path <- paste0(SAVE_DIR, "none_segment_", i, ".png")
  #png(out_path, width = 800, height = 800)
  #par(mar = c(0.8, 0, 0.8, 0))
  
  # NOTE: igraph layout expects (x, y) = (longitude, latitude)
  plot.igraph(
    g_seg,
    directed            = TRUE,
    layout              = cbind(as.numeric(layout$lat),
                                as.numeric(layout$lon)),
    edge.arrow.size     = 0.3,
    #edge.width          = edg_wt,
    edge.color          = "gray60",
    vertex.label        = igraph::V(g_seg)$names,
    vertex.size         = 2,
    vertex.color        = "darkblue",#igraph::V(G[[i]])$color_top,
    vertex.label.dist   = 1.5,
    vertex.frame.color  = NULL,
    vertex.label.font   = 3,
    vertex.label.cex   = 1,
    main                = Date_Range[i]
  )
  
  #dev.off()
  #message(sprintf("Saved: %s", out_path))
}
