# 06_threshold_derivation.R — Context-aware threshold derivation
#
# Three strategies:
#   1. Decision-utility mapping: threshold functions T(prevalence, cost_ratio)
#   2. ROC-style optimization: cut points maximizing Youden's J per stratum
#   3. Empirical clustering: k-means on coefficient vectors

# =====================================================================
# Strategy 1: Decision-utility mapping
# =====================================================================
# For each scenario, we know ground-truth quality (from Se/Sp) and observed
# metrics (from simulation). We learn: at what PABAK (or kappa) value does
# a scenario transition from "acceptable" to "unacceptable" under a given
# cost model?

derive_utility_thresholds <- function(summaries, cost_ratios, loss_threshold = 0.10) {
  # For each cost ratio and decision rule, classify scenarios as acceptable/unacceptable
  # based on whether expected loss < loss_threshold.
  # Then find the PABAK (and kappa) value that best separates the two classes
  # within each prevalence stratum.

  threshold_tables <- list()

  for (cr in cost_ratios) {
    cr_label <- gsub("\\.", "p", as.character(cr))

    for (rule in c("require", "any")) {
      loss_col <- paste0("loss_", rule, "_cr", cr_label)
      if (!loss_col %in% names(summaries)) next

      summaries[, acceptable := get(loss_col) < loss_threshold]

      # Within each prevalence level, find optimal PABAK threshold
      prevs <- sort(unique(summaries$prevalence))

      for (p in prevs) {
        sub <- summaries[prevalence == p, ]
        if (length(unique(sub$acceptable)) < 2) {
          # All same class — record but no threshold derivable
          threshold_tables[[length(threshold_tables) + 1]] <- data.frame(
            cost_ratio = cr,
            rule = rule,
            prevalence = p,
            metric = "PABAK",
            threshold = NA_real_,
            youden_j = NA_real_,
            n_acceptable = sum(sub$acceptable),
            n_unacceptable = sum(!sub$acceptable),
            sensitivity = NA_real_,
            specificity = NA_real_,
            stringsAsFactors = FALSE
          )
          next
        }

        # Find PABAK threshold via Youden's J
        pabak_vals <- sort(unique(sub$PABAK_mean))
        best_j <- -Inf
        best_t <- NA
        best_se <- NA
        best_sp <- NA

        for (t in pabak_vals) {
          pred_accept <- sub$PABAK_mean >= t
          tp <- sum(pred_accept & sub$acceptable)
          fn <- sum(!pred_accept & sub$acceptable)
          fp <- sum(pred_accept & !sub$acceptable)
          tn <- sum(!pred_accept & !sub$acceptable)
          se <- if ((tp + fn) > 0) tp / (tp + fn) else 0
          sp <- if ((tn + fp) > 0) tn / (tn + fp) else 0
          j <- se + sp - 1
          if (j > best_j) {
            best_j <- j
            best_t <- t
            best_se <- se
            best_sp <- sp
          }
        }

        # Same for kappa
        kappa_vals <- sort(unique(sub$kappa_mean))
        best_j_k <- -Inf
        best_t_k <- NA
        best_se_k <- NA
        best_sp_k <- NA

        for (t in kappa_vals) {
          pred_accept <- sub$kappa_mean >= t
          tp <- sum(pred_accept & sub$acceptable)
          fn <- sum(!pred_accept & sub$acceptable)
          fp <- sum(pred_accept & !sub$acceptable)
          tn <- sum(!pred_accept & !sub$acceptable)
          se <- if ((tp + fn) > 0) tp / (tp + fn) else 0
          sp <- if ((tn + fp) > 0) tn / (tn + fp) else 0
          j <- se + sp - 1
          if (j > best_j_k) {
            best_j_k <- j
            best_t_k <- t
            best_se_k <- se
            best_sp_k <- sp
          }
        }

        threshold_tables[[length(threshold_tables) + 1]] <- data.frame(
          cost_ratio = cr,
          rule = rule,
          prevalence = p,
          metric = "PABAK",
          threshold = best_t,
          youden_j = best_j,
          n_acceptable = sum(sub$acceptable),
          n_unacceptable = sum(!sub$acceptable),
          sensitivity = best_se,
          specificity = best_sp,
          stringsAsFactors = FALSE
        )

        threshold_tables[[length(threshold_tables) + 1]] <- data.frame(
          cost_ratio = cr,
          rule = rule,
          prevalence = p,
          metric = "kappa",
          threshold = best_t_k,
          youden_j = best_j_k,
          n_acceptable = sum(sub$acceptable),
          n_unacceptable = sum(!sub$acceptable),
          sensitivity = best_se_k,
          specificity = best_sp_k,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  data.table::rbindlist(threshold_tables)
}


# =====================================================================
# Strategy 2: ROC-style thresholds for quality band classification
# =====================================================================
# Treat quality_band as the ground truth. For each prevalence stratum,
# compute ROC-optimal thresholds for classifying "high" vs "not high"
# and "low" vs "not low".

derive_roc_thresholds <- function(summaries) {
  prevs <- sort(unique(summaries$prevalence))
  roc_results <- list()

  for (p in prevs) {
    sub <- summaries[prevalence == p, ]

    for (target in c("high", "low")) {
      sub[, target_class := quality_band == target]

      if (length(unique(sub$target_class)) < 2) {
        roc_results[[length(roc_results) + 1]] <- data.frame(
          prevalence = p, target = target,
          metric = "PABAK", threshold = NA, youden_j = NA,
          sensitivity = NA, specificity = NA,
          stringsAsFactors = FALSE
        )
        next
      }

      for (metric_name in c("PABAK_mean", "kappa_mean")) {
        vals <- sort(unique(sub[[metric_name]]))
        best_j <- -Inf; best_t <- NA; best_se <- NA; best_sp <- NA

        for (t in vals) {
          if (target == "high") {
            pred <- sub[[metric_name]] >= t
          } else {
            pred <- sub[[metric_name]] <= t
          }
          tp <- sum(pred & sub$target_class)
          fn <- sum(!pred & sub$target_class)
          fp <- sum(pred & !sub$target_class)
          tn <- sum(!pred & !sub$target_class)
          se <- if ((tp + fn) > 0) tp / (tp + fn) else 0
          sp <- if ((tn + fp) > 0) tn / (tn + fp) else 0
          j <- se + sp - 1
          if (j > best_j) {
            best_j <- j; best_t <- t; best_se <- se; best_sp <- sp
          }
        }

        roc_results[[length(roc_results) + 1]] <- data.frame(
          prevalence = p, target = target,
          metric = gsub("_mean$", "", metric_name),
          threshold = best_t, youden_j = best_j,
          sensitivity = best_se, specificity = best_sp,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  data.table::rbindlist(roc_results)
}


# =====================================================================
# Strategy 3: Empirical clustering
# =====================================================================
# Within each prevalence stratum, cluster scenario-level metric vectors
# and assess correspondence with ground-truth quality bands.

derive_clustering_thresholds <- function(summaries, k = 3,
                                          method = c("kmeans", "hclust")) {
  method <- match.arg(method)
  prevs <- sort(unique(summaries$prevalence))
  cluster_results <- list()

  feature_cols <- c("P0_mean", "kappa_mean", "PABAK_mean",
                    "pos_agree_mean", "neg_agree_mean")

  for (p in prevs) {
    sub <- summaries[prevalence == p, ]

    # Extract and standardize features
    feat <- as.data.frame(sub[, ..feature_cols])

    # Handle NAs: replace with column median
    for (col in names(feat)) {
      na_idx <- is.na(feat[[col]])
      if (any(na_idx)) feat[na_idx, col] <- median(feat[[col]], na.rm = TRUE)
    }

    feat_scaled <- scale(feat)

    # Remove zero-variance columns
    keep <- apply(feat_scaled, 2, function(x) sd(x, na.rm = TRUE) > 0)
    feat_scaled <- feat_scaled[, keep, drop = FALSE]

    if (ncol(feat_scaled) < 2 || nrow(feat_scaled) < k) {
      cluster_results[[length(cluster_results) + 1]] <- data.frame(
        prevalence = p, method = method,
        ari = NA, purity = NA, n_scenarios = nrow(sub),
        stringsAsFactors = FALSE
      )
      next
    }

    if (method == "kmeans") {
      set.seed(42)
      cl <- kmeans(feat_scaled, centers = k, nstart = 25)
      labels <- cl$cluster
    } else {
      d <- dist(feat_scaled)
      hc <- hclust(d, method = "ward.D2")
      labels <- cutree(hc, k = k)
    }

    # Map quality bands to integers for comparison
    truth <- as.integer(factor(sub$quality_band, levels = c("low", "medium", "high")))

    # Purity: fraction of each cluster's majority class
    purity <- 0
    for (ci in seq_len(k)) {
      members <- truth[labels == ci]
      if (length(members) > 0) {
        purity <- purity + max(table(members))
      }
    }
    purity <- purity / length(truth)

    # Adjusted Rand Index (simple implementation)
    ari <- compute_ari(labels, truth)

    cluster_results[[length(cluster_results) + 1]] <- data.frame(
      prevalence = p, method = method,
      ari = ari, purity = purity, n_scenarios = nrow(sub),
      stringsAsFactors = FALSE
    )
  }

  data.table::rbindlist(cluster_results)
}

# Simple ARI computation
compute_ari <- function(labels1, labels2) {
  tab <- table(labels1, labels2)
  sum_comb <- function(x) sum(choose(x, 2))

  n <- length(labels1)
  sum_ij <- sum(choose(tab, 2))
  sum_a <- sum_comb(rowSums(tab))
  sum_b <- sum_comb(colSums(tab))
  total_comb <- choose(n, 2)

  expected <- sum_a * sum_b / total_comb
  max_index <- 0.5 * (sum_a + sum_b)

  if (max_index == expected) return(1)
  (sum_ij - expected) / (max_index - expected)
}
