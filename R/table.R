# Build a 2x2 count table from aligned 0/1 integer vectors, forcing all four
# cells to be present (0-count cells included).

build_table <- function(r1, r2) {
  if (length(r1) == 0 || length(r2) == 0) {
    stop("Cannot build 2x2 table from empty rating vectors.", call. = FALSE)
  }
  if (length(r1) != length(r2)) {
    stop("Rating vectors must have equal length.", call. = FALSE)
  }
  if (any(!r1 %in% c(0L, 1L)) || any(!r2 %in% c(0L, 1L))) {
    stop("Rating vectors must contain only 0 and 1 after coercion.", call. = FALSE)
  }
  tab <- table(factor(r1, levels = c(0, 1)),
               factor(r2, levels = c(0, 1)))
  # Convert to plain integer matrix for downstream speed.
  m <- matrix(as.integer(tab), nrow = 2, ncol = 2,
              dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
  m
}
