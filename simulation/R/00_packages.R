# 00_packages.R — Package loading for PABAK simulation
# All packages required by the simulation pipeline.

required_packages <- c(
  "yaml",        # config parsing
  "future",      # parallel backend
  "furrr",       # future-powered purrr
  "data.table",  # fast rbindlist for inner loop
  "ggplot2",     # visualization
  "patchwork"    # multi-panel figures
)

missing <- required_packages[!vapply(required_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(yaml)
  library(future)
  library(furrr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
