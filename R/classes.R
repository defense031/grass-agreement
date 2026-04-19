# S3 constructors for grass objects.

new_grass_metrics <- function(values, n, table, positive_level,
                              n_dropped = 0L, call = NULL) {
  structure(
    list(values = values,
         n = n,
         table = table,
         positive_level = positive_level,
         n_dropped = n_dropped,
         call = call),
    class = "grass_metrics"
  )
}

new_grass_reference <- function(prevalence, quality, reference,
                                call = NULL) {
  structure(
    list(prevalence = prevalence,
         quality = quality,
         reference = reference,
         call = call),
    class = "grass_reference"
  )
}

new_grass_result <- function(metrics, prevalence, prevalence_source,
                             regime, regime_note,
                             spec = NULL,
                             reference = NULL, distance = NULL,
                             call = NULL) {
  structure(
    list(metrics = metrics,
         spec = spec,
         prevalence = prevalence,
         prevalence_source = prevalence_source,
         regime = regime,
         regime_note = regime_note,
         reference = reference,
         distance = distance,
         call = call),
    class = "grass_result"
  )
}
