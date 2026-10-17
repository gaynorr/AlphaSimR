# Rscript DemographicModelsBenchmark.R [output.csv]
#
# Each invocation starts a fresh simulation; no input data files are needed.
# The optional CSV path defaults to segregating-sites-benchmark.csv in the
# working directory. One summary CSV with site counts and timings is written
# beside it, replacing prior output. Timing statistics are in seconds.
# Each model uses 100 replicates at 10, 100, and 1,000 individuals, and
# one replicate at 10,000 individuals. Each replicate runs one chromosome
# on one thread, recording both its segregating-site count and elapsed time.

library(AlphaSimR)

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else "segregating-sites-benchmark.csv"
folder <- dirname(output)
dir.create(folder, recursive = TRUE, showWarnings = FALSE)

models <- c("GENERIC", "CATTLE", "WHEAT", "MAIZE")
sampleSizes <- c(10L, 100L, 1000L, 10000L)
plan <- rbind(
  expand.grid(
    model = models,
    nInd = sampleSizes[1:3],
    replicate = seq_len(100L),
    stringsAsFactors = FALSE
  ),
  expand.grid(
    model = models,
    nInd = 10000L,
    replicate = 1L,
    stringsAsFactors = FALSE
  )
)
plan$key <- with(plan, paste(model, nInd, replicate, sep = "_"))
results <- data.frame()
for (i in seq_len(nrow(plan))) {
  model <- plan$model[i]
  nInd <- plan$nInd[i]
  replicate <- plan$replicate[i]
  key <- plan$key[i]
  gc()
  cat("START", key, format(Sys.time(), tz = "UTC"), "\n")
  flush.console()
  timing <- system.time({
    pop <- runMacs(
      nInd = nInd,
      nChr = 1L,
      species = model,
      ploidy = 2L,
      inbred = FALSE,
      segSites = NULL,
      nThreads = 1L
    )
  })
  row <- data.frame(
    key = key,
    model = model,
    nInd = nInd,
    replicate = replicate,
    nSites = pop@nLoci[1],
    elapsed = timing[["elapsed"]],
    user = timing[["user.self"]],
    system = timing[["sys.self"]],
    nThreads = 1L,
    version = as.character(packageVersion("AlphaSimR"))
  )
  results <- rbind(results, row)
  write.csv(results, paste0(output, ".tmp"), row.names = FALSE)
  if (!file.rename(paste0(output, ".tmp"), output)) {
    stop("Could not save checkpoint")
  }
  cat("DONE", key, "sites", row$nSites, "elapsed", row$elapsed, "\n")
  flush.console()
  rm(pop)
}

benchmark <- results

summary <- do.call(
  rbind,
  lapply(models, function(model) {
    do.call(
      rbind,
      lapply(sampleSizes, function(n) {
        rows <- benchmark$model == model & benchmark$nInd == n
        x <- benchmark$nSites[rows]
        elapsed <- benchmark$elapsed[rows]
        stopifnot(length(x) == if (n == 10000) 1L else 100L)
        data.frame(
          Model = model,
          Individuals = n,
          Replicates = length(x),
          Mean = mean(x),
          SD = if (length(x) > 1) sd(x) else NA_real_,
          Minimum = if (length(x) > 1) min(x) else NA_real_,
          Maximum = if (length(x) > 1) max(x) else NA_real_,
          MeanTime = mean(elapsed),
          MedianTime = median(elapsed),
          MinTime = min(elapsed),
          MaxTime = max(elapsed)
        )
      })
    )
  })
)
write.csv(
  summary,
  file.path(folder, "segregating-sites-summary.csv"),
  row.names = FALSE
)

# Reproduce the vignette table
display <- summary[c(
  "Model", "Individuals", "Replicates", "Mean", "SD", "Minimum", "Maximum",
  "MeanTime"
)]
for (column in c(
  "Individuals",
  "Replicates",
  "Mean",
  "SD",
  "Minimum",
  "Maximum"
)) {
  values <- display[[column]]
  if (column %in% c("Mean", "SD", "Minimum", "Maximum")) {
    values <- round(values, digits = -2)
  }
  display[[column]] <- ifelse(
    is.na(values),
    "—",
    format(values, big.mark = ",", scientific = FALSE, trim = TRUE)
  )
}
display$MeanTime <- sprintf("%.3g", summary$MeanTime)

print(display)
