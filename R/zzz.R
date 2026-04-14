.AlphaSimRStartupMessage = function() {
  version = as.character(utils::packageVersion("AlphaSimR"))
  threads = getNumThreads()
  thread_label = if (threads == 1L) "thread" else "threads"

  if (isOpenMPAvailable()) {
    return(sprintf(
      paste0(
        "AlphaSimR %s with OpenMP support; using %d %s by default.\n",
        "See vignette(\"parallelization\", package=\"AlphaSimR\") for usage details."
      ),
      version,
      threads,
      thread_label
    ))
  }

  sprintf(
    paste0(
      "AlphaSimR %s without OpenMP support; running in single-threaded mode.\n",
      "See vignette(\"parallelization\", package=\"AlphaSimR\") for setup details."
    ),
    version
  )
}

.onAttach = function(libname, pkgname) {
  if (!interactive() || isTRUE(getOption("AlphaSimR.quiet"))) {
    return(invisible(NULL))
  }

  packageStartupMessage(.AlphaSimRStartupMessage())
  invisible(NULL)
}
