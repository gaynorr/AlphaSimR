#!/usr/bin/env Rscript

findRcppTskitLib <- function() {
  libdir <- system.file("libs", package = "RcppTskit")
  if (!nzchar(libdir)) {
    stop("Unable to locate the RcppTskit shared library directory")
  }

  libdirs <- libdir
  if (.Platform$OS.type == "windows") {
    r_arch <- sub("^/", "", .Platform$r_arch)
    if (nzchar(r_arch)) {
      libdirs <- c(libdirs, file.path(libdir, r_arch))
    } else {
      arch_dirs <- c("x64", "i386")
      arch_dirs <- arch_dirs[dir.exists(file.path(libdir, arch_dirs))]
      libdirs <- c(libdirs, file.path(libdir, arch_dirs))
    }
    libdirs <- unique(libdirs)
  }

  candidates <- if (.Platform$OS.type == "unix") {
    c("RcppTskit.so", "RcppTskit.dylib")
  } else {
    c("RcppTskit.dll.a", "RcppTskit.lib", "RcppTskit.dll")
  }
  libpaths <- unlist(lapply(libdirs, file.path, candidates), use.names = FALSE)
  libfile <- libpaths[file.exists(libpaths)][1L]
  if (is.na(libfile) || !nzchar(libfile)) {
    stop("Unable to locate the RcppTskit shared library")
  }

  if (.Platform$OS.type == "unix") {
    sprintf("-Wl,-rpath,%s %s", shQuote(dirname(libfile)), shQuote(libfile))
  } else {
    shQuote(libfile)
  }
}

renderMakevars <- function(template, output) {
  lines <- readLines(template)
  lines <- gsub(
    "@RCPPTSKIT_LIB@",
    findRcppTskitLib(),
    lines,
    fixed = TRUE
  )
  writeLines(lines, output)
}

if (.Platform$OS.type == "unix") {
  renderMakevars("src/Makevars.in", "src/Makevars")
} else {
  renderMakevars("src/Makevars.win.in", "src/Makevars.win")
}
