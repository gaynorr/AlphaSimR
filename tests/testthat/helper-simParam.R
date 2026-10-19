# Supply a deliberately different global SimParam
# while testing an explicit local one.
# NULL is a real binding here;
# use withoutGlobalSimParam to remove SP.
withGlobalSimParam <- function(simParam, code) {
  hadSP <- exists("SP", envir = .GlobalEnv, inherits = FALSE)
  old <- get0("SP", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (hadSP) {
      assign("SP", old, envir = .GlobalEnv)
    } else if (exists("SP", envir = .GlobalEnv, inherits = FALSE)) {
      rm("SP", envir = .GlobalEnv)
    }
  })
  assign("SP", simParam, envir = .GlobalEnv)
  force(code)
}

withoutGlobalSimParam <- function(code) {
  withGlobalSimParam(NULL, {
    rm("SP", envir = .GlobalEnv)
    force(code)
  })
}

# This specialization catches fallback on metadata paths only. For other paths,
# use withGlobalSimParam with conflicting settings and assert the local result.
withWrongGlobalMiscSP <- function(sp, code) {
  wrong <- sp$clone()
  fail <- function(...) stop("Used the wrong global SP")
  wrong$subsetMisc <- wrong$subsetMiscPop <- fail
  wrong$mergeMisc <- wrong$mergeMiscPop <- fail
  withGlobalSimParam(wrong, code)
}
