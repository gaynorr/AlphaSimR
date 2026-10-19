# Run from the package root: Rscript dev/audit_simParam_example.R
# Demonstrates the output and checks detection and CLI exit statuses.
# These small source examples are parsed by the auditor, never executed.
# No installed AlphaSimR package or compiled code is needed.

simParamAuditExamples = function(){
  c(
    'mergePops = function(popList, simParam=NULL) NULL',
    'isPop = function(x) TRUE',
    'knownHelpers = function(xs, simParam=NULL){',
    '  localHelper = function(x) x',
    '  lapply(xs, isPop)',
    '  lapply(xs, localHelper)',
    '  lapply(xs, length)',
    '  lapply(xs, base::is.null)',
    '  do.call("rbind", xs)',
    '}',
    'unknownCallback = function(xs, FUN, simParam=NULL){',
    '  lapply(xs, FUN)',
    '}',
    'shadowedBase = function(xs, length, simParam=NULL){',
    '  lapply(xs, length)',
    '}',
    'noContext = function(popList){',
    '  mergePops(popList)',
    '}',
    'indirectCombine = function(xs, simParam=NULL){',
    '  do.call("c", xs)',
    '}',
    'goodIndirectCombine = function(a, b, simParam=NULL){',
    '  do.call("c", list(a, b, simParam=simParam))',
    '}',
    'badMerge = function(popList, simParam=NULL){',
    '  mergePops(popList)',
    '}',
    'goodMerge = function(popList, simParam=NULL){',
    '  mergePops(popList, simParam=simParam)',
    '}',
    'goodPositional = function(popList, simParam=NULL){',
    '  mergePops(popList, simParam)',
    '}',
    'badDoCall = function(popList, simParam=NULL){',
    '  do.call(mergePops, list(popList))',
    '}',
    'goodDoCall = function(popList, simParam=NULL){',
    '  do.call(mergePops, list(popList, simParam=simParam))',
    '}',
    'differentContext = function(popList, simParam=NULL){',
    '  mergePops(popList, simParam=otherSP)',
    '}',
    'mapped = function(groups, simParam=NULL){',
    '  lapply(groups, mergePops)',
    '}',
    'goodMapped = function(groups, simParam=NULL){',
    '  lapply(groups, mergePops, simParam=simParam)',
    '}',
    'badMap = function(groups, simParam=NULL){',
    '  Map(mergePops, groups)',
    '}',
    'goodMap = function(groups, params, simParam=NULL){',
    '  Map(mergePops, groups, params)',
    '}',
    'goodMapply = function(groups, simParam=NULL){',
    '  mapply(mergePops, groups, MoreArgs=list(simParam=simParam))',
    '}',
    'goodMappedPositional = function(groups, simParam=NULL){',
    '  lapply(groups, mergePops, simParam)',
    '}',
    'badSubset = function(pop, i, simParam=NULL){',
    '  pop[i]',
    '}',
    'indexOnly = function(i, simParam=NULL){',
    '  i[1]',
    '}',
    'badCombine = function(a, b, simParam=NULL){',
    '  c(a, b)',
    '}',
    'goodCombine = function(a, b, simParam=NULL){',
    '  c(a, b, simParam=simParam)',
    '}',
    'goodSubset = function(pop, i, simParam=NULL){',
    '  pop[i, simParam=simParam]',
    '}'
  )
}

runSimParamAuditExample = function(){
  source("dev/audit_simParam.R", local=TRUE)
  path = tempfile("simParam-audit-test-")
  dir.create(file.path(path, "R"), recursive=TRUE)
  on.exit(unlink(path, recursive=TRUE))
  fixture = file.path(path, "R", "example.R")
  writeLines(simParamAuditExamples(), fixture)
  result = audit_simParam(path)
  missing = startsWith(result$issue, "missing:")
  stopifnot(setequal(result$caller[missing], c("badMerge", "badDoCall", "mapped", "noContext", "badMap")),
            all(missing[seq_len(sum(missing))]),
            all(missing),
            !any(startsWith(result$caller, "good")))
  operators = audit_simParam(path, includeOperators=TRUE)
  stopifnot(setequal(setdiff(operators$caller, result$caller),
                    c("badSubset", "indexOnly", "badCombine", "indirectCombine")))

  show = function(result){
    print(result[, c("caller", "issue", "call")], row.names=FALSE)
  }
  cat("Default report: calls that omit an accepted simParam argument.\n")
  show(result)
  cat("\nWith includeOperators=TRUE (CLI: --operators):\n")
  show(operators)
  cat("\n--strict fails for badMerge, badDoCall, mapped, badMap, and noContext.\n",
      "A caller without simParam is still reported when its callee accepts it.\n",
      "Explicit values, known helpers without simParam, and unknown callbacks are omitted.\n",
      "Operator reviews are optional and never fail strict mode.\n",
      "An omission can be intentional: this is not proof that the call is wrong.\n\n",
      sep="")

  # Exercise the actual CLI exit status, not a duplicate of its condition.
  cli = function(flags=character()){
    suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c(shQuote("dev/audit_simParam.R"), shQuote(path), flags),
      stdout=file.path(path, "stdout.txt"), stderr=file.path(path, "stderr.txt")))
  }
  stopifnot(cli()==0L, cli("--strict")==1L,
            cli(c("--strict", "--operators"))==1L)

  # A report containing only review items must not fail strict mode.
  writeLines(c(
    'mergePops = function(popList, simParam=NULL) NULL',
    'reviewOnly = function(popList, i, simParam=NULL){',
    '  mergePops(popList, simParam=otherSP)',
    '  i[1]',
    '}'
  ), fixture)
  stopifnot(cli("--strict")==0L,
            cli(c("--strict", "--operators"))==0L)

  # A caller without simParam still omits an accepted callee argument.
  writeLines(c(
    'mergePops = function(popList, simParam=NULL) NULL',
    'noContext = function(popList) mergePops(popList)'
  ), fixture)
  stopifnot(nrow(audit_simParam(path))==1L, cli("--strict")==1L)

  # A correctly forwarded call produces an empty report.
  writeLines(c(
    'mergePops = function(popList, simParam=NULL) NULL',
    'caller = function(popList, simParam=NULL){',
    '  mergePops(popList, simParam=simParam)',
    '}'
  ), fixture)
  stopifnot(nrow(audit_simParam(path))==0L, cli("--strict")==0L,
            cli("--unknown")!=0L)
  cat("All simParam auditor tests passed.\n")
}

if(sys.nframe()==0L) runSimParamAuditExample()
