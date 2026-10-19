context("misc functions (hooks)")

test_that("SUBSET functions preserve order, repetitions, dimensions and source data",{
  f = miscFunctionFixture()
  original = f$pop
  indices = list(c(4L, 2L, 2L), -2L, c(TRUE, FALSE), integer(),
                 f$pop@id[c(3, 1)])
  for(i in indices){
    positions = if(is.character(i)) match(i, original@id) else i
    out = original[i, simParam=f$sp]
    expect_identical(out@misc$value, original@misc$value[positions])
    expect_identical(out@misc$mat, original@misc$mat[positions,,drop=FALSE])
    expect_identical(out@miscPop$total, sum(original@misc$value[positions]))
    expect_identical(out@miscPop$ids, out@id)
    expect_true(validObject(out))
  }
  expect_identical(original[, simParam=f$sp], original)
  expect_identical(f$pop, original)
  expect_error(original["missing"], "invalid individuals")
  selected = selectInd(original, nInd=2, use=function(x, ...) x@misc$value,
                       simParam=f$sp)
  expect_identical(selected@miscPop$total, 14)
  expect_identical(selected@miscPop$ids, selected@id)
})

test_that("MERGE functions work through c, nested MultiPop, and empty inputs",{
  f = miscFunctionFixture()
  pop = f$pop
  a = pop[c(4, 2), simParam=f$sp]
  b = pop[1, simParam=f$sp]
  expected = pop[c(4, 2, 1), simParam=f$sp]
  expect_identical(c(a, b, simParam=f$sp), expected)
  expect_identical(mergePops(list(a, NULL, b), simParam=f$sp), expected)
  expect_identical(mergePops(newMultiPop(a, newMultiPop(b)), simParam=f$sp),
                   expected)
  expect_identical(mergePops(list(a), simParam=f$sp), a)
  expect_identical(c(pop[integer()], a, simParam=f$sp), a)
  expect_identical(c(pop[integer()], pop[integer()], simParam=f$sp),
                   pop[integer(), simParam=f$sp])
})

test_that("HOOKS are read from explicit or global simulation parameters",{
  f = miscFunctionFixture()
  pop = f$pop
  expect_identical(pop[1, simParam=f$sp]@miscPop$total, 2)
  f$sp$subsetMiscPop = NULL
  expect_identical(pop[1, simParam=f$sp]@miscPop, list())
  # Reading parameters from local environment - not ideal,
  # but would have to pass ... to every function in the call stack to avoid it
  scale = 3
  f$sp$subsetMisc = function(pop, i, simParam){
    expect_identical(simParam, f$sp)
    list(value=pop@misc$value[i]*scale)
  }
  expect_identical(pop[2, simParam=f$sp]@misc$value, 12)
  scale = 1
  old = get0("SP", envir=.GlobalEnv, inherits=FALSE)
  hadSP = exists("SP", envir=.GlobalEnv, inherits=FALSE)
  on.exit(if(hadSP) assign("SP", old, envir=.GlobalEnv) else
    if(exists("SP", envir=.GlobalEnv, inherits=FALSE)) rm("SP", envir=.GlobalEnv))
  assign("SP", f$sp, envir=.GlobalEnv)
  expect_identical(pop[2]@misc$value, 4)
  expect_identical(c(pop, pop)@miscPop$total, 40)
  expect_identical(mergePops(list(pop, pop))@miscPop$total, 40)
  other = SimParam$new(f$founder)
  expect_identical(pop[2, simParam=other]@miscPop, list())
  assign("SP", 1, envir=.GlobalEnv)
  expect_identical(pop[2]@miscPop, list())
  rm("SP", envir=.GlobalEnv)
  expect_identical(pop[2]@miscPop, list())
  expect_error(pop[1, simParam=1], "simParam must be")
  f$sp$subsetMisc = 1
  expect_error(pop[1, simParam=f$sp], "must be a function or NULL")
  f$sp$mergeMisc = 1
  expect_error(c(pop, pop, simParam=f$sp), "must be a function or NULL")
})

test_that("HOOK results and failures are checked",{
  f = miscFunctionFixture()
  pop = f$pop
  for(value in list(NULL, 1, "invalid")){
    f$sp$subsetMisc = function(pop, i, simParam) value
    f$sp$mergeMisc = function(popList, simParam) value
    expect_error(pop[1, simParam=f$sp], "must return a list")
    expect_error(c(pop, pop, simParam=f$sp), "must return a list")
  }
  f$sp$subsetMisc = function(pop, i, simParam){
    list(wrong=1:10)
  }
  expect_error(pop[1, simParam=f$sp], "nInd")
  f$sp$mergeMisc = function(popList, simParam){
    list(wrong=1:10)
  }
  expect_error(c(pop, pop, simParam=f$sp), "nInd")
  f$sp$subsetMisc = function(pop, i, simParam) stop("custom failure")
  expect_error(pop[1, simParam=f$sp], "custom failure")
  expect_identical(pop@misc, f$pop@misc)
})

test_that("CUSTOM merging reconciles schemas and uses captured settings",{
  f = miscFunctionFixture()
  a = f$pop[1:2]
  b = f$pop[3:4]
  b@misc = list(other=b@misc$value)
  # Reading parameters from local environment - not ideal,
  # but would have to pass ... to every function in the call stack to avoid it
  scale = 1
  f$sp$mergeMisc = function(popList, simParam){
    values = unlist(lapply(popList, function(x) x@misc[[1]]), use.names=FALSE)
    list(value=values*scale)
  }
  # NA since we expect no warning.
  expect_warning({
    out = c(a, b, simParam=f$sp)
  }, NA)
  expect_identical(out@misc$value, f$pop@misc$value)
  expect_identical(out@miscPop$total, 20)
  scale = 2
  expect_identical(mergePops(list(a, b), simParam=f$sp)@misc$value,
                   f$pop@misc$value*2)
  scale = 1
  mp = newMultiPop(a, newMultiPop(b))
  expect_identical(mergeMultiPops(mp, simParam=f$sp)@miscPop$total, 20)
  expect_identical(mergePops(mp, simParam=f$sp)@miscPop$total, 20)
  expect_identical(subsetMisc(f$pop, f$pop@id[2])$value, 4)
})

test_that("SLOTS have independent functions receiving unchanged source populations",{
  f = miscFunctionFixture()
  sp = SimParam$new(f$founder)
  pop = f$pop
  sp$subsetMisc = function(pop, i, simParam){
    list(value=pop@misc$value[i]*10)
  }
  sp$mergeMisc = function(popList, simParam){
    list(value=unlist(lapply(popList, function(pop) pop@misc$value))*10)
  }
  expect_identical(pop[2, simParam=sp]@misc$value, 40)
  expect_identical(pop[2, simParam=sp]@miscPop, list())
  expect_identical(c(pop, pop, simParam=sp)@misc$value, rep(pop@misc$value, 2)*10)
  expect_identical(c(pop, pop, simParam=sp)@miscPop, list())
  scale = 2
  sp$subsetMiscPop = function(pop, i, simParam){
    expect_identical(simParam, sp)
    list(total=sum(pop@misc$value[i])*scale, original=pop@miscPop$total)
  }
  sp$mergeMiscPop = function(popList, simParam){
    expect_identical(simParam, sp)
    list(total=sum(vapply(popList, function(pop) sum(pop@misc$value),
                          numeric(1)))*scale)
  }
  out = pop[2, simParam=sp]
  expect_identical(out@misc$value, 40)
  expect_identical(out@miscPop, list(total=8, original=20))
  out = mergePops(list(pop, pop), simParam=sp)
  expect_identical(out@misc$value, rep(pop@misc$value, 2)*10)
  expect_identical(out@miscPop$total, 80)
  sp$subsetMisc = NULL
  sp$mergeMisc = NULL
  expect_identical(pop[2, simParam=sp]@misc, subsetMisc(pop, 2))
  expect_identical(c(pop, pop, simParam=sp)@misc, mergeMisc(list(pop, pop)))
  expect_identical(subsetMiscPop(pop, 2), list())
  expect_identical(mergeMiscPop(list(pop, pop)), list())
  expect_identical(pop, f$pop)
})

test_that("POPULATION functions validate configurations and list results independently",{
  f = miscFunctionFixture()
  for(name in c("subsetMiscPop", "mergeMiscPop")){
    for(value in list(NULL, 1, "invalid")){
      f$sp[[name]] = function(...) value
      if(name=="subsetMiscPop"){
        expect_error(f$pop[1, simParam=f$sp], "subsetMiscPop must return a list")
      }else{
        expect_error(c(f$pop, f$pop, simParam=f$sp),
                     "mergeMiscPop must return a list")
      }
    }
    f$sp[[name]] = 1
    if(name=="subsetMiscPop"){
      expect_error(f$pop[1, simParam=f$sp], "subsetMiscPop must be a function")
    }else{
      expect_error(c(f$pop, f$pop, simParam=f$sp),
                   "mergeMiscPop must be a function")
    }
    f$sp[[name]] = function(...) list()
  }
})

test_that("NESTED misc populations receive the same simulation parameters",{
  f = miscFunctionFixture()
  pop = f$pop
  pop@misc$nested = f$pop
  pop@misc$containers = newMultiPop(f$pop, f$pop, f$pop, f$pop)
  withWrongGlobalMiscSP(f$sp, {
    out = pop[c(4, 2), simParam=f$sp]
    expect_identical(out@misc$nested@miscPop$total, 12)
    expect_identical(out@misc$nested@id, pop@id[c(4, 2)])
    expect_identical(out@misc$containers[[1]], f$pop)
    out = c(pop, pop, simParam=f$sp)
    expect_identical(out@misc$nested@miscPop$total, 40)
    expect_identical(out@misc$nested@nInd, 8L)
    expect_length(out@misc$containers, 8)
    expect_identical(out@misc$containers[[8]], f$pop)
  })
})

test_that("MULTIPOP containers accept simParam without running individual functions",{
  f = miscFunctionFixture()
  pop = f$pop
  mp = newMultiPop(pop, newMultiPop(pop))
  withWrongGlobalMiscSP(f$sp, {
    # This SP deliberately has functions that fail; container operations need none.
    wrong = get("SP", envir=.GlobalEnv)
    expect_identical(c(mp, pop, simParam=wrong), c(mp, pop))
    expect_identical(c(pop, mp, simParam=wrong), c(pop, mp))
    expect_identical(c(mp, mp, simParam=wrong), c(mp, mp))
    expect_identical(mp[1], newMultiPop(pop))
    expect_identical(mp[[1]], pop)
    merged = mergeMultiPops(pop, mp, simParam=f$sp)
    expect_identical(merged@miscPop$total, 60)
    nested = mergeMultiPops(mp, level=1, simParam=f$sp)
    expect_identical(nested[[2]]@miscPop$total, 20)
  })
})

test_that("DEFAULT misc functions initialized on each simulation",{
  founder = quickHaplo(nInd=4, nChr=1, segSites=10)
  sp = SimParam$new(founder)
  sp$nThreads = 1L
  defaults = list(subsetMisc=subsetMisc, subsetMiscPop=subsetMiscPop,
                  mergeMisc=mergeMisc, mergeMiscPop=mergeMiscPop)
  for(name in names(defaults)) expect_identical(sp[[name]], defaults[[name]])
  pop = newPop(founder, simParam=sp)
  pop@misc = list(value=1:4)
  pop@miscPop = list(label="original")
  withoutGlobalSimParam({
    expect_identical(pop[2:3, simParam=sp], pop[2:3])
    expect_identical(c(pop, pop, simParam=sp), c(pop, pop))
    expect_identical(pop[2:3, simParam=sp]@misc$value, 2:3)
    expect_identical(c(pop, pop, simParam=sp)@misc$value, rep(1:4, 2))
    expect_identical(pop[2:3, simParam=sp]@miscPop, list())
    expect_identical(c(pop, pop, simParam=sp)@miscPop, list())
  })
})
