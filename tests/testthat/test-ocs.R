context("OCS")

smallOcsPop = function(sexes="no"){
  set.seed(6023)
  f = quickHaplo(6,1,20)
  sp = SimParam$new(f)
  sp$nThreads = 1L
  sp$addSnpChip(10)
  sp$setSexes(sexes)
  list(pop=newPop(f,simParam=sp),sp=sp)
}

ocsBoundaryCall = function(s,...){
  suppressMessages(selectOC(s$pop,simParam=s$sp,returnPop=FALSE,
    trait=identity,use=function(pop) seq_len(pop@nInd),
    inbrTarget=1,inbrType="absolute",...))
}

test_that("OCS lower-tail selection and global parameters match explicit calls", {
  s = smallOcsPop()
  expect_equal(sort(ocsBoundaryCall(s,nInd=2,selectTop=FALSE)),1:2)
  expected = ocsBoundaryCall(s,nInd=2)
  hadSP = exists("SP",envir=.GlobalEnv,inherits=FALSE)
  if(hadSP) oldSP = get("SP",envir=.GlobalEnv)
  on.exit({if(hadSP) assign("SP",oldSP,envir=.GlobalEnv)
           else rm("SP",envir=.GlobalEnv)},add=TRUE)
  assign("SP",s$sp,envir=.GlobalEnv)
  actual = suppressMessages(selectOC(s$pop,nInd=2,returnPop=FALSE,
    trait=identity,use=function(pop) seq_len(pop@nInd),inbrTarget=1,
    inbrType="absolute"))
  expect_equal(actual,expected)
})

test_that("OCS rejects ambiguous responses and missing or invalid counts", {
  s = smallOcsPop()
  expect_error(suppressMessages(selectOC(s$pop,nInd=2,simParam=s$sp,
    trait=identity,use=function(pop) matrix(1,pop@nInd,2))),"single column")
  expect_error(ocsBoundaryCall(s,nInd=0),"nInd must")
  expect_error(ocsBoundaryCall(s),"specify nInd")
  s = smallOcsPop("yes_sys")
  expect_error(ocsBoundaryCall(s),"specify nFemale")
  expect_error(ocsBoundaryCall(s,nFemale=1),"specify nMale")
  expect_error(ocsBoundaryCall(s,nFemale=0,nMale=1),"nFemale must")
  expect_error(ocsBoundaryCall(s,nFemale=1,nMale=0),"nMale must")
})

test_that("OCS caps requests at the available candidates", {
  s = smallOcsPop()
  expect_warning(take <- ocsBoundaryCall(s,nInd=7),"smaller than nInd")
  expect_equal(sort(take),seq_len(s$pop@nInd))
  s = smallOcsPop("yes_sys")
  expect_warning(take <- ocsBoundaryCall(s,nFemale=7,nMale=1),"smaller than nFemale")
  expect_equal(sum(s$pop@sex[take]=="F"),sum(s$pop@sex=="F"))
  expect_equal(sum(s$pop@sex[take]=="M"),1)
  expect_warning(take <- ocsBoundaryCall(s,nFemale=1,nMale=7),"smaller than nMale")
  expect_equal(sum(s$pop@sex[take]=="F"),1)
  expect_equal(sum(s$pop@sex[take]=="M"),sum(s$pop@sex=="M"))
})

test_that("OCS handles constant merit and constant genotype penalty", {
  s = smallOcsPop()
  # Identical copies have no spread in the genotype penalty or merit.
  pop = s$pop[rep(1,4)]
  take = suppressMessages(selectOC(pop,nInd=2,returnPop=FALSE,simParam=s$sp,
    trait=identity,use=function(pop) rep(1,pop@nInd),inbrTarget=1,inbrType="absolute"))
  expect_length(take,2)
  expect_equal(length(unique(take)),2L)
  expect_true(all(take %in% seq_len(pop@nInd)))
})
