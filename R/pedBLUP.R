# Pedigree BLUP ----

# Fits the single trait model y = X*b + Z*u + e with u ~ N(0, A*Vu) and
# e ~ N(0, I*Ve), where A is the numerator relationship matrix of the
# phenotyped individuals and all of their ancestors, X has one column per
# level of fixEff, and an individual with several records has one breeding
# value. The mixed model equations are built in C++ from the pedigree, with
# individuals that have no Mendelian sampling written through their basis
# ancestors, and solved with Matrix's sparse Cholesky factorization.
#
# Vu and Ve are estimated by derivative-free REML unless both are given.
# With Ve profiled out, the likelihood depends only on lambda = Ve/Vu, so
# it is maximized over log(lambda) by optimize(). Every step refactors the
# same pattern, so the symbolic analysis is done once and each step costs
# one numerical factorization.
#
# Returns the variance components, the fixed effects by level, and the
# estimated breeding values of every individual in the pruned pedigree with
# their iids. Individuals outside the pruned pedigree have no phenotyped
# descendants, so their BLUP is the average of their parents'.
.fitPedBLUP = function(y, iid, fixEff, ploidy, Vu=NULL, Ve=NULL,
                       simParam){
  y = as.numeric(y)
  if((length(iid)!=length(y)) || (length(fixEff)!=length(y))){
    stop("y, iid and fixEff must have the same length", call.=FALSE)
  }
  if(length(y)==0L){
    stop("the training population has no records", call.=FALSE)
  }
  if(anyNA(y)){
    stop("pedigree BLUP needs a record for every individual in the training population, but the responses contain missing values",
         call.=FALSE)
  }
  if(is.null(Vu)!=is.null(Ve)){
    stop("supply both Vu and Ve, or neither", call.=FALSE)
  }
  if(!is.null(Vu)){
    if((length(Vu)!=1L) || (length(Ve)!=1L) || !is.finite(Vu) ||
       !is.finite(Ve) || (Vu<=0) || (Ve<=0)){
      stop("Vu and Ve must each be a single positive value", call.=FALSE)
    }
  }
  ped = .prunePed(iid, ploidy, simParam)
  fix = factor(fixEff)
  N = length(y)
  nFix = nlevels(fix)
  # REML needs residual degrees of freedom. A fit with given variance
  # components only needs a record in every level of the fixed effects,
  # which the levels have by construction.
  if(is.null(Vu) && (N<=nFix)){
    stop("there must be more records than levels of the fixed effects",
         call.=FALSE)
  }

  # Centering the records within each level of the fixed effects changes
  # only the fixed effects, which absorb the level means, and leaves the
  # breeding values and variance components as they are. It keeps y'Py,
  # computed below as y'y less the solutions times the right-hand side,
  # from losing the variation in the records to cancellation when they
  # share a large offset.
  levelMean = as.vector(tapply(y, fix, mean))
  yc = y - levelMean[as.integer(fix)]
  # Records that are equal within every level leave nothing to divide
  # between the two variances, and the likelihood takes the logarithm of
  # zero. They are tested for equality as given, because centering can
  # leave rounding error, and any difference between two records, however
  # small next to their size, is variation that REML can use.
  if(is.null(Vu) &&
     all(tapply(y, fix, function(v) all(v==v[1])))){
    stop("the records do not vary within the levels of the fixed effects, so the variance components cannot be estimated; supply Vu and Ve",
         call.=FALSE)
  }

  mme = calcPedMME(ped$mother, ped$father, ped$isDH, ped$row,
                   as.integer(fix), yc, .pedF(ped, simParam), ped$ploidy)
  n = mme$nFix + mme$nBasis
  C0 = Matrix::sparseMatrix(i=mme$C$i, j=mme$C$j, x=mme$C$x,
                            dims=c(n, n), symmetric=TRUE)
  A0 = Matrix::sparseMatrix(i=mme$A$i, j=mme$A$j, x=mme$A$x,
                            dims=c(n, n), symmetric=TRUE)
  rhs = c(mme$rhs)

  # The pattern of C0+lambda*A0 is the same for every positive lambda, so
  # it is analyzed once and later steps only refactor. The first
  # factorization is made at the first lambda asked for, the given one
  # when Vu and Ve are supplied, rather than at an arbitrary value, at
  # which very small Mendelian sampling variances can leave the matrix
  # indefinite in floating point.
  L = NULL
  solveAt = function(lambda){
    if(is.null(L)){
      L <<- Matrix::Cholesky(C0+lambda*A0, perm=TRUE, LDL=FALSE,
                             super=FALSE)
    }else{
      L <<- Matrix::update(L, C0+lambda*A0)
    }
    as.vector(Matrix::solve(L, rhs, system="A"))
  }

  boundary = FALSE
  if(is.null(Vu)){
    # -2 times the REML log likelihood with Ve profiled out, less a
    # constant. y'Py is y'y less the solutions times the right-hand side.
    m2ll = function(logLambda){
      sol = solveAt(exp(logLambda))
      ypy = mme$yy-sum(sol*rhs)
      (N-nFix)*log(ypy) - mme$nBasis*logLambda +
        as.numeric(Matrix::determinant(L, sqrt=FALSE)$modulus)
    }
    range = log(c(1e-4, 1e4))
    # The likelihood need not have a single peak in lambda, and a single
    # search over the whole range can settle on a local optimum far from
    # the best. It is first evaluated on a grid of one point per factor of
    # 10, which includes both ends of the range, and the search is then
    # confined to the neighbors of the best grid point. Whichever of the
    # two is lower is kept, as the search never evaluates the ends of its
    # interval.
    grid = seq(range[1], range[2], length.out=9L)
    gridVal = vapply(grid, m2ll, numeric(1))
    k = which.min(gridVal)
    opt = optimize(m2ll, interval=c(grid[max(k-1L, 1L)],
                                    grid[min(k+1L, length(grid))]),
                   tol=1e-6)
    logLambda = if(opt$objective<=gridVal[k]) opt$minimum else grid[k]
    if(min(abs(logLambda-range))<1e-3){
      boundary = TRUE
      warning("the REML estimate of the variance ratio is at the edge of the range searched, so the variance components may not be reliable",
              call.=FALSE)
    }
    lambda = exp(logLambda)
    sol = solveAt(lambda)
    Ve = (mme$yy-sum(sol*rhs))/(N-nFix)
    Vu = Ve/lambda
  }else{
    sol = solveAt(Ve/Vu)
  }

  W = Matrix::sparseMatrix(i=mme$W$i, j=mme$W$j, x=mme$W$x,
                           dims=c(mme$W$nrow, mme$W$ncol))
  ebv = as.vector(W %*% sol[-seq_len(nFix)])
  fixEffEst = sol[seq_len(nFix)] + levelMean
  names(fixEffEst) = levels(fix)

  return(list(Vu=Vu, Ve=Ve, fixEff=fixEffEst, iid=ped$iid, ebv=ebv,
              boundary=boundary))
}

#' @title Pedigree BLUP
#'
#' @description
#' Fits a single trait pedigree BLUP model, also known as an animal model,
#' to a training population, using the pedigree recorded in the
#' \code{\link{SimParam}}. The output is used with \code{\link{setEBV}} to
#' set estimated breeding values or genetic values for any individual in
#' the simulation, including individuals created after the model was
#' fitted. The pedigree must have been tracked from the start of the
#' simulation, see \code{\link[=SimParam-setTrackPed]{setTrackPed}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training
#' population. Several generations can be combined with \code{c()}, and an
#' individual that appears more than once has its records treated as
#' repeated records.
#' @param traits an integer indicating the trait to model, a trait name, or
#' a function of the traits returning a single value. Only one trait can be
#' modeled.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param Vu the additive genetic variance. If supplied with Ve, the
#' variance components are not estimated.
#' @param Ve the error variance. If supplied with Vu, the variance
#' components are not estimated.
#' @param varComp an earlier fit of \code{pedBLUP}, whose Vu and Ve are
#' used in place of estimating the variance components. Supply either
#' varComp or Vu and Ve, not both.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param ... additional arguments if using a function for traits
#'
#' @details
#' The model is \eqn{y = Xb + Zu + e}, where \eqn{X} has one column for
#' each level of the population's fixEff slot, \eqn{u \sim N(0, A
#' \sigma^2_u)} with \eqn{A} the numerator relationship matrix of the
#' phenotyped individuals and all of their ancestors, and \eqn{e \sim N(0,
#' I \sigma^2_e)}. The relationships are those of \code{\link{pedigreeA}},
#' including its treatment of doubled haploids, inbred founders and
#' polyploids. An individual with several records has a single breeding
#' value, and no permanent environmental effect is fitted, so any
#' non-additive genetic variance among repeated records is not separated
#' from the additive variance.
#'
#' The mixed model equations are built from the inverse of \eqn{A}, which
#' is sparse, and solved with the sparse Cholesky factorization of the
#' \pkg{Matrix} package. Individuals with no Mendelian sampling variance,
#' such as doubled haploids of inbred parents and F1s of two inbred
#' parents, make \eqn{A} singular, and are written as exact functions of
#' their parents rather than fitted separately.
#'
#' Unless they are supplied, the variance components are estimated by
#' restricted maximum likelihood (REML) with the derivative-free approach
#' of \insertCite{graser1987Derivative;textual}{AlphaSimR} and
#' \insertCite{meyer1989Restricted;textual}{AlphaSimR}: the error variance
#' is profiled out of the likelihood, which is then maximized over the
#' ratio of the error variance to the additive genetic variance. A warning
#' is given if the estimate is at the edge of the range searched, which
#' covers heritabilities from about 0.0001 to 0.9999.
#'
#' An individual outside the fitted pedigree has no phenotyped descendants
#' in the training population, so its best linear unbiased prediction is
#' the average of its parents', which is what \code{\link{setEBV}} uses.
#'
#' @return Returns an object of \code{\link{PedSol-class}}
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=20, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setTrackPed(TRUE)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Train on two generations
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses=20, nProgeny=5, simParam=SP)
#' ans = pedBLUP(c(pop, pop2), simParam=SP)
#'
#' #Predict a generation that was not in the training population
#' pop3 = randCross(pop2, nCrosses=20, simParam=SP)
#' pop3 = setEBV(pop3, ans, value="bv", simParam=SP)
#' cor(bv(pop3, simParam=SP), ebv(pop3))
#'
#' @export
pedBLUP = function(pop, traits=1, use="pheno", Vu=NULL, Ve=NULL,
                   varComp=NULL, simParam=NULL, ...){
  simParam = getSimParam(simParam)
  .checkPedPop(pop)
  if(length(traits)!=1L){
    stop("pedBLUP models one trait at a time", call.=FALSE)
  }
  if(!is.null(varComp)){
    if(!is.null(Vu) || !is.null(Ve)){
      stop("supply either varComp or Vu and Ve, not both", call.=FALSE)
    }
    if(!is(varComp, "PedSol")){
      stop("varComp must be an earlier fit of pedBLUP", call.=FALSE)
    }
    Vu = c(varComp@Vu)
    Ve = c(varComp@Ve)
  }
  y = as.matrix(getResponse(pop=pop, trait=traits, use=use,
                            simParam=simParam, ...))
  if(ncol(y)!=1L){
    stop("pedBLUP models one trait at a time", call.=FALSE)
  }
  ans = .fitPedBLUP(y[,1], pop@iid, pop@fixEff, pop@ploidy, Vu=Vu, Ve=Ve,
                    simParam=simParam)
  # The mean breeding value of the training records, which centers "bv"
  # estimates on the training population as RRBLUP's do
  center = mean(ans$ebv[match(pop@iid, ans$iid)])
  return(new("PedSol",
             iid=as.integer(ans$iid),
             ebv=ans$ebv,
             fixEff=ans$fixEff,
             center=center,
             trait=convertTraitsToNames(traits, simParam),
             Vu=as.matrix(ans$Vu),
             Ve=as.matrix(ans$Ve),
             model="pedBLUP",
             nPedResets=length(.pedResets(simParam))))
}

# The lastId of every reset of the pedigree, in order. A SimParam saved by
# an older version of AlphaSimR has no record, and asking it for the method
# gives NULL.
.pedResets = function(simParam){
  if(is.null(simParam$pedResets)){
    return(integer(0))
  }
  return(simParam$pedResets())
}

# Estimated breeding values of individuals, by iid, from a PedSol. Those in
# the fitted pedigree take their solution. Everyone else has no phenotyped
# descendants in the training population, so their BLUP is the average of
# their parents', a doubled haploid's being its parent's, with 0 for an
# unknown parent. Parents have smaller iids than their progeny, so a single
# pass in iid order fills in every parent before it is used.
.pedPredict = function(solution, iid, simParam){
  # A reset after the fit that cut the pedigree back below one of its
  # individuals lets new individuals take that iid, and a match by iid
  # would then give them the old individual's estimate
  resets = .pedResets(simParam)
  if(length(resets)>solution@nPedResets){
    since = resets[(solution@nPedResets+1L):length(resets)]
    if(any(since<max(solution@iid))){
      stop("the pedigree has been reset since this model was fitted, discarding some of its individuals, whose iids may now belong to others; fit the model again",
           call.=FALSE)
    }
  }
  iid = as.integer(iid)
  take = match(iid, solution@iid)
  if(!anyNA(take)){
    return(solution@ebv[take])
  }
  if(!simParam$isTrackPed){
    stop("predicting individuals outside the training pedigree needs the trackPed option, see setTrackPed in ?SimParam",
         call.=FALSE)
  }
  ped = simParam$pedigree
  if((nrow(ped)!=simParam$lastId) || (max(solution@iid)>nrow(ped)) ||
     any(is.na(iid) | iid<1L | iid>nrow(ped))){
    stop("the recorded pedigree does not cover these individuals or the fit, which happens when it was tracked late or reset after the fit",
         call.=FALSE)
  }
  have = logical(nrow(ped))
  have[solution@iid] = TRUE
  u = numeric(nrow(ped))
  u[solution@iid] = solution@ebv
  need = logical(nrow(ped))
  cur = unique(iid[is.na(take)])
  while(length(cur)>0L){
    cur = cur[cur>0L]
    cur = cur[!need[cur] & !have[cur]]
    need[cur] = TRUE
    cur = unique(c(ped[cur,1:2]))
  }
  for(i in which(need)){
    s = ped[i,1]
    d = ped[i,2]
    u[i] = 0.5*(if(s>0L) u[s] else 0) + 0.5*(if(d>0L) u[d] else 0)
  }
  return(u[iid])
}

# setEBV for a PedSol, with value "gv" or "bv" read as for an RRsol: "gv"
# adds the estimate for the first level of the fixed effects, and "bv" is
# centered on the training population, or on targetPop when it is given
.setEBVPed = function(pop, solution, value, targetPop, append, simParam,
                      onlyValues){
  if(!is(pop, "Pop")){
    stop("pedigree BLUP estimates can only be set for an object of class Pop",
         call.=FALSE)
  }
  value = tolower(value)
  u = .pedPredict(solution, pop@iid, simParam)
  if(value=="gv"){
    ebv = u + solution@fixEff[1]
    name = paste0("est_GV_", solution@trait)
  }else if(value=="bv"){
    if(is.null(targetPop)){
      center = solution@center
    }else{
      center = mean(.pedPredict(solution, targetPop@iid, simParam))
    }
    ebv = u - center
    name = paste0("est_BV_", solution@trait)
  }else if(value %in% c("female", "male")){
    stop("This pedigree model does not produce GCA estimates.")
  }else{
    stop(paste0("value=", value, " is not a valid option"))
  }
  ebv = matrix(unname(ebv), ncol=1L, dimnames=list(NULL, name))
  if(onlyValues){
    return(ebv)
  }
  if(append){
    pop@ebv = cbind(pop@ebv, ebv)
  }else{
    pop@ebv = ebv
  }
  return(pop)
}
