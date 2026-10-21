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
  if(N<=nFix){
    stop("there must be more records than levels of the fixed effects",
         call.=FALSE)
  }

  mme = calcPedMME(ped$mother, ped$father, ped$isDH, ped$row,
                   as.integer(fix), y, .pedF(ped, simParam), ped$ploidy)
  n = mme$nFix + mme$nBasis
  C0 = Matrix::sparseMatrix(i=mme$C$i, j=mme$C$j, x=mme$C$x,
                            dims=c(n, n), symmetric=TRUE)
  A0 = Matrix::sparseMatrix(i=mme$A$i, j=mme$A$j, x=mme$A$x,
                            dims=c(n, n), symmetric=TRUE)
  rhs = c(mme$rhs)

  # The pattern of C0+lambda*A0 is the same for every positive lambda
  L = Matrix::Cholesky(C0+A0, perm=TRUE, LDL=FALSE, super=FALSE)
  solveAt = function(lambda){
    L <<- Matrix::update(L, C0+lambda*A0)
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
    opt = optimize(m2ll, interval=range, tol=1e-6)
    if(min(abs(opt$minimum-range))<1e-3){
      boundary = TRUE
      warning("the REML estimate of the variance ratio is at the edge of the range searched, so the variance components may not be reliable",
              call.=FALSE)
    }
    lambda = exp(opt$minimum)
    sol = solveAt(lambda)
    Ve = (mme$yy-sum(sol*rhs))/(N-nFix)
    Vu = Ve/lambda
  }else{
    sol = solveAt(Ve/Vu)
  }

  W = Matrix::sparseMatrix(i=mme$W$i, j=mme$W$j, x=mme$W$x,
                           dims=c(mme$W$nrow, mme$W$ncol))
  ebv = as.vector(W %*% sol[-seq_len(nFix)])
  fixEffEst = sol[seq_len(nFix)]
  names(fixEffEst) = levels(fix)

  return(list(Vu=Vu, Ve=Ve, fixEff=fixEffEst, iid=ped$iid, ebv=ebv,
              boundary=boundary))
}
