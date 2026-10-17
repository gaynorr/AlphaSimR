# selectOCS and its helpers, which selectCross(restrInbr=TRUE) also
# uses. They approximate optimal contribution selection with equal
# contributions by truncating on a merit that is penalized for sharing
# alleles with the selected parents. The derivation is given in the
# details section of selectOCS.

#' @title Penalty kernel for restricted inbreeding selection
#'
#' @description
#' Precomputes what the restricted inbreeding helpers need from the
#' genotypes. With \eqn{Z = Y/k}, the penalty on candidate \eqn{i} is
#' \eqn{z_i'u} with \eqn{u = Z'c}, which equals \eqn{(ZZ'c)_i}. When
#' there are no more individuals than loci, \eqn{YY'} is computed once,
#' so each iteration works with individuals rather than loci and is no
#' larger in memory than \eqn{Y}. Otherwise the penalty is computed
#' from \eqn{Y} directly.
#'
#' @param Y a genotype matrix coded as \eqn{2x - k}, where \eqn{x} is
#' the allele dosage and \eqn{k} the ploidy
#' @param ploidy the ploidy, \eqn{k}
#' @param useK should \eqn{YY'} be precomputed
#'
#' @return a list with the number of loci, the penalty for the whole
#' population, each individual's \eqn{z_i'z_i}, and functions for the
#' penalty and the expected fixation of a selected set
#'
#' @keywords internal
.ocsKernel = function(Y, ploidy, useK=nrow(Y)<=ncol(Y)){
  # Every entry of Y is a whole number, so every product and sum below
  # is a whole number held exactly in a double. The results are then
  # the same whatever BLAS performs them and however many threads it
  # uses, so they cannot depend on the machine.
  storage.mode(Y) = "double"
  nInd = nrow(Y)
  nLoci = ncol(Y)
  s = ploidy^2
  if(useK){
    K = tcrossprod(Y)
    rm(Y)
    self = diag(K)/s
    zu = function(take){
      rowSums(K[,take,drop=FALSE])/(s*length(take))
    }
    cross = function(a, b){
      sum(K[a,b,drop=FALSE])
    }
  }else{
    self = rowSums(Y^2)/s
    zu = function(take){
      drop(Y%*%colSums(Y[take,,drop=FALSE]))/(s*length(take))
    }
    cross = function(a, b){
      sum(colSums(Y[a,,drop=FALSE])*colSums(Y[b,,drop=FALSE]))
    }
  }
  # Expected fixation, mean(u^2), for equal contributions from take
  fix = function(take){
    cross(take,take)/(s*nLoci*length(take)^2)
  }
  # Expected fixation when each sex supplies half of the gametes, so
  # u is the average of the female and male allele frequencies
  fixPool = function(takeF, takeM){
    nF = length(takeF)
    nM = length(takeM)
    (cross(takeF,takeF)/nF^2 + 2*cross(takeF,takeM)/(nF*nM) +
       cross(takeM,takeM)/nM^2)/(4*s*nLoci)
  }
  return(list(nLoci=nLoci, zu0=zu(seq_len(nInd)), self=self, zu=zu,
              fix=fix, fixPool=fixPool))
}

#' @title Restricted inbreeding truncation without sexes
#'
#' @description
#' Selects \code{N} candidates for a fixed penalty \code{kappa} by
#' iterating truncation selection on
#' \code{m - kappa*(Z\%*\%u + rowSums(Z^2)/(2*N))} until the penalties,
#' which follow the allele frequencies of the selected set, converge.
#'
#' @param kappa the penalty on shared alleles
#' @param ker a kernel from \code{.ocsKernel}
#' @param m a merit vector, where larger values are better
#' @param N the number of individuals to select
#' @param maxit maximum number of iterations
#' @param tol convergence tolerance on the penalties per locus
#'
#' @return a list with the selected individuals and their expected
#' fixation
#'
#' @keywords internal
.ocsTrunc = function(kappa, ker, m, N, maxit=200L, tol=1e-5){
  # Zu holds Z%*%u, starting from the whole population's allele
  # frequencies
  Zu = ker$zu0
  # Replacing j with i in a set of N moves u by (z_i - z_j)/N, which
  # raises fixation by 2(z_i - z_j)'u/N + |z_i - z_j|^2/N^2. The first
  # part is Zu. The second contains z_i'z_i/N^2 for each candidate on
  # its own, mostly its homozygosity, and that part is included here.
  # The remaining part, from the pair, cannot be split between
  # candidates and so cannot enter a ranking.
  own = ker$self/(2*N)
  for(it in seq_len(maxit)){
    take = order(m - kappa*(Zu+own), decreasing=TRUE)[seq_len(N)]
    # A running average rather than a fixed damping factor, because a
    # fixed factor lets one or two borderline candidates flip in and
    # out of the set indefinitely. Convergence is judged on the
    # penalties rather than on the selected set for the same reason.
    d = (ker$zu(take)-Zu)/(it+1)
    Zu = Zu + d
    # Each z_i lies between -1 and 1, so a candidate's penalty moves by
    # at most nLoci times the largest change in u. This test is
    # therefore never stricter than requiring u to move by less than
    # tol at every locus, and it is what decides the ranking.
    if(max(abs(d))<tol*ker$nLoci) break
  }
  # Fixation is computed exactly for the set chosen, not for the
  # averaged u that guided the choice
  return(list(take=take, F=ker$fix(take)))
}

#' @title Restricted inbreeding truncation with sexes
#'
#' @description
#' Selects \code{nFemale} females and \code{nMale} males for a fixed
#' penalty \code{kappa}. Each sex supplies half of the gametes, so both
#' are penalized against the pooled allele frequencies, but they are
#' ranked separately.
#'
#' @param kappa the penalty on shared alleles
#' @param ker a kernel from \code{.ocsKernel}
#' @param m a merit vector, where larger values are better
#' @param female individuals that are female
#' @param male individuals that are male
#' @param nFemale the number of females to select
#' @param nMale the number of males to select
#' @param maxit maximum number of iterations
#' @param tol convergence tolerance on the penalties per locus
#'
#' @return a list with the selected individuals and their expected
#' fixation
#'
#' @keywords internal
.ocsTruncSex = function(kappa, ker, m, female, male, nFemale, nMale,
                        maxit=200L, tol=1e-5){
  Zu = ker$zu0
  # As in .ocsTrunc, but each selected female contributes 1/(2*nFemale)
  # and each male 1/(2*nMale), so replacing one moves u by half as much
  # per individual and the term for a candidate on its own is scaled
  # by its sex
  ownF = ker$self[female]/(4*nFemale)
  ownM = ker$self[male]/(4*nMale)
  for(it in seq_len(maxit)){
    takeF = female[order(m[female] - kappa*(Zu[female]+ownF),
                         decreasing=TRUE)[seq_len(nFemale)]]
    takeM = male[order(m[male] - kappa*(Zu[male]+ownM),
                       decreasing=TRUE)[seq_len(nMale)]]
    # Running average and convergence test, as in .ocsTrunc
    d = ((ker$zu(takeF)+ker$zu(takeM))/2-Zu)/(it+1)
    Zu = Zu + d
    if(max(abs(d))<tol*ker$nLoci) break
  }
  return(list(take=c(takeF, takeM), F=ker$fixPool(takeF, takeM)))
}

#' @title Find the penalty that meets an expected fixation target
#'
#' @description
#' Finds the smallest penalty for which \code{f} returns a selected set
#' meeting the target, by doubling an initial guess until the target
#' is met and then bisecting.
#'
#' @param f a function of the penalty returning a list with element F
#' @param kappa0 the initial guess for the upper bound of the penalty
#' @param target the maximum allowed expected fixation
#' @param nDouble maximum number of times the upper bound is doubled
#' @param relTol bisection stops once the interval holding the penalty
#' is narrower than this proportion of its upper end
#' @param maxBisect maximum number of bisection steps
#'
#' @return the list returned by \code{f}, with the penalty used and
#' whether the target was met
#'
#' @keywords internal
.ocsSolve = function(f, kappa0, target, nDouble=20L, relTol=1e-4,
                     maxBisect=60L){
  ans = f(0)
  if(ans$F<=target){
    return(c(ans, list(kappa=0, met=TRUE)))
  }
  lo = 0
  hi = kappa0
  ansHi = f(hi)
  k = 0L
  while(ansHi$F>target && k<nDouble){
    # The bound just tried missed the target, so it is a lower bound
    lo = hi
    hi = 2*hi
    ansHi = f(hi)
    k = k + 1L
  }
  if(ansHi$F>target){
    # The caller warns and uses this set, which is the most diverse
    # one found
    return(c(ansHi, list(kappa=hi, met=FALSE)))
  }
  # The set kept always meets the target, so stopping early only costs
  # a little merit, never feasibility
  i = 0L
  while((hi-lo)>relTol*hi && i<maxBisect){
    mid = (lo+hi)/2
    ans = f(mid)
    if(ans$F>target){
      lo = mid
    }else{
      hi = mid
      ansHi = ans
    }
    i = i + 1L
  }
  return(c(ansHi, list(kappa=hi, met=TRUE)))
}

#' @title Select individuals with restricted inbreeding
#'
#' @description
#' Selects individuals by an approximation to optimal contribution
#' selection, maximizing merit while restricting the expected increase
#' in fixation. No crosses are made, so the selected individuals can be
#' crossed with a custom crossing plan, for example with
#' \code{\link{makeCross}}. \code{\link{selectCross}} with
#' \code{restrInbr=TRUE} makes the same selection and crosses at random.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param nInd the number of individuals to select. These individuals
#' are selected without regard to sex and it supersedes values for
#' nFemale and nMale. Thus if the simulation uses sexes, it is likely
#' better to leave this value as NULL and use nFemale and nMale instead.
#' @param nFemale the number of females to select. This value is ignored
#' if nInd is set.
#' @param nMale the number of males to select. This value is ignored
#' if nInd is set.
#' @param trait the trait for selection. Either a number indicating
#' a single trait or a function returning a vector of length nInd.
#' @param use the selection criterion. Either a character
#' (genetic values "gv", estimated breeding values "ebv", breeding values
#' "bv", phenotypes "pheno", or randomly "rand") or a function returning
#' a vector of length nInd.
#' @param selectTop selects highest values if true.
#' Selects lowest values if false.
#' @param inbrTarget the target for expected fixation. Its meaning
#' depends on \code{inbrType}. The default of 0.01 with
#' \code{inbrType="relative"} allows a loss of 1\% of the current
#' heterozygosity per generation, which corresponds to an effective
#' population size of 50 for a diploid.
#' @param inbrType either "relative", where \code{inbrTarget} is the
#' allowed increase in expected fixation as a proportion of the
#' remaining heterozygosity in \code{pop}, or "absolute", where
#' \code{inbrTarget} is the maximum allowed expected fixation.
#' @param snpChip an integer indicating which SNP chip genotypes are
#' used to measure expected fixation
#' @param useQtl should QTL genotypes be used instead of a SNP chip
#' to measure expected fixation. If TRUE, snpChip specifies which
#' trait's QTL to use.
#' @param returnPop should results be returned as a
#' \code{\link{Pop-class}}. If FALSE, only the index of selected
#' individuals is returned.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' \code{trait} or \code{use}
#'
#' @details
#' Individuals are chosen by an approximation to optimal contribution
#' selection in which every selected individual contributes equally to
#' the next generation.
#'
#' Genotypes at the SNP chip, or at a trait's QTL if \code{useQtl=TRUE},
#' are coded as \eqn{z = 2x/k - 1}, where \eqn{x} is the allele dosage
#' and \eqn{k} is the ploidy. For a set of selected individuals,
#' \eqn{u_j = 2\bar{p}_j - 1}{u_j = 2 pbar_j - 1}, where
#' \eqn{\bar{p}_j}{pbar_j} is the allele frequency at locus \eqn{j}
#' among the gametes they supply. When \code{nInd} is used, this is the
#' frequency among the selected individuals. When \code{nFemale} and
#' \code{nMale} are used, each sex supplies half of the gametes, so
#' \eqn{u_j} is the average of the female and male values. The expected
#' fixation of the selected individuals is
#' \deqn{F = \frac{1}{m}\sum_j u_j^2 = 1 - H/H_{max},}{F = (1/m) sum_j u_j^2 = 1 - H/H_max,}
#' where \eqn{m} is the number of loci, \eqn{H} is the mean expected
#' heterozygosity of their gametes and \eqn{H_{max} = 0.5} is its value
#' at an allele frequency of 0.5. For diploids, \eqn{F} equals the
#' average coancestry \eqn{c'Gc/2} of the selected individuals, where
#' \eqn{c} holds their equal contributions and \eqn{G = ZZ'/(m/2)} is a
#' genomic relationship matrix scaled with allele frequencies of 0.5.
#'
#' The target \eqn{F^*} is set by \code{inbrType}. With "relative",
#' \eqn{F^* = F_t + \Delta F(1 - F_t)}, where \eqn{F_t} is the expected
#' fixation of all of \code{pop} and \eqn{\Delta F} is
#' \code{inbrTarget}. This allows a loss of a proportion
#' \eqn{\Delta F} of the current heterozygosity per generation. Because
#' each individual carries \eqn{k} copies of every locus,
#' \eqn{\Delta F = 1/(kN_e)}{Delta F = 1/(k Ne)} for an effective
#' population size \eqn{N_e}{Ne}, so the default of 0.01 corresponds to
#' an effective population size of 50 for a diploid and 25 for a
#' tetraploid. With "absolute", \eqn{F^* =} \code{inbrTarget},
#' independent of \code{pop}. This is useful when \code{pop} has
#' already been preselected, so its own expected fixation is not a
#' suitable reference.
#'
#' Selection on merit \eqn{m_i}, given by \code{trait} and \code{use},
#' subject to \eqn{F \le F^*} is relaxed with a Lagrange multiplier.
#' Each candidate \eqn{i} is then scored as
#' \deqn{s_i = m_i - \kappa\left(z_i'u + \frac{c}{2} z_i'z_i\right),}{s_i = m_i - kappa (z_i'u + (c/2) z_i'z_i),}
#' where \eqn{c} is the contribution of each selected individual:
#' \eqn{1/N} for \eqn{N} individuals, or \eqn{1/(2N_f)} for a female
#' and \eqn{1/(2N_m)} for a male when using sexes. The first term
#' penalizes candidates carrying alleles that are common among the
#' selected individuals and rewards those carrying alleles that are
#' rare among them, pulling allele frequencies toward 0.5. The second
#' is the candidate's own contribution to fixation, which grows with
#' its homozygosity and matters most when few individuals are
#' selected. For a
#' given \eqn{\kappa}, \eqn{u} starts from the allele frequencies of
#' \code{pop}, the top candidates on \eqn{s_i} are selected, and
#' \eqn{u} is updated toward their frequencies with a running
#' average. This repeats until \eqn{u} converges. With \code{nFemale}
#' and \code{nMale}, females and males are ranked separately against
#' the shared \eqn{u}. The expected fixation of the final set is
#' computed exactly, and bisection finds the smallest \eqn{\kappa} that
#' meets the target. If \eqn{\kappa = 0} meets it, the result is
#' ordinary truncation selection. If no \eqn{\kappa} tried meets it, a
#' warning is given and the individuals selected with the largest
#' \eqn{\kappa} tried are returned. A message reports \eqn{\kappa}, the
#' achieved \eqn{F}, the target and \eqn{F_t}.
#'
#' A random set of \eqn{N} individuals of ploidy \eqn{k}, contributing
#' equally, increases expected fixation by about \eqn{(1 - F_t)/(kN)}.
#' When using sexes, \eqn{N} is replaced by \eqn{4N_fN_m/(N_f + N_m)}
#' for \eqn{N_f} females and \eqn{N_m} males. A selected set can do
#' better than a random one by choosing individuals whose alleles
#' balance each other, so the target can be met with fewer than
#' \eqn{1/(k\Delta F)}{1/(k Delta F)} individuals, which is 50 for a
#' diploid with the default. Most of the restriction is then spent on
#' keeping allele frequencies balanced, though, so much less gain in
#' merit is made than with a larger number.
#'
#' The target holds only if the selected individuals are used as
#' assumed: each contributes equally to the next generation and, with
#' \code{nFemale} and \code{nMale}, each sex supplies half of the
#' gametes. A custom crossing plan should therefore give every selected
#' individual the same number of progeny, or the realized fixation will
#' differ from the target. The method controls the allele frequencies
#' of the gametes the selected individuals are expected to supply.
#' Sampling at meiosis, including double reduction in polyploids, adds
#' further drift in the progeny that the method does not account for.
#'
#' The method is a Lagrangian relaxation solved by repeated
#' linearization rather than an exact solution of optimal contribution
#' selection. Swapping one selected individual for another changes
#' fixation through a second order term with three parts: two for the
#' individuals on their own, which the score includes, and one for the
#' relationship between them, which cannot be split between candidates
#' and so is ignored.
#'
#' @return Returns an object of \code{\link{Pop-class}}, or the indices
#' of the selected individuals if \code{returnPop=FALSE}. With
#' \code{nFemale} and \code{nMale}, the selected females come before
#' the selected males.
#'
#' @family selection functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Select 4 individuals while restricting the increase in expected
#' #fixation measured at the SNP chip
#' parents = selectOCS(pop, nInd=4, inbrTarget=0.5, simParam=SP)
#'
#' #Cross them in a circle, so each parent is used equally
#' crossPlan = cbind(1:4, c(2:4, 1))
#' pop2 = makeCross(parents, crossPlan, nProgeny=2, simParam=SP)
#'
#' @export
selectOCS = function(pop, nInd=NULL, nFemale=NULL, nMale=NULL, trait=1,
                     use="pheno", selectTop=TRUE, inbrTarget=0.01,
                     inbrType="relative", snpChip=1, useQtl=FALSE,
                     returnPop=TRUE, simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  # A HybridPop has no genotypes to measure fixation with, and a
  # MultiPop has no single set of allele frequencies to restrict
  if(is(pop,"MultiPop")){
    stop("selectOCS does not support a MultiPop")
  }
  if(is(pop,"HybridPop")){
    stop("selectOCS does not support a HybridPop")
  }
  take = .ocsSelectParents(pop=pop, nInd=nInd, nFemale=nFemale,
                           nMale=nMale, trait=trait, use=use,
                           selectTop=selectTop, inbrTarget=inbrTarget,
                           inbrType=inbrType, snpChip=snpChip,
                           useQtl=useQtl, simParam=simParam,
                           nThreads=nThreads, ...)
  if(returnPop){
    return(pop[take])
  }else{
    return(take)
  }
}

#' @title Select parents with restricted inbreeding
#'
#' @description
#' The implementation behind \code{\link{selectOCS}}. See the details
#' section of \code{\link{selectOCS}} for the method.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param nInd the number of individuals to select, or NULL
#' @param nFemale the number of females to select
#' @param nMale the number of males to select
#' @param trait the trait for selection
#' @param use the selection criterion
#' @param selectTop selects highest values if true
#' @param inbrTarget the target for expected fixation
#' @param inbrType "relative" or "absolute"
#' @param snpChip which SNP chip, or which trait's QTL if \code{useQtl}
#' @param useQtl use QTL genotypes instead of a SNP chip
#' @param simParam an object of class \code{\link{SimParam}}
#' @param nThreads number of threads to use if OpenMP is available
#' @param ... additional arguments if using a function for trait
#'
#' @return indices of the selected individuals
#'
#' @keywords internal
.ocsSelectParents = function(pop, nInd, nFemale, nMale, trait, use,
                             selectTop, inbrTarget, inbrType, snpChip,
                             useQtl, simParam, nThreads, ...){
  inbrType = tolower(inbrType)
  if(length(inbrType)!=1 || !(inbrType%in%c("relative","absolute"))){
    stop("inbrType must be \"relative\" or \"absolute\"")
  }
  if(length(inbrTarget)!=1 || !is.numeric(inbrTarget) ||
     is.na(inbrTarget) || inbrTarget<0 || inbrTarget>1){
    stop("inbrTarget must be a single value between 0 and 1")
  }

  # Merit, oriented so that larger is always better
  m = getResponse(pop=pop, trait=trait, use=use, simParam=simParam,
                  nThreads=nThreads, ...)
  if(is.matrix(m)){
    if(ncol(m)!=1) stop("response must have a single column")
  }
  m = c(m)
  if(!selectTop){
    m = -m
  }

  # Genotypes coded as 2x-k, so Z = Y/k runs from -1 to 1 and a column
  # mean of Z is 2p-1 for allele frequency p at any ploidy. Y is kept
  # in whole numbers so that .ocsKernel can work exactly.
  if(useQtl){
    Y = pullQtlGeno(pop, trait=snpChip, simParam=simParam,
                    nThreads=nThreads)
  }else{
    Y = pullSnpGeno(pop, snpChip=snpChip, simParam=simParam,
                    nThreads=nThreads)
  }
  Y = 2*Y - pop@ploidy
  ker = .ocsKernel(Y, pop@ploidy)
  rm(Y)

  currentF = ker$fix(seq_len(pop@nInd))
  if(inbrType=="relative"){
    target = currentF + inbrTarget*(1-currentF)
  }else{
    target = inbrTarget
  }

  # The penalty is in units of merit per unit of z'u, so the ratio of
  # their spreads is a starting point of the right magnitude. The
  # doubling in .ocsSolve corrects it if it is too small.
  sdM = sd(m)
  sdZu = sd(ker$zu0)
  if(!is.finite(sdM) || sdM<=0) sdM = 1
  if(!is.finite(sdZu) || sdZu<=0) sdZu = 1
  kappa0 = sdM/sdZu

  if(!is.null(nInd)){
    if(nInd<1) stop("nInd must be >= 1")
    if(pop@nInd<nInd){
      nInd = pop@nInd
      warning("Suitable candidates smaller than nInd, returning ",nInd," individuals")
    }
    f = function(kappa){
      .ocsTrunc(kappa=kappa, ker=ker, m=m, N=nInd)
    }
  }else{
    if(simParam$sexes=="no")
      stop("You must specify nInd when simParam$sexes is `no`")
    if(is.null(nFemale))
      stop("You must specify nFemale if nInd is NULL")
    if(is.null(nMale))
      stop("You must specify nMale if nInd is NULL")
    if(nFemale<1) stop("nFemale must be >= 1")
    if(nMale<1) stop("nMale must be >= 1")
    female = checkSexes(pop=pop, sex="F", simParam=simParam)
    male = checkSexes(pop=pop, sex="M", simParam=simParam)
    if(length(female)<nFemale){
      nFemale = length(female)
      warning("Suitable candidates smaller than nFemale, returning ",nFemale," individuals")
    }
    if(length(male)<nMale){
      nMale = length(male)
      warning("Suitable candidates smaller than nMale, returning ",nMale," individuals")
    }
    f = function(kappa){
      .ocsTruncSex(kappa=kappa, ker=ker, m=m, female=female, male=male,
                   nFemale=nFemale, nMale=nMale)
    }
  }

  ans = .ocsSolve(f=f, kappa0=kappa0, target=target)
  if(!ans$met){
    warning("selectOCS target was not met; using the individuals ",
            "selected with the largest penalty tried")
  }
  message("selectOCS: kappa = ", signif(ans$kappa, 4),
          ", F = ", signif(ans$F, 4),
          " (target ", signif(target, 4),
          ", current ", signif(currentF, 4), ")")
  return(ans$take)
}
