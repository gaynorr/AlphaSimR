context("Segregating-site moments")

# Use the default solver tolerance; request tighter integration only where needed.
predictMoments = function(n,Ne=100,bp=1e8,mutRate=1e-8,
                          histNe=c(500,1500,6000,12000,100000),
                          histGen=c(100,1000,10000,100000,1000000),
                          tolerance=1e-5){
  predictSegSitesMoments(n,Ne,bp,mutRate,histNe,histGen,tolerance)
}

test_that("Constant-size moments match the exact coalescent formulas", {
  for(n in c(2,20,200)){
    ne = 123
    theta = 4.5
    h1 = sum(1/seq_len(n-1))
    h2 = sum(1/seq_len(n-1)^2)
    z = predictMoments(n,Ne=ne,bp=1,mutRate=theta/(4*ne),histNe=NULL,histGen=NULL)
    expect_equal(z$expectedTc,4*ne*h1,tolerance=1e-10)
    expect_equal(z$varianceTc,(4*ne)^2*h2,tolerance=1e-10)
    expect_equal(z$expectedS,theta*h1,tolerance=1e-10)
    expect_equal(z$varianceS,theta*h1+theta^2*h2,tolerance=1e-10)
  }
})

test_that("Two-lineage moments match analytic survival integrals across epochs", {
  for(sizes in list(c(100,1000,50),c(1000,50,100))){
    boundaries = c(30,700)
    starts = c(0,boundaries)
    durations = c(diff(starts),Inf)
    survival = 1
    meanTime = secondTime = 0
    for(i in seq_along(sizes)){
      rate = 1/(2*sizes[i])
      d = durations[i]
      area = if(is.infinite(d)) 1/rate else -expm1(-rate*d)/rate
      first = if(is.infinite(d)) 1/rate^2 else
        (1-exp(-rate*d)*(1+rate*d))/rate^2
      meanTime = meanTime+survival*area
      secondTime = secondTime+2*survival*(starts[i]*area+first)
      survival = survival*exp(-rate*d)
    }
    z = predictMoments(2,Ne=sizes[1],bp=1,mutRate=8/(4*sizes[1]),
                        histNe=sizes[-1],histGen=boundaries,tolerance=1e-8)
    expect_equal(z$expectedTc,2*meanTime,tolerance=2e-5)
    expect_equal(z$varianceTc,4*(secondTime-meanTime^2),tolerance=2e-5)
  }
})

test_that("Redundant boundaries and generation rescaling preserve moments", {
  a = predictMoments(40,Ne=100,bp=1,mutRate=6/400,histNe=500,histGen=70)
  b = predictMoments(40,Ne=100,bp=1,mutRate=6/400,
                     histNe=c(100,500,500),histGen=c(20,70,500))
  expect_equal(a,b,tolerance=2e-4)
  c = predictMoments(40,Ne=1000,bp=1,mutRate=6/4000,histNe=5000,histGen=700)
  expect_equal(c$expectedTc,10*a$expectedTc,tolerance=1e-9)
  expect_equal(c$varianceTc,100*a$varianceTc,tolerance=1e-9)
  expect_equal(c$expectedS,a$expectedS)
  expect_equal(c$varianceS,a$varianceS)
  zero = predictMoments(40,Ne=100,bp=1,mutRate=0,histNe=500,histGen=70)
  expect_equal(zero$expectedS,0)
  expect_equal(zero$varianceS,0)
})

test_that("Malformed population histories and mutation parameters are rejected", {
  for(n in list(1,2.5,NA,Inf,c(2,3),'2')){
    expect_error(predictMoments(n),"n must")
  }
  for(x in list(0,-1,NA,Inf,c(1,2),'1')){
    expect_error(predictMoments(2,Ne=x),"Ne must")
    expect_error(predictMoments(2,bp=x),"bp must")
  }
  for(x in list(-1,NA,Inf,c(1,2),'1')){
    expect_error(predictMoments(2,mutRate=x),"mutRate must")
  }
  for(sizes in list(0,-1,NA,Inf,'100')){
    expect_error(predictMoments(2,histNe=sizes),"histNe must")
  }
  for(times in list(numeric(0),c(0,2),c(2,1),c(1,1),c(1,NA),c(1,Inf))){
    expect_error(predictMoments(2,histNe=c(200,300),histGen=times),"histGen must")
  }
})

test_that("Relaxed integration tolerance retains accurate moments", {
  # A moderate sample exercises integration without a large-state benchmark.
  precise = predictMoments(200,tolerance=1e-7)
  fast = predictSegSitesMoments(200,Ne=100,bp=1e8,mutRate=1e-8,
    histNe=c(500,1500,6000,12000,100000),
    histGen=c(100,1000,10000,100000,1000000))
  expect_equal(fast,predictMoments(200,tolerance=1e-5))
  expect_equal(fast$expectedS,precise$expectedS,tolerance=1e-3)
  expect_equal(fast$varianceS,precise$varianceS,tolerance=1e-3)
  # Constant populations use exact terminal moments at either tolerance.
  expect_equal(predictMoments(200,histNe=NULL,histGen=NULL,tolerance=1e-5),
               predictMoments(200,histNe=NULL,histGen=NULL,tolerance=1e-7))
  for(x in list(0,-1,1,Inf,NA,c(1e-5,1e-6),"fast")){
    expect_error(predictMoments(20,tolerance=x),"tolerance must")
  }
})

test_that("Finite inputs that overflow the coalescent scale are rejected", {
  # Individually finite inputs can still overflow when multiplied or rescaled.
  expect_error(predictMoments(2,Ne=1e200,bp=1,mutRate=1e200,
                              histNe=NULL,histGen=NULL),
               "mutation parameter .* must be finite")
  expect_error(predictMoments(2,Ne=1e-200,bp=1,mutRate=0,
                              histNe=1e200,histGen=1),
               "cannot be represented on the reference scale")
  expect_error(predictMoments(2,Ne=1e-200,bp=1,mutRate=0,
                              histNe=1,histGen=1e200),
               "cannot be represented on the reference scale")
  # A positive historical size can also become zero after rescaling.
  expect_error(predictMoments(2,Ne=1e200,bp=1,mutRate=0,
                              histNe=1e-200,histGen=1),
               "cannot be represented on the reference scale")
})

test_that("Numerically unrepresentable integration fails explicitly", {
  # A nonzero epoch can be too short to divide into representable steps.
  tiny = .Machine$double.xmin * .Machine$double.eps
  expect_error(predictMoments(2,Ne=0.5,bp=1,mutRate=0,
                              histNe=1,histGen=tiny),
               "integration step underflow")
  # An unachievable tolerance must not silently return inaccurate moments.
  expect_error(predictMoments(2,Ne=100,bp=1,mutRate=0,
                              histNe=1000,histGen=100,tolerance=tiny),
               "Non-finite tree moments during integration")
  # Raw second moments overflow before the corresponding mean does.
  expect_error(predictMoments(2,Ne=1,bp=1,mutRate=0,
                              histNe=1e200,histGen=1),
               "Invalid calculated tree moments")
})

test_that("Mutation rate changes site moments without changing ancestry", {
  a = predictMoments(20,bp=100,mutRate=1e-3,histNe=500,histGen=70)
  b = predictMoments(20,bp=100,mutRate=2e-3,histNe=500,histGen=70)
  expect_equal(b$expectedTc,a$expectedTc)
  expect_equal(b$varianceTc,a$varianceTc)
  expect_equal(b$expectedS,2*a$expectedS)
  # Conditional mutation counts contribute a linear Poisson term;
  # variation in ancestry contributes the remaining quadratic term.
  expect_equal(b$varianceS-b$expectedS,4*(a$varianceS-a$expectedS))
  expect_equal(b,predictMoments(20,bp=200,mutRate=1e-3,
                                histNe=500,histGen=70))
})
