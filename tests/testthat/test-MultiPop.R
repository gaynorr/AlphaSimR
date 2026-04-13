
test_that("MultiPop assignment and replace methods", {

  founderPop = quickHaplo(nInd=100, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
  pop = newPop(founderPop, simParam=SP)

  # Append elements to a MultiPop using `[[<-` and defined names
  mpA = newMultiPop()
  mpA[["Group1"]] = pop[1:5]
  mpA[["Group2"]] = pop[6:10]
  expect_identical(names(mpA), c("Group1", "Group2"))
  expect_identical(mpA[["Group1"]], pop[1:5])

  # Append elements to a MultiPop using `[<-` and defined names
  mpB = newEmptyMultiPop()
  mpB[c('Group1','Group2')] = list(pop[1:5], pop[6:10])
  expect_identical(names(mpB), c("Group1", "Group2"))
  expect_identical(mpB["Group1"], newMultiPop(Group1 = pop[1:5]))

  # Append elements to a MultiPop using `$<-` and defined names
  mpC = newEmptyMultiPop()

  mpC$Group1 = pop[1:5]
  mpC$Group2 = pop[6:10]
  expect_identical(names(mpC), c("Group1", "Group2"))
  expect_identical(mpC$Group2, pop[6:10])

  # Test that the MultiPops created with different methods are identical
  expect_identical(mpA, mpB)
  expect_identical(mpA, mpC)
  expect_identical(mpB, mpC)

  # Replace elements in a MultiPop using `[[<-`, `[<-`, and `$<-`
  mpA[["Group1"]] = pop[11:15]
  mpA[[2]] = pop[16:20]
  mpB[c('Group1')] = list(pop[11:15])
  mpC$Group1 = pop[11:15]
  expect_identical(mpA[[1]], pop[11:15])
  expect_identical(mpA[["Group2"]], pop[16:20])
  expect_identical(mpB["Group1"], newMultiPop(Group1 = pop[11:15]))
  expect_identical(mpC$Group1, pop[11:15])
  mpC[2:3] = list(pop[16:20], pop[21:25])
  expect_identical(names(mpC), c("Group1", "Group2", ""))
  mpB[2:3] = newMultiPop(Group2 = pop[16:20], pop[21:25])
  expect_identical(mpB, mpC)

  # Delete elements from a MultiPop using `[[<-`, `[<-`, and `$<-`
  mpA[["Group2"]] = NULL
  mpB[3] = NULL
  mpB$Group2 = NULL
  mpC[c(2,3)] = NULL
  expect_identical(mpA, mpB)
  expect_identical(mpA, mpC)
  expect_identical(mpB, mpC)

  # Test that accessing a non-existent name returns NULL
  expect_null(mpA[["does_not_exist"]])
  expect_null(mpC$does_not_exist)

  # Assign, replace, and remove names of a MultiPop
  mpD = newMultiPop(pop[1:5], pop[6:10], pop[11:15], pop[16:20], pop[21:25])
  expect_null(names(mpD))
  names(mpD) = c("Group1", "Group2", "Group3", "Group4", "Group5")
  expect_identical(names(mpD), c("Group1", "Group2", "Group3", "Group4", "Group5"))

  names(mpD)[1:3] = c("GroupX", "GroupY", "GroupZ")
  expect_identical(names(mpD), c("GroupX", "GroupY", "GroupZ", "Group4", "Group5"))

  names(mpD)[1:3] = ""
  expect_identical(names(mpD), c("", "", "", "Group4", "Group5"))
  names(mpD) = NULL
  expect_null(names(mpD))

  # Multipop with 3 levels of nesting and direct name assignment
  mp2 = newMultiPop(pop1 = pop[1:20], pop2 = pop[21:40],
                    mpA = newMultiPop(pop3 = pop[41:60],
                                      mpB = newMultiPop(pop4 = pop[61:80], pop5 = pop[81:100])))
  mp2copy = mp2

  # `names` returns names of the top-level elements in a MultiPop, even if it has nested structure
  expect_identical(names(mp2), c("pop1", "pop2", "mpA"))

  # Access nested MultiPop with [[
  expect_identical(mp2[["mpA"]][['mpB']][['pop5']], pop[81:100])

  # Append elements to a MultiPop with nested structure
  mp2[["pop6"]] = pop[1:20]
  expect_identical(mp2[["pop1"]], mp2[["pop6"]])


  # Append and remove elements from a MultiPop with nested structure
  mp2[['mpA']][['mpB']]['pop6'] = mp2[["pop6"]]
  expect_identical(mp2[['mpA']][['mpB']]['pop6'], newMultiPop(pop6 = pop[1:20]))
  expect_identical(mp2[['mpA']][['mpB']][['pop6']], pop[1:20])
  mp2[['mpA']][['mpB']]['pop6'] = NULL
  expect_length(mp2[['mpA']][['mpB']], 2)

  mp2[['mpA']][['mpB']][[3]] = mp2["pop6"]
  expect_identical(mp2[['mpA']][['mpB']][[3]], newMultiPop(pop6 = pop[1:20]))
  mp2[['mpA']][['mpB']][[3]] = NULL
  expect_length(mp2[['mpA']][['mpB']], 2)

  # Delete additional element and check that the MultiPop is unchanged
  mp2[['pop6']] = NULL
  expect_identical(mp2, mp2copy)

  # Delete non-existent name should not change the MultiPop
  mp2[['popX']] = NULL
  mp2['pop6'] = NULL
  expect_identical(mp2, mp2copy)

  # A single-level MultiPop for testing error paths
  mp <- newEmptyMultiPop()
  mp[["A"]] <- pop[1:5]

  # $ getter
  dollar_get = methods::getMethod("$", "MultiPop")

  expect_error(
    dollar_get(mp, quote(c("A", "B"))),
    "$ requires a single name",
    fixed = TRUE
  )

  # [<- invalid value paths
  expect_error(
    {
      mp[1] <- 1L
    },
    "value must be a list or a Pop/MultiPop",
    fixed = TRUE
  )
  expect_error(
    {
      mp[1] <- list(1L)
    },
    "All elements of list must be Pop, MultiPop, or NULL",
    fixed = TRUE
  )

  # [[<- invalid index/value paths
  expect_error(
    {
      mp[[c(TRUE, TRUE)]] <- pop[1:5]
    },
    "logical index must select exactly one element",
    fixed = TRUE
  )
  expect_error(
    {
      mp[[c("A", "B")]] <- pop[1:5]
    },
    "only single character index allowed",
    fixed = TRUE
  )
  expect_error(
    {
      mp[[0]] <- pop[1:5]
    },
    "invalid numeric index",
    fixed = TRUE
  )
  expect_error(
    {
      mp[[c(1, 2)]] <- pop[1:5]
    },
    "invalid numeric index",
    fixed = TRUE
  )
  expect_error(
    {
      mp[[list(1)]] <- pop[1:5]
    },
    "index must be numeric, character, or logical",
    fixed = TRUE
  )
  expect_error(
    {
      mp[[1]] <- 1L
    },
    "value must be a Pop or MultiPop",
    fixed = TRUE
  )

  # $<- invalid name/value paths
  dollar_set = methods::getMethod("$<-", "MultiPop")
  expect_error(
    dollar_set(mp, quote(c("A", "B")), pop[6:10]),
    "$ requires a single name",
    fixed = TRUE
  )

  expect_error(
    dollar_set(mp, "A", 1L),
    "value must be a Pop or MultiPop",
    fixed = TRUE
  )
})

test_that("MultiPop show", {
  founderPop <- quickHaplo(nInd = 40, nChr = 1, segSites = 10)
  SP <- SimParam$new(founderPop)
  SP$addTraitA(10)
  pop <- newPop(founderPop, simParam = SP)

  # Empty MultiPop
  out0 <- capture.output(show(newEmptyMultiPop()))
  expect_true(any(grepl('An object of class "MultiPop" with 0 item\\(s\\)', out0)))

  # Nested + named MultiPop
  mp <- newMultiPop(
    top1 = pop[1:5],
    top2 = newMultiPop(
      mid1 = pop[6:10],
      mid2 = newMultiPop(low1 = pop[11:15],
                         low2 = pop[16:18])
    )
  )
  out1 <- capture.output(show(mp))
  expect_true(any(grepl("Level 1:", out1, fixed = TRUE)))
  expect_true(any(grepl('An object of class "Pop" with 5 individual\\(s\\)', out1)))
  expect_true(any(grepl('An object of class "MultiPop" with 2 item\\(s\\)', out1)))
  expect_true(any(grepl('An object of class "Pop" with 3 individual\\(s\\)', out1)))

  # Name suppression branch for "" and NA
  names(mp) <- c("", NA_character_)
  out2 <- capture.output(show(mp))
  expect_false(any(grepl('"NA" - ', out2, fixed = TRUE)))
  expect_false(any(grepl('"" - ', out2, fixed = TRUE)))
})

test_that("Multipop unnameMultiPop", {
  founderPop <- quickHaplo(nInd = 100, nChr = 1, segSites = 10)
  SP <- SimParam$new(founderPop)
  SP$addTraitA(10)
  pop <- newPop(founderPop, simParam = SP)

  mp2 = newMultiPop(pop1 = pop[1:20], pop2 = pop[21:40],
                    mpA = newMultiPop(pop3 = pop[41:60],
                                      mpB = newMultiPop(pop4 = pop[61:80], pop5 = pop[81:100])))

  # Test unnameMultiPop on MultiPop with nested structure
  expect_error(unnameMultiPop(mp2, level = 0), "level must be a positive integer or Inf", fixed = TRUE)

  u1mp2 = unnameMultiPop(mp2, level = 1)
  expect_null(names(u1mp2))
  expect_identical(names(u1mp2[[3]]), c("pop3", "mpB"))
  expect_identical(names(u1mp2[[3]][['mpB']]), c("pop4", "pop5"))

  u2mp2 = unnameMultiPop(mp2, level = 2)
  expect_identical(names(u2mp2), c("pop1", "pop2", "mpA"))
  expect_null(names(u2mp2[[3]]))
  expect_identical(names(u2mp2[['mpA']][[2]]), c("pop4", "pop5"))

  u3mp2 = unnameMultiPop(mp2, level = 3)
  expect_identical(names(u3mp2), c("pop1", "pop2", "mpA"))
  expect_identical(names(u3mp2[['mpA']]), c("pop3", "mpB"))
  expect_null(names(u3mp2[['mpA']][['mpB']]))

  u4mp2 = unnameMultiPop(mp2, level = 1:2)
  expect_null(names(u4mp2))
  expect_null(names(u4mp2[[3]]))
  expect_identical(names(u2mp2[[3]][[2]]), c("pop4", "pop5"))

  u5mp2 = unnameMultiPop(mp2, level = Inf)
  expect_null(names(u5mp2))
  expect_null(names(u5mp2[[3]]))
  expect_null(names(u5mp2[[3]][[2]]))

  # Test error paths for invalid level arguments
  mp <- newMultiPop(
    l1 = pop[1:5],
    l2 = newMultiPop(
      l2a = pop[6:10],
      l2b = newMultiPop(l3a = pop[11:15])
    )
  )
  # depth = 3

  # length(level) == 1 branch
  expect_error(unnameMultiPop(mp, level = 0), "level must be a positive integer or Inf", fixed = TRUE)
  expect_error(unnameMultiPop(mp, level = NA_real_), "level must be a positive integer or Inf", fixed = TRUE)
  expect_error(unnameMultiPop(mp, level = 1.5), "level must be a positive integer or Inf", fixed = TRUE)
  expect_error(unnameMultiPop(mp, level = 4), "requested level exceed max depth of x (3)", fixed = TRUE)

  # length(level) > 1 branch
  expect_error(
    unnameMultiPop(mp, level = c(1, NA_real_)),
    "levels must be a numeric vector of positive integers (no NA)",
    fixed = TRUE
  )
  expect_error(
    unnameMultiPop(mp, level = c(1, 0)),
    "levels must be a numeric vector of positive integers (no NA)",
    fixed = TRUE
  )
  expect_error(
    unnameMultiPop(mp, level = c(1, 2.5)),
    "levels must be a numeric vector of positive integers (no NA)",
    fixed = TRUE
  )
  expect_error(
    unnameMultiPop(mp, level = c(1, Inf)),
    "cannot mix Inf with integer levels",
    fixed = TRUE
  )
  expect_error(
    unnameMultiPop(mp, level = c(1, 4)),
    "requested level(s) exceed max depth of x (3)",
    fixed = TRUE
  )
})
