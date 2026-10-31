context("seededResults")

# ---- 1. Stored results from seeded simulations ----

# A simulation run from a fixed seed must give the same results in every
# version of AlphaSimR unless a change is made on purpose. The scenarios in
# helper-seededResults.R are run here and compared with the values they
# gave when they were last updated. test-reproducibility.R compares runs
# within one build; these tests compare against earlier builds, so they
# fail when a change moves the random number stream, for example by
# drawing numbers in a different order or making individuals in batches of
# a different size.
#
# If a change does this on purpose, describe it under "Changes to results
# from a seed" in NEWS.md and regenerate the values below as described in
# helper-seededResults.R.
#
# Checksums, identifiers and sexes are exact on every platform. Effects,
# genetic values and phenotypes are real numbers whose last digits can
# depend on the platform's math library, so they are compared with a
# tolerance.

seededExpected = list(
  founders = list(
    haplo = 3780115, qtl = c(81448, 75280, 73381, 75099), snp = 75705,
    sex = "FFFFFMFFFMMMMFFFFMFF",
    addEff = c(2.4459315482661208, -0.1135551194038511, 0.17381439701319368,
               1.8069845421934694),
    gv = c(1.5986251909961029, 0.35996379414406582, -0.42466240095015068,
           -0.69661341200403781, -1.0117370355955082, 1.150701446913561,
           1.527486603803349, 0.82599378434384341, 0.69462439535625486,
           -0.96801163968615156, 0.29972968154089713, 2.3158980453344635),
    pheno = c(4.6668312245155601, -2.5443428971952553, 5.806950964623014,
              7.6663035344448636)),
  crossing = list(
    cross = 35070527,
    crossParents = paste0(
      "2 13,2 13,10 15,10 15,20 1,20 1,16 9,16 9,14 5,14 5,18 17,18 17,",
      "6 19,6 19,8 7,8 7,4 11,4 11,12 3,12 3,2 15,2 15,10 1,10 1,20 9,",
      "20 9,16 5,16 5,14 17,14 17,18 19,18 19,6 7,6 7,8 11,8 11,4 3,4 3,",
      "12 13,12 13,2 1,2 1,10 9,10 9,20 5,20 5,16 17,16 17,14 19,14 19,",
      "18 7,18 7,6 11,6 11,8 3,8 3,4 13,4 13,12 15,12 15"),
    selected = 3918834,
    selectedParents = paste0(
      "78 53,60 55,30 29,70 35,54 43,48 55,26 29,28 35,24 43,50 53,78 29,",
      "60 35,30 43,70 53,54 55,48 35,26 43,28 53,24 55,50 29"),
    dh = 1035330, selfed = 1038971, kept = "22,26",
    pheno = -1.0043428219809156),
  pedigree = list(mapDrop = 937686, popDrop = 945771),
  tracking = list(ibd = 12528626, mutated = 650284,
                  addEff = 1.767693000817629, domEff = 1.083219799295067),
  polyploid = list(cross = 9981662, reduced = 168515, doubled = 673660)
)

# Observations that are real numbers rather than exact
seededReal = c("addEff", "domEff", "gv", "pheno")

for (scenario in names(seededExpected)) {
  test_that(paste("SEEDED", scenario, "results match the stored values"), {
    observed = seededScenarios[[scenario]]()
    expected = seededExpected[[scenario]]
    expect_identical(sort(names(observed)), sort(names(expected)))
    for (item in names(expected)) {
      message = paste0("Seeded result '", scenario, "$", item, "' changed. ",
                       "If this is intended, list it under \"Changes to ",
                       "results from a seed\" in NEWS.md and update the ",
                       "stored values (see helper-seededResults.R).")
      if (item %in% seededReal) {
        expect_equal(observed[[item]], expected[[item]], tolerance = 1e-8,
                     info = message)
      } else {
        expect_identical(observed[[item]], expected[[item]], info = message)
      }
    }
  })
}
