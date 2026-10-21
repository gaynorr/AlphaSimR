context("SimParamDocs")

# Each user-facing method of SimParam has a help page of its own, written
# in R/SimParam-methods.R. Those pages cannot have a real \usage section,
# so R CMD check does not compare their arguments with the code, and a
# method that gained or renamed an argument, or a new method without a
# page, would otherwise go unnoticed. These tests take that check's place.

# Methods for AlphaSimR's own use, which are deliberately undocumented,
# and the constructor, which is documented on the SimParam page itself
undocumented = c("addToRec", "ibdHaplo", "updateLastId", "addToPed",
                 "pedInbreeding",
                 "clone", "initialize")

# Read the help pages from the source tree when the tests run from one,
# so that they check the pages just generated rather than an older
# installed copy, and from the installed package otherwise, as under
# R CMD check.
simParamRd = function(){
  root = test_path("..", "..")
  if(file.exists(file.path(root, "DESCRIPTION")) &&
     dir.exists(file.path(root, "man"))){
    files = list.files(file.path(root, "man"), pattern="^SimParam.*\\.Rd$",
                       full.names=TRUE)
    db = lapply(files, tools::parse_Rd)
    names(db) = basename(files)
  }else{
    db = tools::Rd_db("AlphaSimR")
    db = db[grepl("^SimParam.*\\.Rd$", names(db))]
  }
  db
}

rdArgNames = function(rd){
  tags = vapply(rd, attr, character(1), "Rd_tag")
  if(!any(tags=="\\arguments")) return(character())
  args = rd[[which(tags=="\\arguments")]]
  items = Filter(function(x) identical(attr(x, "Rd_tag"), "\\item"), args)
  nms = vapply(items, function(x) paste(unlist(x[[1]]), collapse=""),
               character(1))
  trimws(unlist(strsplit(nms, ",")))
}

# Topics linked from a parsed Rd page. A link written as
# \link[=topic]{text} keeps its target in the Rd_option attribute.
rdLinks = function(x){
  out = character()
  if(identical(attr(x, "Rd_tag"), "\\link")){
    opt = attr(x, "Rd_option")
    if(is.null(opt)){
      out = paste(unlist(x), collapse="")
    }else{
      out = sub("^=", "", paste(unlist(opt), collapse=""))
    }
  }
  if(is.list(x)){
    for(el in x) out = c(out, rdLinks(el))
  }
  out
}

test_that("PAGES every user-facing SimParam method has a help page", {
  db = simParamRd()
  methods = setdiff(names(SimParam$public_methods), undocumented)
  expect_true(length(methods) > 0L)
  missing = methods[!paste0("SimParam-", methods, ".Rd") %in% names(db)]
  expect_identical(missing, character())
})

test_that("ARGS each method page documents exactly the method's arguments", {
  db = simParamRd()
  methods = setdiff(names(SimParam$public_methods), undocumented)
  for(m in methods){
    rd = db[[paste0("SimParam-", m, ".Rd")]]
    if(is.null(rd)) next # reported by the test above
    expect_setequal(rdArgNames(rd), names(formals(SimParam$public_methods[[m]])))
  }
  expect_setequal(rdArgNames(db[["SimParam.Rd"]]),
                  names(formals(SimParam$public_methods$initialize)))
})

test_that("OVERVIEW the SimParam page links to every method page", {
  db = simParamRd()
  methods = setdiff(names(SimParam$public_methods), undocumented)
  links = rdLinks(db[["SimParam.Rd"]])
  expect_identical(setdiff(paste0("SimParam-", methods), links), character())
})
