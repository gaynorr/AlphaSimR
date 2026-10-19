# Development-only, dependency-free review aid; never evaluates package code.
# From the package root: Rscript dev/audit_simParam.R [package-root] [--strict] [--operators]
# --strict exits with status 1 for omitted simParam arguments.
# An omission may be intentional; this checks calls, not whether context is needed.
# --operators includes uncertain [ and c calls; these are hidden by default.
# Or source this file and call audit_simParam(path) to get a data.frame.
#
# Records function signatures, then checks simulation-dependent direct calls and
# lapply/sapply/vapply/Map/mapply/do.call targets. The default report contains
# omissions only; unresolved targets and explicit replacement contexts are omitted. Operators [ and c are deliberately conservative:
# R's dynamic dispatch prevents this script from knowing the receiver's class.
# Aliases, external functions, computed targets and argument lists need manual
# review. A clean report is not proof of correct context propagation.

audit_simParam = function(path=".", includeOperators=FALSE){
  files = list.files(file.path(path, "R"), "\\.[Rr]$", full.names=TRUE)
  if(!length(files)) stop("No R source files found under ", path)
  nodes = list()
  definitions = list()
  nameOf = function(x){
    if(is.symbol(x) || (is.character(x) && length(x)==1L)) return(as.character(x))
    if(is.call(x) && as.character(x[[1]])[1] %in% c("::", ":::"))
      return(as.character(x[[3]]))
    ""
  }
  for(file in files){
    message("Auditing ", basename(file))
    src = readLines(file, warn=FALSE)
    pd = getParseData(parse(file, keep.source=TRUE))
    children = split(pd$id, pd$parent)
    rowIndex = setNames(seq_len(nrow(pd)), pd$id)
    parents = setNames(pd$parent, pd$id)
    extract = function(row){
      lines = src[seq.int(row$line1, row$line2)]
      lines[length(lines)] = substr(lines[length(lines)], 1, row$col2)
      lines[1] = substring(lines[1], row$col1)
      paste(lines, collapse="\n")
    }
    # Parse individual expressions to retain exact call locations. The parser
    # supplies the ancestry used to identify enclosing functions and quotations.
    exprs = list()
    for(id in pd$id[pd$token %in% c("expr", "expr_or_assign_or_help")]){
      value = tryCatch(parse(text=extract(pd[rowIndex[[as.character(id)]], ]))[[1]],
                       error=function(e) NULL)
      if(is.call(value)) exprs[[as.character(id)]] = value
    }
    functionNames = list()
    for(id in names(exprs)){
      value = exprs[[id]]
      if(!identical(value[[1]], as.name("function"))) next
      parent = parents[[id]]
      siblings = pd[match(children[[as.character(parent)]], pd$id), , drop=FALSE]
      binding = siblings$text[siblings$token=="SYMBOL_SUB"]
      assignment = exprs[[as.character(parent)]]
      if(!is.null(assignment) && nameOf(assignment[[1]]) %in% c("=", "<-"))
        binding = nameOf(assignment[[2]])
      label = if(length(binding)) tail(binding, 1) else paste0("<function@", pd[rowIndex[[id]],"line1"], ">")
      functionNames[[id]] = label
      # Also record functions without simParam so they are not mistaken for
      # unresolved callbacks. Keep ambiguous signatures for review.
      definitions[[label]] = c(definitions[[label]], list(value[[2]]))
    }
    for(id in names(exprs)){
      value = exprs[[id]]
      if(identical(value[[1]], as.name("function"))) next
      ancestors = character()
      parent = parents[[id]]
      while(parent!=0L){
        ancestors = c(ancestors, as.character(parent))
        parent = parents[[as.character(parent)]]
      }
      quoted = vapply(exprs[intersect(ancestors, names(exprs))], function(x)
        nameOf(x[[1]]) %in% c("quote", "substitute", "expression", "alist", "bquote"), logical(1))
      if(any(quoted)) next
      enclosing = intersect(ancestors, names(functionNames))
      # Only executable function bodies are relevant; top-level examples and
      # registration calls do not have a calling simulation to propagate.
      if(!length(enclosing)) next
      context = any(vapply(exprs[enclosing], function(x)
        "simParam" %in% names(x[[2]]), logical(1)))
      nodes[[length(nodes)+1L]] = list(file=basename(file), line=pd[rowIndex[[id]],"line1"],
        column=pd[rowIndex[[id]],"col1"], caller=functionNames[[enclosing[1]]],
        context=context, parameters=unique(unlist(lapply(exprs[enclosing],
          function(x) names(x[[2]])))), value=value)
    }
  }
  needsSimParam = function(name){
    any(vapply(definitions[[name]], function(x) "simParam" %in% names(x), logical(1)))
  }
  knownWithoutSimParam = function(name){
    if(name %in% names(definitions)) return(!needsSimParam(name))
    # Inspect only base bindings, never evaluate source or resolve user callbacks.
    value = get0(name, envir=baseenv(), inherits=FALSE)
    is.function(value) && !"simParam" %in% names(formals(value))
  }
  findings = list()
  add = function(node, callee, issue){
    # The default report answers one question: does a known callee accept
    # simParam while this call omits it? Unknown calls are not omissions.
    operatorReview = identical(issue,
      "review: operator may dispatch to a population method")
    if(!startsWith(issue, "missing:") && !(includeOperators && operatorReview)) return()

    findings[[length(findings)+1L]] <<- data.frame(file=paste0("R/",node$file),
      line=node$line, column=node$column, caller=node$caller, callee=callee,
      context=if(node$context) "simParam in scope" else "no simParam formal",
      issue=issue, call=paste(deparse(node$value, width.cutoff=120), collapse=" "),
      stringsAsFactors=FALSE)
  }
  for(node in nodes){
    call = node$value
    callee = nameOf(call[[1]])
    indirect = callee %in% c("lapply", "sapply", "vapply", "Map", "mapply", "do.call")
    if(indirect){
      wrapper = callee
      args = as.list(call)[-1]
      targetName = if(wrapper=="do.call") "what" else "FUN"
      targetIndex = match(targetName, names(args))
      if(is.na(targetIndex)) targetIndex = if(wrapper %in% c("Map", "mapply", "do.call")) 1L else 2L
      if(length(args)<targetIndex) next
      target = args[[targetIndex]]
      callee = nameOf(target)
      # A callback parameter can shadow a known function, including base names.
      if(is.symbol(target) && callee %in% node$parameters){
        if(node$context) add(node, callee, "review: indirect target not resolved")
        next
      }
      if(callee %in% c("[", "c")){
        if(includeOperators && node$context){
          supplied = "simParam" %in% names(args)
          if(wrapper=="do.call"){
            argIndex = match("args", names(args))
            if(is.na(argIndex)) argIndex = 2L
            bundle = if(length(args)>=argIndex) args[[argIndex]] else NULL
            supplied = is.call(bundle) && nameOf(bundle[[1]])=="list" &&
              "simParam" %in% names(as.list(bundle))
          }
          if(!supplied) add(node, callee,
            "review: operator may dispatch to a population method")
        }
        next
      }
      if(knownWithoutSimParam(callee)) next
      if(!callee %in% names(definitions)){
        inline = is.call(target) && identical(target[[1]], as.name("function"))
        if(node$context && !inline)
          add(node, callee, "review: indirect target not resolved")
        next
      }
      if(wrapper=="do.call"){
        argIndex = match("args", names(args))
        if(is.na(argIndex)) argIndex = 2L
        bundle = if(length(args)>=argIndex) args[[argIndex]] else NULL
        if(!is.call(bundle) || nameOf(bundle[[1]])!="list"){
          add(node, callee, "review: computed do.call arguments")
          next
        }
        call = as.call(c(list(as.name(callee)), as.list(bundle)[-1]))
      }else{
        # These wrappers supply one mapped argument, followed by their dots.
        if(wrapper %in% c("lapply", "sapply", "vapply")){
          matchedWrapper = tryCatch(match.call(get(wrapper, envir=baseenv()),
            node$value, expand.dots=FALSE), error=function(e) NULL)
          if(is.null(matchedWrapper)) next
          forwarded = as.list(matchedWrapper[["..."]])
          call = as.call(c(list(as.name(callee), as.name(".mappedElement")),
                           forwarded))
        }else{
          matchedWrapper = tryCatch(match.call(get(wrapper, envir=baseenv()),
            node$value, expand.dots=FALSE), error=function(e) NULL)
          if(is.null(matchedWrapper)) next
          mapped = as.list(matchedWrapper[["..."]])
          extra = matchedWrapper[["MoreArgs"]]
          if(!is.null(extra)){
            if(!is.call(extra) || nameOf(extra[[1]])!="list") next
            extra = as.list(extra)[-1]
          }
          call = as.call(c(list(as.name(callee)), mapped, extra))
        }
      }
    }
    if(callee %in% c("[", "c")){
      if(includeOperators && node$context && !"simParam" %in% names(as.list(call)))
        add(node, callee, "review: operator may dispatch to a population method")
      next
    }
    if(!callee %in% names(definitions) || !needsSimParam(callee)) next
    if("simParam" %in% names(call)){
      if(!identical(call[["simParam"]], as.name("simParam")))
        add(node, callee, "review: explicit context differs from simParam")
      next
    }
    signatures = unique(definitions[[callee]])
    if(length(signatures)!=1L){
      add(node, callee, "review: multiple function signatures")
      next
    }
    # Construct a function object solely for argument matching, never execution.
    fn = function() NULL
    formals(fn) = signatures[[1]]
    matched = tryCatch(match.call(fn, call, expand.dots=FALSE), error=function(e) NULL)
    if(is.null(matched)){
      add(node, callee, "review: arguments could not be matched")
    }else if(!"simParam" %in% names(matched)){
      hasDots = any(vapply(as.list(call)[-1], function(x)
        identical(x, as.name("...")), logical(1)))
      add(node, callee, if(hasDots) "review: simParam may be supplied through ..." else
        "missing: simParam argument")
    }else if(!identical(matched[["simParam"]], as.name("simParam"))){
      add(node, callee, "review: explicit context differs from simParam")
    }
  }
  if(!length(findings)) return(data.frame(file=character(), line=integer(),
    column=integer(), caller=character(), callee=character(), context=character(),
    issue=character(), call=character()))
  result = do.call(rbind, findings)
  result[order(!startsWith(result$issue, "missing:"),
               result$file, result$line, result$column), , drop=FALSE]
}

if(sys.nframe()==0L){
  args = commandArgs(trailingOnly=TRUE)
  if(any(args %in% c("-h", "--help"))){
    cat("Usage: Rscript dev/audit_simParam.R [package-root] [--strict] [--operators]\n",
        "Reports calls to known functions that accept simParam but omit it.\n",
        "--strict fails on omissions, including callers without a simParam formal.\n",
        "Omissions may be intentional. --operators adds uncertain [ and c calls.\n",
        "Source is parsed, never evaluated.\n")
  }else{
    strict = "--strict" %in% args
    includeOperators = "--operators" %in% args
    unknown = args[startsWith(args, "--") &
                     !args %in% c("--strict", "--operators")]
    if(length(unknown)) stop("Unknown option: ", paste(unknown, collapse=", "))
    paths = args[!args %in% c("--strict", "--operators")]
    if(length(paths)>1L) stop("Expected at most one package root")
    result = audit_simParam(if(length(paths)) paths else ".",
                            includeOperators=includeOperators)
    write.table(result, stdout(), sep="\t", row.names=FALSE, quote=TRUE)
    nMissing = sum(startsWith(result$issue, "missing:"))
    message(nMissing, " missing argument(s); ", nrow(result)-nMissing,
            " review item(s). Review items do not fail --strict.")
    if(!includeOperators) message("Uncertain operator calls omitted; use --operators to include them.")
    if(strict && nMissing>0L) quit(status=1L)
  }
}
