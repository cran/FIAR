lavaan3_merge <- function(input='') {

       # replace semicolons by newlines
    input <- gsub(";","\n", input)

    # break up in lines 
    model <- unlist( strsplit(input, "\n") )

    # remove comments starting with '#'
    model <- gsub("#.*","", model)

    # remove comments starting with '!'
    model <- gsub("!.*","", model)

    # strip all white space
    model <- gsub("[[:space:]]+", "", model)

    # remove empty lines
    idx <- which(nzchar(model))
    model <- model[idx]

    # check for multi-line formulas: they contain no "~" or "=" character
    # but before we do that, we remove everything within round brackets 
    # to avoid confusion with for example equal("f1=~x1") statements
    model.simple <- gsub("\\(.*\\)", "()", model)

    start.idx <- grep("[~=]", model.simple)
    end.idx <- c( start.idx[-1]-1, length(model) )
    model.orig    <- model
    model <- character( length(start.idx) )
    for(i in 1:length(start.idx)) {
        model[i] <- paste(model.orig[start.idx[i]:end.idx[i]], collapse="")
    }

    # ok, in all remaining lines, we should have a '~' operator
    idx.wrong <- which(!grepl("[~]", model))
    if(length(idx.wrong) > 0) {
        cat("lavaan: missing '~' operator in formula(s):\n")
        print(model[idx.wrong])
        stop("syntax error in lavaan model syntax")
    }

    # which ones are mm?
    model.simple <- gsub("\\(.*\\)", "()", model)
    mm.idx <- grep("=~", model.simple)

    # extract lv names
    mlist <- strsplit(model.simple, "=")
    list.names <- rep("", length(model))
    list.names[mm.idx] <- unlist( lapply(mlist[mm.idx], "[", 1) )

    # create formulaList
    formulaList <- lapply(as.list(model), as.formula)
    names(formulaList) <- list.names

    # decompose formulaList in 4 lists: 
    #  1. mm: definition of latent variables
    #  2. eqs:               regressions
    #  3. vc:                variance-covariance parameters
    #  4. int:               intercepts

    # 1. which lines are 'lv definitions'? one-sided formula's only!
    mm.idx <- which(lapply(formulaList, length) == 2)
    mm <- formulaList[mm.idx]
    if(length(mm) > 0 && is.null(names(mm))) {
        stop("unable to extract the latent variable `names` from the mm elements in the formulaList argument")
    }
    if(length(mm.idx) > 0) {
        formulaList <- formulaList[-mm.idx]
    }

    # 2. variance-covariances: double 'tildes'
    vc.idx <- which(lapply(formulaList, function(x) { length(x[[3]]) }) == 2)
    vc <- formulaList[vc.idx]
    if(length(vc.idx) > 0) {
        formulaList <- formulaList[-vc.idx]
    }

    # 3. all the others are single-tilde equations; but the intercept formulas
    # should only contain 1 single variable name
    int.idx <- which(lapply(formulaList, function(x) {length(all.vars(x))}) == 1)
    int <- formulaList[int.idx]
    if(length(int.idx) > 0) {
        formulaList <- formulaList[-int.idx]
    }

    # 4. everything else: regressions
    eqs <- formulaList

    # check syntax, and merge formulas with a common left-hand side variables
    if(length(mm) > 1) {
        mm <- check.mm(mm)
    }
    if(length(eqs) > 1) {
        eqs <- check.eqs(eqs)
    }
    if(length(vc) > 1) {
        vc <- check.vc(vc)
    }
    if(length(int) > 1) {
        int <- check.int(int)
    }

    out <- c(mm, eqs, vc, int)
    out
}

check.mm <- function(mm, warn=FALSE) {

    N <- length(mm)
    if(length( unique(names(mm)) ) == N) {
        return(mm)
    } else {
        # warn, and try to re-assemble
        if(warn) {
            cat("\nlavaan WARNING:\n")
            cat("Found several =~ formulas with the same latent variable name;\n")
            cat("These formulas have been merged\n")
            cat("To avoid this warning, make sure the latent variable names are unique\n")
            cat("\n")
        }
        lv.names <- unique(names(mm))
        lv.idx <- match(names(mm), lv.names)
        N2 <- length(lv.names)
        mm.new <- vector("list", length=N2)
        for(i in 1:N2) {
            mm.new[[i]] <- as.formula(paste("~",
                                      paste(lapply(mm[lv.idx == i], "[[", 2),
                                      collapse=" + ")))
        }
        names(mm.new) <- lv.names
    }

    mm.new
}


check.eqs <- function(eqs, warn=FALSE) {

    N <- length(eqs)
    eqs.y.names <- unlist(lapply(eqs, function(x) all.vars(x[[2]])))
    if(length( unique(eqs.y.names) ) == N) {
        return(eqs)
    } else {
        # warn, and try to re-assemble
        if(warn) {
            cat("\nlavaan WARNING:\n")
            cat("Found several ~ formulas with the same name for the dependent variable;\n")
            cat("These formulas have been merged\n")
            cat("To avoid this warning, make sure there is only one formula for each dependent variable\n")
            cat("\n")
        }
        y.names <- unique(eqs.y.names)
        y.idx <- match(eqs.y.names, y.names)
        N2 <- length(y.names)
        eqs.new <- vector("list", length=N2)
        for(i in 1:N2) {
            eqs.new[[i]] <- as.formula(paste(y.names[i], "~",
                                       paste(lapply(eqs[y.idx == i], "[[", 3),
                                       collapse=" + ")))
        }
    }

    eqs.new
}


check.vc <- function(vc, warn=FALSE) {

    N <- length(vc)
    vc.l.names <- unlist(lapply(vc, function(x) all.vars(x[[2]])))
    if(length( unique(vc.l.names) ) == N) {
        return(vc)
    } else {
        # warn, and try to re-assemble
        if(warn) {
            cat("\nlavaan WARNING:\n")
            cat("Found several ~~ formulas with the same name on the left side of the ~~ operator\n")
            cat("These formulas have been merged\n")
            cat("To avoid this warning, make sure there is only one formula for each variable\n")
            cat("\n")
        }
        vc.names <- unique(vc.l.names)
        vc.idx <- match(vc.l.names, vc.names)
        N2 <- length(vc.names)
        vc.new <- vector("list", length=N2)
        for(i in 1:N2) {
            vc.new[[i]] <- as.formula(paste(vc.names[i], "~~",
                                      paste(lapply(vc[vc.idx == i],
                                            function(x) x[[3]][[2]]),
                                      collapse=" + ")))
        }
    }

    # FIXME
    # check for symmetric entries (eg x1 ~~ x2 and x2 ~~ x1)

    vc.new
}


check.int <- function(int, warn=FALSE) {

    N <- length(int)
    int.y.names <- unlist(lapply(int, function(x) all.vars(x[[2]])))
    if(length( unique(int.y.names) ) == N) {
        return(int)
    } else {
        # warn, and try to re-assemble
        if(warn) {
            cat("\n"); cat("LAVAAN WARNING:"); cat("\n")
            cat(
"Found several intercept (~1) formulas with the same name for the dependent\n")
            cat(
"variable; these formulas have been removed.\n")
            cat(
"To avoid this warning, make sure there is only one intercept formula for\n")
            cat(
"each dependent variable.\n")
            cat("\n")
        }
        idx.duplicated <- which(duplicated(int.y.names))
        int.new <- int[-idx.duplicated]
    }

    int.new
}

