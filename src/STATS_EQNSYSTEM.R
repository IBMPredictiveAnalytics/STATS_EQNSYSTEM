#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.1"

# History
# 12-Jun-2014 Original Version


helptext="Estimate a system of linear equations

STATS EQNSYSTEM
   eqnname: dep = indep indep ... eqnname: dep = indep indep ...
/INSTRUMENTS 
    VARIABLES=list of instrumental variables or AUTOMATIC=YES
/OPTIONS 
    METHOD = OLS or TWOSLS or SUR or THREESLS
    COVMETHOD = GEOMEAN or MAX or THEIL or NODFCOR
    MAXITER = number of iterations
    TOL = iteration convergence criterion
/SAVE 
    DATASET=dataset name ID=id variable
/PRINT 
    RESIDCOV=YES or NO RESIDCOR=YES or NO.

Example:  Three-stage least squares applied to a pair of equations
with consump and price endogenous.

STATS EQNSYSTEM 
demand: consump = price income 
supply: consump=price farmPrice trend
/INSTRUMENTS VARIABLES=income farmPrice trend
/OPTIONS METHOD=THREESLS.

The main subcommand specifies one or more equations to estimate.
Each equation is written in the form
equation name: dependent variable = list of regressors.
The equation name should follow the same rules as variable names,
and must not be duplicated.
Include 0 in the regressor list to suppress the constant term.

Instrumental variables must be provided for 2SLS and 3SLS.  The
instruments can be listed with the VARIABLES keyword.  
Alternatively, if appropriate, specify AUTOMATIC=YES 
to use all the right hand side variables that do not appear
on the left hand side of any equation as the instruments.
The same set of instruments is used for all equations.

METHOD specifies the estimation method.  There is no default.
OLS is ordinary least squares.
TWOSLS is two stage least squares.
SUR is seemingly unrelated regression.
THREESLS is three stage least squares.

COVMETHOD specifies how the residual covariance matrix
is calculated.  The default is GEOMEAN.  Details can be
found in
Judge, George G.; W. E. Griffiths; R. Carter Hill; Helmut Luetkepohl 
and Tsoung-Chao Lee (1985) The Theory and Practice of Econometrics, 
Second Edition, Wiley.

MAXITER specifies how many iterations may be carried out.
It does not apply to OLS or 2SLS.  The default is 1, i.e.,
no iteration.

TOL is the iteration convergence criterion.  The default is
1e-5.

The SAVE command is used to create a dataset of residuals.
There will be one column for each equation.  Specify a name
for the dataset.  It must not already be in use.
ID can specify an ID variable that will be the first
column of the dataset.

The PRINT subcommand specifies whether the residual
covariance and correlation matrices between equations
are displayed.  By default they are not displayed.

STATS EQNSYSTEM /HELP prints this help and does nothing else.
"

methods = list(sur="SUR", threesls="3SLS", threeslsgls="3SLS",
    threeslsiv="3SLS", threeslsgmm="3SLS", threeslsschmidt="3SLS",
    twosls="2SLS", ols="OLS")
submethods3sls=list(threesls="GLS", threeslsiv="IV",
    threeslsgmm="GMM", threeslsschmidt="Schmidt")

### MAIN ROUTINE ###
dosystem = function(eqns, method, insts=NULL, autoinst=FALSE, covmethod="geomean",
    maxiter=1, tol=1e-5, residcov=FALSE, residcor=FALSE,
    dataset=NULL, id=NULL, ignore=TRUE) {
    # Estimate system of equations
    
    setuplocalization("STATS_EQNSYSTEM")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Equation System")
    warningsprocname = gtxt("Equation System: Warnings")
    omsid="STATSEQNSYSTEM"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(systemfit), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "systemfit"),dostop=TRUE)
        }
    )

    if (!is.null(dataset)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
                dostop=TRUE)
        }
        if (dataset %in% alldatasets) {
            warns$warn(gtxt("The output dataset name must not already be in use"),
                dostop=TRUE)
        }
    }
#     item 1 is list of equation names (possibly corrected to be valid as variable names)
#     item 2 is list of left hand side variables
#     item 3 is list of right hand side
#     item 4 is list of formulas

    eqnspec = parseeqns(eqns, warns)

    allvars = union(eqnspec[[2]], union(eqnspec[[3]], insts))

    # automatic instruments specify all right hand variables not appearing on any left hand side
    if (autoinst) {
        if (!is.null(insts)) {
            warns$warn(gtxt("Cannot specify both automatic instruments and a list of instruments"),
                       dostop=TRUE)
        }
        insts = setdiff(eqnspec[[3]], eqnspec[[2]])
    }
    allargs = as.list(environment())
    # Equation variables are not guaranteed to exist, because they have not been
    # validated in the syntax declaration, so problems will be caught here
    dta = tryCatch(spssdata.GetDataFromSPSS(allvars, row.label=id, missingValueToNA=TRUE,
        factorMode="levels"),
        error=function(e) {warns$warn(e$message, dostop=TRUE)}
    )
    # submethods for 3SLS do not apply to other estimation methods
    if (covmethod == "theil") {
        covmethod = "Theil"
    }
    if (covmethod == "nodfcor") {
        covmethod = "noDfCor"
    }
    if (method %in% names(submethods3sls)) {
        ctl = systemfit.control(maxiter=maxiter, tol=tol,
            methodResidCov=covmethod, method3sls=submethods3sls[[method]],
            model=FALSE, x=FALSE, y=FALSE, z=FALSE)
    } else {
        ctl = systemfit.control(maxiter=maxiter, tol=tol,
                methodResidCov=covmethod,
                model=FALSE, x=FALSE, y=FALSE, z=FALSE)
    }

    if (!is.null(insts)) {
        instf = formula(paste("~", paste(insts, collapse="+"), collapse=""))
        res = tryCatch(systemfit(eqnspec[[4]], method=methods[[method]],
            inst=instf, control=ctl, data=dta),
            error = function (e) {warns$warn(e$message, dostop=TRUE)}
        )
    } else {
        res = tryCatch(systemfit(eqnspec[[4]], method=methods[[method]],
                control=ctl, data=dta),
                error = function (e) {warns$warn(e$message, dostop=TRUE)}
        )
    }
    ressum = summary(res)
    displayresults(allargs, res, ressum, warns)
    if (!is.null(dataset)) {
        savepred(allargs, res, ressum, dta, warns)
    }
}

parseeqns = function(eqns, warns) {
    # return list containing list of equation names, list of lhs variables,
    # list of rhs variables, followed by equation formulas
    
    # eqns has the form eqn name, ":", lhs var, "=", rhs vars
    # Equation names should be valid Statistics variable names and be
    # unique.  If they are not, an attempt is made to fix them
    
    lists = list()
    eqnnames = list()
    lhsvars = list()
    rhsvars = list()
    equations = list()
    numeqns = 0
    tokennum = 1
    numtokens = length(eqns)
    eqns[[numtokens+1]] = "<STOP>"
    eqns[[numtokens+2]] = ":"

    while (tokennum < numtokens) {
        tt = eqns[[tokennum]]
        if (tt == ":") {
            warns$warn(gtxt("Invalid syntax or duplicate equation name"), dostop=TRUE)
        }
        numeqns = numeqns + 1
        eqnnames[numeqns] = tt
        tokennum = tokennum + 1
        if (eqns[[tokennum]] != ":") {
            warns$warn(gtxt("Invalid syntax: Expected \":\""), dostop=TRUE)
        }
        tokennum = tokennum + 1
        # left hand side variable
        lhsvar = eqns[[tokennum]]
        lhsvars[[numeqns]] = lhsvar
        if (eqns[[tokennum+1]] != "=") {
            warns$warn(gtxt("Invalid syntax: Expected \"=\""), dostop=TRUE)
        }
        tokennum = tokennum + 2
        # right hand side variables
        rhsvar = list()
        numrhs = 0
        while (eqns[[tokennum + 1]] != ":") {
            numrhs = numrhs + 1
            rhsvar[[numrhs]] = eqns[[tokennum]]
            tokennum = tokennum + 1
        }
        rhsvars = union(rhsvars, rhsvar)
        # not a complete check but should catch typical errors
        if (any(c(":", "=", "+") %in% union(lhsvar, rhsvar))) {
            warns$warn(gtxt("One or more invalid variable names were found in the equations."), dostop=TRUE)
        }
        frml = paste(lhsvar, paste(rhsvar, collapse="+"), sep="~")
        equations[[numeqns]] = formula(frml)
    }

    if (length(eqnnames) != length(union(tolower(eqnnames), tolower(eqnnames)))) {
        warns$warn(gtxt("Duplicate equation names were found.  Names must be unique, ignoring case."), dostop=TRUE)
    }

    lists[[1]] = fixnames(eqnnames)
    lists[[2]] = unlist(lhsvars)
    # screen out 0 (noconstant) from variable list
    lists[[3]] = unlist(rhsvars[-match('0', rhsvars, nomatch=9999)])
    lists[[4]] = equations
    return(lists)
}
displayresults = function(allargs, res, ressum, warns) {
    # display results
    # allargs is the list of input parameters
    # res is the estimated model
    # ressum is summary(res)

    
    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    lbls = c(gtxt("Number of Equations"),
             gtxt("Estimation Method"),
             gtxt("3SLS Method"),
             gtxt("Residual Covariance Method"),
             gtxt("Instrumental Variables"),
             gtxt("Maximum Iterations"),
             gtxt("Tolerance"),
             gtxt("Actual Iterations"),
             gtxt("Convergence Achieved"),
             gtxt("System D. F."),
             gtxt("Output Dataset")
    )

    if (allargs[["maxiter"]] == 1 || (res$method %in% c("OLS", "2SLS"))) {
        convergence = gtxt("--NA--")
    } else if (res$iter < res$control$maxiter) {
        convergence = gtxt("Yes")
    } else {
        convergence = gtxt("NO")
    }
    values = c(
            length(res$eq),
            res$method,
            ifelse(res$method == "3SLS", res$control$method3sls, gtxt("--NA--")),
            res$control$methodResidCov,
            ifelse(!is.null(allargs$insts), paste(allargs$insts, collapse=" "), gtxt("--NA--")),
            res$control$maxiter,
            res$control$tol,
            ifelse(!is.null(res$iter), res$iter, gtxt("--NA--")),
            convergence,
            res$df.residual,
            ifelse(!is.null(allargs[["dataset"]]), allargs[["dataset"]], gtxt("--NA--"))
    )
    names(values) = lbls
    summarydf = data.frame(cbind(values))
    colnames(summarydf) = gtxt("Values")

    spsspivottable.Display(summarydf, title=gtxt("Equation System Summary"), 
                           templateName="EQNSYSTEMSUMMARY",
                           caption=gtxt("Results computed by R systemfit package"),
                           isSplit=FALSE
    )
    
    # iterate over equation results
    for (eqn in 1:length(res[['eq']])) {
        eq = ressum[['eq']][[eqn]]
        coefs = data.frame(coef(eq))
        # NaN values may occur in the table in some cases
        names(coefs) = c(gtxt("Coefficient"), gtxt("Std. Error"), "t", gtxt("Sig."))
        intercept = match("(Intercept)", row.names(coefs))
        if (!is.na(intercept)) {
            row.names(coefs)[[intercept]] = gtxt("(Constant)")
        }
        caption = c(
            gtxtf("R-squared: %.4f Adjusted R-squared: %.4f", eq$r.squared, eq$adj.r.squared),
            gtxtf("Residual std. error: %.4f, D.f.: %s", eq$sigma, eq$df[[2]]),
            gtxtf("SSR: %.4f, MSE: %.4f, Root MSE: %.4f, N: %s", 
                eq$ssr, eq$ssr/eq$df[[2]], sqrt(eq$ssr/eq$df[[2]]), eq$df[[1]] + eq$df[[2]])
        )
        if (NaN %in% coef(eq)) {
            caption[[length(caption)+1]] = 
                gtxt("NaN values indicate an estimation problem.  Possibly too few instruments were given")
        }
        caption = paste(caption, collapse="\n")
        spsspivottable.Display(coefs, 
            title=gtxtf("Equation: %s, Dependent Variable: %s", 
                eq$eqnLabel, allargs$eqnspec[[2]][[eqn]]),
            templateName="EQNSYSTEMSCOEFS",
            caption=caption,
            isSplit=FALSE
        ) 
    }
    if (allargs[["residcov"]]) {
        spsspivottable.Display(data.frame(res[["residCov"]]), title=gtxt("Residual Covariance Matrix"),
            templateName="EQNSYSTEMRESIDCOV", isSplit=FALSE)
    }
    if (allargs[["residcor"]]) {
        spsspivottable.Display(data.frame(ressum$residCor), title=gtxt("Residual Correlation Matrix"),
                               templateName="EQNSYSTEMRESIDCOR", isSplit=FALSE)
    }

    spsspkg.EndProcedure()
}

savepred = function(allargs, res, ressum, dta, warns) {
    # save residuals
    
    # dataset structure will have an id column and one
    # column for each equation
    
    dict = list()
    
    residdf = data.frame(ressum$residuals)

    # row names will always be nominal strings
    # guessing that we need length tripling for Unicode expansion
    rnlength = max(nchar(row.names(residdf)), na.rm=TRUE) * 3
    eqnnames = allargs$eqnspec[[1]]
    # make sure ID name does not conflict
    idname = "ID"
    eqnnamesuc = toupper(eqnnames)
    while (idname %in% eqnnamesuc) {
        idname = paste("ID", rpois(1, lambda=100), sep="")
    }
    
    dict[[1]] = c(idname, allargs[['id']], rnlength, paste("A", rnlength, sep=""), "nominal")
    for (name in 1:length(eqnnames)) {
        dict[[name+1]] = c(eqnnames[[name]], "", 0, "F8.2", "scale")
    }
    dict = spssdictionary.CreateSPSSDictionary(dict)
    spssdictionary.SetDictionaryToSPSS(allargs[["dataset"]], dict)
    tryCatch(spssdata.SetDataToSPSS(allargs[["dataset"]], data.frame(row.names(residdf), residdf)),
        error=function(e) {warns$warn(e$message, dostop=TRUE)}
    )
    spssdictionary.EndDataStep()
}

reserved = c('all','and', 'by','eq','ge','gt','le','lt','ne','not','or','to','with')
fixnames = function(thenames) {
    # fix equation names to be valid as SPSS names

    thenamescopy = thenames
    for (v in 1:length(thenames)) {
        # invalid characters (most of them)
        thenames[v] = gsub("[][+/*() \"\t'),=-]", "_", thenames[v])
        # reserved words
        if (!is.na(match(tolower(thenames[v]), reserved))) {
            thenames[v] = paste(substring(thenames[v], 1, 1), "_", substring(thenames[v], 2), sep="")
        }
        # regular first character
        if (substring(thenames[v],1,1) %in% list("$", "#", "@", "_")) {
            thenames[v] = paste('x', thenames[v], sep="")
        }
        # unique names ignoring case
        name = thenames[v]
        lowernames = lapply(thenames, tolower)
        while (v > 1 && tolower(name) %in% lowernames[1:(v-1)]) {
            name = paste(thenames[v], rpois(1,100), sep="")
        }
        thenames[v] = name
    }

    return(thenames)
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_EQNSYSTEM"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_EQNSYSTEM"))
}


Run = function(args) {
    #Execute the STATS SVM command

    cmdname = args[[1]]
    args = args[[2]]

    oobj = spsspkg.Syntax(list(
        spsspkg.Template("", subc="", ktype="literal", var="eqns", islist=TRUE),
        
        spsspkg.Template("VARIABLES", subc="INSTRUMENTS", ktype="existingvarlist", 
            var="insts", islist=TRUE),
        spsspkg.Template("AUTOMATIC", subc="INSTRUMENTS", ktype="bool", var="autoinst"),
        # all the three stage methods are actually the same currently since
        # the same instruments are used for all equations
        spsspkg.Template("METHOD", subc="OPTIONS", ktype="str", var="method",
            vallist=list("sur", "threesls", "threeslsgls", "threeslsiv", "threeslsgmm",
            "threeslsschmidt", "twosls", "ols")),
        spsspkg.Template("COVMETHOD", subc="OPTIONS", ktype="str", var="covmethod",
            vallist=list("geomean", "max", "theil", "nodfcor")),
        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", var="maxiter",
            vallist=list(1)),
        spsspkg.Template("TOL", subc="OPTIONS", ktype="float", var="tol",
            vallist=list(1e-10)),
        
        spsspkg.Template("RESIDCOV", subc="PRINT", ktype="bool", var="residcov"),
        spsspkg.Template("RESIDCOR", subc="PRINT", ktype="bool", var="residcor"),
        
        spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
        spsspkg.Template("ID", subc="SAVE", ktype="existingvarlist", var="id"),
        spsspkg.Template("RESIDUALS", subc="SAVE", ktype="bool", var="ignore"),
        spsspkg.Template("PRED", subc="SAVE", ktype="bool", var="pred")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dosystem")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
