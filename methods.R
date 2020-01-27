library(testit)
library(ChemmineR)
library(edgeR)
library(caret)
library(caretEnsemble)

library(annotate)
library(hgu133plus2.db)

PATH_TO_AOD8_SDF <- 'data/singledrug/AOD8_Structures.sdf'
REMOTE_URL_AOD8_SDF <- 'https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf' # 'https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf?version=1&modificationDate=1510673147932&api=v2'

argmin <- function(x){
  #' Returns the index of the minimum value in the vector x
  which(x == min(x))
}

argmax <- function(x){
  #' Returns the index of the maximum value in the vector x
  which(x == max(x))
}

retryOnError <- function(func, args, n.tries=Inf) {
  #' Repeatedly call a function until the call does not throw an error (use case is network requests that return 404 and throw an error).
  #' @return The result of \code{do.call(func, args)}
  #' @param func A character string representation of a function name to be called with \code{args}.
  #' @param {args} A list of arguments to be passed to func. Named lists are accepted wherein name specify parameters. 
  #' @param {n.tries} An integer specifying the number of times do.call(func, args) should be allowed to throw an error before giving up and returning FALSE.
  #' @examples
  #' retryOnError(getGEO, list('GSE32474'))
  #' \dontrun{
  #' # This will cause retryOnError to enter an infinite loop.
  #' retryOnError(do.call, list(c,1)) # The second argument to do.call must be a list. 
  #' }

  tempenv <- new.env()
  tempenv$triesLeft <- n.tries
  tempenv$success <- T
  errorCallBack <- function(e){
    print(e); warning('\nRetrying...') 
    tempenv$triesLeft <- tempenv$triesLeft-1
    tempenv$success <- F
    return(F)
  }
  repeat{
    res <- tryCatch(do.call(func, args), error=errorCallBack)
    if(tempenv$success | tempenv$triesLeft==0){
      if(tempenv$success) {
        message('Success.')
      }else {
        warning(sprintf('Maximum tries (%s) reached.', n.tries))
      }
      break
    }
  }
  rm(tempenv)
  return(res)
}


extractDrugName <- function(nscid, sdfSet=PATH_TO_AOD8_SDF){
  #' Get the corresponding drug name for a given NSC number
}


extractDrugNSC <- function(drugname, sdfSet=PATH_TO_AOD8_SDF){
  # TODO:: 
  #   - In download block, dont write to temp file. See comment.
  #   - Too many sdfSet variables in nested tryCatch? If it works it works.
  #   - Make search more robust
  #       - Its pretty greedy with the grepl statement, but works fine for now. 
  #       - Check Chemminr methods for this. 
  #' Match a drug name with its NSC identifier by searching for matches in each SDF's @datablock slot
  #' @param sdfSet An SDFSet object or the path to AOD8_Structures.sdf (available for download at https://wiki.nci.nih.gov/display/NCIDTPdata/).
  #' @param drugname The name of å drug for which an NSC identifier should be found.
  #' @examples
  #' extractDrugNSC('Gefitinib')
  #' extractDrugNSC('cisplatin')
  #' extractDrugNSC(715055)
  
  if(class(sdfSet) != 'SDFset'){
    # Find the AOD8_Structures.sdf file and read it
    sdfSet <- tryCatch(read.SDFset(sdfSet), error=function(e){
      warning(sprintf('%s\n sdfSet not found. Trying local AOD8_Structures.sdf...', e))
      sdfSet <- tryCatch(readSDFset(PATH_TO_AOD8_SDF), error=function(e) {
        warning(sprintf('%s\nLocal AOD8_Structures.sdf not found. Trying to download it...', e))
        x <- tryCatch(getURL(REMOTE_URL_AOD8_SDF), error=function(e) {
          stop(sprintf('%s\nCould not download AOD8_Structures.sdf. Please download it from https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf and pass its path as the argument to sdfSet in extractDrugNSC(drugname, sdfSet).',e))
        })
        tmpf <- tempfile()
        writeLines(x, tempf)
        sdfSet <- readSDFset(tmpf)
        file.remove(tmpf)
        return(sdfSet)
      })
      return(sdfSet)
    })
  } 
  
  # Search the sdfSet for any instances of drugname using grepl()
  out <- c()
  for(sdf in sdfSet@SDF){
    dblock <- sdf@datablock
    found <- sapply(dblock, function(x){
      grepl(tolower(drugname), tolower(x))
    })
    if(any(found)){
      out$nsc <- as.numeric(dblock[['NSC']])
      names(out) <- as.character(drugname)
      # names(out) <- dblock['dtpname:Brand Name']
    }
  }
  
  # Handle case of no NSC ID found
  if(length(out)==0){
    if(is.numeric(drugname)){
      # TODO::
      #   - Use extractDrugName() and add the returned name to names(out)
      out <- drugname # The user is not expected to pass a number, but ust return the number in-case it is the NSC ID itself.
    }else{
      stop(sprintf('An NSC identifier for %s couldn\'t be found. If you know this drug\'s NSC ID, please use that instead' , drugname))  
    }
  }
  
  return(out)
}


grepCellLine <- function(phenoTitle, reformat=T){
  # TODO:
  #   - Loop via sapply() is not necessary. Just do out <- gsub('','',out)
  #' [Not Exported] Lambda function applied to a list of cell names, such as that found in pData(eset)$title
  #' @return An (optionally) reformatted cell name extracted from pData()$title 
  #' @param phenoTitle A character string from which to extract a cell name
  #' @param remormat Should the cell name be reformatted (T)? Or should it simply be extracted and returned (F)? The default is highly recommended. There is almost no reason not to reformat.
  #' @example sapply(pData(getGEO('GSE32474')[[1]])$title, grepCellLine)
  out <- regmatches(phenoTitle, regexpr('(?<=:).*(?=\\s)', phenoTitle, perl=T))
  if(reformat) {
    out <- sapply(out, function(x) {
      x <- gsub('_', '-', x)       # Replace all occurences of "_" with "-"
      x <- gsub('/.*', '\\1', x)   # Select the first name of a pair separeted by / eg. "Name1/Name2" ==> "Name1"
      x <- gsub('\\(.*', '\\1', x) # Remove any characters in parentheses and the parentheses themselves eg. "HL-60(TB)" ==> "HL-60"
      x <- gsub('\\s.*', '\\1', x) # Remove characters following spaces and the space itself eg. "T47D NFkB15" ==> "T47D"
      return(trimws(x))
    })
  }
  return(out)
}


findCellAttr <- function(eset) {
  # Search through all the items in pData(eset) for a list of GI50 cell names
  # This is probably unecessary, because the $title attribute of GSE32474 will always exist in this package. However, other GSEs lack this attribute.
  # See extractCellsFromEset() for use case
  # Another use case exists in cellNames()
}

extractCellsFromEset <- function(eset) {
  #' Reformat the cell names from GSE32474 to match those in the GI50 dataset. The returned value is expected to be used to set pData(GSE32474)$cell <- returned_value for downstream funtionality. If that operation has already been done, the function just returns pData(GSE32474)$cell
  #' @return A list of the reformatted cell names
  #' @param eset An ExpressionSet object containing a \code{@phenoData} slot with a \code{$title} attribute holding a list of cell names
  #' @examples 
  #' gse32474 <- getGEO('GSE32474')
  #' pData(gse32474)$cell <- extractCellsFromEset(gse32474)

  if(is.list(eset)) eset <- eset[[1]]
  if(is.null(pData(eset)$cell) | length(pData(eset)$cell) == 0){
    if(is.null(pData(eset)$title)) {
      # TODO::
      #   - Use findCellAttr() to search the eset for cell names
      #       - Again, this is not necessary with the GSE32474 dataset
      stop('Unable to extract cell names because pData(eset)$title does not exist.')
    }
    return(sapply(pData(eset)$title, grepCellLine))
  }else{
    return(pData(eset)$cell)
  }
}


cellNames <- function(eset) {
  # TODO:: 
  #   - More descriptive warning messages
  #' A case-insensitive wrapper for the pData(eset)$cell attribute (useful if $Cell or $CELL exists but $cell doesn't).
  #' @examples 
  #' gse32474 <- getGEO('GSE32474')
  #' cellNames(gse32747)
  pD.names <- names(pData(eset))
  idx <- which(tolower(pD.names) == 'cell')
  if(length(idx) > 1){
    warning('More than one')
  }else if(length(idx) < 1){
    warning('Not found')
    found.which <- which(grepl('cell', tolower(pD.names)))
    if(length(found.which) == 1){
      idx <- found.which
    }else{
      # TODO::
      #   - Use findCellAttr() to search the eset for cell names and set idx <- findCellAttr()
      idx <- 'title'
    }
  }else{
    out <- pData(eset)[[idx]]  
  }
  return(out)
}


extractResponseByNSC <- function(nsc, gi50, col.out = c('NSC','NLOGGI50', 'CELL', 'LCONC', 'STDDEV')) {
  # TODO::
  #   - Make the cell name reformatting block a separate function
  # Extract samples from the GI50 dataset treated with a given drug. This function also reformats the items in the gi50$CELL column.
  #' @return A data.frame holding a subset of the GI50 dataset corresponding to all cells treated with the drug specified in \code{nsc}
  #' @param ncs Integer representing a single NCS identifier
  #' @param gi50 The GI50 dataset as a data.frame or matrix. NSC identifiers are expected to be in the gi50$NSC column.
  #' @param col.out A list of gi50 columns to return. Downstream methods require at least c('NSC', 'NLOGGI50').
  if(!is.numeric(nsc)) nsc <- extractDrugNSC(nsc)
  out <- gi50[which(gi50[,'NSC']==nsc), col.out]
  if(!is.null(out$CELL)){
    out$CELL <- sapply(out$CELL, function(x) {
      x <- sub(' ', '', x)                       # Remove the first occurence of a space
      # Single outlier names exist that would otherwise be thrown away without the conditional blocks below
      if(grepl('ADR-RES', x)) x <- 'NCI-ADR-RES' # An outlier "NCI/ADR-RES" ==> "NCI-ADR-RES"
      if(grepl('T-47D', x)) x <- 'T47D'          # Another outlier
      if(grepl('RXF393', x)) x <- 'RXF-393'      # Another outlier
      # General rules that can be used for reformatting across all names 
      x <- gsub('/.*', '\\1', x)                 # Select the first name of a pair separated by /
      x <- gsub('\\(.*', '\\1', x)               # Remove parenthesis eg."HL-60(TB)" ==> "HL-60"
      x <- gsub('\\s.*', '\\1', x)               # Remove names following a space eg. "T47D NFkB15" ==> "T47D"
      return(trimws(x))
    })  
  }else{
    stop('Can\'t truncate cell names because column gi50$CELL does not exist or was not included in col.out .')
  }
  return(out)
}


matchResponseToGSMByName <- function(responses, eset, LCONC=NULL) {
  # TODO:
  #   - LCONC = 'all' is not supported
  #       - A method for incorporating all the LCONC values in a regression model is needed.
  #       - Currently, passing LCONC='all' is treated as LCONC=NULL (the LCONC value that yields the greatest number of samples).
  #' @param responses A subset of the GI50 dataset corresponding to the drug of interest.
  #' @param eset An ExpressionSet object (i.e. GSE32474)
  #' @param LCONC Either a numeric value representing the desired log(concentration) of drug exposure that will be used to fit the model, or one of \code{'min'}, \code{'mid'}, or \code{'max'} used to select the minimun, median, or max of all LCONC values found in the GI50 dataset subsetted by drug. If a number is specified, and that number is not found, the function will try to find the nearest LCONC value among those existing in the subsetted dataset. If left NULL, the LCONC value yielding the greatest number of samples is used.
  #' @examples
  #' data(gi50)
  #' gse32747 <- getGEO('GSE32474')
  #' drug <- extractDrugNSC('Paclitaxel')
  #' responses <- extractResponseByNSC(drug, GI50)
  #' eset <- gse32474
  #' matchResponseToGSMByName(responses, gse32747)

  pData(eset)$cell <- extractCellsFromEset(eset)
  responses.cellNames <- responses$CELL
  eset.cellNames <- cellNames(eset)
  eset.cellNames.responseMatches <- match(eset.cellNames, responses.cellNames)
  
  all(eset.cellNames == responses[eset.cellNames.responseMatches, 'CELL']) 
  all(eset.cellNames == responses.cellNames[eset.cellNames.responseMatches])
  all(responses$CELL[eset.cellNames.responseMatches] == eset.cellNames)
  
  # Convenience variables to handle LCONC
  responses.lconc <- responses[eset.cellNames.responseMatches, 'LCONC']
  responses.lconc.levels <- as.numeric(levels(as.factor(responses.lconc)))
  responses.lconc.LCONC <- which(responses.lconc == LCONC) 
  
  # For the message() statement
  notFound <- list(b=FALSE, val=numeric())
  
  # Handles the user's LCONC argument by reassigning LCONC if a string option is passed or an error is thrown. Otherwise LCONC is left unchanged.
  if(!is.null(LCONC)){
    if(!suppressWarnings(is.na(as.numeric(LCONC)))) LCONC <- as.numeric(LCONC) # Make sure a valid LCONC number was not accidentally pass asa string. If so, convert it to numeric and carry on.
    if(is.character(LCONC)){
      LCONC <- tolower(LCONC)
      if(LCONC == 'all'){ 
        samples.counts <- sapply(responses.lconc.levels, function(x){
          # Count the frequency of each factor in the subsetted dataset
          sum(responses.lconc == x)
        })
        LCONC <- responses.lconc.levels[argmax(samples.counts)] # Select responses using the LCONC value with the highest frequency 
        warning(sprintf('The argument LCONC = "all" is not currently supported. The LCONC value that yields the maximum number of samples (LCONC=%s) will be used instead.', LCONC))
      }else if(LCONC == 'mid'){
        # Pick median LCON value from a sorted, unique list. The lowest value is selected to break a tie.
        idx <- floor((length(responses.lconc.levels)+1)/2) # The variable responses.lconc.levels is already sorted ascending, so just pick the middle number. Using median() will not work well for decimals, hence this solution.
        LCONC <- responses.lconc.levels[idx] 
      }else{
        # The user passed min or max, use the arg to call min() or max()
        LCONC <- do.call(LCONC, list(responses.lconc.levels))
      }
    }else if(is.numeric(LCONC) & length(responses.lconc.LCONC) == 0){ 
      LCONC.user <- LCONC
      if(-LCONC %in% responses.lconc){
        warning(sprintf('LCONC value %s not found for drug %s. Found -%s value. Using -%s instead.', LCONC.user, drug, LCONC, LCONC))
        LCONC <- -LCONC
      }else{ 
        # Use the nearest LCONC to the users argument that exists in the subsetted datset
        LCONC <- responses.lconc.levels[argmin(sqrt((LCONC - responses.lconc.levels)^2))]
        warning(sprintf('LCONC value %s not found for drug %s. The nearest LCONC value avilable in the dataset (%s) has been selected instead.', LCONC.previous, drug, LCONC))
      }
      notFound$b <- T; notFound$val <- LCONC.user # Only used for the message
    }else{
      # The user's LCONC argument is valid
      TRUE
    }
  }else{ 
    # Choose the LCONC value that will yield the most samples
    samples.counts <- sapply(responses.lconc.levels, function(x){ 
      # Count the frequency of each factor in the subsetted dataset
      sum(responses.lconc == x)
    })
    # Select the LCONC that yields the most samples
    LCONC <- responses.lconc.levels[argmax(samples.counts)] 
  }
  
  # These values will be returned
  responses.match <- responses[eset.cellNames.responseMatches, ]
  responses.match.LCONC <- responses.match[responses.match$LCONC == LCONC,]
  eset.cellNames.LCONC <- eset.cellNames[responses.match$LCONC == LCONC]
  x.col <- responses.match$LCONC == LCONC
  
  # Assertions
  all(responses.match$CELL == eset.cellNames)
  all(eset.cellNames[x.col] == responses.match.LCONC$CELL)
  
  # Count the samples eliminated after subsetting by LCONC?
  n.eliminated <- nrow(responses.match) - nrow(responses.match.LCONC)
  perc.eliminated <- round(n.eliminated / nrow(responses.match) * 100, digits=2)
  drug <- levels(as.factor(responses.match[,'NSC']))
  if(n.eliminated == 0){
    if(notFound$b == F){
      msg <- sprintf('Only one LCONC value (%s) for drug %s was found. This LCONC value yielded %s samples.', LCONC, drug, nrow(responses.match))
    }else{
      msg <- sprintf('Subsetting by LCONC value %s yielded the maximum number of samples available (%s).', LCONC, nrow(responses.match), drug)
    }
    message(msg)
  }else{
    message(sprintf('Subsetting by LCONC value %s eliminated %s of %s total samples (%s%%). If you\'d like to use the maximum number of samples possible, please set LCONC=NULL .', LCONC, n.eliminated, nrow(responses.match), perc.eliminated))  
  }
  
  # Return the dataset subsetted by LCONC along with the corresponding columns of eset
  list(responses=responses.match.LCONC, x.col=x.col, LCONC=LCONC)
}
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC=LCONC)
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='all')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='min')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='max')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='mid')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='1')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='-6')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC=1.0)


makeDataset <- function(responses, eset, LCONC=NULL, preProFunc=function(x){return(x)}) {
  # TODO:
  #   - Validate the dims of X and Y, or just rely on the downstream check within trainEnsemble()?
  #   - Expand preProFunc.extra by appending it to the returned list 
  #' @param preProFunc A lambda function of the form \code{function(x){return(x)}} or \code{function(x,y){return(list(x=x,y=y))}} used to preprocess the dataset. The \code{x} parameter is a count matrix and the \{y} parameter is a data.frame of the GI50 dataset. If using the y parameter, do not return the full GI50 dataframe, return a vector of response referenced as the $y attribute of your returned list.
  #' @examples 
  #' data(gi50)
  #' eset <- getGEO('GSE32474')
  #' drug <- extractDrugNSC('Doxorubicin')
  #' responses <- extractResponseByNSC(drug, gi50)
  #' preProFunc <- function(x){ list(x=x, extraStuff=list('hello')) }
  #' makeDataset(responses, eset, NULL, preProFunc)
  
  responses.matched <- matchResponseToGSMByName(responses, eset, LCONC=LCONC)
  x.col <- responses.matched$x.col # Use these columns of the count matrix to be sure they match the Y-values created by matchResponseToGSMByName()
  x <- t(Biobase::exprs(eset)[,x.col]) # The count matrix from the ExpressionSet object [features, samples] is the transpose of the canonical format expected by caret [samples, features]
  y <- responses.matched$responses # The Y-values created by matchResponseToGSMByName() are stored in $responses
  LCONC <- responses.matched$LCONC # The LCONC value selected by the user may not be the same as that returned by matchResponseToGSMByName() due to validation done by matchResponseToGSMByName()
  preProFunc.extra <- NULL # Used to hold any additional attributes returned by preProFunc
  
  # Assertions
  all(eset$cell[x.col] == y$CELL) # Make sure samples in X and Y match i.e. they both refer to the same cells
  
  # Handle preProFunc
  if(class(preProFunc) == 'function'){
    preProFunc.nargs <- length(as.list(args(preProFunc)))
    if(preProFunc.nargs == 2){
      # The user expects only x to be passed to preProFunc (the default)
      out.ppf <- preProFunc(x)
      if(is.list(out.ppf)){
        # The user returned a list from preProFunc as expected
        names(out.ppf) <- tolower(names(out.ppf))
        x <- out.ppf$x
        if(length(names(out.ppf)) > 1){
          # Assumes the returned list contains more than one name, and one of those names is "x"
          preProFunc.extra <- out.ppf[names(out.ppf) != 'x']
        }
      }else{
        # Assumes the user only returned x from preProFunc (not a list containing x)
        x <- out.ppf
      }
      y <- y$NLOGGI50
    }else if(preProFunc.nargs == 3) {
      # The user expects both x and y to be passed to preProFunc 
      out.ppf <- preProFunc(x,y)
      names(out.ppf) <- tolower(names(out.ppf))
      x <- out.ppf$x
      y <- out.ppf$y
      preProFunc.extra <- out.ppf[names(out.ppf) != c('x', 'y')]
    }else{
      # The user passed an object that was not a function, or their function did not have either 1 or 2 parameters
      warning('Your preProFunc was not applied to the dataset, please check its parameters. The argument to preProFunc should be a function of the form function(x){return(x)} or function(x,y){return(list(x=x,y=y))}).')
      warning('Skipping preprocessing...')
      x <- x
      y <- y$NLOGGI50
    }
  }else{
    warning('Your preProFunc was not applied to the dataset, because it was not recognized as a function. The argument to preProFunc should be a function of the form function(x){return(x)} or function(x,y){return(list(x=x,y=y))}).')
    warning('Skipping preprocessing...')
    x <- x
    y <- y$NLOGGI50
  }
  
  # Return the dataset and any extra data generated by the user
  list(x=x, y=y, LCONC=LCONC, preProFunc.extra=preProFunc.extra)
}
# makeDataset(extractResponseByNSC(extractDrugNSC('Paclitaxel'), gi50), gse32474, LCONC=NULL)
# makeDataset(extractResponseByNSC(extractDrugNSC('DOESNT_EXIST'), gi50), gse32474, LCONC=NULL)
# makeDataset(extractResponseByNSC(extractDrugNSC('Gefitinib'), gi50), gse32474, LCONC='100')


trainRFE <- function(drug, eset, gi50, LCONC=NULL, preProFunc=function(x){return(x)},
                     rfeControl=rfeControl(), trainControl=trainControl(),
                     method='svmLinear', multiFolds=2, trainPercent=.85, 
                     sizes=seq(nrow(eset)/3, nrow(eset)-1, length.out=10), 
                     ...) 
{
  # TODO::
  #   - Move some of these params to trainControl or rfeControl?
  #   - Default sizes should be coerced to integer
  
  #' @examples
  #' caretFuncs <- list(summary=function(...){c(defaultSummary(...))})
  #' svm.RFE <- trainRFE(613327, gse32474, data.gi50,
  #'                     function(x) {x},
  #'                     rfeControl = rfeControl(
  #'                       method = "repeatedcv",
  #'                       repeats = 2,
  #'                       number = 5,
  #'                       verbose = TRUE,
  #'                       functions = caretFuncs,
  #'                       seeds = set.seed(123),
  #'                       saveDetails = TRUE
  #'                     ),
  #'                     trainControl = trainControl(
  #'                       method = "LGOCV",
  #'                       number = 1,
  #'                       p = 0.8,
  #'                       seeds = set.seed(321),
  #'                       savePredictions = TRUE,
  #'                       classProbs = FALSE,
  #'                       verboseIter = TRUE
  #'                     ),
  #'                     multiFolds = 2, 
  #'                     sizes = seq(nrow(gse32474)/2, nrow(gse32474)-1, 10000), trainPercent = 0.85
  #'                   )
  #' pred <- svm.RFE$pred
  #' result <- svmRFE$results
  #' fit <- svmRFE$fit

  if(is.list(eset)) eset <- eset[[1]]
  if(!is.numeric(drug)) drug <- extractDrugNSC(drug) 
  if(is.null(pData(eset)$cell)) pData(eset)$cell <- extractCellsFromEset(eset) 
  responses <- extractResponseByNSC(drug, gi50)
  data <- makeDataset(responses, eset, LCONC=LCONC, preProFunc=preProFunc)
  LCONC <- data$LCONC
  
  if(length(data$y) != nrow(data$x)) stop('X and Y lengths of dataset differ.')
  
  partition.idx <- createDataPartition(y=data$y, times=1, p=trainPercent)
  train.rows <- partition.idx$Resample1
  test.rows <- !(1:length(data$y) %in% partition.idx$Resample1)
  
  data.train <- data.frame(y=data$y[train.rows], data$x[train.rows,])
  data.test <- data.frame(y=data$y[test.rows], data$x[test.rows,])
  
  rfeControl$index <- createMultiFolds(data.train[, 1], times=multiFolds)
  rfe.out <- rfe(x = data.train[,-(1)],
                 y = data.train$y,
                 method = method,
                 sizes = sizes,
                 rfeControl = rfeControl,
                 trainControl = trainControl, 
                 ...)
  
  return(rfe.out)
}


trainEnsemble <- function(drug=NULL, eset=NULL, gi50=NULL, dataset = NULL, LCONC=NULL, preProFunc=function(x) {return(x)},
                          esbl.trainControl=trainControl(), esbl.tuneList=list(svm=caretModelSpec(method='svm', tuneLength=3)), esbl.methodList=c(), 
                          stack.trainControl=trainControl(), stack.method='glm', stack.metric='RMSE', train.percent=0.85,
                          do.plot = T, ...) 
{
  # TODO:
  #   - Why partition the data before passing to caret? Can't you just get caret to split the dataset? That way the user's resampling indices will work. I think samples are currently being witheld from training that don't need to be held out.
  #   - Don't put the Rsq and RMSE in the plot title. Use text() or something (except text() throws an error in Rmarkdon (Error:plot.new() has not been called)).
  #' @param drug An integer specifying the NSC ID of the drug for which an ensemble will be trained. The ID can be found for a given drug name using extractDrugNSC().
  #' @param eset An ExpressionSet containing a count matrix and a \code{@phenoData} slot with a \code{$title} attribute holding a list of cell names e.g. \code{getGEO(GSE32474)}. The count matrix will be passed as the argument to parameter \code{x} in \code{preProFunc}.
  #' @param gi50 A matrix or dataframe of the GI50 dataset with column names of at least c('NLOGGI50', 'CELL').
  #' @param dataset A user defined dataset. Should be left NULL if drug, eset, and gi50 arguments are passed. \code{preProFunc} is not applied to this dataset.
  #' @param preProFunc A lambda function for preprocessing. It should be of the form \code{function(x){return(x)}} or \code{function(x,y){return(list(x=x,y=y))}} where parameter \code{x} is a count matrix and parameter \code{y} is the GI50 dataset (a dataframe). If using y, return the full dataframe.
  #' @param [esbl.trainControl] <list> The list created by a call to trainControl() which is passed as the trControl argument to caretEnsemble::caretList
  #' @param [esbl.methodList] <list> Alternative to esbl.tuneList passed as the methodList argument to caretEnsemble::caretList(). Can be left NULL if esbl.trainControl argument is passed
  #' @examlpes
  #' data(gi50)
  #' drug <- 613327
  #' eset <- getGEO('GSE32474')
  #' preProFunc <- function(x,y) {
  #'   list(x=edgeR::cpm(x), y=y[,'NLOGGI50'])
  #' }
  #' esbl.trainControl = trainControl(
  #'   method = "cv",
  #'   number=5,
  #'   p=0.9,
  #'   classProbs=FALSE,
  #'   returnData=FALSE,
  #'   verboseIter=TRUE
  #' )
  #' esbl.tuneList=list(
  #'   knn=caretModelSpec('knn', tuneGrid=data.frame(k=c(4,3,2))),
  #'   svmRadial=caretModelSpec('svmRadial', tuneGrid=data.frame(sigma=seq(1e-03, 1e-04, length.out=8), C=seq(3,1, length.out=8))),
  #'   glmnet=caretModelSpec('glmnet', tuneGrid=data.frame(alpha=seq(0.05, 0.001, length.out=5), lambda=seq(0.01, 0.001, length.out=5)))
  #' )
  #' stack.trainControl = trainControl(
  #'   method="cv",
  #'   number=5,
  #'   savePredictions="final",
  #'   classProbs=FALSE,
  #'   summaryFunction=defaultSummary,
  #'   returnData=FALSE
  #' )
  #' stack.method = 'glm'
  #' stack.metric = 'RSME'
  #' ensemble <- trainEnsemble(drug, eset, gi50)
  #' modelCor(resamples(ensemble$model.list))
  #' xyplot(resamples(ensemble$model.list))
  
  if(is.null(dataset)){
    if(is.list(eset)) eset <- eset[[1]]
    if(!is.numeric(drug)) drug <- extractDrugNSC(drug) 
    if(is.null(pData(eset)$cell)) pData(eset)$cell <- extractCellsFromEset(eset) 
    responses <- extractResponseByNSC(drug, gi50)
    data <- makeDataset(responses, eset, LCONC=LCONC, preProFunc=preProFunc)
    LCONC <- data$LCONC
    preProFunc.extra <- data$preProFunc.extra
  }else{
    data <- dataset
  }
  
  # Assertions
  if(length(data$y) != nrow(data$x)) stop(sprintf('For drug:%s x and y lengths differ.', drug))
  
  # Get the train and test rows
  train.rows <- sample(1:length(data$y), size=length(data$y)*train.percent, replace=FALSE)
  test.rows <- (1:length(data$y))[!(1:length(data$y) %in% train.rows)]
  
  assert('The test set does not contain any samples found in the train set.', !any(test.rows %in% train.rows))
  
  # Partition the dataset
  data.train <- list(y=data$y[train.rows], x=data$x[train.rows,])
  data.test <- list(y=data$y[test.rows], x=data$x[test.rows,])
  
  remove(data) # Save some space? Please stop crashing 5 drugs deep...
  
  # Ensemble parameters
  esbl.list.args <- list(
    x = data.train$x,
    y = data.train$y,
    trControl=esbl.trainControl,
    continue_on_fail = T
  )
  
  # Does the user want to use methodList? Default is tuneList.
  if(length(esbl.methodList)>0){ 
    esbl.list.args$methodList <- esbl.methodList
  }else{
    esbl.list.args$tuneList <- esbl.tuneList
  }
  
  # Train the ensemble
  message(sprintf('Training ensemble with %s samples and %s features.', nrow(esbl.list.args$x), ncol(esbl.list.args$x)))
  model.list <- do.call(caretList, esbl.list.args) # Train the model
  
  # Train the stack
  model.stack <- caretStack(
    model.list,
    method = 'glm', # TODO:: Hardcoded - add a new param
    metric = 'RMSE', # TODO:: Hardcoded - add a new param or delete. I think caretEnsemble selects ROC or RMSE depending on classification or regression.
    trControl=stack.trainControl
  )
  
  # Post processing
  predictions <- predict(model.stack, newdata=data.test$x)
  rmse <- RMSE(predictions, data.test$y)
  rsq <- R2(predictions, data.test$y)
  
  # Plotting
  if(do.plot){
    plot(data.test$y, predictions, main=sprintf('NSC%s LCONC%s Rsq=%s RMSE=%s', drug, LCONC, round(rsq, digits=3), round(rmse, digits=3)))
    abline(0,1)
  }
  
  list(
    LCONC = LCONC,
    model.list = model.list,
    model.stack = model.stack, 
    model.stack.testRMSE = rmse,
    model.stack.rsq = rsq,
    model.stack.dataTest.y = data.test$y,
    model.stack.predictions = predictions,
    preProFunc.extra = preProFunc.extra
  )
}


varImpByModel <- function(model.list) 
{
  # !!! INCOMPLETE
  # TODO:: 
  #   - The importance is not averaged like you thought. Need to implement by hand
  #' Calculate variable importance for each model in the ensemble
  #' @example getVarImp(model.ensemble$model.list)
  sapply(model.list, function(mod){
    if(length(mod)>0){
      tryCatch(varImp(mod$bestTune), error=warning)
      tryCatch(varImp(mod$finalModel), error=warning)
    }else{
      warning(sprintf('Model %s is empty.', '<curent model name>'))
    }
  })
}


topVarsByModel <- function(model.list, n=NULL) {
  # !!! INCOMPLETE
  # TODO::
  #   - Incomplete
  #   - Handle/remove the "X" at the beginning of the probeID
  #     - This should have been done while preparing the dataset
  #' Get the top N most important variables for each model in the ensemble
  vI.list <- varImpByModel(model.list)
  vI.list[!sapply(vI.list, is.null)]
  vI.list.rownames <- sapply(vI.list, rownames)
  Reduce(intersect, vI.list.rownames)
}

meanWeightedVarImp <- function(model.output) {
  # !!! INCOMPLETE
  # TODO::
  #   - Average the weighted variable importance across models 
  #     - Use the GLM coefficients to weight varImp of each model
  #     - Sum the weighted vars, and divide by the number of models
  ens <- model.output$model.list
  stack <- model.output$model.stack
  W <- varImp(model.ensemble$model.stack$ens_model)$importance
  # W <- Use the GLM coefficients instead by exctracting with coef().
  X <- varImpByModel(model.ensemble$model.list)
  N <- length(unlist(W))
  U <- sum(W*X)/N
  sapply(rownames(W), function(mod){
    W[mod,'Overall'] * X[[mod]]
  })
}


goEnrichment <- function(model.list, db='hgu133plus2.db') {
  # !!! INCOMPLETE
  do.call(library, list(db))
  strsplit(xx[[1]],'X')[[1]][2]
  xx <- sapply(xx, function(x){
    strsplit(x,'X')[[1]][2]
  })
  annot.out <- select(hgu133plus2.db, xx, c('ENTREZID')) # "SYMBOL","GENENAME", 
  results <- goana(annot.out$ENTREZID) 
  # results <- goana(results, sort= <the $Overall column from varImp> )
  topGO(results)
}
