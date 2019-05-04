library(testit)
library(ChemmineR)
library(edgeR)
library(caret)
library(caretEnsemble)

PATH_TO_AOD8_SDF <- 'data/singledrug/AOD8_Structures.sdf'
REMOTE_URL_AOD8_SDF <- 'https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf' # 'https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf?version=1&modificationDate=1510673147932&api=v2'

# TODO: 
#   - ! Eliminate the globals PATH_TO_AOD8_SDF and REMOTE_URL_AOD8_SDF. Just include AOD8_Structures.sdf as an .Rds object in the data directory of your package or create a data package (the size of datasets might be over the limit for a single package).
#   - ! There is no need to subset the columns of the GI50 dataset in extractResponseByNSC(). You are removing information that the user might want. Change the default col.out argument to be col.out = colnames(gi50)
#   - LCONC = NULL does not work when methods are called within function scopes. To fix this, I've added if(is.null(LCOC)) LCONC <- NULL to every function that has an LCONC parameter. This is less than ideal, but works. Needs a fix.
#   - Think of a way to use all LCONC values for regression. For classification, this could be done by including LCONC as a feature (quick thought).
#   - Make tests to include in package that confirm there are no mistakes in the data parsing
#   - extractDrugNSC() should write the downloaded file to the package data/ directory and and name it the same as PATH_TO_AOD8_SDF (but what if someone changes PATH_TO_AOD_SDF?).
#   - Include the LCONC value that is being used for each model in the output of trainEnsemble, because it matters.
#   - Redundant error checking in makeDataset
#   - Many copies of data in trainEnsemble, find better solution
#   - Hard-coded model.list arguments in trainEnsemble
#   - Should you just delete the POST PROC block in trainEnsemble? The user can do that themself.
#   - Refactor hasty/hacky solutions
#     - tryCatch in extractDrugNSC ==> check if the user passed an sdfSet argument, otherwise the "Trying local file..." will run twice
#     - createDataPartition and test.rows in trainEnsemble

argmin <- function(x){
  # Returns the index of the minimum value in the vector x
  which(x == min(x))
}

argmax <- function(x){
  # Returns the index of the maximum value in the vector x
  which(x == max(x))
}

retryOnError <- function(func, args, n.tries=Inf) {
  # TODO:
  #   - Is tempenv destroyed after the function terminates? Should I even worry about removing it?
  # Repeatedly call a function until the call does not throw an error (use case is network requests that return 404 and throw an error).
  # Returns the result of do.call(func, args)
  # @ param {func} A function definition to be called using the argument to args.
  # @ param {args} A list of arguments to be passed to func. Named lists are accepted wherein name specify parameters. 
  # @ param {n.tries} An integer specifying the number of times do.call(func, args) should be allowed to throw an error before giving up and returning FALSE.
  # @ example retryOnError(do.call, list(c,1)) # NOT RUN Executing do.call(c,1) will throw an error, because the second argument to do.call must be a list. This will cause retryOnError to enter an infinite loop, as expected.
  # @ example retryOnError(getGEO, list('GSE32474'))
  tempenv <- new.env()
  tempenv$triesLeft <- n.tries
  tempenv$success <- T
  errorCallBack <- function(e){
    print(e); print('Retrying...') 
    tempenv$triesLeft <- tempenv$triesLeft-1
    tempenv$success <- F
    return(F)
  }
  repeat{
    res <- tryCatch(do.call(func, args), error=errorCallBack)
    if(tempenv$success | tempenv$triesLeft==0){
      if(tempenv$success) {
        print('Success.')
      }else {
        print(sprintf('Maximum tries (%s) reached.', n.tries))
      }
      break
    }
  }
  remove(tempenv)
  return(res)
}

# DEL
# handleAOD8 <- function(sdfSet=NULL, localURL=PATH_TO_AOD8_SDF, remoteURL=REMOTE_URL_AOD8_SDF, downloadDest=NULL) {
#   # if(!is.character(sdfSet)) warning('The sdfSet argument should be a path to a .sdf file (AOD8_Structures.sdf) or an SDFset object.')
#   sdfSet <- tryCatch(read.SDFset(sdfSet), error=function(e){
#     warning(sprintf('%s\n sdfSet not found. Trying local AOD8_Structures.sdf ...', e))
#     sdfSet <- tryCatch(readSDFset(PATH_TO_AOD8_SDF), error=function(e) {
#       warning(sprintf('%s\nLocal AOD8_Structures.sdf not found. Trying to download it...', e))
#       x <- tryCatch(getURL(REMOTE_URL_AOD8_SDF), error=function(e) {
#         # stop(sprintf('%s Could not download AOD8_Structures.sdf. Check REMOTE_URL_AOD8_SDF and make sure the file exists at the specified URL. Exiting.',e))
#         stop(sprintf('%s\nCould not download AOD8_Structures.sdf. Please download it from https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf and pass its path as the argument to sdfSet in extractDrugNSC(drugname, sdfSet).',e))
#       })
#       tmpf <- tempfile() # This is not a very elegant solution, because /var will get bloated with .sdf files if the user runs this multiple times. Would be better to write to a static directory where the file is expected by this function to exist. http://r-pkgs.had.co.nz/data.html
#       writeLines(x, tempf)
#       sdfSet <- readSDFset(tmpf)
#       return(sdfSet)
#     })
#     return(sdfSet)
#   })
#   return(sdfSet)
# }
# END DEL

extractDrugNSC <- function(drugname, sdfSet=PATH_TO_AOD8_SDF){
  # TODO:: 
  #   - In download block, dont write to temp file. See comment.
  #   - Too many sdfSet variables in nested tryCatch? If it works it works. 
  #   - Change param sdfSet to SDFset to match the object class name?
  #   - Consider handling the reverse situation??? (if argument is NSC ID find a corresponding name).
  #   - Make search more robust. Its pretty greedy with the grepl statement, but it works fine for now. Check Chemminr methods for this. 
  # Match a drug name with its NSC identifier by searching for matches in each SDF's @datablock slot
  # @ param [sdfSet]   <character> The path to the AOD8_Structures.sdf file available for download at https://wiki.nci.nih.gov/display/NCIDTPdata/ . Alternatively, SDFset objects are also accepted.
  # @ param [drugname] <character> The name of the drug for which an NSC identifier should be found
  # @ example extractDrugNSC('cisplatin', PATH_TO_AOD8_SDF)
  
  if(any(class(sdfSet) == 'SDFset')){
    # The user passed an SDFset object.
    TRUE # Break
  }else{
    # We need to find the AOD8_Structures.sdf file and read it
    sdfSet <- tryCatch(read.SDFset(sdfSet), error=function(e){
      warning(sprintf('%s\n sdfSet not found. Trying local AOD8_Structures.sdf...', e))
      sdfSet <- tryCatch(readSDFset(PATH_TO_AOD8_SDF), error=function(e) {
        warning(sprintf('%s\nLocal AOD8_Structures.sdf not found. Trying to download it...', e))
        x <- tryCatch(getURL(REMOTE_URL_AOD8_SDF), error=function(e) {
          stop(sprintf('%s\nCould not download AOD8_Structures.sdf. Please download it from https://wiki.nci.nih.gov/download/attachments/160989212/AOD8_structures_corrected.sdf and pass its path as the argument to sdfSet in extractDrugNSC(drugname, sdfSet).',e))
        })
        tmpf <- tempfile() # This is not a very elegant solution, because /var will get bloated with .sdf files if the user runs this multiple times. Would be better to write to a static directory where the file is expected by this function to exist. http://r-pkgs.had.co.nz/data.html   is there a system.file() for writing files to packages?
        writeLines(x, tempf)
        sdfSet <- readSDFset(tmpf)
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
      # The user is not expected to pass a number. 
      # Just return the number in-case it is the NSC ID itself.
      # Could do a reverse lookup for the drug name here and add it as names(out)...
      out <- drugname
    }else{
      stop(sprintf('An NSC identifier for %s couldn\'t be found. If you know this drug\'s NSC ID, please use that instead' , drugname))  
    }
  }
  
  return(out)
}
# extractDrugNSC('Gefitinib')
# extractDrugNSC(715055)

grepCellLine <- function(phenoTitle, reformat=T){
  # TODO:
  #   - I dont think sapply is necessary. I think you can just do out <- gsub('','',out)
  # Extracts cell name and returns a reformatted version. This is intended to be a lambda function applied to a list of cell names, such as that found in pData(eset)$title
  # @ param [phenoTitle] <character> A cell name to be reformatted
  # @ param [remormat]   <logical>   Should the cell name be reformatted (T)? Or should it simply be extracted and returned (F)? The default is highly recommended. There is almost no reason not to reformat.
  # @ example sapply(pData(getGEO('GSE32474')[[1]])$title, grepCellLine)
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
  # NOT RUN
  # I thought I might search through all the lists in pData(eset) and see if any contain matches to cell names in GI50, but this seems unecessary, especially because the user doesn't really need to worry about this.
  # The $title attribute for GSE32474 will always exist in this package.
  # Maybe I should rename this to findTitleAttribute(). Really, I'm just looking for any names(pData(eset)) that I can use to extract cell names.
  # See extractCellsFromEset() for use case
  # Another use case exists in cellNames()
}

extractCellsFromEset <- function(eset) {
  # The $title attribute in the @phenoData slot of GSE32474 contains cell names that don't exactly match those in the GI50 dataset. This function reformats the GSE32474 cell names to match the GI50 dataset and returns a list of the reformatted cell names. The returned value is expected to be used to set pData(GSE32474)$cell <- returned_value for downstream funtionality. If that operation has already been done, the function just returns pData(GSE32474)$cell
  # @ param [eset] <ExpressionSet> An object containing a @phenoData slot with a $title attribute that holds a list of cell names
  # @ example gse32474 <- getGEO('GSE32474'); pData(gse32474)$cell <- extractCellsFromEset(gse32474);
  if(is.list(eset)) eset <- eset[[1]]
  if(is.null(pData(eset)$cell) | length(pData(eset)$cell) == 0){
    if(is.null(pData(eset)$title)) {
      # Here is the findCellAttr() use case
      stop('Unable to extract cell names because pData(eset)$title does not exist.')
    }
    return(sapply(pData(eset)$title, grepCellLine))
  }else{
    return(pData(eset)$cell)
  }
}

cellNames <- function(eset) {
  # TODO:: 
  #   - Write more descriptive warning messages
  # Same as pData(eset)$cell. More robust matching in case $Cell or $CELL exists but $cell doesn't
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
      # FALSE
      idx <- 'title'
      # idx <- findCellAttr()
    }
  }else{
    out <- pData(eset)[[idx]]  
  }
  return(out)
}

extractResponseByNSC <- function(nsc, gi50, col.out = c('NSC','NLOGGI50', 'CELL', 'LCONC', 'STDDEV')) {
  # Given an NCS drug identifier, return the corresponding rows of the GI50 dataset for all cells receiving that drug. The function also reformats the items in the column gi50$CELL
  # @ param [ncs]      <integer>  A single NCS identifier
  # @ param [gi50]     <matrix>   A matrix or dataframe of the GI50 dataset, the first column of which contains NCS identifiers
  # @ param [col.out]  <vector>   A list of gi50 columns to return. Downstream methods require at least c('NSC', 'NLOGGI50') 
  if(!is.numeric(nsc)) nsc <- extractDrugNSC(nsc)
  out <- gi50[which(gi50[,'NSC']==nsc), col.out]
  if(!is.null(out$CELL)){
    out$CELL <- sapply(out$CELL, function(x) {
      x <- sub(' ', '', x)                       # Remove the first occurence of a space
      # There are outlier names that would otherwise be thrown away without the conditional blocks below
      if(grepl('ADR-RES', x)) x <- 'NCI-ADR-RES' # An outlier "NCI/ADR-RES" ==> "NCI-ADR-RES"
      if(grepl('T-47D', x)) x <- 'T47D'          # Another outlier
      if(grepl('RXF393', x)) x <- 'RXF-393'      # Another outlier
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
  #   - Currently LCONC = 'all' is not supported because I have not come up with an idea for implementing all the LCONC values in a regression model. Currently, if the user passes LCONC='all' the funciton treats this as if LCONC=NULL (the LCONC value that yields the greatest number of samples).
  # @ param [responses] The *subsetted* responses by drug. Note, this should not be the GI50 datset. Instead, it should be the GI50 dataset subsetted by the drug of interest.
  # @ param [eset] The complete ExpressionSet object (GSE32474)
  # @ param [LCONC] Either a numeric value representing the desired log(concentration) of drug exposure that will be used to fit the model, or a string being either 'min' 'mid' 'max' to select the minimun, median, or max value of all LCONC values found in the GI50 dataset subsetted by drug. If a number is specified, and that number is not found, it will try to find the number with smallest RSE between the given value and those that exist in the subsetted dataset. If left null, the highest concentration of drug available is used.
  # DEL
    # drug <- extractDrugNSC('Paclitaxel')
    # responses <- extractResponseByNSC(drug, data.gi50)
    # eset <- gse32474
    # LCONC <- NULL
  # END

  pData(eset)$cell <- extractCellsFromEset(eset)
  responses.cellNames <- responses$CELL
  eset.cellNames <- cellNames(eset) # Same as pData(eset)$cell. More robust matching in case the $Cell or $CELL exists but $cell doesn't
  eset.cellNames.responseMatches <- match(eset.cellNames, responses.cellNames)
  
  # Assertions
  all(eset.cellNames == responses[eset.cellNames.responseMatches, 'CELL']) 
  all(eset.cellNames == responses.cellNames[eset.cellNames.responseMatches])
  all(responses$CELL[eset.cellNames.responseMatches] == eset.cellNames)
  
  # Convenience variables to handle LCONC
  responses.lconc <- responses[eset.cellNames.responseMatches, 'LCONC']
  responses.lconc.levels <- as.numeric(levels(as.factor(responses.lconc)))
  responses.lconc.LCONC <- which(responses.lconc == LCONC) 
  
  # For the message() statement
  notFound <- list(b=FALSE, val=numeric())
  
  # This conditional block handles the user's LCONC argument. Several checks are performed as well as handling the string options. The result is reassignment of LCONC if the checks catch error or a string option is passed., otherwise LCONC is not changed.
  if(!is.null(LCONC)){ # Handle the user's argument to LCONC
    if(!suppressWarnings(is.na(as.numeric(LCONC)))) LCONC <- as.numeric(LCONC) # Make sure they didn't accidentally pass a valid LCONC number as a string. If they did, convert it for them and carry on.
    if(is.character(LCONC)){
      LCONC <- tolower(LCONC) # In case the user capitalized something on accident, we know what they mean.
      if(LCONC == 'all'){ 
        # Handle select all, this requires some refactoring to the project as a whole. Currently not implemented. 
        samples.counts <- sapply(responses.lconc.levels, function(x){
          # Count the number of instances of each factor level in the subsetted dataset
          sum(responses.lconc == x)
        })
        LCONC <- responses.lconc.levels[argmax(samples.counts)] # Select responses using the LCONC value with the highest frequency 
        warning(sprintf('The argument LCONC = "all" is not currently supported. The LCONC value that yields the maximum number of samples (LCONC=%s).', LCONC))
      }else if(LCONC == 'mid'){
        # Pick the middle number ie. the median in a sorted, unique list. I will err on the side of lower concentration in case there are an even number of options.
        which.idx <- floor((length(responses.lconc.levels)+1)/2) # The variable responses.lconc.levels is already sorted ascending. Just pick the middle number. Using median() will not work well for decimals, hence this solution.
        LCONC <- responses.lconc.levels[which.idx] 
      }else{
        # The user passed min or max, can just use the user's argument to call the function (already converted tolower)
        LCONC <- do.call(LCONC, list(responses.lconc.levels))
      }
    }else if(is.numeric(LCONC) & length(responses.lconc.LCONC) == 0){ 
      # Couldnt find that LCONC value in the subsetted data
      LCONC.user<- LCONC
      if(-LCONC %in% responses.lconc){
        warning(sprintf('LCONC value %s not found for drug %s. Found -%s value. Using -%s instead.', LCONC, drug, LCONC, LCONC))
        LCONC <- -LCONC
      }else{ 
        # Use the closest LCONC to the users argument that exists in the subsetted datset
        LCONC <- responses.lconc.levels[argmin(sqrt((LCONC - responses.lconc.levels)^2))]
        warning(sprintf('LCONC value %s not found for drug %s. The nearest LCONC value avilable in the dataset (%s) has been selected instead.', LCONC.previous, drug, LCONC))
      }
      notFound$b <- T; notFound$val <- LCONC.user # Only used for the message
    }else{
      # Do nothing, the user's LCONC argument is sufficient to subset the data 
      TRUE 
    }
  }else{ 
    # Choose the LCONC value that will yield the most samples
    samples.counts <- sapply(responses.lconc.levels, function(x){ 
      # Count the frequency of each factor in the subsetted dataset
      sum(responses.lconc == x)
    })
    # Select LCONC using the factor level with the most instances (the LCONC that yields the most samples)
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
  
  # How many samples were eliminated after subsetting by LCONC?
  n.eliminated <- nrow(responses.match) - nrow(responses.match.LCONC)
  perc.eliminated <- round(n.eliminated / nrow(responses.match) * 100, digits=2)
  drug <- as.factor(levels(responses.match[,'NSC']))
  if(n.eliminated == 0){
    if(notFound$b == F){
      msg <- sprintf('Subsetting the dataset by LCONC value %s yielded the maximum number of samples available (%s) ie. only one LCONC value was available for drug %s.', LCONC, nrow(responses.match), drug)
    }else{
      msg <- sprintf('Subsetting the dataset by LCONC value %s yielded the maximum number of samples available (%s).', LCONC, nrow(responses.match), drug)
    }
    message(msg)
  }else{
    message(sprintf('Subsetting the dataset by LCONC value %s eliminated %s%% of the dataset (%s samples eliminated out of %s total samples available). If you would like to use the maximum number of samples possible, please set LCONC=NULL .', LCONC, perc.eliminated, n.eliminated, nrow(responses.match)))  
  }
  
  # Return the dataset subsetted by LCONC and the corresponding columns of the eset to use
  return(list(responses = responses.match.LCONC, x.col=x.col, LCONC=LCONC))
  
}

# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC=LCONC)
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='all')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='min')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='max')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='mid')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='1')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='-6')
# matchResponseToGSMByName(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC=1.0)

makeDataset <- function(responses, eset, LCONC=NULL, preProFunc=function(x){return(x)}) {
  # TODO:
  #   - Should you do an error check for the dims of X and Y? Or should you leave as-is and just rely on the downstream check within trainEnsemble(), thereby allowing the user to create a dataset with dims that don't match?
  # @ param [preProFunc] <function> A lambda function for preprocessing with parameters x and y, where x is a count matrix and y is the GI50 dataset (a dataframe). It should be of the form function(x){return(x)} or function(x,y){return(list(x=x,y=y))}. If using the y parameter, do not return the full GI50 dataframe, only return a vector of response variables as the $y attribute in your returned list.
  
  # DEL
    # drug <- extractDrugNSC('Doxorubicin')
    # gi50 <- data.gi50
    # responses <- extractResponseByNSC(drug, gi50)
    # eset <- gse32474
    # preProFunc <- function(x){return(list(x=x,extraStuff=list('hello')))}
  # END
  
  responses.matched <- matchResponseToGSMByName(responses, eset, LCONC=LCONC)
  x.col <- responses.matched$x.col # Use these columns of the count matrix to be sure they match the Y-values created by matchResponseToGSMByName()
  x <- t(Biobase::exprs(eset)[,x.col]) # The count matrix in the ExpressionSet object [features,samples] is the transpose of the format expected by caret [samples,features]
  y <- responses.matched$responses # The Y-values created by matchResponseToGSMByName() are stored in $responses
  LCONC <- responses.matched$LCONC # The LCONC value selected by the user may not be the same as that returned by matchResponseToGSMByName() due to error checking
  preProFunc.extra <- NULL # This var will hold the additional information returned by the user's preProFunc
  
  # Assertions
  all(eset$cell[x.col] == y$CELL) # Make sure X and Y will match up ie. they both refer to the same cells
  
  # Handle preProFunc
  if(class(preProFunc) == 'function'){
    preProFunc.nargs <- length(as.list(args(preProFunc)))
    if(preProFunc.nargs == 2){
      # The user expects only x to be passed to preProFunc (the default)
      out.ppf <- preProFunc(x)
      if(is.list(out.ppf)){
        # The user returned a list from preProFunc
        names(out.ppf) <- tolower(names(out.ppf))
        x <- out.ppf$x
        if(length(names(out.ppf)) > 1){
          # Assumes the returned list contains more than one name, and that one of those names is "x"
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
      warning('Your preProFunc was not applied to the dataset, please check its parameters. The argumen to preProFunc should be a function with 1 or 2 parameters of the form function(x){return(x)} or function(x,y){return(list(x=x,y=y))}).')
      warning('Skipping preprocessing...')
      x <- x
      y <- y$NLOGGI50
    }
  }else{
    warning('Your preProFunc was not applied to the dataset. The argument to preProFunc should be a function with 1 or 2 parameters of the form function(x){return(x)} or function(x,y){return(list(x=x,y=y))}).')
    warning('Skipping preprocessing...')
    x <- x
    y <- y$NLOGGI50
  }
  return(list(x=x, y=y, LCONC=LCONC, preProFunc.extra=preProFunc.extra))
}
# makeDataset(extractResponseByNSC(extractDrugNSC('Paclitaxel'), data.gi50), gse32474, LCONC=NULL)
# makeDataset(extractResponseByNSC(extractDrugNSC('DOESNT_EXIST'), data.gi50), gse32474, LCONC=NULL)
# makeDataset(extractResponseByNSC(extractDrugNSC('Gefitinib'), data.gi50), gse32474, LCONC='100')


trainRFE <- function(drug, eset, gi50, preProFunc=function(x) {return(x)},
                     rfeControl=rfeControl(), trainControl=trainControl(),
                     method='svmLinear', multiFolds=2, trainPercent=.85, 
                     sizes=seq(nrow(eset)/3, nrow(eset)-1, length.out = 10), 
                     ...) 
{
  # TODO::
  #   Can I move some of these params to trainControl or rfeControl?
  #   Default sizes should be converted to integer for good measure
  
  if(!is.numeric(drug)) {
    if(file.exists(PATH_TO_AOD8_SDF)){
      drug <- extractDrugNSC(drug, read.SDFset(PATH_TO_AOD8_SDF))  
    }else{
      # TODO: consider a remote query if an SDF set is not saved locally, or consider downloading it. 
      stop('Argument drug should be an NSC drug identifier (integer). Please use extractDrugNSC() to find your drug\'s NSC identifier.')
    }
  }
  if(is.list(eset)) eset <- eset[[1]]
  if(!is.integer(drug)) 
    if(is.null(pData(eset)$cell)) pData(eset)$cell <- extractCellsFromEset(eset)
    
    responses.nsc <- extractResponseByNSC(drug, gi50)
    
    data <- makeDataset(responses.nsc, eset)
    if(length(data$y) != nrow(data$x)) stop('x and y lengths of dataset differ.')
    
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

# ###########
# # EXAMPLE #
# ###########
# caretFuncs$summary <- function(...) c(defaultSummary(...))
# svm.RFE <- trainRFE(613327, gse32474, data.gi50,
#                     function(x) {return(x)},
#                     rfeControl = rfeControl(
#                       method = "repeatedcv",
#                       repeats = 2,
#                       number = 5,
#                       verbose = TRUE,
#                       functions = caretFuncs,
#                       seeds = set.seed(123),
#                       saveDetails = TRUE
#                     ),
#                     trainControl = trainControl(
#                       method = "LGOCV",
#                       number = 1,
#                       p = 0.8,
#                       seeds = set.seed(321),
#                       savePredictions = TRUE,
#                       classProbs = FALSE,
#                       verboseIter = TRUE
#                     ),
#                     multiFolds = 2, 
#                     sizes = seq(nrow(gse32474)/2, nrow(gse32474)-1, 10000), trainPercent = 0.85
#                   )
# pred <- svm.RFE$pred
# result <- svmRFE$results
# fit <- svmRFE$fit
# ###############
# # END EXAMPLE #
# ###############

trainEnsemble <- function(drug=NULL, eset=NULL, gi50=NULL, dataset = NULL, LCONC=NULL, preProFunc=function(x) {return(x)},
                          esbl.trainControl=trainControl(), esbl.tuneList=list(svm=caretModelSpec(method='svm', tuneLength=3)), esbl.methodList=c(), 
                          stack.trainControl=trainControl(), stack.method='glm', stack.metric='RMSE', train.percent=0.85,
                          do.plot = T, ...) 
{
  # TODO:
  #   - Do you need to partition the data before passing to caret? Can't you just get caret to split the dataset and return the data? That way the user's resampling indices will work. I think currently samples might be witheld from training that don't need to be held out...? 
  #   - Don't put the Rsq and RMSE in the plot title. Use text() or something (except this throws error in Rmarkdon (Error:plot.new() has not been called)). Find a fix.
  # @ param [drug]       <numeric>  An integer specifying the NSC ID of the drug for which an ensemble will be trained. The ID can be found for a given drug name using extractDrugNSC().
  # @ param [eset] <ExpressionSet> An ExpressionSet containing a count matrix and a @phenoData slot with a $title attribute that holds a list of cell names. eg. getGEO(GSE32474). The count matrix will be passed as parameter x to the lambda function preProFunc=function(x){...}.
  # @ param [gi50] <matrix> A matrix or dataframe of the GI50 dataset with column names of at least c('NLOGGI50', 'CELL').
  # @ param [dataset] <matrix> A user defined dataset. Should be left NULL if drug, eset, and gi50 arguments are passed. preProFunc is not applied to this dataset.
  # @ param [preProFunc] <function> A lambda function for preprocessing. It should be of the form function(x){return(x)} or function(x,y){return(list(x=x,y=y))} where x is a count matrix and y is the GI50 dataset (a dataframe). If using y, return the full dataframe.
  # @ param [esbl.trainControl] <list> The list created by a call to trainControl() which is passed as the trControl argument to caretEnsemble::caretList
  # @ param [esbl.methodList] <list> Alternative to esbl.tuneList passed as the methodList argument to caretEnsemble::caretList(). Can be left NULL if esbl.trainControl argument is passed
  
  # DEL
  # drug=613327
  # eset=gse32474
  # gi50=data.gi50
  # dataset=NULL
  # preProFunc = function(x,y) {
  #   out <- list(x=edgeR::cpm(x), y=y[,'NLOGGI50'])
  #   return(out)
  # }
  # esbl.trainControl = trainControl(
  #   method = "cv",
  #   number=5,
  #   # repeats=2,
  #   p = 0.90,
  #   # index = resampleIndex,
  #   savePredictions = 'final',
  #   classProbs = FALSE,
  #   verboseIter = TRUE,
  #   returnData=FALSE
  # )
  # esbl.tuneList=list(
  #   # treebag=caretModelSpec('treebag'), CRASHES
  #   knn=caretModelSpec('knn', tuneGrid=data.frame(k=c(4,3,2))),
  #   # gbm=caretModelSpec('gbm', tuneGrid=data.frame(n.trees=c(100,200,500,1000), interaction.depth=c(2,4,6,8), shrinkage=shrinkage=c(0.01,0.05,0.1,0.2), n.minobsinnode=c(2,5,10,20))),
  #   svmRadial=caretModelSpec('svmRadial', tuneGrid=data.frame(sigma=seq(1e-03, 1e-04, length.out=8), C=seq(3,1, length.out=8))),
  #   glmnet=caretModelSpec(method='glmnet', tuneGrid=data.frame(alpha=seq(0.05, 0.001, length.out=5), lambda=seq(0.01, 0.001, length.out=5)))
  # )
  # stack.trainControl = trainControl(
  #   method="boot",
  #   number=1,
  #   # tuneLength = 5,
  #   savePredictions="final",
  #   classProbs=FALSE,
  #   summaryFunction=defaultSummary,
  #   returnData=FALSE
  # )
  # stack.method = 'glm'
  # stack.metric = 'RSME'
  # END
  
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
  
  # Assertions
  assert('The test set does not contain any samples found in the train set.', !any(test.rows %in% train.rows))
  
  # Partition the dataset
  data.train <- list(y=data$y[train.rows], x=data$x[train.rows,])
  data.test <- list(y=data$y[test.rows], x=data$x[test.rows,])
  
  # remove(data) # Save some space? Please stop crashing 5 drugs deep...
  
  
  # ENSEMBLE PARAMETERS
  esbl.list.args <- list(
    x = data.train$x,
    y = data.train$y,
    trControl=esbl.trainControl,
    continue_on_fail = T
  )
  if(length(esbl.methodList)>0){ # Does the user want to use methodList
    esbl.list.args$methodList <- esbl.methodList
  }else{ # Or tune list
    esbl.list.args$tuneList <- esbl.tuneList
  }
  
  # TRAIN ENSEMBLE
  message(sprintf('Training ensemble with %s samples and %s features.', nrow(esbl.list.args$x), ncol(esbl.list.args$x)))
  model.list <- do.call(caretList, esbl.list.args) # Train the model
  
  # TRAIN STACK
  model.stack <- caretStack(
    model.list,
    method = 'glm', # Delete or add param
    metric = 'RMSE', # Delete or add param, I think caretEnsemble selects ROC or RMSE depending on classification or regression
    trControl=stack.trainControl
  )
  
  # POST PROC
  predictions <- predict(model.stack, newdata=data.test$x)
  rmse <- RMSE(predictions, data.test$y)
  rsq <- R2(predictions, data.test$y)
  
  
  if(do.plot){
    plot(data.test$y, predictions, main=sprintf('NSC%s LCONC%s Rsq=%s RMSE=%s', drug, LCONC, round(rsq, digits=3), round(rmse, digits=3)))
    abline(0,1)
  }
  
  out <- list(
    LCONC = LCONC,
    model.list = model.list,
    model.stack = model.stack, 
    model.stack.testRMSE = rmse,
    model.stack.rsq = rsq,
    model.stack.dataTest.y = data.test$y,
    model.stack.predictions = predictions,
    preProFunc.extra = preProFunc.extra
  )
  
  return(out)
}

varImpByModel <- function(model.list) 
{
  # The importance is not averaged like you thought. Need to implement by hand
  # @ example getVarImp(model.ensemble$model.list)
  # lapply(model.list, function(mod){
  #   tryCatch(varImp(mod$finalModel), error=function(e){
  #     # This tryCtch is slow, do something else?
  #     return(NULL)
  #   })
  # })
  sapply(model.list, function(mod){
    if(length(model)>0){
      return(varImp(mod))
    }else{
      return(0)
    }
  })
}

topVarsByModel <- function(model.list, n=NULL) {
  # Not complete, also need to handle the X at the beginning of the probeID
  #DEL
  # model.list <- model.ensemble$model.list
  # n <- 9
  #END
  vI.list <- varImpByModel(model.list)
  vI.list[!sapply(vI.list, is.null)]
  vI.list.rownames <- sapply(vI.list, rownames)
  Reduce(intersect, vI.list.rownames)
}

meanWeightedVarImp <- function(model.output) {
  # Broken, how to I get the coefficients of the GLM?
  # Use the glmcoefficients to weight the varImp of each model, sum the weighted vars, and divide by the number of models
  ens <- model.output$model.list
  stack <- model.output$model.stack
  
  model.ensemble$model.stack$models$glm$finalModel
  
  # Wut? 
  # ccc <- coefficients(model.ensemble$model.stack$models$glm$finalModel)
  # sapply(model.ensemble$model.stack$models, function(mod) coefficients(mod$finalModel))
  # dim(ccc)
  
  
  W <- varImp(model.ensemble$model.stack$ens_model)$importance
  # W <- Use the GLM coefficients instead (can't find them, though).
  X <- varImpByModel(model.ensemble$model.list)
  N <- length(unlist(W))
  
  # sum(W*X)/N
  rownames(W)
  X[3]
  # W['svmRadial','Overall'] * X[['svmRadial']]
  sapply(rownames(W), function(mod){
    W[mod,'Overall'] * X[[mod]]
  })
  

}

goEnrichment <- function(model.list, db='hgu133plus2.db') {
  source("https://bioconductor.org/biocLite.R")
  biocLite("hgu133plus2.db")
  library(annotate)
  library(hgu133plus2.db)
  do.call(library, list(db))
  strsplit(xx[[1]],'X')[[1]][2]
  xx <- sapply(xx, function(x){
    strsplit(x,'X')[[1]][2]
  })
  
  annot.out <- select(hgu133plus2.db, xx, c('ENTREZID')) # "SYMBOL","GENENAME", 
  results <- goana(annot.out$ENTREZID) # do goana(results, sort={the Overall column from varImp})
  topGO(results)
}


# ###########
# # EXAMPLE #
# ###########
# resampleIndex <- createResample(ncol(gse32474), 10)
# trainEnsemble(drug='Paclitaxel', eset=gse32474, gi50=data.gi50, LCONC=-1,
#               preProFunc = function(x,y) {
#                 # Feature elimination/extraction here
#                 
#                 out <- list(x=x, y=y[,'NLOGGI50'], foo='Hello World!')
#                 # out <- list(x=t(cpm(t(x))), y=y[,'NLOGGI50'])
#                 return(out)
#               },
#               esbl.trainControl = trainControl(
#                 method = "cv",
#                 number=3,
#                 p = 0.90,
#                 savePredictions = 'final',
#                 classProbs = FALSE,
#                 verboseIter = TRUE,
#                 returnData=FALSE
#               ),
#               esbl.tuneList = list(
#                 glmnet = caretModelSpec(method='glmnet', tuneGrid=data.frame(alpha=seq(0.2, 0.05, length.out=10), lambda=seq(0.2, 0.05, length.out=10))),
#                 knn = caretModelSpec('knn', tuneGrid=data.frame(k=c(30,25,22,20,18,16,12,10,8,7,6,5))),
#                 gbm = caretModelSpec('gbm', tuneGrid=data.frame(n.trees=c(300,500,600,800,1000,1200), interaction.depth=c(3,4,6,8,10,12), shrinkage=seq(0.05,0.5, length.out=6), n.minobsinnode=c(2,4,6,8,10,12))), #  Fitting n.trees = 600, interaction.depth = 6, shrinkage = 0.064, n.minobsinnode = 8 on full training set
#                 svmRadial = caretModelSpec('svmRadial', tuneGrid=data.frame(sigma=seq(8e-04, 2e-04, length.out=20), C=seq(2.5,0.5, length.out=20))), # Fitting sigma = 0.000579, C = 1.76 on full training set
#                 xgbLinear = caretModelSpec('xgbLinear', tuneGrid=data.frame(lambda=seq(0.01,0.5,length.out=10), alpha=seq(0.01,0.5, length.out=10), eta=seq(0.1,0.9,length.out=10), nrounds=seq.int(100,1000,length.out=10))) # Best nrounds=810, lambda=039476 alpha = 0.158, eta=0.732
#                 # Broken xgbTree = caretModelSpec('xgbTree',  tuneGrid=data.frame(eta=seq(0.1,0.9,length.out=5), max_depth=c(1,2,3,4,5), gamma=seq(0.00001,0.9,length.out=5), colsample_bytree=seq(0.5,0.9, length.out=5), min_child_weight=seq(0.5,2,length.out=5), subsample=seq(0.5,0.9,length.out=5))) # eta=0.3, max_depth= 1, gamma=0, colsample_bytree=0.6, min_child_weight=1, subsample=seq(0.5,1,length.out=5, nrounds=500,1000,length.out=5 .... nrounds 200,  max_depth 1, eta 0.3, gamma 0, colsample_bytree 0.8 min_child_weight 1 subsample 0.76
#                 # Works  xgbTree = caretModelSpec('xgbTree',  tuneLength=10)
#                 # svmSpectrumString = caretModelSpec('svmSpectrumString', tuneLength=10),
#                 # svmLinear = caretModelSpec('svmLinear', tuneGrid=data.frame(C=seq(0.25, 2, length.out=10))),
#                 # svmPoly = caretModelSpec('svmPoly', tuneGrid=data.frame(degree=c(1,2,3,4), scale=seq(1,1e-04,length.out=5), C=c(1,seq(.5,2.5,length.out=4)))),
#                 # svmRadialCost = caretModelSpec('svmRadialCost', tuneGrid=data.frame(C=seq(.5,2.5,length.out=5)))
#               ),
#               stack.trainControl = trainControl(
#                 method="cv",
#                 number=5,
#                 # tuneLength = 5,
#                 savePredictions="final",
#                 classProbs=FALSE,
#                 summaryFunction=defaultSummary,
#                 returnData=FALSE
#               ),
#               stack.method = 'glm',
#               stack.metric = 'RSME'
# )
# 
# xyplot(resamples(model.ensemble$model.list))
# modelmodelCor(resamples(model.ensemble$model.list))

# ###############
# # END EXAMPLE #
# ###############
