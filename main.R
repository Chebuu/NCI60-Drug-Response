library(GEOquery)
library(edgeR)
library(caret)
library(caretEnsemble)
source('methods.R')

# TODO:
#   - Consider normalizing with SCAN package from Bioconductor in preProFunc
#   - Consider some feature extration/dimensionality reduction
#   - Consider making an email beeper using mailR to send the results of the run when execution finishes using mailR

USE_LOCAL_DATA = T

data.gi50 <- read.csv('data/singledrug/CANCER60GI50.csv', header=T)
if(USE_LOCAL_DATA){
  gse32474 <- getGEO(filename='data/GSE32474/GSE32474_series_matrix.txt.gz', getGPL = F)
}else{
  gse32474 <- retryOnError(getGEO, list('GSE32474'))[[1]]
} 


drugNames <- cbind('Gefitinib','Gemcitabine','Cisplatin','Carboplatin','Doxorubicin','Docetaxel','Paclitaxel')
drugs <- sapply(drugNames, function(drugName) extractDrugNSC(drugName))
model.ensemble <-
  trainEnsemble(drug='Cisplatin', eset=gse32474, gi50=data.gi50, LCONC=NULL, dataset=NULL,
                preProFunc = function(x,y) {
                  # Do feature engineering here if desired
                  # @ param {x} <Biobase::ExpressionSet> The *transposed* probe count matrix. The matrix has been transposed to match the canonical data table format where ROW x COLUMN = SAMPLES x FEATURES.
                  out <- list(x=x, y=y[,'NLOGGI50'], foo='Hello World!')
                  # out <- list(x=t(cpm(t(x))), y=y[,'NLOGGI50'])
                  return(out)
                },
                esbl.trainControl = trainControl(
                  method = "cv",
                  number=3,
                  p = 0.90,
                  savePredictions = 'final',
                  classProbs = FALSE,
                  verboseIter = TRUE,
                  returnData=FALSE
                ),
                esbl.tuneList = list(
                  glmnet = caretModelSpec(method='glmnet', tuneGrid=expand.grid(alpha=seq(0.2, 0.05, length.out=10), lambda=seq(0.2, 0.05, length.out=10))),
                  knn = caretModelSpec('knn', tuneGrid=expand.grid(k=c(30,25,22,20,18,16,12,10,8,7,6,5))),
                  gbm = caretModelSpec('gbm', tuneGrid=expand.grid(n.trees=c(300,500,600,800,1000,1200), interaction.depth=c(3,4,6,8,10,12), shrinkage=seq(0.05,0.5, length.out=6), n.minobsinnode=c(2,4,6,8,10,12))), #  Fitting n.trees = 600, interaction.depth = 6, shrinkage = 0.064, n.minobsinnode = 8 on full training set
                  svmRadial = caretModelSpec('svmRadial', tuneGrid=expand.grid(sigma=seq(8e-04, 2e-04, length.out=20), C=seq(2.5,0.5, length.out=20))), # Fitting sigma = 0.000579, C = 1.76 on full training set
                  xgbLinear = caretModelSpec('xgbLinear', tuneGrid=expand.grid(lambda=seq(0.01,0.5,length.out=10), alpha=seq(0.01,0.5, length.out=10), eta=seq(0.1,0.9,length.out=10), nrounds=seq.int(100,1000,length.out=10))) # Best nrounds=810, lambda=039476 alpha = 0.158, eta=0.732
                ),
                stack.trainControl = trainControl(
                  method="cv",
                  number=5,
                  # tuneLength = 5,
                  savePredictions="final",
                  classProbs=FALSE,
                  summaryFunction=defaultSummary,
                  returnData=FALSE
                ),
                stack.method = 'glm',
                stack.metric = 'RSME'
  )

saveRDS(model.ensemble, file='model.cisplatin.Rds')

model.ensemble$model.list
lapply(model.ensemble$model.list, length)

varImp(model.ensemble$model.list)
varImpGenes <- varImpByModel(model.ensemble$model.list)

## NO -- See topVarsByModel & meanWeightedVarImp
varImpGenes.ana <- lapply(varImpGenes, function(modl){
  goana(rownames(varImpGenes[[modl]]))  
})
do.call(rbind, varImpGenes.ana)
varImpGenes.ana.unique <- varImpGenes[varImpGenes.ana[1,] == unique(varImpGenes.ana[1,]),]
varImpGenes.ana.unique



# model.ensemble$model.stack$error
# modelCor(resamples(model.ensemble$model.list))
# summary(model.ensemble$model.stack)
# 
# for(model.name in names(model.ensemble$model.list)){
#   model <- model.ensemble$model.list[[model.name]]
#   print(model.name)
#   print(model$bestTune)
# }
# model.ensemble$model.list$glmnet$bestTune
# xyplot(resamples(model.ensemble$model.list))
# modelCor(resamples(model.ensemble$model.list))
# model.ensemble$model.stack$error # RMSE 0.20, Rsq .90
# varImp(model.ensemble$model.list)
# rownames(model.ensemble$model.stack.testData)
# plot(model.ensemble$model.stack.testPredictions,model.ensemble$model.stack.testData$y)

# timestamp()
# models.ensemble.drugs <- sapply(drugNames, function(drugName) {
#   trainEnsemble(drug=drugs[[drugName]], eset=gse32474, gi50=data.gi50, dataset=NULL,
#                 preProFunc = function(x,y) {
#                   out <- list(x=edgeR::cpm(x), y=y[,'NLOGGI50'])
#                   return(out)
#                 },
#                 esbl.trainControl = trainControl(
#                   method = "cv",
#                   number=5,
#                   p = 0.90,
#                   savePredictions = 'final',
#                   classProbs = FALSE,
#                   verboseIter = TRUE,
#                   returnData=FALSE
#                 ),
#                 esbl.tuneList = list(
#                     glmnet = caretModelSpec(method='glmnet', tuneGrid=data.frame(alpha=seq(0.2, 0.05, length.out=10), lambda=seq(0.2, 0.05, length.out=10))),
#                     knn = caretModelSpec('knn', tuneGrid=data.frame(k=c(30,25,22,20,18,16,12,10,8,7,6,5))),
#                     gbm = caretModelSpec('gbm', tuneGrid=data.frame(n.trees=c(300,500,600,800,1000,1200), interaction.depth=c(3,4,6,8,10,12), shrinkage=seq(0.05,0.5, length.out=6), n.minobsinnode=c(2,4,6,8,10,12))), #  Fitting n.trees = 600, interaction.depth = 6, shrinkage = 0.064, n.minobsinnode = 8 on full training set
#                     svmRadial = caretModelSpec('svmRadial', tuneGrid=data.frame(sigma=seq(8e-04, 2e-04, length.out=20), C=seq(2.5,0.5, length.out=20))), # Fitting sigma = 0.000579, C = 1.76 on full training set
#                     xgbLinear = caretModelSpec('xgbLinear', tuneGrid=data.frame(lambda=seq(0.01,0.5,length.out=10), alpha=seq(0.01,0.5, length.out=10), eta=seq(0.1,0.9,length.out=10), nrounds=seq.int(100,1000,length.out=10))) # Best nrounds=810, lambda=039476 alpha = 0.158, eta=0.732
#                     # Broken xgbTree = caretModelSpec('xgbTree',  tuneGrid=data.frame(eta=seq(0.1,0.9,length.out=5), max_depth=c(1,2,3,4,5), gamma=seq(0.00001,0.9,length.out=5), colsample_bytree=seq(0.5,0.9, length.out=5), min_child_weight=seq(0.5,2,length.out=5), subsample=seq(0.5,0.9,length.out=5))) # eta=0.3, max_depth= 1, gamma=0, colsample_bytree=0.6, min_child_weight=1, subsample=seq(0.5,1,length.out=5, nrounds=500,1000,length.out=5 .... nrounds 200,  max_depth 1, eta 0.3, gamma 0, colsample_bytree 0.8 min_child_weight 1 subsample 0.76
#                     # Works  xgbTree = caretModelSpec('xgbTree',  tuneLength=10)
#                     # svmSpectrumString = caretModelSpec('svmSpectrumString', tuneLength=10), 
#                     # svmLinear = caretModelSpec('svmLinear', tuneGrid=data.frame(C=seq(0.25, 2, length.out=10))), 
#                     # svmPoly = caretModelSpec('svmPoly', tuneGrid=data.frame(degree=c(1,2,3,4), scale=seq(1,1e-04,length.out=5), C=c(1,seq(.5,2.5,length.out=4)))), 
#                     # svmRadialCost = caretModelSpec('svmRadialCost', tuneGrid=data.frame(C=seq(.5,2.5,length.out=5)))   
#                   ),
#                 stack.trainControl = trainControl(
#                   method="cv",
#                   number=20,
#                   # tuneLength = 5,
#                   savePredictions="final",
#                   classProbs=FALSE,
#                   summaryFunction=defaultSummary,
#                   returnData=FALSE
#                 ),
#                 stack.method = 'glm',
#                 stack.metric = 'RSME'
#   )
# })
# timestamp()
# for(drugName in drugNames[3:6]){
#   print(drugName)
#   print(models.ensemble.drugs[,drugName]$model.stack$error[,-(1)])
#   for(model in models.ensemble.drugs[,drugName]$model.list) {
#     print(model$bestTune)
#   }
# }
# 
# models.ensemble.drugs[,'Docetaxel']$model.stack$models$glmnet$results
