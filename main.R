library(GEOquery)
library(edgeR)
library(caret)
library(caretEnsemble)

source('methods.R')

USE_LOCAL_DATA = T

# gi50 <- read.csv('data/singledrug/CANCER60GI50.csv', header=T)
data(gi50)

if(USE_LOCAL_DATA){
  gse32474 <- getGEO(filename='data/GSE32474/GSE32474_series_matrix.txt.gz', getGPL = F)
}else{
  gse32474 <- retryOnError(getGEO, list('GSE32474'))[[1]]
} 

ensemble <-
  trainEnsemble(drug='cisplatin', eset=gse32474, gi50=gi50, LCONC=NULL, dataset=NULL,
                preProFunc = function(x,y) {
                  list(x=t(cpm(t(x))), y=y[,'NLOGGI50'])
                },
                esbl.trainControl = trainControl(
                  method = 'cv',
                  number=5,
                  p = 0.90,
                  savePredictions = 'final',
                  classProbs = FALSE,
                  verboseIter = TRUE,
                  returnData=FALSE
                ),
                esbl.tuneList = list(
                  knn = caretModelSpec('knn', tuneGrid=expand.grid(k=c(30,25,22,20,18,16,12,10,8,7,6,5))),
                  glmnet = caretModelSpec(method='glmnet', tuneGrid=expand.grid(alpha=seq(0.2, 0.05, length.out=10), lambda=seq(0.2, 0.05, length.out=10))),
                  gbm = caretModelSpec('gbm', tuneGrid=expand.grid(n.trees=c(300,500,600,800,1000,1200), interaction.depth=c(3,4,6,8,10,12), shrinkage=seq(0.05,0.5, length.out=6), n.minobsinnode=c(2,4,6,8,10,12))), 
                  svmRadial = caretModelSpec('svmRadial', tuneGrid=expand.grid(sigma=seq(8e-04, 2e-04, length.out=20), C=seq(2.5,0.5, length.out=20))),
                  xgbLinear = caretModelSpec('xgbLinear', tuneGrid=expand.grid(lambda=seq(0.01,0.5,length.out=10), alpha=seq(0.01,0.5, length.out=10), eta=seq(0.1,0.9,length.out=10), nrounds=seq.int(100,1000,length.out=10))),
                  glmnet = caretModelSpec(method='glmnet', tuneGrid=data.frame(alpha=seq(0.2, 0.05, length.out=10), lambda=seq(0.2, 0.05, length.out=10)))
                ),
                stack.trainControl = trainControl(
                  method="cv",
                  number=10,
                  savePredictions='final',
                  classProbs=FALSE,
                  summaryFunction=defaultSummary,
                  returnData=FALSE
                ),
                stack.method = 'glm',
                stack.metric = 'RSME'
  )

summary(ensemble)
modelCor(resamples(ensemble))
xyplot(resamples(ensemble))

# GO enrichment 
# TODO:: Use topVarsByModel() & meanWeightedVarImp() instead
ensemble$impVars <- varImpByModel(ensemble$model.list)
ensemble$impGenes <- sapply(ensemble$varImp, function(imps){
  goana(rownames(imps))
})

saveRDS(ensemble, file='ensemble.cisplatin.Rds')

