# ######################################################################
#
##  Meta-Analysis Model Building, Complete pipeline
#
# ######################################################################

# Packages
# Packages
library("tidyverse")
library("SIAMCAT")
library("yaml")

# ##############################################################################
# general 
cat('Starting model building script\n')
start.time <- proc.time()[1]

set.seed(2018)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]


parameters <- yaml.load_file('../parameters.yaml')
# extract parameters
norm.method <- parameters$model.building$norm.method
n.p <- list(log.n0=ifelse(tag %in% c('species', 'genus'),
                          as.numeric(parameters$model.building$log.n0),
                          as.numeric(parameters$model.building$log.n0.func)),
            sd.min.q=as.numeric(parameters$model.building$sd.min.q),
            n.p=as.numeric(parameters$model.building$n.p),
            norm.margin=as.numeric(parameters$model.building$norm.margin))
num.folds <- as.numeric(parameters$model.building$num.folds)
num.resample <- as.numeric(parameters$model.building$num.resample)
ml.method <- parameters$model.building$ml.method
min.nonzero.coeff <- as.numeric(parameters$model.building$min.nonzero.coeff)
modsel.crit <- list(parameters$model.building$modsel.crit)
perform.fs <- FALSE
param.fs <- list()
if (!tag %in% c('species', 'genus')){
  perform.fs <- TRUE
  param.fs.ss <- 
    list(thres.fs = as.numeric(
      parameters$model.building$feature.selection$cutoff),
         method.fs = parameters$model.building$feature.selection$type)
  param.fs.loso <- 
    list(thres.fs = 3200,
      method.fs = parameters$model.building$feature.selection$type)
}

# overwrite ml method with command line argument
if (length(args) == 2){
  ml.method <- args[2]
}

# ##############################################################################
# Get Data
fn.feat <- paste0('../data/', tag, '/feat_rel_crc.tsv')
feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE, 
                                 stringsAsFactors = FALSE, 
                                 check.names = FALSE, quote=''))

meta <- read_tsv('../data/meta/meta_crc.tsv')
stopifnot(all(meta$Sample_ID %in% colnames(feat.all)))

# ##############################################################################
# Model Building
models <- list()
for (study in parameters$ref.studies){
  # single study model
  meta.train <- meta %>%
    filter(Study == study)

  feat.train <- feat.all[,meta.train %>% pull(Sample_ID)]

  meta.train <- data.frame(meta.train)
  rownames(meta.train) <- meta.train$Sample_ID

  siamcat <- siamcat(feat=feat.train, meta=meta.train,
                     label = 'Group', case='CRC')
  siamcat <- normalize.features(siamcat, norm.method = norm.method,
                                norm.param = n.p, feature.type = 'original',
                                verbose=3)
  siamcat <- create.data.split(siamcat, num.folds = num.folds,
                               num.resample = num.resample)
  siamcat <- train.model(siamcat,
                         method = ml.method,
                         modsel.crit=modsel.crit,
                         min.nonzero.coeff = min.nonzero.coeff,
                         perform.fs = perform.fs,
                         param.fs = param.fs.ss)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)
  models[[study]] <- siamcat
  save(siamcat, file=paste0('../models/',tag,'/', study, '_', 
                            ml.method ,'_model.RData'))
  cat("Successfully trained a single study model for study", study, '\n')

  # LOSO models
  meta.train <- meta %>%
    filter(Study != study)
  
  feat.train <- feat.all[,meta.train %>% pull(Sample_ID)]
  
  meta.train <- data.frame(meta.train)
  rownames(meta.train) <- meta.train$Sample_ID
  
  siamcat <- siamcat(feat=feat.train, meta=meta.train, 
                     label = 'Group', case='CRC')
  siamcat <- normalize.features(siamcat, norm.method = norm.method, 
                                norm.param = n.p, feature.type = 'original', 
                                verbose=3)
  siamcat <- create.data.split(siamcat, num.folds = num.folds, 
                               num.resample = num.resample)
  siamcat <- train.model(siamcat,
                         method = ml.method,
                         modsel.crit=modsel.crit,
                         min.nonzero.coeff = min.nonzero.coeff,
                         perform.fs = perform.fs,
                         param.fs = param.fs.loso)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)
  models[[paste0(study, '_LOSO')]] <- siamcat
  save(siamcat, file=paste0('../models/',tag,'/', study, '_loso_', 
                            ml.method, '_model.RData'))
  cat("Successfully trained a LOSO model for study", study, '\n')
}


# ##############################################################################
# make Predictions
pred.matrix <- matrix(NA, nrow=nrow(meta), 
                      ncol=length(parameters$ref.studies)+1, 
                      dimnames = list(meta$Sample_ID, 
                                      c(parameters$ref.studies, 'LOSO')))

for (study in parameters$ref.studies){
  
  # load model
  siamcat <- models[[study]]
  temp <- rowMeans(pred_matrix(siamcat))
  pred.matrix[names(temp), study] <- temp
  
  # predict other studies
  for (study_ext in setdiff(parameters$ref.studies, study)){
    
    meta.test <- meta %>%
      filter(Study == study_ext)
    
    feat.test <- feat.all[,meta.test %>% pull(Sample_ID)]
    
    meta.test <- data.frame(meta.test)
    rownames(meta.test) <- meta.test$Sample_ID
    
    siamcat.test <- siamcat(feat=feat.test)
    
    siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
    
    temp <- rowMeans(pred_matrix(siamcat.test))
    pred.matrix[names(temp), study] <- temp
    
  }
  
}

# ##############################################################################
# make LOSO Predictions
for (study in parameters$ref.studies){
  
  # load model
  siamcat <- models[[paste0(study, '_LOSO')]]
  
  meta.test <- meta %>%
    filter(Study == study)
  
  feat.test <- feat.all[,meta.test %>% pull(Sample_ID)]
  
  meta.test <- data.frame(meta.test)
  rownames(meta.test) <- meta.test$Sample_ID
  
  siamcat.test <- siamcat(feat=feat.test)
  
  siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
  
  temp <- rowMeans(pred_matrix(siamcat.test))
  pred.matrix[names(temp), 'LOSO'] <- temp
}

# ##############################################################################
# save predictions
write.table(pred.matrix, file=paste0('../files/',tag,'/predictions_', 
                                     ml.method, '.tsv'), 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

cat('Successfully nuild models in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################
