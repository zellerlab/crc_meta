# ##############################################################################
#
##  External validation, Classification
#
# ##############################################################################

# Packages
library("tidyverse")
library("SIAMCAT")
library("pROC")
library("yaml")
library("matrixStats")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]

# ##############################################################################
parameters <- yaml.load_file('../parameters.yaml')
ml.method <- parameters$model.building$ml.method
if (length(args) == 2){
  ml.method <- args[2]
}
ref.studies <- parameters$ref.studies
ext.studies <- parameters$ext.crc.studies

# ##############################################################################
# load metadata
meta <- read_tsv('../data/meta/meta_crc_ext.tsv')

# load features
feat.rel <- read.table(paste0('../data/', tag, '/feat_rel_crc_ext.tsv'), 
                       sep='\t',  quote='',
                       stringsAsFactors = FALSE, check.names = FALSE)

# ##############################################################################
# reference studies
siamcat.list <- list()
siamcat.list[['CV']] <- list()
siamcat.list[['LOSO']] <- list()
for (study.ref in ref.studies){
  siamcat.list[['CV']][[study.ref]] <- get(load(
    paste0('../models/', tag, '/', study.ref, '_', ml.method, '_model.RData')))
  siamcat.list[['LOSO']][[study.ref]] <- get(load(
    paste0('../models/', tag, '/', study.ref, '_loso_', 
           ml.method, '_model.RData')))
}


# ##############################################################################
# loop through studies
# meta
df.plot <- tibble()
for (ext.study in ext.studies){
  cat('############################################\n', ext.study, '\n')
  meta.red <- meta %>%
    filter(Study == ext.study) %>%
    filter(Group %in% c('CTR', 'CRC'))

  # feat
  stopifnot(all(meta.red$Sample_ID %in% colnames(feat.rel)))
  feat.red <- feat.rel[,meta.red$Sample_ID]

  # get f.idx from pipeline_data_cleaning
  meta.red <- data.frame(meta.red)
  rownames(meta.red) <- meta.red$Sample_ID
  siamcat.holdout <- siamcat(feat=feat.red, meta=meta.red,
                             label='Group', case='CRC', verbose=0)
  pred.mat <- matrix(NA, nrow=ncol(feat.red), ncol=2*length(ref.studies),
                     dimnames=list(colnames(feat.red),
                                   paste0(rep(ref.studies, each=2),
                                          c('-CV', '-LOSO'))))
  for (study.ref in ref.studies){
    # CV
    siamcat <- siamcat.list[['CV']][[study.ref]]
    siamcat.holdout <- make.predictions(siamcat, siamcat.holdout, verbose=0)
    siamcat.holdout <- evaluate.predictions(siamcat.holdout, verbose=0)
    df.plot <- bind_rows(
      df.plot, tibble(type='CV', model=study.ref,
                      auroc=eval_data(siamcat.holdout)$auroc,
                      ext.study=ext.study))

    pred.mat[,paste0(study.ref, '-CV')] <-
      rowMeans(pred_matrix(siamcat.holdout))
    # LOSO
    siamcat <- siamcat.list[['LOSO']][[study.ref]]
    siamcat.holdout <- make.predictions(siamcat, siamcat.holdout, verbose=0)
    siamcat.holdout <- evaluate.predictions(siamcat.holdout, verbose=0)
    df.plot <- bind_rows(
      df.plot, tibble(type='LOSO', model=study.ref,
                      auroc=eval_data(siamcat.holdout)$auroc,
                      ext.study=ext.study))
    pred.mat[,paste0(study.ref, '-LOSO')] <-
      rowMeans(pred_matrix(siamcat.holdout))

  }
  # all LOSO
  temp <- rowMeans(pred.mat[,grep('LOSO', colnames(pred.mat))])
  df.plot <- bind_rows(
    df.plot, tibble(type='LOSO-all', model='all',
                    auroc=roc(predictor=temp, response=meta.red$Group)$auc,
                    ext.study=ext.study)
  )
}

# ##############################################################################
# plot results
temp <- df.plot %>% 
  filter(type!='LOSO') %>%
  group_by(type, ext.study) %>%
  summarize(mean=mean(auroc), sd=sd(auroc))

g <- df.plot %>%
  filter(type!='LOSO') %>%
  ggplot(aes(x=ext.study, y=auroc, fill=type)) +
    geom_bar(stat='summary', fun.y='mean', position = position_dodge(), 
             col='black') + 
    geom_point(position = position_dodge(0.9)) + 
    geom_errorbar(data=temp, inherit.aes = FALSE, 
                  aes(x=ext.study, ymin=mean-sd, ymax=mean+sd),
                  position = position_dodge(0.9),
                  width=0.2) + 
    coord_cartesian(ylim = c(0.5, 1)) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(colour='grey')) +
    scale_fill_manual(values=c('white', 'darkgrey'), name='Model Type') +
    xlab('') + ylab('AUROC')
ggsave(g, filename = paste0('../figures/', tag,
                            '/external_classification_', ml.method, '.pdf'),
       width = 6, height = 5)
