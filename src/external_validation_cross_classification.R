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

# ##############################################################################
parameters <- yaml.load_file('../parameters.yaml')
ml.method <- parameters$model.building$ml.method
args = commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  ml.method <- args[1]
}
ref.studies <- parameters$ref.studies
ext.studies <- parameters$ext.studies

# ##############################################################################
# load metadata
meta <- read_tsv('../data/meta/meta_ext.tsv')
meta.crc <- read_tsv('../data/meta/meta_crc.tsv')

# reference predictions
pred.ref <- read.table(
  paste0('../files/species/predictions_', ml.method, '.tsv'), 
  sep='\t', quote='', stringsAsFactors = FALSE, check.names = FALSE)

# load features
feat <- read.table(file = '../data/species/feat_rel_ext.tsv', sep='\t', 
                   quote='', stringsAsFactors = FALSE, check.names = FALSE)

# ##############################################################################
# reference studies
siamcat.list <- list()
siamcat.list[['CV']] <- list()
siamcat.list[['LOSO']] <- list()
for (study.ref in ref.studies){
  siamcat.list[['CV']][[study.ref]] <- get(load(
    paste0('../models/species/', study.ref, '_', ml.method, '_model.RData')))
  siamcat.list[['LOSO']][[study.ref]] <- get(load(
    paste0('../models/species/', study.ref, '_loso_', 
           ml.method, '_model.RData')))
}

# ##############################################################################
# predict external studies
# meta

df.pred <- tibble()
for (ext.study in ext.studies){
  
  cat('############################################\n', ext.study, '\n')
  meta.red <- meta %>%
    filter(Study == ext.study)
    
  stopifnot(all(meta.red$Sample_ID %in% colnames(feat)))
  feat.red <- feat[,meta.red$Sample_ID]

  
  if (any(is.na(colSums(feat.red)))){
    feat.red <- feat.red[,!is.na(colSums(feat.red))]
    meta.red <- meta.red %>% filter(Sample_ID %in% colnames(feat.red))
  }
  
  # get f.idx from pipeline_data_cleaning
  meta.red <- data.frame(meta.red)
  rownames(meta.red) <- meta.red$Sample_ID
  siamcat.holdout <- siamcat(feat=feat.red, meta=meta.red,
                             label='Group', case='CTR', verbose=0)
  pred.mat <- tibble(Sample_ID=meta.red$Sample_ID, Group=meta.red$Group)
  for (study.ref in ref.studies){
    # CV
    siamcat <- siamcat.list[['CV']][[study.ref]]
    siamcat.holdout <- make.predictions(siamcat, siamcat.holdout, verbose=0)
    pred.mat[[paste0(study.ref, '-CV')]] <- 
      rowMeans(pred_matrix(siamcat.holdout))
  
    # LOSO
    siamcat <- siamcat.list[['LOSO']][[study.ref]]
    siamcat.holdout <- make.predictions(siamcat, siamcat.holdout, verbose=0)
    pred.mat[[paste0(study.ref, '-LOSO')]] <- 
      rowMeans(pred_matrix(siamcat.holdout))
    cat('finished predictions for study: ', study.ref, '...\n', sep='')
  }
  
  # all LOSO
  pred.mat[['LOSO']] <- rowMeans(pred.mat[,grep('LOSO', colnames(pred.mat))])
  df.pred <- bind_rows(df.pred, pred.mat)
}

# ##############################################################################
# compute fpr
df.pred.red <- df.pred %>% 
  filter(!Group %in% c('CTR', 'IGT'))

df.plot <- tibble()
for (study.ref in ref.studies){
  roc.ref <- roc(cases=pred.ref[meta.crc %>% 
                                  filter(Study==study.ref) %>% 
                                  filter(Group=='CRC') %>% 
                                  pull(Sample_ID), study.ref],
                 controls=pred.ref[meta.crc %>% 
                                     filter(Study==study.ref) %>% 
                                     filter(Group=='CTR') %>% 
                                     pull(Sample_ID), study.ref])
  threshold <- roc.ref$thresholds[which(roc.ref$specificities > .9)[1]]
  temp <- tibble(Group=df.pred.red$Group, 
                 value=df.pred.red[[paste0(study.ref, '-CV')]]) %>% 
    group_by(Group) %>% 
    summarize(fpr=sum(value>threshold)/n()) %>% 
    mutate(ref=study.ref, type='CV')
  df.plot <- bind_rows(df.plot, temp)
}

# loso fpr
roc.ref <- roc(cases=pred.ref[meta.crc %>% 
                                filter(Group=='CRC') %>% 
                                pull(Sample_ID), 'LOSO'],
               controls=pred.ref[meta.crc %>% 
                                   filter(Group=='CTR') %>% 
                                   pull(Sample_ID), 'LOSO'])
threshold <- roc.ref$thresholds[which(roc.ref$specificities > .9)[1]]
temp <- tibble(Group=df.pred.red$Group, 
               value=df.pred.red[['LOSO']]) %>% 
  group_by(Group) %>% 
  summarize(fpr=sum(value>threshold)/n()) %>% 
  mutate(ref='all', type='LOSO')
temp <- bind_rows(
  temp, tibble(Group='CTR', 
               fpr=sum(pred.ref[meta.crc %>% 
                                  filter(Group=='CTR') %>% 
                                  pull(Sample_ID), 'LOSO'] > threshold)/
                 meta.crc %>% filter(Group=='CTR') %>% nrow, 
                         ref='all', type='LOSO'))
df.plot <- bind_rows(df.plot, temp)

# ##############################################################################
# plot
g <- df.plot %>% 
  group_by(Group, type) %>% 
  summarise(m=mean(fpr), sd=sd(fpr)) %>% 
  ungroup() %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'T2D', 'PD', 'UC', 'CD'))) %>% 
  ggplot(aes(x=Group, y=m, fill=type)) +
    geom_bar(stat='identity', position = position_dodge(), colour='black') +
    geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.25, 
                  position = position_dodge(width = 0.9)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) +  
    scale_fill_manual(values=c('white', 'grey')) + 
    ylab('False positive rate (FPR)')
    
ggsave(g, filename = paste0('../figures/species/figure_cross_classification_', 
                            ml.method, '.pdf'), width = 6, height = 4)
