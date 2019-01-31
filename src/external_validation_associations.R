# ##############################################################################
#
##  External Validation, Associations
#
# ##############################################################################

# Packages
library("tidyverse")
library("pROC")
library("yaml")

# ##############################################################################
parameters <- yaml.load_file('../parameters.yaml')
ref.studies <- parameters$ref.studies
ext.studies <- parameters$ext.crc.studies
log.n0 <- as.numeric(parameters$associations$log.n0)

# ##############################################################################
# load data
meta <- read_tsv('../data/meta/meta_crc_ext.tsv') %>% 
  filter(Group %in% c('CTR', 'CRC'))

# load features
feat.rel <- read.table(file = '../data/species/feat_rel_crc_ext.tsv', sep='\t', 
                       quote='', stringsAsFactors = FALSE, check.names = FALSE)

# load reference pvalues
p.adj <- read.table('../files/species/p_adj.tsv', sep='\t',
                    stringsAsFactors = FALSE, check.names = FALSE,
                    quote='')
p.val <- read.table('../files/species/p_val.tsv', sep='\t',
                    stringsAsFactors = FALSE, check.names = FALSE,
                    quote = '')
auroc.ref <- read.table('../files/species/aucs.tsv', sep='\t',
                        stringsAsFactors = FALSE, check.names = FALSE,
                        quote = '')
fc.mat <- read.table('../files/species/fc.tsv', sep='\t', quote='',
                     stringsAsFactors = FALSE, check.names = FALSE)
marker.set <- rownames(p.adj)[
  p.adj$all < as.numeric(parameters$associations$alpha.meta)]
marker.set.extended <- rownames(p.adj)[
  p.adj$all < parameters$associations$alpha.single.study]


# ##############################################################################
# single study assocations
for (ext.study in ext.studies){
  print(ext.study)
  # meta
  meta.red <- meta %>% 
    filter(Study==ext.study)
  # feat
  stopifnot(all(meta.red$Sample_ID %in% colnames(feat.rel)))
  feat.red <- feat.rel[,meta.red$Sample_ID]
  label <- meta.red$Group
  
  # compute pvals
  p.vals <- apply(feat.red, 1, FUN=function(x){
    wilcox.test(x~label, exact=FALSE)$p.value
  })
  p.vals.adj <- p.adjust(p.vals, method='fdr')
  p.adj[[ext.study]] <- p.vals[rownames(p.adj)]

  # compute auroc
  auroc <- apply(feat.red, 1, FUN=function(x, label){
    cases=x[label=='CRC']
    controls=x[label=='CTR']
    roc(cases=cases, controls=controls)$auc
  }, label=label)
  
  auroc.ref[[ext.study]] <- auroc
  
  # compute FCs
  fcs <- c()
  for (i in rownames(feat.red)){
    x <- feat.red[i, meta.red %>% filter(Group=='CRC') %>% pull(Sample_ID)]
    y <- feat.red[i, meta.red %>% filter(Group=='CTR') %>% pull(Sample_ID)]
    q.p <- quantile(log10(x+log.n0), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+log.n0), probs=seq(.1, .9, .05))
    fcs <- c(fcs, sum(q.p - q.n)/length(q.p))
  }
  fc.mat[[ext.study]] <- fcs
}

# ##############################################################################
# plot PR curves
df.plot <- tibble()
for (study in ext.studies){
  df.temp <- tibble(
    predictor=-log10(p.adj[[study]]),
    fc.study=ifelse(sign(fc.mat[[study]]) != 0, sign(fc.mat$all), 1),
    fc.all=ifelse(sign(fc.mat$all) != 0, sign(fc.mat$all), 1),
    response1=as.numeric(rownames(p.adj) %in% marker.set),
    response2=as.numeric(rownames(p.adj) %in% marker.set.extended)
  )
  # remove NAs
  df.temp <- df.temp %>% 
    filter(!is.na(df.temp$predictor)) %>% 
    # set positive cases with opposite FCs to negative cases
    mutate(response1=case_when(fc.study!=fc.all ~ 0,
                               fc.study==fc.all ~ response1),
           response2=case_when(fc.study!=fc.all ~ 0,
                               fc.study==fc.all ~ response2))
  
  # marker set
  temp <- SIAMCAT:::evaluate.classifier(predictions = df.temp$predictor,
                                        test.label=df.temp$response1, 
                                        label=list(info=c('test'=1, 'neg'=0)))
  temp.2 <- SIAMCAT:::evaluate.get.pr(temp)
  threshold.0.1 <- abs(temp$thresholds + log10(0.1))
  threshold.0.01 <- abs(temp$thresholds + log10(0.01))
  
  df.plot <- bind_rows(
    df.plot, tibble(
      recall=temp.2$recall, precision=temp.2$precision, 
      threshold=temp$thresholds,
      study=study, set='Meta-FDR 1e-05',
      set.0.1=threshold.0.1 == min(threshold.0.1),
      set.0.01=threshold.0.01 == min(threshold.0.01)
    ))
  
  # extended marker set
  temp <- SIAMCAT:::evaluate.classifier(predictions = df.temp$predictor,
                                        test.label=df.temp$response2, 
                                        label=list(info=c('test'=1, 'neg'=0)))
  temp.2 <- SIAMCAT:::evaluate.get.pr(temp)
  threshold.0.1 <- abs(temp$thresholds + log10(0.1))
  threshold.0.01 <- abs(temp$thresholds + log10(0.01))
  
  df.plot <- bind_rows(
    df.plot, tibble(
      recall=temp.2$recall, precision=temp.2$precision, 
      threshold=temp$thresholds,
      study=study, set='Meta-FDR 0.05',
      set.0.1=threshold.0.1 == min(threshold.0.1),
      set.0.01=threshold.0.01 == min(threshold.0.01)
    ))
}

g1 <- df.plot %>% 
  ggplot(aes(x=recall, y=precision, col=study)) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        panel.grid.minor = element_blank()) + 
  geom_line() +
  facet_grid(.~set) +
  # scale_color_manual(values=study.cols, guide=FALSE) + 
  xlab('Recall') + ylab('Precision')

ggsave(g1, filename = '../figures/species/pr_curves_external.pdf', 
       width = 10, height = 5, useDingbats=FALSE)

# ##############################################################################
# plot ranks for marker set
temp <- apply(p.adj, 2, rank, ties.method='min')
df.plot <- as_tibble(temp) %>% 
  mutate(species=rownames(temp)) %>% 
  filter(species %in% marker.set) %>% 
  arrange(all) %>%
  mutate(species=factor(species, levels = .$species)) %>% 
  gather(key=study, value=rank, -species)

g <- df.plot %>% 
  mutate(study=factor(study, levels = c('all', 
                                        ref.studies, 
                                        ext.studies))) %>% 
  ggplot(aes(x=species, y=1, fill=rank)) + 
    geom_tile() +
    facet_grid(study~.) + 
    theme_bw() + 
    theme(axis.text.y=element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          panel.grid = element_blank(), 
          panel.border = element_rect(fill=NA, colour=NA),
          axis.ticks.x=element_blank(),
          strip.background = element_blank()) + 
    scale_fill_gradientn(colours=paste0('grey', seq(20, 90, 5)),
                         limits=c(1, 100), na.value = 'red') + 
    xlab('') + ylab('') + 
    geom_text(aes(label=rank), colour='white')

ggsave(g, filename = '../figures/species/rank_heatmap_external.pdf',
       width = 12, height = 10, useDingbats=FALSE)
