# ######################################################################
#
##  Feature Weight Plot
#
# ######################################################################

# packages
library("tidyverse")
library("SIAMCAT")
library("cowplot")
library("yaml")

# ######################################################################
# parameters
parameters <- yaml.load_file('../parameters.yaml')
ml.method <- parameters$model.building$ml.method
# ######################################################################
# get data

feat <- read.table('../data/species/feat_rel_crc.tsv', sep='\t', quote='',
                   stringsAsFactors = FALSE, check.names = FALSE)

auc.mat <- read.table('../files/species/aucs.tsv', sep='\t', quote='',
                      stringsAsFactors = FALSE, check.names = FALSE)

# load models
siamcat.list <- list()
siamcat.list[['CV']] <- list()
siamcat.list[['LOSO']] <- list()
for (study.ref in parameters$ref.studies){
  siamcat.list[['CV']][[study.ref]] <- get(load(
    paste0('../models/species/', study.ref, '_', 
           ml.method,'_model.RData')))
  siamcat.list[['LOSO']][[study.ref]] <- get(load(
    paste0('../models/species/', study.ref, '_loso_', 
           ml.method,'_model.RData')))
}

# ######################################################################
# calculate relative weights for CV models

df.plot <- tibble()
for (study.ref in parameters$ref.studies){
  temp <- weight_matrix(siamcat.list[['CV']][[study.ref]])
  # normalize
  temp.norm <- apply(temp, 2, FUN=function(i){i/sum(abs(i), na.rm = TRUE)})
  
  df.temp <- tibble(species=rownames(temp.norm), 
                    weight=rowMeans(temp.norm),
                    study=study.ref,
                    type='CV',
                    selected=((rowMeans(temp.norm != 0) >= 0.5) & 
                                (rowSums(abs(temp.norm) > .1) >= 10)),
                    AUROC=auc.mat[rownames(temp.norm),'all'])
  df.plot <- bind_rows(df.plot, df.temp)
  temp <- weight_matrix(siamcat.list[['LOSO']][[study.ref]])
  # normalize
  temp.norm <- apply(temp, 2, FUN=function(i){i/sum(abs(i), na.rm = TRUE)})
  
  df.temp <- tibble(species=rownames(temp.norm), 
                    weight=rowMeans(temp.norm),
                    study=study.ref,
                    type='LOSO',
                    selected=((rowMeans(temp.norm != 0) >= 0.5) & 
                                (rowSums(abs(temp.norm) > .025) >= 10)),
                    AUROC=auc.mat[rownames(temp.norm),'all'])
  df.plot <- bind_rows(df.plot, df.temp)
}


g1 <- df.plot %>% 
  filter(type=='CV') %>% 
  sample_frac() %>% 
  ggplot(aes(x=-weight, y=AUROC, col=study)) + 
  geom_hline(yintercept = 0.5, colour='grey') +
  geom_vline(xintercept=0, colour='grey') +
    geom_point() + 
    scale_color_manual(values=unlist(parameters$plotting$study.cols), 
                       guide=FALSE) +
  geom_hline(aes(yintercept=AUROC), data=. %>% filter(selected), 
             linetype=2, colour='darkgrey') +
  geom_text(data=. %>% filter(selected), 
           aes(y=AUROC, label=species, x=0.14), hjust=0, col='black') +
  ylab('Single feature AUROC') + xlab('Mean coefficient weight')

g2 <- df.plot %>% 
  filter(type=='LOSO') %>% 
  sample_frac() %>% 
  ggplot(aes(x=-weight, y=AUROC, col=study)) + 
  geom_hline(yintercept = 0.5, colour='grey') +
  geom_vline(xintercept=0, colour='grey') +
  geom_point() + 
  scale_color_manual(values=unlist(parameters$plotting$study.cols), 
                     guide=FALSE) +
  geom_hline(aes(yintercept=AUROC), data=. %>% filter(selected), 
             linetype=2, colour='darkgrey') +
  geom_text(data=. %>% filter(selected), 
            aes(y=AUROC, label=species, x=0.05), hjust=0, col='black') +
  ylab('Single feature AUROC') + xlab('Mean coefficient weight')

# ######################################################################
# number of nonzero coefficients shared
df.plot.2 <- tibble(CV=c(df.plot %>% filter(type=='CV') %>% 
                      group_by(species) %>% summarize(n=sum(weight!=0)) %>% 
                      ungroup %>% pull(n) %>% table)) %>% 
  mutate(LOSO=c(df.plot %>% filter(type=='LOSO') %>% 
                 group_by(species) %>% summarize(n=sum(weight!=0)) %>% 
                 ungroup %>% pull(n) %>% table)) %>% 
  mutate(count=c(1:nrow(.)-1)) %>% 
  gather(key=key, value=value, -count) %>% 
  filter(count!=0)
g3 <- ggplot(df.plot.2, aes(x=as.factor(count), y=value, fill=key)) +
  geom_bar(stat='identity', position = position_dodge(), col='black') + 
  xlab('Number of models sharing coefficients') + 
  ylab('Number of non-zero coefficients') + 
  theme_classic() + 
  scale_fill_manual(values=c('white', '#515151'), guide=FALSE)

# ######################################################################
# mean difference between study models
df.plot.3 <- df.plot %>%
  group_by(type, species) %>% 
  summarise(m=median(as.numeric(dist(weight)))) %>% 
  ungroup %>% 
  spread(key=type, value=m) %>% 
  mutate(alpha=CV>0.02)

g4 <- df.plot.3 %>% 
  ggplot(aes(x=CV, y=LOSO, alpha=alpha)) + 
  geom_point() +
  xlim(0, 0.08) + 
  ylim(0, 0.08) + 
  xlab('Median difference between coefficients for CV models') + 
  ylab('Median difference between coefficients for LOSO models') +
  geom_text(aes(label=species), data=. %>% filter(alpha), 
            angle=90, hjust=0, nudge_y=0.01) + 
  scale_alpha_manual(values=c(0.2, 0.8), guide=FALSE)

# ######################################################################
# plot everything together
pdf('../figures/species/figure_weights.pdf', width = 8.3, height = 11, 
    useDingbats = FALSE)
plot_grid(g1, g2, g3, g4, labels=c('a', 'b', 'd', 'e'), align='hv')
dev.off()


# ######################################################################
# distribution of nonzero coefficients
df.plot.4 <- tibble()
for (study.ref in parameters$ref.studies){
  temp <- weight_matrix(siamcat.list[['CV']][[study.ref]])
  # normalize
  temp.norm <- apply(temp, 2, FUN=function(i){i/sum(abs(i), na.rm = TRUE)})
  counts <- apply(temp.norm, 2, FUN=function(x){
    sum(cumsum(-sort(-abs(x))) < 0.8)})
  
  df.temp <- tibble(counts=counts,
                    study=study.ref)
  df.plot.4 <- bind_rows(df.plot.4, df.temp)
}

g5 <- df.plot.4 %>% 
  ggplot(aes(x=counts)) + 
    geom_density() + 
    xlab('Number of nonzero coefficients')
ggsave(g5, filename = '../figures/species/distribution_nonzero_weights.pdf',
       width = 4, height = 3)
