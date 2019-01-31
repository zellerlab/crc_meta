# ##############################################################################
#
##  Compute clustering of marker species
#
# ##############################################################################

# Packages
library("tidyverse")
library("gdata")
library("vegan")
library("ComplexHeatmap")
library('RColorBrewer')
library("yaml")

parameters <- yaml.load_file("../parameters.yaml")

# ##############################################################################
# get data

# meta
meta <- read_tsv('../data/meta/meta_crc.tsv') %>% 
  mutate(Stage=case_when(Group=='CTR'~'CTR',
                         AJCC_stage %in% c('0','I', 'II')~'early',
                         AJCC_stage %in% c('III', 'IV')~'late')) %>% 
  mutate(Localization=ifelse(Localization=='TODO', NA, Localization)) 

# feat
feat <- read.table('../data/species/feat_rel_crc.tsv', sep='\t', quote='',
                   stringsAsFactors = FALSE, check.names = FALSE)

# pvals
p.adj <- read.table('../files/species/p_adj.tsv', sep='\t', quote='',
                    stringsAsFactors = FALSE, check.names = FALSE)
marker.set <- rownames(p.adj[p.adj$all < 1e-05,])

# ##############################################################################
# calculate positivity matrix
pos.mat <- apply(feat[marker.set,], 1, FUN=function(x){
  cutoff <- quantile(x[meta %>% filter(Group=='CTR') %>% pull(Sample_ID)], 
                     probs = c(1.0-0.05))
  ifelse(x > cutoff, 1, 0)
})
rownames(pos.mat) <- meta$Sample_ID

# ##############################################################################
# cluster species
jaccard.dist = vegdist(feat[marker.set,meta %>% 
                              filter(Group=='CRC') %>% 
                              pull(Sample_ID)], 
                       method='jaccard', binary=TRUE)
jacc.dend <- hclust(jaccard.dist, method='ward.D2')
clustering <- cutree(jacc.dend, k=4)
save(clustering, file='../files/species/clustering.RData')

# ##############################################################################
# plot jaccard distance between marker species
hm.mat <- as.matrix(jaccard.dist)
hm.mat <- 1 - hm.mat
ht.jaccard <- Heatmap(hm.mat, na_col = 'white',
                      col = c('white', rev(brewer.pal(11, 'RdBu')[1:6])),
                      cluster_rows = jacc.dend, cluster_columns = jacc.dend,
                      width = .8, row_names_gp = gpar(cex=.6), 
                      column_names_gp = gpar(cex=.6), 
                      row_names_max_width = unit(3, 'inches'),
                      name='Jaccard Distance')

pdf(paste0('../figures/species/jaccard_distance.pdf'), width = 8, height = 5)
draw(ht.jaccard)
dev.off()

# ##############################################################################
# plot jaccard distance withint each cluster versus the background
hm.mat.2 <- hm.mat
diag(hm.mat.2) <- NA
upperTriangle(hm.mat.2) <- NA

g.inset <- as_tibble(hm.mat.2) %>% 
  mutate(species=rownames(hm.mat.2)) %>% 
  mutate(cluster=clustering[rownames(hm.mat.2)]) %>% 
  gather(key=species.2, value=jaccard, -species, -cluster) %>% 
  mutate(cluster.2=clustering[species.2]) %>% 
  mutate(background=
           case_when(cluster!=cluster.2~'bg',
                     TRUE~paste0('Cluster ', as.character(cluster)))) %>% 
  filter(complete.cases(.)) %>% 
  ggplot(aes(x=background, y=jaccard)) +
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(width = .2, size=3, alpha=.6, pch=16) + 
    theme_classic() + 
    xlab('') + ylab('Jaccard Similarity')
ggsave(g.inset, filename = '../figures/species/cluster_similarity.pdf',
       height = 5, width = 6)

# ##############################################################################
# plot heatmap
hm.data <- t(pos.mat)[,meta %>% filter(Group == 'CRC') %>% pull(Sample_ID)]
hm.data <- hm.data[,sample(1:ncol(hm.data))]

# annotation
h.col.pred <- 
  HeatmapAnnotation(
    Pos= anno_barplot(colSums(hm.data)/nrow(hm.data), 
                      ylim=c(0, .8),border=TRUE, axis=TRUE, 
                      gp=gpar(fill='grey20', col='grey20', lwd=1.5, 
                              xlim=c(0, 0.4))))

h.col.study <- 
  HeatmapAnnotation(
    Study=meta %>% 
      arrange(match(Sample_ID, colnames(hm.data))) %>% 
      filter(Sample_ID %in% colnames(hm.data)) %>% 
      pull(Study),
    Stage=meta %>% 
      arrange(match(Sample_ID, colnames(hm.data))) %>% 
      filter(Sample_ID %in% colnames(hm.data)) %>% 
      pull(Stage),
    Loc=meta %>% 
      arrange(match(Sample_ID, colnames(hm.data))) %>% 
      filter(Sample_ID %in% colnames(hm.data)) %>% 
      pull(Localization),
    Gender=meta %>% 
      arrange(match(Sample_ID, colnames(hm.data))) %>% 
      filter(Sample_ID %in% colnames(hm.data)) %>% 
      pull(Gender),
    col=list("Study"=unlist(parameters$plotting$study.cols),
             'Stage'=unlist(parameters$plotting$stage.cols),
             'Gender'=unlist(parameters$plotting$sex.cols),
             'Loc'=unlist(parameters$plotting$loc.cols)))

h.row.pf <- 
  rowAnnotation(
    Pos = row_anno_barplot(rowSums(hm.data)/ncol(hm.data), border=TRUE, 
                           ylim=c(0,.4),
                           gp=gpar(fill='grey20', col='grey20', lwd=1.5), 
                           axis_gp=gpar(ylim=c(0,0.4), fontsize=8),
                           axis=TRUE), width=unit(1, 'inches'))

ht.main <- Heatmap(hm.data, cluster_rows = jacc.dend,
                   cluster_columns = FALSE,
                   top_annotation = h.col.pred,
                   bottom_annotation = h.col.study,
                   column_order = order(colSums(hm.data)), 
                   show_column_names = FALSE,
                   show_row_names = FALSE,
                   col=c('0'='grey95', '1'='grey20'), width=2,
                   name='Positivity Heatmap')

pdf(paste0('../figures/species/positivity_heatmap.pdf'), width = 12, height = 7)
draw(ht.main + h.row.pf, heatmap_legend_side='bottom', 
     annotation_legend_side='bottom')
dev.off()

# ##############################################################################
# plot all barplots
source('./utils_cluster_associations.R')

meta.crc <- meta %>% 
  filter(Group == 'CRC') %>% 
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>% 
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))
  

pdf('../figures/species/cluster_associations.pdf', width = 5, height = 6)

# stage
g.stage <- cluster.associations(features=t(pos.mat), 
                                meta.data=meta.crc,
                                conf='Stage',
                                block='block',
                                clustering = clustering)
g.stage + scale_fill_manual(values=unlist(parameters$plotting$stage.cols))

# localization
g.loc <- cluster.associations(features=t(pos.mat), 
                              meta.data=meta.crc,
                              conf='Localization',
                              block='block',
                              clustering = clustering)
g.loc + scale_fill_manual(values=unlist(parameters$plotting$loc.cols))

# gender
g.sex <- cluster.associations(features=t(pos.mat),
                              meta.data=meta.crc,
                              conf='Gender',
                              block='block',
                              clustering=clustering)
g.sex + scale_fill_manual(values=unlist(parameters$plotting$sex.cols))

# bmi
g.bmi <- cluster.associations(
  features=t(pos.mat),
  meta.data=meta.crc %>% 
    mutate(bmi_factor=cut(meta.crc$BMI, breaks = c(0, 25, 100), 
                          labels = c('lean', 'overweight/obese'))) %>% 
    filter(Group == 'CRC'),
  conf='bmi_factor',
  block='block',
  clustering=clustering
)
g.bmi + scale_fill_manual(values=unlist(parameters$plotting$bmi.cols))

# age
g.age <- cluster.associations(
  features=t(pos.mat),
  meta.data=meta.crc %>% 
    mutate(age_factor=cut(meta.crc$Age, breaks=quantile(meta.crc$Age), 
                          labels = c('A', 'B', 'C', 'D'),
                          include.lowest = TRUE)) %>% 
    filter(age_factor %in% c('A', 'D')) %>% 
    filter(Group == 'CRC'),
  conf='age_factor',
  block='block',
  clustering=clustering
)
g.age + scale_fill_manual(values=unlist(parameters$plotting$age.cols))

dev.off()
