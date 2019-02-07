# ##############################################################################
#
##  Connect functional and taxonomic profiles
#
# ##############################################################################

# Packages
library("tidyverse")
library("mlr")
library("SIAMCAT")
library("matrixStats")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
library("cowplot")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')
log.n0 <- as.numeric(parameters$model.building$log.n0)
log.n0.func <- as.numeric(parameters$model.building$log.n0.func)

# ##############################################################################
# Get Data

# meta
meta <- read_tsv('../data/meta/meta_crc.tsv')

gmm.module.info <- read_tsv('../files/GMM/GMM_info.tsv')

# species
feat.species <- read.table('../data/species/feat_rel_crc.tsv', sep='\t',
                           stringsAsFactors = FALSE, check.names = FALSE,
                           quote = '')
feat.species <- as.matrix(feat.species)

p.adj <- read.table('../files/species/p_adj.tsv', sep='\t', quote='',
                    stringsAsFactors = FALSE, check.names = FALSE)
fc.mat <- read.table('../files/species/fc.tsv', sep='\t', quote='',
                    stringsAsFactors = FALSE, check.names = FALSE)
# GMM
feat.gmm <- read.table('../data/GMM/feat_rel_crc.tsv', sep='\t',
                       stringsAsFactors = FALSE, check.names = FALSE,
                       quote = '')
feat.gmm <- as.matrix(feat.gmm)

# species clustering
load('../files/species/clustering.RData')

# ##############################################################################
# preprocess species data
stopifnot(all(rownames(p.adj) == rownames(feat.species)))

# restrict to associated species
feat.species <- feat.species[p.adj$all < 0.005,]

# preprocess GMM data as well
p.gmm <- read.table('../files/GMM/p_adj.tsv', sep='\t', quote = '',
                    stringsAsFactors = FALSE, check.names = FALSE)
feat.gmm.red <- feat.gmm[rownames(p.gmm[p.gmm$all < 0.01,]),]

# ##############################################################################
# correlation plot
cor.mat <- cor(t(log10(feat.species + log.n0)), 
               t(log10(feat.gmm.red + log.n0.func)))
cor.mat[is.na(cor.mat)] <- 0
# annotate the heatmap
annot.df <- data.frame(gmm.module.info %>% select(Module, HL1))
rownames(annot.df) <- annot.df$Module
annot.df$Module <- NULL
annot.df <- annot.df[rownames(feat.gmm.red),,drop=FALSE]
annot.df$HL1[!annot.df$HL1 %in% names(parameters$plotting$gmm.cols)] <- NA
bottom.annot <- HeatmapAnnotation(
  df = annot.df, 
  col=list('HL1' = unlist(parameters$plotting$gmm.cols))
  )
side.annot.df <- data.frame(cluster=clustering[rownames(cor.mat)],
                            fc=fc.mat[rownames(cor.mat), 'all'])
rownames(side.annot.df) <- rownames(cor.mat)
side.annot.df$cluster [side.annot.df$fc < 0] <- '-1'
side.annot.df$fc <- NULL
side.annot <- HeatmapAnnotation(
  df = side.annot.df,
  col=list('cluster' = c('1'='green', '2'='orange','3'='blue', 
                         '4'='red', '-1'='slategrey')),
  which = 'row'
)
temp <- Heatmap(cor.mat, 
                show_column_names = TRUE,
                bottom_annotation = bottom.annot,
                show_row_names = TRUE,
                clustering_method_rows = 'mcquitty',
                clustering_method_columns = 'mcquitty', 
                col=colorRamp2(seq(from=-0.4, to=0.4, length.out = 19),
                  c(rev(brewer.pal(9, 'Blues')), 'white', 
                    brewer.pal(9, 'Reds'))),
                column_title = 'GMMs',
                row_title = 'Associated gut microbial species (FDR 0.005)')
pdf('../figures/GMM/correlation_plot_FDR_005.pdf', height = 16, width = 12)
draw(side.annot + temp, heatmap_legend_side='right', 
     annotation_legend_side='right')
dev.off()

# significance heatmap
cor.mat.signif <- cor.mat
for (i in rownames(cor.mat)){
  for (j in colnames(cor.mat)){
    cor.mat.signif[i, j] <- 
      cor.test(log10(feat.species[i,] + log.n0), 
               log10(feat.gmm.red[j,] + log.n0.func))$p.value
  }
}
cor.mat.signif <- cor.mat.signif < 0.005
cor.mat.signif <- apply(cor.mat.signif, 2, as.numeric)
rownames(cor.mat.signif) <- rownames(cor.mat)
temp2 <- Heatmap(cor.mat.signif,
                 row_order = labels(row_dend(temp)[[1]]),
                 column_order = labels(column_dend(temp)),
                 cluster_rows = FALSE, cluster_columns = FALSE)
pdf('../figures/GMM/correlation_plot_FDR_005_significance.pdf', 
    height = 16, width = 12)
draw(temp2)
dev.off()

# ##############################################################################
# link species and GMMs by machine learning
# normalize as log.std
normalize.matrix <- function(matrix, log.count){
  matrix.log <- log10(matrix + log.count)
  m <- rowMeans(matrix.log)
  s <- rowSds(matrix.log)
  q <- quantile(s, 0, names = FALSE)
  matrix.norm <- (matrix.log - m) / (s + q)
  return(matrix.norm)
}

# preprocess species
feat.species.markers <- feat.species[rownames(p.adj[p.adj$all < 1e-5,]),]
feat.species.markers.ext <- feat.species[rownames(p.adj[p.adj$all < 0.005,]),]
feat.species.norm <- normalize.matrix(feat.species, log.n0)
feat.species.markers.norm <- normalize.matrix(feat.species.markers, log.n0)
feat.species.markers.ext.norm <- normalize.matrix(feat.species.markers.ext, 
                                                  log.n0)

# preprocess GMM data as well
feat.gmm.red <- feat.gmm[rownames(p.gmm[p.gmm$all < 0.01,]),]
feat.gmm.norm <- normalize.matrix(feat.gmm.red, log.n0.func)

# ##############################################################################
# function for ML
predict.gmms <- function(feat.mat, label.mat, alpha=NULL){
  
  # mlr stuff
  if (is.null(alpha)){
    lrn <- makeLearner(cl='regr.cvglmnet', predict.type = 'response')
  } else {
    lrn <- makeLearner(cl='regr.cvglmnet', predict.type = 'response',
                       alpha=alpha)
  }
  rdesc <- makeResampleDesc('CV', iters=10)
  
  df.temp <- tibble()
  # loop through the labels
  for (gmm in rownames(label.mat)){
    feat.task <- data.frame(t(feat.mat))
    feat.task$gmm <- label.mat[gmm,]
    task <- makeRegrTask(target='gmm', data=feat.task)
    
    model <- resample(learner = lrn, task=task, resampling = rdesc, 
                      show.info = FALSE, models = FALSE)
    df.temp <- bind_rows(df.temp, tibble(module=gmm, 
                                         mse=mean(model$measures.test$mse)))
    cat('Finished modeeling for:', gmm, '\n')
  }
  
  return(df.temp)
  
}

# ##############################################################################
# function for enet
test <- predict.gmms(feat.mat=feat.species.markers.norm, 
                     label.mat=feat.gmm.norm)
test2 <- predict.gmms(feat.mat=feat.species.norm,
                      label.mat=feat.gmm.norm)
test3 <- predict.gmms(feat.mat=feat.species.markers.ext.norm,
                      label.mat=feat.gmm.norm)

df.plot <- full_join(test2, test, by='module')
colnames(df.plot) <- c('module', 'All species', 'Marker set')
g.enet <- df.plot %>% 
  mutate(module=factor(module, labels(column_dend(temp)))) %>% 
  gather(key=type, value=mse, -module) %>% 
  ggplot(aes(x=module, y=mse, shape=type)) +
  geom_point() + 
  theme(panel.grid.major.y = element_line(colour = 'lightgrey')) + 
  scale_shape_manual(values=c(15, 2)) + 
  xlab('') +
  theme(axis.text.x = element_text(angle=90),
        axis.ticks.x = element_blank())
ggsave(g.enet, filename = '../figures/revision/enet_link_gmm_species.pdf',
       width = 7, height = 4)
wilcox.test(test$mse, test2$mse, paired = TRUE)
#p-value = 2.328e-10