# ##############################################################################
#
## bai clusters, correlation with species, and abundance plots
#
# ##############################################################################

library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("cowplot")
library("coin")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')
data.loc <- parameters$setup$data.location
v.tag <- parameters$setup$tag

# ##############################################################################
# get data
meta <- read_tsv('../data/meta/meta_crc.tsv') %>% 
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

bai.hits <- read_tsv('../files/genes/hits/bai_hits.tsv') %>% 
  filter(!is.na(cluster)) %>% 
  mutate(cluster = as.character(cluster))

bai.mat <- read.table('../data/genes/bai_genes_crc.tsv', sep='\t', quote='',
                      stringsAsFactors = FALSE, check.names = FALSE)
bai.mat <- bai.mat[,meta$Sample_ID]
# load original features (since C.scindens is filtered out in the cleaning step)
# features
fn.feat <- paste0(data.loc, 'mOTU_profiles/species_profiles_',
                  v.tag, '_motus2.0.0.tsv')
feat.ab.all <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                          header = TRUE, check.names = FALSE, row.names = 1,
                          quote='')
feat.ab.all <- as.matrix(feat.ab.all)
feat.ab.crc <- feat.ab.all[,meta$Sample_ID]
feat.rel.crc <- prop.table(feat.ab.crc, 2)


# ##############################################################################
# plot clusters

bai.cor.mat <- as.matrix(cor(log10(t(bai.mat)+1e-09)))

# calculate abundance, correlation, completenes for each cluster
df.cluster <- tibble(cluster=character(0),
                     correlation=double(0),
                     completeness=integer(0),
                     identity=double(0),
                     abundance=double(0))

for (i in unique(bai.hits$cluster)){
  temp <- bai.hits %>% 
    filter(cluster==i)
  comp <- temp %>% pull(query) %>% unique %>% length
  comp <- ifelse(comp <= 5, 5, comp)
  id <- temp %>% pull(id) %>% max
  if (nrow(temp) == 1){
    avg.cor <- 1
    avg.ab <- log10(sum(bai.mat[temp$gene_id,]))
  } else {
    temp.cor <- bai.cor.mat[temp %>% pull(gene_id),
                            temp %>% pull(gene_id)]
    diag(temp.cor) <- NA
    avg.cor <- mean(temp.cor, na.rm=TRUE)
    avg.ab <- mean(log10(rowSums(bai.mat[temp$gene_id,]) + 1e-09))
  }
  
  
  df.cluster <- df.cluster %>% 
    add_row(cluster=i,
            correlation=avg.cor,
            completeness=comp,
            identity=id,
            abundance=avg.ab)
}


g1 <- df.cluster %>% 
  mutate(complteness = factor(completeness)) %>% 
  mutate(included=cluster%in%c(bai.hits %>% filter(included) %>% 
           pull(cluster))) %>% 
  ggplot(aes(x=abundance, fill=correlation, size=completeness, y=identity,
             col=included)) + 
    geom_point(pch=21) + 
    scale_fill_gradientn(colours=c(rev(brewer.pal(9, 'Blues')), 
                                   brewer.pal(9, 'Reds')), 
                         limits=c(-1, 1), name='Mean cluster\ncorrelation') +
    scale_size_continuous(range=c(2.2,6.5), name='Operon\ncompleteness') +
    scale_colour_manual(values=c('grey90', 'black'), guide=FALSE) + 
    xlab('Average cluster abundance') + 
    ylab('Average cluster score') + 
    geom_hline(yintercept = 75, linetype=2, colour='grey')


# ##############################################################################
# plot correlation

# summarize clusters
bai.hits.red <- bai.hits %>% 
  filter(included)

# make a matrix for each query
feat.cluster <- matrix(NA, ncol=ncol(bai.mat), 
                     nrow=length(unique(bai.hits.red$cluster)),
                     dimnames=list(as.character(unique(bai.hits.red$cluster)), 
                                   colnames(bai.mat)))
for (q in as.character(unique(bai.hits.red$cluster))){
  feat.cluster[q,] <- colSums(bai.mat[bai.hits.red %>% 
                                        filter(cluster==q) %>% 
                                        pull(gene_id),])
}

# correlation
cor.mat.species <- cor(t(log10(feat.cluster + 1e-08)), 
                       t(log10(feat.rel.crc + 1e-05)))

# prepare plots
df.plot.cor <- tibble(species=character(0), 
                      cor=double(0),
                      cluster=character(0))

for (i in as.character(unique(bai.hits.red$cluster))){
  temp <- cor.mat.species[i,]
  temp <- -head(sort(-temp), n=4)
  df.temp <- tibble(species=paste0(names(temp), '_', i),
                    cor=temp, cluster=i)
  df.plot.cor <- bind_rows(df.plot.cor, df.temp)
}

g2 <- df.plot.cor %>% 
  mutate(species=factor(species, levels=species)) %>% 
  ggplot(aes(x=species, y=cor)) + 
    geom_bar(stat='identity') + 
    facet_grid(~cluster, scales='free_x') + 
    theme(strip.background = element_blank(), strip.text = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab('') + 
    ylab('Correlation') +
    theme(axis.text.x = element_text(size=8))

# ##############################################################################
# plot abundance

df.plot.ab <- meta %>% 
  select(Group, Sample_ID, Study, block) %>% 
  mutate(bai=log10(colSums(bai.mat[bai.hits %>% 
                                 filter(included) %>% 
                                 pull(gene_id), Sample_ID])+1e-08)+3) %>% 
  # 3 is added to have functions and species on the same scale...
  mutate(scindens=log10(feat.rel.crc[
    which(str_detect(rownames(feat.rel.crc), 'scindens')), 
    Sample_ID] + 1e-05)) %>% 
  mutate(hylemonae=log10(feat.rel.crc[
    which(str_detect(rownames(feat.rel.crc), 'hylemonae')), 
    Sample_ID] + 1e-05)) %>% 
  mutate(`6259`=log10(feat.rel.crc[
    which(str_detect(rownames(feat.rel.crc), 'meta_mOTU_v2_6259')), 
    Sample_ID] + 1e-05)) %>% 
  mutate(`5569`=log10(feat.rel.crc[
    which(str_detect(rownames(feat.rel.crc), 'meta_mOTU_v2_5569')), 
    Sample_ID] + 1e-05))


df.plot.ab <- df.plot.ab %>% 
  select(-Study, -Sample_ID) %>% 
  gather(key=key, value=value, -block, -Group) %>% 
  mutate(type=ifelse(key=="bai", 'bai', 'species')) %>% 
  mutate(Group=factor(Group, levels = names(parameters$plotting$group.cols)))

g3 <- df.plot.ab %>% 
  ggplot(aes(x=key, y=value)) + 
    geom_boxplot(aes(fill=Group), outlier.shape = NA) +
    geom_jitter(aes(col=Group), 
                position = position_jitterdodge(jitter.width = 0.3)) + 
    facet_grid(~type, scales = 'free_x', space='free') + 
    xlab('') + ylab('log10(normalized abundance)') + 
    scale_colour_manual(values=unlist(parameters$plotting$group.cols), 
                        guide=FALSE) + 
    scale_fill_manual(values=alpha(unlist(parameters$plotting$group.cols), 
                                   alpha = 0.75), guide=FALSE)
  
# p-values
for (i in unique(df.plot.ab$key)){
  print(i)
  t <- wilcox_test(value~Group|block, data=
                     df.plot.ab %>% filter(key==i) %>% 
                     mutate(block=factor(block)))
  print(pvalue(t))
}

# ##############################################################################
# plug everything together
pdf('../figures/genes/bai_extended.pdf', width = 7, height = 5, 
    useDingbats = FALSE)
print(g1)
print(g2)
print(g3)
dev.off()
