# ##############################################################################
#
##  External validation, Gene profiles
#
# ##############################################################################

# Packages
library("tidyverse")
library("ggpubr")
library("coin")
library("matrixStats")
library("cowplot")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')

# ##############################################################################
# load metadata
meta <- read_tsv('../data/meta/meta_crc_ext.tsv')

# load features
# bai
bai.ext <- as.matrix(
  read.table('../data/genes/bai_genes_ext.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))
# fada
fada.ext <- as.matrix(
  read.table('../data/genes/fada_genes_ext.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))
# bft
bft.ext <- as.matrix(
  read.table('../data/genes/bft_genes_ext.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))
# clb
clb.ext <- as.matrix(
  read.table('../data/genes/clb_genes_ext.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))

clb.hits <- read_tsv('../files/genes/hits/clb_hits.tsv')
bai.hits <- read_tsv('../files/genes/hits/bai_hits.tsv')

# ##############################################################################
# prepare plot
df.plot <- meta %>% 
  filter(Group != 'ADA') %>% 
  select(Group, Sample_ID, Study) %>% 
  mutate(bai=log10(colSums(bai.ext[bai.hits %>% 
                                     filter(included) %>% 
                                     pull(gene_id), Sample_ID])+1e-08)) %>% 
  mutate(clb=log10(colSums(clb.ext[clb.hits %>% 
                                     filter(included) %>% 
                                     pull(gene_id), Sample_ID])+1e-08)) %>% 
  mutate(fada=log10(colSums(fada.ext[,Sample_ID]) + 1e-08)) %>% 
  mutate(bft=log10(colSums(bft.ext[,Sample_ID]) + 1e-08))

# ##############################################################################
# plot
g1 <- df.plot %>% 
  gather(key=gene, value=value, -Group, -Sample_ID, -Study) %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>% 
  ggplot(aes(x=gene, y=value, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
                aes(col=Group, shape=Study)) + 
    scale_fill_manual(values=alpha(unlist(parameters$plotting$group.cols), 
                                   alpha=0.75), guide=FALSE) + 
    scale_colour_manual(values=unlist(parameters$plotting$group.cols), 
                        guide=FALSE) + 
    xlab('') + ylab('log10(abundance)')

ggsave(g1, filename = '../figures/genes/external_validation_functional.pdf', 
       width = 6, height = 5)

# ##############################################################################
# test
# bai
wilcox_test(fada~Group|Study, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    Study=factor(Study)),
            alternative='greater')

wilcox_test(bft~Group|Study, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    Study=factor(Study)),
            alternative='greater')

wilcox_test(clb~Group|Study, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    Study=factor(Study)),
            alternative='greater')

wilcox_test(bai~Group|Study, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    Study=factor(Study)),
            alternative='greater')
