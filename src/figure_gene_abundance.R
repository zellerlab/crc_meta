# ##############################################################################
#
##  Figure 4c: Gene abundance
#
# ##############################################################################

# Packages
library("tidyverse")
library("coin")
library("matrixStats")
library("cowplot")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')

# ##############################################################################
# load metadata
meta <- read_tsv('../data/meta/meta_crc.tsv') %>% 
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

# load features
# bai
bai <- as.matrix(
  read.table('../data/genes/bai_genes_crc.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))
# fada
fada <- as.matrix(
  read.table('../data/genes/fada_genes_crc.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))
# bft
bft <- as.matrix(
  read.table('../data/genes/bft_genes_crc.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))
# clb
clb <- as.matrix(
  read.table('../data/genes/clb_genes_crc.tsv', sep='\t', 
             stringsAsFactors = FALSE, check.names = FALSE))

clb.hits <- read_tsv('../files/genes/hits/clb_hits.tsv')
bai.hits <- read_tsv('../files/genes/hits/bai_hits.tsv')
fada.hits <- read_tsv('../files/genes/hits/fada_hits.tsv')
bft.hits <- read_tsv('../files/genes/hits/bft_hits.tsv')

stopifnot(all(fada.hits$gene %in% rownames(fada)))
stopifnot(all(bft.hits$gene_id %in% rownames(bft)))

# ##############################################################################
# prepare plot
df.plot <- meta %>% 
  select(Group, Sample_ID, Study, block) %>% 
  mutate(bai=log10(colSums(bai[bai.hits %>% 
                                     filter(included) %>% 
                                     pull(gene_id), Sample_ID])+1e-08)) %>% 
  mutate(clb=log10(colSums(clb[clb.hits %>% 
                                     filter(included) %>% 
                                     pull(gene_id), Sample_ID])+1e-08)) %>% 
  mutate(fada=log10(colSums(fada[,Sample_ID]) + 1e-08)) %>% 
  mutate(bft=log10(colSums(bft[,Sample_ID]) + 1e-08))

# ##############################################################################
# plot
g1 <- df.plot %>% 
  gather(key=gene, value=value, -Group, -Sample_ID, -Study, -block) %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>%
  ggplot(aes(x=gene, y=value, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), 
                aes(col=Group)) + 
    scale_fill_manual(values=alpha(unlist(parameters$plotting$group.cols), 
                                   alpha=0.75), guide=FALSE) + 
    scale_colour_manual(values=unlist(parameters$plotting$group.cols), 
                        guide=FALSE) + 
    xlab('') + ylab('log10(normalized abundance)')

ggsave(g1, filename = '../figures/genes/gene_boxplot.pdf', 
       width = 6, height = 5)

# ##############################################################################
# test

wilcox_test(fada~Group|block, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    block=factor(block)))

wilcox_test(bft~Group|block, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    block=factor(block)))

wilcox_test(clb~Group|block, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    block=factor(block)))

wilcox_test(bai~Group|block, 
            data=df.plot %>% mutate(Group=factor(Group),
                                    block=factor(block)))
