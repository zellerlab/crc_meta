# ##############################################################################
#
## baiF expression/qPCR quantification plots
##  and bai heatmap/pvalues
#
# ##############################################################################

library("tidyverse")
library("RColorBrewer")
library("ggpubr")
library("cowplot")
library("pROC")
library("coin")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')

# ##############################################################################
# get data
meta <- read_tsv('../data/meta/meta_crc.tsv') %>% 
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

bai.hits <- read_tsv('../files/genes/hits/bai_hits.tsv')

bai.mat <- read.table('../data/genes/bai_genes_crc.tsv', sep='\t', quote='',
                      stringsAsFactors = FALSE, check.names = FALSE)

qPCR <- read_tsv('../files/qPCR/qPCR_values.tsv')

# ##############################################################################
# process primer CT values


df.plot <- qPCR %>% 
  group_by(GeneCoreID, Bacteria, type) %>% 
  summarise(mean.ct=mean(CT.value)) %>% 
  mutate(delta.ct=0) %>% 
  ungroup()

for (p in unique(df.plot$GeneCoreID)){
  df.plot <- df.plot %>% 
    mutate(delta.ct=ifelse(GeneCoreID != p, 
                           delta.ct, 
                           df.plot %>% 
                             filter(GeneCoreID == p, 
                                    Bacteria == '16S gene') %>% 
                             pull(mean.ct) - mean.ct))
}

# ##############################################################################
# correlate with metaG
df.plot.2 <- meta %>% 
  select(Sample_ID, External_ID, Group, Study) %>% 
  mutate(bai=log10(colSums(bai.mat[bai.hits %>% 
                                     filter(included) %>% 
                                     filter(query=='baiF') %>%  
                                     pull(gene_id),
                                   Sample_ID])+1e-08)) %>% 
  filter(Study=='DE-CRC')

df.plot.2 <- inner_join(df.plot.2, 
                        df.plot %>% filter(type=='RNA', 
                                           Bacteria == 'baiFsciF2') %>% 
                          mutate(External_ID=as.character(GeneCoreID)) %>% 
                          select(External_ID, delta.ct), by="External_ID")
df.plot.2 <- inner_join(df.plot.2,
                        df.plot %>% filter(type=='DNA', 
                                           Bacteria == 'baiFsciF2') %>% 
                          mutate(External_ID=as.character(GeneCoreID)) %>% 
                          select(External_ID, delta.ct), by="External_ID")


a <- df.plot.2 %>% 
  ggplot(aes(x=bai, y=delta.ct.y, col=Group)) + 
  geom_point() + 
  xlab('metaG quantification') + 
  ylab('gDNA qPCR quantification') +
  theme(text = element_text(size=7),
        axis.text = element_text(size=6)) +
  scale_colour_manual(values=unlist(parameters$plotting$group.cols),
                      guide=FALSE)
b <- df.plot.2 %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>% 
  ggplot(aes(x=Group, y=delta.ct.y, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col=Group), width = 0.1) +
  stat_compare_means(method.args=list(alternative='greater')) + 
  ylab('') + xlab('') +
  theme(text = element_text(size=7),
        axis.text = element_text(size=6)) +
  scale_colour_manual(values=unlist(parameters$plotting$group.cols),
                      guide=FALSE)+
  scale_fill_manual(values=alpha(unlist(parameters$plotting$group.cols), 
                                 alpha=0.75), guide=FALSE)


# ##############################################################################
# correlate with expression
c <- df.plot.2 %>% 
  ggplot(aes(x=delta.ct.y, y=delta.ct.x, col=Group)) + 
  geom_point() + 
  xlab('qPCR genomic DNA') + 
  ylab('qPCR metaT') +
  theme(text = element_text(size=7),
        axis.text = element_text(size=6)) +
  scale_colour_manual(values=unlist(parameters$plotting$group.cols),
                      guide=FALSE)

d<- df.plot.2 %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>% 
  ggplot(aes(x=Group, y=delta.ct.x, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col=Group), width = 0.1) +
  stat_compare_means(method.args=list(alternative='greater')) + 
  ylab('') + xlab('') +
  theme(text = element_text(size=7),
        axis.text = element_text(size=6)) +
  scale_colour_manual(values=unlist(parameters$plotting$group.cols),
                      guide=FALSE)+
  scale_fill_manual(values=alpha(unlist(parameters$plotting$group.cols), 
                                 alpha=0.75), guide=FALSE)

# ##############################################################################
# plot everything together

pdf('../figures/genes/qPCR_baiF.pdf',width = 10, height = 2.5, 
    useDingbats = FALSE)
plot_grid(a, b, c, d, nrow = 1, rel_widths = c(1, 0.5, 1, 0.5))
dev.off()

# calculate correlation
print(cor(df.plot.2$bai, df.plot.2$delta.ct.y))
print(cor(df.plot.2$delta.ct.y, df.plot.2$delta.ct.x))

# ##############################################################################
# plot roc plots for baiF
roc.bai.qPCR <- roc(cases=df.plot.2 %>% 
                      filter(Group=='CRC') %>% 
                      pull(delta.ct.y),
                    controls=df.plot.2 %>% 
                      filter(Group=='CTR') %>% 
                      pull(delta.ct.y), ci=TRUE)
ci.ob <- ci(roc.bai.qPCR, of='se', boot.n=500)

# roc plot for metaG bai
pdf('../figures/genes/baiF_genomic_qPCR_roc.pdf', width = 5, height = 5)
plot(roc.bai.qPCR)
plot(ci.ob, type = 'shape', border='white')
dev.off()


# ##############################################################################
# bai heatmap
bai.hits.red <- bai.hits %>% 
  filter(included) %>% 
  filter(query != 'baiK') # remove the single baiK gene

# make a matrix for each query
feat.query <- matrix(NA, ncol=ncol(bai.mat), 
                     nrow=length(unique(bai.hits.red$query)),
                     dimnames=list(unique(bai.hits.red$query), 
                                   colnames(bai.mat)))
for (q in unique(bai.hits.red$query)){
  feat.query[q,] <- colSums(bai.mat[bai.hits.red %>% 
                                      filter(query==q) %>% 
                                      pull(gene_id),])
}

# calculate FC each study and blocked p.val
df.plot.heatmap <- tibble(query=character(0), study=character(0), fc=double(0))
df.plot.pval <- tibble(query=character(0), pval=double(0))
for (q in unique(bai.hits.red$query)){
  for (s in parameters$ref.studies){
    
    x <- log10(feat.query[q, meta %>% filter(Study==s, Group=='CRC') %>% 
                            pull(Sample_ID)] + 1e-09)
    y <- log10(feat.query[q, meta %>% filter(Study==s, Group=='CTR') %>% 
                            pull(Sample_ID)] + 1e-09)
    q.p <- quantile(x, probs = seq(.1, .9, length.out = 9))
    q.n <- quantile(y, probs = seq(.1, .9, length.out = 9))
    fc <- sum(q.p - q.n)/9
    df.plot.heatmap <- df.plot.heatmap %>% 
      add_row(query=q, study=s, fc=fc)
  }
  temp <- data.frame(x=feat.query[q,meta$Sample_ID],
                     group=meta$Group, 
                     block=meta$block)
  t <- wilcox_test(x~group|block, data=temp)
  df.plot.pval <- df.plot.pval %>% 
    add_row(query=q, pval=pvalue(t))
}


# plot as heatmap
bai.levels <- c('baiB', 'baiCD', 'baiE', 'baiA', 'baiF', 
                'baiG', 'baiH', 'baiI', 'baiJ', 'baiL')
mx <- round(max(df.plot.heatmap$fc), digits = 1)
brs = seq(-mx, mx, by=0.05)
num.col.steps = length(brs)-1
n = floor(0.45*num.col.steps)
col.hm = c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
           rep('#FFFFFF', num.col.steps-2*n),
           colorRampPalette(brewer.pal(9, 'Reds'))(n))

g <- df.plot.heatmap %>% 
  mutate(study=factor(study, levels = rev(parameters$ref.studies))) %>% 
  mutate(query=factor(query, levels = bai.levels)) %>% 
  ggplot(aes(x=query, y=study, fill=fc)) + 
    geom_tile() + 
    theme_minimal() + 
    theme(panel.grid = element_blank()) + 
    scale_fill_gradientn(colours = c(col.hm), limits=c(-0.9, 0.9))

# p value barplot
g2 <- df.plot.pval %>% 
  mutate(log.p.val=-log10(pval)) %>% 
  mutate(query=factor(query, levels = bai.levels)) %>% 
  ggplot(aes(x=query, y=log.p.val)) +
    geom_bar(stat='identity') +
    ylim(0, 10)

pdf('../figures/genes/bai_fc_heatmap.pdf', width = 6, height = 6)
plot_grid(g2, g, ncol=1)
dev.off()
