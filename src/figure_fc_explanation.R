# ##############################################################################
#
##  Figure traditional FC versus new FC
#
# ##############################################################################

# Packages
library("tidyverse")
library("ggExtra")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')

# ##############################################################################
# get data
meta <- read_tsv('../data/meta/meta_crc.tsv')
feat <- read.table('../data/species/feat_rel_crc.tsv', sep='\t', quote='',
                   stringsAsFactors = FALSE, check.names = FALSE)
feat.log <- as.matrix(log10(feat + as.numeric(parameters$associations$log.n0)))

fc.new <- read.table('../files/species/fc.tsv', sep='\t', quote='',
                     stringsAsFactors = FALSE, check.names = FALSE)
p.adj <- read.table('../files/species/p_adj.tsv', sep='\t', quote='',
                    stringsAsFactors = FALSE, check.names = FALSE)
auroc <- read.table('../files/species/aucs.tsv', sep='\t', quote='',
                    stringsAsFactors = FALSE, check.names = FALSE)

# ##############################################################################
# calculate traditional FC
fc.trad <- c()
for (i in rownames(feat.log)){
  temp.fc <- c()
  for (s in parameters$ref.studies){
    temp.ctr <- median(feat.log[i, meta %>% 
                                  filter(Group=='CTR', Study==s) %>% 
                                  pull(Sample_ID)])
    temp.crc <- median(feat.log[i, meta %>% 
                                  filter(Group=='CRC', Study==s) %>% 
                                  pull(Sample_ID)])
    temp.fc <- c(temp.fc, temp.crc-temp.ctr)
  }
  fc.trad <- c(fc.trad, mean(temp.fc))
}

# calculate prevalence shift

# calculate prev shift
prev.shift <- c()
for (f in rownames(feat)){
  prev.crc <- sum(feat[f, meta %>% 
                         filter(Group=='CRC') %>% 
                         pull(Sample_ID)] == 0)/
    meta %>% filter(Group=='CRC') %>% nrow
  prev.ctr <- sum(feat[f, meta %>% 
                         filter(Group=='CTR') %>% 
                         pull(Sample_ID)] == 0)/
    meta %>% filter(Group=='CTR') %>% nrow
  prev.shift <- c(prev.shift, prev.ctr - prev.crc)
}

df.plot <- tibble(old=fc.trad, 
                  new=fc.new$all, 
                  p.adj=p.adj$all, 
                  signif=p.adj < 
                    as.numeric(parameters$associations$alpha.meta),
                  AUROC=auroc$all,
                  prev=prev.shift)

# ##############################################################################
# plot correlation with traditional FC
g <- df.plot %>% 
  filter(!signif) %>% 
  ggplot(aes(x=old, y=new)) + 
    geom_point(pch=16, size=2, alpha=.4, col='darkgrey') +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) + 
    xlab('Median fold change (FC)') + 
    ylab('Generalized fold change (gFC)') + 
    geom_point(data=df.plot %>% filter(signif), col='#d9544d', size=2.5)
g.marg <- ggMarginal(g, type = "histogram") 

# ######################################################################
# correlation with AUROC

g2 <- df.plot %>% 
  filter(!signif) %>% 
  ggplot(aes(x=old, y=AUROC)) + 
    geom_point(pch=16, size=2, alpha=.4, col='darkgrey') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab('Median fold change (FC)') + ylab('AUROC') +
    geom_point(data=df.plot %>% filter(signif), col='#d9544d', size=2.5) + 
    ylim(0.25, 0.75)
g.marg.2 <- ggMarginal(g2, type='histogram') 

g3 <- df.plot %>% 
  filter(!signif) %>% 
  ggplot(aes(x=new, y=AUROC)) + 
    geom_point(pch=16, size=2, alpha=.4, col='darkgrey') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab('Generalized fold change (gFC)') + ylab('AUROC') +
    geom_point(data=df.plot %>% filter(signif), col='#d9544d', size=2.5) + 
    ylim(0.25, 0.75)
g.marg.3 <- ggMarginal(g3, type='histogram') 

# ######################################################################
# correlation with prevalence shift

g4 <- df.plot %>% 
  filter(!signif) %>% 
  ggplot(aes(x=old, y=prev)) + 
    geom_point(pch=16, size=2, alpha=.4, col='darkgrey') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    geom_point(data=df.plot %>% filter(signif), col='#d9544d', size=2.5) + 
    xlab('Median fold change (FC)') + ylab('Prevalence shift')
g.marg.4 <- ggExtra::ggMarginal(g4, type='histogram') 

g5 <- df.plot %>% 
  filter(!signif) %>% 
  ggplot(aes(x=new, y=prev)) + 
    geom_point(pch=16, size=2, alpha=.4, col='darkgrey') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    geom_point(data=df.plot %>% filter(signif), col='#d9544d', size=2.5) + 
    xlab('Generalized fold change (gFC)') + ylab('Prevalence shift')
g.marg.5 <- ggExtra::ggMarginal(g5, type='histogram') 


# ######################################################################
# save everything
pdf(paste0('../figures/species/fold_change_explanation.pdf'), 
    width = 4, height = 3.5)
print(g.marg)
plot.new()
print(g.marg.2)
plot.new()
print(g.marg.3)
plot.new()
print(g.marg.4)
plot.new()
print(g.marg.5)
dev.off()

# #######################
# End of script
# #######################
