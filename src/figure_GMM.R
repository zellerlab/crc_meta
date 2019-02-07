# ##############################################################################
#
##  Plot GMM figure
#
# ##############################################################################

# Packages
library("tidyverse")
library("RColorBrewer")
library("yaml")
library("coin")

# ##############################################################################
# general 
start.time <- proc.time()[1]

set.seed(2018)
parameters <- yaml.load_file('../parameters.yaml')
study.cols <- unlist(parameters$plotting$study.cols)
gmm.cols <- unlist(parameters$plotting$gmm.cols)
group.cols <- unlist(parameters$plotting$group.cols)

# ##############################################################################
# Get Data
# meta
meta <- read_tsv('../data/meta/meta_crc.tsv') %>% 
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

# feat
fn.feat <- '../data/GMM/feat_rel_crc.tsv'
feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE, 
                                 stringsAsFactors = FALSE, 
                                 check.names = FALSE, quote=''))
feat.all <- feat.all[,meta$Sample_ID]

# pvals
fn.pvals <- '../files/GMM/p_adj.tsv'
pvals <- read.table(fn.pvals, sep='\t', header=TRUE, 
                    stringsAsFactors = FALSE, 
                    check.names = FALSE, quote='')
# FC
fn.fc <- '../files/GMM/fc.tsv'
fc <- read.table(fn.fc, sep='\t', header=TRUE, 
                    stringsAsFactors = FALSE, 
                    check.names = FALSE, quote='')

stopifnot(all(colnames(feat.all) == meta$Sample_ID))
stopifnot(all(rownames(pvals) %in% rownames(feat.all)))
stopifnot(all(rownames(fc) %in% rownames(feat.all)))

# GMM info
module.info <- read_tsv('../files/GMM/GMM_info.tsv')

# ##############################################################################
# plot heatmap

# ##############################################################################
# order stuff
pvals <- pvals[!is.na(pvals$all),]
fc <- fc[rownames(pvals),]
df.plot <- tibble(log.p.val=-log10(pvals$all)*
                    sign(fc$all),
                  module=rownames(pvals),
                  text=module.info$Name[match(rownames(pvals), 
                                              module.info$Module)],
                  type=module.info$HL1[match(rownames(pvals), 
                                              module.info$Module)]) %>% 
  arrange(log.p.val) %>% 
  filter(abs(log.p.val) > 2)
fc.mat <- fc[df.plot$module,-6]


pdf('../figures/GMM/GMM_heatmap.pdf', height = 8, width = 8)

  # layout
  layout(matrix(c(2,1,3,5,4,6), nrow=2, byrow = TRUE), 
         heights = c(0.9, 0.1), widths = c(0.3, 0.45, 0.25))

  # plot heatmap
  par(mar=c(0,0,2,.5))
  mx <- max(abs(range(fc.mat, na.rm=TRUE)))
  mx <- ifelse(round(mx, digits = 1) < mx, round(mx, digits = 1) + 0.1, 
               round(mx, digits = 1))
  mx <- 0.4
  fc.mat[fc.mat > mx] <- mx
  brs = seq(-mx, mx, by=0.05)
  num.col.steps = length(brs)-1
  n = floor(0.45*num.col.steps)
  col.hm = c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
             rep('#FFFFFF', num.col.steps-2*n),
             colorRampPalette(brewer.pal(9, 'Reds'))(n))
  image(t(as.matrix(fc.mat)),
         col=col.hm,
        breaks=brs,
        xaxt='n', yaxt='n',
        ylab='', xlab='')
  
  # plot functional category
  par(mar=c(0,14,2,.5))
  plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(0, 1), ylim=c(0, nrow(fc.mat)), type='n')
  barplot(as.matrix(rep(1, nrow(fc.mat))), 
          col=unlist(ifelse(df.plot$type %in% names(gmm.cols),
                     gmm.cols[df.plot$type],
                     'white')),
          width=1, space=0, border=NA, axes=FALSE, add=TRUE, names.arg=FALSE)
  axis(side=2, labels = df.plot$text, 
       at=seq(1, nrow(fc.mat))-0.5, tick=FALSE, las=2)
  
  # plot significance as barplot
  par(mar=c(0,1,2,2))
  plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(0, ceiling(max(abs(df.plot$log.p.val)))), 
       ylim=c(0, nrow(fc.mat)+0.25), type='n')
  grid(NULL, NA)
  barplot(t(as.matrix(abs(df.plot$log.p.val))), col='darkgrey', horiz=TRUE, 
          width=.8, space=1/4, border=NA, axes=FALSE, add=TRUE)
  axis(side=1)
  mtext('-log10(q-Value)', side=1, line=2.5)
  box()
  
  # plot color scheme for studies
  par(mar=c(1,0,2,.5))
  plot(NULL, xlab='', ylab='', xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(0, 1), ylim=c(0, 1), type='n')
  barplot(as.matrix(rep(1/ncol(fc.mat),ncol(fc.mat))), col=study.cols,
          horiz=TRUE, ylab='', xlab='', border=0, xaxt='n', add=TRUE)
  text(sub(names(study.cols), replacement = '', pattern = '-CRC'), 
       x=seq(0.1,0.9,length.out = ncol(fc.mat)), y=0.7)
  
  # plot heatmap color scheme
  par(mar=c(3,2,2,2))
  barplot(as.matrix(rep(1,length(col.hm))), col = col.hm,
          space=0,
          horiz=TRUE, border=0, ylab='', axes=FALSE)
  key.ticks <- format(seq(-mx, mx, length.out=7), digits=1)
  key.label <- 'Feature fold change over controls'
  axis(side=1, at=seq(0, length(col.hm), length.out=7), labels=key.ticks)

dev.off()


# ##############################################################################
# plot boxplot

# get data
module.info$pval <- pvals[match(module.info$Module, rownames(pvals)), 'all']

amino.acid.modules <- module.info %>% 
  filter(HL1=='amino acid degradation') %>% 
  filter(pval < 0.01) %>% 
  pull(Module)
carbohydrate.degradation <- module.info %>% 
  filter(HL1=='carbohydrate degradation') %>% 
  filter(pval < 0.01) %>% 
  pull(Module)
organic.acid.metabolism <- module.info %>% 
  filter(HL1=='organic acid metabolism') %>% 
  filter(pval < 0.01) %>% 
  pull(Module)
glycoprotein.degradation <- module.info %>% 
  filter(HL1=='glycoprotein degradation') %>% 
  filter(pval < 0.01) %>% 
  pull(Module)


df.plot <- tibble(
  `amino acid degradation` = 
    colMeans(log10(feat.all[amino.acid.modules,] + 1e-08), na.rm=TRUE),
  `glycoprotein degradation` = 
    colMeans(log10(feat.all[glycoprotein.degradation,,drop=FALSE] + 1e-08), 
             na.rm=TRUE),
  `carbohydrate degradation` = 
    colMeans(log10(feat.all[carbohydrate.degradation,] + 1e-08), na.rm=TRUE),
  `organic acid metabolism` =
    colMeans(log10(feat.all[organic.acid.metabolism,] + 1e-08), na.rm=TRUE),
  Group=meta$Group,
  Study=meta$block)

g <- df.plot %>% 
  gather(key=group, value=value, -Group, -Study) %>% 
  mutate(Group=factor(Group, levels=names(group.cols))) %>% 
  mutate(group=factor(group, levels=names(gmm.cols))) %>% 
  ggplot(aes(x=group, fill=Group, y=value)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(col=Group), position = position_jitterdodge(), 
                size=0.5, shape=20) +
    theme_classic() + 
    xlab('') + ylab('log(normalized Abundance)') + 
    scale_fill_manual(values=alpha(group.cols, 0.75), guide=FALSE) + 
    scale_colour_manual(values=group.cols, guide=FALSE) +
    theme(axis.text = element_text(size=7), axis.title = element_text(size=7))
ggsave(g, filename = '../figures/GMM/GMM_boxplots.pdf', 
       useDingbats = FALSE, width = 4, height = 4)

# calculate p values
df.temp <- df.plot %>% 
  gather(key=key, value=value, -Group, -Study)
for (i in colnames(df.plot)[-c(5,6)]){
  t <- wilcox_test(value~as.factor(Group)|as.factor(Study), data=df.temp %>% 
                filter(key==i))
  print(i)
  print(pvalue(t))
}
