# ##############################################################################
#
##  Marker Heatmap
#
# ##############################################################################

# Packages
library("tidyverse")
library("RColorBrewer")
library("cowplot")
library("yaml")

theme_set(theme_gray())

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]


parameters <- yaml.load_file('../parameters.yaml')
ref.studies <- parameters$ref.studies
start.time <- proc.time()[1]

# ##############################################################################
# get data

p.vals <- read.table(paste0('../files/', tag, '/p_adj.tsv'), 
                     sep='\t', quote='',
                     stringsAsFactors = FALSE, check.names = FALSE)

fc.mat <- read.table(paste0('../files/', tag, '/fc.tsv'), 
                     sep='\t', quote='',
                     stringsAsFactors = FALSE, check.names = FALSE)

auc.mat <- read.table(paste0('../files/', tag, '/aucs.tsv'), 
                      sep='\t', quote='',
                      stringsAsFactors = FALSE, check.names = FALSE)

# ##############################################################################
# select and order
species.heatmap <- 
  rownames(p.vals)[which(p.vals$all < 
                           parameters$associations$alpha.single.study)]
fc.sign <- sign(fc.mat)
fc.sign[fc.sign == 0] <- 1

p.val.signed <- -log10(p.vals[species.heatmap,"all", drop=FALSE])*
  fc.sign[species.heatmap, 'all']
top.markers <- rownames(p.val.signed[is.infinite(p.val.signed$all),,drop=FALSE])
p.val.signed[top.markers, 'all'] <- 100 + auc.mat[top.markers, 'all']

species.heatmap.orderd <- rownames(p.val.signed[order(p.val.signed$all),,
                                                drop=FALSE])

# take only those
fc.mat.plot <- fc.mat[species.heatmap.orderd,]
p.vals.plot <- p.vals[species.heatmap.orderd,]

# ##############################################################################
# prepare plotting

# colorscheme for fc heatmap 
mx <- max(abs(range(fc.mat.plot, na.rm=TRUE)))
mx <- ifelse(round(mx, digits = 1) < mx, 
             round(mx, digits = 1) + 0.1, 
             round(mx, digits = 1))
brs = seq(-mx, mx, by=0.05)
num.col.steps = length(brs)-1
n = floor(0.45*num.col.steps)
col.hm = c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
           rep('#FFFFFF', num.col.steps-2*n),
           colorRampPalette(brewer.pal(9, 'Reds'))(n))
# color scheme for pval heatmap
alpha.breaks=c(1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01)
p.vals.bin <- data.frame(apply(p.vals.plot, 2, FUN=.bincode, 
                               breaks = c(0, alpha.breaks, 1), 
                               include.lowest = TRUE),
                         check.names = FALSE)
p.val.greys <-  c(paste0('grey', 
                         round(seq(from=10, to=80, 
                                   length.out = length(alpha.breaks)))), 
                  'white')
names(p.val.greys) <- as.character(1:7)

# function to plot both into a grid
plot.single.study.heatmap <- function(x){
  df.plot <- tibble(species=factor(rownames(p.vals.plot), 
                                   levels=rev(rownames(p.vals.plot))),
                    p.vals=as.factor(p.vals.bin[[x]]),
                    fc=fc.mat.plot[[x]])
  
  g1 <- df.plot %>% 
    ggplot(aes(x=species, y=1, fill=fc)) + 
      geom_tile() + theme_minimal() + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(), 
            panel.grid = element_blank(),
            panel.background = element_rect(fill=NULL, colour='black'),
            plot.margin = unit(c(0, 0, 0, 0), 'cm')) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), guide=FALSE)
  
  g2 <- df.plot %>% 
    ggplot(aes(x=species, y=1, fill=p.vals)) +
      geom_tile() + theme_minimal() + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(), 
            panel.grid = element_blank(),
            panel.background = element_rect(fill=NULL, colour='black'),
            plot.margin = unit(c(0, 0, 0, 0), 'cm')) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_manual(values=p.val.greys, na.value='white', guide=FALSE)
  g.return <- plot_grid(g2, g1, ncol = 1, rel_heights = c(0.25, 0.75))
}

# ##############################################################################
# plot

# p.value histogram
g1 <- tibble(species=factor(rownames(p.vals.plot), 
                            levels=rev(rownames(p.vals.plot))),
             p.vals=-log10(p.vals.plot$all),
             colour=p.vals > 5) %>% 
  ggplot(aes(x=species, y=p.vals, fill=colour)) + 
    geom_bar(stat='identity') + 
    theme_classic() + 
    xlab('Gut microbial species') + 
    ylab('-log10(q-value)') + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill=NULL, colour = 'black')) + 
    scale_y_continuous(limits=c(0, 15), expand = c(0, 0)) + 
    scale_x_discrete(position='top') + 
    scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide=FALSE)

g.lst <- lapply(ref.studies, plot.single.study.heatmap)

pdf(paste0('../figures/', tag, '/fold_change_heatmap.pdf'), 
    width = 10, height = 6, useDingbats = FALSE)
plot_grid(g1, g.lst[[1]], g.lst[[2]], g.lst[[3]], g.lst[[4]],g.lst[[5]], 
          ncol=1, align = 'v', rel_heights = c(0.3,rep(0.12, 5)))
dev.off()

cat('Successfully plotted marker heatmap in',
    proc.time()[1]-start.time, 'second...\n')


# ##############################################################################
# Forest plot
# ##############################################################################

# load aucs
load(paste0('../files/', tag, '/aucs_all.RData'))

# ##############################################################################
# select and order
marker.set <- rownames(p.val.signed)[
  abs(p.val.signed$all) > 
    -log10(as.numeric(parameters$associations$alpha.meta))]
p.val.signed.red <- p.val.signed[marker.set, ,drop=FALSE]
marker.set.orderd <- rev(rownames(p.val.signed.red[order(p.val.signed.red$all),,
                                               drop=FALSE]))

# ##############################################################################
# extract those from the auc list
df.plot <- tibble()
for (i in marker.set.orderd){
  for (s in ref.studies){
    temp <- aucs.all[[i]][[s]]
    df.plot <- bind_rows(df.plot, tibble(
      species=i, study=s,
      low=temp[1], auc=temp[2], high=temp[3]
    ))
  }
}

df.plot <- df.plot %>% 
  mutate(species=factor(species, levels = marker.set.orderd)) %>% 
  mutate(study=factor(study, levels = ref.studies))

# ##############################################################################
# plot everything

g <- df.plot %>% 
  ggplot(aes(x=study, y=auc)) + 
    geom_linerange(aes(ymin=low, ymax=high), colour='lightgrey') + 
    geom_point(pch=23, aes(fill=study)) + 
    facet_grid(~species, scales = 'free_x', space='free') + 
    theme_minimal() + 
    scale_y_continuous(limits=c(0, 1)) + 
    theme(panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text = element_text(angle=90, hjust=0)) + 
    scale_fill_manual(values=unlist(parameters$plotting$study.cols), 
                      guide=FALSE) + 
    ylab('AUROC') + xlab('Gut microbial species')

ggsave(g, filename = paste0('../figures/', tag, '/forest_plot.pdf'), 
       width = 10, height = 10, useDingbats=FALSE)

cat('Successfully plotted forest plot in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################