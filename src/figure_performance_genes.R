# ##############################################################################
#
##  Figure showing how well the classifier can be transferred between 
##    different datasets, based on the results from Paul with the 
##    complete gene catalogue
#
# ##############################################################################

# ##############################################################################
# Packages
library("tidyverse")
library("cowplot")
library("pROC")
library("yaml")

# for cowplot
theme_set(theme_bw())

parameters <- yaml.load_file('../parameters.yaml')

# ##############################################################################
# general 
auroc.values <- read_tsv('../files/genes/ml_results/auroc_values.tsv')

#LOSO values from Paul

loso <- c('AT-CRC'= 0.795767196, 'CN-CRC' = 0.861136479, 'DE-CRC'= 0.898694444,
          'FR-CRC'= 0.820507269, 'US-CRC'= 0.750115856)

# ##############################################################################
# AUROC heatmap

col.scheme.heatmap <- parameters$plotting$peformance.cols
plot.levels <- parameters$ref.studies

g <- auroc.values %>% 
  mutate(Test=str_replace(Test, '_', '-')) %>% 
  mutate(Training=str_replace(Training, '_', '-')) %>% 
  mutate(study.test=factor(Test, levels=plot.levels)) %>% 
  mutate(study.train=factor(Training, levels=rev(plot.levels))) %>% 
  mutate(CV=study.train == study.test) %>%
  mutate(AUC=Mean) %>% 
  ggplot(aes(y=study.train, x=study.test, fill=AUC, size=CV, color=CV)) +
    geom_tile() + theme_bw() +
    # test in tiles
    geom_text(aes_string(label="format(AUC, digits=2)"), col='white', size=7)+
    # color scheme
    scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
    # axis position/remove boxes/ticks/facet background/etc.
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.x.top = element_text(angle=45, hjust=.1), 
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    xlab('Test Set') + ylab('Training Set') + 
    scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
    scale_size_manual(values=c(0, 3), guide=FALSE)

# model average
g2 <- auroc.values %>% 
  filter(Test != Training) %>% 
  group_by(Training) %>% 
  summarise(AUROC=mean(Mean)) %>% 
  mutate(Training=str_replace(Training, '_', '-')) %>% 
  mutate(study.train=factor(Training, levels=rev(plot.levels))) %>% 
  ggplot(aes(y=study.train, x=1, fill=AUROC)) + 
    geom_tile() + theme_bw() +
    geom_text(aes_string(label="format(AUROC, digits=2)"), col='white', size=7)+
    scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.5, 1), 
                         guide=FALSE) + 
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    xlab('Model Average') + ylab('')

pdf(paste0('../figures/genes/performance_heatmap.pdf'), 
    width = 12, height = 7.5, useDingbats = FALSE)
plot_grid(g, g2, rel_widths = c(5/6, 2/6), align = 'h')
dev.off()

# ##############################################################################
# LOSO bars
df.plot.loso <- auroc.values %>% 
  filter(Training != Test) %>% 
  mutate(AUC=Mean) %>% 
  select(Test, AUC) %>% 
  mutate(type='Single Study') %>% 
  mutate(Test=str_replace(Test, '_', '-'))

df.plot.loso <- bind_rows(df.plot.loso, 
                          tibble(Test=names(loso), AUC=loso, type='LOSO'))


temp <- df.plot.loso %>% 
  group_by(Test, type) %>% 
  summarise(AUROC=mean(AUC), sd=sd(AUC)) %>% 
  ungroup() %>% 
  mutate(study.test = factor(Test, levels = plot.levels))

g3 <- df.plot.loso %>% 
  mutate(study.test = factor(Test, levels = plot.levels)) %>% 
  mutate(type=factor(type, levels=c('Single Study', 'LOSO'))) %>% 
  ggplot(aes(x=study.test, y=AUC, fill=type)) + 
    geom_bar(stat='summary', position=position_dodge(), fun.y='mean', 
             size=1, colour='black') + 
    geom_point(position = position_dodge(width = 0.9)) + 
    scale_y_continuous(breaks=seq(0, 1, by=0.1)) + 
    ylab('AUROC') + xlab('') + 
    coord_cartesian(ylim=c(0.5, 1), expand=FALSE) + 
    theme(panel.grid.major.x=element_blank()) + 
    theme_classic() + 
    scale_fill_manual(values=c('white', 'darkgrey'), 
                      labels=c('Dataset average', 'LOSO'), name='') +
    geom_errorbar(data=temp, inherit.aes = FALSE, 
                  aes(x=study.test, ymin=AUROC-sd, ymax=AUROC+sd), 
                  position = position_dodge(width = 0.9), width=0.2)
    
    
ggsave(g3, filename = '../figures/genes/loso_performance.pdf',
       width = 5.5, height = 3, useDingbats=FALSE)

cat('Successfully plotted performance heatmap and LOSO barplots in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################