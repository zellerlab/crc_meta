# ##############################################################################
#
##  Figure showing how well the classifier can be transferred between 
##    different datasets
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
cat('Plot performance heatmap\n')
start.time <- proc.time()[1]

set.seed(2018)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]

ml.method <- parameters$model.building$ml.method
if (length(args) == 2){
  ml.method <- args[2]
}

fn.pred <- paste0('../files/',tag,'/predictions_', ml.method, '.tsv')
if (!file.exists(fn.pred)){
  stop('The prediciton file is not available. Exiting...\n')
}

# ##############################################################################
# Get Data
meta <- read_tsv('../data/meta/meta_crc.tsv')

studies <- meta %>% pull(Study) %>% unique

pred.matrix <- read.table(fn.pred, 
                          sep='\t', check.names = FALSE)
pred.matrix$Sample_ID <- rownames(pred.matrix)
pred.matrix <- as_tibble(pred.matrix)

df.all <- inner_join(meta, pred.matrix, by='Sample_ID')

# ##############################################################################
# Calculate AUROCs
auroc.all <- tibble()

for (study.train in studies){
  for (study.test in studies){
    predictor <- df.all %>%
      filter(Study == study.test) %>% 
      pull(study.train) 
    response <- df.all %>%
      filter(Study == study.test) %>% 
      pull(Group)                  
    temp <- roc(predictor=predictor, response = response, ci=TRUE)
    
    auroc.all <- bind_rows(auroc.all, 
                           tibble(study.train=study.train, 
                                  study.test=study.test,
                                  AUC=c(temp$auc)))
    
  }
}

# ##############################################################################
# AUROC heatmap

col.scheme.heatmap <- parameters$plotting$peformance.cols
plot.levels <- parameters$ref.studies

g <- auroc.all %>% 
  mutate(study.test=factor(study.test, levels=plot.levels)) %>% 
  mutate(study.train=factor(study.train, levels=rev(plot.levels))) %>% 
  mutate(CV=study.train == study.test) %>%
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
g2 <- auroc.all %>% 
  filter(study.test != study.train) %>% 
  group_by(study.train) %>% 
  summarise(AUROC=mean(AUC)) %>% 
  mutate(study.train=factor(study.train, levels=rev(plot.levels))) %>% 
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

pdf(paste0('../figures/', tag, '/performance_heatmap_', ml.method, '.pdf'), 
    width = 12, height = 7.5, useDingbats = FALSE)
plot_grid(g, g2, rel_widths = c(5/6, 2/6), align = 'h')
dev.off()

# ##############################################################################
# LOSO bars
df.plot.loso <- auroc.all %>% 
  filter(study.test != study.train) %>% 
  select(study.test, AUC) %>% 
  mutate(type='Single Study')

for (study in studies){
  cases <- df.all %>%
    filter(Study == study) %>% 
    filter(Group=='CRC') %>% 
    pull(LOSO) 
  controls <- df.all %>%
    filter(Study == study) %>% 
    filter(Group=='CTR') %>% 
    pull(LOSO)                  
  temp <- roc(cases=cases, controls=controls)
  
  df.plot.loso <- bind_rows(
    df.plot.loso, tibble(study.test=study, AUC=c(temp$auc), 
                         type='LOSO'))
  
}

temp <- df.plot.loso %>% 
  group_by(study.test, type) %>% 
  summarise(AUROC=mean(AUC), sd=sd(AUC)) %>% 
  ungroup() %>% 
  mutate(study.test = factor(study.test, levels = plot.levels))

g3 <- df.plot.loso %>% 
  mutate(study.test = factor(study.test, levels = plot.levels)) %>% 
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
    
    

ggsave(g3, filename = paste0('../figures/', tag, '/loso_performance_', 
                             ml.method, '.pdf'),
       width = 5.5, height = 3)

cat('Successfully plotted performance heatmap and LOSO barplots in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################