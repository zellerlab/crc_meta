# ##############################################################################
#
##  TPR and Prediction bias for LOSO predictions
#
# ##############################################################################

# Packages
library("tidyverse")
library("pROC")
library("coin")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')

ml.method <- parameters$model.building$ml.method

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]
if (length(args)==2){
  ml.method <- args[2]
}


start.time <- proc.time()[1]

# ##############################################################################
# Get Data

# predictions
predictions <- read.table(paste0('../files/', tag, '/predictions_', 
                                 ml.method,'.tsv'))


# meta
meta <- read_tsv('../data/meta/meta_crc.tsv') %>% 
  mutate(Stage=case_when(Group == 'CTR' ~ 'CTR',
                         AJCC_stage %in% c('II', 'III', 'IV') ~ 
                           paste0('Stage ', AJCC_stage),
                         AJCC_stage == 0 ~ 'Stage 0/I',
                         AJCC_stage == 'I' ~ 'Stage 0/I')) %>% 
  # add predictions
  mutate(Predictions=predictions$LOSO)

# ##############################################################################
# preprocess metadata/predictions for plotting
df.plot <- meta %>% 
  
  # change continuous values to factors
  mutate(Age=as.character(cut(meta$Age, breaks = quantile(meta$Age), 
                              labels=c('Age < 57', 'Age 57 - 64', 
                                       'Age 64 - 71', 'Age > 71')))) %>% 
  mutate(BMI=case_when(BMI > 25 ~ 'BMI > 25',
                       BMI <= 25 ~ 'BMI < 25')) %>% 
  # restrict to cases
  filter(Group=='CRC') %>% 
  mutate(Localization=case_when(Localization == 'Rectum' ~ Localization,
                                Localization == 'RC' ~ 'Right Colon',
                                Localization == 'LC' ~ 'Left Colon')) %>% 
  mutate(Sex=case_when(Gender == 'F' ~ 'Female',
                       Gender == 'M' ~ 'Male')) %>% 
  select(Sex, Localization, Stage, BMI, Age, Predictions, Study) %>% 
  # reorder for plotting
  gather(key=confounder, value=value, -Predictions, -Study) %>% 
  filter(complete.cases(.)) %>% 
  # re-order for plotting
  mutate(confounder = factor(confounder, levels = c('Age', 'Sex', 'BMI', 
                                                    'Localization', 'Stage')),
         value=factor(value, 
                      levels=c('Female', 'Male',
                               'Right Colon', 'Left Colon', 'Rectum', 
                               'Stage 0/I', 'Stage II', 'Stage III', 'Stage IV',
                               'BMI < 25', 'BMI > 25',
                               'Age < 57', 'Age 57 - 64', 'Age 64 - 71', 
                                'Age > 71')))

# ##############################################################################
# Prediction bias for different metadata, species
g <- df.plot %>% 
  ggplot(aes(x=value, y=Predictions)) + 
    facet_grid(~confounder, scales = 'free', space = 'free') +
    xlab('') + 
    theme_classic() +
    theme(strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
          # axis.text.x = element_blank(),
          panel.grid.major.y = element_line(colour='lightgrey')) +
    ylab('Prediction Score') + 
    geom_boxplot()

ggsave(g, filename = paste0('../figures/', tag, 
                            '/figure_prediction_bias_',ml.method,'.pdf'),
       width = 7, height = 4, useDingbats=FALSE)

# ##############################################################################
# compute p-values
for (i in levels(df.plot$confounder)){
  df.temp <- df.plot %>% 
    filter(confounder==i) %>% 
    mutate(Study=factor(Study)) %>% 
    as.data.frame
  if (length(levels(df.temp$value)) != 2){
    t <- kruskal_test(Predictions~value|Study, data=df.temp)
  } else {
    t <- wilcox_test(Predictions~value|Study,data=df.temp)
  }
  print(i)
  print(pvalue(t))
}


# ##############################################################################
# Compute stage-specific TPR
roc.curve <- roc(cases=meta %>% 
                   filter(Group=='CRC') %>% 
                   pull(Predictions),
                 controls=meta %>% 
                   filter(Group=='CTR') %>% 
                   pull(Predictions))
threshold <- roc.curve$thresholds[which(roc.curve$specificities >= 0.9)[1]]
df.plot <- meta %>% 
  group_by(Stage) %>% 
  summarize(tp=sum(Predictions > threshold), n=n()) %>% 
  mutate(tpr=tp/n)

g.a <- df.plot %>% 
  filter(Stage != 'CTR') %>% 
  ggplot(aes(x=Stage, y=tpr, fill=Stage)) +
    geom_bar(stat="identity", col='black') + 
    theme_classic() + 
    scale_y_continuous(limits=c(0, 1)) +
    xlab('') +
    ylab('True positive rate') + 
    scale_fill_manual(
      values=c('#d2d2d2', '#9e9e9e', '#5d5d5d', '#292929'), guide=FALSE) + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          panel.grid.major.y=element_line(colour='lightgrey'))
ggsave(g.a, filename = paste0('../figures/', tag, 
                              '/figure_tpr_stage_',ml.method,'.pdf'),
       width = 3, height = 4)


cat('Successfully plotted prediction bias and TPR in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################