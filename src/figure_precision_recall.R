# ##############################################################################
#
##  Precision Recall Plots for each study with the Meta-Analysis as 
##    "Common Truth"
##  Use two different truth sets with alpha.meta and alpha.single.study
#
# ##############################################################################

# ##############################################################################
# Packages
library("tidyverse")
library("SIAMCAT")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')
study.cols <- unlist(parameters$plotting$study.cols)
alpha.meta <- as.numeric(parameters$associations$alpha.meta)
alpha.single <- as.numeric(parameters$associations$alpha.single.study)

# ##############################################################################
# get data
p.adj <- read.table('../files/species/p_adj.tsv', sep='\t', quote='',
                    check.names = FALSE, stringsAsFactors = FALSE)
p.val <- read.table('../files/species/p_val.tsv', sep='\t', quote='',
                    check.names = FALSE, stringsAsFactors = FALSE)
fc.mat <- read.table('../files/species/fc.tsv', sep='\t', quote='',
                     check.names = FALSE, stringsAsFactors = FALSE)

marker.set <- rownames(p.adj)[p.adj$all < alpha.meta]
marker.set.extended <- rownames(p.adj)[p.adj$all < alpha.single]

# ##############################################################################
# calculate PR curves for all
df.plot <- tibble()
for (study in parameters$ref.studies){
  df.temp <- tibble(
    predictor=-log10(p.val[[study]] + 1e-20),
    fc.study=ifelse(sign(fc.mat[[study]]) != 0, sign(fc.mat$all), 1),
    fc.all=ifelse(sign(fc.mat$all) != 0, sign(fc.mat$all), 1),
    response1=as.numeric(rownames(p.adj) %in% marker.set),
    response2=as.numeric(rownames(p.adj) %in% marker.set.extended)
  )
  # remove NAs
  df.temp <- df.temp %>% 
    filter(!is.na(df.temp$predictor)) %>% 
    # set positive cases with opposite FCs to negative cases
    mutate(response1=case_when(fc.study!=fc.all ~ 0,
                               fc.study==fc.all ~ response1),
           response2=case_when(fc.study!=fc.all ~ 0,
                               fc.study==fc.all ~ response2))
  
  # marker set
  temp <- SIAMCAT:::evaluate.classifier(predictions = df.temp$predictor,
                                        test.label=df.temp$response1, 
                                        label=list(info=c('test'=1, 'neg'=0)))
  temp.2 <- SIAMCAT:::evaluate.get.pr(temp)
  threshold.0.1 <- abs(temp$thresholds + log10(alpha.single))
  threshold.0.01 <- abs(temp$thresholds + log10(alpha.meta))

  df.plot <- bind_rows(
    df.plot, tibble(
      recall=temp.2$recall, precision=temp.2$precision, 
      threshold=temp$thresholds,
      study=study, set='Meta-FDR 1e-05',
      set.0.1=threshold.0.1 == min(threshold.0.1),
      set.0.01=threshold.0.01 == min(threshold.0.01)
    ))
  
  # extended marker set
  temp <- SIAMCAT:::evaluate.classifier(predictions = df.temp$predictor,
                                        test.label=df.temp$response2, 
                                        label=list(info=c('test'=1, 'neg'=0)))
  temp.2 <- SIAMCAT:::evaluate.get.pr(temp)
  threshold.0.1 <- abs(temp$thresholds + log10(alpha.single))
  threshold.0.01 <- abs(temp$thresholds + log10(alpha.meta))
  
  df.plot <- bind_rows(
    df.plot, tibble(
      recall=temp.2$recall, precision=temp.2$precision, 
      threshold=temp$thresholds,
      study=study, set='Meta-FDR 0.005',
      set.0.1=threshold.0.1 == min(threshold.0.1),
      set.0.01=threshold.0.01 == min(threshold.0.01)
    ))
}

# ##############################################################################
# plot
df.plot <- df.plot %>% 
  mutate(study=factor(study, levels = names(study.cols)))

g1 <- df.plot %>% 
  ggplot(aes(x=recall, y=precision, col=study)) + 
    theme_bw() + 
    theme(strip.background = element_blank(), 
          panel.grid.minor = element_blank()) + 
    geom_line() +
    facet_grid(set~study) + 
    scale_color_manual(values=study.cols, guide=FALSE) + 
    xlab('Recall') + ylab('Precision')

# add points for different cutoffs
g1 <- g1 + 
  geom_point(data=df.plot %>% filter(set.0.1), 
             aes(fill=study), colour='black', shape=21, size=2) + 
  geom_point(data=df.plot %>% filter(set.0.01), 
             aes(fill=study), colour='black', shape=23, size=2) + 
  scale_fill_manual(values = study.cols, guide=FALSE) 
  

ggsave(g1, filename = '../figures/species/precision_recall.pdf',
       useDingbats=FALSE, width = 11.69, height = 6)
