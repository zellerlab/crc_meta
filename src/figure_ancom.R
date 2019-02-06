# ##############################################################################
#
## Compare Wilcoxon and ANCOM results
#
# ##############################################################################

library("matrixStats")
library("tidyverse")
library("ggpubr")
library("cowplot")


# ##############################################################################
# load data

# features
feat <- read.table('../data/species/feat_rel_crc.tsv', sep='\t', quote='',
                   stringsAsFactors = FALSE, check.names = FALSE)
# test results
test.results <- read_tsv('../files/species/all_measures.tsv')
# confounders
conf.results <- read_tsv('../files/species/confounder_table.tsv')
# ancom results
load('../files/species/ancom_results.RData')
ancom.threshold <- ancom.results$all %>% filter(detected_0.8) %>% 
  pull(W_stat) %>% min

# ##############################################################################
# plot ANCOM W vs Wilcoxon pvalue

df.plot <- tibble(species=rownames(feat),
                  p.adj=test.results$p.adj,
                  log.p.adj=-log10(test.results$p.adj + 1e-20),
                  ancom.effect=test.results$ancom.w,
                  marker.set=p.adj < 1e-05,
                  extended.marker.set=p.adj < 0.005,
                  ancom.called=test.results$ancom.w > ancom.threshold,
                  study.effect=conf.results$Study,
                  crc.effect=conf.results$disease)

g1 <- df.plot %>%
  filter(complete.cases(.)) %>% 
  mutate(type=case_when(ancom.called & marker.set ~ 'both',
                        ancom.called & !marker.set ~ 'ancom only',
                        !ancom.called & marker.set ~ 'wilcoxon only',
                        !ancom.called & !marker.set ~ 'neither')) %>% 
  mutate(type=factor(type, levels = 
                       c('ancom only', 'wilcoxon only', 'both', 'neither'))) %>% 
  ggplot(aes(x=ancom.effect, y=log.p.adj, col=type)) +
    geom_point() +
    theme(panel.grid.major.y = element_line(colour='lightgrey')) +  
    geom_hline(yintercept = 5, col='darkgrey') +
    geom_vline(xintercept = df.plot %>% 
                 filter(ancom.called == TRUE) %>% 
                 pull(ancom.effect) %>% min, col='darkgrey') +
    xlab('ANCOM - W') + 
    ylab('-log10(q-value)') + 
    scale_colour_manual(values=c('#BDCD00', '#006165', '#CC071E', 'grey'), 
                        name='') + 
    annotate('text', x=400, y=10, label=
               paste0("Spearman's rho: ", 
                      formatC(cor(df.plot$log.p.adj, 
                                  df.plot$ancom.effect, method='spearman', 
                                  use='pairwise.complete.obs')))) + 
    theme(legend.position = c(0.2, 0.8))

# ##############################################################################
# plot summed up relative abundance of differentially abundant makers

stopifnot(all(test.results$species == rownames(feat)))

df.plot.2 <- tibble(
  `Marker Set` = colSums(feat[df.plot %>% 
                                filter(marker.set) %>% 
                                pull(species),]),
  `Extended Marker Set` = colSums(feat[df.plot %>% 
                                         filter(extended.marker.set) %>% 
                                         pull(species),]),
  `ANCOM Set` = colSums(feat[df.plot %>% 
                               filter(ancom.called) %>% 
                               pull(species),]))

g2 <- df.plot.2 %>%  
  gather(key=key, value=value) %>%
  filter(key != 'Extended Marker Set') %>% 
  ggplot(aes(x=key, y=value)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width=0.2, alpha=0.4, stroke=0, size=2) + 
    theme(panel.grid.major.y = element_line(colour='lightgrey')) + 
    xlab('') + 
    ylab('Cumulative relative abundance') + 
    scale_y_continuous(limits=c(0, 1))

# ##############################################################################
# plot variance explained by disease and by study for extended and ANCOM set
df.plot.3 <- df.plot %>% filter(extended.marker.set | ancom.called)

# overlap to be plotted into a Venn diagramm
table(df.plot.3$ancom.called, df.plot.3$extended.marker.set)

g3 <- df.plot.3 %>% 
  mutate(type=case_when(extended.marker.set & ancom.called ~ 'both',
                        extended.marker.set & !ancom.called ~ 'Wilcoxon only',
                        !extended.marker.set & ancom.called ~ 'ANCOM only')) %>% 
  select(type, study.effect, crc.effect, marker.set) %>% 
  gather(key=panel, value=var, -type, -marker.set) %>% 
  ggplot(aes(x=type, y=var)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, aes(colour=marker.set)) + 
    facet_grid(panel~.) + 
    xlab('') + ylab('Variance explained') + 
    theme(panel.grid.major.y = element_line(colour='lightgrey')) + 
    scale_colour_manual(values=c('darkgrey', 'red'), guide=FALSE) + 
    stat_compare_means(comparisons = list(c('Wilcoxon only', 'both'),
                                          c('ANCOM only', 'both'),
                                          c('ANCOM only', 'Wilcoxon only')))


final.plot <- plot_grid(g1, g2, g3, labels='auto', nrow=1)
ggsave(final.plot, filename = '../figures/revision/ancom_figure.pdf', 
       width = 10, height = 6, useDingbats=FALSE)
