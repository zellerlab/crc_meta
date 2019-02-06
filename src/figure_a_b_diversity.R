# ##############################################################################
#
##  Alpha and Beta Diversity
#
# ##############################################################################

# Packages
library("tidyverse")
library("labdsv")
library("coin")
library("vegan")
library("yaml")
library("ggpubr")
library("cowplot")

parameters <- yaml.load_file('../parameters.yaml')
ref.studies <- parameters$ref.studies
study.cols <- unlist(parameters$plotting$study.cols)
group.cols <- unlist(parameters$plotting$group.cols)
start.time <- proc.time()[1]

# ##############################################################################
# Get Data
fn.feat <- '../data/species/feat_rar_crc.tsv'
feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE, 
                                 stringsAsFactors = FALSE, 
                                 check.names = FALSE, quote=''))
feat.all <- prop.table(feat.all)

meta <- read_tsv('../data/meta/meta_crc.tsv')

# ##############################################################################
# Compute Alpha diversity

df.div <- meta %>% 
  ### ### Richness
  # mutate(`All species`=colSums(feat.all != 0)) %>%
  # mutate(`reference mOTUs`=colSums(feat.all[grep(pattern='ref_mOTU', 
  #   x=rownames(feat.all), value=TRUE),] != 0)) %>%
  # mutate(`meta mOTUs`=colSums(feat.all[grep(pattern='meta_mOTU', 
  #   x=rownames(feat.all), value=TRUE),] != 0)) %>%
  ### ### Shannon
  mutate(`All species` = diversity(t(feat.all))) %>%
  mutate(`reference mOTUs` = diversity(t(feat.all[grep(pattern = 'ref_mOTU',
    x=rownames(feat.all), value = TRUE),]))) %>%
  mutate(`meta mOTUs` = diversity(t(feat.all[grep(pattern = 'meta_mOTU',
    x=rownames(feat.all), value = TRUE),]))) %>%
  
  select(Group, Study, `All species`, `reference mOTUs`, `meta mOTUs`) %>% 
  gather(key=type, value=diversity, -c(Study, Group)) %>% 
  mutate(Study=str_remove(Study, pattern='-CRC')) %>%
  mutate(Study=factor(Study, 
                      levels = str_remove(ref.studies, pattern='-CRC'))) %>%
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>% 
  mutate(type=factor(type, levels=c('All species', 
                                      'reference mOTUs', 'meta mOTUs')))

# blocked wilcoxon
df.div %>% 
  filter(type=='All species') %>% 
  wilcox_test(diversity~Group|Study, data=.)
df.div %>% 
  filter(type=='reference mOTUs') %>% 
  wilcox_test(diversity~Group|Study, data=.)
df.div %>% 
  filter(type=='meta mOTUs') %>% 
  wilcox_test(diversity~Group|Study, data=.)

# anova
summary(aov(rank~Group+Study, 
            data=df.div %>% 
              filter(type=='All species') %>% 
              mutate(rank=rank(diversity))))
summary(aov(rank~Group+Study, 
            data=df.div %>% 
              filter(type=='reference mOTUs') %>% 
              mutate(rank=rank(diversity))))
summary(aov(rank~Group+Study, 
            data=df.div %>% 
              filter(type=='meta mOTUs') %>% 
              mutate(rank=rank(diversity))))
# ##############################################################################
# plot

g <- df.div %>% 
  ggplot(aes(x=Study, fill=Group, y=diversity)) +
    geom_boxplot() +
    facet_wrap(~type, scales = 'fixed') + 
    theme_bw() + 
    ylab('Alpha diversity (Shannon)') +
    xlab('') +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank()) + 
    scale_fill_manual(values=group.cols) + 
    stat_compare_means(label = "p.format", method = "wilcox.test", 
                       size=2)
ggsave(g, filename = '../figures/species/alpha_diversity.pdf', 
       width = 8, height = 3, useDingbats=FALSE)    

cat('Successfully plotted alpha diversity in',
    proc.time()[1]-start.time, 'second...\n')

# ##############################################################################
# Beta Diversity
# species
feat.all <- read.table('../data/species/feat_ab_crc.tsv', sep='\t',
                       stringsAsFactors = FALSE, check.names = FALSE,
                       quote='')


# ##############################################################################
# compute PCoA
dist = vegdist(t(feat.all), method = 'bray')
pco.results = pco(dist, k=2)

axis.1.title <- paste('PCo1 [', 
                      round((pco.results$eig[1]/sum(pco.results$eig))*100,1),
                      '%]', sep='')
axis.2.title <- paste('PCo2 [', 
                      round((pco.results$eig[2]/sum(pco.results$eig))*100,1),
                      '%]', sep='')

df.plot <- tibble(Axis1 = -1*pco.results$points[,1],
                  Axis2 = pco.results$points[,2],
                  Sample_ID = rownames(pco.results$points),
                  Group=meta$Group,
                  Study=meta$Study)

# ##############################################################################
# subplots

# main plot
g.main <- df.plot %>% 
  ggplot(aes(x=Axis1, y=Axis2, shape=Group, col=Study)) +
  geom_point() + 
  scale_colour_manual(values=study.cols, guide=FALSE) + 
  scale_shape_manual(values=c(19, 1), guide=FALSE) + 
  scale_x_continuous(position='top') +
  xlab(axis.1.title) + ylab(axis.2.title) +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        axis.ticks=element_blank(), axis.text = element_blank(),
        panel.grid = element_blank())

# study boxplot axis 1
g.s.1 <- df.plot %>% 
  mutate(Study=factor(Study, levels=names(study.cols))) %>% 
  ggplot(aes(y=Axis1, x=Study, fill=Study)) + 
  geom_boxplot() + 
  scale_fill_manual(values=study.cols, guide=FALSE) +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank()) + 
  coord_flip()

# study boxplot axis 2
g.s.2 <- df.plot %>% 
  mutate(Study=factor(Study, levels=names(study.cols))) %>% 
  ggplot(aes(y=Axis2, x=Study, fill=Study)) + 
  geom_boxplot() + 
  scale_fill_manual(values=study.cols, guide=FALSE) +
  scale_x_discrete(position='top') +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

# group plot axis1
g.g.1 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis1, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.cols, guide=FALSE) + 
  ylab(axis.1.title) + 
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) + 
  coord_flip()
# group plot axis2
g.g.2 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis2, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.cols, guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  ylab(axis.2.title) + 
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())

# ##############################################################################
# Plot everyting together
pdf('../figures/species/beta_diversity.pdf', useDingbats = FALSE)
plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL, NULL, g.g.1, NULL, NULL,
          nrow=3,
          rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))
dev.off()

cat('Successfully plotted beta diversity in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################