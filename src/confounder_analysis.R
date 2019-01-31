# ##############################################################################
#
##  Meta-Analysis Confounder Analysis
#
# ##############################################################################

# Packages
library("tidyverse")
library("coin")

# ##############################################################################
# general
cat('Starting confounder testing script\n')
start.time <- proc.time()[1]

set.seed(2018)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]

alpha.meta <- 1e-05


# ##############################################################################
# Get Data
fn.feat <- paste0('../data/', tag, '/feat_rel_crc.tsv')
feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE,
                                 stringsAsFactors = FALSE,
                                 check.names = FALSE, quote=''))

meta <- read_tsv('../data/meta/meta_crc.tsv')
stopifnot(all(meta$Sample_ID %in% colnames(feat.all)))

studies <- meta %>% pull(Study) %>% unique

feat.all <- feat.all[,meta$Sample_ID]

fn.pval <- paste0('../files/', tag, '/p_adj.tsv')
if (!file.exists(fn.pval)){
  stop("Please run the marker analysis script first. Exiting...\n")
}
p.vals.adj <- read.table(fn.pval, sep='\t', check.names = FALSE)

# ##############################################################################
# preprocess confounder variables to test later

meta <- meta %>%
  # age
  mutate(age_factor=as.factor(
    cut(meta$Age, breaks = quantile(meta$Age), labels=c(1,2,3,4)))) %>%
  # bmi
  mutate(bmi_factor=as.factor(
    cut(meta$BMI, breaks = c(0, 25, 30, 100),
        labels=c('lean', 'overweight', 'obese')))) %>%
  # library size
  mutate(lib_size_factor=as.factor(
    cut(meta$Library_Size, breaks = quantile(meta$Library_Size),
        labels=c(1,2,3,4))))

# ##############################################################################
#  variance explained by disease status
ss.disease <- apply(feat.all, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
  return(1-ss.o.i/ss.tot)
}, label=meta %>% pull(Group))

# calculate trimmed mean abundance
t.mean <- apply(feat.all, 1, mean, trim=0.1)

df.plot.all <- tibble(
  species=rownames(feat.all),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val=p.vals.adj[rownames(feat.all), 'all'],
  meta.significance=p.vals.adj[rownames(feat.all), 'all'] < 1e-05)

# ##############################################################################
# Test all possible confounder variables
df.list <- list()
for (meta.var in c('age_factor', 'Gender', 'bmi_factor', 'Study',
                   'lib_size_factor', 'Sampling_rel_to_colonoscopy',
                   'Diabetes', 'Vegetarian', 'Smoking')){

  cat('###############################\n', meta.var, '\n')
  meta.c <- meta %>%
    filter(!is.na(eval(parse(text=meta.var))))

  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$Group, meta.c %>% pull(meta.var)))
  print(table(meta.c$Study))
  feat.red <- feat.all[,meta.c$Sample_ID]

  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var

  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}
# ##############################################################################
# plot

g <- df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when(type=='Gender' ~ 'Sex',
                         type=='age_factor' ~ 'Age',
                         type=='bmi_factor' ~ 'BMI',
                         type=='lib_size_factor' ~ 'Library Size',
                         type=='Sampling_rel_to_colonoscopy' ~ 'Colonoscopy',
                         TRUE ~ type)) %>%
  ggplot(aes(x=disease, y=meta, size=t.mean+1e-08, col=meta.significance)) +
    geom_point(shape=19) +
    xlab('Variance explained by disease status') +
    ylab('Variance explained by metadata variable') +
    theme_bw() +
    facet_wrap(~facet, ncol=3) +
    theme(strip.background = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(from=0, to=0.32, by=0.1)) +
    scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
    scale_colour_manual(values = alpha(c('black', '#CC071E'),
                                       alpha=c(0.1, .75)),
                        name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
    scale_size_area(name='Trimmed mean\nabundance',
               breaks=c(1e-05, 1e-03, 1e-02)) +
    guides( size = "legend", colour='legend')

ggsave(g, filename = paste0('../figures/', tag, '/confounder_plot.pdf'),
       width = 8, height = 8, useDingbats=FALSE)

write_tsv(df.plot.all, path = paste0('../files/', tag, '/confounder_table.tsv'))

# ##############################################################################
# plot only study
df.plot.study <- df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  filter(type=='Study')

g2 <- df.plot.study %>%
  ggplot(aes(x=disease, y=meta)) +
    geom_point(aes(size=t.mean, fill=meta.significance), shape=21,
               col=alpha(c('black'), alpha=0.4)) +
    xlab(paste0('Variance explained by Disease\n',tag,' average: ',
                formatC(mean(df.plot.study$disease)*100, digits=2), '%')) +
    ylab(paste0('Variance explained by Study\n',tag,' average: ',
                formatC(mean(df.plot.study$meta)*100, digits=2), '%')) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(from=0, to=0.2, by=0.1)) +
    scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
    scale_fill_manual(values = alpha(c('grey', '#CC071E'),
                                       alpha=c(0.4, .8)),
                        name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
    scale_size_area(name='Trimmed mean abundance',
                    breaks=c(1e-05, 1e-03, 1e-02)) +
    guides( size = "legend", colour='legend')

ggsave(g2,
       filename = paste0('../figures/', tag, '/confounder_plot_study.pdf'),
       width = 6, height = 6)

cat('Successfully computed confounder effects in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################
