# ##############################################################################
#
##  Meta-Analysis Marker Analysis, Complete pipeline
#
# ##############################################################################

# Packages
library("tidyverse")
library("coin")
library("pROC")
library("yaml")

# ##############################################################################
# general 
cat('Starting association checking script\n')
start.time <- proc.time()[1]

set.seed(2018)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("The analysis tag needs to be provided! Exiting...\n")
}
tag <- args[1]


parameters <- yaml.load_file('../parameters.yaml')
log.n0 <- ifelse(tag=='species',
                 as.numeric(parameters$associations$log.n0),
                 as.numeric(parameters$associations$log.n0.func))
mult.corr <- parameters$associations$mult.corr

# ##############################################################################
# Get Data
fn.feat <- paste0('../data/', tag, '/feat_rel_crc.tsv')
feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE, 
                                 stringsAsFactors = FALSE, 
                                 check.names = FALSE, quote=''))

meta <- read_tsv('../data/meta/meta_crc.tsv')
stopifnot(all(meta$Sample_ID %in% colnames(feat.all)))

studies <- meta %>% pull(Study) %>% unique

# ##############################################################################
# block for colonoscopy and study as well
meta <- meta %>%
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

feat.all <- feat.all[,meta$Sample_ID]

# ##############################################################################
# calculate pval, gFC, AUROC

p.val <- matrix(NA, nrow=nrow(feat.all), ncol=length(studies)+1, 
                dimnames=list(row.names(feat.all), c(studies, 'all')))
fc <- p.val
aucs.mat <- p.val
aucs.all  <- vector('list', nrow(feat.all))

cat("Calculating effect size for every feature...\n")
pb <- txtProgressBar(max=nrow(feat.all), style=3)

# caluclate wilcoxon test and effect size for each feature and study
for (f in row.names(feat.all)) {
  
  # for each study
  for (s in studies) {
    
    x <- feat.all[f, meta %>% filter(Study==s) %>% 
                    filter(Group=='CRC') %>% pull(Sample_ID)]
    y <- feat.all[f, meta %>% filter(Study==s) %>% 
                    filter(Group=='CTR') %>% pull(Sample_ID)]
    
    # Wilcoxon
    p.val[f,s] <- wilcox.test(x, y, exact=FALSE)$p.value
    
    # AUC
    aucs.all[[f]][[s]] <- c(roc(controls=y, cases=x, 
                                direction='<', ci=TRUE, auc=TRUE)$ci)
    aucs.mat[f,s] <- c(roc(controls=y, cases=x, 
                           direction='<', ci=TRUE, auc=TRUE)$ci)[2]
    
    # FC
    q.p <- quantile(log10(x+log.n0), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+log.n0), probs=seq(.1, .9, .05))
    fc[f, s] <- sum(q.p - q.n)/length(q.p)
  }
  
  # calculate effect size for all studies combined
  # Wilcoxon + blocking factor
  d <- data.frame(y=feat.all[f,], 
                  x=meta$Group, block=meta$block)
  p.val[f,'all'] <- pvalue(wilcox_test(y ~ x | block, data=d))
  # other metrics
  x <- feat.all[f, meta %>% filter(Group=='CRC') %>% pull(Sample_ID)]
  y <- feat.all[f, meta %>% filter(Group=='CTR') %>% pull(Sample_ID)]
  # FC
  fc[f, 'all'] <- mean(fc[f, studies])
  # AUC
  aucs.mat[f,'all'] <- c(roc(controls=y, cases=x, 
                             direction='<', ci=TRUE, auc=TRUE)$ci)[2]
  
  # progressbar
  setTxtProgressBar(pb, (pb$getVal()+1))
}
cat('\n')

# multiple hypothesis correction
p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method=mult.corr),
                    check.names = FALSE)

# ##############################################################################
# save results
write.table(fc, file=paste0('../files/', tag, '/fc.tsv'),
            quote=FALSE, sep='\t')
write.table(p.adj, file = paste0('../files/', tag, '/p_adj.tsv'), 
            quote=FALSE, sep='\t')
write.table(p.val, file = paste0('../files/', tag, '/p_val.tsv'), 
            quote=FALSE, sep='\t')
write.table(aucs.mat, file = paste0('../files/', tag, '/aucs.tsv'), 
            quote=FALSE, sep='\t')
save(aucs.all, file=paste0('../files/', tag, '/aucs_all.RData'))

# ##############################################################################
# compute ANCOM's W
if (tag == 'species'){
  source('./ANCOM_updated.R') ## updated ANCOM
  ancom.results <- list()
  ancom.w.mat <- matrix(NA, nrow=nrow(p.val), ncol=ncol(p.val)+1,
                        dimnames = list(rownames(p.val), 
                                        c(colnames(p.val), 'all-blocked')))
  fn.feat <- '../data/species/feat_ab_all.tsv' 
  ## ancom seems to need absolute abudnances
  feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE, 
                                   stringsAsFactors = FALSE, 
                                   check.names = FALSE, quote=''))
  for (s in studies){
    cat('Calculating ANCOM W for study', s, '\n')
    meta.red <- meta %>% filter(Study==s)
    feat.red <- t(feat.all[,meta.red$Sample_ID])

    data.test <- data.frame(feat.red)
    sample.names <- data.frame(Sample.ID=rownames(feat.red))
    data.test <- cbind(sample.names, data.test)

    meta.test <- cbind(sample.names, meta.red %>%
                         select(Age, Gender, BMI, Study, Group) %>%
                         as.data.frame)
    temp <- ANCOM.main(OTUdat = data.test, Vardat = meta.test, adjusted=FALSE,
                       repeated = FALSE, main.var='Group', adj.formula = NULL,
                       repeat.var=NULL, longitudinal = FALSE,
                       random.formula = NULL, multcorr = 2, sig=0.1,
                       prev.cut = 0.99)

    ancom.results[[s]] <- temp$W.taxa
    ancom.w.mat[match(temp$W.taxa$otu.names,
                      make.names(rownames(ancom.w.mat))),s] <-
      temp$W.taxa$W_stat
    cat('\nFinished calculations for study', s, '\n')
  }

  feat.red <- t(feat.all[, meta$Sample_ID])
  data.test <- data.frame(feat.red)
  sample.names <- data.frame(Sample.ID=rownames(feat.red))
  data.test <- cbind(sample.names, data.test)
  
  meta.test <- cbind(sample.names, meta %>% 
                       select(Age, Gender, BMI, Study, Group, block) %>% 
                       as.data.frame)
  # blocked all
  temp <- ANCOM.main(OTUdat = data.test, Vardat = meta.test, adjusted=TRUE, 
                     repeated = FALSE, main.var='Group', adj.formula = 'block', 
                     repeat.var=NULL, longitudinal = FALSE, 
                     random.formula = NULL, multcorr = 2, sig=0.05, 
                     prev.cut = 0.99)
  ancom.results[['all']] <- temp$W.taxa
  ancom.w.mat[match(temp$W.taxa$otu.names, 
                    make.names(rownames(ancom.w.mat))),'all'] <- 
    temp$W.taxa$W_stat
  
  # save ancom results
  write.table(ancom.w.mat, file=paste0('../files/', tag, '/ancom_w.tsv'),
              quote=FALSE, sep='\t')
  save(ancom.results, file = paste0('../files/', tag, '/ancom_results.RData'))
}

# ##############################################################################
# save results for meta-analysis in supplementary table
df.all <- tibble(species=rownames(feat.all),
                 fc=fc[,'all'],
                 auroc=aucs.mat[,'all'],
                 p.val=p.val[,'all'],
                 p.adj=p.adj[,'all'])
if (tag == 'species'){
  df.all <- df.all %>% 
    mutate(ancom.w=ancom.w.mat[,'all'])
  
  # add everything together for supplementary material
  auc <- aucs.mat
  ancom <- ancom.w.mat
  df.everything <- tibble(species=rownames(feat.all))
  for (metric in c('p.val', 'p.adj', 'fc', 'auc', 'ancom')){
    temp <- eval(parse(text=paste0(metric,' %>% as.data.frame()'))) 
    colnames(temp) <- paste0(metric, '.', colnames(temp))
    df.everything <- bind_cols(df.everything, temp)
  }
  write_tsv(df.everything, path='../files/species/all_measures_all.tsv')
}

write_tsv(df.all, path=paste0('../files/', tag, '/all_measures.tsv'))

cat('Successfully computed associations measures in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################
