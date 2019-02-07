# ##############################################################################
#
##  Compute GMMs
#
# ##############################################################################

# Packages
library("tidyverse")
library("coin")
library("yaml")

parameters <- yaml.load_file('../parameters.yaml')
data.loc <- parameters$setup$dataset
ref.studies <- parameters$ref.studies
log.n0 <- as.numeric(parameters$associations$log.n0.func)

module.info <- read.table('../files/GMM/GMM.hierarchy.v1.07.tsv', 
                          sep='\t', header = TRUE)
# ##############################################################################
# Get Data

# meta
meta <- read_tsv('../data/meta/meta_crc.tsv')
  # filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  # mutate(block=ifelse(Study!='CN-CRC', Study, 
                      # paste0(Study, '_', Sampling_rel_to_colonoscopy)))

# KEGG
feat.kegg <- read.table('../data/KEGG/feat_rel_crc.tsv', sep='\t', quote='',
                       stringsAsFactors = FALSE, check.names = FALSE)

# ##############################################################################
# compute GMMs

# MAJOR HACK! the gmm tool does not seem to work otherwise...
#             and the first KO is not in any GMM, so it should be fine
feat.kegg[1,] <- colnames(feat.kegg)
write.table(feat.kegg, 
            file = '../files/GMM/temp/test.tsv',  # any temporary file
            sep='\t', quote=FALSE, col.names = FALSE)

# GMM call, to be run with the GMM tool in the command line
## java -jar gmms.jar -a 2 -d ./GMMs.v1.07.txt 
##    -i PATH_TO_FOLDER/files/GMM/temp/test.tsv 
##    -o PATH_TO_FOLDER/files/GMM/gmm_output/

# ##############################################################################
# get GMM profiles into one matrix
sample.files <- list.files('../files/GMM/gmm_output/', pattern='modules')

module.mat <- matrix(NA, nrow=nrow(module.info), ncol=length(sample.files),
                     dimnames=list(module.info$Module, 
                                   str_replace(sample.files, '.modules', '')))
for (s in sample.files){
  sample <- str_replace(s, '.modules', '')
  temp <- read.table(paste0('../files/GMM/gmm_output/', s), sep='\t', 
                     header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  module.mat[rownames(temp),sample] <- temp$Value
}

# ##############################################################################
# compute pvalues and FCs
GMM <- module.mat[,meta$Sample_ID]

p.val <- matrix(NA, nrow=nrow(GMM), ncol=length(ref.studies)+1, 
                dimnames=list(row.names(GMM), c(ref.studies, 'all')))
fc <- p.val

# caluclate wilcoxon test and effect size for each feature and study
for (f in row.names(GMM)) {
  
  if (sum(!is.na(GMM[f,])) < 10){
    next()
  }
  
  # for each study
  for (s in ref.studies) {
    
    x <- GMM[f, meta %>% filter(Study==s) %>% 
                    filter(Group=='CRC') %>% pull(Sample_ID)]
    y <- GMM[f, meta %>% filter(Study==s) %>% 
                    filter(Group=='CTR') %>% pull(Sample_ID)]
    
    # Wilcoxon
    p.val[f,s] <- wilcox.test(x, y, exact=FALSE)$p.value
    
    # FC
    q.p <- quantile(log10(x+log.n0), probs=seq(.1, .9, .05), na.rm=TRUE)
    q.n <- quantile(log10(y+log.n0), probs=seq(.1, .9, .05), na.rm=TRUE)
    fc[f, s] <- sum(q.p - q.n)/length(q.p)
  }
  
  # calculate effect size for all studies combined
  # Wilcoxon + blocking factor
  d <- data.frame(y=GMM[f,], 
                  x=meta$Group, block=meta$Study)
  p.val[f,'all'] <- pvalue(wilcox_test(y ~ x | block, data=d))
  
  # FC
  fc[f, 'all'] <- mean(fc[f, ref.studies])
}

# multiple hypothesis correction
p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method='fdr'),
                    check.names = FALSE)
p.adj <- p.adj[!is.na(p.adj$all),]
fc <- fc[rownames(p.adj),]

# ##############################################################################
# save pvalues and FCs
write.table(fc, file = '../files/GMM/fc.tsv', sep='\t', quote=FALSE,
            row.names = TRUE, col.names = TRUE)
write.table(p.adj, file = '../files/GMM/p_adj.tsv', sep='\t', quote=FALSE,
            row.names = TRUE, col.names = TRUE)

# ##############################################################################
# save GMM matrix

module.mat[is.na(module.mat)] <- 0
GMM <- module.mat[,meta$Sample_ID]


write.table(GMM, file = '../data/GMM/feat_rel_crc.tsv', sep='\t',
            quote=FALSE, row.names=TRUE, col.names=TRUE)
