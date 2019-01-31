# ##############################################################################
#
##  Data preparation script
#
# ##############################################################################

# Packages
library('vegan')
library('tidyverse')
library('matrixStats')
library('yaml')

cat('Starting data cleaning script\n')
start.time <- proc.time()[1]
# ##############################################################################
# Parameters
set.seed(2018)

parameters <- yaml.load_file('../parameters.yaml')
data.loc <- parameters$setup$data.location
v.tag <- parameters$setup$tag
ab.cutoff <- as.numeric(parameters$setup$ab.cutoff)
crc.studies <- parameters$ref.studies

# ##############################################################################
# read in raw data

# metadata
fn.meta <- paste0(data.loc, 'metadata/meta_all.tsv')
meta.all <- read_tsv(fn.meta)

meta.crc <- meta.all %>%
  filter(Study %in% crc.studies) %>%
  filter(Group %in% c('CTR', 'CRC'))

# features
fn.feat <- paste0(data.loc, 'mOTU_profiles/species_profiles_',
                  v.tag, '_motus2.0.0.tsv')
feat.ab.all <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                        header = TRUE, check.names = FALSE, row.names = 1,
                        quote='')
feat.ab.all = as.matrix(feat.ab.all)
cat('Successfully loaded data...\n')
cat('Feature dimensions:', dim(feat.ab.all), '\n')

# checks
stopifnot(all(meta.crc$Sample_ID %in% colnames(feat.ab.all)))
feat.ab.crc <- feat.ab.all[,meta.crc$Sample_ID]

print(table(meta.crc$Group))
# CRC CTR
# 285 290
print(table(meta.crc$Study))
# AT-CRC CN-CRC DE-CRC FR-CRC US-CRC
#    109    128    120    114    104
print(table(meta.crc$Study, meta.crc$Group))
#        CRC CTR
# AT-CRC  46  63
# CN-CRC  74  54
# DE-CRC  60  60
# FR-CRC  53  61
# US-CRC  52  52

# ##############################################################################
# convert to relative abundance
feat.rel.crc = prop.table(feat.ab.crc, 2)
lib.size = 5000
feat.rar.crc = t(rrarefy(floor(t(feat.ab.crc)), sample=lib.size))

temp.max.ab = t(sapply(row.names(feat.rel.crc),
  FUN=function(marker){sapply(unique(meta.crc$Study),
    FUN=function(study, marker){
      max.ab = max(feat.rel.crc[marker, which(meta.crc$Study == study)])
    },
  marker=marker)}))

### low abundance filter
f.idx = rowSums(temp.max.ab >= ab.cutoff) >= 3 &
  row.names(feat.rel.crc) != '-1'

feat.ab.red = feat.ab.crc[f.idx,]
feat.rel.red = feat.rel.crc[f.idx,]
feat.rar.red = feat.rar.crc[f.idx,]
cat('Retaining', sum(f.idx), 'features after low-abundance filtering...\n')

# ##############################################################################
# save cleaned data
fn.tax.ab <- '../data/species/feat_ab_crc.tsv'
write.table(feat.ab.red, file=fn.tax.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

fn.tax.rel.ab <- '../data/species/feat_rel_crc.tsv'
write.table(feat.rel.red, file=fn.tax.rel.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

fn.tax.rar.ab <- '../data/species/feat_rar_crc.tsv'
write.table(feat.rar.red, file=fn.tax.rar.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)


# ##############################################################################
# meta
write_tsv(meta.crc, path = '../data/meta/meta_crc.tsv')

# save other metdata as well
crc.ext.studies <- parameters$ext.crc.studies
meta.crc.ext <- meta.all %>% filter(Study %in% crc.ext.studies)
write_tsv(meta.crc.ext, path='../data/meta/meta_crc_ext.tsv')


ext.studies <- parameters$ext.studies
meta.ext <- meta.all %>% filter(Study %in% ext.studies)

# adjust Sample IDs 
# a bit Hacky....
meta.ext <- meta.ext %>% 
  mutate(Sample_ID = str_replace(Sample_ID, 
                                 pattern = '\\.0$', replacement = '-0')) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, 
                                 pattern = '\\.1$', replacement = '-1')) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, 
                                 pattern = '\\.2$', replacement = '-2')) %>% 
  mutate(Sample_ID = str_replace(Sample_ID, 
                                 pattern = '\\.4$', replacement = '-4')) %>% 
  filter(Sample_ID %in% colnames(feat.ab.all))
write_tsv(meta.ext, path='../data/meta/meta_ext.tsv')

# ##############################################################################
# save cleaned data for other studies

# external CRC studies
feat.rel.crc.ext <- prop.table(feat.ab.all, 2)[rownames(feat.ab.red), 
                                               meta.crc.ext$Sample_ID]
fn.tax.rel <- '../data/species/feat_rel_crc_ext.tsv'
write.table(feat.rel.crc.ext, file=fn.tax.rel, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

# external other studies
feat.rel.ext <- prop.table(feat.ab.all, 2)[rownames(feat.ab.red), 
                                           meta.ext$Sample_ID]
fn.tax.rel <- '../data/species/feat_rel_ext.tsv'
write.table(feat.rel.ext, file=fn.tax.rel, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

# ##############################################################################
# label
### save data in appropriate format
fn.label ='../data/label/label_crc.tsv'
label.header = '#BINARY:1=CRC;-1=CTR'
label = rep(-1, nrow(meta.crc))
names(label) = meta.crc$Sample_ID
label[meta.crc$Group=='CRC'] = +1
write(label.header, file=fn.label, append=FALSE)
suppressWarnings(write.table(t(label), file=fn.label, quote=FALSE, sep='\t',
                             row.names=FALSE, col.names=TRUE, append=TRUE))


# ##############################################################################
# same for functional data
## KEGG
fn.kegg <- paste0(data.loc, 'functional_profiles/KEGG_profiles_crc_meta.tsv')
kegg.all <- as.matrix(read.csv(fn.kegg, sep='\t', stringsAsFactors = FALSE,
                               header = TRUE, check.names = FALSE,
                               row.names = 1))
stopifnot(all(meta.crc$Sample_ID %in% colnames(kegg.all)))
kegg.all <- kegg.all[,meta.crc$Sample_ID]

temp.max.ab <- sapply(unique(meta.crc$Study), FUN=function(s){
  rowMaxs(kegg.all[,meta.crc %>% filter(Study == s) %>% pull(Sample_ID)])
})

### low abundance filter
f.idx = rowSums(temp.max.ab >= 1e-06) >= 3
kegg.red <- kegg.all[f.idx,]
cat('Retaining', sum(f.idx), 'KEGG features after low-abundance filtering...\n')


fn.kegg.rel.ab <- '../data/KEGG/feat_rel_all.tsv'
write.table(kegg.red, file=fn.kegg.rel.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

## eggNOG
fn.eggnog <- paste0(data.loc,
                    'functional_profiles/eggNOG_profiles_crc_meta.tsv')
eggnog.all <- as.matrix(read.csv(fn.eggnog, sep='\t',
                                 stringsAsFactors = FALSE,
                                 header = TRUE, check.names = FALSE,
                                 row.names = 1))
stopifnot(all(meta.crc$Sample_ID %in% colnames(eggnog.all)))
eggnog.all <- eggnog.all[,meta.crc$Sample_ID]

temp.max.ab <- sapply(unique(meta.crc$Study), FUN=function(s){
  rowMaxs(eggnog.all[,meta.crc %>% filter(Study == s) %>% pull(Sample_ID)])
})

### low abundance filter
f.idx = rowSums(temp.max.ab >= 1e-06) >= 3
eggnog.red <- eggnog.all[f.idx,]
cat('Retaining', sum(f.idx), 
    'eggNOG features after low-abundance filtering...\n')

fn.eggnog.rel.ab <- '../data/eggNOG/feat_rel_all.tsv'
write.table(eggnog.red, file=fn.eggnog.rel.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

# ##############################################################################
# download and save the functional profiles of external studies as well

# KEGG
fn.feat <- paste0(data.loc, 
                  'functional_profiles/KEGG_profiles_crc_external.tsv')
kegg.ext <- read.csv(file = fn.feat, sep='\t', stringsAsFactors = FALSE,
                     row.names = 1, check.names = FALSE)
stopifnot(all(meta.crc.ext$Sample_ID %in% colnames(kegg.ext)))
kegg.ext.red <- kegg.ext[rownames(kegg.red), meta.crc.ext$Sample_ID]
fn.kegg.rel.ab <- '../data/KEGG/feat_rel_crc_ext.tsv'
write.table(kegg.ext.red, file=fn.kegg.rel.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

# eggNOG
fn.feat <- paste0(data.loc, 
                  'functional_profiles/eggNOG_profiles_crc_external.tsv')
eggnog.ext <- read.csv(file = fn.feat, sep='\t', stringsAsFactors = FALSE,
                       row.names = 1, check.names = FALSE)
stopifnot(all(meta.crc.ext$Sample_ID %in% colnames(eggnog.ext)))
eggnog.ext.red <- eggnog.ext[rownames(eggnog.red), meta.crc.ext$Sample_ID]
fn.eggnog.rel.ab <- '../data/eggNOG/feat_rel_crc_ext.tsv'
write.table(eggnog.ext.red, file=fn.eggnog.rel.ab, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

cat('Successfully cleaned and saved data in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################
