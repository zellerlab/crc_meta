# ######################################################################
#
## Analyze assocations between clusters and meta-variables
#
# ######################################################################

require("coin")

cluster.associations <- function(features, meta.data, conf, 
                                 block=NULL, clustering=NULL){
  
  # assert compatability
  stopifnot(all(meta.data$Sample_ID %in% colnames(features)))
  meta.data <- meta.data %>% filter(Sample_ID %in% colnames(features))
  stopifnot(unique(c(features)) == c(0, 1))
  
  # take out all without conf data
  meta.data <- meta.data %>% 
    filter_(paste0("!is.na(", conf, ")"))
  if (!is.null(block)) {
    print(table(meta.data[[block]], meta.data[[conf]]))
    stopifnot(all(!is.na(meta.data[[block]])))
  } else {
    print('No blocking factor chosen!!!...')
  }
  features.red <- features[,meta.data$Sample_ID]
  
  if (!is.null(clustering)){
    features.red <- t(sapply(unique(clustering), FUN=function(x){
      as.numeric(colSums(features.red[which(clustering == x),]) >= 1)
    }))
    rownames(features.red) <- c('LP Markers', 'Porphyromonas', 
                                'HP Markers', 'Clostridiales')
  }

  if (is.null(block)){
    block.vec <- NULL
  } else {
    block.vec <- meta.data[[block]]
  }
  # compute p-vals
  if (meta.data %>% pull(conf) %>% unique() %>% length < 6){
    
    p.vals <- apply(features.red, 1, FUN=function(x, c, b){
      d <- data.frame(y=as.factor(x),
                      x=as.factor(c))
      if (is.null(b)){
        t <- chisq_test(y~x, data=d)
        return(pvalue(t))
      } else {
        d$block <- as.factor(b)
        t <- cmh_test(y~x|block, data=d)
        return(pvalue(t))
      }
    }, 
    c=meta.data[[conf]], 
    b=block.vec)
  } else {
    p.vals <- apply(features.red, 1, FUN=function(x, c, b){
      d <- data.frame(y=as.numeric(x),
                      x=as.numeric(c))
      if (is.null(b)){
        t <- cor.test(d$x, d$y, method='spearman', exact = FALSE)
        return(t$p.value)
      } else {
        d$block <- as.factor(b)
        t <- spearman_test(y ~x | block, data=d)
        return(pvalue(t))
      }
    }, c=meta.data[[conf]], b=block.vec)
  }
  print(p.vals)
  
  # prepare plot
  df.plot <- as_tibble(t(features.red)) %>% 
    mutate(conf=meta.data[[conf]]) %>% 
    gather(key=key, value=group, -conf) %>% 
    group_by(conf, key) %>% 
    summarize(value=sum(group)/n()) %>% 
    mutate(pval=p.vals[key])
  
  # plot
  g <- df.plot %>% 
    ggplot(aes(x=key, y=value, fill=conf)) + 
      geom_bar(stat='identity', position = position_dodge()) +
      theme_classic() + xlab('') + 
      ylab('Positive fraction') + scale_y_continuous(limits = c(0, 1))
  if (any(df.plot$pval < 0.1)){
    g <- g + geom_text(aes(y=.9),
                       label=formatC(df.plot$pval, digits = 2, 
                                     format='E'), angle=90)
  }
  return(g)
}
