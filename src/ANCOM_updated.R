### 
# Downloaded from https://sites.google.com/site/siddharthamandal1985/research
# on 2019-01-28
# I commented out the examples so that the code can be sourced from 
#   other scripts, otherwise the code is unchanged - Jakob Wirbel
#
# EDIT: Added progressbar for better estimation of remaining run-time
###


library(exactRankTests)
library(nlme)
library(ggplot2)

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
   # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
    if( repeated==FALSE & adjusted == FALSE){
       if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
          tfun <- exactRankTests::wilcox.exact
       } else{
          tfun <- stats::kruskal.test
    }
    }else if( repeated==FALSE & adjusted == TRUE){
      tfun <- stats::aov
    }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
      tfun <- stats::friedman.test
    }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
      tfun <- nlme::lme
    }else if( repeated== TRUE & adjusted == TRUE){
      tfun <- nlme::lme
    }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  ### add progressbar
  pb <- txtProgressBar(max=((n_otu^2)/2), style=3)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      setTxtProgressBar(pb, pb$getVal()+1)
    }
    
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
              geom_histogram(binwidth=0.1,colour="black",fill="white") + 
              xlab("Proportion of zeroes") + ylab("Number of taxa") +
              theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}




# 
# ## Example
# ys_data=read.delim("~/Dropbox/Work at NIEHS/sid_reference_data/ys_genera.txt",header=T)
# 
# otu_test=ys_data[,c(1,5:665)]
# map_test=ys_data[,c(1:4)]
#   
# 
# comparison_test=ANCOM.main(OTUdat=otu_test,
#                            Vardat=map_test[which(map_test$SEX!="unknown"&map_test$COUNTRY!="GAZ:Venezuela"),],
#                            adjusted=TRUE,
#                            repeated=F,
#                            main.var="COUNTRY",
#                            adj.formula="SEX+AGE",
#                            repeat.var=NULL,
#                            multcorr=2,
#                            sig=0.05,
#                            prev.cut=0.90)
# 
# comparison_test$W.taxa
# 
# 
# 
# 
# ####Longitudinal test
# OTU_ANCOM_Test <- read.delim("~/Dropbox/Work at NIEHS/OTU_Ancom_Test_Counts.txt")
# 
# Mapping_File_ANCOM_Test <- read.delim("~/Dropbox/Work at NIEHS/OTU_Ancom_Test_Mapping.txt")
# 
# #data_test=merge(OTU_ANCOM_Test,Mapping_File_ANCOM_Test,by="Sample.ID")
# 
# longitudinal_comparison_test=ANCOM.main(OTUdat=OTU_ANCOM_Test,
#                                         Vardat=Mapping_File_ANCOM_Test[Mapping_File_ANCOM_Test$group_info%in%c("2FL","3SL","pHMO"),],
#                                         adjusted=F,
#                                         repeated=T,
#                                         main.var="group_info",
#                                         adj.formula=NULL,
#                                         repeat.var="day_of_treatment",
#                                         longitudinal=T,
#                                         random.formula="~1|host_subject_id",
#                                         multcorr=2,
#                                         sig=0.05,
#                                         prev.cut=0.90)
# 
# longitudinal_comparison_test$W.taxa
# 
# 
# longitudinal_comparison_test=ANCOM.main(OTUdat=OTU_ANCOM_Test,
#                                         Vardat=Mapping_File_ANCOM_Test[Mapping_File_ANCOM_Test$group_info%in%c("2FL","3SL","pHMO"),],
#                                         adjusted=F,
#                                         repeated=T,
#                                         main.var="group_info",
#                                         adj.formula=NULL,
#                                         repeat.var="day_of_treatment",
#                                         longitudinal=T,
#                                         random.formula="~1|host_subject_id",
#                                         multcorr=2,
#                                         sig=0.05,
#                                         prev.cut=0.90)
# 
# longitudinal_comparison_test$W.taxa
# 
# 
# 
# 
# longitudinal_comparison_test_HC_RA=ANCOM.main(OTUdat=otu_test,
#                                         Vardat=map_test[which(map_test$Group%in%c("HC","RA")),],
#                                         adjusted=F,
#                                         repeated=T,
#                                         main.var="Group",
#                                         adj.formula=NULL,
#                                         repeat.var="Visitorder",
#                                         longitudinal=T,
#                                         random.formula="~1|individual",
#                                         multcorr=2,
#                                         sig=0.05,
#                                         prev.cut=0.50)
# 
# longitudinal_comparison_test_HC_RA$W.taxa
# 
# 
# 
# longitudinal_comparison_test_HC_RR=ANCOM.main(OTUdat=otu_test,
#                                               Vardat=map_test[which(map_test$Group%in%c("HC","RR")),],
#                                               adjusted=F,
#                                               repeated=T,
#                                               main.var="Group",
#                                               adj.formula=NULL,
#                                               repeat.var="Visitorder",
#                                               longitudinal=T,
#                                               random.formula="~1|individual",
#                                               multcorr=2,
#                                               sig=0.05,
#                                               prev.cut=0.50)
# 
# longitudinal_comparison_test_HC_RR$W.taxa
# 
# longitudinal_comparison_test_RA_RR=ANCOM.main(OTUdat=otu_test,
#                                               Vardat=map_test[which(map_test$Group%in%c("RA","RR")),],
#                                               adjusted=F,
#                                               repeated=T,
#                                               main.var="Group",
#                                               adj.formula=NULL,
#                                               repeat.var="Visitorder",
#                                               longitudinal=T,
#                                               random.formula="~1|individual",
#                                               multcorr=2,
#                                               sig=0.05,
#                                               prev.cut=0.50)
# 
# longitudinal_comparison_test_RA_RR$W.taxa
# 
# 
# # plots_for_ancom=function(otu_list,otu_data,metadata,var){
# #   data_sub=otu_data[,which(colnames(otu_data)%in%c("Sample.ID",check))]
# #   meta_sub=metadata[,which(colnames(metadata)%in%c("Sample.ID",var))]
# #   
# #   data_for_plot=merge(data_sub,meta_sub,by="Sample.ID",all.x=T)
# #   
# #   
# # }
# #   
# # check=c("OTU_564806","OTU_360015","OTU_1062061","OTU_581079")
# # 
# 
# 
