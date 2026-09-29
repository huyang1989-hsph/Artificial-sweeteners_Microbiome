homepath = "~project/AS/PR"
setwd(homepath)
source("~project/microb_metabo/MLVSfunctions.R")
source("~project/microb_metabo/myfunctions.R")
source("~project/food_biomarker/myfunction_2.R")
library(data.table)
library(ggplot2)
library(ggcorrplot)
library(nlme)
# library(egg)
library(patchwork)
library(rlist)
library(gdata)
# library(xlsx)
library(dplyr)
library(purrr)
library(tidyverse)
library(chanmetab)
# library(gee)
library(readxl)
# library(UpSetR)
# library(VennDiagram)
# library(ComplexHeatmap)
library(vegan)
library(caret)
# library(CCA)
# library(CCP)
# library(gdata)
# library(tabula)
options(bitmapType='cairo')
###################################################################################
###################################################################################
lnames1 = load(file="Taxon_use_AVG.RData")

lnames2 = load(file="TAXON_ONLY.RData")

AS_var=c("acesk_","aspart_","sach_","sucral_","tag_","eryth_",
         "maltitol_","mani_",
         # "pini_",
         "sorb_",
         "xyli_")
########################################################################################
count = taxon.avg[,Annot_Taxon$microName]
count[count>=0.0001]=1
count[count<0.0001]=0

abd = data.frame(abundence=colMeans(taxon.avg[,Annot_Taxon$microName]), 
                 detection=colSums(count)/dim(count)[1])
abd$microName = rownames(abd)

abd = merge(abd,Annot_Taxon,by="microName")
head(abd)
dim(abd)

use = taxon.avg[,as.character(Annot_Taxon[which(Annot_Taxon$type=="species"),"microName"])]
# use = as.matrix(use)
head(use)[,1:10]
dim(use) # 923 139
########################
use_ast<-asin(sqrt(use))
use_ast_scale<-as.data.frame(scale(use_ast))
use_scale<-as.data.frame(scale(use))
########################
INT<-function(dta,X)
{
  
  for (i in 1:length(X)){
    
    dta[,paste(X[i],"_INT",sep="")]<-qnorm((rank(dta[,X[i]],na.last="keep")-0.5)/sum(!is.na(dta[,X[i]])))
    
  }
  
  return(dta)
  
}

#Log(X+ minimum non-zero value)
LOG.mini<-function(dta,X)
{
  for (i in 1:length(X)){
    
    mini<-min(dta[dta[X[i]]>0,X[i]])
    dta[,paste(X[i],"_LOG",sep="")]<-log(dta[X[i]]+mini)
    
  }
  
  return(dta)  
  
}

taxon.avg<-INT(taxon.avg,paste(AS_var, "avgall"   ,sep=""))
taxon.avg<-LOG.mini(taxon.avg,paste(AS_var, "avgall"   ,sep=""))
########################################################################
eln.AS<-function(data,
                 a,
                 filt,
                 X,
                 Y){
  
  # set.seed(8888)
  
  sp.unique<-NULL
  co_rr<-NULL
  tune.para<-NULL
  sp.id<-list()
  result<-list()
  
grid <- expand.grid(
    alpha = a,
    lambda = 10^seq(-5,1,length=1000)
  )

X<-X[,abd[which(abd$detection>filt & abd$type=="species"),"microName"]]
print(dim(X))
  
  for(i in 1:length(Y)){
    
    x<-as.data.frame(X)
    y<-data[,Y[i]]
    z<-cbind(x,y)
  
    print(dim(z))
    
    model <- train(
      y ~., 
      data = z, 
      method = "glmnet",
      # preProcess="scale",
      trControl = trainControl("LOOCV"),
      # tuneLength=25
      tuneGrid=grid
    )

    result<-list.append(result,model$result)

    names(result)[length(result)]<-paste(Y[i],sep = "")
    
    elnet<-coef(model$finalModel, model$bestTune$lambda)
    
    tune.para<-rbind(tune.para,model$bestTune)
    
    rownames(tune.para)[i]<-Y[i]
    
    elnet_coef<-as.data.frame(elnet[which(elnet[,1]!=0),1][-1])
    
    coef.out=as.data.frame(t(elnet_coef))
    
    rownames(coef.out)<-"Beta"
    
    CORR<-cor(as.matrix(data[colnames(coef.out)])%*%t(coef.out),data[Y[i]])
    
    co_rr<-cbind(co_rr,CORR)
    
    rownames(co_rr)<-"r"
    colnames(co_rr)[i]<-Y[i]
    
    sp.unique<-c(sp.unique,rownames(elnet_coef))
    
    sp.id<-list.append(sp.id,coef.out)
    
    names(sp.id)[length(sp.id)]<-paste(Y[i],sep = "")
    
  }
  
  out=list(tune.para,co_rr,sp.id,unique(sp.unique),result)
  
  names(out)<-c("tuning_para",
                "cor_score",
                "sp_id",
                "sp_unique",
                "training")
  
  return(out)
  
}
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.1.3  <-eln.AS(taxon.avg,  
                                   0.1,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.1.3,file="INT_alpha_0.1_0.3.RData")
################################################################################################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.2.3  <-eln.AS(taxon.avg,  
                                   0.2,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.2.3,file="INT_alpha_0.2_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.3.3  <-eln.AS(taxon.avg,  
                                   0.3,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.3.3,file="INT_alpha_0.3_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.4.3  <-eln.AS(taxon.avg,  
                                   0.4,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.4.3,file="INT_alpha_0.4_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.5.3  <-eln.AS(taxon.avg,  
                                   0.5,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.5.3,file="INT_alpha_0.5_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.6.3  <-eln.AS(taxon.avg,  
                                   0.6,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.6.3,file="INT_alpha_0.6_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.7.3  <-eln.AS(taxon.avg,  
                                   0.7,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.7.3,file="INT_alpha_0.7_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.8.3  <-eln.AS(taxon.avg,  
                                   0.8,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.8.3,file="INT_alpha_0.8_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT.9.3  <-eln.AS(taxon.avg,  
                                   0.9,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT.9.3,file="INT_alpha_0.9_0.3.RData")
########################################################################
start.time <- Sys.time()

ast_scale_avgall_INT1.3  <-eln.AS(taxon.avg,  
                                   1,
                                   0.3,
                                   use_ast_scale, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(ast_scale_avgall_INT1.3,file="INT_alpha_1_0.3.RData")