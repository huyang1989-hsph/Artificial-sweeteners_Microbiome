library(data.table)
library(ggplot2)
library(ggh4x)
library(ggnewscale)
library(ggforce)
library(ggcorrplot)
library(ggrepel)
library(ggalluvial)
library(ggtext)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(hrbrthemes)
library(nlme)
library(egg)
library(patchwork)
library(rlist)
library(gdata)
library(forcats)
# library(xlsx)
library(dplyr)
library(purrr)
library(tidyverse)
library(gee)
library(readxl)
library(UpSetR)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(caret)
library(glmnet)
library(CCA)
library(CCP)
library(tabula)
library(haven)
library(PMA)
library(openxlsx)
library(grid)
library(patchwork)
library(cowplot)
###################################################################################
setwd("~/taxon/0.01%")
###################################################################################
####Prevalence>0.3
load(file="INT_alpha_0.1_0.3.RData")
load(file="INT_alpha_0.2_0.3.RData")
load(file="INT_alpha_0.3_0.3.RData")
load(file="INT_alpha_0.4_0.3.RData")
load(file="INT_alpha_0.5_0.3.RData")
load(file="INT_alpha_0.6_0.3.RData")
load(file="INT_alpha_0.7_0.3.RData")
load(file="INT_alpha_0.8_0.3.RData")
load(file="INT_alpha_0.9_0.3.RData")
load(file="INT_alpha_1_0.3.RData")

load(file="~/Taxon_use_AVG.RData")
load(file="~/TAXON_ONLY.RData")

load(file="~/MET_pm.RData")
load(file="~/MET_fea.RData")
###################################################################################
AS_var=c("acesk_","aspart_","sach_","sucral_",
         # "tag_",
         "eryth_",
         "maltitol_","mani_",
         # "pini_",
         "sorb_",
         "xyli_")

##################################################################################
acesk_INT.3<-Reduce(intersect,
                    list(
                      names(ast_scale_avgall_INT.1.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.2.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.3.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.4.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.5.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.6.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.7.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.8.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT.9.3$sp_id$acesk_avgall_INT),
                      names(ast_scale_avgall_INT1.3$sp_id$acesk_avgall_INT)
                    ) 
                    ,accumulate = FALSE)

acesk_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.2.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.3.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.4.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.5.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.6.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.7.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.8.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT.9.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         ast_scale_avgall_INT1.3$sp_id$acesk_avgall_INT[,acesk_INT.3])

acesk_INT_beta_final.3<-apply(acesk_INT_beta.3,2,mean)
acesk_INT_beta_final.3<-as.data.frame(acesk_INT_beta_final.3)
acesk_INT_beta_final.3$microName<-rownames(acesk_INT_beta_final.3)
##################################################################################
aspart_INT.3<-Reduce(intersect,
                     list(
                       names(ast_scale_avgall_INT.1.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.2.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.3.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.4.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.5.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.6.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.7.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.8.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT.9.3$sp_id$aspart_avgall_INT),
                       names(ast_scale_avgall_INT1.3$sp_id$aspart_avgall_INT)
                     ) 
                     ,accumulate = FALSE)

aspart_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.2.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.3.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.4.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.5.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.6.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.7.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.8.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT.9.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          ast_scale_avgall_INT1.3$sp_id$aspart_avgall_INT[,aspart_INT.3])

aspart_INT_beta_final.3<-apply(aspart_INT_beta.3,2,mean)
aspart_INT_beta_final.3<-as.data.frame(aspart_INT_beta_final.3)
aspart_INT_beta_final.3$microName<-rownames(aspart_INT_beta_final.3)
##################################################################################
sach_INT.3<-Reduce(intersect,
                   list(
                     names(ast_scale_avgall_INT.1.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.2.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.3.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.4.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.5.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.6.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.7.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.8.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT.9.3$sp_id$sach_avgall_INT),
                     names(ast_scale_avgall_INT1.3$sp_id$sach_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

sach_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.2.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.3.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.4.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.5.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.6.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.7.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.8.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT.9.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        ast_scale_avgall_INT1.3$sp_id$sach_avgall_INT[,sach_INT.3])

sach_INT_beta_final.3<-apply(sach_INT_beta.3,2,mean)
sach_INT_beta_final.3<-as.data.frame(sach_INT_beta_final.3)
sach_INT_beta_final.3$microName<-rownames(sach_INT_beta_final.3)
##################################################################################
sucral_INT.3<-Reduce(intersect,
                     list(
                       names(ast_scale_avgall_INT.1.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.2.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.3.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.4.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.5.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.6.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.7.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.8.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT.9.3$sp_id$sucral_avgall_INT),
                       names(ast_scale_avgall_INT1.3$sp_id$sucral_avgall_INT)
                     ) 
                     ,accumulate = FALSE)

sucral_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.2.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.3.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.4.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.5.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.6.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.7.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.8.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT.9.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          ast_scale_avgall_INT1.3$sp_id$sucral_avgall_INT[,sucral_INT.3])

sucral_INT_beta_final.3<-apply(sucral_INT_beta.3,2,mean)
sucral_INT_beta_final.3<-as.data.frame(sucral_INT_beta_final.3)
sucral_INT_beta_final.3$microName<-rownames(sucral_INT_beta_final.3)
##################################################################################
eryth_INT.3<-Reduce(intersect,
                    list(
                      names(ast_scale_avgall_INT.1.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.2.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.3.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.4.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.5.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.6.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.7.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.8.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT.9.3$sp_id$eryth_avgall_INT),
                      names(ast_scale_avgall_INT1.3$sp_id$eryth_avgall_INT)
                    ) 
                    ,accumulate = FALSE)

eryth_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.2.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.3.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.4.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.5.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.6.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.7.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.8.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT.9.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         ast_scale_avgall_INT1.3$sp_id$eryth_avgall_INT[,eryth_INT.3])

eryth_INT_beta_final.3<-apply(eryth_INT_beta.3,2,mean)
eryth_INT_beta_final.3<-as.data.frame(eryth_INT_beta_final.3)
eryth_INT_beta_final.3$microName<-rownames(eryth_INT_beta_final.3)
##################################################################################
maltitol_INT.3<-Reduce(intersect,
                       list(
                         names(ast_scale_avgall_INT.1.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.2.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.3.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.4.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.5.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.6.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.7.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.8.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT.9.3$sp_id$maltitol_avgall_INT),
                         names(ast_scale_avgall_INT1.3$sp_id$maltitol_avgall_INT)
                       ) 
                       ,accumulate = FALSE)

maltitol_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.2.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.3.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.4.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.5.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.6.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.7.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.8.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT.9.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            ast_scale_avgall_INT1.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3])

maltitol_INT_beta_final.3<-apply(maltitol_INT_beta.3,2,mean)
maltitol_INT_beta_final.3<-as.data.frame(maltitol_INT_beta_final.3)
maltitol_INT_beta_final.3$microName<-rownames(maltitol_INT_beta_final.3)
##################################################################################
mani_INT.3<-Reduce(intersect,
                   list(
                     names(ast_scale_avgall_INT.1.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.2.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.3.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.4.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.5.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.6.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.7.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.8.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT.9.3$sp_id$mani_avgall_INT),
                     names(ast_scale_avgall_INT1.3$sp_id$mani_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

mani_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.2.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.3.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.4.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.5.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.6.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.7.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.8.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT.9.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        ast_scale_avgall_INT1.3$sp_id$mani_avgall_INT[,mani_INT.3])

mani_INT_beta_final.3<-apply(mani_INT_beta.3,2,mean)
mani_INT_beta_final.3<-as.data.frame(mani_INT_beta_final.3)
mani_INT_beta_final.3$microName<-rownames(mani_INT_beta_final.3)
##################################################################################
sorb_INT.3<-Reduce(intersect,
                   list(
                     names(ast_scale_avgall_INT.1.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.2.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.3.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.4.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.5.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.6.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.7.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.8.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT.9.3$sp_id$sorb_avgall_INT),
                     names(ast_scale_avgall_INT1.3$sp_id$sorb_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

sorb_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.2.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.3.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.4.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.5.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.6.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.7.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.8.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT.9.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        ast_scale_avgall_INT1.3$sp_id$sorb_avgall_INT[,sorb_INT.3])

sorb_INT_beta_final.3<-apply(sorb_INT_beta.3,2,mean)
sorb_INT_beta_final.3<-as.data.frame(sorb_INT_beta_final.3)
sorb_INT_beta_final.3$microName<-rownames(sorb_INT_beta_final.3)
##################################################################################
xyli_INT.3<-Reduce(intersect,
                   list(
                     names(ast_scale_avgall_INT.1.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.2.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.3.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.4.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.5.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.6.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.7.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.8.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT.9.3$sp_id$xyli_avgall_INT),
                     names(ast_scale_avgall_INT1.3$sp_id$xyli_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

xyli_INT_beta.3<- rbind(ast_scale_avgall_INT.1.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.2.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.3.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.4.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.5.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.6.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.7.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.8.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT.9.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        ast_scale_avgall_INT1.3$sp_id$xyli_avgall_INT[,xyli_INT.3])

xyli_INT_beta_final.3<-apply(xyli_INT_beta.3,2,mean)
xyli_INT_beta_final.3<-as.data.frame(xyli_INT_beta_final.3)
xyli_INT_beta_final.3$microName<-rownames(xyli_INT_beta_final.3)
####################################################################################
sp.3<-unique(c(acesk_INT_beta_final.3$microName,
               aspart_INT_beta_final.3$microName,
               sach_INT_beta_final.3$microName,
               sucral_INT_beta_final.3$microName,
               eryth_INT_beta_final.3$microName,  
               #tag_INT_beta_final$microName,
               maltitol_INT_beta_final.3$microName,
               mani_INT_beta_final.3$microName,
               sorb_INT_beta_final.3$microName,
               xyli_INT_beta_final.3$microName
)
)
length(sp.3)
###########################
OUT.LONG<-rbind(as.matrix(cbind(TYPE="A",merge(acesk_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="B",merge(aspart_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="C",merge(sach_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="D",merge(sucral_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                
                as.matrix(cbind(TYPE="E",merge(eryth_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="F",merge(maltitol_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="G",merge(mani_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="H",merge(sorb_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                as.matrix(cbind(TYPE="I",merge(xyli_INT_beta_final.3,Annot_Taxon[,c("microName","species")],by="microName",all.x=T))),
                
                as.matrix(cbind(TYPE="J",M.score[which(M.score$M1 !=0),c("microName","M1","species")])),
                as.matrix(cbind(TYPE="K",M.score[which(M.score$M2 !=0),c("microName","M2","species")]))
)

OUT.LONG<-as.data.frame(OUT.LONG)

OUT.LONG$microName<-NULL
OUT.LONG$phylum<-NULL

names(OUT.LONG)[2]<-"Betas"

write.csv(OUT.LONG,"~/AS.long_final.csv",row.names = F)
########################
INT<-function(dta,X)
{
  
  for (i in 1:length(X)){
    
    dta[,paste(X[i],"_INT",sep="")]<-qnorm((rank(dta[,X[i]],na.last="keep")-0.5)/sum(!is.na(dta[,X[i]])))
    
  }
  
  return(dta)
  
}
########################
taxon.avg<-INT(taxon.avg,paste(AS_var, "avgall"   ,sep=""))
########################
acesk.score.3<- as.matrix(taxon.avg[acesk_INT_beta_final.3$microName])  %*% acesk_INT_beta_final.3$acesk_INT_beta_final.3  ;colnames(acesk.score.3)<-"acesk.score.3" 
aspart.score.3<-as.matrix(taxon.avg[aspart_INT_beta_final.3$microName]) %*% aspart_INT_beta_final.3$aspart_INT_beta_final.3 ;colnames(aspart.score.3)<-"aspart.score.3" 
sach.score.3<-  as.matrix(taxon.avg[sach_INT_beta_final.3$microName])   %*% sach_INT_beta_final.3$sach_INT_beta_final.3   ;colnames(sach.score.3)<-"sach.score.3" 
sucrl.score.3<- as.matrix(taxon.avg[sucral_INT_beta_final.3$microName]) %*% sucral_INT_beta_final.3$sucral_INT_beta_final.3 ;colnames(sucrl.score.3)<-"sucrl.score.3"  

eryth.score.3<-as.matrix(taxon.avg[eryth_INT_beta_final.3$microName])    %*% eryth_INT_beta_final.3$eryth_INT_beta_final.3    ;colnames(eryth.score.3)<-"eryth.score.3" 
malt.score.3<- as.matrix(taxon.avg[maltitol_INT_beta_final.3$microName]) %*% maltitol_INT_beta_final.3$maltitol_INT_beta_final.3 ;colnames(malt.score.3)<-"malt.score.3" 
mani.score.3<- as.matrix(taxon.avg[mani_INT_beta_final.3$microName])     %*% mani_INT_beta_final.3$mani_INT_beta_final.3     ;colnames(mani.score.3)<-"mani.score.3" 
sorb.score.3<- as.matrix(taxon.avg[sorb_INT_beta_final.3$microName])     %*% sorb_INT_beta_final.3$sorb_INT_beta_final.3     ;colnames(sorb.score.3)<-"sorb.score.3" 
xyli.score.3<- as.matrix(taxon.avg[xyli_INT_beta_final.3$microName])     %*% xyli_INT_beta_final.3$xyli_INT_beta_final.3     ;colnames(xyli.score.3)<-"xyli.score.3" 

AS.corr.3<-as.data.frame(rbind(cor(acesk.score.3,taxon.avg$acesk_avgall_INT),
                               cor(aspart.score.3,taxon.avg$aspart_avgall_INT),
                               cor(sach.score.3,taxon.avg$sach_avgall_INT),
                               cor(sucrl.score.3,taxon.avg$sucral_avgall_INT),
                               cor(eryth.score.3,taxon.avg$eryth_avgall_INT),
                               cor(malt.score.3,taxon.avg$maltitol_avgall_INT),
                               cor(mani.score.3,taxon.avg$mani_avgall_INT),
                               cor(sorb.score.3,taxon.avg$sorb_avgall_INT),
                               cor(xyli.score.3,taxon.avg$xyli_avgall_INT)
)
)
##########################################################
### Heatmap to show ELN results ###
##########################################################
annot_sp<-Annot_Taxon[which(Annot_Taxon$microName %in% sp.3),]
annot_sp<-annot_sp[,c("microName","phylum","species")]

annot_sp_sig_FDR<-list(annot_sp, 
                       abd[,c("microName","abundance")],
                       acesk_INT_beta_final.3,
                       aspart_INT_beta_final.3,
                       sach_INT_beta_final.3,
                       sucral_INT_beta_final.3,
                       eryth_INT_beta_final.3,
                       maltitol_INT_beta_final.3,
                       mani_INT_beta_final.3,
                       sorb_INT_beta_final.3,
                       xyli_INT_beta_final.3) %>% reduce(left_join, by='microName') 

annot_sp_sig_FDR<-annot_sp_sig_FDR[order(annot_sp_sig_FDR$phylum),]


htmap_coef<-annot_sp_sig_FDR[,c("species",
                                "acesk_INT_beta_final.3",
                                "aspart_INT_beta_final.3",
                                "sach_INT_beta_final.3",
                                "sucral_INT_beta_final.3",  
                                "eryth_INT_beta_final.3", 
                                "maltitol_INT_beta_final.3", 
                                "mani_INT_beta_final.3", 
                                "sorb_INT_beta_final.3", 
                                "xyli_INT_beta_final.3"
)
]

rownames(htmap_coef)<-htmap_coef$species
htmap_coef<-htmap_coef[,-1]
colnames(htmap_coef)<-c("Acesulfame K",
                        "Aspartame",
                        "Saccharin",
                        "Sucralose",
                        "Erythritol",
                        "Maltitol",
                        "Mannitol",
                        "Sorbitol",
                        "Xylitol")

rownames(htmap_coef)<-gsub("^s__", "", rownames(htmap_coef))

write.csv(htmap_coef,"ring.csv",row.names=T)

htmap_coef<-scale(htmap_coef)

row_ha1 <-columnAnnotation("Abundance" = anno_barplot(-log(annot_sp_sig_FDR$abundance),border=FALSE,
                                                      gp = gpar(fill = "#4C0099")
),
annotation_name_side = "left",
annotation_name_rot=0,
annotation_label = NULL)

row_ha2<-columnAnnotation(
  foo = tolower((substr(annot_sp_sig_FDR$phylum, 4 , 10))),
  col = list(foo = c("actinob" = "#F26A0F",
                     "bacteri" = "#000000", 
                     "bactero" = "#00FF80", 
                     "firmicu" = "#7231E2",
                     "proteob" = "#00FFFF",
                     "thaumar" = "#FF007F",
                     "verruco" = "#CCCC00"
  )
  ),
  height = unit(13, "cm"),
  gap = unit(8, "cm"),
  
  # gp = gpar(col = "black"),
  annotation_name_side = "left",annotation_name_rot=45,
  annotation_legend_param = list(foo = list(title = "Phyla",
                                            labels = c("Actinobacteria",
                                                       "Bacteria_unclassified",
                                                       "Bacteroidetes", 
                                                       "Firmicutes",
                                                       "Proteobacteria",
                                                       "Thaumarchaeota",
                                                       "Verrucomicrobia"
                                            )
  )
  ),
  annotation_label = "Phyla")

dff<-AS.corr.3
dff$lbl=NULL
dff$group=NULL
dff<-round(dff,2)

dff<-scale(dff,center=F,scale=F)

a=default_axis_param(t(dff))
a$labels<-AS.corr.3$lbl

a$side<-"top"

.hist = anno_barplot(dff[,1], 
                     width = unit(5, "cm"),
                     which="row",
                     gp = gpar(fill = "#FFFF00"),
                     border=FALSE,
                     add_numbers = TRUE,
                     numbers_gp=gpar(fontface="bold",fontsize = 9),
                     numbers_rot = 0
                     
                     
)

ha_mix_top = rowAnnotation(
  hist = .hist,
  annotation_label="Pearson \nr",
  annotation_name_rot=0,
  annotation_name_side="bottom",
  annotation_name_gp=gpar(fontface="bold",fontsize=8)
)

col_fun = colorRamp2(c(min(htmap_coef,na.rm=T), 0, max(htmap_coef,na.rm=T)), c("darkblue", "aliceblue", "darkred"))

ht<-Heatmap(t(htmap_coef), 
            col=col_fun,
            na_col="#FFFFFF",
            rect_gp = gpar(col = "black", lwd = 1),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 10),  
            show_row_dend = FALSE,
            row_order = order(as.numeric(gsub("row", "", rownames(t(htmap_coef))))),
            row_split = c(rep("A",4),rep("B",5)),
            row_title = NULL,
            row_gap = unit(5, "mm"),
            width = unit(50, "cm"),
            column_names_side = "bottom", 
            column_names_rot = 65,
            column_names_gp = gpar(fontsize = 12),
            column_gap = unit(3, "mm"),
            show_column_dend = FALSE,
            column_order = order(as.numeric(gsub("column", "", colnames(t(htmap_coef))))),
            column_split = rep(LETTERS[1:length(table(annot_sp_sig_FDR$phylum))],table(annot_sp_sig_FDR$phylum)[1:length(table(annot_sp_sig_FDR$phylum))]),
            column_title = NULL,
            # left_annotation = row_ha1,
            right_annotation = ha_mix_top,
            top_annotation = row_ha1,
            bottom_annotation = row_ha2,
            heatmap_legend_param = list(title = "Beta coefficients",
                                        at=c(-3,-2,-1,0,1,2,3),
                                        labels=c(-3,-2, -1,0,1,2,3),
                                        legend_direction="horizontal")
            
            # cell_fun = function(j, i, x, y, w, h, fill) {
            #   if(htmap_adjp[i, j] < 0.05) {
            #     grid.text("**", x, y)
            #   } else if(htmap_pval[i, j] < 0.05) {
            #     grid.text("*", x, y)
            #   }
            # }
)

png(paste("~/fig2_rot.png",sep=''),res=400,width=10000, height=4000)
draw(ht,annotation_legend_side = "top",heatmap_legend_side = "top",
     padding = unit(c(2.4, 0, 1, 0), "cm"))
dev.off()
##########################################################
count = taxon.avg[,Annot_Taxon$microName]
count[count>=0.0001]=1
count[count<0.0001]=0

abd = data.frame(abundance=colMeans(taxon.avg[,Annot_Taxon$microName]), 
                 detection=colSums(count)/dim(count)[1])
abd$microName = rownames(abd)

abd = merge(abd,Annot_Taxon,by="microName")
head(abd)
dim(abd)

use = taxon.avg[,as.character(Annot_Taxon[which(Annot_Taxon$type=="species"),"microName"])]
# use = as.matrix(use)
rownames(use) = taxon.avg$id
head(use)[,1:10]
dim(use) # 305 2201

#Prevalence >0.3
use = use[,abd[which(abd$detection>0.3 & abd$type=="species"),"microName"]]
# use = as.matrix(use)
head(use)[,1:10]
dim(use) # 305 169

d = vegan::vegdist(use,distance="bray",na.rm=TRUE) 
mds = cmdscale(d,k=2,eig=T, add=T) # isoMDS
p=as.data.frame(mds$points)
names(p) = paste("Pco",1:2,sep='')
p$id = rownames(p)
head(p)
dim(p) # 305 3

T_use = merge(taxon.avg,p,by="id")
head(T_use)[,1:10]
dim(T_use) # 305 9213

traiti = paste(AS_var, "avgall_INT",sep="")

covs = c("age_fec","bmi_bld","totMETs_paq","calor_fo_dr_wtavg","smoke_bld", "lt_ahei_dm",
         "probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec",
         "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")

set.seed(2025)
pnv<-c()
for (i in 1:length(traiti)){
  
  modeli = paste("~", paste(c(colnames(use),traiti[i],covs),collapse=" + "))
  dati = model.matrix(as.formula(modeli),data=T_use)
  fml<-paste("dati[,colnames(use)]~",traiti[i],"+ age_fec + bmi_bld + totMETs_paq + calor_fo_dr_wtavg + lt_ahei_dm + probio_2m_fec + antibio_12m_fec + stooltype_fec.1 + stooltype_fec.2 + stooltype_fec.3 + stooltype_fec.4 + stooltype_fec.5 + stooltype_fec.6",sep="")

  xx<-vegan::adonis(as.formula(fml), data=as.data.frame(dati), 
                    permutations = 999, method = "bray", strata=as.data.frame(dati)$id)
  pnv<-rbind(pnv,xx$aov.tab[1,])
}

pnv.plot<-as.data.frame(pnv[,c("R2","Pr(>F)")])
pnv.plot$AS<-c("Acesulfame K", 
               "Aspartame", 
               "Saccharin", 
               "Sucralose",
               "Erythritol",
               "Maltitol",
               "Mannitol",
               "Sorbitol",
               "Xylitol")
names(pnv.plot)[2]<-"P"
pnv.plot$AS <- factor(pnv.plot$AS, levels = pnv.plot$AS)
pnv.plot$significance <- ifelse(pnv.plot$P < 0.05, "p < 0.05", "ns")

png(paste("~/PERMANOVA.png",sep=''),res=400,width=7000, height=3500)
ggplot(pnv.plot, aes(x = AS, y = R2, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("p < 0.05" = "#E74C3C", "ns" = "#95A5A6")) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = NULL,
    y = expression(R^2),
    fill = NULL,
    title = " PERMANOVA test for overall associations"
  ) +
  theme(
    axis.text.y = element_text(size=18, face = "bold", color="black"),
    axis.text.x = element_text(size=18, face = "bold", color="black",hjust = -0.05, vjust = 1),
    axis.title.x = element_text(size=18, face = "bold", color="black"),
    axis.title.y = element_text(size=18, face = "bold", color="black"),
    legend.text = element_text(size=15, face = "bold", color="black"),
    panel.grid.major = element_blank(),
 
    panel.background = element_blank(),
    plot.title = element_text(size=25,face = "bold", hjust = 0.5)
  )
dev.off()
#############################################################
fit <- vegan::envfit(mds, use, permutations = 0)
loadings <- as.data.frame(fit$vectors$arrows)
loadings$microName <- rownames(loadings)

loadings<-merge(loadings,Annot_Taxon[,c("microName","species")],by="microName")
loadings$microName<-NULL
names(loadings)[3]<-"feature"
loadings$feature<-gsub("^s__", "", loadings$feature)

colnames(loadings)[1:2] <- c("Pco1", "Pco2")
loadings$Pco1 <- loadings$Pco1 * fit$vectors$r
loadings$Pco2 <- loadings$Pco2 * fit$vectors$r

loadings$arrow_length <- sqrt(loadings$Pco1^2 + loadings$Pco2^2)
loadings_filtered <- loadings %>%
  arrange(desc(arrow_length)) %>%
  slice(1:10)

dati = T_use[,c(paste(AS_var, "avgall_INT"   ,sep=""),"Pco1","Pco2")]

dati_long <- dati %>%
  pivot_longer(
    cols = ends_with("_INT"),  # or specify: acesk_avgall_INT:sucral_avgall_INT
    names_to = "exposure",
    values_to = "value"
  )

dati_long$exposure <- factor(dati_long$exposure,
                             levels = c("acesk_avgall_INT", 
                                        "aspart_avgall_INT", 
                                        "sach_avgall_INT", 
                                        "sucral_avgall_INT",
                                        "eryth_avgall_INT",
                                        "maltitol_avgall_INT",
                                        "mani_avgall_INT",
                                        "sorb_avgall_INT",
                                        "xyli_avgall_INT"),
                             labels = c("Acesulfame K", 
                                        "Aspartame", 
                                        "Saccharin", 
                                        "Sucralose",
                                        "Erythritol",
                                        "Maltitol",
                                        "Mannitol",
                                        "Sorbitol",
                                        "Xylitol"))
png(paste("~/pco_ss.png",sep=''),res=400,width=7500, height=5000)
ggplot(dati_long %>%
         filter(exposure %in% c("Acesulfame K", 
                                "Aspartame", 
                                "Saccharin", 
                                "Sucralose"))
       , aes(x = Pco1, y = Pco2, color = value)) +
  geom_point(shape = 16, size = 2) + 
  scale_color_gradientn(colors = c("blue", "white", "red"),name = "Normalized intake") + 
  facet_wrap(~ exposure, ncol = 2) +
  # Arrows for features
  geom_segment(data = loadings_filtered, 
               aes(x = 0, y = 0, xend = Pco1, yend = Pco2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  
  # Text labels for feature arrows
  geom_text(data = loadings_filtered, 
            mapping = aes(x = Pco1, y = Pco2, label = feature),
            inherit.aes = FALSE,
            size = 3.5, fontface = "bold.italic" ,hjust = 1, vjust = 1.5) +
  
  xlab("Pco1") +
  ylab("Pco2") +
  
  ggtitle("") +
  
  theme(strip.background = element_rect(fill = "pink", color = "black"),
        strip.text = element_text(size = 16, face = "bold", family = "serif"),
        axis.text.y = element_text(size=11, face = "bold", color="black"),
        axis.text.x = element_text(size=11, face = "bold", color="black",hjust = -0.05, vjust = 1),
        axis.title.x = element_text(size=15, face = "bold", color="black"),
        axis.title.y = element_text(size=15, face = "bold", color="black"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size=15, face = "bold", color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

png(paste("~/pco_sa.png",sep=''),res=400,width=7500, height=5000)
ggplot(dati_long %>%
         filter(exposure %in% c("Erythritol",
                                "Maltitol",
                                "Mannitol",
                                "Sorbitol",
                                "Xylitol"))
       , aes(x = Pco1, y = Pco2, color = value)) +
  geom_point(shape = 16, size = 2) + 
  scale_color_gradientn(colors = c("blue", "white", "red"),name = "Normalized intake") + 
  facet_wrap(~ exposure, ncol = 2) +
  # Arrows for features
  geom_segment(data = loadings_filtered, 
               aes(x = 0, y = 0, xend = Pco1, yend = Pco2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  
  # Text labels for feature arrows
  geom_text(data = loadings_filtered, 
            mapping = aes(x = Pco1, y = Pco2, label = feature),
            inherit.aes = FALSE,
            size = 3.5, fontface = "bold.italic" ,hjust = 1, vjust = 1.5) +
  
  xlab("Pco1") +
  ylab("Pco2") +
  
  ggtitle("") +
  
  theme(strip.background = element_rect(fill = "lightblue", color = "black"),
        strip.text = element_text(size = 16, face = "bold", family = "serif"),
        axis.text.y = element_text(size=11, face = "bold", color="black"),
        axis.text.x = element_text(size=11, face = "bold", color="black",hjust = -0.05, vjust = 1),
        axis.title.x = element_text(size=15, face = "bold", color="black"),
        axis.title.y = element_text(size=15, face = "bold", color="black"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size=15, face = "bold", color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()
########################
use_ast<-asin(sqrt(use))
use_ast_scale<-as.data.frame(scale(use_ast))
use_scale<-as.data.frame(scale(use))

sp.corr<-use_ast[,sp.3]
sp.label<-Annot_Taxon[,c("microName","species")]

sp.corr <- sp.corr %>%
  rename_with(~ sp.label$species[match(.x, sp.label$microName)])

corr_mat <- cor(
  sp.corr,
  use = "pairwise.complete.obs",
  method = "spearman"   # use "pearson" if appropriate
)
corr_long <- reshape2::melt(corr_mat)
colnames(corr_long) <- c("Species1", "Species2", "Correlation")

hc <- hclust(as.dist(1 - corr_mat))

order <- hc$labels[hc$order]

corr_long$Species1 <- factor(corr_long$Species1, levels = order)
corr_long$Species2 <- factor(corr_long$Species2, levels = order)

# corr_long <- corr_long[
#   as.numeric(corr_long$Species1) <= as.numeric(corr_long$Species2), ]

png(paste("~/SP_CORR.png",sep=''),res=400,width=5000, height=5000)
ggplot(corr_long, aes(x = Species1, y = Species2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nCorrelation"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5, face="bold"),
    axis.text.y = element_text(size = 5, face="bold"),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
dev.off()
corr_vals <- corr_mat[upper.tri(corr_mat)]
pos_corr <- corr_vals[corr_vals > 0]
neg_corr <- corr_vals[corr_vals < 0]

# Positive correlations
pos_med <- median(pos_corr, na.rm = TRUE)
pos_q <- quantile(pos_corr, probs = c(0.25, 0.75), na.rm = TRUE)

# Negative correlations
neg_med <- median(neg_corr, na.rm = TRUE)
neg_q <- quantile(neg_corr, probs = c(0.25, 0.75), na.rm = TRUE)

pos_fmt <- sprintf(
  "%.2f (%.2f, %.2f)",
  pos_med, pos_q[1], pos_q[2]
)

neg_fmt <- sprintf(
  "%.2f (%.2f, %.2f)",
  neg_med, neg_q[1], neg_q[2]
)

pos_fmt
# [1] "0.09 (0.04, 0.15)"
neg_fmt
# [1] "-0.06 (-0.12, -0.03)"
########################
##    sCCA   ##
########################

X<-use_ast_scale[,sp.3]
Y<-taxon.avg[,paste(AS_var, "avgall_INT",sep="")]

#################Parameter tuning grid search################# 
start.time <- Sys.time()

xs<-seq(0.2,0.6,length=21)
ys<-seq(0.4,0.8,length=21)
x1y1=matrix(nrow = length(xs), ncol = length(ys), dimnames=list(xs,ys))
x2y2=matrix(nrow = length(xs), ncol = length(ys), dimnames=list(xs,ys))
y12=matrix(nrow = length(xs), ncol = length(ys), dimnames=list(xs,ys))
x12=matrix(nrow = length(xs), ncol = length(ys), dimnames=list(xs,ys))

num_samples <- nrow(X)

score.x<-matrix(nrow = num_samples, ncol = 2)
score.y<-matrix(nrow = num_samples, ncol = 2)

for (i in 1:length(xs))
{
  for (j in 1:length(ys))
  {
    for (k in 1:num_samples)
    {
      out <- CCA(X[-k,],Y[-k,],typex="standard",typez="standard",K=2,
                 penaltyx=xs[i],penaltyz=ys[j])
      
      score.x[k,] <- t(scale(t(X[k,]))) %*% out$u
      score.y[k,] <- t(scale(t(Y[k,]))) %*% out$v
      
    }
    
    print(c(i,j))
    
    x1y1[i,j]<-cor(score.x[,1],score.y[,1])
    x2y2[i,j]<-cor(score.x[,2],score.y[,2])
    
    x12[i,j]<-cor(score.x[,1],score.x[,2])
    y12[i,j]<-cor(score.y[,1],score.y[,2])
  }
}

end.time <- Sys.time()
end.time - start.time

round(x1y1,2)
#############################################################
out <- CCA(X,Y,typex="standard",typez="standard",K=9,
           penaltyx=0.39,penaltyz=0.78)

score.x <- scale(X) %*% out$u
score.y <- scale(Y) %*% out$v
diag(cor(score.x, score.y))[1:2]
# [1] 0.5001866 0.5047993

cor(scale(Y),score.x)[,1:2]
# [,1]        [,2]
# acesk_avgall_INT    -0.09791970  0.37374819
# aspart_avgall_INT   -0.30098724  0.33309440
# sach_avgall_INT     -0.13886946  0.35207647
# sucral_avgall_INT   -0.13219196  0.22211190
# eryth_avgall_INT     0.31776893 -0.09533604
# maltitol_avgall_INT -0.04153871  0.02923521
# mani_avgall_INT      0.30201257 -0.12149137
# sorb_avgall_INT      0.22828845 -0.18807658
# xyli_avgall_INT      0.28327662 -0.24809416

# taxon.avg[,c("Y1","Y2")]<-score.y[,1:2]

corr_pval <- c()
corr_r <- c()
for(j in 1:9){
  corr <- cor.test(score.x[,j],score.y[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}
corr_pval
# [1] 1.030311e-20 3.990587e-21 4.289869e-19 1.214993e-10 1.087584e-13 1.354114e-12 1.118828e-12 1.496084e-11
# [9] 6.412014e-18
p.adjust(corr_pval, method = "BH")
# [1] 4.636401e-20 3.591528e-20 1.286961e-18 1.214993e-10 1.957651e-13 1.741004e-12 1.678241e-12 1.683094e-11
# [9] 1.442703e-17

M.score<-as.data.frame(out$u)
names(M.score)<-paste("M",seq(1,9,1),sep="")
M.score$microName<-colnames(X)
rownames(M.score)<-colnames(X)
M.score<-M.score[which(M.score$M1 !=0 | M.score$M2!=0),]
M.score<-merge(M.score,Annot_Taxon[,c("microName","species")],by="microName")
dim(M.score)
# [1] 55 11
##################################################################################
### Coefficients for each species by 2 variates ###
##################################################################################
scca.X<-as.data.frame(out$u)
scca.X$microName<-colnames(X)
scca.X<-scca.X[which(scca.X$V1 !=0 | scca.X$V2!=0),]

#Plots to show betas for each species
scca.X.p<-scca.X
names(scca.X.p)[1:2]<-c("M1","M2")
scca.X.p<-merge(scca.X.p,Annot_Taxon[,c("microName","phylum","species")],by="microName")
scca.X.p<-scca.X.p[order(scca.X.p$M2,decreasing = TRUE),]
scca.X.p<-scca.X.p[order(scca.X.p$phylum),]
type.label<-unique(sub("p__","",scca.X.p$phylum))

scca.X.p<-scca.X.p[,c("species","M1","M2")]
scca.X.p$species<-gsub("^s__", "", scca.X.p$species)
scca.X.p$species<-factor(scca.X.p$species,levels=unique(scca.X.p$species))
scca.X.p<-reshape2::melt(scca.X.p,id="species")

png(paste("~/scca_coef.png",sep=''),res=400,width=5000, height=7000)
ggplot(scca.X.p, aes(x=species, y=value)) +
  geom_segment(aes(y=0, yend=value, color=variable), size=0.5, alpha=0.9) +
  geom_point(data = subset(scca.X.p, value != 0), aes(color = variable), size = 3, alpha=0.7)+ 
  scale_color_manual(name="Species patterns",labels = c("Sugar alochols", "Synthetic sweeteners"),values = c("red",  "blue"))+
  geom_hline(yintercept = 0)+
  xlab("")+
  ylab("Beta coefficients from sCCA")+
  # scale_color_discrete(name="sCCA variates")+
  theme_classic()+
  scale_x_discrete(limits=rev)+
  theme(legend.text = element_text(size=18,face="bold",color = "black"),
        legend.title = element_text(size=18,face="bold",color = "black"),
        axis.ticks.y = element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(size=12,face="bold",color = "black"),
        axis.title=element_text(size=12,face="bold",color = "black"),
        axis.text.y=element_text(size=12,face="bold",color = "black"))+
  
  geom_rect(aes(xmin=50,xmax=55,ymin=-0.5,ymax=-0.46,fill=type.label[1]))+
  geom_rect(aes(xmin=0.5,xmax=50,ymin=-0.5,ymax=-0.46,fill=type.label[2]))+
  scale_fill_discrete(name="Phyla")+
  
  coord_flip()
dev.off()
##################################################################################
### CCA representation plot ###
##################################################################################
corr.X.xscores<-cor(scale(X),score.x)
corr.Y.xscores<-cor(scale(Y),score.x)
Ynames<-rownames(corr.Y.xscores)

scca.plot<-as.data.frame(cbind(Ynames,corr.Y.xscores))
names(scca.plot)<-c("AS","Variate1", "Variate2")
scca.plot$Variate1<-as.numeric(scca.plot$Variate1)
scca.plot$Variate2<-as.numeric(scca.plot$Variate2)
scca.plot$type<-c(rep("Synthetic artificial sweeteners",4),
                  rep("Sugar alcohols",5))
scca.plot$type<-factor(scca.plot$type,levels=c("Synthetic artificial sweeteners","Sugar alcohols"))
scca.plot$AS<-factor(scca.plot$AS,levels=scca.plot$AS)

write.csv(scca.plot,"scca.plot.csv",row.names = T)

scca.p<-
  ggplot(scca.plot, aes(x=Variate1, y=Variate2)) +
  geom_point(size=12, stroke=2, aes(color=type, shape=AS)) +  # Bold markers
  ylim(-0.5, 0.5) + xlim(-0.5, 0.5) + 
  expand_limits(x=c(-1.1, 1.1), y=c(-1.1, 1.1)) +
  xlab("Dimension 1") + ylab("Dimension 2") +
  
  scale_shape_manual(name="Artificial Sweeteners",
                     labels=c("Acesulfame K", "Aspartame", "Saccharin", 
                              "Sucralose", "Erythritol", "Maltitol", 
                              "Mannitol", "Sorbitol", "Xylitol"),
                     values=seq(4,12)) +  # All shapes are unique
  scale_color_manual(name="Types",
                     labels=c("Synthetic sweeteners", "Sugar alcohols"),
                     values=c("#D55E00", "#0072B2")) +
  
  geom_hline(yintercept=0, color="gray50", linetype="dashed") +
  geom_vline(xintercept=0, color="gray50", linetype="dashed") +
  # geom_circle(aes(x0=0, y0=0, r=1), col="black", size=0.8, inherit.aes=T) +
  geom_circle(aes(x0=0, y0=0, r=0.5), col="black", size=0.8, inherit.aes=T) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.title.x=element_text(size=15, face="bold", color="black"),
    axis.text.x=element_text(size=13, face="bold", color="black"),
    axis.title.y=element_text(size=15, face="bold", color="black"),
    axis.text.y=element_text(size=13, face="bold", color="black"),
    legend.title=element_text(size=15, face="bold", color="black"),
    legend.text=element_text(size=13, face="bold", color="black"),
    legend.position="right",  # Legend on the right
    panel.grid.major=element_line(size=0.5, linetype='dotted', color='grey'),
    legend.key.height = unit(3, "lines")
  )

png(paste("~/scca_repre.png",sep=''),res=400,width=4000, height=4000)
scca.p
dev.off()

ggsave(
  filename = "~/scca_repre.svg",
  plot = scca.p,
  device = "svg",
  width = 10,
  height = 10,
  units = "in"
)
##################################################################################
### Correlation between 2 Y variates and 2 Y ###
##################################################################################
out.y.corr<-as.data.frame(cbind(t(cor(score.y[,1],Y)),
                                t(cor(score.y[,2],Y)))
)

names(out.y.corr)<-c("Sweetener pattern 1","Sweetener pattern 2")
out.y.corr$name<-c("Acesulfame K",
                   "Aspartame",
                   "Saccharin",
                   "Sucralose",
                   "Erythritol",
                   "Maltitol",
                   "Mannitol",
                   "Sorbitol",
                   "Xylitol")
out.y.corr$name<-factor(out.y.corr$name,levels=out.y.corr$name)

write.csv(out.y.corr,"scca.cor.Y.p.csv",row.names = T)

data_melted <- reshape2::melt(out.y.corr)

scca.cor.Y.p<-ggplot(data_melted, aes(x=variable, y=name, fill=value)) + 
  geom_tile(color="white", size=0.5) +  # Adds tile borders
  xlab("") + ylab("") +  # Add meaningful axis labels
  guides(fill = guide_colorbar(title = "Pearson Correlation")) +  # Colorbar legend
  scale_fill_gradient2(
    low = "#4575b4",  # Blue for low values
    mid = "#f7f7f7",  # White/gray for midpoint
    high = "#d73027",  # Red for high values
    midpoint = 0, limits=c(-1, 1),
    breaks=c(-1, -0.5, 0, 0.5, 1),  # 5 breaks for consistent intervals
    labels=c("-1", "-0.5", "0", "0.5", "1")  # Matching labels
  ) +
  geom_text(aes(label=round(value, 2)), color="black", size=6) +  # Add annotations
  theme_minimal() +  # Use a clean theme
  theme(
    legend.position = "bottom",
    legend.text = element_text(face="bold"),
    legend.title = element_text(face="bold"),
    axis.text.x = element_text(face="bold", size=14, color="black"),  # Rotated x-axis
    axis.text.y = element_text(face="bold", size=14, color="black"),  # Bold y-axis
    panel.grid = element_blank(),  # Remove gridlines
    panel.border = element_blank()  # Remove panel borders
  )
png(paste("~/scca.cor.Y.png",sep=''),res=400,width=2500, height=3500)
scca.cor.Y.p
dev.off()
##################################################################################
### Correlation between 2 X variates and 2 X ###
##################################################################################
out.x.corr<-as.data.frame(cbind(t(cor(score.x[,1],X)),t(cor(score.x[,2],X))))
names(out.x.corr)<-c("Species pattern 1","Species pattern 2")
out.x.corr$microName<-rownames(out.x.corr)
out.x.corr<-merge(out.x.corr,Annot_Taxon[,c("microName","phylum","species")],by="microName")
names(out.x.corr)[4:5]<-c("type","name")
out.x.corr$microName<-NULL

out.x.corr$name<-gsub("^s__", "", out.x.corr$name)

out.x.corr<-out.x.corr[order(out.x.corr$type),]
out.x.corr$type<-factor(out.x.corr$type,levels=unique(out.x.corr$type))
out.x.corr$name<-factor(out.x.corr$name,levels=unique(out.x.corr$name))



write.csv(out.x.corr,"scca.cor.X.p.csv",row.names = T)

out.x.corr.long<-reshape2::melt(out.x.corr,id=c("type","name"))
out.x.corr.long<-out.x.corr.long[order(out.x.corr.long$type),]

out.x.corr$type<-sub("p__", "", out.x.corr$type)
table(out.x.corr$type)
# Actinobacteria Bacteria_unclassified         Bacteroidetes            Firmicutes        Proteobacteria        Thaumarchaeota       Verrucomicrobia 
# 7                     1                    23                    92                     5                     1                     1

bin.col<-brewer.pal(n = 7, name = "Set2")
#Locate position to draw rect

type.label<-str_to_title(unique(out.x.corr$type))
scca.cor.X.p<-ggplot(out.x.corr.long, aes(x=variable, y=name, fill= value)) + 
  geom_tile() +
  xlab("")+ylab("")+
  guides(fill = guide_legend(title = "Pearson correlations")) +
  scale_fill_gradient2(low = "blue",
                       mid = "grey",    
                       high = "red", 
                       midpoint = 0,
                       limits=c(-1,1),
                       breaks=c(-0.6,-0.3,0,0.3,0.6)) +
  scale_y_discrete(limits=rev)+
  new_scale("fill")+
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=123.5,ymax=130.5,fill=type.label[1])) +
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=122.5,ymax=123.5,fill=type.label[2])) +
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=100,ymax=122.5,fill=type.label[3])) +
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=8,ymax=100,fill=type.label[4])) +
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=2.1,ymax=8,fill=type.label[5])) +  
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=1.3,ymax=2.1,fill=type.label[6])) +  
  geom_rect(aes(xmin=0.2,xmax=0.48,ymin=0.5,ymax=1.3,fill=type.label[7])) + 
  scale_fill_discrete(name="Phyla")+
  theme(
    axis.text.y = element_text(colour = "black", face = "bold", size = 6),
    axis.text.x = element_text(colour = "black", face = "bold", size = 10),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face="bold"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

png(paste("~/scca.cor.X.png",sep=''),res=400,width=3000, height=4000)
scca.cor.X.p
dev.off()
##################################################################################
### Scatter plots to show correlation between AS intake and species score ###
##################################################################################
M1.scca.score<-score.x[,1]
M2.scca.score<-score.x[,2]

A1.scca.score<-score.y[,1]
A2.scca.score<-score.y[,2]

cc.score = as.data.frame(cbind(M1.scca.score,M2.scca.score,
                               A1.scca.score,A2.scca.score))  
names(cc.score)=c("M1.scca.score","M2.scca.score",  
                  "A1.scca.score","A2.scca.score")
new.dat<-cbind(taxon.avg,
               acesk.score.3 ,
               aspart.score.3,
               sach.score.3  ,
               sucrl.score.3 ,
               eryth.score.3 ,
               malt.score.3  ,
               mani.score.3  ,
               sorb.score.3  ,
               xyli.score.3  ,
               cc.score)

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(acesk.score.3), !is.na(acesk_avgall_INT)) %>%
  mutate(is_min = acesk_avgall_INT == min(acesk_avgall_INT))

acesk<-ggplot(new.dat2,aes(x=acesk.score.3,y=acesk_avgall_INT))+

  # Main sphere (conditional color)
  geom_point(aes(fill = is_min,
                 color = is_min),
             size = 6, shape = 21,
             stroke = 1.5) +
  
  scale_fill_manual(
    values = c(`TRUE` = "grey60",
               `FALSE` = "#ff3b3b"),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(`TRUE` = "grey60",   # grey border
               `FALSE` = "darkred"),  # red border
    guide = "none"
  ) +
  
  geom_smooth(
    data = new.dat2 %>% filter(!is_min),
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Acesulfame K")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-1,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(aspart.score.3), !is.na(aspart_avgall_INT)) %>%
  mutate(is_min = aspart_avgall_INT == min(aspart_avgall_INT))

aspart<-ggplot(new.dat2,aes(x=aspart.score.3,y=aspart_avgall_INT))+
  # Main sphere (conditional color)
  geom_point(aes(fill = is_min,
                 color = is_min),
             size = 6, shape = 21,
             stroke = 1.5) +
  
  scale_fill_manual(
    values = c(`TRUE` = "grey60",
               `FALSE` = "#ff3b3b"),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(`TRUE` = "grey60",   # grey border
               `FALSE` = "darkred"),  # red border
    guide = "none"
  ) +
  
  geom_smooth(
    data = new.dat2 %>% filter(!is_min),
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Aspartame")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-1,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(sach.score.3), !is.na(sach_avgall_INT)) %>%
  mutate(is_min = sach_avgall_INT == min(sach_avgall_INT))

sach<-ggplot(new.dat2,aes(x=sach.score.3,y=sach_avgall_INT))+
  # Main sphere (conditional color)
  geom_point(aes(fill = is_min,
                 color = is_min),
             size = 6, shape = 21,
             stroke = 1.5) +
  
  scale_fill_manual(
    values = c(`TRUE` = "grey60",
               `FALSE` = "#ff3b3b"),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(`TRUE` = "grey60",   # grey border
               `FALSE` = "darkred"),  # red border
    guide = "none"
  ) +
  
  geom_smooth(
    data = new.dat2 %>% filter(!is_min),
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Saccharin")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-1,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(sucrl.score.3), !is.na(sucral_avgall_INT)) %>%
  mutate(is_min = sucral_avgall_INT == min(sucral_avgall_INT))

sucrl<-ggplot(new.dat2,aes(x=sucrl.score.3,y=sucral_avgall_INT))+
  # Main sphere (conditional color)
  geom_point(aes(fill = is_min,
                 color = is_min),
             size = 6, shape = 21,
             stroke = 1.5) +
  
  scale_fill_manual(
    values = c(`TRUE` = "grey60",
               `FALSE` = "#ff3b3b"),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(`TRUE` = "grey60",   # grey border
               `FALSE` = "darkred"),  # red border
    guide = "none"
  ) +
  
  geom_smooth(
    data = new.dat2 %>% filter(!is_min),
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Sucralose")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-1,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(eryth.score.3), !is.na(eryth_avgall_INT)) %>%
  mutate(is_min = eryth_avgall_INT == min(eryth_avgall_INT))

eryth<-ggplot(new.dat2,aes(x=eryth.score.3,y=eryth_avgall_INT))+
  # Main sphere (conditional color)
  geom_point(aes(fill = is_min,
                 color = is_min),
             size = 6, shape = 21,
             stroke = 1.5) +
  
  scale_fill_manual(
    values = c(`TRUE` = "grey60",
               `FALSE` = "#ff3b3b"),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(`TRUE` = "grey60",   # grey border
               `FALSE` = "darkred"),  # red border
    guide = "none"
  ) +
  
  geom_smooth(
    data = new.dat2 %>% filter(!is_min),
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Erythritol")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-3,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(malt.score.3), !is.na(maltitol_avgall_INT)) %>%
  mutate(is_min = maltitol_avgall_INT == min(maltitol_avgall_INT))

malti<-ggplot(new.dat2,aes(x=malt.score.3,y=maltitol_avgall_INT))+
  # Main sphere (conditional color)
  geom_point(aes(fill = is_min,
                 color = is_min),
             size = 6, shape = 21,
             stroke = 1.5) +
  
  scale_fill_manual(
    values = c(`TRUE` = "grey60",
               `FALSE` = "#ff3b3b"),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(`TRUE` = "grey60",   # grey border
               `FALSE` = "darkred"),  # red border
    guide = "none"
  ) +
  
  geom_smooth(
    data = new.dat2 %>% filter(!is_min),
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Maltitol")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-3,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

new.dat2 <- new.dat %>%
  dplyr::filter(!is.na(mani.score.3), !is.na(mani_avgall_INT)) %>%
  mutate(is_min = mani_avgall_INT == min(mani_avgall_INT))

mani<-ggplot(new.dat,aes(x=mani.score.3,y=mani_avgall_INT))+
  
  geom_point(size = 6, shape = 21, fill = "#ff3b3b", color = "darkred", stroke = 1.5)+
  
  geom_smooth(
    data = new.dat,
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Mannitol")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-3,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

sorb<-ggplot(new.dat2,aes(x=sorb.score.3,y=sorb_avgall_INT))+
  
  geom_point(size = 6, shape = 21, fill = "#ff3b3b", color = "darkred", stroke = 1.5)+
  
  geom_smooth(
    data = new.dat,
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  ggtitle("Sorbitol")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-3,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

xyli<-ggplot(new.dat2,aes(x=xyli.score.3,y=xyli_avgall_INT))+
  geom_point(size = 6, shape = 21, fill = "#ff3b3b", color = "darkred", stroke = 1.5)+
  
  geom_smooth(
    data = new.dat,
    method = "lm", se = FALSE, color = "blue", linewidth = 2
  ) +
  
  ggtitle("Xylitol")+xlab("Species score")+ylab("Normalized intake")+
  ylim(-3,3)+
  theme_classic()+
  theme(plot.title = element_text(size=20,face="bold",color = "black",hjust = 0.5),
        axis.title.x =  element_text(size=18,face="bold",color = "black"),
        axis.title.y =  element_text(size=18,face="bold",color = "black"),
        axis.text.x=element_text(size=13,face="bold",color = "black"),
        axis.text.y=element_text(size=13,face="bold",color = "black"))

bin.col<-brewer.pal(n = 7, name = "Set2")

png("~/scatter.all.png",res=350,width=6500, height=6000)
ggpubr::ggarrange(ggarrange(acesk,aspart,sach,sucrl,ncol=4),
                  ggarrange(eryth,malti,mani,sorb,xyli,ncol=5),
                  nrow=2)+
  theme(plot.margin = margin(0.1,0.1,0.1,1, "cm"))+
  geom_rect(aes(xmin=-0.1,xmax=0.003,ymin=0,ymax=0.52,fill=bin.col[1],alpha=0.3))+
  annotate('text',
           x=-0.015,
           y=0.25,
           label='Sugar alcohols',
           color='black',
           family='serif',
           size=8,
           fontface="bold",
           angle=90)+
  geom_rect(aes(xmin=-0.1,xmax=0.003,ymin=0.525,ymax=1,fill=bin.col[1],alpha=0.3))+
  annotate('text',
           x=-0.015,
           y=0.75,
           label='Synthetic sweeteners',
           color='black',
           family='serif',
           size=8,
           fontface="bold",
           angle=90)
dev.off()
###########################################################################
### Metabolites analysis ###
###########################################################################
new.dat.met<-merge(new.dat,lvs_hpfs_pm,by="id")

met_exclude<-names(which(apply(is.na(lvs_hpfs_pm), 2, sum)>100))[-1]#Remove missing >100

met<-rownames(lvs_hpfs_fea)

met<-setdiff(met,met_exclude)

met_name<-lvs_hpfs_fea[which(!(rownames(lvs_hpfs_fea) %in% met_exclude)),"metabolite_name"]

AS_inv<-paste(AS_var, "avgall_INT",sep="")
score<-c(AS_inv,names(new.dat.met)[grep(c("score"),names(new.dat.met))])

#Scale to ensure same unit for beta
new.dat.met[,score]<-scale(new.dat.met[,score])

AdjVars = c("age_fec","bmi_bld","totMETs_paq","lt_ahei_dm","calor_fo_dr_wtavg",
            "probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec",
            "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")
#################################################################################
regress<-function(dta,
                  Y,
                  Y_name,
                  X,
                  X_name,
                  cov,
                  fdr,
                  out_fdr=T){
  k<-list()
  k_fdr<-list()
  
  for (i in 1:length(X)){
    
    coef_out=NULL
    
    for (j in 1:length(Y)){
      
      print(j)
      
      model=as.formula(paste(paste(Y[j],"~",sep=""),paste(c(X[i],cov),sep="",collapse="+"),sep=""))
      
      print(model)
      
      output=lm(model, data = dta)
      
      output_summary=summary(output)
      
      coef=data.frame(Estimate=output_summary$coefficients[X[i],"Estimate"],
                      P=output_summary$coefficients[X[i],"Pr(>|t|)"],
                      met.id=Y[j],
                      nmiss=sum(is.na(dta[,Y[j]])))
      rownames(coef)=Y_name[j]
      
      coef_out=rbind(coef_out,coef)
      
    }
    coef_out$P.adj=p.adjust(coef_out$P, method = "BH")
    k<-list.append(k,coef_out)
    names(k)[length(k)]<-X_name[i]   
    
    coef_out_fdr=coef_out[which(coef_out$P.adj<fdr),]
    
    if (dim(coef_out_fdr)[1]>0)
    {
      k_fdr<-list.append(k_fdr,coef_out_fdr)
      names(k_fdr)[length(k_fdr)]<-X_name[i]   
    }
    
  }
  
  if (out_fdr){return(k)}
  else {return(k_fdr)}
  
  
}

metab<-regress(new.dat.met,
                met,
                met_name,
                score,
                score,
                AdjVars,
                0.05,
                F)

write.csv(metab$M1.scca.score,"M1_met.csv")
write.csv(metab$M2.scca.score,"M2_met.csv")
##################################################################################
### Regress each of 55 species with each MET ###
AS_SP<-regress(new.dat.met,
               unique(c(metab$M1.scca.score$met.id,
                        metab$M2.scca.score$met.id)),
               unique(c(rownames(metab$M1.scca.score),
                        rownames(metab$M2.scca.score))),
               paste("scale(","asin(","sqrt(",M.score$microName,")",")",")",sep=""),
               M.score$microName,
               AdjVars,
               0.05,
               T)

est_list <- lapply(AS_SP, function(df) df[, "Estimate", drop = FALSE])
fdr_list <- lapply(AS_SP, function(df) df[, "P.adj", drop = FALSE])

est_cbind=do.call(cbind, est_list)
names(est_cbind)<-Annot_Taxon[which(Annot_Taxon$microName %in% names(AS_SP)),"species"]
names(est_cbind)<-sub("^s__", "", names(est_cbind))
rownames(est_cbind)[8]<-"Ketoleucine"
write.csv(est_cbind,"est_cbind.csv",row.names = T)

fdr_cbind=do.call(cbind, fdr_list)
names(fdr_cbind)<-Annot_Taxon[which(Annot_Taxon$microName %in% names(AS_SP)),"species"]
names(fdr_cbind)<-sub("^s__", "", names(fdr_cbind))
write.csv(fdr_cbind,"fdr_cbind.csv",row.names = T)

data_matrix <- as.matrix(est_cbind)

ht<-Heatmap(data_matrix, 
            name = "values", 
            cluster_rows = F, 
            cluster_columns = F, 
            show_row_names = F, 
            show_column_names = TRUE,
            column_names_gp = gpar(fontface = "bold"),
            row_names_side = "right",
            column_names_rot = 45,
            width = unit(55.5, "cm") ,
            height =unit(41.2, "cm") ,
            heatmap_legend_param = list(title = "Beta",
                                        at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),
                                        labels=c(-0.3,-0.2, -0.1,0,0.1,0.2,0.3),
                                        title_gp = gpar(fontsize = 25),  # Increase title fontsize
                                        labels_gp = gpar(fontsize = 20),  # Increase labels fontsize
                                        title_position="topcenter",
                                        legend_direction="horizontal"),
            
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(fdr_cbind[i, j] < 0.05) {
                grid.text("**", x, y)
              } 
            }
)

ht_draw<-draw(ht, heatmap_legend_side = "top")

png(paste("~/ht.png",sep=''),
    res=400,width=10000, height=8100)
print(ht_draw)
legend_text <- "FDR < 0.05"
legend_grob <- legendGrob(
  labels = paste0("** ", legend_text),
  gp = gpar(fontsize = 25)
)

# Draw the custom legend
pushViewport(viewport(x = unit(0.35, "npc"), y = unit(1.47, "npc"), 
                      just = c("center", "top")))
grid.draw(legend_grob)
popViewport()
dev.off()

#####################################################
ind_hr=read.csv("indv_met_hr.csv",header=T)

#Only keep HR
ind_hr$Per.SD<-as.numeric(trim(gsub("\\(.*", "", ind_hr$Per.SD)))
ind_hr<-ind_hr[-c(1,2,20,21),]
ind_hr<-rbind(ind_hr[1:2,],
              c("T2D","alpha-tocopherol",NA),
              ind_hr[3:23,],
              c("CHD","alpha-tocopherol",NA),
              ind_hr[24:34,])

ind_hr$Outcome<-factor(ind_hr$Outcome,levels=c("T2D","CHD"))
ind_hr$Metabolites<-factor(ind_hr$Metabolites,levels=unique(ind_hr$Metabolites))
ind_hr$Per.SD<-as.numeric(ind_hr$Per.SD)


ind_hr_alluv<-ggplot(ind_hr, aes(axis1 = Metabolites, axis2 = Outcome)) +
  geom_alluvium(aes(fill = Per.SD), alpha = 0.8) +  # Adjust transparency if needed
  geom_stratum(width = 1/2, fill = c(rep("white",18),"lightblue","lightblue"), color = "gray") +  # Control width here
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 15) +
  scale_x_discrete(limits = c("Metabolites", "Cardiometabolic disease"), expand = c(.05, .05)) +
  
  # Applying a color gradient for continuous estimate values
  scale_fill_gradient(low = "green", high = "red", name = "Hazard Ratios",na.value = NA) +
  
  theme_minimal() + 
  guides(
    fill = guide_colourbar(title.position = "top")  # Applies to continuous color scales
  )+
  
  theme(
    axis.text.x = element_text(size = 50),
    axis.title.y = element_blank(),   # remove y-axis title
    axis.text.y = element_blank(),    # remove y-axis text
    axis.ticks.y = element_blank(),   # remove y-axis ticks
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),   # remove minor grid lines
    legend.text = element_text(size = 20),         # Increase legend text size
    legend.title = element_text(size = 23),        # Increase legend title size
    legend.position = c(0.5,0.98),        # Move legend to the top
    legend.direction = "horizontal",
    axis.text = element_text(colour = "black", size = 30,face = "bold"), 
    legend.justification = "center" # Center the legend
  )+
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.vjust = 0.5, 
      barwidth = 20,        # Increase the width of the color bar
      barheight = 1         # Optionally, adjust the height to keep it proportional
    )
  )

png(paste("~/ht_alluv.png",sep=''),
    res=400,width=10000, height=8000)
ind_hr_alluv
dev.off()
# 
# 
# # Convert heatmap to a grid object
# ht_grob <- grid.grabExpr(draw(ht, heatmap_legend_side = "top"))
# 
# # Create a dummy plot with sufficient height to manipulate layout
# dummy_plot <- ggplot() + theme_void()
# 
# ind_hr_alluv_adjusted <- ind_hr_alluv + 
#   theme(plot.margin = unit(c(1, 1, 2, 5), "lines"))  # Top, right, bottom, left
# 
# 
# # Align plots, here using axis = "tblr" ensures all four margins are considered for alignment
# combined_plot <- plot_grid(
#   plot_grid(dummy_plot, ht_grob, nrow = 2, rel_heights = c(0.1,0.36), align = "v", axis = "tb"),
#   ind_hr_alluv_adjusted,       # ggplot object
#   ncol = 2,           # Number of columns
#   align = "h",        # Aligns plots vertically along top
#   axis = "r"         # Align based on the top and bottom margins
# )
# 
# # Display the combined plot
# print(combined_plot)
# 
# legend_text <- "FDR < 0.05"
# legend_grob <- legendGrob(
#   labels = paste0("* ", legend_text),
#   gp = gpar(fontsize = 12)
# )
# 
# # Draw the custom legend
# pushViewport(viewport(x = unit(0.35, "npc"), y = unit(1.47, "npc"), 
#                       just = c("center", "top")))
# grid.draw(legend_grob)
# popViewport()
# 
# png(paste("~/ht_alluv.png",sep=''),res=400,width=12000, height=4000)
# print(combined_plot)
# 
# legend_text <- "FDR < 0.05"
# legend_grob <- legendGrob(
#   labels = paste0("* ", legend_text),
#   gp = gpar(fontsize = 12)
# )
# 
# # Draw the custom legend
# pushViewport(viewport(x = unit(0.35, "npc"), y = unit(1.47, "npc"), 
#                       just = c("center", "top")))
# grid.draw(legend_grob)
# popViewport()
# dev.off()




#################################################################################
### Bubble plots to show coefficients of significant metabolites for each species score ###
#################################################################################
mutual_MET<-unique(unlist(lapply(metab,rownames)))

k<-regress(new.dat.met,
           met,
           met_name,
           score,
           score,
           AdjVars,
           0.05,
           T)

df.mutual<-lapply(k,function(x) x[mutual_MET,])

df.beta<-cbind(
               # df.mutual$A1.scca.score$Estimate,
               # df.mutual$A2.scca.score$Estimate,
               df.mutual$M1.scca.score$Estimate,
               df.mutual$M2.scca.score$Estimate,               
               
               df.mutual$acesk.score.3$Estimate,
               df.mutual$aspart.score.3$Estimate,
               df.mutual$sach.score.3$Estimate,
               df.mutual$sucrl.score.3$Estimate,
               
               df.mutual$eryth.score.3$Estimate,
               df.mutual$malt.score.3$Estimate,
               df.mutual$mani.score.3$Estimate,
               df.mutual$sorb.score.3$Estimate,
               df.mutual$xyli.score.3$Estimate
)

df.beta<-as.data.frame(t(df.beta))
names(df.beta)<-paste(mutual_MET)
df.beta$AS<-
  c(
    # "Sugar alcohol\nenriched intake",
    # "Synthetic sweetener\nenriched intake",
    "Sugar alcohols",
    "Synthetic sweeteners",
    "Acesulfame K",
    "Aspartame",
    "Saccharin",
    "Sucralose",
    "Erythritol",
    "Maltitol",
    "Mannitol",
    "Sorbitol",
    "Xylitol"
  )

df.beta<-df.beta[,c("AS",mutual_MET)]
df.beta$type<-c(rep("CCA variates",2),
                rep("Synthetic sweeteners",4),
                rep("Sugar alcohols",5))
################################################
df.fdr <-cbind(
               # df.mutual$A1.scca.score$P.adj,
               # df.mutual$A2.scca.score$P.adj,
               df.mutual$M1.scca.score$P.adj,
               df.mutual$M2.scca.score$P.adj,
               
               df.mutual$acesk.score$P.adj,
               df.mutual$aspart.score$P.adj,
               df.mutual$sach.score$P.adj,
               df.mutual$sucrl.score$P.adj,
               
               df.mutual$eryth.score$P.adj,
               df.mutual$malt.score$P.adj,
               df.mutual$mani.score$P.adj,
               df.mutual$sorb.score$P.adj,
               df.mutual$xyli.score$P.adj
)
df.fdr<-as.data.frame(t(df.fdr))
names(df.fdr)<-paste(mutual_MET)
df.fdr$AS<-
  c(
    # "Sugar alcohol\nenriched intake",
    # "Synthetic sweetener\nenriched intake",
    "Sugar alcohols",
    "Synthetic sweeteners",
    "Acesulfame K",
    "Aspartame",
    "Saccharin",
    "Sucralose",
    "Erythritol",
    "Maltitol",
    "Mannitol",
    "Sorbitol",
    "Xylitol"
  )
df.fdr<-df.fdr[,c("AS",mutual_MET)]
df.fdr$type<-c(rep("CCA variates",2),
               rep("Synthetic sweeteners",4),
               rep("Sugar alcohols",5))
################################################
pcm.beta = reshape2::melt(df.beta, id = c("AS","type"))
names(pcm.beta)[4]<-"Beta"
pcm.fdr  = reshape2::melt(df.fdr, id = c("AS","type"))
names(pcm.fdr)[4]<-"FDR.P"
pcm.fdr$FDR.P<- -log10(pcm.fdr$FDR.P)

pcm<-with(pcm.fdr,cbind(pcm.beta,FDR.P))  

pcm$AS <- factor(pcm$AS,levels=c(
                                 # "Sugar alcohol\nenriched intake",
                                 # "Synthetic sweetener\nenriched intake",
                                 "Sugar alcohols",
                                 "Synthetic sweeteners",
                                 "Acesulfame K",
                                 "Aspartame",
                                 "Saccharin",
                                 "Sucralose",
                                 "Erythritol",
                                 "Maltitol",
                                 "Mannitol",
                                 "Sorbitol",
                                 "Xylitol"
)
)

pcm$type<-factor(pcm$type,levels= c("CCA variates","Synthetic sweeteners","Sugar alcohols"))

met_name_class<-lvs_hpfs_fea[,c("metabolite_name","super_class_metabolon")]
names(met_name_class)<-c("variable","class_met")

pcm<-merge(pcm,met_name_class,by="variable")
pcm[which(pcm$variable=="C20:1 LPC"),"class_met"]<-"Lysophosphatidylcholines"
pcm[which(pcm$variable=="C22:4 LPC"),"class_met"]<-"Lysophosphatidylcholines"
pcm[which(pcm$variable=="LPE(20:0)"),"class_met"]<-"Lysophosphatidylcholines"
# pcm[which(pcm$variable=="3-Methyl-2-oxovaleric acid or Ketoleucine"),"variable"]<-factor("3-Methyl-2-oxovaleric acid")
pcm[which(pcm$variable=="3-Methyl-2-oxovaleric acid or Ketoleucine"),"class_met"]<-"Organic acids and derivatives"

pcm$class_met<-factor(pcm$class_met,levels=unique(pcm$class_met))

pcm<-pcm[order(pcm$class_met),]

pcm$variable<-factor(pcm$variable,levels=unique(pcm$variable))

write.csv(pcm,"CCA_MET.csv",row.names = F)

bubble.p<-ggplot(pcm, aes(x = AS, y = variable)) + 
  geom_point(aes(size = FDR.P, fill = Beta), alpha = 0.75, shape = 21,stroke=1) + 
  scale_size_continuous(limits = c(0, 8), range = c(1,20), breaks = c(0,2,4,6)) +
  labs( x= "Species score", y = "Plasma metabolites", size = "-Log10 (q-value)", fill = "Beta coefficients")  + 
  theme(legend.key=element_blank(), 
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 14, face ="bold", colour ="black"), 
        legend.title = element_text(size = 17, face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = "NA", linewidth = 1.2), 
        panel.grid.major.y = element_line(color = "pink",linewidth = 0.5,linetype = 2),
        legend.position = "right") +  
  scale_fill_distiller(palette = "RdPu", direction= 1)+
  scale_y_discrete(limits = rev(levels(pcm$variable)))+
  new_scale("fill")+ 
  geom_tile(aes(fill=class_met),data=pcm,alpha=0.3)+
  scale_fill_viridis_d(option = "H", name = "Metabolites Class")+
  labs(fill="") 
png(paste("~/CCA_MET.png",sep=''),res=400,width=12000, height=7500)
bubble.p
dev.off()


pcm_split<-pcm[!pcm$AS %in% c("Sugar alcohols", "Synthetic sweeteners"), ]

bubble.p.split<-ggplot(pcm_split, aes(x = AS, y = variable)) + 
  geom_point(aes(size = FDR.P, fill = Beta), alpha = 0.75, shape = 21,stroke=1) + 
  scale_size_continuous(limits = c(0, 8), range = c(1,20), breaks = c(0,2,4,6)) +
  labs( x= "", y = "Plasma metabolites", size = "-Log10 (q-value)", fill = "Beta coefficients")  + 
  theme(legend.key=element_blank(), 
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 14, face ="bold", colour ="black"), 
        legend.title = element_text(size = 17, face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = "NA", linewidth = 1.2), 
        panel.grid.major.y = element_line(color = "pink",linewidth = 0.5,linetype = 2),
        legend.position = "right") +  
  scale_fill_distiller(palette = "RdPu", direction= 1)+
  scale_y_discrete(limits = rev(levels(pcm_split$variable)))+
  new_scale("fill")+ 
  geom_tile(aes(fill=class_met),data=pcm_split,alpha=0.3)+
  scale_fill_viridis_d(option = "H", name = "Metabolites Class")+
  labs(fill="") 

ggsave(
  filename = "~/CCA_MET.svg",
  plot = bubble.p.split,
  device = "svg",
  width = 30,
  height = 18.75,
  units = "in"
)
############################################################################
### Figure for food contributors ###
############################################################################
library(ggalluvial)
library(grid)

comp<-read.csv("~//food con.csv",header=T) #csv from Pengfei
comp<-comp %>% filter(r> 0.1)

comp$AS<-factor(comp$AS,levels=unique(comp$AS)) 
comp$Foods<-factor(comp$Foods,levels=unique(comp$Foods))
comp$Groups<-factor(comp$Groups,levels=unique(comp$Groups)) 
comp$Type<-factor(comp$Type,levels=unique(comp$Type)) 


comp$adjusted_r <- comp$r  # Copy the original 'r'

# Add a constant value to create space between groups
spacing_factor <- 0.6  # Adjust as needed
comp$adjusted_r <- comp$r + spacing_factor * as.numeric(comp$Groups)

# Plot with adjusted y-values
allu <- ggplot(as.data.frame(comp),
               aes(y = adjusted_r, axis1 = AS, axis2 = Groups, axis3 = Foods)) +
  geom_alluvium(aes(fill = r), width = 1/3) +
  geom_stratum(
    width = 5/12,
    fill = c(rep("lightgreen", 5), 
             rep("orange", 4), 
             rep("pink", length(unique(comp$Groups))), 
             rep("grey", length(unique(comp$Foods))))
  ) +
  geom_text(
    stat = "stratum", 
    aes(label = after_stat(stratum)),
    color = "black",
    size = c(rep(7, 16),
             rep(5, 25))
  ) +
  scale_x_discrete(
    limits = c("Artificial sweeteners", "Food groups", "Individual foods"), 
    expand = c(.05, .05)
  ) +
  scale_fill_gradientn(
    colours = rev(brewer.pal(10, "Spectral"))
  ) +
  labs(fill = 'Correlation coefficients') +
  ylab("") +
  ggtitle("") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.margin = margin(3, 3, 3, 3),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.position = "right",
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_blank(),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(size = 15)
  )

# Create and position the legend
legend.grob1 <- legendGrob(
  labels = c("Synthetic sweeteners", "Sugar alcohols"),
  pch = 16,
  gp = gpar(col = c("orange", "lightgreen"))
)


png(paste("~/food.con.png",sep=''),res=400,width=9500, height=6500)
allu
pushViewport(viewport(x = 0.85, y = 0.4, width = 0.2, height = 0.2, just = c("left", "top")))
grid.draw(legend.grob1)
dev.off()

library(networkD3)
library(htmltools)
library(scales)
library(RColorBrewer)
library(htmlwidgets)

comp <- read.csv(
  "~/food con.csv",
  header = TRUE,
  stringsAsFactors = FALSE) 
# %>%
# filter(r > 0.1)

# Clean whitespace & NBSPs in the key columns
comp <- comp %>%
  mutate(
    AS     = str_squish(gsub("\u00A0", " ", AS)),
    Groups = str_squish(gsub("\u00A0", " ", Groups)),
    Foods  = str_squish(gsub("\u00A0", " ", Foods))
  )

# Build links
edges1 <- comp %>% transmute(source = AS,     target = Groups, r)
edges2 <- comp %>% transmute(source = Groups, target = Foods,  r)
links <- bind_rows(edges1, edges2) %>% mutate(value = r)

# ---- Explicitly set the AS order you want ----
# Replace with the exact order you want; here I put Sucralose before Erythritol.
desired_as_order <- c(
  "Acesulfame K", "Aspartame", "Saccharin",
  "Sucralose", "Erythritol"    # <-- ensure this exact order
)

# If some AS levels are not in desired list, append remaining ones (preserve their original order)
remaining_as <- setdiff(unique(comp$AS), desired_as_order)
as_levels <- c(desired_as_order, remaining_as)

# Groups and Foods - preserve their first appearance order
group_levels <- unique(comp$Groups)
food_levels  <- unique(comp$Foods)

# Build final node order: all AS, then Groups, then Foods
node_order <- unique(c(as_levels, group_levels, food_levels))

# Create nodes data.frame in that order
nodes <- data.frame(name = node_order, stringsAsFactors = FALSE)

# Assign node groups
nodes$group <- case_when(
  nodes$name %in% as_levels    ~ "AS",
  nodes$name %in% group_levels ~ "Group",
  nodes$name %in% food_levels  ~ "Food",
  TRUE ~ "Other"
)

# Rebuild 0-based ids based on the ordered nodes
name_to_id <- setNames(seq_len(nrow(nodes)) - 1, nodes$name)
links$source_id <- unname(name_to_id[links$source])
links$target_id <- unname(name_to_id[links$target])

# If any NA ids produced, inspect mismatches
if (any(is.na(links$source_id)) || any(is.na(links$target_id))) {
  stop("Some link names did not match the nodes. Check for typos or hidden characters.")
}

# Color mapping (same approach as you had)
r_scaled  <- scales::rescale(links$r, to = c(0, 1))
n_colors  <- 100
palette   <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(n_colors)
r_bins    <- cut(r_scaled, breaks = n_colors, labels = FALSE, include.lowest = TRUE)
links$group <- paste0("r", r_bins)

link_colors <- setNames(palette, paste0("r", 1:n_colors))
node_colors <- c(AS = "lightgreen", Group = "orange", Food = "cyan", Other = "grey")
all_colors  <- c(node_colors, link_colors)

color_scale <- JS(sprintf(
  "d3.scaleOrdinal().domain(%s).range(%s)",
  jsonlite::toJSON(names(all_colors)),
  jsonlite::toJSON(unname(all_colors))
))

# Create sankey with iterations = 0 to limit layout relaxation (helps preserve order)
sankey <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source_id",
  Target = "target_id",
  Value  = "value",
  NodeID = "name",
  NodeGroup = "group",
  LinkGroup = "group",
  fontSize = 15,
  nodeWidth = 20,
  sinksRight = FALSE,
  colourScale = color_scale,
  iterations = 0   # important: avoid extra relaxation that can reshuffle nodes
)

legend_html <- tags$div(
  style = "margin-bottom: 20px; text-align: center;",  # center everything
  tags$b("Correlation coefficient:"),
  tags$div(
    style = paste(
      "height: 20px;",
      "width: 400px;",
      sprintf("background: linear-gradient(to right, %s);",
              paste(palette, collapse = ", ")),
      "border: 1px solid #ccc;",
      "margin: 8px auto;"  # auto margins = horizontal centering
    )
  ),
  tags$div(
    style = "display: flex; justify-content: space-between; width: 400px; margin: 0 auto;",
    tags$span(sprintf("%.2f", min(links$r))),
    tags$span(sprintf("%.2f", max(links$r)))
  )
)

sankey_container <- tags$div(
  style = "position: relative; width: 960px; margin: auto;",
  sankey
)

final_plot <- browsable(tagList(
  legend_html,       # top color bar
  sankey_container  # sankey wrapped for layout
  
))
########################################################################################################
df.beta.M1<-read.csv("~//M1_met.csv", header=T)
names(df.beta.M1)[1:2]<-c("Metabolites","Beta")
df.beta.M1$type<-"Sugar alcohols"
df.beta.M2<-read.csv("~//M2_met.csv", header=T)
names(df.beta.M2)[1:2]<-c("Metabolites","Beta")
df.beta.M2$type<-"Synthetic sweeteners"

df.beta.M12<-rbind(df.beta.M1,df.beta.M2)
########################################################################################################
df <- read.csv("~//indv_met_hr_rev.csv", header=T)

df <- df %>%
  separate(
    col = Per.SD,
    into = c("HR", "CI"),    # First, split into HR and CI
    sep = " \\(",            # Split by " (" (space + open parenthesis)
    remove = TRUE
  ) %>%
  separate(
    col = CI,
    into = c("Lower_CI", "Upper_CI"), # Then split CI into Lower_CI and Upper_CI
    sep = ",|\\)",           # Split by "," or ")"
  )
df$Outcome<-factor(df$Outcome,levels = unique(df$Outcome))
df$Metabolites<-factor(df$Metabolites,levels = unique(df$Metabolites))

df$HR <- as.numeric(as.character(df$HR))
df$Lower_CI <- as.numeric(as.character(df$Lower_CI))
df$Upper_CI <- as.numeric(as.character(df$Upper_CI))

df.beta.T2D<-left_join(df.beta.M12,df[which(df$Outcome=="T2D"),],by="Metabolites")
df.beta.T2D<-df.beta.T2D[,c("Metabolites","Beta","Outcome","HR","Lower_CI", "Upper_CI", "type")]
df.beta.T2D<-df.beta.T2D[!is.na(df.beta.T2D$HR),]

df.beta.CHD<-left_join(df.beta.M12,df[which(df$Outcome=="CHD"),],by="Metabolites")
df.beta.CHD<-df.beta.CHD[,c("Metabolites","Beta","Outcome","HR","Lower_CI", "Upper_CI", "type")]
df.beta.CHD<-df.beta.CHD[!is.na(df.beta.CHD$HR),]

png(paste("~/quarant_T2D_notext.png",sep=''),res=400,width=5000, height=3500)
ggplot(df.beta.T2D, aes(x = Beta, y = HR, color = type)) +
  # error bars dodged
  geom_errorbar(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    width = 0.01
  ) +
  scale_color_manual(name="Species score-associated metabolites",
                     labels=c("Synthetic sweeteners", "Sugar alcohols"),
                     values=c("#D55E00", "#0072B2")) +
  # points dodged
  geom_point(size = 2) +
  # quadrant lines
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # labels also dodged so they follow their points
  # geom_text(aes(label = Metabolites, y = Lower_CI),
  #           nudge_y = -0.02,   # move slightly below the lower CI
  #           size = 3,
  #           show.legend = FALSE) +
  labs(
    x = "Beta coefficients",
    y = "Hazard Ratio of T2D"
  ) +
  theme_classic(base_size = 12)+
  theme(legend.position = "top")
dev.off()

png(paste("~/quarant_CHD_notext.png",sep=''),res=400,width=5000, height=3500)
ggplot(df.beta.CHD, aes(x = Beta, y = HR, color = type)) +
  # error bars dodged
  geom_errorbar(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    width = 0.01
  ) +
  scale_color_manual(name="Species score-associated metabolites",
                     labels=c("Synthetic sweeteners", "Sugar alcohols"),
                     values=c("#D55E00", "#0072B2")) +
  # points dodged
  geom_point(size = 2) +
  # quadrant lines
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # labels also dodged so they follow their points
  geom_text(aes(label = Metabolites, y = Lower_CI),
            nudge_y = -0.01,   # move slightly below the lower CI

            size = 3,
            show.legend = FALSE) +

  labs(
    x = "Beta coefficients",
    y = "Hazard Ratio of CHD"
  ) +
  theme_classic(base_size = 12)+
  theme(legend.position = "top")
dev.off()


png(paste("~/metabolites.png",sep=''),res=400,width=6000, height=3500)
ggplot(df, aes(x = Metabolites, y = HR)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) + 
  facet_wrap2(~ Outcome, ncol = 1,
              strip = strip_themed(
                background_x = list(
                  "Outcome1" = element_rect(fill = "lightblue"),
                  "Outcome2" = element_rect(fill = "lightgreen")
                )
              )) + # Create separate plots for each outcome
  labs(
    x = "Metabolites",
    y = "Hazard Ratio (Per SD)",
    title = 
      ""
  ) +
  scale_y_continuous(
    breaks = seq(0.6, 1.8, by = 0.4), # Custom y-axis breaks (from 0 to 2, step 0.2)
    limits = c(0.6, 1.8)              # Set y-axis range (optional)
  ) +
  theme_classic() + # Minimalistic theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
    strip.text = element_text(size = 16, face = "bold", family = "serif")
  )
dev.off()

png(paste("~/metabolite_forest.png",sep=''),res=400,width=5000, height=6000)
ggplot(df, aes(y = Metabolites, x = HR)) +
  geom_point(size = 3, color = "orange") +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, color = "black") +  # horizontal error bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) + 
  facet_wrap2(~ Outcome, ncol = 1,
              strip = strip_themed(
                background_x = list(
                  "Outcome1" = element_rect(fill = "lightblue"),
                  "Outcome2" = element_rect(fill = "lightgreen")
                )
              )) +
  labs(
    y = "",
    x = "Hazard Ratio (Per SD)",
    title = ""
  ) +
  coord_cartesian(xlim = c(0.6, 2))+
  theme_classic() +
  scale_x_log10() +
  theme(
    axis.text.y = element_text(size=20, face = "bold", color="black"),
    axis.text.x = element_text(size=20, face = "bold", color="black"),
    axis.title.x = element_text(size=22, face = "bold", color="black"),
    strip.text = element_text(size = 24, face = "bold", family = "serif")
  )
dev.off()
################################################################################################
path.all<-read.csv("~//path_all.csv", header=T)

df_path <- path.all

df_path_long <- df_path %>%
  pivot_longer(
    cols = ends_with("_Est"),
    names_to = "Sweetener",
    names_pattern = "(.*)_Est",
    values_to = "Effect"
  ) %>%
  mutate(Sweetener = factor(Sweetener, levels = unique(Sweetener)) ) %>%
  left_join(
    df_path %>%
      pivot_longer(
        cols = ends_with("_FDR"),
        names_to = "Sweetener",
        names_pattern = "(.*)_FDR",
        values_to = "FDR"
      ),
    by = c("Sweetener", "pathway")
  ) %>%
  filter(!is.na(FDR), 
         !is.na(Effect),
         FDR<0.25) %>%
  mutate(logFDR = -log10(FDR),
         # pathway = str_replace(pathway, ".*[.|_]s__", ""),
         Sweetener = factor(Sweetener,
                            levels = c(
                              "acesk_avgall", "aspart_avgall", "sach_avgall", "sucral_avgall",
                              "eryth_avgall", "maltitol_avgall", "mani_avgall", "sorb_avgall"
                            ),
                            labels = c(
                              "Acesulfame K", "Aspartame", "Saccharin", "Sucralose",
                              "Erythritol", "Maltitol", "Mannitol", "Sorbitol"
                            )
         ))

png(paste("~/path.png",sep=''),res=400,width=8000, height=6000)
ggplot(df_path_long, aes(x = Sweetener, y = pathway)) +
  geom_point(aes(color = Effect, size = logFDR)) +
  scale_color_gradient2(
    low = "#2166AC",   # deep blue
    mid = "white",
    high = "#B2182B",  # dark red
    midpoint = 0
  ) +
  geom_text(
    data = df_path_long %>% filter(FDR < 0.05),
    aes(label = "*"),
    color = "black",
    size = 5,
    hjust = -1.5,        # Push text to the right
    vjust = 0,         # Vertically centered
    fontface = "bold"
  )+
  scale_size_continuous(range = c(3, 12)) +
  labs(
    title = "",
    x = "Artificial Sweeteners",
    y = "",
    color = "Beta coefficients",
    size = "-Log10 (q-value)"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        fill = NA,           # Make legend circle hollow
        shape = 21,          # Circle with outline
        color = "black",     # Border color
        stroke = 1           # Border thickness
      )
    )
  )+
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 17,color="black", face="bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1, face="bold", color="black"),
    axis.text.y = element_text(size = 15, face="bold", color="black"),
    legend.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12,face = "bold"),
    panel.grid.major = element_line(color = "gray80"),
    panel.background = element_rect(fill = "#F0F0F0", color = NA)
  )
dev.off()
########################################################################################################################################
enz.all<-read.csv("~//enz_all.csv", header=T)

df_enz <- enz.all

df_enz_long <- df_enz %>%
  pivot_longer(
    cols = ends_with("_Est"),
    names_to = "Sweetener",
    names_pattern = "(.*)_Est",
    values_to = "Effect"
  ) %>%
  mutate(Sweetener = factor(Sweetener, levels = unique(Sweetener)) ) %>%
  left_join(
    df_enz %>%
      pivot_longer(
        cols = ends_with("_FDR"),
        names_to = "Sweetener",
        names_pattern = "(.*)_FDR",
        values_to = "FDR"
      ),
    by = c("Sweetener", "enzyme")
  ) %>%
  filter(!is.na(FDR), 
         !is.na(Effect),
         FDR<0.25) %>%
  mutate(logFDR = -log10(FDR),
         # enzyme = str_replace(enzyme, ".*[.|_]s__", ""),
         Sweetener = factor(Sweetener,
                            levels = c(
                              "acesk_avgall", "aspart_avgall", "sach_avgall", "sucral_avgall",
                              "eryth_avgall", "maltitol_avgall", "mani_avgall", "sorb_avgall"
                            ),
                            labels = c(
                              "Acesulfame K", "Aspartame", "Saccharin", "Sucralose",
                              "Erythritol", "Maltitol", "Mannitol", "Sorbitol"
                            )
         ))
png(paste("~/enz.png",sep=''),res=400,width=8000, height=6000)
ggplot(df_enz_long, aes(x = Sweetener, y = enzyme)) +
  geom_point(aes(color = Effect, size = logFDR)) +
  scale_color_gradient2(
    low = "#2166AC",   # deep blue
    mid = "white",
    high = "#B2182B",  # dark red
    midpoint = 0
  ) +
  geom_text(
    data = df_enz_long %>% filter(FDR < 0.05),
    aes(label = "*"),
    color = "black",
    size = 5,
    hjust = -0.5,        # Push text to the right
    vjust = 0,         # Vertically centered
    fontface = "bold"
  )+
  scale_size_continuous(range = c(3, 12)) +
  labs(
    title = "",
    x = "Artificial Sweeteners",
    y = "",
    color = "Beta coefficients",
    size = "-Log10 (q-value)"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        fill = NA,           # Make legend circle hollow
        shape = 21,          # Circle with outline
        color = "black",     # Border color
        stroke = 1           # Border thickness
      )
    )
  )+
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 17,color="black", face="bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1, face="bold", color="black"),
    axis.text.y = element_text(size = 15, face="bold", color="black"),
    legend.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12,face = "bold"),
    panel.grid.major = element_line(color = "gray80"),
    panel.background = element_rect(fill = "#F0F0F0", color = NA)
  )
dev.off()








