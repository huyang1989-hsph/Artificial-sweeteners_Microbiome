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
library(zCompositions) # For compositional zero replacement
###################################################################################
setwd("~/NM_rev")
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

load(file="Taxon_use_AVG.RData")
load(file="TAXON_ONLY.RData")

load(file="~/data/MET_pm.RData")
load(file="~/data/MET_fea.RData")
###################################################################################
AS_var=c("acesk_","aspart_","sach_","sucral_",
         # "tag_",
         "eryth_",
         "maltitol_","mani_",
         # "pini_",
         "sorb_",
         "xyli_")

taxon.avg<-taxon.avg[-215,]
###################################################################################
acesk_INT.3<-Reduce(intersect,
                    list(
                      names(clr_scale_avgall_INT.1.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.2.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.3.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.4.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.5.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.6.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.7.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.8.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT.9.3$sp_id$acesk_avgall_INT),
                      names(clr_scale_avgall_INT1.3$sp_id$acesk_avgall_INT)
                    ) 
                    ,accumulate = FALSE)

acesk_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.2.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.3.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.4.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.5.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.6.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.7.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.8.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT.9.3$sp_id$acesk_avgall_INT[,acesk_INT.3],
                         clr_scale_avgall_INT1.3$sp_id$acesk_avgall_INT[,acesk_INT.3])

acesk_INT_beta_final.3<-apply(acesk_INT_beta.3,2,mean)
acesk_INT_beta_final.3<-as.data.frame(acesk_INT_beta_final.3)
acesk_INT_beta_final.3$microName<-rownames(acesk_INT_beta_final.3)
##################################################################################
aspart_INT.3<-Reduce(intersect,
                     list(
                       names(clr_scale_avgall_INT.1.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.2.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.3.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.4.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.5.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.6.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.7.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.8.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT.9.3$sp_id$aspart_avgall_INT),
                       names(clr_scale_avgall_INT1.3$sp_id$aspart_avgall_INT)
                     ) 
                     ,accumulate = FALSE)

aspart_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.2.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          # clr_scale_avgall_INT.3.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.4.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.5.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.6.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.7.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.8.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.9.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT1.3$sp_id$aspart_avgall_INT[,aspart_INT.3])

aspart_INT_beta_final.3<-apply(aspart_INT_beta.3,2,mean)
aspart_INT_beta_final.3<-as.data.frame(aspart_INT_beta_final.3)
aspart_INT_beta_final.3$microName<-rownames(aspart_INT_beta_final.3)
##################################################################################
sach_INT.3<-Reduce(intersect,
                   list(
                     names(clr_scale_avgall_INT.1.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.2.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.3.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.4.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.5.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.6.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.7.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.8.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT.9.3$sp_id$sach_avgall_INT),
                     names(clr_scale_avgall_INT1.3$sp_id$sach_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

sach_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.2.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.3.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.4.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.5.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.6.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.7.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.8.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT.9.3$sp_id$sach_avgall_INT[,sach_INT.3],
                        clr_scale_avgall_INT1.3$sp_id$sach_avgall_INT[,sach_INT.3])

sach_INT_beta_final.3<-apply(sach_INT_beta.3,2,mean)
sach_INT_beta_final.3<-as.data.frame(sach_INT_beta_final.3)
sach_INT_beta_final.3$microName<-rownames(sach_INT_beta_final.3)
##################################################################################
sucral_INT.3<-Reduce(intersect,
                     list(
                       names(clr_scale_avgall_INT.1.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.2.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.3.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.4.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.5.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.6.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.7.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.8.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT.9.3$sp_id$sucral_avgall_INT),
                       names(clr_scale_avgall_INT1.3$sp_id$sucral_avgall_INT)
                     ) 
                     ,accumulate = FALSE)

sucral_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.2.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.3.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.4.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.5.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.6.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.7.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.8.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT.9.3$sp_id$sucral_avgall_INT[,sucral_INT.3],
                          clr_scale_avgall_INT1.3$sp_id$sucral_avgall_INT[,sucral_INT.3])

sucral_INT_beta_final.3<-apply(sucral_INT_beta.3,2,mean)
sucral_INT_beta_final.3<-as.data.frame(sucral_INT_beta_final.3)
sucral_INT_beta_final.3$microName<-rownames(sucral_INT_beta_final.3)
##################################################################################
eryth_INT.3<-Reduce(intersect,
                    list(
                      names(clr_scale_avgall_INT.1.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.2.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.3.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.4.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.5.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.6.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.7.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.8.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT.9.3$sp_id$eryth_avgall_INT),
                      names(clr_scale_avgall_INT1.3$sp_id$eryth_avgall_INT)
                    ) 
                    ,accumulate = FALSE)

eryth_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.2.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.3.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.4.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.5.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.6.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.7.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.8.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT.9.3$sp_id$eryth_avgall_INT[,eryth_INT.3],
                         clr_scale_avgall_INT1.3$sp_id$eryth_avgall_INT[,eryth_INT.3])

eryth_INT_beta_final.3<-apply(eryth_INT_beta.3,2,mean)
eryth_INT_beta_final.3<-as.data.frame(eryth_INT_beta_final.3)
eryth_INT_beta_final.3$microName<-rownames(eryth_INT_beta_final.3)
##################################################################################
maltitol_INT.3<-Reduce(intersect,
                       list(
                         names(clr_scale_avgall_INT.1.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.2.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.3.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.4.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.5.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.6.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.7.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.8.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT.9.3$sp_id$maltitol_avgall_INT),
                         names(clr_scale_avgall_INT1.3$sp_id$maltitol_avgall_INT)
                       ) 
                       ,accumulate = FALSE)

maltitol_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.2.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.3.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.4.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.5.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.6.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.7.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.8.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT.9.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3],
                            clr_scale_avgall_INT1.3$sp_id$maltitol_avgall_INT[,maltitol_INT.3])

maltitol_INT_beta_final.3<-apply(maltitol_INT_beta.3,2,mean)
maltitol_INT_beta_final.3<-as.data.frame(maltitol_INT_beta_final.3)
maltitol_INT_beta_final.3$microName<-rownames(maltitol_INT_beta_final.3)
##################################################################################
mani_INT.3<-Reduce(intersect,
                   list(
                     names(clr_scale_avgall_INT.1.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.2.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.3.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.4.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.5.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.6.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.7.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.8.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT.9.3$sp_id$mani_avgall_INT),
                     names(clr_scale_avgall_INT1.3$sp_id$mani_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

mani_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.2.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.3.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.4.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.5.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.6.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.7.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.8.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT.9.3$sp_id$mani_avgall_INT[,mani_INT.3],
                        clr_scale_avgall_INT1.3$sp_id$mani_avgall_INT[,mani_INT.3])

mani_INT_beta_final.3<-apply(mani_INT_beta.3,2,mean)
mani_INT_beta_final.3<-as.data.frame(mani_INT_beta_final.3)
mani_INT_beta_final.3$microName<-rownames(mani_INT_beta_final.3)
##################################################################################
sorb_INT.3<-Reduce(intersect,
                   list(
                     names(clr_scale_avgall_INT.1.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.2.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.3.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.4.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.5.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.6.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.7.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.8.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT.9.3$sp_id$sorb_avgall_INT),
                     names(clr_scale_avgall_INT1.3$sp_id$sorb_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

sorb_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.2.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.3.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.4.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.5.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.6.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.7.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.8.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.9.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT1.3$sp_id$sorb_avgall_INT[,sorb_INT.3])

sorb_INT_beta_final.3<-apply(sorb_INT_beta.3,2,mean)
sorb_INT_beta_final.3<-as.data.frame(sorb_INT_beta_final.3)
sorb_INT_beta_final.3$microName<-rownames(sorb_INT_beta_final.3)
##################################################################################
xyli_INT.3<-Reduce(intersect,
                   list(
                     names(clr_scale_avgall_INT.1.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.2.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.3.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.4.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.5.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.6.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.7.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.8.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT.9.3$sp_id$xyli_avgall_INT),
                     names(clr_scale_avgall_INT1.3$sp_id$xyli_avgall_INT)
                   ) 
                   ,accumulate = FALSE)

xyli_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.2.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.3.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.4.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.5.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.6.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.7.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.8.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT.9.3$sp_id$xyli_avgall_INT[,xyli_INT.3],
                        clr_scale_avgall_INT1.3$sp_id$xyli_avgall_INT[,xyli_INT.3])

xyli_INT_beta_final.3<-apply(xyli_INT_beta.3,2,mean)
xyli_INT_beta_final.3<-as.data.frame(xyli_INT_beta_final.3)
xyli_INT_beta_final.3$microName<-rownames(xyli_INT_beta_final.3)
####################################################################################
sp.3<-unique(c(acesk_INT_beta_final.3$microName,
               aspart_INT_beta_final.3$microName,
               sach_INT_beta_final.3$microName,
               sucral_INT_beta_final.3$microName,
               eryth_INT_beta_final.3$microName,  
               maltitol_INT_beta_final.3$microName,
               mani_INT_beta_final.3$microName,
               sorb_INT_beta_final.3$microName,
               xyli_INT_beta_final.3$microName
)
)
length(sp.3)

annot_sp<-Annot_Taxon[which(Annot_Taxon$microName %in% sp.3),]
annot_sp<-annot_sp[,c("microName","phylum","species")]

annot_sp_sig_FDR<-list(annot_sp, 
                       # abd[,c("microName","abundance")],
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

write.csv(htmap_coef,"ring_CLR.csv",row.names=T)

ring_AST<-read.csv("~/data/taxon/0.01%/ring.csv",header=T)

a<-setdiff(ring_AST$X,rownames(htmap_coef))
b<-setdiff(rownames(htmap_coef),ring_AST$X)

intersect(rownames(htmap_coef),ring_AST$X)


###############################################################
INT<-function(dta,X)
{
  
  for (i in 1:length(X)){
    
    dta[,paste(X[i],"_INT",sep="")]<-qnorm((rank(dta[,X[i]],na.last="keep")-0.5)/sum(!is.na(dta[,X[i]])))
    
  }
  
  return(dta)
  
}
taxon.avg<-INT(taxon.avg,paste(AS_var, "avgall"   ,sep=""))
###############################################################
########################################################################################
# 1. Calculate abundance and detection rate
count = taxon.avg[,Annot_Taxon$microName]
count[count>=0.0001]=1
count[count<0.0001]=0

abd = data.frame(abundence=colMeans(taxon.avg[,Annot_Taxon$microName]), 
                 detection=colSums(count)/dim(count)[1])
abd$microName = rownames(abd)

abd = merge(abd,Annot_Taxon,by="microName")
head(abd)
dim(abd)

########################################################################################
# 2. Filter species by type AND detection threshold BEFORE transformation
# Change 0.1 to your desired detection threshold (e.g., >10% prevalence)
detection_cutoff = 0.3 

species_keep = abd$microName[abd$type == "species" & abd$detection >= detection_cutoff]

# Subset taxon.avg to filtered species only
use = taxon.avg[, species_keep]
head(use)[, 1:10]
dim(use)

########################################################################################
# 3. Re-normalize to sum to 1 (Closure) after filtering rare taxa
use = use / rowSums(use)

# 4. Zero Imputation (cmultRepl handles zeros prior to log ratio calculation)
# Note: cmultRepl expects samples as rows and features as columns
# Calculate the proportion of zeros per row
zero_prop <- rowMeans(use == 0)

table(zero_prop <= 0.80)

# Keep only rows with <= 80% zeros
use_clean <- use[zero_prop <= 0.80, ]

# Run CZM on the dynamically cleaned matrix
use_imp <- cmultRepl(use_clean, label = 0, method = "CZM")

# 5. Perform Centered Log-Ratio (CLR) Transformation
# CLR formula: log(x) - mean(log(x)) for each sample row
use_clr = t(apply(use_imp, 1, function(x) log(x) - mean(log(x))))
use_clr = as.data.frame(use_clr)


# use_clr_scale = as.data.frame(scale(use_clr))


acesk.score.3<- as.matrix(scale(use_clr[acesk_INT_beta_final.3$microName]))%*% acesk_INT_beta_final.3$acesk_INT_beta_final.3  ;colnames(acesk.score.3)<-"acesk.score.3" 
aspart.score.3<-as.matrix(scale(use_clr[aspart_INT_beta_final.3$microName]))%*% aspart_INT_beta_final.3$aspart_INT_beta_final.3 ;colnames(aspart.score.3)<-"aspart.score.3" 
sach.score.3<-  as.matrix(scale(use_clr[sach_INT_beta_final.3$microName]))%*% sach_INT_beta_final.3$sach_INT_beta_final.3   ;colnames(sach.score.3)<-"sach.score.3" 
sucrl.score.3<- as.matrix(scale(use_clr[sucral_INT_beta_final.3$microName]))%*% sucral_INT_beta_final.3$sucral_INT_beta_final.3 ;colnames(sucrl.score.3)<-"sucrl.score.3"  

eryth.score.3<-as.matrix(scale(use_clr[eryth_INT_beta_final.3$microName]))%*% eryth_INT_beta_final.3$eryth_INT_beta_final.3    ;colnames(eryth.score.3)<-"eryth.score.3" 
malt.score.3<- as.matrix(scale(use_clr[maltitol_INT_beta_final.3$microName]))%*% maltitol_INT_beta_final.3$maltitol_INT_beta_final.3 ;colnames(malt.score.3)<-"malt.score.3" 
mani.score.3<- as.matrix(scale(use_clr[mani_INT_beta_final.3$microName]))%*% mani_INT_beta_final.3$mani_INT_beta_final.3     ;colnames(mani.score.3)<-"mani.score.3" 
sorb.score.3<- as.matrix(scale(use_clr[sorb_INT_beta_final.3$microName]))%*% sorb_INT_beta_final.3$sorb_INT_beta_final.3     ;colnames(sorb.score.3)<-"sorb.score.3" 
xyli.score.3<- as.matrix(scale(use_clr[xyli_INT_beta_final.3$microName]))%*% xyli_INT_beta_final.3$xyli_INT_beta_final.3     ;colnames(xyli.score.3)<-"xyli.score.3" 

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
########################
##    sCCA   ##
########################
X<-use_clr[,sp.3]
Y<-taxon.avg[,paste(AS_var, "avgall_INT",sep="")]

grid <- expand.grid(
  penaltyx = seq(0.3, 0.5, 0.05),
  penaltyz = seq(0.8, 1, 0.05)
)

set.seed(2026)

tune <- CCA.permute(
  x = X,
  z = Y,
  typex = "standard",
  typez = "standard",
  penaltyxs = grid$penaltyx,
  penaltyzs = grid$penaltyz,
  nperms = 100,
  niter = 10,
  standardize = TRUE
)

best_penaltyx <- tune$bestpenaltyx
best_penaltyz <- tune$bestpenaltyz

out <- CCA(
  X, Y,
  typex = "standard",
  typez = "standard",
  K = 9,
  penaltyx = best_penaltyx,
  penaltyz = best_penaltyz,
  standardize = TRUE
)

out <- CCA(X,Y,typex="standard",typez="standard",K=9,
           penaltyx=0.3,penaltyz=0.85)#0.35,0.85

score.x <- scale(X) %*% out$u
score.y <- scale(Y) %*% out$v
diag(cor(score.x, score.y))[1:2]

cor(scale(Y),score.x)[,1:2]

M.score<-as.data.frame(out$u)
names(M.score)<-paste("M",seq(1,9,1),sep="")
M.score$microName<-colnames(X)
rownames(M.score)<-colnames(X)
M.score<-M.score[which(M.score$M1 !=0 | M.score$M2!=0),]
M.score<-merge(M.score,Annot_Taxon[,c("microName","species")],by="microName")
dim(M.score)

write.csv(M.score,"CLR_ELN_coef.csv",row.names=T)

ast<-read.csv("AST_ELN_coef.csv",header = T)


setdiff(M.score$species,ast$species)
length(intersect(ast$species,M.score$species))/nrow(M.score)
###############################################################
##################################################################################
### CCA representation plot ###
##################################################################################
corr.Y.xscores<-cor(scale(Y),score.x)
Ynames<-rownames(corr.Y.xscores)

scca.plot<-as.data.frame(cbind(Ynames,corr.Y.xscores))
names(scca.plot)<-c("AS","Variate2", "Variate1")
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

png(paste("~/NM_rev/scca_repre.png",sep=''),res=400,width=4000, height=4000)
scca.p
dev.off()
########################################################################################
M2.scca.score<-score.x[,1]
M1.scca.score<-score.x[,2]

A2.scca.score<-score.y[,1]
A1.scca.score<-score.y[,2]

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

new.dat.met<-merge(new.dat,lvs_hpfs_pm,by="id")
met_exclude<-names(which(apply(is.na(lvs_hpfs_pm), 2, sum)>100))[-1]#Remove missing >100

met<-rownames(lvs_hpfs_fea)

met<-setdiff(met,met_exclude)

met_name<-lvs_hpfs_fea[which(!(rownames(lvs_hpfs_fea) %in% met_exclude)),"metabolite_name"]

AS_inv<-paste(AS_var, "avgall_INT",sep="")

score<-c("M1.scca.score","M2.scca.score")

new.dat.met[,score]<-scale(new.dat.met[,score])

AdjVars = c("age_fec",
            "bmi_bld",
            "totMETs_paq",
            "lt_ahei_dm",
            "calor_fo_dr_wtavg",
            "probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec",
            "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")
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
k<-regress(new.dat.met,
               met,
               met_name,
               score,
               score,
               AdjVars,
               0.05,
               T)

metab$M1.scca.score
# Estimate            P      met.id nmiss      P.adj
# 3-Indolepropionic acid 0.2148628 6.502558e-05 HMDB0002302     0 0.02347424

metab$M2.scca.score
# Estimate            P      met.id nmiss        P.adj
# vitamin A               0.2497616 3.967389e-05 HMDB0000305     0 3.580568e-03
# Hippuric acid          -0.2925425 3.466809e-07 HMDB0000714     0 4.171727e-05
# alpha-tocopherol        0.2603245 9.575177e-05 HMDB0001893     0 6.913278e-03
# 3-Indolepropionic acid -0.2881377 1.519055e-07 HMDB0002302     0 2.741894e-05
# cinnamoylglycine       -0.2893654 1.439813e-07 HMDB0011621     0 2.741894e-05
# N-acetylleucine         0.2508759 3.760670e-04 HMDB0011756     0 2.262670e-02


