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
load(file="zero_0.90/INT_alpha_0.1_0.3.RData")
load(file="zero_0.90/INT_alpha_0.2_0.3.RData")
load(file="zero_0.90/INT_alpha_0.3_0.3.RData")
load(file="zero_0.90/INT_alpha_0.4_0.3.RData")
load(file="zero_0.90/INT_alpha_0.5_0.3.RData")
load(file="zero_0.90/INT_alpha_0.6_0.3.RData")
load(file="zero_0.90/INT_alpha_0.7_0.3.RData")
load(file="zero_0.90/INT_alpha_0.8_0.3.RData")
load(file="zero_0.90/INT_alpha_0.9_0.3.RData")
load(file="zero_0.90/INT_alpha_1_0.3.RData")

load(file="Taxon_use_AVG.RData")
load(file="TAXON_ONLY.RData")

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
acesk_INT.3.union<-Reduce(union,
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
length(acesk_INT.3)/length(acesk_INT.3.union)
# [1] 0.8135593

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
aspart_INT.3.union<-Reduce(union,
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
length(aspart_INT.3)/length(aspart_INT.3.union)
# [1] 0.4444444

aspart_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.2.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
                          clr_scale_avgall_INT.3.3$sp_id$aspart_avgall_INT[,aspart_INT.3],
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
sach_INT.3.union<-Reduce(union,
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
length(sach_INT.3)/length(sach_INT.3.union)
# [1] 0.4102564

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
sucral_INT.3.union<-Reduce(union,
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
length(sucral_INT.3)/length(sucral_INT.3.union)
# [1] 0.6216216

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
eryth_INT.3.union<-Reduce(union,
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
length(eryth_INT.3)/length(eryth_INT.3.union)
# [1] 0.4107143

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
maltitol_INT.3.union<-Reduce(union,
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
length(maltitol_INT.3)/length(maltitol_INT.3.union)
# [1] 0.7222222

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
mani_INT.3.union<-Reduce(union,
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
length(mani_INT.3)/length(mani_INT.3.union)
# [1] 0.78125

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
sorb_INT.3.union<-Reduce(union,
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
length(sorb_INT.3)/length(sorb_INT.3.union)
# [1] 0.6097561

sorb_INT_beta.3<- rbind(clr_scale_avgall_INT.1.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.2.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.3.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.4.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.5.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.6.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.7.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.8.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT.9.3$sp_id$sorb_avgall_INT[,sorb_INT.3],
                        clr_scale_avgall_INT1.3$sp_id$sorb_avgall_INT[,sorb_INT.3]
)

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
xyli_INT.3.union<-Reduce(union,
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
length(xyli_INT.3)/length(xyli_INT.3.union)
# [1] 0.7948718

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

write.csv(htmap_coef,"ELN_CLR_0.90.csv",row.names=T)
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
# use_clean: filtered relative-abundance matrix
# Rows = samples; columns = species
X <- as.matrix(use_clean)

stopifnot(
  all(is.finite(X)),
  all(X >= 0),
  all(rowSums(X) > 0),
  all(colSums(X > 0) > 0)
)

# Ensure each sample sums to one
X <- X / rowSums(X)

# Species-specific minimum positive abundance
dl <- apply(X, 2, function(x) min(x[x > 0]))

# Set to 0.40, 0.65, or 0.90 for the three analyses
replacement_fraction <- 0.9

# Multiplicative zero replacement
use_imp <- zCompositions::multRepl(
  X,
  label = 0,
  dl = dl,
  frac = replacement_fraction,
  z.warning = 0.80,
  z.delete = FALSE
)

# Verify that all values are suitable for log transformation
use_imp <- as.matrix(use_imp)
stopifnot(
  identical(dim(use_imp), dim(X)),
  all(is.finite(use_imp)),
  all(use_imp > 0),
  all(abs(rowSums(use_imp) - 1) < 1e-8)
)

# CLR transformation within each sample
log_imp <- log(use_imp)
use_clr <- as.data.frame(
  sweep(log_imp, 1, rowMeans(log_imp), FUN = "-")
)

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

AS.corr.3


########################
##    sCCA   ##
########################
X<-use_clr[,sp.3]
Y<-taxon.avg[,paste(AS_var, "avgall_INT",sep="")]
rownames(Y)<-taxon.avg$id

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

write.csv(M.score,"CLR_sCCA_coef_0.90.csv",row.names=T)



