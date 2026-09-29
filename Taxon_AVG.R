setwd("~/PR")
source("~/MLVSfunctions.R")
library(Maaslin2)
library(data.table)
library(ggplot2)
library(ggpubr)
library(sqldf)
library(plyr)
library(dplyr)
library(MASS)

options(bitmapType='cairo')
#--------------------------------------------------------------
# load data for analysis 
#--------------------------------------------------------------
lnames = load(file="TAXON_ONLY.RData")

Microbiome<-Annot_Taxon$microName

ffq_var<-names(Taxon_use)[grep(c("_ffq"),names(Taxon_use))]
ddr_var<-names(Taxon_use)[grep(c("avg|wtavg|_fo_dr"),names(Taxon_use))]
bio_var<-c("adj_choline",
           "adj_tmao",
           "adj_carnitine",
           "tc_plasma",
           "hdlc_plasma",
           "tg_plasma",
           "adj_crp",
           "tc.hdl",
           "hba1cp",
           "choline_std_bld",
           "tmao_std_bld",
           "carnitine_std_bld",
           "tc_std_bld",
           "hdlc_std_bld",     
           "tg_std_bld",
           "tc.hdl_std_bld",
           "crp_std_bld",
           "hba1c_std_bld")
basic_var<-c(
             "AnalySet",
             "SampleID",
             "IDTime",
             "fast_bld",
             "smoke_bld",
             "antibio_12m_fec",
             "colsc_2m_fec",
             "probio_2m_fec",
             "acid_2m_fec",
             "bile_2m_fec",
             "stooltype_fec.1",
             "stooltype_fec.2",
             "stooltype_fec.3",
             "stooltype_fec.4",
             "stooltype_fec.5",
             "stooltype_fec.6",
             "age_bld",
             "age_fec",
             "bmi_bld",
             "totMETs_paq",
             "lt_fiber",
             "lt_whgrn",         
             "lt_ahei",
             "lt_ahei_na",
             "lt_ahei_dm",
             "lt_ahei_nowgr",
             "lt_readmt",
             "lt_trypt")
ToBeUsed.Taxon<-Taxon_use[,c("id",
                             # "d1","d2","d3","d4","d5","d6","d7",
                             "flag1","flag7",Microbiome,basic_var,bio_var,ffq_var,ddr_var)]

taxon.avg = as.data.table(ToBeUsed.Taxon[,!names(ToBeUsed.Taxon) %in% names(ToBeUsed.Taxon)[which(sapply(ToBeUsed.Taxon, is.numeric) == FALSE)]])[, lapply(.SD,mean), by=id]

taxon.avg=as.data.frame(taxon.avg)

taxon.avg[which(taxon.avg[,"stooltype_fec.1"] > 0), "stooltype_fec.1"] = 1
taxon.avg[which(taxon.avg[,"stooltype_fec.2"] > 0), "stooltype_fec.2"] = 1
taxon.avg[which(taxon.avg[,"stooltype_fec.3"] > 0), "stooltype_fec.3"] = 1
taxon.avg[which(taxon.avg[,"stooltype_fec.4"] > 0), "stooltype_fec.4"] = 1
taxon.avg[which(taxon.avg[,"stooltype_fec.5"] > 0), "stooltype_fec.5"] = 1
taxon.avg[which(taxon.avg[,"stooltype_fec.6"] > 0), "stooltype_fec.6"] = 1
taxon.avg[which(taxon.avg[,"smoke_bld"] > 0), "smoke_bld"] = 1
taxon.avg[which(taxon.avg[,"colsc_2m_fec"] > 0), "colsc_2m_fec"] = 1
taxon.avg[which(taxon.avg[,"antibio_12m_fec"] > 0), "antibio_12m_fec"] = 1
taxon.avg[which(taxon.avg[,"probio_2m_fec"] > 0), "probio_2m_fec"] = 1
taxon.avg[which(taxon.avg[,"acid_2m_fec"] > 0), "acid_2m_fec"] = 1
taxon.avg[which(taxon.avg[,"bile_2m_fec"] > 0), "bile_2m_fec"] = 1

as_var=c("acesk_","aspart_","sach_","sucral_",#Synthetic sweetener
         # "isom_","lactl_",
         "eryth_","maltitol_","mani_","sorb_","xyli_","tag_")#Sugar alcohol 

pdi_var=c("pdi", "hpdi", "updi", "fhpdi", "yhpdi", "phpdi", "fpyhpdi", 
          "whlsvg", "frusvg", "vegsvg", "nutsvg", "legsvg", "tcfsvg", "frjsvg", "refsvg", "potsvg", 
          "ssbsvg", "swtsvg", "daisvg", "dmysvg", "yogsvg", "eggsvg", "fshsvg", "mtsvg", "rmtsvg", "pltsvg", "mscsvg", "oilsvg", "fatsvg")
d_use=c(
  "avg1d",
  "avg2d",
  "avg3d",
  "avg4d",
  "avg5d",
  "avg6d",
  "avg7d",
  "avgall",
  "fo_dr_wtavg")

var_use<-NULL

for (i in 1:length(as_var)){
  
  var<-paste(as_var[i], d_use,sep="")
  
  var_use<-c(var_use,var)
  
}
#############################################################################################################################
##Find how many AS has too much 0 intake, isomalt and lactitol should be exlcuded
most.0=sapply(taxon.avg[,var_use],function(x) table(x)[1]/length(x))
# sink("AS_too_many_0.csv")
# print(a)
# sink()
exclude.0=substr(names(most.0[most.0>0.99]),1,nchar(names(most.0[most.0>0.99]))-2)
# [1] "isom_avg1d"  "isom_avg2d"  "isom_avg3d"  "isom_avg4d"  "isom_avg5d"  "isom_avg6d"  "isom_avg7d"  "lactl_avg1d" "lactl_avg2d" "lactl_avg3d" "lactl_avg4d" "lactl_avg5d" "lactl_avg6d"
# [14] "lactl_avg7d"
########
#1. Isomalt and Lactitol are exlcuded due to mostly 0, hopeless
########
var_use<-var_use[! var_use %in% exclude.0]

AdjVars = c("age_fec","bmi_bld","totMETs_paq", 
            # "smoke_bld", 
            "lt_ahei_dm","calor_fo_dr_wtavg",
            "probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec",
            "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")

# table(is.na(taxon.avg[,var_use]))


summary(taxon.avg[,AdjVars])


# taxon.avg[1:10,c("eryth_avg1d","eryth_avg3d","eryth_fo_dr_w2d4","eryth_fo_dr_w2d5","eryth_fo_dr_w2d6","eryth_fo_dr_w1d7")]

dim(taxon.avg)
##################################################
Phenotypes = setdiff(names(taxon.avg),Microbiome)

save(taxon.avg,file="Taxon_use_AVG.RData")