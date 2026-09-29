library(survival) 
library(survminer)
library(mice)
library(tibble)
library(dplyr)
M1.met <- read.csv("~/M1_met.csv",header=T)
M2.met <- read.csv("~/M2_met.csv",header=T)

indiv_met1<-M1.met$met.id
indiv_met2<-M2.met$met.id

unique(c(indiv_met1,indiv_met2))
# [1] "HMDB0000305" "HMDB0000714" "HMDB0001893" "HMDB0002302" "HMDB0011621" "HMDB0011756" "HMDB0000019" "HMDB0000491" "HMDB0000885" "HMDB0001563" "HMDB0006736" "HMDB0007098" "HMDB0007100" "HMDB0007102" "HMDB0007103"
# [16] "HMDB0007199" "HMDB0007216" "HMDB0007218"
met.list<-unique(c(indiv_met1,indiv_met2))
met.list<-setdiff(met.list,"HMDB0001893") # HMDB0001893 was not measured in any of 3 cohorts

pool<-read.csv("~/AS_outcome.csv",header=T)

summary(pool[,met.list])

met.name<-c("vitA",
            "hipp_acid",
            "ipa",
            "c_glycine",
            "N_A_leucine",
            "k_isov",
            "k_leuc",
            "CE160",
            "m_guanosine",
            "CE20_3",
            "DG32_0",
            "DG34_0",
            "DG34_1",
            "DG34_2",
            "DG38_5",
            "DG36_1",
            "DG36_2")

pool[,met.name]<-pool[,met.list]

#Assign median value in each quantile
assign_med <- function(dta,       
                       var,      #Continuous variable
                       var.rank, #Quantile variable
                       prefix = "m",
                       n_quantile = 3) #Number of quantiles
{ # Add a parameter to define the number of quantiles
  
  # Create new column names
  mdq <- paste(prefix, var.rank, sep = "")
  
  # Loop through variables
  for (i in seq_along(var)) {
    
    # Calculate the median for each rank
    md <- tapply(dta[, var[i]], dta[, var.rank[i]], median, na.rm = TRUE)
    
    # Initialize the new column with NA
    dta[, mdq[i]] <- NA
    
    # Assign median values to the appropriate quantile group dynamically
    for (q in 1:n_quantile) {
      dta[, mdq[i]] <- ifelse(dta[, var.rank[i]] == q, md[q], dta[, mdq[i]])
    }
  }
  
  return(dta)
}

pool<-assign_med(pool,
                 met.name,
                 paste(met.name,"t",sep=""))

pool <- pool %>%
  mutate(vegqq = ntile(veg, 4),         
         vegqq = ifelse(is.na(veg), 1, vegqq),
         fruitqq= ntile(veg, 4),
         fruitqq = ifelse(is.na(veg), 1, fruitqq))   
#########################################################################################################################################################
covar<-"ageyr+factor(race)+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq)"

exposure<-met.name
exposure_q<-paste(met.name,"t",sep="")
exposure_trend<-paste("m",met.name,"t",sep="")

outcome<-c("dbd","db","tchd","chd"
           # ,"tstr","str","tcvd","cvd"
)
outcome_name<-c("T2D","CHD"
                # ,"STR","CVD"
)

exposure_f<-paste("factor(",exposure_q,")",sep="") 
exposure_sd<-paste("scale(",exposure,")",sep="") 

x<-NULL
for (i in 1: length(exposure_f)){
  xx<-paste(exposure_f[i],"+",covar)
  x<-c(x,xx)
}

x.trend<-NULL
for (i in 1: length(exposure_f)){
  xx.trend<-paste(exposure_trend[i],"+",covar)
  x.trend<-c(x.trend,xx.trend)
}

x.sd<-NULL
for (i in 1: length(exposure)){
  xx.sd<-paste(exposure_sd[i],"+",covar)
  x.sd<-c(x.sd,xx.sd)
}

y<-NULL
for (i in (seq(2,length(outcome),2)-1)){
  yy<-paste('Surv',"(",outcome[i],",",outcome[i+1],")","~",sep="")
  y<-c(y,yy)
}

a<-as.list(y)

b<-sapply(y,function(xx)(paste(xx,x)))
b.trend<-sapply(y,function(xx.trend)(paste(xx.trend,x.trend)))
b.sd<-sapply(y,function(xx.sd)(paste(xx.sd,x.sd)))

c <- sapply(b, as.formula)
c.trend <- sapply(b.trend, as.formula)
c.sd <- sapply(b.sd, as.formula)

d <- lapply(c, function(x) coxph(x, data = pool)) 
d.trend <- lapply(c.trend, function(x) coxph(x, data = pool)) 
d.sd <- lapply(c.sd, function(x) coxph(x, data = pool)) 

e<-lapply(d, summary) 
names(e)<-as.character(seq(1,length(e),1))
f<-lapply(as.character(seq(1,length(e),1)), function(x) as.data.frame(get("conf.int",get(x,e)))[1:2,c("exp(coef)","lower .95","upper .95")]) 
names(f)<-sapply(outcome_name,function(x)(paste(x,"_",exposure,sep="")))

e.trend<-lapply(d.trend, summary) 
names(e.trend)<-as.character(seq(1,length(e.trend),1))
f.trend<-lapply(as.character(seq(1,length(e.trend),1)), function(x) as.data.frame(get("coefficients",get(x,e.trend)))[1,5]) 
names(f.trend)<-sapply(outcome_name,function(x)(paste(x,"_",exposure,"_","trend",sep="")))

e.sd<-lapply(d.sd, summary) 
names(e.sd)<-as.character(seq(1,length(e.sd),1))
f.sd<-lapply(as.character(seq(1,length(e.sd),1)), function(x) as.data.frame(get("conf.int",get(x,e.sd)))[1,c("exp(coef)","lower .95","upper .95")]) 
names(f.sd)<-sapply(outcome_name,function(x)(paste(x,"_",exposure,"_","sd",sep="")))

results<-NULL
for (i in 1:(length(exposure)*length(outcome_name))){
  
  res<-as.data.frame(get(names(f)[i],f))
  res.trend<-round(get(names(f.trend)[i],f.trend),3)
  res.sd<-as.data.frame(get(names(f.sd)[i],f.sd))
  res$HR<-paste(round(res$`exp(coef)`,2)," (", round(res$`lower .95`,2), ",", round(res$`upper .95`,2), ")", sep="")
  res.sd$HR<-paste(round(res.sd$`exp(coef)`,2)," (", round(res.sd$`lower .95`,2), ",", round(res.sd$`upper .95`,2), ")", sep="")
  event<-names(f[i])
  ref<-"1.00(Reference)"
  res<-t(rbind(event,ref,res,res.trend,res.sd))[4,]
  if(as.numeric(res[5]<0.001)){res[5]<-"<0.001"}
  
  results<-rbind(results,res)
}
colnames(results)<-c("Outcome",paste("T",seq(1,3,1),sep=""),"P trend","Per SD")

results<-as.data.frame(results)

results
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# Strategy to handle missing, IPA should stay #
pool.IPA<- pool %>% 
  filter(if_all(c(HMDB0002302,HMDB0011621), ~ !is.na(.)))


table(pool.IPA$cohort)
# 1    2    3 
# 2234 2198  635
by(pool.IPA[,met.list], pool.IPA$cohort, summary)

## Remove HMDB0000305,HMDB0000019,HMDB0000491
met.list.rd<-setdiff(met.list,c("HMDB0000305","HMDB0000019","HMDB0000491"))

summary(pool.IPA[,met.list.rd])

# HMDB0000714        HMDB0002302        HMDB0011621        HMDB0011756       HMDB0000885        HMDB0001563         HMDB0006736         HMDB0007098       HMDB0007100       HMDB0007102        HMDB0007103      
# Min.   :-3.69913   Min.   :-3.61601   Min.   :-6.07283   Min.   :-4.2885   Min.   :-7.84315   Min.   :-13.77183   Min.   :-10.21306   Min.   :-3.9964   Min.   :-4.2367   Min.   :-3.20762   Min.   :-3.53290  
# 1st Qu.:-0.61204   1st Qu.:-0.54043   1st Qu.:-0.40240   1st Qu.:-0.5370   1st Qu.:-0.39693   1st Qu.: -0.47148   1st Qu.: -0.45866   1st Qu.:-0.7651   1st Qu.:-0.7620   1st Qu.:-0.73808   1st Qu.:-0.74906  
# Median : 0.03487   Median : 0.06864   Median : 0.11974   Median :-0.1349   Median : 0.20891   Median :  0.05842   Median :  0.13072   Median :-0.0653   Median :-0.1420   Median :-0.08794   Median :-0.03748  
# Mean   : 0.03702   Mean   : 0.00020   Mean   : 0.03251   Mean   :-0.0252   Mean   : 0.03667   Mean   : -0.02816   Mean   :  0.01516   Mean   :-0.0400   Mean   :-0.0327   Mean   :-0.05522   Mean   :-0.05298  
# 3rd Qu.: 0.71453   3rd Qu.: 0.62277   3rd Qu.: 0.62726   3rd Qu.: 0.3400   3rd Qu.: 0.64333   3rd Qu.:  0.48063   3rd Qu.:  0.64875   3rd Qu.: 0.6435   3rd Qu.: 0.5933   3rd Qu.: 0.62798   3rd Qu.: 0.63875  
# Max.   : 2.92773   Max.   : 3.66063   Max.   : 3.82625   Max.   : 8.9634   Max.   : 9.10829   Max.   :  3.24739   Max.   :  2.66652   Max.   : 3.2473   Max.   : 4.1653   Max.   : 3.22446   Max.   : 3.17408  
# NA's   :1458      NA's   :3                              NA's   :3           NA's   :660       NA's   :660       NA's   :3          NA's   :3         
#   HMDB0007199        HMDB0007216       HMDB0007218      
#  Min.   :-3.30794   Min.   :-4.0307   Min.   :-3.70748  
#  1st Qu.:-0.71342   1st Qu.:-0.7012   1st Qu.:-0.74602  
#  Median :-0.05308   Median :-0.0083   Median :-0.05198  
#  Mean   :-0.03770   Mean   :-0.0395   Mean   :-0.04460  
#  3rd Qu.: 0.60912   3rd Qu.: 0.6366   3rd Qu.: 0.62875  
#  Max.   : 4.00774   Max.   : 3.3872   Max.   : 3.76416  
#  NA's   :3          NA's   :660       NA's   :3     

# Impute missing with median
pool.IPA[met.list.rd] <- lapply(pool.IPA[met.list.rd], function(col) {
  ifelse(is.na(col), median(col, na.rm = TRUE), col)
})

summary(pool.IPA[,met.list.rd])

met.name<-c("hipp_acid",
            "ipa",
            "c_glycine",
            "N_A_leucine",
            "CE160",
            "m_guanosine",
            "CE20_3",
            "DG32_0",
            "DG34_0",
            "DG34_1",
            "DG34_2",
            "DG38_5",
            "DG36_1",
            "DG36_2")

pool.IPA[,met.name]<-pool.IPA[,met.list.rd]

pool.IPA <- pool.IPA %>%
  mutate(
    across(
      all_of(met.name), # Specify the columns to rank
      ~ ntile(., 3), # Divide values into 4 quartiles
      .names = "{col}_tt" # New column names
    )
  )

pool.IPA <- pool.IPA %>%
  mutate(
    across(
      all_of(met.name), # Specify the columns to rank
      ~ ntile(-., 3), # Multiply by -1 to reverse the ranking
      .names = "r_{col}_tt" # New column names
    )
  )

pool.IPA$M1.score<-pool.IPA$hipp_acid_tt+
  pool.IPA$ipa_tt+
  pool.IPA$c_glycine_tt+
  pool.IPA$r_N_A_leucine_tt #Reverse

# pool.IPA$M1.score.sd <- (pool.IPA$M1.score - mean(pool.IPA$M1.score, na.rm = TRUE)) / sd(pool.IPA$M1.score, na.rm = TRUE)

pool.IPA$M2.score<-pool.IPA$r_CE160_tt+        #Reverse
  pool.IPA$r_m_guanosine_tt+  #Reverse
  pool.IPA$r_ipa_tt+          #Reverse
  pool.IPA$r_CE20_3_tt+       #Reverse
  pool.IPA$DG32_0_tt+
  pool.IPA$DG34_0_tt+
  pool.IPA$DG34_1_tt+
  pool.IPA$DG34_2_tt+
  pool.IPA$DG38_5_tt+
  pool.IPA$DG36_1_tt+
  pool.IPA$DG36_2_tt+
  pool.IPA$c_glycine_tt

# pool.IPA$M2.score.sd <- (pool.IPA$M2.score - mean(pool.IPA$M2.score, na.rm = TRUE)) / sd(pool.IPA$M2.score, na.rm = TRUE)

#########################################################################################################################################################
covar<-"ageyr+factor(race)+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq)"

exposure<-c(met.name,"M1.score","M2.score")

outcome<-c("dbd","db","tchd","chd"
           # ,"tstr","str","tcvd","cvd"
)
outcome_name<-c("T2D","CHD"
                # ,"STR","CVD"
)

exposure_sd<-paste("scale(",exposure,")",sep="") 


x.sd<-NULL
for (i in 1: length(exposure)){
  xx.sd<-paste(exposure_sd[i],"+",covar)
  x.sd<-c(x.sd,xx.sd)
}

y<-NULL
for (i in (seq(2,length(outcome),2)-1)){
  yy<-paste('Surv',"(",outcome[i],",",outcome[i+1],")","~",sep="")
  y<-c(y,yy)
}

a<-as.list(y)

b.sd<-sapply(y,function(xx.sd)(paste(xx.sd,x.sd)))

c.sd <- sapply(b.sd, as.formula)

d.sd <- lapply(c.sd, function(x) coxph(x, data = pool.IPA)) 

e.sd<-lapply(d.sd, summary) 
names(e.sd)<-as.character(seq(1,length(e.sd),1))
f.sd<-lapply(as.character(seq(1,length(e.sd),1)), function(x) as.data.frame(get("conf.int",get(x,e.sd)))[1,c("exp(coef)","lower .95","upper .95")]) 
names(f.sd)<-sapply(outcome_name,function(x)(paste(x,"_",exposure,"_","sd",sep="")))

results<-NULL
for (i in 1:(length(exposure)*length(outcome_name))){
  
  res.sd<-as.data.frame(get(names(f.sd)[i],f.sd))
  res.sd$HR<-paste(round(res.sd$`exp(coef)`,2)," (", round(res.sd$`lower .95`,2), ",", round(res.sd$`upper .95`,2), ")", sep="")
  event<-names(f.sd[i])
  ref<-"1.00(Reference)"
  res<-t(rbind(event,res.sd))[4,]
  
  results<-rbind(results,res)
}
colnames(results)<-c("Outcome","Per SD")

results<-as.data.frame(results)


results

covar<-"ageyr+factor(race)+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq)"
####################################################################################
pool.IPA <- pool.IPA %>%
  mutate(q.M1.score = ntile(M1.score, 4),
         q.M2.score = ntile(M2.score, 4))

table(pool.IPA$q.M1.score)
table(pool.IPA$q.M2.score)

pool.IPA$white<-ifelse(pool.IPA$race=="White",1,0)

surv.ipw<-function(dta,exp){
  data<-dta[!is.na(dta[,exp]),]
  data$gg <- as.numeric(data[,exp])
  
  lfit1 <- glm(I(gg==1) ~ ageyr+white+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq), data=data,family="binomial")
  lfit2 <- glm(I(gg==2) ~ ageyr+white+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq), data=data,family="binomial")
  lfit3 <- glm(I(gg==3) ~ ageyr+white+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq), data=data,family="binomial")
  lfit4 <- glm(I(gg==4) ~ ageyr+white+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq), data=data,family="binomial")
  # lfit5 <- glm(I(gg==5) ~ ageyr+white+factor(cohort)+factor(actqq)+factor(alcoqq)+factor(ahei_noalqq)+smoke2+smoke3+bmic2+bmic3+bmicm+factor(calorqq), data=data,family="binomial")
  
  data$sbw <- ifelse(data$gg==1, 0.25/predict(lfit1, type='response'),
                     ifelse(data$gg==2, 0.25/predict(lfit2, type='response'),
                            ifelse(data$gg==3, 0.25/predict(lfit3, type='response'),
                                                     0.25/predict(lfit4, type='response'))))
  return(data)
}

M1.KM<-surv.ipw(pool.IPA,"q.M1.score")

M1.KM.fit.t2d<-survfit(Surv(dbd, db) ~ q.M1.score, weight=sbw, data=M1.KM )
M1.gg.t2d<-ggsurvplot(M1.KM.fit.t2d, fun="cumhaz", 
                      data = M1.KM, 
                      risk.table = F,
                      # risk.table.title="Baseline plasma IPA",
                      # break.time.by = 20,
                      censor=F,
                      ylim=c(0,0.5),
                      xlim=c(0,300), 
                      xlab="Months",
                      ylab="Cumulative incidence of T2D",
                      font.x = c(30, "bold", "black"),
                      font.y = c(30, "bold", "black"),
                      font.tickslab = c(30, "bold", "black"),
                      # pval="P trend = 0.002",
                      # pval.coord=c(320,0.95), 
                      legend.title="Metabolite composite score of sugar alcohols",legend.lab=c("Q1","Q2","Q3","Q4"),font.legend=c(28, "bold", "black"),
                      ggtheme = theme_classic() + 
                                theme(legend.direction = "horizontal"))
M1.gg.t2d$plot<-M1.gg.t2d$plot+guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

summary(coxph(Surv(dbd, db) ~ factor(q.M1.score),weight=sbw,ties="efron",data = M1.KM))

M1.KM.fit.chd<-survfit(Surv(tchd, chd) ~ q.M1.score, weight=sbw, data=M1.KM )
M1.gg.chd<-ggsurvplot(M1.KM.fit.chd, fun="cumhaz", 
                      data = M1.KM, 
                      risk.table = F,
                      # risk.table.title="Baseline plasma IPA",
                      # break.time.by = 20,
                      censor=F,
                      # ylim=c(0.75,1),
                      xlim=c(0,350), 
                      xlab="Months",
                      ylab="Cumulative incidence of CHD",
                      font.x = c(30, "bold", "black"),
                      font.y = c(30, "bold", "black"),
                      font.tickslab = c(30, "bold", "black"),
                      # pval="P trend = 0.002",
                      # pval.coord=c(320,0.95), 
                      legend.title="Metabolite composite score of sugar alcohols",legend.lab=c("Q1","Q2","Q3","Q4"),font.legend=c(28, "bold", "black"),
                      ggtheme = theme_classic() + 
                                theme(legend.direction = "horizontal"))
M1.gg.chd$plot<-M1.gg.chd$plot+guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

summary(coxph(Surv(tchd, chd) ~ factor(q.M1.score),weight=sbw,ties="efron",data = M1.KM))


M2.KM<-surv.ipw(pool.IPA,"q.M2.score")

M2.KM.fit.t2d<-survfit(Surv(dbd, db) ~ q.M2.score, weight=sbw, data=M2.KM )
M2.gg.t2d<-ggsurvplot(M2.KM.fit.t2d, fun="cumhaz", 
                      data = M2.KM, 
                      risk.table = F,
                      # risk.table.title="Baseline plasma IPA",
                      # break.time.by = 20,
                      censor=F,
                      ylim=c(0,0.5),
                      xlim=c(0,300), 
                      xlab="Months",
                      ylab="Cumulative incidence of T2D",
                      font.x = c(30, "bold", "black"),
                      font.y = c(30, "bold", "black"),
                      font.tickslab = c(30, "bold", "black"),
                      # pval="P trend = 0.002",
                      # pval.coord=c(320,0.95), 
                      legend.title="Metabolite composite score of synthetic sweeteners",legend.lab=c("Q1","Q2","Q3","Q4"),font.legend=c(28, "bold", "black"),
                      ggtheme = theme_classic() + 
                                theme(legend.direction = "horizontal"))
M2.gg.t2d$plot<-M2.gg.t2d$plot+guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

summary(coxph(Surv(dbd, db) ~ factor(q.M2.score),weight=sbw,ties="efron",data = M2.KM))

M2.KM.fit.chd<-survfit(Surv(tchd, chd) ~ q.M2.score, weight=sbw, data=M2.KM )
M2.gg.chd<-ggsurvplot(M2.KM.fit.chd, fun="cumhaz", 
                      data = M2.KM, 
                      risk.table = F,
                      # risk.table.title="Baseline plasma IPA",
                      # break.time.by = 20,
                      censor=F,
                      # ylim=c(0.75,1),
                      xlim=c(0,400), 
                      xlab="Months",
                      ylab="Cumulative incidence of CHD",
                      font.x = c(30, "bold", "black"),
                      font.y = c(30, "bold", "black"),
                      font.tickslab = c(30, "bold", "black"),
                      # pval="P trend = 0.002",
                      # pval.coord=c(320,0.95), 
                      legend.title="Metabolite composite score of synthetic sweeteners",legend.lab=c("Q1","Q2","Q3","Q4"),font.legend=c(28, "bold", "black"),
                      ggtheme = theme_classic() + 
                                theme(legend.direction = "horizontal"))
M2.gg.chd$plot<-M2.gg.chd$plot+guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

summary(coxph(Surv(tchd, chd) ~ factor(q.M2.score),weight=sbw,ties="efron",data = M2.KM))



options(bitmapType='cairo')
png(paste("~/sbw_KM.png",sep=''),res=400,width=8800, height=9000)
ggarrange(M1.gg.t2d$plot,M2.gg.t2d$plot,
          M1.gg.chd$plot,M2.gg.chd$plot,
          ncol=2,nrow=2)
dev.off()
















