homepath = "~/NM_rev"
setwd(homepath)

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
# library(chanmetab)
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
library(glmnet)
library(zCompositions) # For compositional zero replacement
options(bitmapType='cairo')
###################################################################################
###################################################################################
lnames1 = load(file="~/Taxon_use_AVG.RData")

lnames2 = load(file="~/TAXON_ONLY.RData")

AS_var=c("acesk_","aspart_","sach_","sucral_","tag_","eryth_",
         "maltitol_","mani_",
         # "pini_",
         "sorb_",
         "xyli_")

taxon.avg<-taxon.avg[-215,]
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
########################
INT<-function(dta,X)
{
  
  for (i in 1:length(X)){
    
    dta[,paste(X[i],"_INT",sep="")]<-qnorm((rank(dta[,X[i]],na.last="keep")-0.5)/sum(!is.na(dta[,X[i]])))
    
  }
  
  return(dta)
  
}
taxon.avg<-INT(taxon.avg,paste(AS_var, "avgall"   ,sep=""))
########################################################################
eln.AS <- function(data, a, X, Y) {
  
  # Initialize output storage
  sp.unique <- NULL
  co_rr <- NULL
  tune.para <- NULL
  sp.id <- list()
  result <- list()
  
  grid <- expand.grid(
    alpha = a,
    lambda = 10^seq(-5,1,length=1000)
  )
  
  X_df <- as.data.frame(X)
  print(dim(X_df))
  
  for (i in seq_along(Y)) {
    y_var <- Y[i]
    
    # 1. Assemble dataset for current outcome (prevents data leakage)
    z <- X_df
    z$y <- data[[y_var]]
    
    # 2. Train model via caret
    model <- train(
      y ~ ., 
      data = z, 
      method = "glmnet",
      preProcess = c("center", "scale"),
      trControl = trainControl("LOOCV"),
      # tuneLength=25
      tuneGrid=grid
    )
    
    # 3. Store caret training results
    result[[y_var]] <- model$results                # Fixed: 'model$results' (plural)
    
    # 4. Extract best tuning parameters
    tune.para <- rbind(tune.para, model$bestTune)
    rownames(tune.para)[nrow(tune.para)] <- y_var
    
    # 5. Safely extract non-zero coefficients
    elnet_raw <- coef(model$finalModel, s = model$bestTune$lambda)
    coef_matrix <- as.matrix(elnet_raw)
    
    # Filter out Intercept and zero-weighted features
    non_zero_names <- rownames(coef_matrix)[coef_matrix[, 1] != 0]
    non_zero_names <- setdiff(non_zero_names, "(Intercept)")
    
    # 6. Calculate correlation if features were selected
    if (length(non_zero_names) > 0) {
      feat_coefs <- coef_matrix[non_zero_names, 1, drop = FALSE]
      coef.out <- as.data.frame(t(feat_coefs))
      rownames(coef.out) <- "Beta"
      
      # Multiply features from X (NOT data) by coefficients
      x_selected <- as.matrix(X_df[, non_zero_names, drop = FALSE])
      pred_y <- x_selected %*% t(coef.out)
      CORR <- cor(pred_y, data[[y_var]], use = "complete.obs")
      
      sp.unique <- c(sp.unique, non_zero_names)
      sp.id[[y_var]] <- coef.out
    } else {
      # Fallback for models where all features were penalized to zero
      CORR <- NA
      sp.id[[y_var]] <- "No features selected"
    }
    
    # Append correlation result
    corr_df <- matrix(CORR, nrow = 1, dimnames = list("r", y_var))
    co_rr <- cbind(co_rr, corr_df)
  }
  
  out <- list(
    tuning_para = tune.para,
    cor_score   = co_rr,
    sp_id       = sp.id,
    sp_unique   = unique(sp.unique),
    training    = result
  )
  
  return(out)
}

# 
# y_problem <- taxon.avg$sorb_avgall_INT # Replace with your target exposure
# 
# # 1. Fit LOOCV ELN for a mid-range alpha (e.g., 0.5)
# cv_loocv <- cv.glmnet(
#   x = scale(use_clr), 
#   y = y_problem, 
#   alpha = 0.5, 
#   grouped = FALSE,            # Important for LOOCV in glmnet!
#   nfolds = length(y_problem)   # LOOCV
# )
# 
# # 2. Plot the CV Error Curve
# plot(cv_loocv)


# ########################################################################
# clr_scale_avgall_INT  <-eln.AS(taxon.avg,  
#                                    use_clr, 
#                                    paste(AS_var, "avgall_INT"   ,sep=""))
# 
# save(clr_scale_avgall_INT,file="CLR_cv_5.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.1.3  <-eln.AS(taxon.avg,  
                                   0.1,
                                   
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.1.3,file="INT_alpha_0.1_0.3.RData")
################################################################################################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.2.3  <-eln.AS(taxon.avg,  
                                   0.2,
                                   
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.2.3,file="INT_alpha_0.2_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.3.3  <-eln.AS(taxon.avg,  
                                   0.3,
                                    
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.3.3,file="INT_alpha_0.3_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.4.3  <-eln.AS(taxon.avg,  
                                   0.4,
                                    
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.4.3,file="INT_alpha_0.4_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.5.3  <-eln.AS(taxon.avg,  
                                   0.5,
                                    
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.5.3,file="INT_alpha_0.5_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.6.3  <-eln.AS(taxon.avg,  
                                   0.6,
                                   
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.6.3,file="INT_alpha_0.6_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.7.3  <-eln.AS(taxon.avg,  
                                   0.7,
                                   
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.7.3,file="INT_alpha_0.7_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.8.3  <-eln.AS(taxon.avg,  
                                   0.8,
                                    
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.8.3,file="INT_alpha_0.8_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT.9.3  <-eln.AS(taxon.avg,  
                                   0.9,
                                    
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT.9.3,file="INT_alpha_0.9_0.3.RData")
########################################################################
start.time <- Sys.time()

clr_scale_avgall_INT1.3  <-eln.AS(taxon.avg,  
                                   1,
                                    
                                   use_clr, 
                                   paste(AS_var, "avgall_INT"   ,sep=""))

end.time <- Sys.time()
end.time - start.time

save(clr_scale_avgall_INT1.3,file="INT_alpha_1_0.3.RData")