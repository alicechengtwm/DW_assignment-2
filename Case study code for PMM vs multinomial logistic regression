
library(mice)
library(nnet) 
library(readr)

################################################################################
# Read in data
################################################################################

#the dataset is in github repo
effect1_dataset <- read_csv("synthetic_observed_dataset_with_missingness.csv") 

effect1.df <- data.frame(effect1_dataset)
effect1.df$smokhx <- as.factor(effect1.df$smokhx)

################################################################################
# Multiple imputation using MICE: Use PMM for smoking history.
################################################################################

meth <- make.method(effect1.df)
meth["smokhx"] <- "pmm"
nimp <- round(100 * nic(effect1.df)/nrow(effect1.df),0)

effect1.imp.pmm <- mice(effect1.df,
                        method=meth,
                        m=nimp,
                        maxit=20,
                        seed=12022025, 
                        printFlag=F)
remove(meth,nimp)

################################################################################
# Multiple imputation using MICE: use multinomial logistic regression for smoking history.
################################################################################

meth <- make.method(effect1.df)
nimp <- round(100 * nic(effect1.df)/nrow(effect1.df),0)

effect1.imp.mlr <- mice(effect1.df,
                        method=meth,
                        m=nimp,
                        maxit=20,
                        seed=12022025,
                        printFlag=F)
remove(meth,nimp)

################################################################################
# Fit logistic regression model in imputed datasets
# Adjust for smoking history and all other covariates.
################################################################################

outcome.model.imp.mlr <- with(effect1.imp.mlr, 
                              glm(death1yr ~ female + age + acpulmed + 
  cshock + diabetes + highbp + cva + dyslip + famhxcad + angina + cancer + 
  dementia + pud + prevmi + asthma + depres + perartdis + chf + hyprthyr + as +
  prevrevasc + sysbp + diasbp + hrtrate + resp + hgb + wbc + sod + pot + 
  glucose + urea + cr + smokhx, family="binomial"))

outcome.model.mlr <- summary(pool(outcome.model.imp.mlr),conf.int=TRUE)


outcome.model.imp.pmm <- with(effect1.imp.pmm, 
                              glm(death1yr ~ female + age + acpulmed + 
  cshock + diabetes + highbp + cva + dyslip + famhxcad + angina + cancer + 
  dementia + pud + prevmi + asthma + depres + perartdis + chf + hyprthyr + as + 
  prevrevasc + sysbp + diasbp + hrtrate + resp + hgb + wbc + sod + pot + 
  glucose + urea + cr + smokhx, family="binomial"))

outcome.model.pmm <- summary(pool(outcome.model.imp.pmm),conf.int=TRUE)

################################################################################
# Extract odds ratios and 95% CIs.
################################################################################

Variables <- outcome.model.mlr[,"term"]

OR.mlr <- round(exp(outcome.model.mlr[,"estimate"]),3)
CI.lower.mlr <- round(exp(outcome.model.mlr[,"conf.low"]),3)
CI.upper.mlr <- round(exp(outcome.model.mlr[,"conf.high"]),3)

OR.pmm <- round(exp(outcome.model.pmm[,"estimate"]),3)
CI.lower.pmm <- round(exp(outcome.model.pmm[,"conf.low"]),3)
CI.upper.pmm <- round(exp(outcome.model.pmm[,"conf.high"]),3)

MLR.summary <- paste0(OR.mlr," (",CI.lower.mlr,",",CI.upper.mlr,")")
PMM.summary <- paste0(OR.pmm," (",CI.lower.pmm,",",CI.upper.pmm,")")

Table1.df <- data.frame(Variables,MLR.summary,PMM.summary)

OUTPUT_DIR <- "case_study_output" #change to your output directory
write.csv(Table1.df, file.path(OUTPUT_DIR, "table1.csv"), row.names = FALSE)

