library(tidyverse)

male_df_18_19_h <- read.table(file = "./data/raw_data/2018_2019_h_male_seed_counts.tsv",
                           header = TRUE,
                           sep = '\t')

vc_high_18_19_h<- male_df_18_19_h[male_df_18_19_h$expression_category == "vegetative_cell_high", ] #NOTE CHANGE: instead of "expression_category_a, "expression_category" is used"
sc_high_18_19_h_18_19_h <- male_df_18_19_h[male_df_18_19_h$expression_category == "sperm_cell_high", ]
seedling_only_18_19_h <- male_df_18_19_h[male_df_18_19_h$expression_category == "seedling_only", ]
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

#####################################################################################################################################
###### Generalized linear models (GLM) ######
# Here we create generalized linear models with a quasibinomial distribution 
# and a logit link function

### GLM for Vegetative Cell category ###
glm_vc_quasi_18_19_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                      as.factor(allele) + as.factor(year),
                    data = vc_high_18_19_h,
                    family = quasibinomial(link = "logit"))
summary(glm_vc_quasi_18_19_h)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)

S = vcov(glm_vc_quasi_18_19_h)
p_vc_18_19_h = coef(summary(glm_vc_quasi_18_19_h))[1,4]
for (i in 2:nrow(coef(summary(glm_vc_quasi_18_19_h)))){
  beta = coef(summary(glm_vc_quasi_18_19_h))[1,1] + coef(summary(glm_vc_quasi_18_19_h))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_vc_18_19_h = c(p_vc_18_19_h, 2*(1-pnorm(abs(beta/sigma))))
}
p_vc_18_19_h
names(p_vc_18_19_h) = sort(unique(vc_high_18_19_h$allele))

p_vc_18_19_h
names(p_vc_18_19_h)[length(p_vc_18_19_h)] = 2019 #added this as we needed to name last p value but we know it's 2019
p_vc_18_19_h
# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_vc_adj_18_19_h = p.adjust(p_vc_18_19_h,method = "BH")

p_results_vc_18_19_h=data.frame(raw_p = p_vc_18_19_h, adjusted_p = p_vc_adj_18_19_h)
p_results_vc_18_19_h
p_results_vc_18_19_h = arrange(p_results_vc_18_19_h, p_vc_adj_18_19_h)
write.table(p_results_vc_18_19_h,"clipboard",sep="\t")



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

### GLM for Sperm Cell category ###
glm_sc_quasi_18_19_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                      as.factor(allele) + as.factor(year),
                    data = sc_high_18_19_h,
                    family = quasibinomial(link = "logit"))
summary(glm_sc_quasi_18_19_h)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_sc_quasi_18_19_h)
p_sc_18_19_h = coef(summary(glm_sc_quasi_18_19_h))[1,4]
for (i in 2:nrow(coef(summary(glm_sc_quasi_18_19_h)))){
  beta = coef(summary(glm_sc_quasi_18_19_h))[1,1] + coef(summary(glm_sc_quasi_18_19_h))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_sc_18_19_h = c(p_sc_18_19_h, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_sc_18_19_h) = sort(unique(sc_high_18_19_h$allele))

p_sc_18_19_h
names(p_sc_18_19_h)[length(p_sc_18_19_h)] = 2019 #added this as we needed to name last p value but we know it's 2019
p_sc_18_19_h

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_sc_18_19_h_adj = p.adjust(p_sc_18_19_h,method="BH")
p_results_sc_18_19_h=data.frame(raw_p = p_sc_18_19_h, adjusted_p = p_sc_18_19_h_adj)
p_results_sc_18_19_h
p_results_sc = arrange(p_results_sc_18_19_h, p_sc_18_19_h_adj)
write.table(p_results_sc_18_19_h,"clipboard",sep="\t")


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

### GLM for Seedling Only category ###
glm_seed_quasi_18_19_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                        as.factor(allele) + as.factor(year),
                      data = seedling_only_18_19_h,
                      family = quasibinomial(link = "logit"))
summary(glm_seed_quasi_18_19_h)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_seed_quasi_18_19_h)
p_seed_18_19_h = coef(summary(glm_seed_quasi_18_19_h))[1,4]
for (i in 2:nrow(coef(summary(glm_seed_quasi_18_19_h)))){
  beta = coef(summary(glm_seed_quasi_18_19_h))[1,1] + coef(summary(glm_seed_quasi_18_19_h))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_seed_18_19_h = c(p_seed_18_19_h, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_seed_18_19_h) = sort(unique(seedling_only_18_19_h$allele))

p_seed_18_19_h
names(p_seed_18_19_h)[length(p_seed_18_19_h)] = 2019 #added this as we needed to name last p value but we know it's 2019
p_seed_18_19_h


# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_seed_18_19_h_adj = p.adjust(p_seed_18_19_h,method="BH")
p_results_seed_18_19_h=data.frame(raw_p = p_seed_18_19_h, adjusted_p = p_seed_18_19_h_adj)
p_results_seed_18_19_h
p_results_seed_18_19_h = arrange(p_results_seed_18_19_h, p_seed_18_19_h_adj)
write.table(p_results_seed_18_19_h,"clipboard",sep="\t")


#Checking which ones pass below .05

library(ggplot2)
library(dplyr)

vc_signals_2018 = p_results_vc[p_results_vc$adjusted_p < .05,]
sc_signals_2018 = p_results_sc[p_results_sc$adjusted_p < .05,]
seed_signals_2018 = p_results_seed[p_results_seed$adjusted_p < .05,]

#2019 signals

vc_signals_18_19_h = p_results_vc_18_19_h[p_results_vc_18_19_h$adjusted_p < .05,]
sc_signals_18_19_h = p_results_sc_18_19_h[p_results_sc_18_19_h$adjusted_p < .05,]
seed_signals_18_19_h = p_results_seed_18_19_h[p_results_seed_18_19_h$adjusted_p < .05,]

vc_signals_18_19_h$allele = rownames(vc_signals_18_19_h)
sc_signals_18_19_h$allele = rownames(sc_signals_18_19_h)
seed_signals_18_19_h$allele = rownames(seed_signals_18_19_h)




