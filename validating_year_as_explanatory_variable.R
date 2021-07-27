#checking dispersion in categories between 2018 and 2019

#a = sum(residuals(glm_vc_quasi, type = "deviance")^2)
#pchisq(a,lower = F, glm_vc_quasi$df.residual, lower = F)
#pchisq(summary(glm_vc_quasi_2019)$dispersion * glm_vc_quasi$df.residual, glm_vc_quasi$df.residual, lower = F)

#pchisq(summary(glm_sc_quasi_2019)$dispersion * glm_sc_quasi$df.residual, glm_sc_quasi$df.residual, lower = F)

#pchisq(summary(glm_seed_quasi_2019)$dispersion * glm_seed_quasi$df.residual, glm_seed_quasi$df.residual, lower = F)

#Plotting signals for 2018 and 2019, to check if there is overlap. 

#2018 signals
library(ggplot2)
library(dplyr)

vc_signals_2018 = p_results_vc[p_results_vc$adjusted_p < .05,]
sc_signals_2018 = p_results_sc[p_results_sc$adjusted_p < .05,]
seed_signals_2018 = p_results_seed[p_results_seed$adjusted_p < .05,]

#2019 signals

vc_signals_2019 = p_results_vc_2019[p_results_vc_2019$adjusted_p < .05,]
sc_signals_2019 = p_results_sc_2019[p_results_sc_2019$adjusted_p < .05,]
seed_signals_2019 = p_results_seed_2019[p_results_seed_2019$adjusted_p < .05,]

vc_signals_2019$allele = rownames(vc_signals_2019)
sc_signals_2019$allele = rownames(sc_signals_2019)
seed_signals_2019$allele = rownames(seed_signals_2019)


vc_signals_2018$allele = rownames(vc_signals_2018)
sc_signals_2018$allele = rownames(sc_signals_2018)
seed_signals_2018$allele = rownames(seed_signals_2018)

#2018 signals that don't overlap in 2019:
p_results_vc_2019$allele = rownames(p_results_vc_2019)
p_results_vc$allele = rownames(p_results_vc)

miss_signal_2018_vc = anti_join(vc_signals_2018, vc_signals_2019, by = "allele")
miss_signal_2019_vc = semi_join(p_results_vc_2019, miss_signal_2018_vc, by = "allele")




#VC
#2019
ggplot(vc_signals_2019, aes(x = log(adjusted_p), y = 1)) +
  geom_point() + xlim(-10,0) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)
#2018
ggplot(vc_signals_2018, aes(x = log(adjusted_p), y = 1)) +
  geom_point() +xlim(-10,0) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)

#SC
#2019
ggplot(sc_signals_2019, aes(x = log(adjusted_p), y = 1)) +
  geom_point() + xlim(-50,1) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)

#2018
ggplot(sc_signals_2018, aes(x = log(adjusted_p), y = 1)) +
  geom_point() + xlim(-50,1) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)

#SEEDLING ONLY
#2019
ggplot(seed_signals_2019, aes(x = log(adjusted_p), y = 1)) +
  geom_point() + xlim(-50,1) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)

#2018
ggplot(seed_signals_2018, aes(x = log(adjusted_p), y = 1)) +
  geom_point() + xlim(-50,1) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)


#Now adding in H data for 2019, and checking the same results...

male_df_2019_h <- read.table(file = "./data/raw_data/2019_male_seed_counts_h.tsv",
                           header = TRUE,
                           sep = '\t')

all_exp_2019_h <- male_df_2019_h
vc_high_2019_h <- male_df_2019_h[male_df_2019_h$expression_category == "vegetative_cell_high", ] 
sc_high_2019_h <- male_df_2019_h[male_df_2019_h$expression_category == "sperm_cell_high", ]
seedling_only_2019_h <- male_df_2019_h[male_df_2019_h$expression_category == "seedling_only", ]


### GLM for Vegetative Cell category ###
glm_vc_quasi_2019_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                           as.factor(allele),
                         data = vc_high_2019_h,
                         family = quasibinomial(link = "logit"))
summary(glm_vc_quasi_2019_h)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_vc_quasi_2019_h)
p_vc_2019_h = coef(summary(glm_vc_quasi_2019_h))[1,4] #gives p value of intercept
for (i in 2:nrow(coef(summary(glm_vc_quasi_2019_h)))){ #loops through all B coeffecients besides the intercept
  beta = coef(summary(glm_vc_quasi_2019_h))[1,1] + coef(summary(glm_vc_quasi_2019_h))[i,1] #beta intercept and currenta beta
  sigma = sqrt(sum(S[c(1,i),c(1,i)])) #standard deviation is the square root of variance given by variance matrix 
  p_vc_2019_h = c(p_vc_2019_h, 2*(1-pnorm(abs(beta/sigma)))) #keeps adding the p value of the new intercept, z test one sided
}
names(p_vc_2019_h) = sort(unique(vc_high_2019_h$allele)) #sorts it

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_vc_adj_2019_h = p.adjust(p_vc_2019_h,method = "BH")
p_results_vc_2019_h=data.frame(raw_p = p_vc_2019_h, adjusted_p = p_vc_adj_2019_h)
p_results_vc_2019_h
p_results_vc_2019_h = arrange(p_results_vc_2019_h, p_vc_adj_2019_h)
write.table(p_results_vc_2019_h,"clipboard",sep="\t")



### GLM for Sperm Cell category ###
glm_sc_quasi_2019_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                           as.factor(allele),
                         data = sc_high_2019_h,
                         family = quasibinomial(link = "logit"))
summary(glm_sc_quasi_2019_h)

# Getting p-values for each line using a quasi-likelihood test (comparing to 
# Mendelian, 50% transmission)
S = vcov(glm_sc_quasi_2019_h)
p_sc_2019_h = coef(summary(glm_sc_quasi_2019_h))[1,4]
for (i in 2:nrow(coef(summary(glm_sc_quasi_2019_h)))){
  beta = coef(summary(glm_sc_quasi_2019_h))[1,1] + coef(summary(glm_sc_quasi_2019_h))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  p_sc_2019_h = c(p_sc_2019_h, 2*(1-pnorm(abs(beta/sigma))))
}
names(p_sc_2019_h) = sort(unique(sc_high_2019_h$allele))

# Multiple testing correction for the p-values using Benjamini-Hochberg method
p_sc_adj_2019_h = p.adjust(p_sc_2019_h,method="BH")
p_results_sc_2019_h=data.frame(raw_p = p_sc_2019_h, adjusted_p = p_sc_adj_2019_h)
p_results_sc_2019_h
p_results_sc_2019_h = arrange(p_results_sc_2019_h, p_sc_adj_2019_h)
write.table(p_results_sc_2019_h,"clipboard",sep="\t")

#Only looking at VC and SC since Seedling only typically had no signals anyway.... (Look at seedling only later)
#Now we will limit it to a .05 threshold...

vc_signals_2019_h = p_results_vc_2019_h[p_results_vc_2019_h$adjusted_p < .05,]
sc_signals_2019_h = p_results_sc_2019_h[p_results_sc_2019_h$adjusted_p < .05,]
#seed_signals_2019_h = p_results_seed_2019_h[p_results_seed_2019_h$adjusted_p < .05,]

vc_signals_2019_h$allele = rownames(vc_signals_2019_h)
sc_signals_2019_h$allele = rownames(sc_signals_2019_h)
#seed_signals_2019_h$allele = rownames(seed_signals_2019_h)

#VC
#2019_h
ggplot(vc_signals_2019_h, aes(x = log(adjusted_p), y = 1)) +
  geom_point() +xlim(-20,0) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)

#SC
#2019_h
ggplot(sc_signals_2019_h, aes(x = log(adjusted_p), y = 1)) +
  geom_point() + xlim(-100,1) + ylim(.9999,1.0001) + geom_text(aes(label = allele), hjust=0, vjust=-1, size=3)


#2018 signals that don't overlap in 2019:
p_results_vc_2019_h$allele = rownames(p_results_vc_2019_h)
#p_results_vc$allele = rownames(p_results_vc), already done above. 

miss_signal_2018_vc_h = anti_join(vc_signals_2018, vc_signals_2019_h, by = "allele")
miss_signal_2019_vc_h = semi_join(p_results_vc_2019_h, miss_signal_2018_vc_h, by = "allele")

