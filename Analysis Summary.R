#### R Code for Analysis Summary ####
# This code takes in a dataset, and then runs an analysis with the given dataset. It then copys the information, and you can paste it on excel.
# In the output, the rows are representing the outputs for a single allele. This includes the name, logit of the percentage,
# estimated transmission rate, lower and upper confidence interval, raw p value, adjusted p value, expression category, and number of crosses.
# The analysis done is a generalized linear model, with a logit function. We also took an extra binomial variation approach to the analysis. 
# More details provided on the analsis in "Analysis Summary" page on Onenote

#install.packages("tidyverse")
# ^ above line is commented out as my computer already has the library installed
# If not installed, remove "#" and run the line above
library(tidyverse)
# tidyverse is a highly useful library, and is used throughout the code.

#### Loading Data ####

male_df <- read.table(file = "./data/raw_data/2018_2019_2020.tsv",
                              header = TRUE,
                              sep = '\t')
# This simply saves the data from the file desired. Here it's very simply, if you add a new tsv file that you want analyzed, just change
# the naming in the line after "file =" and in the quotations "Analysis Summary", put the location of your new file.


vc_high<- male_df[male_df$expression_category == "vegetative_cell_high", ] 
sc_high <- male_df[male_df$expression_category == "sperm_cell_high", ]
seedling_only <- male_df[male_df$expression_category == "seedling_only", ]
#The lines above seperate the dataset into subsets based on the expression category.

#### Analysis ####

#Vegetative Cell:

rm(S)
rm(raw_p)
rm(beta_list)
rm(estimated_coeffecient)
rm(low_int)
rm(up_int)
rm(beta)
rm(sigma)
# These variables are being removed as they are used within each categories code to generate output
# I just thought it would be cleaner to clear the data, and then create it again. Although, this is not necessary as the code will work either way.

#GLM
glm_vc_quasi_18_19_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~ # using the ratio between kernels as the output
                              as.factor(allele), # allele as explanatory variable, factored and there are multiple crosses per allele
                            data = vc_high, # data from vc subset, changes depending on expression category
                            family = quasibinomial(link = "logit")) # this is where the logit function is specificed, and quasibinomial approach

#This analysis is done automatically by rstudio, however the analysis uses a maximum likelihood approach which assumes that the data is independant.
#We are not sure if it's independant, so it's best to assume that it's not and therefore a different approach was taken to calculate the p values.
#The alternate approach is to use standard deviation in place of standard error. The standard deviation is calculated by the variance.
#It uses the formula of the variance of 2 random variables
#Var(X + Y) = Var(X) + Var(Y) + 2Cov(X,Y)
#This is easily calcualted using the variance covariance matrix below.

summary(glm_vc_quasi_18_19_h)
# Viewing the output from the rstudio GLM
#END GLM

#Calculating p values by hand, because of our methodology:
S = vcov(glm_vc_quasi_18_19_h) #variance covariance matrix
# this puts the variance of alleles, and the covariance of allele interactions in a matrix where each row and column corresponds to an allele
raw_p = coef(summary(glm_vc_quasi_18_19_h))[1,4] #raw p value of intercept
beta_list = coef(summary(glm_vc_quasi_18_19_h))[1,1] #intercept coeffecient value
estimated_coeffecient = (exp(beta_list)/(1+exp(beta_list))) #intercept estimated percentage between kernels based on beta, inverse of the logit function
low_int = exp(beta_list - qnorm(.975)*coef(summary(glm_vc_quasi_18_19_h))[1,2])/(1+exp(beta_list - qnorm(.975)*coef(summary(glm_vc_quasi_18_19_h))[1,2])) #Getting the lower confidence interval of the estimate
up_int = exp(beta_list + qnorm(.975)*coef(summary(glm_vc_quasi_18_19_h))[1,2])/(1+exp(beta_list + qnorm(.975)*coef(summary(glm_vc_quasi_18_19_h))[1,2]))
# calculating the confidence interval using the estimated coeffecient, more details in notes
# It only appears very long because it's collecting the initial values, the for loop belows shows how it's done with a loop. 

#This adds the rest of the alleles to the data, then combines the vectors. 
for (i in 2:nrow(coef(summary(glm_vc_quasi_18_19_h)))){ #The 2 simple indicates our starting index, since we already took our initial values above we start at 2. 
  beta = coef(summary(glm_vc_quasi_18_19_h))[1,1] + coef(summary(glm_vc_quasi_18_19_h))[i,1] #Collecting Intercept + allele coeffecients to calculate transmission rate
  sigma = sqrt(sum(S[c(1,i),c(1,i)])) #This is calculating the standard deviation by using the variance covariance matrix and having the location of alleles with the index
  raw_p = c(raw_p, 2*(1-pnorm(abs(beta/sigma)))) #Calculated the raw p values using a 2 sided z test
  beta_list = c(beta_list, beta) #Collecting intercept + allele coeffecients for analysis output
  estimated_coeffecient = c(estimated_coeffecient, (exp(beta)/(1+exp(beta)))) #calculating the estimated coeffecient, and adding it to a vector
  low_int = c(low_int, exp(beta-qnorm(.975)*sigma)/(1+exp(beta-qnorm(.975)*sigma))) #calculating the lower conf interval, and adding it to a vector
  up_int = c(up_int, exp(beta+qnorm(.975)*sigma)/(1+exp(beta+qnorm(.975)*sigma))) #same as above but for upper interval
                          
}

Alleles = sort(unique(vc_high$allele)) #creating unique allele vector
adjusted_p = p.adjust(raw_p, method = "BH") #calculating adjusted p's from raw p's using BH method

names(raw_p) = sort(unique(vc_high$allele))
names(beta_list) = sort(unique(vc_high$allele))
names(estimated_coeffecient) = sort(unique(vc_high$allele))
names(Alleles) = sort(unique(vc_high$allele))
names(adjusted_p) = sort(unique(vc_high$allele))
names(low_int) = sort(unique(vc_high$allele))
names(up_int) = sort(unique(vc_high$allele))
#This helps sort the data correctly in the table with respect to the allele names

analysis_summary_vc = data.frame(allele = Alleles, logit_of_percentage = beta_list, estimated_tranmission_percentage = estimated_coeffecient, lower_confidence_interval = low_int, upper_confidence_interval = up_int, raw_p_value = raw_p, adj_p_value = adjusted_p)
analysis_summary_vc$expression_category = "Vegetative Cell"
#Combining all the organized vectors into one large data table for vegetative cell analysis output.
#This is done for each expression category and is combines at the very end to have all the analysis output
#Finished creating output for VC


###########################################################################################################################################################

#Sperm Cell:

rm(S)
rm(raw_p)
rm(beta_list)
rm(estimated_coeffecient)
rm(low_int)
rm(up_int)
rm(beta)
rm(sigma)

#GLM
glm_sc_quasi_18_19_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                              as.factor(allele),
                            data = sc_high,
                            family = quasibinomial(link = "logit"))
summary(glm_sc_quasi_18_19_h)
#END GLM

#Creating Sperm Cell data frame

S = vcov(glm_sc_quasi_18_19_h) #variance covariance matrix
raw_p = coef(summary(glm_sc_quasi_18_19_h))[1,4] #raw p value of intercept
beta_list = coef(summary(glm_sc_quasi_18_19_h))[1,1] #intercept coeffecient value
estimated_coeffecient = (exp(beta_list)/(1+exp(beta_list))) #estimated percentage between kernels based on beta
low_int = exp(beta_list - qnorm(.975)*coef(summary(glm_sc_quasi_18_19_h))[1,2])/(1+exp(beta_list - qnorm(.975)*coef(summary(glm_sc_quasi_18_19_h))[1,2])) #Getting the lower confidence interval of the estimate
up_int = exp(beta_list + qnorm(.975)*coef(summary(glm_sc_quasi_18_19_h))[1,2])/(1+exp(beta_list + qnorm(.975)*coef(summary(glm_sc_quasi_18_19_h))[1,2]))
#This adds the rest of the alleles to the data, then combines the vectors. 
for (i in 2:nrow(coef(summary(glm_sc_quasi_18_19_h)))){
  beta = coef(summary(glm_sc_quasi_18_19_h))[1,1] + coef(summary(glm_sc_quasi_18_19_h))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  raw_p = c(raw_p, 2*(1-pnorm(abs(beta/sigma))))
  beta_list = c(beta_list, beta)
  estimated_coeffecient = c(estimated_coeffecient, (exp(beta)/(1+exp(beta))))
  low_int = c(low_int, exp(beta-qnorm(.975)*sigma)/(1+exp(beta-qnorm(.975)*sigma)))
  up_int = c(up_int, exp(beta+qnorm(.975)*sigma)/(1+exp(beta+qnorm(.975)*sigma)))
  
}

Alleles = sort(unique(sc_high$allele)) #creating unique allele vector
adjusted_p = p.adjust(raw_p, method = "BH") #calculating adjusted p's from raw p's

names(raw_p) = sort(unique(sc_high$allele))
names(beta_list) = sort(unique(sc_high$allele))
names(estimated_coeffecient) = sort(unique(sc_high$allele))
names(Alleles) = sort(unique(sc_high$allele))
names(adjusted_p) = sort(unique(sc_high$allele))
names(low_int) = sort(unique(sc_high$allele))
names(up_int) = sort(unique(sc_high$allele))
analysis_summary_sc = data.frame(allele = Alleles, logit_of_percentage = beta_list, estimated_tranmission_percentage = estimated_coeffecient, lower_confidence_interval = low_int, upper_confidence_interval = up_int, raw_p_value = raw_p, adj_p_value = adjusted_p)
analysis_summary_sc$expression_category = "Sperm Cell"

write.table(analysis_summary_sc,"clipboard",sep="\t")

###########################################################################################################################################################

#Seedline Only:

rm(S)
rm(raw_p)
rm(beta_list)
rm(estimated_coeffecient)
rm(low_int)
rm(up_int)
rm(beta)
rm(sigma)

#GLM
glm_seedling_only_quasi_18_19_h <- glm(cbind(number_of_GFP_kernels, number_of_WT_kernels) ~
                              as.factor(allele),
                            data = seedling_only,
                            family = quasibinomial(link = "logit"))
summary(glm_seedling_only_quasi_18_19_h)
#END GLM

#Creating Sperm Cell data frame

S = vcov(glm_seedling_only_quasi_18_19_h) #variance covariance matrix
raw_p = coef(summary(glm_seedling_only_quasi_18_19_h))[1,4] #raw p value of intercept
beta_list = coef(summary(glm_seedling_only_quasi_18_19_h))[1,1] #intercept coeffecient value
estimated_coeffecient = (exp(beta_list)/(1+exp(beta_list))) #estimated percentage between kernels based on beta
low_int = exp(beta_list - qnorm(.975)*coef(summary(glm_seedling_only_quasi_18_19_h))[1,2])/(1+exp(beta_list - qnorm(.975)*coef(summary(glm_seedling_only_quasi_18_19_h))[1,2])) #Getting the lower confidence interval of the estimate
up_int = exp(beta_list + qnorm(.975)*coef(summary(glm_seedling_only_quasi_18_19_h))[1,2])/(1+exp(beta_list + qnorm(.975)*coef(summary(glm_seedling_only_quasi_18_19_h))[1,2]))
#This adds the rest of the alleles to the data, then combines the vectors. 
for (i in 2:nrow(coef(summary(glm_seedling_only_quasi_18_19_h)))){
  beta = coef(summary(glm_seedling_only_quasi_18_19_h))[1,1] + coef(summary(glm_seedling_only_quasi_18_19_h))[i,1]
  sigma = sqrt(sum(S[c(1,i),c(1,i)]))
  raw_p = c(raw_p, 2*(1-pnorm(abs(beta/sigma))))
  beta_list = c(beta_list, beta)
  estimated_coeffecient = c(estimated_coeffecient, (exp(beta)/(1+exp(beta))))
  low_int = c(low_int, exp(beta-qnorm(.975)*sigma)/(1+exp(beta-qnorm(.975)*sigma)))
  up_int = c(up_int, exp(beta+qnorm(.975)*sigma)/(1+exp(beta+qnorm(.975)*sigma)))
  
}

Alleles = sort(unique(seedling_only$allele)) #creating unique allele vector
adjusted_p = p.adjust(raw_p, method = "BH") #calculating adjusted p's from raw p's

names(raw_p) = sort(unique(seedling_only$allele))
names(beta_list) = sort(unique(seedling_only$allele))
names(estimated_coeffecient) = sort(unique(seedling_only$allele))
names(Alleles) = sort(unique(seedling_only$allele))
names(adjusted_p) = sort(unique(seedling_only$allele))
names(low_int) = sort(unique(seedling_only$allele))
names(up_int) = sort(unique(seedling_only$allele))
analysis_summary_seedling_only = data.frame(allele = Alleles, logit_of_percentage = beta_list, estimated_tranmission_percentage = estimated_coeffecient, lower_confidence_interval = low_int, upper_confidence_interval = up_int, raw_p_value = raw_p, adj_p_value = adjusted_p)
analysis_summary_seedling_only$expression_category = "Seedling Only"

#Full Table of results:

rm(analysis_summary) #clearing any previous saved data

analysis_summary = rbind(analysis_summary_sc,analysis_summary_seedling_only, analysis_summary_vc) #combining the expression category outputs

analysis_summary$number_of_crosses = 0
  
num_crosses = table(unlist(male_df$allele))
names(num_crosses) = sort(unique(analysis_summary$allele))

for (i in 1:length(analysis_summary$allele)){
  analysis_summary$number_of_crosses[i] = num_crosses[analysis_summary$allele[i]]
}
#Adding one last column, and that's the number of crosses. This is calculated by the simple program above. 


write.table(analysis_summary, "clipboard", sep = "\t") #Copies the table to your computer, you can now paste it anywhere.

#End of file