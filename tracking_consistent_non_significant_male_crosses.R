library(dplyr)
library(tidyverse)
library(stats)
library(ggplot2)
library(xlsx)
#2018_h

male_df_2018_h <- read.table(file = "./data/raw_data/2018_h_male_seed_counts.tsv",
                      header = TRUE,
                      sep = '\t')

unique_male_alleles_18_h = data.frame(unique(male_df_2018_h[c("allele")]))
unique_male_alleles_18_h

unique_male_alleles_18_h$number_of_crosses = NA

counter = 0

length(male_df_2018_h$allele)

for (val in 1:length(unique_male_alleles_18_h$allele)) 
{
  for (val2 in 1:length(male_df_2018_h$allele)) 
  {
    if (unique_male_alleles_18_h$allele[val] == male_df_2018_h$allele[val2])
    {
      counter = counter + 1
    }
  }
  unique_male_alleles_18_h$number_of_crosses[val] = counter
  counter = 0

}

#Looking at the tables, we already know that 

ggplot(unique_male_alleles_18_h, aes(x = number_of_crosses)) + geom_histogram(bins = 8)

#2019_h


male_df_2019_h <- read.table(file = "./data/raw_data/2019_male_seed_counts.tsv",
                             header = TRUE,
                             sep = '\t')
unique_male_alleles_19_h = data.frame(unique(male_df_2019_h[c("allele")]))
unique_male_alleles_19_h

unique_male_alleles_19_h$number_of_crosses = NA

counter = 0

length(male_df_2019_h$allele)

for (val in 1:length(unique_male_alleles_19_h$allele)) 
{
  for (val2 in 1:length(male_df_2019_h$allele)) 
  {
    if (unique_male_alleles_19_h$allele[val] == male_df_2019_h$allele[val2])
    {
      counter = counter + 1
    }
  }
  unique_male_alleles_19_h$number_of_crosses[val] = counter
  counter = 0
  
}

#Looking at the tables, we already know that 

ggplot(unique_male_alleles_19_h, aes(x = number_of_crosses)) + geom_histogram(bins = 10)


#2018 + 2019 + h

#male_df_18_19_h already in environment...

unique_male_alleles_18_19_h = data.frame(unique(male_df_18_19_h[c("allele")]))
unique_male_alleles_18_19_h
length(unique_male_alleles_18_19_h$allele)

unique_male_alleles_18_19_h$number_of_crosses = NA

counter = 0

length(male_df_18_19_h$allele)

for (val in 1:length(unique_male_alleles_18_19_h$allele)) 
{
  for (val2 in 1:length(male_df_18_19_h$allele)) 
  {
    if (unique_male_alleles_18_19_h$allele[val] == male_df_18_19_h$allele[val2])
    {
      counter = counter + 1
    }
  }
  unique_male_alleles_18_19_h$number_of_crosses[val] = counter
  counter = 0
  
}

#Looking at the tables, we already know that 

print(unique_male_alleles_18_19_h)

ggplot(unique_male_alleles_18_19_h, aes(x = number_of_crosses)) + geom_histogram(bins = 15, color = "black", size = 1) + ggtitle("2018 + 2019 + h w/ signals")

#Now let's exclude all alleles that attained an adjusted p of < .05

unique_male_alleles_no_s = unique_male_alleles_18_19_h#

#filtering out vc signals...
for (val in 1:length(unique_male_alleles_no_s$allele))
{
  for (val2 in 1:length(vc_signals_18_19_h$allele))
  {
    if (unique_male_alleles_no_s$allele[val] == vc_signals_18_19_h$allele[val2])
    {
      print(unique_male_alleles_no_s$number_of_crosses[val])
      unique_male_alleles_no_s = unique_male_alleles_no_s[-c(val), ]
      break

    }
  }
  
}

#filtering out sc signals...
for (val in 1:length(unique_male_alleles_no_s$allele))
{
  for (val2 in 1:length(sc_signals_18_19_h$allele))
  {
    if (unique_male_alleles_no_s$allele[val] == sc_signals_18_19_h$allele[val2])
    {
      print(unique_male_alleles_no_s$number_of_crosses[val])
      unique_male_alleles_no_s = unique_male_alleles_no_s[-c(val), ]
      break
    }
  }
  
}


#no need to check seedling as they are all non signals. 

ggplot(unique_male_alleles_no_s, aes(x = number_of_crosses)) + geom_histogram(bins = 15, color = "black", size = 1) + ggtitle("2018 + 2019 + h no signals")

median(unique_male_alleles_no_s$number_of_crosses)

values_to_check = unique_male_alleles_no_s[unique_male_alleles_no_s$number_of_crosses > 9, ] #looking at a random threshold, looking above the median in this case

edge_p_vc = p_results_vc_18_19_h[ p_results_vc_18_19_h$adjusted_p < .25, ]
edge_p_vc = edge_p_vc[edge_p_vc$adjusted_p > .05, ]

edge_p_sc = p_results_sc_18_19_h[ p_results_sc_18_19_h$adjusted_p < .25, ]
edge_p_sc = edge_p_sc[edge_p_sc$adjusted_p > .05, ]

edge_p = union(edge_p_vc, edge_p_sc)

#filtering out edge cases...
#from edge_p, we remove the alleles...

values_to_check = values_to_check[values_to_check$allele != "R103E04",]
values_to_check = values_to_check[values_to_check$allele != "R108A02",]
values_to_check = values_to_check[values_to_check$allele != "R02A05",]
values_to_check = values_to_check[values_to_check$allele != "R96B12",]

values_to_check
length(values_to_check$allele)

a = length(p_results_vc_18_19_h$raw_p)
b = length(p_results_sc_18_19_h$raw_p)
c = length(p_results_seed_18_19_h$raw_p)
print(a)
print(b)
print(c)

list = row.names(p_results_vc_18_19_h)
row.names(p_results_vc_18_19_h) = NULL
p_results_vc_18_19_h = cbind(p_results_vc_18_19_h, list)

list = row.names(p_results_sc_18_19_h)
row.names(p_results_sc_18_19_h) = NULL
p_results_sc_18_19_h = cbind(p_results_sc_18_19_h, list)

list = row.names(p_results_seed_18_19_h)
row.names(p_results_seed_18_19_h) = NULL
p_results_seed_18_19_h = cbind(p_results_seed_18_19_h, list)

p1 = p_results_vc_18_19_h[p_results_vc_18_19_h$list != "2019",]
p2 = p_results_sc_18_19_h[p_results_sc_18_19_h$list != "2019",]
p3 = p_results_seed_18_19_h[p_results_seed_18_19_h$list != "2019",]


p_results_all = bind_rows(p1,p2,p3)
p_results_all = p_results_all[c("raw_p", "adjusted_p","list...3")]

colnames(p_results_all)[3] = "allele"


p_results_all_below_005 = p_results_all[p_results_all$raw_p < .05,]
p_results_all_below_010 = p_results_all[p_results_all$raw_p < .10,]
p_results_all_below_015 = p_results_all[p_results_all$raw_p < .15,]
p_results_all_below_020 = p_results_all[p_results_all$raw_p < .20,]
p_results_all_below_025 = p_results_all[p_results_all$raw_p < .25,]

p_results_all_above_025 = p_results_all[p_results_all$raw_p >= .25, ]
p_results_all_above_020 = p_results_all[p_results_all$raw_p >= .20, ]
p_results_all_above_015 = p_results_all[p_results_all$raw_p >= .15, ]
p_results_all_above_010 = p_results_all[p_results_all$raw_p >= .10, ]
p_results_all_above_005 = p_results_all[p_results_all$raw_p >= .05, ]

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}
