#explore the dataset

# June 23
# KD

library(tidyverse)
library(janitor)

df <- read_csv("data/mac-means.csv") %>% 
  clean_names()
head(df)
names(df)

#conversion efficiency has individual and population level
#resource growth rate has individual and population level
# carrying capacity is pop level

df %>% group_by(simple_parameter, general_taxon_grouping_3, heterotroph_or_autotroph) %>% tally() %>% view()
#most of the data for EA estimates comes from heterotrophs across the board, EXCEPT resource growth rate, which comes largely from algae 