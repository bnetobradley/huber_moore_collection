## Index Herbariorum stats on Herbarium Staff

## load libraries
library(jsonlite)
library(tidyverse)

## get the herbarium staff data
url <- "http://sweetgum.nybg.org/science/api/v1/staff"
df <- fromJSON(url)$data

## summarise the number of distinct staff codes per herbarium
df_sum <- df  %>% group_by(code) %>% summarise(n_distinct(irn))

## count the number of herbaria with 1 staff member
one_staff <- length(which(df_sum$`n_distinct(irn)` == 1))

## get the number of active herbaria data
url <- "http://sweetgum.nybg.org/science/api/v1/institutions"
df <- fromJSON(url)$data

n_herbaria <- df %>% filter(currentStatus == "Active") %>% summarise(n_distinct(irn))

## calculate percentage of herbaria with only one staff member
one_staff/n_herbaria