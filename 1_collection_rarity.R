## collection rarity
## january 3rd 2024

## load libraries
library(rgbif)
library(ggthemes)
library(scales)
library(tidyverse)

## load data
hm <- readRDS("tidy_hm_data.rds")

## sort out a list of taxon names & taxon keys
hm_taxa <- unique(hm$species)
hm_keys <- c(unique(hm$usageKey,unique(hm$acceptedUsageKey)))
here <- ifelse(is.na(hm$acceptedUsageKey), hm$usageKey, hm$acceptedUsageKey)
hm_check_keys <- unique(here)
nothere <- ifelse(is.na(hm$acceptedUsageKey), hm$acceptedUsageKey,hm$usageKey)
nothere <- unique(nothere)

## get all herbarium records for hm_taxa
#hm_gbif_data <- occ_data(taxonKey = hm_keys, basisOfRecord = "PRESERVED_SPECIMEN", limit = 10000)
#saveRDS(hm_gbif_data, "gbif_data_hm_03jan2024.rds") # saved January 7th 2024
hm_gbif_data <- readRDS("gbif_data_hm_03jan2024.rds")

## get a list of accepted plant species from canadensys
canadensys <- read.delim("DWC-bdd48e6f-7abe-4ece-b100-c8b84b356a34/taxon.txt")
can_species <- canadensys %>% filter(taxonRank == "species" | taxonRank == "variety" | taxonRank == "subspecies") %>% filter(taxonomicStatus == "accepted") %>% select(scientificName)

#get taxon keys for canadian plant species (canadensys list filtered for accepted species, varieties and subspecies)
#aa <- c()
#bb <- c()
#for (i in 1:length(can_species$scientificName)) {
#  aa <- name_backbone(can_species$scientificName[i])
#  aa$sciname <- can_species$scientificName[i]
#  bb <- bind_rows(bb,aa)
#}
#saveRDS(bb,"bb_file.rds")
bb <- readRDS("bb_file.rds")
bb_keys <- unique(if_else(is.na(bb$acceptedUsageKey),bb$usageKey,bb$acceptedUsageKey))
not_bb_keys <- if_else(is.na(bb$acceptedUsageKey),bb$acceptedUsageKey,bb$usageKey)
bb_keys <- bb_keys[-which(bb_keys %in% not_bb_keys)]

## get number of records in canada
#d <- c()
#e <- c()
#for (j in 1:length(bb_keys)) {
#  d <- as.data.frame(occ_count(taxonKey = bb_keys[j], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA"))
#  d$taxonKey <- bb_keys[j]
#  e <- bind_rows(e,d)
#}
#saveRDS(e,"e_file.rds")
e <- readRDS("e_file.rds")
e <- e[-which(! e$taxonKey %in% bb_keys),]
## GET NUMBER OF RECORDS FOR HM COLLECTED TAXA
#f <- c()
#g <- c()
#for (k in 1:length(hm_check_keys)) {
#  f <- as.data.frame(occ_count(taxonKey = hm_check_keys[k], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA"))
#  f$taxonkey <- hm_check_keys[k]
#  g <- bind_rows(g,f)
#}
#saveRDS(g, "g_file.rds")
g <- readRDS("g_file.rds")

# remove outliers at wrong taxon level (carex genus taxon ID & "NA")
e <-  e %>% filter(`occ_count(taxonKey = bb_keys[j], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")` < 10000)
g <- g[-which(g$taxonkey %in% nothere),]
g <- g %>% filter(! is.na(taxonkey)) 

# remove bryophytes from hm data
bryos <- c("2679715","4279498","4279479","2680238","7990504")
g <- g[-which(g$taxonkey %in% bryos),]

#calc mean # records per canadian species
l1 <- mean(e$`occ_count(taxonKey = bb_keys[j], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")`, na.rm = T)
sd(e$`occ_count(taxonKey = bb_keys[j], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")`)
l2 <- mean(g$`occ_count(taxonKey = hm_keys[k], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")`)
sd(g$`occ_count(taxonKey = hm_keys[k], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")`, na.rm = T)
l2 <- as.numeric(l2[1])

## plot of number of records per plant species
 ggplot() + geom_histogram(aes(e$`occ_count(taxonKey = bb_keys[j], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")`+1), fill = "#117733", alpha = 0.8) + geom_vline(aes(xintercept = l1), colour = "#117733",  size = 5, alpha = 0.5) + geom_histogram(aes(g$`occ_count(taxonKey = hm_keys[k], basisOfRecord = "PRESERVED_SPECIMEN", country = "CA")`+1), fill ="#023d54", alpha = 0.95) +  theme_wsj(color = "white")   + geom_vline(aes(xintercept = l2), colour = "#023d54",  size = 5, alpha = 0.5) + scale_y_continuous(position = 'right') + theme( panel.grid.major.x = element_line(), panel.grid.major.y = element_blank(), axis.title = element_text(size=12, family = "Helvetica", face = "bold")) + labs(y = "Number of Plant Species", x = "Specimens collected in Canada") + annotate(x = 9000, y= 350, geom = "text",fontface = 3,size = 2.5, label = "commonly\ncollected") + annotate(x = 2, y= 350, geom = "text",fontface = 3, size = 2.5, label = "rarely\ncollected") + annotate(x = l2, y= 180, geom = "text", label = "HM collection mean", color = "white") + annotate(x = l1, y= 175, geom = "text", label = "country level mean",color = "white") + scale_x_continuous(trans = c("log10", "reverse"), breaks=c(1,2,10,100,1000,10000),labels=c("0","1","10","100","1000","10000"))  +  coord_flip()

ggsave("F1_rarity_plot.jpg", width = 3.4, height = 7)
