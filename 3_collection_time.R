## load libraries
library(rgbif)
library(dplyr)

## skip to line 70 if not donwloading GBIF data

# Function to query GBIF for observations and specimens within a geographic range
get_gbif_data <- function(species_key, lat, lon, radius_km, start_date = "1800-01-01", end_date = Sys.Date()) {
  # Define GBIF query parameters
  radius_deg <- radius_km / 111.32  # Convert radius to degrees (approximation)
  
  # Query for human observations
  human_obs <- occ_search(
    taxonKey = species_key,
    basisOfRecord = "HUMAN_OBSERVATION",
    decimalLatitude = paste(lat - radius_deg, lat + radius_deg, sep = ","),
    decimalLongitude = paste(lon - radius_deg, lon + radius_deg, sep = ","),
    limit = 300,
    eventDate = paste(start_date, end_date, sep = ","),
    fields = c("gbifID", "decimalLatitude", "decimalLongitude", "species", "eventDate", "basisOfRecord")
  )$data
  
  # Query for preserved specimens
  preserved_specimens <- occ_search(
    taxonKey = species_key,
    basisOfRecord = "PRESERVED_SPECIMEN",
    decimalLatitude = paste(lat - radius_deg, lat + radius_deg, sep = ","),
    decimalLongitude = paste(lon - radius_deg, lon + radius_deg, sep = ","),
    limit = 300,
    eventDate = paste(start_date, end_date, sep = ","),
    fields = c("gbifID", "decimalLatitude", "decimalLongitude", "species", "eventDate", "basisOfRecord")
  )$data
  
  # Combine the results
  combined_data <- bind_rows(human_obs, preserved_specimens)
  return(combined_data)
}

# Loop over the close_neighbors_results to retrieve GBIF data for each record
results_with_gbif <- list()

for (i in 1:nrow(close_neighbors_results)) {
  ref_id <- close_neighbors_results$ReferenceID[i]
  species_key <- plants$key_to_run[plants$ID == ref_id]
  lat <- plants$Latitude[plants$ID == ref_id]
  lon <- plants$Longitude[plants$ID == ref_id]
  
  # Skip if coordinates are missing
  if (is.na(lat) || is.na(lon)) next
  
  cat("Querying GBIF for Reference ID:", ref_id, "\n")
  
  gbif_data <- get_gbif_data(
    species_key = species_key,
    lat = lat,
    lon = lon,
    radius_km = 1  # Adjust radius as needed (e.g., 1 km)
  )
  
  # Add the reference ID to the retrieved data
  gbif_data <- gbif_data %>% mutate(ReferenceID = ref_id)
  results_with_gbif[[length(results_with_gbif) + 1]] <- gbif_data
}

# Combine all results into a single data frame
all_gbif_results <- bind_rows(results_with_gbif)
# Save the results to a CSV file
write.csv(all_gbif_results, "gbif_observations_and_specimens.csv", row.names = FALSE)

## START HERE IF SKIPPING GBIF DOWNLOAD
## load data
hm <- readRDS("tidy_hm_data.rds")
all_gbif_results <- read.csv("gbif_observations_and_specimens.csv")

# sort dates
all_gbif_results$dia <- as.Date(all_gbif_results$eventDate)

#find HM records
temporalnewplot <- hm %>% filter(Accession.Number %in% all_gbif_results$ReferenceID)
temporalnewplot$dates <- as.Date(temporalnewplot$Date.Collected, tryFormats = c("%d-%b-%y", "%Y"))
temporalnewplot$dates[3] <- as.Date("2001-01-01")

ggplot() + geom_point(aes(x =all_gbif_results$species, y =all_gbif_results$dia, color = all_gbif_results$basisOfRecord)) + geom_point(aes(x = temporalnewplot$species,y=temporalnewplot$dates), shape = 15, size = 2, color = "#023D54") + coord_flip() + theme_wsj(color = "white") + scale_color_manual(values = c("#D17a22","#023D54")) + theme(panel.grid.major = element_line(color = "gray"), legend.title = element_blank())
#ggsave("F3_time_plot.pdf", width = 7.16, height = 3.5)

y1 <- all_gbif_results %>% group_by(species) %>% arrange(all_gbif_results$dia) %>% summarise(first(dia), last(dia))
y2 <- temporalnewplot %>% group_by(species) %>% arrange(temporalnewplot$dates) %>% summarise(last(dates))

colnames(y2) <- c("species","hm_date")
colnames(y1) <- c("species","first","last") 

yy <- inner_join(y1,y2,by = "species")
yy$diff_start <- yy$first-yy$hm_date
yy$diff_end <- yy$hm_date-yy$last
