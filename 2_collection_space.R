## load libraries
library(sp)        
library(geosphere) 
library(dplyr)     

## load data
hm <- readRDS("tidy_hm_data.rds")

## filter for georef'd specimens
plants <- hm %>% filter(!is.na(Latitude) & !is.na(Longitude)) %>% select(ID = `Accession.Number`, Latitude, Longitude, Species = canonicalName, usageKey, acceptedUsageKey)
plants$key_to_run <- ifelse(is.na(plants$acceptedUsageKey), plants$usageKey, plants$acceptedUsageKey)

## load herbarium records for taxa in hm collection
gbif_data <- readRDS("gbif_data_hm_03jan2024.rds")

## function to calculate the next closest specimen and track neighbors within 1km
calc_next_closest_by_id <- function(ref_plants, target_plants, distance_threshold = 1000) {
  results <- list()  
  close_neighbors <- list()  
  
  for (i in seq_len(nrow(ref_plants))) {
    # Get coordinates for the reference plant
    ref_point <- c(ref_plants$Longitude[i], ref_plants$Latitude[i])
    ref_id <- ref_plants$ID[i]
    ref_species <- ref_plants$key_to_run[i]
    sppp <- ref_plants$Species[i] 
    
    if (any(is.na(ref_point)) || length(ref_point) != 2) {
      # If coordinates are invalid, skip
      results[[i]] <- list(
        ReferenceID = ref_id,
        SecondNearestID = NA,
        SecondNearestDistance = NA,
        SecondNearestBearing = NA
      )
      next
    }
    
    # Filter target plants to match the species of the reference plant
    matching_targets <- target_plants 
    
    # Skip if no matching target plants are found
    if (nrow(matching_targets) == 0) {
      results[[i]] <- list(
        ReferenceID = ref_id,
        SecondNearestID = NA,
        SecondNearestDistance = NA,
        SecondNearestBearing = NA
      )
      next
    }
    
    # Calculate distances to all matching target plants
    distances <- distm(ref_point, matching_targets[, c("Longitude", "Latitude")], fun = distHaversine)
    
    # Check for neighbors within the threshold
    neighbors_within_threshold <- which(distances <= distance_threshold)
    if (length(neighbors_within_threshold) > 0) {
      close_neighbors[[length(close_neighbors) + 1]] <- list(
        ReferenceID = ref_id,
        NeighborID = matching_targets$ID[neighbors_within_threshold],
        NeighborDistance = distances[neighbors_within_threshold]
      )
    }
    
    # Sort distances and find the second closest neighbor
    sorted_idx <- order(distances)
    if (length(sorted_idx) < 1) {
      # If less than two neighbors, no second closest
      results[[i]] <- list(
        ReferenceID = ref_id,
        SecondNearestID = NA,
        SecondNearestDistance = NA,
        SecondNearestBearing = NA
      )
      next
    }
    
    second_closest_idx <- sorted_idx[1]
    second_nearest_distance <- distances[second_closest_idx]
    
    # Calculate bearing (from the second nearest target plant to the reference plant)
    target_point <- c(matching_targets$Longitude[second_closest_idx], matching_targets$Latitude[second_closest_idx])
    second_nearest_bearing <- bearing(target_point, ref_point)  # Reversed coordinates for bearing
    
    # Save results
    results[[i]] <- list(
      ReferenceID = ref_id,
      SecondNearestID = matching_targets$ID[second_closest_idx],
      SecondNearestDistance = second_nearest_distance,
      SecondNearestBearing = (second_nearest_bearing + 360) %% 360  # Normalize bearing to 0-360 degrees
    )
  }
  
  # Convert results and close neighbors lists to data frames
  results_df <- bind_rows(results)
  close_neighbors_df <- bind_rows(close_neighbors)
  
  return(list(results = results_df, close_neighbors = close_neighbors_df))
}

## create results dataframe
next_closest_results <- data.frame()
close_neighbors_results <- data.frame()

## process each unique reference ID
mini_lookup <- plants %>% select(key_to_run,Species)
plants <- plants[!duplicated(plants),]
unique_ids <- unique(plants$ID)

for (id in unique_ids) {
  cat("Processing specimen ID:", id, "\n")
  
  # Filter the current reference plant
  ref_plants <- plants %>% filter(ID == id)
  
  # Safeguard for empty or missing data
  if (nrow(ref_plants) == 0) {
    cat("  No data for specimen ID:", id, "- Skipping.\n")
    next
  }
  
  # Extract species of the current reference plant
  ref_species <- ref_plants$key_to_run[1]
  
  # Filter target plants to include only matching species
  target_plants <- as.data.frame(bind_rows(lapply(gbif_data, function(x) x$data))) %>%
    filter(!is.na(decimalLatitude) & !is.na(decimalLongitude) & acceptedTaxonKey == ref_species) %>% select(ID = gbifID, Latitude = decimalLatitude, Longitude = decimalLongitude, Species = species)
  
  # Skip if no valid target plants for the species
  if (nrow(target_plants) == 0) {
    cat("  No valid target plants for specimen ID:", id, "- Skipping.\n")
    next
  }
  
  # Calculate the next closest specimen and neighbors within the threshold
  id_results <- calc_next_closest_by_id(ref_plants, target_plants)
  
  # Append results
  next_closest_results <- bind_rows(next_closest_results, id_results$results)
  close_neighbors_results <- bind_rows(close_neighbors_results, id_results$close_neighbors)
}

# Save results
saveRDS(next_closest_results, "HM_next_closest_by_id_results_species_filtered.rds")
saveRDS(close_neighbors_results, "HM_close_neighbors_by_id_results_species_filtered.rds")


next_closest_results <- readRDS("HM_next_closest_by_id_results_species_filtered.rds")
# Create the circle plot with lines from the center to each point
ggplot(next_closest_results, aes(x = SecondNearestBearing, y = log(SecondNearestDistance, 10)))  +  
  #geom_text(aes(y = 3, x = 235, label = "1km"), angle = 310 ,size = 2, color = "gray") +  geom_text(aes(y = 4, x = 235, label = "10 km"), angle = 310, size = 2, color = "gray") +  geom_text(aes(y = 5, x = 235, label = "100 km"), angle = 310,size = 2, color = "gray") +  geom_text(aes(y = 6, x = 233, label = "1000 km"), angle = 310,size = 2, color = "gray") +
  geom_point(color = "#023d54",size = 2, shape = 20) + 
  geom_segment(aes(x = SecondNearestBearing, xend = SecondNearestBearing, y = 2, yend = log(SecondNearestDistance, 10)), color = "#023d54", alpha = 0.3,size = 1) +  scale_x_continuous(limits = c(0,360),breaks = c(0,45, 90, 135, 180, 225, 270, 315),  labels = c("N","NE", "E", "SE", "S", "SW", "W", "NW")) +   coord_polar() +   theme_wsj(color = "white") + theme(axis.text.y = element_blank(),  legend.position = "none", axis.line.x.bottom = element_blank(), panel.grid.major = element_line(color = "gray")) +ylim(c(2,6)) 
ggsave("F2_space_plot.png",  width = 3.4,  height = 3.4)

mean(next_closest_results$SecondNearestDistance,na.rm = T)
sd(next_closest_results$SecondNearestDistance,na.rm = T)

## calculate mean bearing
calculate_mean_bearing <- function(bearings) {
  # Convert degrees to radians
  radians <- bearings * pi / 180
  
  # Compute x and y components
  x <- cos(radians)
  y <- sin(radians)
  
  # Compute mean x and y
  mean_x <- mean(x, na.rm = TRUE)
  mean_y <- mean(y, na.rm = TRUE)
  
  # Compute mean bearing in radians
  mean_bearing_radians <- atan2(mean_y, mean_x)
  
  # Convert to degrees
  mean_bearing_degrees <- mean_bearing_radians * 180 / pi
  
  # Normalize to [0, 360)
  mean_bearing_degrees <- (mean_bearing_degrees + 360) %% 360
  
  return(mean_bearing_degrees)
}


bearings <- next_closest_results$SecondNearestBearing
mean_bearing <- calculate_mean_bearing(bearings)
cat("Mean Bearing:", mean_bearing, "degrees\n")
## looks north/north-northeast

tocheck <- plants[! (plants$ID %in% next_closest_results$ReferenceID), ]
