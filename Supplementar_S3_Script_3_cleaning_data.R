#*****************************************************************************************
# Manuscript: "Social interactions do not  impact infection status in griffon vultures"
# Elvira D’Bastiani1*, Nili Anglister2, Inna Lysynyansky3, Inna Mikula3, Marta Acácio2,
# Gideon Vaadia2, Kaija Gahm1, Orr Spiegel2, Noa Pinter-Wollman1
# 1.Department of Ecology and Evolutionary Biology, University of California Los Angeles, 
# Los Angeles, California, USA
# 2.School of Zoology, Faculty of Life Sciences, Tel Aviv University, Tel Aviv, Israel.
# 3.Department of Avian Diseases, Professor Kimron Str 1, Kimron Veterinary Institute (KVI),
# POB12, Beit Dagan, Israel.
# *Corresponding author: Elvira D’Bastiani e-mail: elviradbastiani@gmail.com
#*****************************************************************************************

############################################################################
#Script S3. Cleaning the spatial data
############################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
setwd()

# Step 3: List files in the directory
dir()

# Step 4: Load necessary packages
library(vultureUtils) # For vulture data cleaning and preliminary analyses
library(tidyverse)    # For data wrangling
library(move)         # For logging into Movebank
library(dplyr)        # For data manipulation
library(lubridate)    # For working with dates
library(geosphere)    # For spatial and geographic coordinates
library(maptools)     # For handling spatial objects
library(car)          # For Generalized Linear Models (GLM)
library(sf)           # For working with spatial data
library(devtools)     # For package development
library(mapview)      # For interactive maps in R

# Step 5: Load supporting data files
# CutOffRegion.kml is a polygon that roughly outlines Israel. We use this to exclude vulture points
# that fall outside of Israel (such as when the birds go on long-range forays).
mask <- sf::st_read("data/CutOffRegion.kml")
roostPolygons <- sf::st_read("data/new_roosts.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")

# Step 6: Load the dataset from movebank
load("final_data_movebank_2021.rda")

# Step 7: Display unique identifiers and their count
length(unique(final_data_movebank_2021$local_identifier))
sort(unique(final_data_movebank_2021$local_identifier), decreasing = TRUE)
length(sort(unique(final_data_movebank_2021$local_identifier), decreasing = TRUE))

# Step 8: Assign the loaded dataset to my_dataset
my_dataset <- final_data_movebank_2021

# Step 9: Convert specific columns to numeric
my_dataset <- my_dataset %>% mutate(across(c("location_long",       
                                             "location_lat",
                                             "external_temperature",
                                             "ground_speed",
                                             "barometric_height",
                                             "gps_time_to_fix",
                                             "heading",
                                             "gps_satellite_count",
                                             "height_above_msl"), as.numeric))

# Step 10: Define column names
colnames(my_dataset)[1] <- "idCol"
colnames(my_dataset)[4] <- "longCol"
colnames(my_dataset)[5] <- "latCol"

# Step 11: Argument checks
checkmate::assertClass(mask, "sf")
checkmate::assertDataFrame(my_dataset)
checkmate::assertChoice("gps_time_to_fix", names(my_dataset))
checkmate::assertChoice("heading", names(my_dataset))
checkmate::assertChoice("gps_satellite_count", names(my_dataset))
checkmate::assertChoice("ground_speed", names(my_dataset))
checkmate::assertChoice("external_temperature", names(my_dataset))
checkmate::assertChoice("barometric_height", names(my_dataset))
checkmate::assertClass(my_dataset$timestamp, "POSIXct")

# Step 12: For checking as we go along and getting a report: a little function to calculate rows, columns, and individuals.
getStats <- function(df) {
  rows <- nrow(df)
  cols <- ncol(df)
  indivs <- length(unique(df$idCol))
  return(c("rows" = rows, "cols" = cols, "indivs" = indivs))
  }
init <- getStats(my_dataset) # Get an initial baseline from the input data.

#-------------------------------------------------------------------------------
# Remove outliers based on speed
#-------------------------------------------------------------------------------

# Step 13: Remove outlier points based on zero
my_dataset <- my_dataset %>%
  dplyr::mutate(outlier = ifelse(.data$external_temperature == 0 & .data$barometric_height == 0 & .data$ground_speed == 0, 1, 0))

my_dataset <- my_dataset %>%
  dplyr::filter(is.na(outlier) | outlier == 0) %>%
  dplyr::select(-"outlier")
outliers <- getStats(my_dataset)

# Step 14: Filter out bad gps data
my_dataset <- my_dataset %>%
  dplyr::filter(.data$gps_time_to_fix <= 89)
badTimeToFix <- getStats(my_dataset)

# Step 15: Filter out bad heading data
my_dataset <- my_dataset %>%
  dplyr::filter(.data$heading <= 360 & .data$heading >= 0) # Only reasonable headings, between 0 and 360.
badHeading <- getStats(my_dataset) 

# Step 16: Only take locs that have at least 3 satellites
my_dataset <- my_dataset %>%
  dplyr::filter(.data$gps_satellite_count >= 3) # Must have at least 3 satellites in order to triangulate.
badSatellites <- getStats(my_dataset) 

# Step 17: Remove unrealistic "spiky" speeds
df <- vultureUtils::calcSpeeds(my_dataset, grpCol = "idCol", longCol = "longCol", latCol = "latCol")

# Step 18: Remove those that are for sure outliers: lead + lag > 180km/h
df2 <- df %>%
  dplyr::filter(lead_speed_m_s <= 50 & abs(lag_speed_m_s) <= 50) %>%
  dplyr::select(-c("lead_hour_diff_sec", "lead_dist_m", "lead_speed_m_s",
                   "lag_hour_diff_sec", "lag_dist_m", "lag_speed_m_s"))

# Step 19: Re-calculate speeds (because we removed some observations before)
df2 <- vultureUtils::calcSpeeds(df2, grpCol = "idCol", longCol = "longCol", latCol = "latCol")

# Step 20: Filter outliers based on lead_speed_m_s (observation -> Unfortunately the previous step did not get rid of all the outliers. 
#So we'll use only the lead to remove the other outliers)
df3 <- df2 %>%
  dplyr::filter(lead_speed_m_s <= 50)

# Step 21:see the maximum of lead_speed_m_s to confirm if to remove the outliers
max(df3$lead_speed_m_s)

# Step 22: This also does not get rid of all the outliers...But most of them are at night,
# which because of the reduced schedule, does not seem like such a large speed 
spikySpeeds <- getStats(df3)

# Step 23: So now we have to calculate if the fix is during the day or the night. 
# So, now determine if the data is during the day or night
times <- suncalc::getSunlightTimes(date = unique(lubridate::date(df3$timestamp)),
                                   lat = 31.434306, lon = 34.991889,
                                   keep = c("sunrise", "sunset")) %>%
  dplyr::select("dateOnly" = date, sunrise, sunset)

# Step 24: Manipulate and add daylight information to df3
df4 <- df3 %>%
  dplyr::mutate(dateOnly = lubridate::ymd(dateOnly)) %>%
  dplyr::left_join(times, by = "dateOnly") %>%
  dplyr::mutate(daylight = ifelse(timestamp >= sunrise & timestamp <= sunset, "day", "night")) %>%
  dplyr::select(-c(sunrise, sunset))

# Step 25: Manipulate and add daylight information to df3
df4 <- vultureUtils::calcSpeeds(df4, grpCol = "idCol", longCol = "longCol", latCol = "latCol")

# Step 26: Exclude if during the night the distance between two locations are 
# more than 10km apart, so filter out night outliers based on distance
df5 <- df4 %>%
  dplyr::mutate(
    day_diff = as.numeric(difftime(dplyr::lead(lubridate::date(timestamp)),
                                   lubridate::date(timestamp), units = "days")),
    night_outlier = ifelse(daylight == "night" &
                             day_diff %in% c(0, 1) &
                             dplyr::lead(daylight) == "night" &
                             lead_dist_m > 10000, T, F)
  ) %>%
  dplyr::filter(!night_outlier) %>%
  dplyr::select(-c("lead_hour_diff_sec", "lag_hour_diff_sec", "lead_dist_m",
                   "lag_dist_m", "lead_speed_m_s", "lag_speed_m_s"))

# Step 27: Calculate altitude "speeds" using vultureUtils package
nightDistance <- getStats(df5)

# Step 28: Create a copy of the original dataset
my_dataset <- df5

# Step 29: Remove unrealistic altitude values and calculate altitude "speeds"
dfAlt <- vultureUtils::calcSpeedsVert(df = my_dataset, grpCol = "idCol", altCol = "height_above_msl")

# Step 30: Rough filtering based on altitude speeds
dfAlt$height_above_msl[abs(dfAlt$lead_speed_m_s) > 2 |
                         abs(dfAlt$lag_speed_m_s) > 2] <- NA

dfAlt <- dfAlt %>%
  dplyr::select(-c("lead_hour_diff_sec", "lag_hour_diff_sec", "lead_dist_mV",
                   "lag_dist_mV", "lead_speed_m_s", "lag_speed_m_s"))

# Step 31: Count NA values after filtering
nAltitudesToNA <- sum(is.na(dfAlt$height_above_msl)) - sum(is.na(my_dataset$height_above_msl))
my_dataset <- dfAlt

# Step 32: Loading the data and calculating statistics
n_my_dataset <- getStats(my_dataset)
rawData <- my_dataset

# Step 33: Data cleaning - removing fixes with temperature == 0, height == 0, etc.
rawData$timestamp <- as.POSIXct(rawData$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
rawData$date <- as.Date(rawData$timestamp)

#-------------------------------------------------------------------------------
# Capture dates
#-------------------------------------------------------------------------------

# Step 34:  Classify roosts using get_roosts() function
roosts <- get_roosts(id = rawData$idCol, timestamp = as.POSIXct(rawData$timestamp),
                     x = rawData$longCol, y = rawData$latCol,
                     ground_speed = rawData$ground_speed, speed_units = "m/s",
                     buffer = 1, twilight = 61,
                     morning_hours = c(0:12), night_hours = c(13:23))

# Step 35: Parameters
start.day <- 01
start.month <- 01
end.day <- 31
end.month <- 12
distance <- 500 # distance is in meters, to calculate from the cage

# Step 36: Load the capture sites
sites <- read.csv("capture_sites.csv")

# Step 37: Subset the dataset for the start and end dates
sub.roosts <- roosts %>%
  mutate(start_date = as.Date(paste(start.day, start.month, year(date), sep="-"),
                              format="%d-%m-%Y"),
         end_date =as.Date(paste(end.day, end.month, year(date), sep="-"),
                           format="%d-%m-%Y")) %>%
  filter(date>=start_date&date<=end_date)

# Step 38: Print unique months in the subset
unique(month(sub.roosts$date))

# Step 39: Calculate roost distance to each capture cage and keep lines within 500m
crds <- matrix(c(sub.roosts$location_long, sub.roosts$location_lat),
               nrow = nrow(sub.roosts), ncol = 2)

DistanceMat <- matrix(ncol = nrow(sites), nrow = nrow(crds))
colnames(DistanceMat) <- unique(sites$name)

# Step 40: Calculate distances and print progress
for (i in 1:nrow(crds)) {
  DistanceMat[i,] <- round(distm(crds[i,], sites[,c(3,2)]),2)
  print(i)
}

# Step 41: ID of the closest capture site and distance from it
ClosestCaptureSite<-colnames(DistanceMat)[apply(DistanceMat,1,which.min)] 
ClosestCaptureDist<-apply(DistanceMat,1,min) 

# Step 42: Add columns to subset indicating capture status
sub.roosts <- cbind(sub.roosts, ClosestCaptureSite, ClosestCaptureDist)
sub.roosts$Captured <- ifelse(sub.roosts$ClosestCaptureDist <= distance, "Yes", "No")

# Step 43: Get a list of dates when the bird was inside the cage
sub.captured <- subset(sub.roosts, Captured == "Yes") 

# Step 44: Subset columns for captured dates
sub.captured.dates <- subset(sub.captured,
                             select=c("id", "date", "ClosestCaptureSite",
                                      "ClosestCaptureDist", "Captured"))

# Step 45: This does not work for the Carmel, because the roosts are within 500m of the cage and very often
# the birds roost on top of the cage without actually being inside it. So for the Carmel captures,
# we will use another protocol: if the birds were sleeping within 50m of the cage and
# it was a capture day (or 3 days before the release day), we consider the birds were captured and we remove the data
sub.captured.no.carmel <- subset(sub.captured.dates, ClosestCaptureSite != "Carmel")
sub.captured.carmel <- subset(sub.captured.dates, ClosestCaptureSite == "Carmel")

# Step 46: Read Carmel capture dates and adjust format
AllCarmelDates <- read.csv("all_captures_carmel_2010-2021.csv")
AllCarmelDates$Date <- as.Date(AllCarmelDates$Date, format = "%d/%m/%Y")

# Step 47: Create data frames for dates before captured dates
AllCarmelDates.1 <- data.frame(Date = as.Date(paste(AllCarmelDates$Date-1)))
AllCarmelDates.2 <- data.frame(Date = as.Date(paste(AllCarmelDates$Date-2)))
AllCarmelDates.3 <- data.frame(Date = as.Date(paste(AllCarmelDates$Date-3)))

# Step 48: Combine all Carmel dates for comparison
AllCarmelDates.all <- rbind(AllCarmelDates, AllCarmelDates.1, AllCarmelDates.2, AllCarmelDates.3)

# Step 49: Update captured status for Carmel captures
sub.captured.carmel <- sub.captured.carmel %>%
  mutate(known_capture = ifelse(date %in% AllCarmelDates.all$Date, 1, 0),
         Captured = ifelse(known_capture == 1 & ClosestCaptureDist <= 50, "yes", "no")) %>%
  filter(Captured == "yes") %>%
  dplyr::select(-c(known_capture))

names(sub.captured.no.carmel)
names(sub.captured.carmel)

# Step 50: Combine captured dates excluding Carmel and Carmel captured dates
sub.captured.dates <- rbind(sub.captured.no.carmel, sub.captured.carmel)

# Step 51: Exclude the day after the bird was captured
sub.captured.dates.1 <- sub.captured.dates
sub.captured.dates.1$date <- sub.captured.dates.1$date+1

# Step 52: Combine the datasets
sub.captured.dates <- rbind(sub.captured.dates, sub.captured.dates.1)
sub.captured.dates <- sub.captured.dates %>%
                      distinct(id, date, .keep_all = T)

# Step 53: Subset the dataset to exclude the capture dates
rawData.noout.5<-rawData
colnames(rawData.noout.5)[8]<-"id"
captured.rawData <- merge(rawData.noout.5, sub.captured.dates,
                          by=c("id","date"),
                          all.x=T)

# Step 54: Filter out non-captured rows and select relevant columns
captured.rawData1 <- captured.rawData %>%
  filter (is.na(Captured)) %>%
  dplyr::select(-c("daylight","ClosestCaptureSite", "Captured"))

# Step 55: Load the movebank data without cage dates
movebank_dataset<-captured.rawData1

# Step 56: Output the count of unique IDs in the movebank dataset
length(unique(movebank_dataset$id))

#-------------------------------------------------------------------------------
# Restrict the dataset to southern individuals based on a cutoff latitude
#-------------------------------------------------------------------------------
# Step 57: Based on previous investigations for the 2022 breeding and non-breeding seasons, 
# have found that a good cutoff for southern vs. non-southern is 3550000 (in ITM)
# So transform to SF object, so we can get centroids
sf <- movebank_dataset %>%
  sf::st_as_sf(coords = c("longCol","latCol"), remove = F) %>%
  sf::st_set_crs("WGS84") %>%
  sf::st_transform(32636)

# Step 58: Get centroids to determine southern individuals
centroids <- sf %>% group_by(id) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_centroid()

# Step 59: Examine a histogram of centroid latitudes 
hist(st_coordinates(centroids[,2])) # looks like 3550000 is generally a good cutoff point here

# Step 60: Get southern individuals based on cutoff latitude
southernIndividuals <- centroids %>%
  filter(st_coordinates(.)[,2] <= 3550000) %>% 
  pull(id)

# Step 61: Remove individuals from the dataset that are not in the south.
# Note that even if a "southern individual" has some points that
# are not in the south, the individual will still be kept and the
# points will also be kept.
Movebank_2021<- sf %>% filter(id %in% southernIndividuals)
#mapview(Movebank_2021)

# Step 62: Save to .rda file
save(Movebank_2021, file = "Movebank_2021.rda")

#**********************************************************************************************************
#end
#**********************************************************************************************************
