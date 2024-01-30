#*****************************************************************************************
# Manuscript: "Social interactions do not  impact infection status in griffon vultures"
# Elvira D'Bastiani1*, Nili Anglister2, Inna Lysynyansky3, Inna Mikula3, Marta Acácio2,
# Gideon Vaadia2, Kaija Gahm1, Orr Spiegel2, Noa Pinter-Wollman1
# 1. Department of Ecology and Evolutionary Biology, University of California Los Angeles, 
# Los Angeles, California, USA
# 2. School of Zoology, Faculty of Life Sciences, Tel Aviv University, Tel Aviv, Israel.
# 3. Mycoplasma unit, Department of Avian Diseases, Professor Kimron Str 1, Kimron Veterinary Institute (KVI),
# POB12, Beit Dagan, Israel.
# *Corresponding author: Elvira D'Bastiani e-mail: elviradbastiani@gmail.com
#*****************************************************************************************

Citation: D’Bastiani, E.*, Anglister, N., Lysynyansky, I., Mikula, I., Acácio, M., Vaadia, G., Gahm, K., Spiegel, O., Pinter-Wollman, N. "Social Interactions and Mycoplasma Infection in Griffon Vultures." Journal of [Journal Name], [Year], [Volume(Issue)], [Page Range]. DOI: [DOI Number].

###############################################################################################################################
###############################################################################################################################
# Script to creating networks to co-roosting interactions
###############################################################################################################################
###############################################################################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
#setwd()

# Step 3: List files in the directory
dir()

# Step 4: Install and load necessary packages
#install.packages()
library(tidyverse) # for data wrangling
library(dplyr)#for data manipulation
library(lubridate) #for dates 
library(tidyverse) # for data wrangling
library(dplyr)#for data manipulation
library(ggplot2)#for plots
library(geosphere)# for spatial and geographic coordinates.
library(maptools) #for handling spatial objects
library(sf)
library(assortnet)
library(knitr)
library(devtools)
library(mapview)
library(ggpubr)

# Step 5: Loading the sampled data for mycoplasma sampling data of individuals used in 2021 #OBSERVATION: To test to 2022, simply change the year throughout the entire script.
load("output_infection_and_attributes_2021.rda")

# Step 6: See the dimension
dim(output_infection_and_attributes_2021)

# Step 7: Choose and sort the list_myDates only the sampling dates only from 2021 ex: "2021-09-14"
list_myDates <- sort(unique(as.Date(output_infection_and_attributes_2021$sample_event_date_edb)), decreasing = TRUE)
class(list_myDates)  # Displaying the class of list_myDates
list_myDates  # Displaying the content of list_myDates

# Step 8: Define the distThreshold in meters
distance <- c(25)  # You can also include c(10, 50, 75, 100)

# Step 9: Define the time lag
days <- c(14)  # You can also consider other values like c(7, 21, 28) 

for (i in 1:length(list_myDates)) {
  for (dist in 1:length(distance)) {
    for (da in 1:length(days)) {
      
  #***************************************************************************************************************************
  #SELECT THE WINDOW DAYS FROM MOVEBANK
  #***************************************************************************************************************************
  # Step 10: Select the data from sampling data
  load(file=paste0("Movebank_2021_day_",i,"_meters_",distance[dist],"_window_",days[da],".rda"))
  head(mydata_day)
  colnames(mydata_day)
  dim(mydata_day)
  sort(unique(mydata_day$id),F)
  table(mydata_day$id)
  length(unique(mydata_day$id))
  
  colnames(mydata_day)[6]<-"location_long"
  colnames(mydata_day)[7]<-"location_lat"
  
  #***************************************************************************************************************************
  #ANALYSES CO-ROOSTING TO DIRECT AND INDIRECT NETWORK
  #***************************************************************************************************************************
  # Step 11: Group the data by "IndividualName"
  grouped_data <- mydata_day %>%
    group_by(id)
  # Check the dimensions of the resulting dataset
  dim(grouped_data)
  
  # Step 12: Group the data by Individual and Date
  grouped_data <- mydata_day %>%
    group_by(id, date)
  # Check the dimensions of the resulting dataset
  dim(grouped_data)
  
  # Step 13: Sample one day for each individual
  selected_data <- grouped_data %>%
    distinct(id, .keep_all = TRUE)
  # Check the dimensions of the resulting dataset
  dim(selected_data)
  
  # Step 14: Count how many days each individual occurs
  day_counts <- selected_data %>%
    group_by(id) %>%
    summarise(DaysOccurred = n_distinct(date))
  # Check the dimensions of the resulting dataset
  dim(day_counts)
  day_counts
  
  # Step 15: Update the dataset with the selected data
  mydata_day <- selected_data
  dim(mydata_day)
  
  #---------------------------------------------------------------------------------
  # ---- CREATE THE DIRECT NETWORK TO ROOST INTERACTIONS ----
  #---------------------------------------------------------------------------------
  
  #For co-roosting direct run all_time_dist only once for each sampling time,
  #save it and load it for testing different windows rather than 
  #re-compute it every time because it will be very large and take very long to compute
  
  # Step 16: Create a SpatialPoints object with longitude and latitude coordinates
  coords_direct <- SpatialPoints(cbind(mydata_day$location_long, mydata_day$location_lat), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Step 17: Define the desired projection system for converting to meters (e.g., UTM) 
  # for the Israel zone IS36
  utm_proj_direct <- CRS("+proj=utm +zone=36 +datum=WGS84")
  
  # Step 18: Convert longitude and latitude to meters in UTM
  coords_utm_direct <- spTransform(coords_direct, utm_proj_direct)
  coords_utm_direct@coords
  
  # Step 19: Calculate all possible spatial distances between vultures
  xy_coords_direct <- coords_utm_direct@coords
  
  # Step 20: convert to a matrix format
  all_xy_dist_roosting_direct <- as.matrix(round(dist(xy_coords_direct)))
  
  # Step 21: Calculate all possible time distances between vultures based on the 'date' information
  all_time_dist_roosting_direct <- as.matrix(round(dist(mydata_day$date)))
  
  # Step 22: Select distances and times that fit specific conditions
  select_by_time_direct <- all_time_dist_roosting_direct <= 1  #Co-roosting considered less or equal a night
  # windows()
  # image(select_by_time_direct)
  select_by_dist_direct <- all_xy_dist_roosting_direct < distance
  # windows()
  # image(select_by_dist_direct)
  
  # Step 23: Combine time and distance thresholds
  comb_roosting_direct = select_by_time_direct * select_by_dist_direct
  diag(comb_roosting_direct) = 0  # Get rid of values on the diagonal (vulture can't interact with itself)
  ix = which(comb_roosting_direct > 0, arr.ind = TRUE)
  
  # Step 24: Find which vultures are interacting
  vul1 = mydata_day$id[ix[,1]]
  vul2 = mydata_day$id[ix[,2]]
  
  # Step 25: Create a data frame with interacting vultures
  edge_df_roosting_direct = as.data.frame(cbind(as.character(vul1), as.character(vul2)))
  net_roosting_direct<-graph_from_data_frame(edge_df_roosting_direct,directed = FALSE)
  
  # Step 26: Set up 'dummy' weights for edges
  E(net_roosting_direct)$weight = rep(1, length(E(net_roosting_direct)))
  
  # Step 27: Simplify the network, removing multiple edges and consolidating weights
  net_roosting_direct<-igraph::simplify(net_roosting_direct, remove.multiple = T,
                                        remove.loops = T,  edge.attr.comb = "sum")
  
  # Step 28: Access the edge weights of the simplified network
  E(net_roosting_direct)$weight
  direct_network_coroosting<-net_roosting_direct
  
  # Step 29: Write the simplified network to a file
  write_graph(direct_network_coroosting,file=paste0("direct_network_coroosting_2021_day_",
                                                    i,"_meters_",distance[dist],"_window_",
                                                    days[da],".txt"),format="ncol")
  
  #---------------------------------------------------------------------------------
  # ---- CREATE THE INDIRECT NETWORK TO ROOST INTERACTIONS ----
  #---------------------------------------------------------------------------------
  
  #For co-roosting indirect run all_time_dist only once for each sampling time,
  #save it and load it for testing different windows rather than 
  #re-compute it every time because it will be very large and take very long to compute
  
  # Step 30: Select only the distances and times that fit the conditions we are interested in:
  # interactions larger than 1 day so that we won't include direct interactions
  select_by_time_indirect = all_time_dist_roosting_direct > 1  # To co-roosting consider more than 1 day
  # windows()
  # image(select_by_time_indirect)
  
  # Step 31: distances smaller than the distance threshold so that we are choosing vultures that interacted with each other
  select_by_dist_indirect = all_xy_dist_roosting_direct < distance
  # windows()
  # image(select_by_dist_indirect)
  
  # Step 32: combine time and dist thresholds:
  comb_roosting_indirect = select_by_time_indirect * select_by_dist_indirect
  diag(comb_roosting_indirect) = 0  # Get rid of values on diagonal - a vulture can't interact with itself
  ix = which(comb_roosting_indirect > 0, arr.ind = TRUE)
  
  # Step 33: Find which vultures are interacting:
  vul1 = mydata_day$id[ix[,1]]
  vul2 = mydata_day$id[ix[,2]]
  
  # Step 34: Create a data frame with interacting vultures
  edge_df_roosting_indirect = as.data.frame(cbind(as.character(vul1), as.character(vul2)))
  net_roosting_indirect<-graph_from_data_frame(edge_df_roosting_indirect, directed = FALSE)
  
  # Step 35: Set up 'dummy' weights for edges
  E(net_roosting_indirect)$weight = rep(1, length(E(net_roosting_indirect)))
  
  # Step 36: Simplify the network, removing multiple edges and consolidating weights
  net_roosting_indirect<-igraph::simplify(net_roosting_indirect, remove.multiple = T,
                                          remove.loops = T,  edge.attr.comb = "sum")  # remove self loops -> remove.loops = T # consolidate multiple edges into one edge with a ->  edge.attr.comb = "sum"
  
  # Step 37: Access the edge weights of the simplified network
  E(net_roosting_indirect)$weight
  indirect_network_coroosting<-net_roosting_indirect
  
  # Step 38: Write the simplified network to a file
  write_graph(indirect_network_coroosting,file=paste0("indirect_network_coroosting_2021_day_",i,"_meters_",distance[dist],"_window_",days[da],".txt"),format="ncol")

    }
  }
}
