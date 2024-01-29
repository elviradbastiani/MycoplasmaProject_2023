#*****************************************************************************************
# Manuscript: "Social interactions do not  impact infection status in griffon vultures"
# Elvira D'Bastiani1*, Nili Anglister2, Inna Lysynyansky3, Inna Mikula3, Marta Acácio2,
# Gideon Vaadia2, Kaija Gahm1, Orr Spiegel2, Noa Pinter-Wollman1
# 1.Department of Ecology and Evolutionary Biology, University of California Los Angeles, 
# Los Angeles, California, USA
# 2.School of Zoology, Faculty of Life Sciences, Tel Aviv University, Tel Aviv, Israel.
# 3.Department of Avian Diseases, Professor Kimron Str 1, Kimron Veterinary Institute (KVI),
# POB12, Beit Dagan, Israel.
# *Corresponding author: Elvira D'Bastiani e-mail: elviradbastiani@gmail.com
#*****************************************************************************************

###############################################################################################################################
###############################################################################################################################
#S2 Script: Separating the window days of movement ecology data for griffon vultures prior to Mycoplasma spp. sampling in 2021.
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
load("output_infection_and_traits_2021.rda")

# Step 6: Load feed polygons
feedPolygons <- sf::st_read("roosts_to_feeding_situation.kml", quiet = TRUE) %>%
  sf::st_transform("WGS84")
feedPolygons <- feedPolygons %>% st_transform(32636)

# Step 7: Choose and sort the list_myDates only the sampling dates only from 2021 ex: "2021-09-14"
list_myDates <- sort(unique(as.Date(output_infection_and_traits_2021$sample_event_date_edb)), decreasing = TRUE)
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
  #ANALYSES CO-FEEDING TO DIRECT AND INDIRECT NETWORK
  #***************************************************************************************************************************
  # Step 11: Parameters
  feedBuffer=50
  consecThreshold=2
  speedThreshUpper=5
  speedThreshLower=NULL
  timeThreshold="10 minutes"
  daytimeOnly = T
  quiet = T
  includeAllVertices=T

  # Step 12: Filter points for speed <5 meters for feeding  (restrict interactions based on ground speed)
  filteredData <- vultureUtils::filterLocs(df = mydata_day,
                                           speedThreshUpper = 5,
                                           speedThreshLower = NULL, 
                                           speedCol = "ground_speed")
  
  # Step 13: Remove points that fall in or very near the feed site polygons,
  # by applying a buffer to those polygons. If feed polygons were provided, 
  # use them to filter out data
  if(!is.null(feedPolygons)){
    # Buffer the feed polygons
    if(!is.null(feedBuffer)){
      feedPolygons <- convertAndBuffer(feedPolygons, dist = feedBuffer)
    }
    # Exclude any points that fall within a (buffered) roost polygon
    points <- filteredData[lengths(sf::st_intersects(filteredData, feedPolygons)) == 0,]
  }else{
    message("No roost polygons provided; points will not be filtered by spatial intersection.")
    points <- filteredData
  }
  
  dim(points)
  head(points)
  
  # Step 14: Get sunrise and sunset times for each date in the mydata_day (restrict based on daylight)
  times <- suncalc::getSunlightTimes(date = unique(lubridate::date(points$timestamp)), 
                                     lat = 31.434306, lon = 34.991889,
                                     keep = c("sunrise", "sunset")) %>%
           dplyr::select(date, sunrise, sunset) 
          #The coordinates I'm using here are from the centroid of Israel calculated 
          #here: https://rona.sh/centroid.
    
  dim(times)
  
  # Step 15: Attach those sunrise and sunset times to the main mydata_day by left_join. 
  # Then, using the conditional case_when() function, create a new variable 
  # called daytime, which is set to TRUE for timestamps between sunrise and sunset, 
  # and FALSE otherwise.
  # remove leftover sunrise/sunset cols just in case
  points <- points %>%
    {if("sunrise" %in% names(.)) dplyr::select(., -sunrise) else .}%>%
    {if("sunset" %in% names(.)) dplyr::select(., -sunset) else .}%>%
    dplyr::left_join(times, by = c("dateOnly" = "date")) %>%
    dplyr::mutate(daytime = dplyr::case_when(timestamp > .data[["sunrise"]] &
                                               timestamp < .data[["sunset"]] ~ T, TRUE ~ F))
  dim(points)
  
  # Step 16: remove the nighttime points and filter out nighttimes
  points <- points %>%dplyr::filter(daytime == T)
  
  dim(points)
  
  # Step 17: Final mydata_day file
  mydata_1<- points
  
  # Step 18: Create a SpatialPoints object with the longitude and latitude coordinates
  coords <- SpatialPoints(cbind(mydata_1$location_long, mydata_1$location_lat), 
                          proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  ######!!!!! ##################################################################
  #The UTM (Universal Transverse Mercator) coordinate system
  #Within each zone, coordinates are measured as northings and eastings in meters
  ######!!!!! ##################################################################
  
  # Step 19: Define the desired projection system for converting to meters (e.g., UTM) 
  # Zone 36 is the Israel zone
  # WGS84 is a coordinate system, is the spatial reference system of GPS satellites with 
  #an error of less than 2 centimeters 
  utm_proj <- CRS("+proj=utm +zone=36 +datum=WGS84")
  
  # Step 20: Convert longitude and latitude to meters in UTM
  coords_utm <- spTransform(coords, utm_proj)
  coords_utm@coords
  
  dim(coords_utm)
  
  # Step 21: Get all possible space distances between all vultures in the data set
  xy_coords=coords_utm@coords
  all_xy_dist_feeding_direct= as.matrix(dist(xy_coords))
  #You can see the data here:
  #windows()
  #image(all_xy_dist_feeding_direct)
   
  # Step 22: Get all possible time distances between all vultures in the data set 
  # from the day information
  all_time_dist_feeding_direct = as.matrix(dist(mydata_1$timestamp))
  #You can see the data here:
  #windows()
  #image(all_time_dist_feeding_direct)
  
  
  #---------------------------------------------------------------------------------
  # ---- CREATE THE DIRECT NETWORK TO FEED INTERACTIONS ----
  #---------------------------------------------------------------------------------
  
  #For co-feeding direct run all_time_dist only once for each sampling time,
  #save it and load it for testing different windows rather than 
  #re-compute it every time because it will be very large and take very long to compute
  
  ## select only the distances and times that fit the conditions we are interested in:
  # Step 23:  Select the direct co-feeding interactions were recorded if vultures were co-feeding (i.e., within 25m of each other) within 0-30 minutes.
  select_by_time_feeding_direct = all_time_dist_feeding_direct<=1800 #30'*60" = 1800"
  
  # Step 24: Select the distances smaller than the distance threshold so that we are choosing vultures that interacted with each other
  select_by_dist_feeding_direct = all_xy_dist_feeding_direct<distance
  
  # Step 25: Combine time and distance thresholds:
  comb_feeding_direct = select_by_time_feeding_direct*select_by_dist_feeding_direct
  
  # Step 26: Get rid of values on diagonal (because a vulture can't interact with itself)
  diag(comb_feeding_direct)=0 
  ix = which(comb_feeding_direct>0, arr.ind = TRUE)
  
  # Step 27: Find which vultures are interacting:
  vul1 = mydata_1$id[ix[,1]]
  vul2 = mydata_1$id[ix[,2]]
  
  # Step 28: Creates a data frame combining the vultures that interact
  #`as.character(vul1)`: Converts the elements of the vector `vul1` to character type.
  #`as.character(vul2)`: Converts the elements of the vector `vul2` to character type.
  #`cbind(...)`: Combines `vul1` and `vul2` column-wise, creating a matrix.
  #`as.data.frame(...)`: Converts the matrix into a data frame, where each column corresponds to `vul1` and `vul2`.
  edge_df = as.data.frame(cbind(as.character(vul1), as.character(vul2)))
  
  # Step 29: Creates a graph object for the list of vultures that interact.
  net<-graph_from_data_frame(edge_df,directed = FALSE) 
  
  # Step 30: Accesses the edge weights. Creates a vector of 1s with the same length as the number of edges 
  # in the graph, setting the weight of each edge to 1.
  E(net)$weight = rep(1, length(E(net))) 
  
  #Step 31: Removes multiple edges between the same pair of vertices, and removes self-loops 
  #(edges where the source and target vertices are the same).
  #If there are multiple edges between the same pair of vertices, their weights are summed up.
  net_feeding_direct<-igraph::simplify(net, remove.multiple = T,
                                       remove.loops = T,  
                                       edge.attr.comb = "sum") # remove self loops -> remove.loops = T #consolidate multiple edges into one edge with a ->  edge.attr.comb = "sum"
  #Accesses the edge weights of the simplified graph `net_feeding_direct`.
  E(net_feeding_direct)$weight
  
  #Step 32: Creates a new graph object named 
  direct_network_cofeeding<-net_feeding_direct
  #You  can plot 
  #plot(direct_network_cofeeding, main=paste0("direct day - ", i))
  
  # Step 33: Writes the graph `direct_network_cofeeding` to a file in the ncol format with a specific filename based on the values of `i`, `distance[dist]`, and `days[da]`.
  write_graph(direct_network_cofeeding,file=paste0("direct_network_cofeeding_2021_day_",i,"_meters_",distance[dist],"_window_",days[da],".txt"),format="ncol")

  #---------------------------------------------------------------------------------
  # ---- CREATE THE INDIRECT NETWORK TO FEED INTERACTIONS ----
  #---------------------------------------------------------------------------------
  
  #For co-feeding indirect run all_time_dist only once for each sampling time,
  #save it and load it for testing different windows rather than 
  #re-compute it every time because it will be very large and take very long to compute
  
  # Step 34: Decide on the time lag
  time_lag =  14400#14400 # 4 hours in seconds # 4 horas*60 minutes = 240, 240 minutes*60 seconds = 14400
  
  # Step 35: Select only the distances and times that fit the conditions we are interested in:
  # Select the indirect co-feeding interactions were recorded if vultures were co-feeding (i.e., within 25m of each other) more than 4 hours
  select_by_time_feeding_indirect = all_time_dist_feeding_direct>=time_lag
  
  # Step 36: Distances smaller than the spatial distance threshold so that we are choosing vultures 
  #that interacted with each other
  select_by_dist_feeding_indirect = all_xy_dist_feeding_direct<distance
  
  # Step 37: combine time and dist thresholds:
  comb_feeding_indirect = select_by_time_feeding_indirect*select_by_dist_feeding_indirect
  
  # Step 38: get rid of values on diagonal - a vulture can't interact with itself
  diag(comb_feeding_indirect)=0 
  ix = which(comb_feeding_indirect>0, arr.ind = TRUE)
  
  # Step 39: Find which vultures are interacting:
  vul1 = mydata_1$id[ix[,1]]
  vul2 = mydata_1$id[ix[,2]]
  
  # Step 40: Creates a data frame combining the vultures that interact
  #`as.character(vul1)`: Converts the elements of the vector `vul1` to character type.
  #`as.character(vul2)`: Converts the elements of the vector `vul2` to character type.
  #`cbind(...)`: Combines `vul1` and `vul2` column-wise, creating a matrix.
  #`as.data.frame(...)`: Converts the matrix into a data frame, where each column corresponds to `vul1` and `vul2`.
  edge_df = as.data.frame(cbind(as.character(vul1), as.character(vul2)))
  
  # Step 41: Creates a graph object for the list of vultures that interact.
  net<-graph_from_data_frame(edge_df,directed = FALSE)
  
  # Step 42: Accesses the edge weights. Creates a vector of 1s with the same length as the number of edges 
  # in the graph, setting the weight of each edge to 1.
  E(net)$weight = rep(1, length(E(net))) # set up 'dummy' wights - giving each edge the wight of 1
  
  #Step 43: Removes multiple edges between the same pair of vertices, and removes self-loops 
  #(edges where the source and target vertices are the same).
  #If there are multiple edges between the same pair of vertices, their weights are summed up.
  net_feeding_indirect<-igraph::simplify(net, remove.multiple = T,
                                         remove.loops = T,  edge.attr.comb = "sum") # remove self loops -> remove.loops = T #consolidate multiple edges into one edge with a ->  edge.attr.comb = "sum"
  
  #Accesses the edge weights of the simplified graph `net_feeding_direct`.
  E(net_feeding_indirect)$weight
  
  #Step 44: Creates a new graph object named 
  indirect_network_cofeeding<-net_feeding_indirect
  #You  can plot 
  #plot(indirect_network_cofeeding, main=paste0("indirect day - ", i))
  
  # Step 33: Writes the graph `direct_network_cofeeding` to a file in the ncol format with a specific filename based on the values of `i`, `distance[dist]`, and `days[da]`.
  write_graph(indirect_network_cofeeding,file=paste0("indirect_network_cofeeding_2021_day_",i,"_meters_",distance[dist],"_window_",days[da],".txt"),format="ncol")
  
    }
  }
}
