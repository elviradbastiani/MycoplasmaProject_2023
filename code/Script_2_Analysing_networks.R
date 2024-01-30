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

###############################################################################################################################
###############################################################################################################################
#Script to analyzing the social networks of griffon vultures
###############################################################################################################################
###############################################################################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
#setwd()

# Step 3: List files in the directory
dir()

# Step 4: Install and load necessary packages
library(igraph) # for network analysis


# Step 5: Loading the sampled data for mycoplasma sampling data of individuals used in 2021 #OBSERVATION: To test to 2022, simply change the year throughout the entire script.
load("output_infection_and_attributes_2021.rda")#mycoplasma_data
mycoplasma_sampling_date<-sort(unique(output_infection_and_attributes_2021$sample_event_date), F)
mycoplasma_sampling_date

#checking
output_infection_and_attributes_2021$edb_id
length(output_infection_and_attributes_2021$edb_id)

output_infection_and_attributes_2021$pos_neg
length(output_infection_and_attributes_2021$pos_neg)

output_infection_and_attributes_2021$sex
length(output_infection_and_attributes_2021$sex)

output_infection_and_attributes_2021$age
length(output_infection_and_attributes_2021$age)

# Step 6: Define the distThreshold in meters
distance = c(25)  # You can also include c(10, 25, 50, 75, 100)

# Step 7: Define the time lag
days = c(14)  # You can consider other values like c(7, 14, 21, 28) for comparison

# Step 8: Loop through all sampling dates, distances, and window sizes
for (i in 1:length(mycoplasma_sampling_date)) {
  for (dist in 1:length(distance)) {
    for (da in 1:length(days)) {
      
      # Step 9: Loop through all sampling dates, distances, and window sizes
      measures_results_day_direct<-data.frame(day=as.numeric(),
                                                  sampling_date=as.character(),
                                                  id=as.character(),
                                                  degree=as.numeric(),
                                                  betweenness=as.numeric(),
                                                  strength=as.numeric(),
                                                  infection=as.numeric(),
                                                  age=as.character(),
                                                  sex=as.character(),
                                                  color=as.character())
          
      # Step 10: Filter the mycoplasma data for the current sampling date
      output_infection_day <- output_infection_and_attributes_2021 %>%
                                  filter(sample_event_date == mycoplasma_sampling_date[i])
          
      # Step 11: Load the direct network roost for the current date, distance, and window
      direct_network<-igraph::read_graph(file=paste0("direct_network_cofeeding_2021_day_",i,"_meters_",distance[dist],"_window_",days[da],".txt"),format =  "ncol")
      #plot(direct_network)
      V(direct_network)
      E(direct_network)
      
      # Step 12: Calculate degree for each vulture in the direct network
      D_direct<-igraph::degree(direct_network, mode="total",normalized = T)
      
      # Step 13: Calculate betweenness for each vulture in the direct network
      B_direct<-igraph::betweenness(direct_network, directed=F, normalized = T)
      
      # Step 14: Calculate strength for each vulture in the direct network
      S_direct<-igraph::strength(direct_network, mode="total")
      
      # Step 15: Normalize population strength by individual number
      S_direct <- S_direct/sum(S_direct)
      
      #Individuals attributes
      # Step 16: Set default values
      infection_i<- NA
      age_i<- NA
      sex_i<- NA
      color_i<- "grey"
      cont_j=0
      
      # Step 17: Loop through all individuals in the direct network
      for(j in 1:length(V(direct_network)$name)) {
        name_direct_indj<- V(direct_network)$name[j]
        name_direct_indj
        
        D_direct_j<-D_direct[j]
        D_direct_j
        B_direct_j<-B_direct[j]
        B_direct_j
        S_direct_j<-S_direct[j]
        S_direct_j
        
        # Step 18: Retrieve individual attributes from the mycoplasma data
        attributes_day<-output_infection_day[output_infection_day$edb_id==name_direct_indj,]
        
        # Step 19: Determine infection status, age, sex, and color for the individual
        if (nrow(attributes_day)==0) {
          infection_i<- NA
          age_i<- NA
          sex_i<- NA
          color_i<- "grey"
        } else if (nrow(attributes_day)>=1) {
          
          if (unique(attributes_day$pos_neg)=="neg") {
            infection_i<- 0
            age_i<- unique(attributes_day$age)
            sex_i<- unique(attributes_day$sex)
            color_i<- "blue"
          } else if (unique(attributes_day$pos_neg)=="pos") {
            infection_i<- 1
            age_i<- unique(attributes_day$age)
            sex_i<- unique(attributes_day$sex)
            color_i<- "red"
          }
        }
        
        # Step 20: Update the results data frame with the calculated values
        cont_j=cont_j+1
        measures_results_day_direct[cont_j,]=c(i,
                                               paste0("cofeeding_direct_day_",i,"_meters_",distance[dist],"_window_",days[da]),
                                               name_direct_indj,
                                               D_direct_j,
                                               B_direct_j,
                                               S_direct_j,
                                               infection_i, 
                                               age_i, 
                                               sex_i, 
                                               color_i)
        # Step 21: Write the results to a file
        write.table(measures_results_day_direct, paste0("output_cofeeding_direct_day_",i,"_meters_",distance[dist],"_window_",days[da],".txt"))
        }
      }
    }
  }    


# ---- Indirect measures ----

# Step 22: Loop through all sampling dates, distances, and window sizes for indirect network
for (i in 1:length(mycoplasma_sampling_date)) { 
  for (dist in 1:length(distance)) {
    for (da in 1:length(days)) {
      
      # Step 23: Initialize a data frame to store the results for each day, distance, and window for indirect network
      measures_results_day_indirect <- data.frame(
        day = as.numeric(),
        sampling_date = as.character(),
        id = as.character(),
        degree = as.numeric(),
        betweenness = as.numeric(),
        strength = as.numeric(),
        infection = as.numeric(),
        age = as.character(),
        sex = as.character(),
        color = as.character()
      )
      
      # Step 24: Filter the mycoplasma data for the current sampling date
      output_infection_day <- output_infection_and_attributes_2021 %>%
        filter(sample_event_date == mycoplasma_sampling_date[i])
      
      # Step 25: Load the indirect network roost for the current date, distance, and window
      indirect_network <- igraph::read_graph(file = paste0("indirect_network_cofeeding_2021_day_", i, "_meters_", distance[dist], "_window_", days[da], ".txt"), format = "ncol")
      # plot(indirect_network)
      V(indirect_network)
      E(indirect_network)
      
      # Step 26: Calculate degree for each vulture in the indirect network
      D_indirect <- igraph::degree(indirect_network, mode = "total", normalized = TRUE)
      
      # Step 27: Calculate betweenness for each vulture in the indirect network
      B_indirect <- igraph::betweenness(indirect_network, directed = FALSE, normalized = TRUE)
      
      # Step 28: Calculate strength for each vulture in the indirect network
      S_indirect <- igraph::strength(indirect_network, mode = "total")
      
      # Step 29:  Normalize population strength by individual number
      S_indirect <- S_indirect / sum(S_indirect)
      
      
      #Individuals attributes
      # Step 30: Set default values
      infection_i<- NA
      age_i<- NA
      sex_i<- NA
      color_i<- "grey"
      cont_j=0
      
      # Step 31: Loop through all individuals in the indirect network
      for (j in 1:length(V(indirect_network)$name)) {
        name_indirect_indj <- V(indirect_network)$name[j]
        name_indirect_indj
        
        D_indirect_j <- D_indirect[j]
        D_indirect_j
        B_indirect_j <- B_indirect[j]
        B_indirect_j
        S_indirect_j <- S_indirect[j]
        S_indirect_j
        
        # Step 32: Retrieve individual attributes from the mycoplasma data
        attributes_day <- output_infection_day[output_infection_day$edb_id == name_indirect_indj,]
        
        # Step 33: Determine infection status, age, sex, and color for the individual
        if (nrow(attributes_day) == 0) {
          infection_i <- NA
          age_i <- NA
          sex_i <- NA
          color_i <- "grey"
        } else if (nrow(attributes_day) >= 1) {
          if (unique(attributes_day$pos_neg) == "neg") {
            infection_i <- 0
            age_i <- unique(attributes_day$age)
            sex_i <- unique(attributes_day$sex)
            color_i <- "blue"
          } else if (unique(attributes_day$pos_neg) == "pos") {
            infection_i <- 1
            age_i <- unique(attributes_day$age)
            sex_i <- unique(attributes_day$sex)
            color_i <- "red"
          }
        }
        
        # Step 34: Update the results data frame with the calculated values
        cont_j = cont_j + 1
        measures_results_day_indirect[cont_j,] = c(
          i,
          paste0("cofeeding_indirect_day_", i, "_meters_", distance[dist], "_window_", days[da]),
          name_indirect_indj,
          D_indirect_j,
          B_indirect_j,
          S_indirect_j,
          infection_i, 
          age_i, 
          sex_i, 
          color_i
        )
        
        # Step 35: Write the results to a file
        write.table(measures_results_day_indirect, paste0("output_cofeeding_indirect_day_", i, "_meters_", distance[dist], "_window_", days[da], ".txt"))
        }
      }
    }
  }


################################################
# Direct Network Results
################################################

# Day 1
day_1_direct_meters_25_window_14 <- read.table("output_cofeeding_direct_day_1_meters_25_window_14.txt")
# Day 2
day_2_direct_meters_25_window_14 <- read.table("output_cofeeding_direct_day_2_meters_25_window_14.txt")
# Day 3
day_3_direct_meters_25_window_14 <- read.table("output_cofeeding_direct_day_3_meters_25_window_14.txt")
# Day 4
day_4_direct_meters_25_window_14 <- read.table("output_cofeeding_direct_day_4_meters_25_window_14.txt")
# Day 5
day_5_direct_meters_25_window_14 <- read.table("output_cofeeding_direct_day_5_meters_25_window_14.txt")

# Combine results for all days
output_results_direct_2021 <- data.frame(rbind(
                                        day_1_direct_meters_25_window_14,
                                        day_2_direct_meters_25_window_14,
                                        day_3_direct_meters_25_window_14,
                                        day_4_direct_meters_25_window_14,
                                        day_5_direct_meters_25_window_14
                                        ))

# Check the length of the 'id' column
length(output_results_direct_2021$id)

# Write combined results to a file
write.table(output_results_direct_2021, paste0("output_results_direct_2021_meters_25_window_14.txt"))

# Filter for individuals with infection status 0 or 1
output_infection_direct_2021 <- output_results_direct_2021 %>% filter(infection == 0 | infection == 1)

# Write filtered results to a file
write.table(output_infection_direct_2021, paste0("output_infection_direct_2021_meters_25_window_14.txt"))

################################################
# Indirect Network Results
################################################

# Day 1
day_1_indirect_meters_25_window_14 <- read.table("output_cofeeding_indirect_day_1_meters_25_window_14.txt")
# Day 2
day_2_indirect_meters_25_window_14 <- read.table("output_cofeeding_indirect_day_2_meters_25_window_14.txt")
# Day 3
day_3_indirect_meters_25_window_14 <- read.table("output_cofeeding_indirect_day_3_meters_25_window_14.txt")
# Day 4
day_4_indirect_meters_25_window_14 <- read.table("output_cofeeding_indirect_day_4_meters_25_window_14.txt")
# Day 5
day_5_indirect_meters_25_window_14 <- read.table("output_cofeeding_indirect_day_5_meters_25_window_14.txt")

# Combine results for all days
output_results_indirect_2021 <- data.frame(rbind(
                              day_1_indirect_meters_25_window_14,
                              day_2_indirect_meters_25_window_14,
                              day_3_indirect_meters_25_window_14,
                              day_4_indirect_meters_25_window_14,
                              day_5_indirect_meters_25_window_14
                              ))


# Check the length of the 'id' column
length(output_results_indirect_2021$id)

# Write combined results to a file
write.table(output_results_indirect_2021, paste0("output_results_indirect_2021_meters_25_window_14.txt"))

# Filter for individuals with infection status 0 or 1
output_infection_indirect_2021 <- output_results_indirect_2021 %>% filter(infection == 0 | infection == 1)

# Write filtered results to a file
write.table(output_infection_indirect_2021, paste0("output_infection_indirect_2021_meters_25_window_14.txt"))





















