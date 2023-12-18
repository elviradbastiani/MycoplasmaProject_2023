#*****************************************************************************************
#Manuscript: "Social interactions do not  impact infection status in griffon vultures"
# Elvira D’Bastiani1*, Nili Anglister2, Inna Lysynyansky3, Inna Mikula3, Marta Acácio2,
# Gideon Vaadia2, Kaija Gahm1, Orr Spiegel2, Noa Pinter-Wollman1
# 1.Department of Ecology and Evolutionary Biology, University of California Los Angeles, 
# Los Angeles, California, USA
# 2.School of Zoology, Faculty of Life Sciences, Tel Aviv University, Tel Aviv, Israel.
# 3.Department of Avian Diseases, Professor Kimron Str 1, Kimron Veterinary Institute (KVI),
# POB12, Beit Dagan, Israel.
# *Corresponding author: Elvira D’Bastiani e-mail: elviradbastiani@gmail.com
#*****************************************************************************************

#########################################################################################
#Script S2. Checking individual names from GPS data for Ornitella and Inpa from movebank
#########################################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
setwd()

# Step 3: List files in the directory
dir()

# Step 4: Load packages
library(vultureUtils) #Refer: https://github.com/kaijagahm/vultureUtils 

#########################################################################
# Read and organize the names in the datasets for the year 2021 and 2022
#########################################################################

# Step 5: Loading the dataset from Movebank of ornitella and inpa
load("Download_Data2021_ornitella.Rda")
unique(sort(Download_Data2021_ornitella$local_identifier), decreasing = FALSE)

load("Download_Data2021_inpa.Rda")
unique(sort(Download_Data2021_inpa$individual.local.identifier), decreasing = FALSE)

#-------------------------------------------------------------------------------
# Select individuals from movebank from ornitella
#-------------------------------------------------------------------------------

# Step 6: Createa new object
substitute_names_Data2021_ornitella <- Download_Data2021_ornitella

# Step 7: Loop through rows to substitute names in the ornitella dataset
for (i in 1:nrow(substitute_names_Data2021_ornitella)) {
  if (substitute_names_Data2021_ornitella$local_identifier[i] == "E03") {
    substitute_names_Data2021_ornitella$local_identifier[i] <- "E03w"
    print(substitute_names_Data2021_ornitella$local_identifier[i])
  }
}

# Step 8: Display unique identifiers and their length after substitution
unique(sort(substitute_names_Data2021_ornitella$local_identifier), decreasing = FALSE)
length(unique(substitute_names_Data2021_ornitella$local_identifier))

# Step 9: Update the ornitella dataset with substituted names
Download_Data2021_ornitella <- substitute_names_Data2021_ornitella

# Step 10: List from ornitella of individuals used in 2021
movebank_ornitella <- Download_Data2021_ornitella

# Step 11: Display the new dimensions of the ornitella dataset
dim(Download_Data2021_ornitella)


#-------------------------------------------------------------------------------
# Select individuals from Movebank for inpa in 2021
#-------------------------------------------------------------------------------
# Step 12: Selecting the inpa data and planning to change the names in inpa 2021
T13w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="Y69>T13W", ]
T13w$id<-"T13w"
T23w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="K62>T23 white", ]
T23w$id<-"T23w"
T23b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T23 B (J69>E81>S24>S73)", ]
T23b$id<-"T23b"
A16w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T20 White", ]
A16w$id<-"A16w"
E07w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T22 B", ]
E07w$id<-"E07w"

T08b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T08 B (A73>Y31)", ]
T08b$id<-"T08b"
T42w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T42 W (L05>A92>Y09)", ]
T42w$id<-"T42w"
T91b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T91 B (N25>H66>Y90)", ]
T91b$id<-"T91b"
T98w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="Y11>T98 W", ]
T98w$id<-"T98w"
A50w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A50W (Y40>S70)", ]
A50w$id<-"A50w"

A63w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A63 White", ]
A63w$id<-"A63w"
A64w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A64 White", ]
A64w$id<-"A64w"
A65w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A65 White", ]
A65w$id<-"A65w"
A66w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A66 White", ]
A66w$id<-"A66w"
T40w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T40 White (Y14)", ]
T40w$id<-"T40w"

A68w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A68 White", ]
A68w$id<-"A68w"
A90w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A90 White", ]
A90w$id<-"A90w"
A99w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="S94>A99W", ]
A99w$id<-"A99w"
C53y<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="C53 Yellow", ]
C53y$id<-"C53y"
J79w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J79 White", ]
J79w$id<-"J79w"

T51b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T51 B (E08>H86>K51>S19)", ]
T51b$id<-"T51b"
J22w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J22 White", ]
J22w$id<-"J22w"
J24w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J24 White", ]
J24w$id<-"J24w"
J25w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J25 White", ]
J25w$id<-"J25w"
J29w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J29 White", ]
J29w$id<-"J29w"

J50w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J50 White", ]
J50w$id<-"J50w"
J51w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J51 White", ]
J51w$id<-"J51w"
J52w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J52 White", ]
J52w$id<-"J52w"
J54w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J54 White", ]
J54w$id<-"J54w"
J55w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J55 White", ]
J55w$id<-"J55w"

J56w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J56 White", ]
J56w$id<-"J56w"
J70w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J70 White", ]
J70w$id<-"J70w"
J94w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J94 White", ]
J94w$id<-"J94w"
J99w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="J99 White", ]
J99w$id<-"J99w"
T02w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="A74>T02 white", ]
T02w$id<-"T02w"

T24b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T24 B (J63>Y 65>R47>P05>S44)", ]
T24b$id<-"T24b"
T35w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T35 White", ]
T35w$id<-"T35w"
T36w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T36 White", ]
T36w$id<-"T36w"
T41w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T41 W (L10>L25>L10>A84)", ]
T41w$id<-"T41w"
T44w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T44 W", ]
T44w$id<-"T44w"

T51w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T51 White", ]
T51w$id<-"T51w"
T54w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T54 White", ]
T54w$id<-"T54w"
T55w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T55 White", ]
T55w$id<-"T55w"
T58w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T58 White", ]
T58w$id<-"T58w"
T60w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="Y01>T60 W", ]
T60w$id<-"T60w"

T68w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T68 White", ]
T68w$id<-"T68w"
T72b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T72 B", ]
T72b$id<-"T72b"
T75b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T75 B", ]
T75b$id<-"T75b"
T83b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T83 B", ]
T83b$id<-"T83b"
T93w<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T93 White", ]
T93w$id<-"T93w"

T95b<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="T95 B (N45>Y 78>X46>E05>K76)", ]
T95b$id<-"T95b"
K39<-Download_Data2021_inpa[Download_Data2021_inpa$individual.local.identifier=="K39", ]
K39$id<-"K39"

# Step 13: create a data frame with 51 individuals
movebank_inpa <- data.frame(rbind(T13w,T23w,T23b,A16w,E07w,
                                  
                                  T08b,T42w,T91b,T98w,A50w,
                                  
                                  A63w,A64w,A65w,A66w,T40w,
                                  
                                  A68w,A90w,A99w,C53y,J79w,
                                  
                                  T51b,J22w,J24w,J25w,J29w,
                                  
                                  J50w,J51w,J52w,J54w,J55w,
                                  
                                  J56w,J70w,J94w,J99w,T02w,
                                  
                                  T24b,T35w,T36w,T41w,T44w,
                                  
                                  T51w,T54w,T55w,T58w,T60w,
                                  
                                  T68w,T72b,T75b,T83b,T93w,
                                  
                                  T95b,K39))

#-------------------------------------------------------------------------------
# Select individuals from Movebank for ornitella and inpa in 2021
#-------------------------------------------------------------------------------
# Step 14: Copy the movebank_ornitella dataset into a new object 
movebank_ornitella_2021 <- movebank_ornitella

# Step 15: Display dimensions of the data
dim(movebank_ornitella_2021)

# Step 16: Display the column names 
colnames(movebank_ornitella_2021)

# Step 17: Sort the 'local_identifier' column of movebank_ornitella_2021 in ascending order 
#and remove duplicates. This command sorts the data in ascending order 
#(decreasing = FALSE) and displays unique values.
unique(sort(movebank_ornitella_2021$local_identifier, decreasing = FALSE))

# Step 18: Calculate and display the number of unique values in the 'local_identifier' column 
length(unique(sort(movebank_ornitella_2021$local_identifier, decreasing = FALSE)))
       
# Step 19: Copy the movebank_inpa dataset into a new object.
movebank_inpa_2021 <- movebank_inpa

# Step 20: Display dimensions of the data
dim(movebank_ornitella_2021)

# Step 21: Display the column names 
colnames(movebank_inpa_2021)

# Step 22: Sort the 'id' column of movebank_inpa_2021 in ascending order and remove duplicates
# This command sorts the data in ascending order (decreasing = FALSE) and 
# displays unique values.
unique(sort(movebank_inpa_2021$id, decreasing = FALSE))

# Step 23: Calculate and display the number of unique values in the 'id' column of movebank_inpa_2021
length(unique(sort(movebank_inpa_2021$id, decreasing = FALSE)))
              
# Step 24: Clean and select specific columns for Ornitella dataset
movebank_ornitella_specific <- movebank_ornitella_2021 %>%
dplyr::select(c("trackId", 
                "dateOnly", 
                "timestamp",
                "location_long",
                "location_lat",
                "external_temperature",
                "ground_speed",
                "local_identifier",
                "barometric_height",
                "gps_time_to_fix",
                "heading",
                "gps_satellite_count",
                "height_above_msl"))

# Step 25: Clean and select specific columns for INPA dataset
movebank_inpa_specific <- movebank_inpa_2021 %>%
dplyr::select(c("individual.local.identifier", 
                "date", 
                "timestamp",
                "location.long",
                "location.lat",
                "external.temperature",
                "ground.speed",
                "id",
                "bar.barometric.height",
                "gps.time.to.fix",
                "heading",
                "gps.satellite.count",
                "height.above.msl"))

# Step 26: Rename columns for consistency
colnames(movebank_inpa_specific) <- c("trackId", 
                                      "dateOnly",
                                      "timestamp", 
                                      "location_long", 
                                      "location_lat", 
                                      "external_temperature", 
                                      "ground_speed", 
                                      "local_identifier", 
                                      "barometric_height", 
                                      "gps_time_to_fix", 
                                      "heading", 
                                      "gps_satellite_count", 
                                      "height_above_msl")

# Step 27: Combine data.frame objects
final_data_movebank_2021 <- rbind(movebank_ornitella_specific, movebank_inpa_specific)
class(final_data_movebank_2021)

# Step 28: Save to .rda file
save(final_data_movebank_2021, file = "final_data_movebank_2021.rda")

#**********************************************************************************************************
#end
#**********************************************************************************************************