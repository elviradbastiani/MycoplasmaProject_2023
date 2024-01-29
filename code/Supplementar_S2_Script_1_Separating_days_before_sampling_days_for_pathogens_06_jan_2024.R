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

######################################################################################################################
#Script S4: Separating the 14 days of movement ecology data for griffon vultures prior to sampling to 2021
######################################################################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
#setwd()

# Step 3: List files in the directory
dir()

# Step 4: Install and load necessary packages
#install.packages("lubridate")
library(tidyverse) # for data wrangling
library(dplyr)#for data manipulation
library(lubridate) #for dates 

# Step 5: Load the cleaned data from Movebank for 2021
load("Movebank_2021.rda")

dim(Movebank_2021)
unique(sort(Movebank_2021$id), F)
length(unique(Movebank_2021$id))

# Step 6: Load the sampled data and list general individuals used in 2021
mycoplasma_data <- read.csv2("myco_parasite_results_2021.csv", sep = ",")
dim(mycoplasma_data)
colnames(mycoplasma_data)
unique(sort(mycoplasma_data$edb_id), F)
length(unique(mycoplasma_data$edb_id))
mycoplasma_data$sex

# ---- Define the age and sex of the individuals ----

# Step 7: Define birth date of the individuals
birth_date_edb <- as.Date(paste0(mycoplasma_data$birth_year, "-03-01", sep = ""))
mycoplasma_data_1 <- cbind(mycoplasma_data, birth_date_edb, birth_date_edb)
colnames(mycoplasma_data_1)
colnames(mycoplasma_data_1)[67] <- "age"
colnames(mycoplasma_data_1)

# Step 8: Check the birth date with year, month, and day
mycoplasma_data_1$birth_date_edb

# Step 9: Convert the date that individuals were sampled to as.Date
mycoplasma_data_1$sample_event_date_edb <- as.Date(mycoplasma_data_1$sample_event_date_edb, origin = "%YY-%m-%d")

# Step 10: Define the individuals that have data for infection and the n
individuals <- unique(mycoplasma_data_1$edb_id)
length(unique(mycoplasma_data_1$edb_id))

# Step 11: Loop for define the age of vultures
for (i in 1:length(individuals)) {
  name<-individuals[i]
  name
  ind_select<-mycoplasma_data_1[mycoplasma_data_1$edb_id==name,]
  ind_select
  
  if (nrow(ind_select)>0) {
    age<-unique(as.numeric(round(floor(unique(ind_select$sample_event_date_edb)-unique(ind_select$birth_date_edb))/365.25)))
    age
    
    if (unique(age)<=4) {
      ind_select$age<-"Juvenile"
      print(paste0(unique(age)," age ",ind_select$age, i))
    } else if (unique(age)>4)  {
      ind_select$age<-"Adult"
      print(paste0(unique(age)," age ",ind_select$age, i))
    } 
    
    if (i==1) {
      saida=ind_select
    } else {
      saidab=ind_select
      saida=rbind(saida,saidab)
    }
    
  }
}

output_infection_and_traits_2021 <- saida

# Step 12: Save rda file
save(output_infection_and_traits_2021, file = "output_infection_and_traits_2021.rda")

#*********************************************************************************************************
#Selecting the GPS only for 14 days with the distance between 25 meters of the distance among individuals
#*********************************************************************************************************
# Step 13: Choose and sort the list_myDates only the sampling dates only from 2021 ex: "2021-09-14"
list_myDates<-sort(unique(as.Date(output_infection_and_traits_2021$sample_event_date_edb)), F)
class(list_myDates)
list_myDates

# Step 14: Define the spatial distance threshold of the vultures
distance=25 #You can also include - c(10, 25, 50, 75, 100) meters

# Step 15: Define the windows of movement of vultures before sampling date
window=14 #You can test to other values - c(7,14,21,28) days

# Step 16: Define days that vultures are in the cage to remove
days_cage=3

# Step 17: Separate the data for all windows and for all sampling dates.
for (j in 1:length(window)) {
  for (d in 1:length(distance)) {
    for (ii in 1:length(list_myDates)) {
      
      #***************************************************************************************************************************
      #SELECT THE PREWINDOW DAYS FROM MOVEBANK
      #***************************************************************************************************************************
      # Step 18: Select the date before infection sampling to analyse
      Date_choose<-as_date(list_myDates[ii])
      print(paste0("Date_choose= ",Date_choose, "; Window= ",window[j]," days; Distance= ", distance[d]," meters"))
      
      # Step 19: Define the time less prewindow days
      Prewindow<-(Date_choose - days(days_cage))-days(window[j])
      print(paste0("Prewindow= ",Prewindow, "; Window= ",window[j]," days; Distance= ", distance[d]," meters"))
      
      # Step 20: Define the time less cage date
      Precage<-(Date_choose - days(days_cage))
      print(paste0("Precage= ",Precage, "; Window= ",window[j]," days; Distance= ", distance[d]," meters"))
      
      # Step 21: Individuals in the network during the specific selected time.
      mydata_day<-filter(Movebank_2021, date>=Prewindow&date<Precage)
      mydata_day<-cbind(mydata_day, Date_choose)
      dim(mydata_day)
      table(mydata_day$id)
      unique(mydata_day$id)
      sort(unique(mydata_day$id), F)
      length(unique(mydata_day$id))
      print(sort(unique(mydata_day$id), F))
      print(length(unique(mydata_day$id)))
      
      # Step 22: convert the data to an sf object
      mydata_day <- st_as_sf(mydata_day, crs = "WGS84", coords = c("location_long", "location_lat"), remove = F)
      
      # Step 23: save the dataset for this day
      save(mydata_day, file = paste0("Movebank_2021_day_",ii,"_meters_",distance,"_window_",window,".rda"))
      
    }
  }
}


#*######################################################################################################################
#Script S4: Separating the 14 days of movement ecology data for griffon vultures prior to sampling to 2022
######################################################################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
#setwd()

# Step 3: List files in the directory
dir()

# Step 4: Load necessary packages

# Step 5: Load the cleaned data from Movebank for 2022
load("Movebank_2022.rda")
dim(Movebank_2022)
unique(sort(Movebank_2022$id), F)
length(unique(Movebank_2022$id))

# Step 6: Load the sampled data and list general individuals used in 2022
mycoplasma_data <- read.csv2("myco_parasite_results_2022.csv", sep = ",")
dim(mycoplasma_data)
colnames(mycoplasma_data)
unique(sort(mycoplasma_data$edb_id), F)
length(unique(mycoplasma_data$edb_id))
mycoplasma_data$sex

# ---- Define the age and sex of the individuals ----

# Step 7: Define birth date of the individuals
birth_date_edb <- as.Date(paste0(mycoplasma_data$birth_year, "-03-01", sep = ""))
mycoplasma_data_1 <- cbind(mycoplasma_data, birth_date_edb, birth_date_edb)
colnames(mycoplasma_data_1)
colnames(mycoplasma_data_1)[30] <- "age"
colnames(mycoplasma_data_1)

# Step 8: Check the birth date with year, month, and day
mycoplasma_data_1$birth_date_edb

# Step 9: Convert the date that individuals were sampled to as.Date
mycoplasma_data_1$sample_event_date_edb <- as.Date(mycoplasma_data_1$sample_event_date_edb, origin = "%YY-%m-%d")

# Step 10: Define the individuals that have data for infection and the n
individuals <- unique(mycoplasma_data_1$edb_id)
length(unique(mycoplasma_data_1$edb_id))

# Step 11: Loop for define the age of vultures
for (i in 1:length(individuals)) {
  name<-individuals[i]
  name
  ind_select<-mycoplasma_data_1[mycoplasma_data_1$edb_id==name,]
  ind_select
  
  if (nrow(ind_select)>0) {
    age<-unique(as.numeric(round(floor(unique(ind_select$sample_event_date_edb)-unique(ind_select$birth_date_edb))/365.25)))
    age
    
    if (unique(age)<=4) {
      ind_select$age<-"Juvenile"
      print(paste0(unique(age)," age ",ind_select$age, i))
    } else if (unique(age)>4)  {
      ind_select$age<-"Adult"
      print(paste0(unique(age)," age ",ind_select$age, i))
    } 
    
    if (i==1) {
      saida=ind_select
    } else {
      saidab=ind_select
      saida=rbind(saida,saidab)
    }
    
  }
}

output_infection_and_traits_2022 <- saida

# Step 12: Save to rda file
save(output_infection_and_traits_2022, file = "output_infection_and_traits_2022.rda")

#*********************************************************************************************************
#Selecting the GPS only for 14 days with the distance between 25 meters of the distance among individuals
#*********************************************************************************************************
# Step 13: Choose and sort the list_myDates only the sampling dates only from 2022 ex: "2022-09-14"
list_myDates<-sort(unique(as.Date(output_infection_and_traits_2022$sample_event_date_edb)), F)
class(list_myDates)
list_myDates

# Step 14: Define the spatial distance threshold of the vultures
distance=25 #You can also include - c(10, 25, 50, 75, 100) meters

# Step 15: Define the windows of movement of vultures before sampling date
window=14 #You can test to other values - c(7,14,21,28) days

# Step 16: Define days that vultures are in the cage to remove
days_cage=3

# Step 17: Separate the data for all windows and for all sampling dates.
for (j in 1:length(window)) {
  for (d in 1:length(distance)) {
    for (ii in 1:length(list_myDates)) {
      
      #***************************************************************************************************************************
      #SELECT THE PREWINDOW DAYS FROM MOVEBANK
      #***************************************************************************************************************************
      # Step 18: Select the date before infection sampling to analyse
      Date_choose<-as_date(list_myDates[ii])
      print(paste0("Date_choose= ",Date_choose, "; Window= ",window[j]," days; Distance= ", distance[d]," meters"))
      
      # Step 19: Define the time less prewindow days
      Prewindow<-(Date_choose - days(days_cage))-days(window[j])
      print(paste0("Prewindow= ",Prewindow, "; Window= ",window[j]," days; Distance= ", distance[d]," meters"))
      
      # Step 20: Define the time less cage date
      Precage<-(Date_choose - days(days_cage))
      print(paste0("Precage= ",Precage, "; Window= ",window[j]," days; Distance= ", distance[d]," meters"))
      
      # Step 21: Individuals in the network during the specific selected time.
      mydata_day<-filter(Movebank_2022, date>=Prewindow&date<Precage)
      mydata_day<-cbind(mydata_day, Date_choose)
      dim(mydata_day)
      table(mydata_day$id)
      unique(mydata_day$id)
      sort(unique(mydata_day$id), F)
      length(unique(mydata_day$id))
      print(sort(unique(mydata_day$id), F))
      print(length(unique(mydata_day$id)))
      
      # Step 22: convert the data to an sf object
      mydata_day <- st_as_sf(mydata_day, crs = "WGS84", coords = c("location_long", "location_lat"), remove = F)
      
      # Step 23: save the dataset for this day
      save(mydata_day, file = paste0("Movebank_2022_day_",ii,"_meters_",distance,"_window_",window,".rda"))
      
    }
  }
}

#**********************************************************************************************************
#end
#**********************************************************************************************************
