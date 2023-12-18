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

######################################################################################
#Script S1. Download the data from movebank
######################################################################################

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
setwd()

# Step 3: List files in the directory
dir()

# Step 4: Install the package vultureUtils
devtools::install_github("kaijagahm/vultureUtils") 

# Step 5: Load packages
library(vultureUtils) #Refer: https://github.com/kaijagahm/vultureUtils 

#-------------------------------------------------------------------------------
# Load the Movebank data from INPA 2021
#-------------------------------------------------------------------------------
# Step 6: Read the CSV file
rawData_inpa_2021 <- read.csv("2021_inpa_2021.csv")

# Step 7: Display dimensions of the data
dim(rawData_inpa_2021)

# Step 8: Display column names of the data
colnames(rawData_inpa_2021)

# Step 9: Extract unique timestamps and individual identifiers
unique(rawData_inpa_2021$timestamp)
unique(rawData_inpa_2021$individual.local.identifier)

# Step 10: Extract and format timestamp components
timestamp_i <- format(as.POSIXct(rawData_inpa_2021$timestamp, "%Y-%m-%d %H:%M:%S", tz = "UTC"), format = "%Y-%m-%d %H:%M:%S")
hour <- format(as.POSIXct(rawData_inpa_2021$timestamp, "%H:%M:%S"), format = "%H:%M:%S")
date <- format(as.Date(rawData_inpa_2021$timestamp, "%Y-%m-%d"), format = "%Y-%m-%d")
year <- format(as.Date(rawData_inpa_2021$timestamp, "%Y"), format = "%Y")

# Step 11: Add formatted timestamp components to the data frame
rawData_inpa_2021 <- data.frame(cbind(rawData_inpa_2021, rawData_inpa_2021$timestamp, date, hour, year, rawData_inpa_2021$individual.local.identifier))

# Step 12: Display new dimensions of the data
dim(rawData_inpa_2021)

# Step 13: Display new column names of the data
colnames(rawData_inpa_2021)

# Step 14: Rename a specific column to "id"
colnames(rawData_inpa_2021)[86] <- "id"

# Step 15: Save the downloaded data (.Rda format)
save(Download_Data2021_inpa, file = "Download_Data2021_inpa.Rda")

#-------------------------------------------------------------------------------
# Loading the Movebank data from Ornitella 2021
#-------------------------------------------------------------------------------
# Step 16: Define file paths
pwFilePath <- "pw.Rda"#insert the movebankCredentials
downloadDataPath <- "Download_Data2021_ornitella.Rda"

# Step 17: Load Movebank password from the .Rda file
load(pwFilePath)
MB.LoginObject <- movebankLogin(username = 'insert_name_here', password = pw)

# Step 18: Download data in a specific date range
Download_Data2021_ornitella <- vultureUtils::downloadVultures(
  loginObject = MB.LoginObject,
  extraSensors = FALSE,
  removeDup = TRUE,
  dateTimeStartUTC = as.POSIXct("2021-01-01 00:00"),
  dateTimeEndUTC = as.POSIXct("2021-12-31 00:00"),
  addDateOnly = TRUE,
  dfConvert = TRUE,
  quiet = FALSE
  )

# Step 19: Save the downloaded data (.Rda format)
save(Download_Data2021_ornitella, file = downloadDataPath)

#**********************************************************************************************************
#end
#**********************************************************************************************************