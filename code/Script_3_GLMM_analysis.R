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
#Script to analyzing the relationship between mycoplasma infection status and social position (degree, betweenness, and strength) 
#of griffon vultures during co-feeding and co-roosting.
###############################################################################################################################
###############################################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# CO-FEEDING SITUATION
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
#setwd()

# Step 3: List files in the directory
dir()

# Step 4: Install and load necessary packages
#install.packages()

# Step 5: Load required libraries
library(ggplot2) # is a powerful and flexible plotting system in R
library(tidyverse) # is a collection of R packages designed to work together seamlessly for data analysis
library(dplyr)     # is a package for data manipulation and transformation
library(ggpubr)    # is an extension of ggplot2 that facilitates the creation and customization of publication-ready plots
library(performance) # is used for assessing and comparing the performance of statistical models.
library(DHARMa) #is a package for checking the distributional assumptions of regression models
library(lme4) #is a package for fitting linear mixed-effects models (LMMs) and generalized linear mixed-effects models (GLMMs). It extends the functionality of the base R `lm()` function to handle hierarchical and nested data structures.
library(broom.mixed) # is an extension of the broom package, which provides a consistent approach to tidying model outputs in R. broom.mixed specifically extends this functionality to mixed-effects models (LMMs and GLMMs) fitted with the lme4 package.

#*##############################
# DIRECT AND INDIRECT NETWORK
###############################

# ---- read file ----
dir()
results_DI <- read.table("output_infection_cofeeding_direct.txt", sep = " ")

# Group by 'id' and select the first sampling date
results_direct <- results_DI %>%
  group_by(id) %>%
  arrange(sampling_date) %>%
  slice(1) %>%
  ungroup()

# Print the new data frame
dim(results_direct)
print(results_direct)

# ---- read file ----
results_IND <- read.table("output_infection_cofeeding_indirect.txt", sep = " ")

# Group by 'id' and select the first sampling date
results_indirect <- results_IND %>%
                  group_by(id) %>%
                  arrange(sampling_date) %>%
                  slice(1) %>%
                  ungroup()

# Print the new data frame
dim(results_indirect)
print(results_indirect)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## GLMM MODELS
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Results Direct
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Dataset
dim(results_direct)

# Fit additional models for comparison
mod1_direct <- glmer(infection ~ scale(degree) + scale(strength) +
                       scale(betweenness) + as.factor(age) + (1|sampling_date), data = results_direct, family = binomial())

# Checking model assumptions
check_model(mod1_direct)

# Summary
summary(mod1_direct)

# Tidy the model results
tidy_results <- tidy(mod1_direct, conf.int = TRUE)
tidy_results$term <- c("Intercept",
                       "Degree" ,               
                       "Strength" ,             
                       "Betweenness" ,          
                       "Age",
                       "Sd Intercept")

# Create a forest plot using ggplot2
plot_1 <- ggplot(tidy_results, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(color = factor(term)), size = 3) +
  geom_errorbarh(height = 0.1) +
  labs(title = "GLMM Model Results",
       x = "Estimate",
       y = "Predictor") +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,1.2,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

print(plot_1)

ggsave(paste0("plot_feeding_direct_Estimate.pdf"), dpi = 600,
       width = 35,  height = 22,  units = "cm")

#***********************************************************************
# Cofeeding direct******************************************************
#***********************************************************************

results_direct$betweenness <- as.numeric(results_direct$betweenness)
results_direct$degree <- as.numeric(results_direct$degree)
results_direct$strength <- as.numeric(results_direct$strength)

colnames(results_direct)

# Plot Predicted data and original data points for betweenness, degree and strength
r1_direct <- ggplot(results_direct, aes(x = betweenness, y = as.numeric(infection))) +
  labs(x = "Betweenness", y = "Infection status") +
  geom_point(fill = "#d30000ff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_direct$infection), 
                                max(results_direct$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#d30000ff", se = T, fill = "#d30000ff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,1.2,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r1_direct

r2_direct <- ggplot(results_direct, aes(x = degree, y = as.numeric(infection))) +
  labs(x = "Degree", y = " ") +
  geom_point(fill = "#d30000ff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_direct$infection), 
                                max(results_direct$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#d30000ff", se = T, fill = "#d30000ff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,.8,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r2_direct

r3_direct <- ggplot(results_direct, aes(x = strength, y = as.numeric(infection))) +
  labs(x = "Strength", y = " ") +
  geom_point(fill = "#d30000ff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_direct$infection), 
                                max(results_direct$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#d30000ff", se = T, fill = "#d30000ff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,.8,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r3_direct

plots <- ggarrange(r1_direct, r2_direct, r3_direct,
                 labels = c("a.", "b.", "c."), 
                 ncol = 3, 
                 nrow = 1,
                 font.label = list(size = 30, family = "Helvetica", color = "black"),
                 hjust = -0.5, 
                 common.legend = T)

annotate_figure(plots, left = text_grob(" ", face = "bold", 
                                       color = "black", family = "Helvetica", rot = 90, size = 40),
                top = text_grob(" ", color = "black", family = "Helvetica",
                                face = "bold", size = 41),
                bottom = text_grob(" ", color = "black", family = "Helvetica",
                                   face = "bold", size = 35))

ggsave(paste0("plot_feeding_direct_new.pdf"), dpi = 600,
       width = 65,  height = 22,  units = "cm")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Results Indirect
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Dataset
dim(results_indirect)

# Fit additional models for comparison
mod1_indirect <- glmer(infection ~ scale(degree) + scale(strength) +
                         scale(betweenness) + as.factor(age) + (1|sampling_date), data = results_indirect, family = binomial())

# Checking model assumptions
check_model(mod1_indirect)

# Summary
summary(mod1_indirect)

# Tidy the model results
tidy_results <- tidy(mod1_indirect, conf.int = TRUE)
tidy_results$term <- c("Intercept",
                       "Degree",
                       "Strength",
                       "Betweenness",
                       "Age",
                       "Sd Intercept")

# Create a forest plot using ggplot2
plot_1 <- ggplot(tidy_results, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(color = factor(term)), size = 3) +
  geom_errorbarh(height = 0.1) +
  labs(title = "GLMM Model Results",
       x = "Estimate",
       y = "Predictor") +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, 1.2, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

print(plot_1)

ggsave(paste0("plot_feeding_indirect_Estimate.pdf"), dpi = 600,
       width = 35,  height = 22,  units = "cm")

#***********************************************************************
#Cofeeding Indirect
#***********************************************************************
results_indirect$betweenness <- as.numeric(results_indirect$betweenness)
results_indirect$degree <- as.numeric(results_indirect$degree)
results_indirect$strength <- as.numeric(results_indirect$strength)
colnames(results_indirect)

# Plot Predicted data and original data points for degree
r1_indirect <- ggplot(results_indirect, aes(x = betweenness, y = as.numeric(infection))) +
  labs(x = "Betweenness", y = "Infection status") +
  geom_point(fill = "#fcbdbaff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_indirect$infection), 
                                max(results_indirect$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#fcbdbaff", se = T, fill = "#fcbdbaff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, 1.2, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r1_indirect

r2_indirect <- ggplot(results_indirect, aes(x = degree, y = as.numeric(infection))) +
  labs(x = "Degree", y = " ") +
  geom_point(fill = "#fcbdbaff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_indirect$infection), 
                                max(results_indirect$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#fcbdbaff", se = T, fill = "#fcbdbaff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, .8, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r2_indirect

r3_indirect <- ggplot(results_indirect, aes(x = strength, y = as.numeric(infection))) +
  labs(x = "Strength", y = " ") +
  geom_point(fill = "#fcbdbaff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_indirect$infection), 
                                max(results_indirect$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#fcbdbaff", se = T, fill = "#fcbdbaff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, .8, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r3_indirect

plots <- ggarrange(r1_indirect, r2_indirect, r3_indirect,
                   labels = c("d.", "e.", "f."), 
                   ncol = 3, 
                   nrow = 1,
                   font.label = list(size = 30, family = "Helvetica", color = "black"),
                   hjust = -0.5, 
                   common.legend = T)

annotate_figure(plots, left = text_grob(" ", face = "bold", 
                                        color = "black", family = "Helvetica", rot = 90, size = 40),
                top = text_grob(" ", color = "black", family = "Helvetica",
                                face = "bold", size = 41),
                bottom = text_grob(" ", color = "black", family = "Helvetica",
                                   face = "bold", size = 35))

ggsave(paste0("plot_feeding_indirect_new.pdf"), dpi = 300,
       width = 65,  height = 22,  units = "cm")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# GLM Table
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library("broom")
library("tidyr")
library("dplyr")

# Tidy the coefficient summaries
mod1_direct_table <- tidy(mod1_direct)
mod1_indirect_table <- tidy(mod1_indirect)
class(mod1_direct_table)

# Assuming 'dplyr' package is loaded

# Combine the two data frames
combined_table <- bind_rows(
                  mod1_direct_table %>% mutate(Model = "mod1_direct"),
                  mod1_indirect_table %>% mutate(Model = "mod1_indirect"))

# Print or use the combined data frame
print(combined_table)

# Save the combined summary table to a CSV file
write.csv(combined_table, "combined_table_25m_14days.csv", row.names = FALSE)

# Print or save the combined table
print(combined_table)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# CO-ROOSTING SITUATION
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Step 1: Start by cleaning the desktop
rm(list=ls())

# Step 2: Set the working directory (modify as needed)
#setwd()

# Step 3: List files in the directory
dir()

# Step 4: Install and load necessary packages

# Load required libraries
library(ggplot2) # is a powerful and flexible plotting system in R
library(tidyverse) # is a collection of R packages designed to work together seamlessly for data analysis
library(dplyr)     # is a package for data manipulation and transformation
library(ggpubr)    # is an extension of ggplot2 that facilitates the creation and customization of publication-ready plots
library(performance) # is used for assessing and comparing the performance of statistical models.
library(DHARMa) #is a package for checking the distributional assumptions of regression models
library(lme4) #is a package for fitting linear mixed-effects models (LMMs) and generalized linear mixed-effects models (GLMMs). It extends the functionality of the base R `lm()` function to handle hierarchical and nested data structures.
library(broom.mixed) # is an extension of the broom package, which provides a consistent approach to tidying model outputs in R. broom.mixed specifically extends this functionality to mixed-effects models (LMMs and GLMMs) fitted with the lme4 package.

#*##############################
# DIRECT AND INDIRECT NETWORK
###############################

# ---- read file ----
dir()
results_DI <- read.table("output_infection_coroosting_direct.txt", sep = " ")

# Group by 'id' and select the first sampling date
results_direct <- results_DI %>%
  group_by(id) %>%
  arrange(sampling_date) %>%
  slice(1) %>%
  ungroup()

# Print the new data frame
dim(results_direct)
print(results_direct)

# ---- read file ----
results_IND <- read.table("output_infection_coroosting_indirect.txt", sep = " ")

# Group by 'id' and select the first sampling date
results_indirect <- results_IND %>%
  group_by(id) %>%
  arrange(sampling_date) %>%
  slice(1) %>%
  ungroup()

# Print the new data frame
dim(results_indirect)
print(results_indirect)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## MODELS
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Results Direct
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Dataset
dim(results_direct)

# Fit additional models for comparison
mod1_direct <- glmer(infection ~ scale(degree) + scale(strength) +
                       scale(betweenness) + as.factor(age) + (1|sampling_date), data = results_direct, family = binomial())

# Checking model assumptions
check_model(mod1_direct)

# Summary
summary(mod1_direct)

# Tidy the model results
tidy_results <- tidy(mod1_direct, conf.int = TRUE)
tidy_results$term <- c("Intercept",
                       "Degree" ,               
                       "Strength" ,             
                       "Betweenness" ,          
                       "Age",
                       "Sd Intercept")

# Create a forest plot using ggplot2
plot_1 <- ggplot(tidy_results, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(color = factor(term)), size = 3) +
  geom_errorbarh(height = 0.1) +
  labs(title = "GLMM Model Results",
       x = "Estimate",
       y = "Predictor") +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,1.2,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

print(plot_1)

ggsave(paste0("plot_roosting_direct_Estimate.pdf"), dpi = 600,
       width = 35,  height = 22,  units = "cm")

#***********************************************************************
# Coroosting direct******************************************************
#***********************************************************************

results_direct$betweenness <- as.numeric(results_direct$betweenness)
results_direct$degree <- as.numeric(results_direct$degree)
results_direct$strength <- as.numeric(results_direct$strength)

colnames(results_direct)

# Plot Predicted data and original data points for betweenness, degree and strength
r1_direct <- ggplot(results_direct, aes(x = betweenness, y = as.numeric(infection))) +
  labs(x = "Betweenness", y = "Infection status") +
  geom_point(fill = "#016cd4ff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_direct$infection), 
                                max(results_direct$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#016cd4ff", se = T, fill = "#016cd4ff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,1.2,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r1_direct

r2_direct <- ggplot(results_direct, aes(x = degree, y = as.numeric(infection))) +
  labs(x = "Degree", y = " ") +
  geom_point(fill = "#016cd4ff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_direct$infection), 
                                max(results_direct$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#016cd4ff", se = T, fill = "#016cd4ff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,.8,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r2_direct

r3_direct <- ggplot(results_direct, aes(x = strength, y = as.numeric(infection))) +
  labs(x = "Strength", y = " ") +
  geom_point(fill = "#016cd4ff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_direct$infection), 
                                max(results_direct$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#016cd4ff", se = T, fill = "#016cd4ff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35), 
        axis.title.y = element_text(color = "black", size = 35), 
        axis.text = element_text(color = "black", size = 28), 
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5,.8,.1,.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r3_direct

plots <- ggarrange(r1_direct, r2_direct, r3_direct,
                   labels = c("a.", "b.", "c."), 
                   ncol = 3, 
                   nrow = 1,
                   font.label = list(size = 30, family = "Helvetica", color = "black"),
                   hjust = -0.5, 
                   common.legend = T)

annotate_figure(plots, left = text_grob(" ", face = "bold", 
                                        color = "black", family = "Helvetica", rot = 90, size = 40),
                top = text_grob(" ", color = "black", family = "Helvetica",
                                face = "bold", size = 41),
                bottom = text_grob(" ", color = "black", family = "Helvetica",
                                   face = "bold", size = 35))

ggsave(paste0("plot_roosting_direct_new.pdf"), dpi = 600,
       width = 65,  height = 22,  units = "cm")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Results Indirect
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Dataset
dim(results_indirect)

# Fit additional models for comparison
mod1_indirect <- glmer(infection ~ scale(degree) + scale(strength) +
                         scale(betweenness) + as.factor(age) + (1|sampling_date), data = results_indirect, family = binomial())

# Checking model assumptions
check_model(mod1_indirect)

# Summary
summary(mod1_indirect)

# Tidy the model results
tidy_results <- tidy(mod1_indirect, conf.int = TRUE)
tidy_results$term <- c("Intercept",
                       "Degree",
                       "Strength",
                       "Betweenness",
                       "Age",
                       "Sd Intercept")

# Create a forest plot using ggplot2
plot_1 <- ggplot(tidy_results, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(color = factor(term)), size = 3) +
  geom_errorbarh(height = 0.1) +
  labs(title = "GLMM Model Results",
       x = "Estimate",
       y = "Predictor") +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, 1.2, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

print(plot_1)

ggsave(paste0("plot_roosting_indirect_Estimate.pdf"), dpi = 600,
       width = 35,  height = 22,  units = "cm")

#***********************************************************************
#Coroosting Indirect
#***********************************************************************
results_indirect$betweenness <- as.numeric(results_indirect$betweenness)
results_indirect$degree <- as.numeric(results_indirect$degree)
results_indirect$strength <- as.numeric(results_indirect$strength)
colnames(results_indirect)

# Plot Predicted data and original data points for degree
r1_indirect <- ggplot(results_indirect, aes(x = betweenness, y = as.numeric(infection))) +
  labs(x = "Betweenness", y = "Infection status") +
  geom_point(fill = "#b2dbffff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_indirect$infection), 
                                max(results_indirect$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#b2dbffff", se = T, fill = "#b2dbffff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, 1.2, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r1_indirect

r2_indirect <- ggplot(results_indirect, aes(x = degree, y = as.numeric(infection))) +
  labs(x = "Degree", y = " ") +
  geom_point(fill = "#b2dbffff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_indirect$infection), 
                                max(results_indirect$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#b2dbffff", se = T, fill = "#b2dbffff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, .8, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r2_indirect

r3_indirect <- ggplot(results_indirect, aes(x = strength, y = as.numeric(infection))) +
  labs(x = "Strength", y = " ") +
  geom_point(fill = "#b2dbffff", size = 10, alpha = 5/9, shape = 21, color = "black", stroke = 2) +
  scale_y_continuous(limits = c(min(results_indirect$infection), 
                                max(results_indirect$infection)),
                     breaks = c(0, 1),  labels = c("No", "Yes")) +
  stat_smooth(method = "glm", color = "#b2dbffff", se = T, fill = "#b2dbffff", alpha = 3/9,
              method.args = list(family = binomial), size = 1.5) +
  theme(axis.title.x = element_text(color = "black", size = 35),
        axis.title.y = element_text(color = "black", size = 35),
        axis.text = element_text(color = "black", size = 28),
        legend.position = "none",
        plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(.5, .8, .1, .5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(fill = 'white', colour = "black"))

r3_indirect

plots <- ggarrange(r1_indirect, r2_indirect, r3_indirect,
                   labels = c("d.", "e.", "f."), 
                   ncol = 3, 
                   nrow = 1,
                   font.label = list(size = 30, family = "Helvetica", color = "black"),
                   hjust = -0.5, 
                   common.legend = T)

annotate_figure(plots, left = text_grob(" ", face = "bold", 
                                        color = "black", family = "Helvetica", rot = 90, size = 40),
                top = text_grob(" ", color = "black", family = "Helvetica",
                                face = "bold", size = 41),
                bottom = text_grob(" ", color = "black", family = "Helvetica",
                                   face = "bold", size = 35))

ggsave(paste0("plot_roosting_indirect_new.pdf"), dpi = 300,
       width = 65,  height = 22,  units = "cm")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# GLM Table
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library("broom")
library("tidyr")
library("dplyr")

# Tidy the coefficient summaries
mod1_direct_table <- tidy(mod1_direct)
mod1_indirect_table <- tidy(mod1_indirect)
class(mod1_direct_table)

# Assuming 'dplyr' package is loaded

# Combine the two data frames
combined_table <- bind_rows(
  mod1_direct_table %>% mutate(Model = "mod1_direct"),
  mod1_indirect_table %>% mutate(Model = "mod1_indirect"))

# Print or use the combined data frame
print(combined_table)

# Save the combined summary table to a CSV file
write.csv(combined_table, "combined_table_25m_14days.csv", row.names = FALSE)

# Print or save the combined table
print(combined_table)



