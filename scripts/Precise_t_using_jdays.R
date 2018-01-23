##########################################################################
#More precise time interval for trees from 2008 - 2017 and cored trees
#based on mortality censuses in 2008 and 2017 and trees cored in 2010/11
#For use with SCBI ForestGEO grid
#Author: Ryan Helcoski - 12\/18/2018
##########################################################################

rm(list = ls())

library(dplyr)
library (kimisc)
library(lubridate)


#############################################################
##### creating wd and loading data
# Set up working directory ####
setwd("C:/Users/helcoskiR/Dropbox (Smithsonian)/Tree Cores/literature/effect of coring on mortality/Data/R_project")

# INPUT DATA LOCATION ####
Input_data_location <- "INPUT_FILES/"

# OUTPUT DATA LOCATION ####
Output_data_location <- "OUTPUT_FILES/"

# load data ####
# load 2017 mortality data + trees cored and uncored. ####
#File is csv of 2017 mortality with one extra column, Cored_2009. A 1 indicates a tree was cored a 0 indicates it was not
cmortall <- read.csv(paste0(Input_data_location, "Mortality_Survey_2017+core2009.csv"))
s2008 <- read.csv(paste0(Input_data_location, "scbi.stem1_tori.csv"))
c2010 <- read.csv(paste0(Input_data_location, "SCBI_SIGEO_all_trees_cored.csv"))
#add column for Cored_2009
c2010$Cored_2009 <- 1

#merge files
m1 <- merge(cmortall, c2010, by.x= c("tag", "Cored_2009"), by.y =c("tag", "Cored_2009") , all.x=TRUE)
m2 <- merge(m1, s2008, by.x=c("tag", "stem"), by.y=c("tag", "stem"))

#convert to jdays
m2$survey_2017 <- NA
m2$survey_2017 <- as.Date(m2$date.x, "%m/%d/%Y")
m2$survey_2017 <- as.numeric(m2$survey_2017)
m2$survey_2008 <- NA
m2$survey_2008 <- as.Date(m2$ExactDate, "%m/%d/%Y")
m2$survey_2008 <- as.numeric(m2$survey_2008)
m2$cored_2010 <- NA
m2$cored_2010 <- as.Date(m2$Date_Collected, "%m/%d/%Y")
m2$cored_2010 <- as.numeric(m2$cored_2010)

#determine number of jdays since 2008 survey for all uncored, and jdays since 2010/11 for all cored
m2$t_uncored <- ifelse (m2$Cored_2009 == 0, m2$survey_2017 - m2$survey_2008, NA)
m2$t_cored <- ifelse (m2$Cored_2009 == 1, m2$survey_2017 - m2$cored_2010, NA)

#save csv
write.csv(m2, "Raw_for_t_calculation.csv")

#determine new mean t_cored and t_uncored
cored_t <- mean(m2$t_cored, na.rm=TRUE) / 365
uncored_t <- mean(m2$t_uncored, na.rm=TRUE) / 365
t <-cbind(cored_t,uncored_t)
write.csv(t, "New_t_cored_uncored.csv")

