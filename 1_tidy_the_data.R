################################################################################
##### AUTHOR:    Maggie Harings
##### DATE:      1/4/24; last update: 6/26/24
##### INPUT:     data_raw.csv
##### OUTPUT:    data_tidy.csv
##### PURPOSE:   tidy survey data for analyses and figure generation
################################################################################

# start fresh!
dev.off(dev.list()["RStudioGD"])    # clear plots
cat("\f")                           # clear console
rm(list=ls())                       # clear environment

# paste working directory file pathway between r"(.......)")
(wd_slash <- r"(G:\.shortcut-targets-by-id\1uIdPWWUw_u5QwRXArWtGUTKwv_bsBOPd\Environmental DNA review paper\data_analyses\)")     

# convert \ to / in file pathway; save as wd
(wd <- gsub("\\\\", "/", wd_slash))    
wd

# set working directory
setwd(wd)

# load packages
library(tidyverse)
library(scales)    # converts numbers to Percents_labels with % signs (as characters)

#### LOAD SURVEY DATA + CLEAN ##################################################

data_tem <- read_csv("./raw/data_raw.csv", 
                     trim_ws = TRUE)

# correct sampling_start for entry #45
data_tem$sampling_start[data_tem$entry_order == "45"] <- "10-04"


data_temp <- data_tem %>%
  # remove row for proposed project
  filter(sampling_start != "Proposed") %>% 
         # add sampling_start_year column
  mutate(sampling_start_year = as.factor(paste("20",
                                               str_sub(sampling_start,1,2),
                                               sep="")),
         # modify data types of each column
         date_time = as.POSIXct(date_time, 
                                format = "%m/%d/%Y %H:%M"),
         proj_status = as.factor(proj_status), 
         proj_region = as.factor(proj_region),
         proj_status <- as.factor(proj_status),
         proj_collabs <- as.factor(proj_collabs),
         proj_type <- as.factor(proj_type),
         waters_sampled <- as.factor(waters_sampled),
         target_spp <- as.factor(target_spp),
         filter_type <- as.factor(filter_type),
         filter_pore_size <- as.factor(filter_pore_size),
         blanks <- as.factor(blanks),
         DNA_analyses <- as.factor(DNA_analyses),
         findings_comms <- as.factor(findings_comms),
         barriers <- as.factor(barriers),
         sampling_method <- as.factor(sampling_method),
         stats <- as.factor(stats),
         blanks  <- as.factor(blanks),
         contracted_work  <- as.factor(contracted_work),
         manual_input  <- as.factor(manual_input),
         # get rid of reference to 'including lampreys'
         target_spp = if_else(target_spp == 'Fish (including lampreys)', 'Fishes',target_spp))

# check data types
str(data_temp)

#### REORDER COLUMNS ###########################################################
data_colOrder <- data_temp %>% 
  relocate(c(entry_order,sampling_start,sampling_start_year, sampling_end,proj_region,waters_sampled,    
             proj_type, target_spp, sampling_method, filter_type, 
             filter_pore_size,DNA_analyses, stats,blanks, replicates,     
             findings_comms,barriers, contracted_work,proj_collabs, 
             funding, proj_status, sampling_locations, publication_link,   
             emails, acknowledgements,  manual_input, date_time))

#### TIDY DATA #################################################################

# transform wide --> long based on: proj_region
data_long <- data_colOrder %>% 
  separate_rows(sampling_start_year, sep = ",") %>% 
  separate_rows(proj_region, sep = ",") %>% 
  separate_rows(waters_sampled, sep = ",") %>% 
  separate_rows(proj_type, sep = ",") %>% 
  separate_rows(target_spp, sep = ",") %>% 
  separate_rows(sampling_method, sep = ",") %>% 
  separate_rows(filter_type, sep = ",") %>% 
  separate_rows(filter_pore_size, sep = ",") %>% 
  separate_rows(DNA_analyses, sep = ",") %>% 
  separate_rows(stats, sep = ",") %>% 
  separate_rows(blanks, sep = ",") %>% 
  separate_rows(replicates, sep = ",") %>% 
  separate_rows(findings_comms, sep = ", ") %>% 
  separate_rows(barriers, sep = ",") %>% 
  separate_rows(contracted_work, sep = ",") %>% 
  separate_rows(proj_collabs, sep = ",") %>% 
  separate_rows(funding, sep = ",") %>% 
  separate_rows(proj_status, sep = ",")

# trim white space (do this by column to reduce R run time)
data_long$sampling_start_year   <- trimws(data_long$sampling_start_year, which = "both")
data_long$proj_region   <- trimws(data_long$proj_region, which = "both")
data_long$waters_sampled  <- trimws(data_long$waters_sampled, which = "both")
data_long$proj_type       <- trimws(data_long$proj_type, which = "both")
data_long$target_spp      <- trimws(data_long$target_spp, which = "both")
data_long$filter_type     <- trimws(data_long$filter_type, which = "both")
data_long$filter_pore_size <- trimws(data_long$filter_pore_size, which = "both")
data_long$DNA_analyses     <- trimws(data_long$DNA_analyses, which = "both")
data_long$stats    <- trimws(data_long$stats, which = "both")
data_long$blanks    <- trimws(data_long$blanks, which = "both")
data_long$replicates    <- trimws(data_long$replicates, which = "both")
data_long$findings_comms  <- trimws(data_long$findings_comms, which = "both")
data_long$barriers  <- trimws(data_long$barriers, which = "both")
data_long$contracted_work <- trimws(data_long$contracted_work, which = "both")
data_long$proj_collabs <- trimws(data_long$proj_collabs, which = "both")
data_long$funding <- trimws(data_long$funding, which = "both")
data_long$proj_status  <- trimws(data_long$proj_status, which = "both")

#### REMOVE 2024 PROJECTS ######################################################

# total data_long rows - total data_long rows with just 2024
minus2024 <-count(data_long[data_long$sampling_start_year,]) - count(data_long[data_long$sampling_start_year == "2024",])

# remove sampling_start_year  == "2024
data_long2 <- data_long[data_long$sampling_start_year != "2024",]

# ensure 2024 was removed
count(data_long2[data_long2$sampling_start_year,]) == minus2024

#### CONVERT ATKA ISLAND TO ALEUTIAN ISLANDS  ##################################
data_long2$proj_region[data_long2$proj_region == "Atka Island"] <- "Aleutian Islands"

#### CONVERT BERING CHUKCHI TO BERING, CHUKCHI SEAS  #########@#################
data_long2$proj_region[data_long2$proj_region == "Bering Chukchi"] <- "Bering, Chukchi Seas"

#### CONVEVERT DISPO...WHIRLPAK TO WHIR-PAK ####################################
data_long2$sampling_method[data_long2$sampling_method == "Disposable whirlpak" ] <- "Whirl-Pak"

#### CORRECT NUMBER ############################################################
data_long2$filter_pore_size[data_long2$filter_pore_size == "1.0"] <- "1"

#### RENAME FINAL DATASET ######################################################

data <- data_long2

#### EXPORT TIDY DATASET #######################################################

write.csv(data, "./raw/data_tidy.csv", row.names = FALSE)

dev.off()

