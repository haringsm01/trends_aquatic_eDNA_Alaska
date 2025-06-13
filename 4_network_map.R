################################################################################
##### AUTHOR:    Maggie Harings
##### DATE:      11/5/24; updated 6/3/25
##### INPUT:     data_raw.csv
##### OUTPUT:    network map tally tables
################################################################################
# start fresh!
# dev.off(dev.list()["RStudioGD"])    # clear plots
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
library(ggraph)    # network maps
# install.packages("igraph", type="binary")
library(igraph)    # network maps

#### LOAD SURVEY DATA ##########################################################

data_tem <- read_csv("./raw/data_raw.csv", 
                     trim_ws = TRUE)

# correct sampling_start for entry #45
data_tem$sampling_start[data_tem$entry_order == "45"] <- "10-04"

# remove row for proposed project
data_temp <- data_tem %>%
  filter(sampling_start != "Proposed")

#### ADD COLUMN: SAMPLING START YEAR ###########################################

# add sampling_start_year column
data_temp$sampling_start_year <- as.factor(paste("20",
                                                 str_sub(data_temp$sampling_start,1,2),
                                                 sep=""))

#### CONVERT DATA TYPES ########################################################

# check data types
str(data_temp)

# convert dat_time to POSIXct
data_temp$date_time <- as.POSIXct(data_temp$date_time, 
                                  format = "%m/%d/%Y %H:%M")

# convert proj_status to factor
data_temp$proj_status <- as.factor(data_temp$proj_status)


# convert proj_region, proj_status, proj_collabs, proj_type, waters_sampled, 
# target_spp, filter_type, filter_pore_size, blanks, DNA_analyses, 
# findings_comms, barriers, sampling_method, blanks, blanks, contracted_work, 
# manual_input to factors

data_temp$proj_region <- as.factor(data_temp$proj_region)
data_temp$proj_status <- as.factor(data_temp$proj_status)
data_temp$proj_collabs <- as.factor(data_temp$proj_collabs)
data_temp$proj_type <- as.factor(data_temp$proj_type)
data_temp$waters_sampled <- as.factor(data_temp$waters_sampled)
data_temp$target_spp <- as.factor(data_temp$target_spp)
data_temp$filter_type <- as.factor(data_temp$filter_type)
data_temp$filter_pore_size <- as.factor(data_temp$filter_pore_size)
data_temp$blanks <- as.factor(data_temp$blanks)
data_temp$DNA_analyses <- as.factor(data_temp$DNA_analyses)
data_temp$findings_comms <- as.factor(data_temp$findings_comms)
data_temp$barriers <- as.factor(data_temp$barriers)
data_temp$sampling_method <- as.factor(data_temp$sampling_method)
data_temp$stats <- as.factor(data_temp$stats)
data_temp$blanks  <- as.factor(data_temp$blanks)
data_temp$contracted_work  <- as.factor(data_temp$contracted_work)
data_temp$manual_input  <- as.factor(data_temp$manual_input)

# check structure to ensure columns were modified to correct data types
str(data_temp)

#### REORDER COLUMNS ###########################################################

colnames(data_temp)

data_colOrder <- data_temp %>% 
  relocate(c(entry_order,sampling_start,sampling_start_year, sampling_end,proj_region,waters_sampled,    
             proj_type, target_spp, sampling_method, filter_type, 
             filter_pore_size,DNA_analyses, stats,blanks, replicates,     
             findings_comms,barriers, contracted_work,proj_collabs, 
             funding, proj_status, sampling_locations, publication_link,   
             emails, acknowledgements,  manual_input, date_time))


#### REMOVE 2024 PROJECTS ######################################################

# remove sampling_start_year  == "2024
data <- data_colOrder[data_colOrder$sampling_start_year != "2024",]

#### CREATE NETWORK MAP DATA: PROJECT COLLABORATIONS ###########################

# Function to generate unique pairs or flag "No outside collaborations"
generate_pairs <- function(observations) {
  # Split by comma and trim whitespace
  obs <- str_split(observations, ",\\s*")[[1]]
  obs <- unique(str_trim(obs))  # Remove duplicates
  
  # If only one collaboration exists, return it in COLLAB1 and "No outside collaborations" in COLLAB2
  if (length(obs) == 1) {
    return(data.frame(COLLAB1 = obs, COLLAB2 = "No outside collaborations", stringsAsFactors = FALSE))
  }
  
  # Generate combinations (pairs) in the form "item1-item2"
  pair_combinations <- combn(sort(obs), 2, FUN = function(x) paste(x, collapse = "-"), simplify = TRUE)
  
  # Return a data frame of pairs with COLLAB1 and COLLAB2 as separate columns
  result <- data.frame(COLLAB1 = sapply(str_split(pair_combinations, "-"), `[`, 1),
                       COLLAB2 = sapply(str_split(pair_combinations, "-"), `[`, 2), 
                       stringsAsFactors = FALSE)
  
  return(result)
}


# Apply the function and count each pair
pair_counts <- data %>%
  # Use `mutate` to generate pairs or "No outside collaborations" for each row
  mutate(pairs = map(proj_collabs, generate_pairs)) %>%
  # Unnest pairs to get each pair or "No outside collaborations" on a separate row
  unnest(pairs) %>%
  # Group by pair and count occurrences
  count(COLLAB1, COLLAB2, sort = TRUE)

# Now create the expanded table by repeating rows based on their counts
expanded_table_temp <- pair_counts %>%
  # Use `uncount` to repeat the pairs based on their counts
  uncount(n) %>%
  select(COLLAB1, COLLAB2)

expanded_table <- expanded_table_temp %>%
  left_join(pair_counts, by = c("COLLAB1", "COLLAB2"))


#### GENERATE NETWORK MAP: PROJECT COLLABORATIONS ##############################
# https://r-graph-gallery.com/257-input-formats-for-network-charts.html
# https://r-graph-gallery.com/247-network-chart-layouts.html
# color lines by year, etc: https://ggraph.data-imaginist.com/articles/Layouts.html
#https://kateto.net/network-visualization

# Create the network object with 'n' as edge weights
network <- graph_from_data_frame(d = expanded_table, directed = FALSE)

# Add 'n' as edge weights (to represent the number of collaborations)
E(network)$weight <- expanded_table$n

# wrap labels that are too long so they fit inside the nodes
V(network)$label <- sapply(V(network)$name, function(label) {
  str_wrap(label, width = 8)  # Adjust 'width' based on how many characters you want per line
})

# change color of just 'no outside collaborations' edge 
V(network)$color <- "turquoise" # Set default color for all vertices
# print(V(network)$name) 
V(network)$color[V(network)$name == "No outside collaborations"] <- "gray90"  # Ensure you're using the correct vertex name (not index)

# change color of just 'no outside collaborations' edge frame
V(network)$frame.color <- "turquoise4"
V(network)$frame.color[V(network)$name == "No outside collaborations"] <- "gray50"  # Replace with the actual name or ID you want

tiff("./output/network_maps/network_plot.tiff", width = 13, height = 13, units = "in", res = 300)
plot.igraph(network, 
            layout = layout_in_circle, 
            vertex.frame.width = 3,               # node border width
            vertex.size = 31,                     # node size
            vertex.label.color = "black",         # node text color
            vertex.label.cex = 1.2,               # node text size
            vertex.frame.cex = 100,               # node border thickness
            edge.color = "gray",                  # edge (line) color
            edge.width = E(network)$weight,       # weighted edge widths based on occurrance of relationship
            # edge.label = E(network)$weight,       # Label with weight (number of collaborations)
            # edge.label.cex = 1,                   # Increase label size
            # edge.label.dist = -3,                 # Offset from edge line (try 0.5–1)
            # edge.label.color = "black",           # edge label text color
            # #can use this to modify location of edge label but ended up just doing it manually
            # edge .label.y = c(0, 0, 0, 0, 0,      # 1-5
            #                  0, 0, 0, 0, 0,      # 6-10
            #                  0, 0, 0, 0, 0,      # 11-15
            #                  0, 0, 0, 0, 0,      # 16-20
            #                  0, 0, 0, 0, 0,      # 21-25
            #                  0, 0, 0, 0),        # 26-29
            edge.curved = 0)                      # do not curve edges
dev.off()

