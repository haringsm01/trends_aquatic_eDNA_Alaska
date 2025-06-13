################################################################################
##### AUTHOR:    Maggie Harings
##### DATE:      8/11/24
##### INPUT:     data_tidy.csv
##### OUTPUT:    circle plots of survey results
################################################################################

# start fresh!
# dev.off(dev.list()["RStudioGD"])    # clear plots
cat("\f")                           # clear console
rm(list=ls())                       # clear environment

# paste working directory file pathway between r"(.......)")
(wd_slash <- r"(G:/.shortcut-targets-by-id/1uIdPWWUw_u5QwRXArWtGUTKwv_bsBOPd/Environmental DNA review paper/data_analyses)")     

# convert \ to / in file pathway; save as wd
(wd <- gsub("\\\\", "/", wd_slash))    
wd

# set working directory
setwd(wd)

# load packages
library(tidyverse)

#### LOAD SURVEY DATA ##########################################################
df <- read_csv("./raw/data_tidy.csv",##
                trim_ws = TRUE) 

#### CONVERT DATA TYPES ########################################################

# create vector of column names to convert to factors
factor_cols <- c("entry_order", "sampling_start_year", "proj_region", "proj_status",
                 "proj_collabs", "funding", "replicates", "proj_type", "waters_sampled",
                 "target_spp", "filter_type", "filter_pore_size", "blanks", "DNA_analyses",
                 "findings_comms", "barriers", "sampling_method", "stats", "blanks",
                 "contracted_work", "manual_input")

# loop through each column name and convert it to factor in df
for (col in factor_cols) {
  df[[col]] <- as.factor(df[[col]])
}

# store total number of studies as object
nStudies <- length(unique(df$entry_order))

#### Circle plot: project type vs year #########################################

nStudies <- length(unique(df$entry_order))

# freq table 1: freq of each proj_type by year
table_yr_projTypeFreq <- df %>%
  group_by(sampling_start_year, proj_type) %>%
  summarise(freq_proj_type = n_distinct(entry_order)) %>%
  mutate(freqProjType_div_by_nStudies = freq_proj_type / nStudies)

# freq table 2: sum of total # of entries for proj_type ACROSS years
table_total_projTypes <- table_yr_projTypeFreq %>% 
  group_by(proj_type) %>% 
  summarise(total_proj_type = sum(freq_proj_type))

# table for plotting: add column to freq table 1 (table_yr_projTypeFreq): 
# freq of each proj_type by year/total of proj_type(s) entries for each year 
(table_projType_plot <- table_yr_projTypeFreq %>%
  left_join(table_total_projTypes, by = "proj_type") %>%
  mutate(prcnt_proj_type = freq_proj_type / total_proj_type))

order <- rev(c("Presence/non-detection", "Species richness",
           "Species quantification", "eDNA ecology", 
           "Invasive species", "Rare species assessment",
           "Methods comparison (field)", 
           "Methods comparison (laboratory)", 
           "Assay development/validation", 
           "Sequencing"))

# set plot theme
theme_set(theme_classic())

#  PLOT
(plot_projType_yr <- ggplot(table_projType_plot) +
    # plot circles
    geom_point(aes(x= sampling_start_year, 
                   y= proj_type,
                   # make circle size reflective of percent of region
                   # mentions for that year compared to mentions across all yrs
                   size= prcnt_proj_type)) +
    # adjust the size of the circles
    scale_size_continuous(range = c(1, 8),
                          # set legend breaks from 0-max value of data separated by .1 increments
                          breaks = seq(0, max(table_total_projTypes$total_proj_type), by = .2)) +
    labs(
      # rename x-axis
      x = "Project initiation year",
      # rename y-axis
      y = "Project type",
      # rename legend title
      size = "Annual\nproportion\nof project\ntype\nacross time") +
    # rotate x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # set axes text font size
          axis.title = element_text(size=12),
          # set legend title text font size
          legend.title = element_text(size=10),
          # add spacing between plot and legend (i.e. move legend to the right)
          # legend.box.spacing = margin(0,0,0,40),
          # adjust margin around plot (top, right, bottom, left)
          plot.margin = margin(20, 0, 0, 0, "pt")) +
    scale_y_discrete(limits = order) +
    # add column of total number of instances each region was reported
        # across all years
    geom_text(data = table_total_projTypes,
              aes(x = 13.5,
                  y = proj_type,
                  label = total_proj_type),
              size = 4,
              # fontface = "bold",
              hjust = 1) +
    annotate("text",
             x = 13.3,  # Position on the x-axis after the last year
             y = 0,  # Adjust y-position as needed
             size = 3,
             label = "Total\nmentions",
             angle = 0,  # Rotate the text if needed
             hjust = 0.5,  # Center justify horizontally
             vjust = 1.2,
             color = "gray30") +
    # allow 'Total mentions' label to exist outside of bounds of plot itself
    coord_cartesian(clip = "off") +
    # add vertical dotted line
    geom_segment(x = 12.7,
                 y = -Inf,
                 xend = 12.7,
                 yend = Inf,
                 color = "black",
                 linetype = "dashed",
                 size = 0.5) +
    scale_y_discrete(limits = rev(c(order))))
    
# export figure
ggsave(plot = plot_projType_yr, 
       filename = "./output/circle_plots/Proportional_Project_Type_by_Start_Year.png",
       width = 9,
       height = 4,
       dpi = 300)

dev.off()

#### Circle plot: filter pore size vs year #####################################

# freq table 1: freq of each filter pore size by year
table_yr_filterPoreSize <- df %>%
  group_by(sampling_start_year, filter_pore_size) %>%
  summarise(freq_filterPoreSize = n_distinct(entry_order)) %>%
  mutate(freqFilterPoreSize_div_by_nStudies = freq_filterPoreSize / nStudies)

# freq table 2: sum of total # of entries for filter pore size ACROSS years
table_total_filterPoreSize <- table_yr_filterPoreSize %>% 
  group_by(filter_pore_size) %>% 
  summarise(total_filterPoreSize = sum(freq_filterPoreSize))

# table for plotting: add column to freq table 1 (table_yr_filterPoreSize): 
# freq of each filter_pore_size by year/total of filter_pore_size(s) entries for each year 
(table_filterPoreSize_plot <- table_yr_filterPoreSize %>%
    left_join(table_total_filterPoreSize, by = "filter_pore_size") %>%
    mutate(prcnt_filterPoreSize = freq_filterPoreSize / total_filterPoreSize))

order <- c( "Unknown", "None", "0.1", "0.22","0.4", "0.45", "0.45-1.5", "1","1.2","1.5","3", 
            "5", "7", "10")

# set plot theme
theme_set(theme_classic())

#  PLOT
(plot_filterPoreSize <- ggplot(table_filterPoreSize_plot) +
    # plot circles
    geom_point(aes(x= sampling_start_year, 
                   y= filter_pore_size,
                   # make circle size reflective of percent of filter pore size
                   # mentions for that year compared to mentions across all yrs
                   size = prcnt_filterPoreSize)) +
    # adjust the size of the circles
    scale_size_continuous(range = c(1, 8),
                          # set legend breaks from 0-max value of data separated by .1 increments
                          breaks = seq(0, max(table_total_filterPoreSize$total_filterPoreSize), by = .2)) +
    labs(
      # rename x-axis
      x = "Project initiation year",
      # rename y-axis
      y = "Filter pore size (um)",
      # rename legend title
      size = "Annual\nproportion\nof filter\npore size\nacross time") +
    # rotate x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # set axes text font size
          axis.title = element_text(size=12),
          # set legend title text font size
          legend.title = element_text(size=10),
          # add spacing between plot and legend (i.e. move legend to the right)
          legend.box.spacing = margin(0,0,0,40),
          # adjust margin around plot (top, right, bottom, left)
          plot.margin = margin(20, 0, 0, 0, "pt")) +
    # add column of total number of instances each region was reported
    # across all years
    geom_text(data = table_filterPoreSize_plot,
              aes(x = 13.5,
                  y = filter_pore_size,
                  label = total_filterPoreSize),
              size = 4,
              # fontface = "bold",
              hjust = 1) +
    annotate("text",
             x = 13.3,  # Position on the x-axis after the last year
             y = 0,  # Adjust y-position as needed
             size = 3,
             label = "Total\nmentions",
             angle = 0,  # Rotate the text if needed
             hjust = 0.5,  # Center justify horizontally
             vjust = 1.2,
             color = "gray30") +
    # allow 'Total mentions' label to exist outside of bounds of plot itself
    coord_cartesian(clip = "off") +
    # add vertical dotted line
    geom_segment(x = 12.7, 
                 y = -Inf, 
                 xend = 12.7, 
                 yend = Inf,
                 color = "black",
                 linetype = "dashed",
                 size = 0.5)+
    scale_y_discrete(limits = c(order))
)

# export figure
ggsave(plot = plot_filterPoreSize, 
       filename = "./output/circle_plots/Proportional_Filter_Pore_Size_by_Start_Year.png",
       width = 9,
       height = 4,
       dpi = 300)

dev.off()

#### Circle plot: filter type vs year ##########################################

# freq table 1: freq of each filter type size by year
table_yr_filterType <- df %>%
  group_by(sampling_start_year, filter_type ) %>%
  summarise(freq_filterType = n_distinct(entry_order)) %>%
  mutate(freqFilterType_div_by_nStudies = freq_filterType / nStudies)

# freq table 2: sum of total # of entries for filter types ACROSS years
table_total_filterType <- table_yr_filterType %>% 
  group_by(filter_type) %>% 
  summarise(total_filterType = sum(freq_filterType))

# table for plotting: add column to freq table 1 (table_yr_filterType): 
# freq of each filter_type by year/total of filter_type(s) entries for each year 
(table_filterType_plot <- table_yr_filterType %>%
    left_join(table_total_filterType, by = "filter_type") %>%
    mutate(prcnt_filterType = freq_filterType / total_filterType))

# set plot theme
theme_set(theme_classic())

#  PLOT
(plot_filterType <- ggplot(table_filterType_plot) +
    # plot circles
    geom_point(aes(x= sampling_start_year, 
                   y= filter_type,
                   # make circle size reflective of percent of filter pore size
                   # mentions for that year compared to mentions across all yrs
                   size = prcnt_filterType)) +
    # adjust the size of the circles
    scale_size_continuous(range = c(1, 8),
                          # set legend breaks from 0-max value of data separated by .1 increments
                          breaks = seq(0, max(table_total_filterType$total_filterType), by = .2)) +
    # reverse order of y-axis labels
    scale_y_discrete(limits = rev(unique(table_filterType_plot$filter_type))) +
    labs(
      # rename x-axis
      x = "Project initiation year",
      # rename y-axis
      y = "Filter type",
      # rename legend title
      size = "Annual\nproportion\nof filter\ntypes\nacross time") +
    # rotate x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # set axes text font size
          axis.title = element_text(size=12),
          # set legend title text font size
          legend.title = element_text(size=10),
          # add spacing between plot and legend (i.e. move legend to the right)
          legend.box.spacing = margin(0,0,0,40),
          # adjust margin around plot (top, right, bottom, left)
          plot.margin = margin(20, 0, 0, 0, "pt")) +
    # add column of total number of instances each region was reported
    # across all years
    geom_text(data = table_filterType_plot,
              aes(x = 13.5,
                  y = filter_type,
                  label = total_filterType),
              size = 4,
              # fontface = "bold",
              hjust = 1) +
    annotate("text",
             x = 13.3,  # Position on the x-axis after the last year
             y = 0,  # Adjust y-position as needed
             size = 3,
             label = "Total\nmentions",
             angle = 0,  # Rotate the text if needed
             hjust = 0.5,  # Center justify horizontally
             vjust = 1.2,
             color = "gray30") +
    # allow 'Total mentions' label to exist outside of bounds of plot itself
    coord_cartesian(clip = "off") +
    # add vertical dotted line
    geom_segment(x = 12.7, 
                 y = -Inf, 
                 xend = 12.7, 
                 yend = Inf,
                 color = "black",
                 linetype = "dashed",
                 size = 0.5))

# export figure
ggsave(plot = plot_filterPoreSize, 
       filename = "./output/circle_plots/Proportional_Filter_Pore_Size_by_Start_Year.png",
       width = 9,
       height = 4,
       dpi = 300)

dev.off()

#### Circle plot: blanks vs year #########################################

# freq table 1: freq of each region by year
table_yr_blanksFreq <- df %>%
  group_by(sampling_start_year, blanks) %>%
  summarise(freq_blanks = n_distinct(entry_order)) %>%
  mutate(freqblanks_div_by_nStudies = freq_blanks / nStudies)

# freq table 2: sum of total # of entries for proj_region ACROSS years
table_total_blanks <- table_yr_blanksFreq %>% 
  group_by(blanks) %>% 
  summarise(total_blanks = sum(freq_blanks))

# table for plotting: add column to freq table 1 (table_yr_regionFreq): 
# freq of each proj_region by year/total of proj_region(s) entries for each year 
(table_blanks_plot <- table_yr_blanksFreq %>%
    left_join(table_total_blanks, by = "blanks") %>%
    mutate(prcnt_blanks = freq_blanks / total_blanks))

order <- c("Sample collection", "Filter subsetting","DNA extractions", 
           "DNA analyses", "None", "Unknown")

# set plot theme
theme_set(theme_classic())

#  PLOT
(plot_blanks<- ggplot(table_blanks_plot) +
    # plot circles
    geom_point(aes(x= sampling_start_year, 
                   y= blanks,
                   # make circle size reflective of percent of blanks
                   # mentions for that year compared to mentions across all yrs
                   size= prcnt_blanks)) +
    # adjust the size of the circles
    scale_size_continuous(range = c(1, 8),
                          # set legend breaks from 0-max value of data separated by .1 increments
                          breaks = seq(0, max(table_total_blanks$total_blanks), by = .2)) +
    labs(
      # rename x-axis
      x = "Project initiation year",
      # rename y-axis
      y = "Negative control use",
      # rename legend title
      size = "Annual\nproportion\nof negative\nrcontrol use\nacross time") +
    # rotate x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # set axes text font size
          axis.title = element_text(size=12),
          # set legend title text font size
          legend.title = element_text(size=10),
          # add spacing between plot and legend (i.e. move legend to the right)
          legend.box.spacing = margin(0,0,0,40),
          # adjust margin around plot (top, right, bottom, left)
          plot.margin = margin(20, 0, 0, 0, "pt")) +
    # add column of total number of instances each region was reported
    # across all years
    geom_text(data = table_blanks_plot,
              aes(x = 13.5,
                  y = blanks,
                  label = total_blanks),
              size = 4,
              # fontface = "bold",
              hjust = 1) +
    annotate("text",
             x = 13.3,  # Position on the x-axis after the last year
             y = 0,  # Adjust y-position as needed
             size = 3,
             label = "Total\nmentions",
             angle = 0,  # Rotate the text if needed
             hjust = 0.5,  # Center justify horizontally
             vjust = 1.2,
             color = "gray30") +
    # allow 'Total mentions' label to exist outside of bounds of plot itself
    coord_cartesian(clip = "off") +
    # add vertical dotted line
    geom_segment(x = 12.7, 
                 y = -Inf, 
                 xend = 12.7, 
                 yend = Inf,
                 color = "black",
                 linetype = "dashed",
                 size = 0.5)+
    scale_y_discrete(limits = rev(c(order)))
)


# export figure
ggsave(plot = plot_blanks, 
       filename = "./output/circle_plots/Proportional_Use_of_NegControls_by_Start_Year.png",
       width = 9,
       height = 4,
       dpi = 300)

dev.off()


#### Circle plot: filter pore size vs waterbody type ###########################

# freq table 1: freq of each pSize by waterbody type
table_water_pSize_plot <- df %>%
  group_by(waters_sampled, filter_pore_size) %>%
  summarise(freq_pSize = n_distinct(entry_order)) 

# freq table 2: sum of total # of entries for filter pore size ACROSS all waterbodies
table_total_pSize<- table_water_pSize_plot %>% 
  group_by(filter_pore_size) %>% 
  summarise(total_pSize = sum(freq_pSize))

order_x <- c("Wetland" ,"Lake",  "River/stream", "Estuary", 
             "Tidewater glacier fjords","Ocean (nearshore)",
             "Ocean (mesopelagic)", "Ocean (pelagic)")

order_y <- c( "Unknown", "None", "0.1", "0.22","0.4", "0.45", "0.45-1.5", "1","1.2","1.5","3",
              "5", "7", "10")

# set plot theme
theme_set(theme_classic())

#  PLOT
(plot_water_pSize <- ggplot(table_water_pSize_plot) +
    # plot circles
    geom_point(aes(x= waters_sampled , 
                   y= filter_pore_size ,
                   # make circle size reflective of percent of filter_pore_size 
                   # mentions for that waterbody type
                   size = freq_pSize)) +
    # adjust the size of the circles
    scale_size_continuous(range = c(1, 8),
                          # set legend breaks from 0-max value of data separated by .1 increments
                          breaks = seq(0, max(table_water_pSize_plot$freq_pSize), by = 1)) +
    labs(
      # rename x-axis
      x = "Waterbody type",
      # rename y-axis
      y = "Filter pore size (um)",
      # rename legend title
      size = "Frequency of\nfilter pore size\nused by waterbody\ntype sampled") +
    # rotate x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # set axes text font size
          axis.title = element_text(size=12),
          # set legend title text font size
          legend.title = element_text(size=10),
          # add spacing between plot and legend (i.e. move legend to the right)
          legend.box.spacing = margin(0,0,0,40),
          # adjust margin around plot (top, right, bottom, left)
          plot.margin = margin(20, 0, 0, 0, "pt")) +
    # add column of total number of instances each pore size was reported
    # across each waterbody type
    geom_text(data = table_total_pSize,
              aes(x = 9,
                  y = filter_pore_size,
                  label = total_pSize),
              size = 4,
              # fontface = "bold",
              hjust = 1) +
    annotate("text",
             x = 9,  # Position on the x-axis after the last year
             y = 0,  # Adjust y-position as needed
             size = 3,
             label = "Total\nmentions",
             angle = 0,  # Rotate the text if needed
             hjust = 0.5,  # Center justify horizontally
             vjust = 1.2,
             color = "gray30") +
    # allow 'Total mentions' label to exist outside of bounds of plot itself
    coord_cartesian(clip = "off") +
    # add vertical dotted line
    geom_segment(x = 8.5,
                 y = -Inf,
                 xend = 8.5,
                 yend = Inf,
                 color = "black",
                 linetype = "dashed",
                 size = 0.5)+
    scale_x_discrete(limits = rev(c(order_x))) +
    scale_y_discrete(limits = c(order_y))
)


# export figure
ggsave(plot = plot_water_pSize, 
       filename = "./output/circle_plots/Freq_filtPoreSize_by_waterbody.png",
       width = 10,
       height = 4.2,
       dpi = 300)

dev.off()

#### Circle plot: filter pore size vs species ##################################

# freq table 1: freq of each pSize by species
table_spp_pSize_plot_spp <- df %>%
  group_by(target_spp, filter_pore_size) %>%
  summarise(freq_pSize = n_distinct(entry_order)) 

# freq table 2: sum of total # of entries for filter pore size ACROSS each spp. type
table_total_pSize_spp <- table_spp_pSize_plot_spp %>% 
  group_by(filter_pore_size) %>% 
  summarise(total_pSize = sum(freq_pSize))

order_x <- c("Fishes", "Macroinvertebrates",       
             "Cephalopods", "Pathogens", "Crustaceans", 
             "Mollusks","Mammals (marine)", "Mammals (land-based)",     
             "Microalgae" ,"Plants (macrophytes)", "Chordates","All species")

order_y <- c( "Unknown", "None", "0.1", "0.22","0.4", "0.45", "0.45-1.5", "1","1.2","1.5","3",
              "5", "7", "10")

# set plot theme
theme_set(theme_classic())

#  PLOT
(plot_spp_pSize <- ggplot(table_spp_pSize_plot_spp) +
    # plot circles
    geom_point(aes(x= target_spp , 
                   y= filter_pore_size ,
                   # make circle size reflective of percent of filter_pore_size 
                   # mentions for that waterbody type
                   size = freq_pSize)) +
    # adjust the size of the circles
    scale_size_continuous(range = c(1, 8),
                          # set legend breaks from 0-max value of data separated by .1 increments
                          breaks = seq(0, max(table_spp_pSize_plot_spp$freq_pSize), by = 1)) +
    labs(
      # rename x-axis
      x = "Target species",
      # rename y-axis
      y = "Filter pore size (um)",
      # rename legend title
      size = "Frequency of\nfilter pore size\nused by target\nspecies") +
    # rotate x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # set axes text font size
          axis.title = element_text(size=12),
          # set legend title text font size
          legend.title = element_text(size=10),
          # add spacing between plot and legend (i.e. move legend to the right)
          legend.box.spacing = margin(0,0,0,40),
          # adjust margin around plot (top, right, bottom, left)
          plot.margin = margin(20, 0, 0, 0, "pt")) +
    # add column of total number of instances each region was reported
    # across all years
    geom_text(data = table_total_pSize_spp,
              aes(x = 13,
                  y = filter_pore_size,
                  label = total_pSize),
              size = 4,
              # fontface = "bold",
              hjust = 1) +
    annotate("text",
             x = 13,  # Position on the x-axis after the last year
             y = 0,  # Adjust y-position as needed
             size = 3,
             label = "Total\nmentions",
             angle = 0,  # Rotate the text if needed
             hjust = 0.5,  # Center justify horizontally
             vjust = 1.2,
             color = "gray30") +
    # allow 'Total mentions' label to exist outside of bounds of plot itself
    coord_cartesian(clip = "off") +
    # add vertical dotted line
    geom_segment(x = 12.5,
                 y = -Inf,
                 xend = 12.5,
                 yend = Inf,
                 color = "black",
                 linetype = "dashed",
                 size = 0.5)+
    scale_x_discrete(limits = c(order_x)) +
    scale_y_discrete(limits = c(order_y))
)


# export figure
ggsave(plot = plot_spp_pSize, 
       filename = "./output/circle_plots/Freq_filtPoreSize_by_species.png",
       width = 9.1,
       height = 4.5,
       dpi = 300)

dev.off()


