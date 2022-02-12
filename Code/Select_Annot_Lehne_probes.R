# Intro and pseudocode
############################################

# 1
# Open EPIC annotation file. 
# Pull one region from a single chromosome with probe names - 
# 60 or so probes would be cool
# Generate a distance matrix

# 2
# Use those probe names to filter through the Lenhe file and find in file

############################################
# Packages
############################################
library(tidyverse)
library(readr)
# First try to load the annotation file for
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


# Pull off the location info
# ############################################
data(Locations);force(Locations)

# Pull off a chromosome and region
chroms <-Locations@listData$chr == "chr1"

Locations[chroms,] 
# 82k or so

# Look at density
Locations[chroms,] %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = pos))+
  geom_density()+
  theme_bw()

# Now add positions
pos <-Locations@listData$pos > 1.6e8 & Locations@listData$pos < 1.65e8

Locations[chroms & pos,]
# Now about 2000 probes


# Look at where it is relative to the full chromosome
Locations[chroms,] %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = pos))+
  geom_density()+
  theme_bw()+ 
  geom_rect(aes(xmin = 1.6e8, 
                xmax = 1.65e8,
                ymin = 0, 
                ymax = 1e-8))


# Check it out in just that region
Locations[chroms & pos,] %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = pos))+
  geom_density()+
  theme_bw() 

my_locations <-Locations[chroms & pos,] %>% 
  as.data.frame() %>% 
  rownames_to_column("probe_id")

write_csv(my_locations, here::here("Output/Data", "my_locations.csv"))

my_locations_probes <-my_locations$probe_id

############################################
# Read in Lehne data (NOTE THAT IT"S NOT EPIC ARRAY but 450K array)
# Select a subset of columns (100)
# Then remove the Detection P-value rows
############################################
dat=read_tsv(file = here::here("Data", "GSE55763_normalized_betas.txt"), col_select = 1:100)

my_chr_1_dat <-dat %>% 
  select(-contains("Detection")) %>% 
  filter(ID_REF %in% my_locations_probes)


# Write File 
############################################
write_csv(my_chr_1_dat, here::here("Output/Data", "my_chr_1_dat.csv"))


