# Intro and pseudocode
############################################
# 1
# Open EPIC annotation file. 
# Pull one region from a single chromosome with probe names - 
# 60 or so probes would be cool
# Generate a distance matrix

# 2
# Use those probe names to filter through the Lenhe file and find in file
# Extract those and setup as necessary for my correlation matrix. 
# 

# Packages
############################################
library(tidyverse)



############################################


############################################



# First try to load the annotation file for
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# Pull off the location info
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





newurl <-"https://cumccolumbia-my.sharepoint.com/:t:/r/personal/db3275_cumc_columbia_edu/Documents/Projects/MK/GSE55763_Lehne/GSE55763_normalized_betas.txt?csf=1&web=1&e=ULcgpW"



