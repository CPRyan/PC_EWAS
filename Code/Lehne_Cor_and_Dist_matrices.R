# Extract those and setup as necessary for my correlation matrix. 
library(tidyverse)
library(readr)


# Load data
############################################
my_locations <-read_csv(here::here("Output/Data", "my_locations.csv"))
my_chr_1_dat <-read_csv(here::here("Output/Data", "my_chr_1_dat.csv"))

joined_ordered_df <-inner_join(my_locations, my_chr_1_dat, by = c("probe_id" = "ID_REF"))
# Make distance matrix
############################################
# pairwise distances
# dist_mat <-dist(dist_vect, diag = TRUE)

# First filter my_locations to only those on 450k (in the DNAm dataset)

ten_pos <-joined_ordered_df$pos[1:10]
dist_mat <-dist(ten_pos, diag = TRUE)

# Check 1-2
ten_pos[1];ten_pos[2]
ten_pos[1]-ten_pos[2] 

# Check 4-10
ten_pos[4];ten_pos[10]
ten_pos[4]-ten_pos[10] 

# Ok looks good. 

# Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)


# Make correlation matrix
############################################

# Transpose so that I'm comparing probe ids, not sample ids
joined_ordered_df_transposed <-joined_ordered_df %>% 
  select(-c("chr", "pos", "strand")) %>% 
  pivot_longer(-probe_id, 'variable', 'value') %>%
  pivot_wider(names_from = probe_id, values_from = value)


# Make sure transposed and distance data are aligned
identical ( joined_ordered_df[1:10,]$probe_id, names(joined_ordered_df_transposed[,2:11]) )
# Must be TRUE

joined_ordered_df[1:10,]
joined_ordered_df_transposed[,2:11]
# They are the same.

# Make correlation matrix
cor_mat <-cor(joined_ordered_df_transposed[,2:11], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]


# Plot it out
############################################
plot(dist_mat_flat, cor_mat_flat)
# Ok
# 

############################################ 
# Repeat with larger set
############################################ 

hundred_pos <-joined_ordered_df$pos[1:100]

dist_mat <-dist(hundred_pos, diag = TRUE)

# Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)

############
# Make correlation matrix
cor_mat <-cor(joined_ordered_df_transposed[,2:101], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]


# Plot it out
############################################
plot(dist_mat_flat, cor_mat_flat)
# Ok

############################################ 
# Repeat with full set
############################################ 

all_pos <-joined_ordered_df$pos

dist_mat <-dist(all_pos, diag = TRUE)

# Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)

############
# Make correlation matrix
cor_mat <-cor(joined_ordered_df_transposed[,-1], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]

###########################
plot(dist_mat_flat, cor_mat_flat)
# Ok
# 
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat, y = cor_mat_flat))+
  geom_bin2d(bins = 100)+
  scale_fill_continuous(type = "viridis")

cor(dist_mat_flat, cor_mat_flat)
summary(lm(dist_mat_flat~cor_mat_flat))
# Really low. Weird.

# Look at absolute values
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat, y = abs(cor_mat_flat)))+
  geom_bin2d(bins = 100)+
  scale_fill_continuous(type = "viridis")

summary(lm(dist_mat_flat~abs(cor_mat_flat)))
# Still really low. 