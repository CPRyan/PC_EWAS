#' ---
#' title: "Lehne Cor and Dist matrices"
#' author: "Calen Patrick Ryan"
#' date: "February 14, 2022"
#' ---

# Extract those and setup as necessary for my correlation matrix. 
library(tidyverse)
library(readr)


#' Load data
############################################
my_locations <-read_csv(here::here("Output/Data", "my_locations.csv"))
my_chr_1_dat <-read_csv(here::here("Output/Data", "my_chr_1_dat.csv"))

joined_ordered_df <-inner_join(my_locations, my_chr_1_dat, by = c("probe_id" = "ID_REF"))

#' Make distance matrix (10 observations)
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
ten_pos[9];ten_pos[10]
ten_pos[9]-ten_pos[10] 

# Ok looks good. 

##' Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)


##' Make correlation matrix
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

##' Make correlation matrix
cor_mat <-cor(joined_ordered_df_transposed[,2:11], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]


##' Plot it out
############################################
plot(dist_mat_flat, cor_mat_flat)
# Ok
# 

############################################ 
#' Repeat with larger set
############################################ 

hundred_pos <-joined_ordered_df$pos[1:100]

dist_mat <-dist(hundred_pos, diag = TRUE)

##' Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)

############
##' Make correlation matrix
cor_mat <-cor(joined_ordered_df_transposed[,2:101], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]


##' Plot it out
############################################
plot(dist_mat_flat, cor_mat_flat)
# Ok

############################################ 
#' Repeat with full set
############################################ 

all_pos <-joined_ordered_df$pos

dist_mat <-dist(all_pos, diag = TRUE)

##' Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)

############
##' Make correlation matrix
cor_mat <-cor(joined_ordered_df_transposed[,-1], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]

###########################
plot(dist_mat_flat, cor_mat_flat)
# Ok

##' Look at plot (correlation both - and +)
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat, y = cor_mat_flat))+
  geom_bin2d(bins = 100)+
  scale_fill_continuous(type = "viridis")

cor(dist_mat_flat, cor_mat_flat)
summary(lm(dist_mat_flat~cor_mat_flat))
# Really low. Weird.

##' Look at absolute values
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat, y = abs(cor_mat_flat)))+
  geom_bin2d(bins = 400)+
 # geom_smooth()+
  scale_fill_continuous(type = "viridis")


#' Cumulative density plot of distance
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat)) +
  stat_ecdf(geom = "step")

summary(lm(dist_mat_flat~abs(cor_mat_flat)))
# Still really low. 
 
summary(lm(dist_mat_flat[dist_mat_flat < 1e6]~abs(cor_mat_flat[dist_mat_flat < 1e6])))
summary(lm(dist_mat_flat[dist_mat_flat > 1e6]~abs(cor_mat_flat[dist_mat_flat > 1e6])))

######################################
#' What about methylated sites that show variability?
######################################

##' Calculate variability

# my_chr_1_dat <-read_csv(here::here("Output/Data", "my_chr_1_dat.csv"))
# Should already be loaded

top90 <-apply(my_chr_1_dat[,-1], 1, function(x) quantile(x, probs=.90, na.rm = TRUE))
bottom10 <-apply(my_chr_1_dat[,-1], 1, function(x) quantile(x, probs=.10, na.rm = TRUE))

five_perc <-top90-bottom10 < 0.05

myvar_sites <-my_chr_1_dat[c(five_perc),]
myvar_rows <-c(five_perc)


var_df <-joined_ordered_df[joined_ordered_df$probe_id %in% myvar_sites$ID_REF,]

dist_mat <-dist(var_df$pos, diag = TRUE)


##' Flatten
###########
# I'd like to name these at some point, but not for now. 
dist_mat_flat <-c(dist_mat)


##' Make correlation matrix
############################################

# Transpose so that I'm comparing probe ids, not sample ids
# joined_ordered_df_transposed <-joined_ordered_df %>% 
#   select(-c("chr", "pos", "strand")) %>% 
#   pivot_longer(-probe_id, 'variable', 'value') %>%
#   pivot_wider(names_from = probe_id, values_from = value)
# Should be calculated already

##' Filter to the columns with the names that match the probe ids that are variable...
var_df
var_df_t <-joined_ordered_df_transposed[,c(names(joined_ordered_df_transposed) %in% c('variable', myvar_sites$ID_REF))]

# Make sure transposed and distance data are aligned
identical (var_df$probe_id, names(var_df_t)[-1])
# Must be TRUE

##' Make correlation matrix
cor_mat <-cor(var_df_t[,-1], 
              use = "pairwise.complete.obs")

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]

##' Look at plot of absolute values
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat, y = abs(cor_mat_flat)))+
  geom_bin2d(bins = 400)+
  # geom_smooth()+
  scale_fill_continuous(type = "viridis")


##' Cumulative density plot of distance
cbind(dist_mat_flat, cor_mat_flat) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = dist_mat_flat)) +
  stat_ecdf(geom = "step")
#######################################

#' Look at stats of variable DNAm 

# # Take y -axis of cdf and use as weight


# Create ecdf function with your data:
  
fun.ecdf <- ecdf(dist_mat_flat) # x is a vector of your data

# Now use this "ecdf function" to generate the cumulative probabilities of any vector you feed it, including your original, sorted data:
  
my.ecdf <- fun.ecdf(dist_mat_flat)

summary(lm(dist_mat_flat~abs(cor_mat_flat)))


summary(lm(dist_mat_flat~abs(cor_mat_flat), weights = 1/my.ecdf))

#################################################
# Save files for PCs
#################################################
save(joined_ordered_df, 
     joined_ordered_df_transposed,
     var_df, 
     var_df_t, file = here::here("Output/Data", "joined_ordered.Rdata"))
