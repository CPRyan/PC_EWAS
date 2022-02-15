#' ---
#' title: "Lehne Chr1 region PCs"
#' author: "Calen Patrick Ryan"
#' date: "February 14, 2022"
#' ---

## Load packages
library(visdat)
library(factoextra)
library(tidyverse)
library(readr)


load(file = here::here("Output/Data", "joined_ordered.Rdata"))
vis_miss(joined_ordered_df_transposed)
# Ok I'm missing some for these probes - remove them and my data will be cleaner
# Remember to account for this when you use the distance matrix

remove <-names(which(colSums(is.na(joined_ordered_df_transposed)) > 0))




chr1_pcs <-prcomp(joined_ordered_df_transposed %>% 
                    select(-c(remove, variable)), scale = TRUE)

summary(chr1_pcs)

# Eigenvalues
eig.val <- get_eigenvalue(chr1_pcs)
head(eig.val)

# Results for Variables
res.var <- get_pca_var(chr1_pcs)
head(res.var$coord)       # Coordinates
head(res.var$contrib)       # Contributions to the PCs
head(res.var$cos2)        # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(chr1_pcs)
head(res.ind$coord)        # Coordinates
head(res.ind$contrib )       # Contributions to the PCs
head(res.ind$cos2 )       # Quality of representation 


# Viz
factoextra::fviz_eig(chr1_pcs)


factoextra::fviz_pca_ind(chr1_pcs,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_var(chr1_pcs,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# Ok but I want to look at the MOST informative CpGs for the first 2 dimensions
# 
# Pick top CpGs with highest contributions to Dim 1
dim.1_probes <-res.var$contrib %>%
  as.data.frame() %>% 
  rownames_to_column("probe") %>%
  arrange(desc(Dim.1)) %>% 
  select(probe, Dim.1) %>% 
  top_n(50)
# Pick top CpGs with highest contributions to Dim 2
dim.2_probes <-res.var$contrib %>%
  as.data.frame() %>% 
  rownames_to_column("probe") %>%
  arrange(desc(Dim.2)) %>% 
  select(probe, Dim.2) %>% 
  top_n(50)

dim1.2_union_probes <-union(dim.1_probes$probe, dim.2_probes$probe)


fviz_pca_var(chr1_pcs,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             select.var = list(contrib = 20))


pulled_position <-joined_ordered_df %>% 
  filter(!probe_id %in% c(remove)) %>% 
  pull(pos)

fviz_pca_var(chr1_pcs,
             col.var = (pulled_position-min(pulled_position)), # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             select.var = list(contrib = 50))




# Try with 'variable' sites only
 
vis_miss(var_df_t)
# Ok I'm missing some for these probes - remove them and my data will be cleaner
# Remember to account for this when you use the distance matrix

remove <-names(which(colSums(is.na(var_df_t)) > 0))




chr1_pcs <-prcomp(var_df_t %>% 
                    select(-c(remove, variable)), scale. = TRUE)

summary(chr1_pcs)


# Eigenvalues
eig.val <- get_eigenvalue(chr1_pcs)
head(eig.val)

# Results for Variables
res.var <- get_pca_var(chr1_pcs)
head(res.var$coord)      # Coordinates
head(res.var$contrib)        # Contributions to the PCs
head(res.var$cos2 )         # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(chr1_pcs)
head(res.ind$coord)          # Coordinates
head(res.ind$contrib)        # Contributions to the PCs
head(res.ind$cos2 )      # Quality of representation 


# Viz
factoextra::fviz_eig(chr1_pcs)


factoextra::fviz_pca_ind(chr1_pcs,
                         col.ind = "cos2", # Color by the quality of representation
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE     # Avoid text overlapping
)


fviz_pca_var(chr1_pcs,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
