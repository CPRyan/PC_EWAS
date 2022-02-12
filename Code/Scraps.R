# url <-("https://cumccolumbia-my.sharepoint.com/:u:/r/personal/cpr2139_cumc_columbia_edu/Documents/GSE55763_Lehne/GSE55763_normalized_betas.txt.gz?csf=1&web=1&e=YERZCC")
# 
# readr::read_table(file = here::here("Data", "GSE55763_normalized_betas.txt.gz"), n_max = 10)
# 
# readr::read_table(url, n_max = 10)
# Eventually just downloaded the file as is (gzip) and am running this as is.

short_lehne <-readr::read_table(here::here("Data", "GSE55763_normalized_betas.txt"), n_max = 100, progress = show_progress())
# Proof of reading




# my_locations_probes

short_lehne %>% 
  select(-contains("Detection") & -contains("Pval")) %>% 
  filter(ID_REF %in% my_locations_probes)



short_lehne %>%
  select(-contains("Detection") & -contains("Pval"))  %>% 
  select(1:10)


# zz=gzfile('GSE55763_normalized_betas.txt.gz')  
# dat=read.txt(zz,header=F)

zz=gzfile(here::here("Data", "GSE55763_normalized_betas.txt")) 
read_tsv(file = zz, n_max = 5)
# Ok, I can read the first 5 rows. It has 5423 columns

# Can I only load a subset of the columns? 
zz=gzfile(here::here("Data", "GSE55763_normalized_betas.txt")) 
names(read_tsv(file = zz, n_max = 5) %>% 
        data.frame())


# Try to select the ten names

zz=gzfile(here::here("Data", "GSE55763_normalized_betas.txt")) 
ten_names <-names(read_tsv(file = zz, n_max = 5) %>% 
                    data.frame())[1:10]

zz=gzfile(here::here("Data", "GSE55763_normalized_betas.txt")) 
read_tsv(file = zz, col_select = c(ten_names), n_max = 5)
# Can't do it. 



zz=gzfile(here::here("Data", "GSE55763_normalized_betas.txt")) 
read.delim(zz, nrows = 5)
# Still can't do it using original read.delim

###################################

# Code pillaged to flatten matrix - 
###################################
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
  )
}

ut <-upper.tri(cor_mat) 
data.frame(
  row = rownames(cor_mat),
  colum = rownames(cor_mat),
  cor  =(cor_mat)[ut])