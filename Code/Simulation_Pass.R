source("http://www.sthda.com/upload/rquery_cormat.r")

df <-matrix(runif(80), ncol = 8)
colnames(df)<-paste("cg", 1:8, sep = "")
rownames(df)<-paste("id", 1:10, sep = "")

# Correlation matrix
cor_mat <-cor(df)

# Distance
dist_vect <-sample(10, size = 8, replace = TRUE)
dist_vect <-dist_vect[order(dist_vect)]

dist_vect

# pairwise distances
dist_mat <-dist(dist_vect, diag = TRUE)

plot(dist_mat, cor_mat)
dist_mat*cor_mat
# Ok, we got it and we got some decent calculations from this but the cor_mat is full, I want only bottom

cor_mat_flat <-cor_mat[lower.tri(cor_mat, diag = FALSE)]
dist_mat_flat <-c(dist_mat)

plot(dist_mat_flat, cor_mat_flat)
