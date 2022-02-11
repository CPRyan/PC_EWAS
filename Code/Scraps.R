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