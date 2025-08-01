library(mashr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(assertthat)

setwd("/data/ADRD/glia_across_NDDs")

########################################

res_list <- readRDS("./analysis/astrocytes/differential_expression/astrocytes_diffexp_results_for_mashr.rds")

########################################

# modified mashr functions to return eigenvalues for genes for different PCs

cov_pca = function(data,npc,subset = NULL) {
  assert_that(npc > 1)
  assert_that(npc <= n_conditions(data))
  if (is.null(subset))
    subset = 1:n_effects(data)
  res.svd = svd(data$Bhat[subset,],nv = npc,nu = npc)
  
  # FIXME: we need to think of for the EE case what to use for svd
  # input: Bhat or Bhat/Shat
  f     = res.svd$v
  Ulist = cov_from_factors(t(f),"PCA")
  d     = diag(res.svd$d[1:npc])
  Ulist = c(Ulist,list("tPCA" = f %*% d^2 %*% t(f)/length(subset)))
  for (i in 1:length(Ulist)) {
    rownames(Ulist[[i]]) <- colnames(data$Bhat)
    colnames(Ulist[[i]]) <- colnames(data$Bhat)
  }
  return(list(Ulist = Ulist, 
              svd = res.svd))
}

n_conditions = function(data){ncol(data$Bhat)}

cov_from_factors = function(f, name){
  Ulist = list()
  for(i in 1:nrow(f)){
    Ulist = c(Ulist,list(r1cov(f[i,])))
  }
  names(Ulist) = paste0(name,"_",(1:nrow(f)))
  return(Ulist)
}

r1cov=function(x){x %*% t(x)}

########################################

# make bhat and pval matrices (bhat = z-scores from nebula, pval = un-adj pval from nebula)

bhat_list_named <- lapply(names(res_list), function(name) {
  df <- res_list[[name]]
  df <- data.frame(RowNames = df$gene, Value = df[,11])
  colnames(df)[2] <- name
  df
})

bhat_df <- Reduce(function(x, y) merge(x, y, by = "RowNames", all = TRUE), bhat_list_named)
rownames(bhat_df) <- bhat_df$RowNames
bhat_df <- bhat_df[, -1] %>% 
  as.matrix()

###

pval_list_named <- lapply(names(res_list), function(name) {
  df <- res_list[[name]]
  df <- data.frame(RowNames = df$gene, Value = df[,4])
  colnames(df)[2] <- name  # Rename column to match dataframe name
  df
})

pval_df <- Reduce(function(x, y) merge(x, y, by = "RowNames", all = TRUE), pval_list_named)
rownames(pval_df) <- pval_df$RowNames
pval_df <- pval_df[, -1] %>% 
  as.matrix()

########################################

# set up & mash canonical covariance matrix

mash_data = mash_set_data(Bhat = bhat_df, pval = pval_df)

U.c = cov_canonical(mash_data)

m.c = mash(mash_data, U.c)

########################################

# set up & mash data-driven covariance matrix

m.1by1 = mash_1by1(mash_data)

strong = get_significant_results(m.1by1, thresh = 1e-4)
strong_subset_genes <- names(strong)

# run cov_pca on ALL PCs, then check how many PCs to use based on variance explained -- only use PCs that explain >5% of variance
U.pca = cov_pca(mash_data, npc = n_conditions(mash_data), subset = strong)
svd = U.pca$svd
plot(svd$d^2/sum(svd$d^2))
t <- as.data.frame(svd$d^2/sum(svd$d^2))


# re-run cov_pca w/ desired number of PCs
U.pca.2 = cov_pca(mash_data, npc = 5, subset = strong)

Ulist = U.pca.2$Ulist

U.ed = cov_ed(mash_data, Ulist, subset = strong)

m.ed = mash(mash_data, c(U.c, U.ed))

########################################

# save outputs

saveRDS(m.ed, file = "./analysis/astrocytes/differential_expression/astrocytes_diffexp_mashr_object.rds")
saveRDS(strong_subset_genes, file = "./analysis/astrocytes/differential_expression/astrocytes_strong_subset_genes.rds")
