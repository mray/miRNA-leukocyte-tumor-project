# import mRNA leukocyte reference matrix
sig_matrix = read.delim("LM22.txt", header=T, row.names=1)

# import mRNA mixture file
mixture_file_all = read.delim("OV_HT_HG-U133A.txt.txt", header=T, row.names=1)

# remove samples in which >75% of values are NA
mixture_file = mixture_file_all[ , colSums(is.na(mixture_file_all)) < 0.75*dim(mixture_file_all)[1] ]

# replace "NA" with zero
mixture_file[is.na(mixture_file)] = 0

# transpose mixture file so that samples are rows and genes are columns
mixture_file = t(mixture_file)
mixture_file = as.data.frame(mixture_file)

# remove non-tumor samples
mixture_file = mixture_file[- grep(".11", rownames(mixture_file)),]

# standardize sample names
rownames(mixture_file) = gsub("[.]","",rownames(mixture_file))
rownames(mixture_file) = strtrim(rownames(mixture_file), 12)

# import miR mixture file
mir_df_all = read.delim("OV_miRNA_HiSeq_gene_t.txt", header=TRUE, row.names = 1)

# remove miRNAs in which >75% of values are NA
mir_df = mir_df_all[ , colSums(is.na(mir_df_all)) < 0.75*dim(mir_df_all)[1] ] 

# replace "NA" with zero
mir_df[is.na(mir_df)] = 0

# remove non-tumor samples
mir_df = mir_df[- grep("-11", rownames(mir_df)),]

# standardize sample name
rownames(mir_df) = gsub("[-]","",rownames(mir_df))
rownames(mir_df) = strtrim(rownames(mir_df), 12)

# subset mRNA and miRNA files using matched IDs
mixture_subset_df = mixture_file[rownames(mixture_file) %in% rownames(mir_df), ]
mir_subset_df = mir_df[rownames(mir_df) %in% rownames(mixture_subset_df), ]

# order
mixture_subset_df = mixture_subset_df[order(rownames(mixture_subset_df)),]
mir_subset_df = mir_subset_df[order(rownames(mir_subset_df)),]

# export subset files
write.table(mixture_subset_df, "OV_mRNA_mixture_subset.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(mir_subset_df, "OV_miR_mixture_subset.txt", sep="\t", row.names = TRUE, col.names = TRUE)

# transpose for CIBERSORT
mixture_subset_df_t = t(mixture_subset_df)
mixture_subset_df_t = as.data.frame(mixture_subset_df_t)

mir_subset_df_t = t(mir_subset_df)
mir_subset_df_t = as.data.frame(mir_subset_df_t)


### mRNA CIBERSORT ###

# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the ouput
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt

#dependencies
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm
CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#do permutations
doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE, absolute=FALSE, abs_method='sig.score'){
  
  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  #read in data
  X <- sig_matrix
  Y <- mixture_file
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  #if(max(Y) < 50) 
  Y <- 2^Y
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) out <- c(out, sum(w))
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  write.table(rbind(header,output), file="OV_mRNA_CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
  else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
  obj
}

results_mRNA <- CIBERSORT(sig_matrix, mixture_subset_df_t, perm=0, QN=FALSE, absolute=TRUE, abs_method='no.sumto1')

# SIGNATURE MATRIX CONSTRUCTION #

# import libraries
library(stringdist)

# import CIBERSORT results
cibersort.df = read.delim("OV_mRNA_CIBERSORT-Results.txt", header=TRUE)

# remove samples in which >75% of samples = 0
i.cibersort.df = cibersort.df[,colSums(cibersort.df != 0) > 0.75*dim(cibersort.df)[1]]

# delete last four columns
i.cibersort.df = i.cibersort.df[-c(12:dim(i.cibersort.df)[2])]

# create A and B matrices
A.matrix.df = i.cibersort.df
rownames(A.matrix.df) = A.matrix.df$Mixture
A.matrix.df = A.matrix.df[-c(1)]
B.matrix.df = mir_subset_df[rownames(mir_subset_df) %in% rownames(A.matrix.df), ]

A.matrix.final.df = A.matrix.df[rownames(A.matrix.df) %in% rownames(B.matrix.df), ]
B.matrix.final.df = B.matrix.df[rownames(B.matrix.df) %in% rownames(A.matrix.df), ]

# order by IDs
A.matrix.final.df = A.matrix.final.df[order(rownames(A.matrix.final.df)),]
B.matrix.final.df = B.matrix.final.df[order(rownames(B.matrix.final.df)),]

A.matrix.final.df[is.na(A.matrix.final.df)] = 0
B.matrix.final.df[is.na(B.matrix.final.df)] = 0

A.matrix = as.matrix(A.matrix.final.df)
B.matrix = as.matrix(B.matrix.final.df)

write.table(A.matrix, "OV_A.matrix.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(B.matrix, "OV_B.matrix.txt", sep="\t", row.names = TRUE, col.names = TRUE)


# Unconstrained least-squares fit using lm
model = lm(B.matrix ~ A.matrix + 0)
signature_matrix = model$coefficients
cat('Dimensions of signature matrix:\n')
print(dim(signature_matrix))

# Should have A matrix * Signature Matrix = B matrix
B.matrix_reconstructed = A.matrix %*% signature_matrix

# Calculate reconstruction difference
reconstruction_difference = norm(B.matrix_reconstructed - B.matrix, 'F')
cat('Norm of difference between reconstructed B.matrix and real matrix:\n')
print(reconstruction_difference)

# Initialize "empty" signature matrix
signature_matrix_nonnegative = matrix(0, nrow=dim(signature_matrix)[1], ncol=dim(signature_matrix)[2])

# Optim with non-negative constraint 
for (col_index in 1:dim(signature_matrix_nonnegative)[2]) {
  initial_vector = rep(1, dim(signature_matrix_nonnegative)[1])
  matrix_product_norm = function(vector, col_index) {
    reconstructed = A.matrix %*% vector
    norm(B.matrix[,col_index] - reconstructed, "F")
  }
  
  positivity_lower_bound = rep(0, dim(signature_matrix_nonnegative)[1])
  
  optim_result = optim(initial_vector, matrix_product_norm, col_index=col_index, lower=positivity_lower_bound, method='L-BFGS-B')
  signature_vector = optim_result$par
  signature_matrix_nonnegative[, col_index] = signature_vector
}


# Should have A matrix * Signature Matrix = B matrix
B.matrix_nonnegative_reconstructed = A.matrix %*% signature_matrix_nonnegative
# Should be a slightly worse reconstruction than using the unconstrained signature matrix
# since allowing negative values was as good a reconstruction as possible
# and adding constraints will make the reconstruction worse, but hopefully not by much
reconstruction_difference_nonnegative = norm(B.matrix_nonnegative_reconstructed - B.matrix, 'F')
cat('Norm of difference between non-negative reconstructed B.matrix and real matrix:\n')
print(reconstruction_difference_nonnegative)

# The value of these norms isn't meaningful on its own, 
# interpret the amount of increase in constrained compared to unconstrained
increase_proportion = (reconstruction_difference_nonnegative - reconstruction_difference) / reconstruction_difference
cat('Reconstruction error proportional increase:\n')
print(increase_proportion)

# Export unconstrained signature matrix as a tab-deliminated table
write.table(signature_matrix, "OV_signature_matrix_unconstrained.txt", sep="\t", row.names = TRUE, col.names = TRUE)

# Add row and column names to constrained signature matrix
row.names(signature_matrix_nonnegative) = row.names(signature_matrix)
colnames(signature_matrix_nonnegative) = colnames(signature_matrix)

# Transpose
t_sig_matrix_nonneg = t(signature_matrix_nonnegative)

# Export signature matrix nonnegative as a tab-deliminated table
write.table(signature_matrix_nonnegative, "OV_signature_matrix_nonnegative.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(t_sig_matrix_nonneg, "OV_signature_matrix_nonnegative_t.txt", sep="\t", row.names = TRUE, col.names = TRUE)


# REFINE SIGNATURE MATRIX USING DIFFERENTIAL EXPRESSION ANALYSIS #

i.sig.df = as.data.frame(t_sig_matrix_nonneg)

# create row means function
r.means = function(df,col.to.eliminate){
  apply(df[-c(col.to.eliminate)],1,mean)
}

# create row sds function
r.sds = function(df,col.to.eliminate){
  apply(df[-c(col.to.eliminate)],1,sd)
}

# initialize results table
z.scores = matrix(0, nrow=dim(i.sig.df)[1], ncol=dim(i.sig.df)[2])

# initialize vector of number of immune cell types
vector = seq(1,dim(i.sig.df)[2])

# 2-tailed pnorm results for all immune cell/miRNAs
for (i in 1:dim(i.sig.df)[2]) {
  rmeans_i = r.means(i.sig.df,vector[i])
  rsds_i = r.sds(i.sig.df,vector[i])
  for (j in 1:dim(i.sig.df)[1]) {
    z = (i.sig.df[j,i] - rmeans_i[j]) / rsds_i[j]
    z.scores[j,i] = z
  }
}

T.full.probability.results = 2 * pnorm(abs(z.scores), lower.tail = FALSE)

# name rows and columns in pnorm results matrix
colnames(T.full.probability.results) = colnames(i.sig.df)
rownames(T.full.probability.results) = rownames(i.sig.df)

# export tab-delim tables
write.table(T.full.probability.results, "OV_pnorm_results.txt", sep="\t", row.names = TRUE, col.names = TRUE)

# compile miR names+each immune cell type into a list
pList <- list()
for (i in 1:dim(T.full.probability.results)[2]) {
  p = T.full.probability.results[,i,drop = FALSE]
  pList[[i]] = p
}

# take out miRs = 0
nonzero_pList = lapply(pList,function(x) x[ x != 0,,drop = FALSE])

# order list in decending order of miR expression
ordered.pList = lapply(nonzero_pList, function(x) x[order(x, decreasing = FALSE),,drop = FALSE])

# output ranked lists of miRs for each immune cell type
for (i in 1:dim(T.full.probability.results)[2]){
  ranked = do.call(rbind.data.frame, ordered.pList[i])
  filename = paste("20180531_OV_", "ranked_", i, ".txt", sep="")
  write.table(ranked, filename, sep="\t", row.names = TRUE, col.names = TRUE)
}

# select top (n) miRs from each list
# create complete list of top (n) from each immune cell type
# subset full matrix by unique miR names from complete list
# calculate condition number

# initialize condition number results matrix
cond.num.results = NULL

for (i in 10:100){
  full.miR.list = lapply(ordered.pList, function(x) rownames(head(x,i)))
  full.miR.list = unlist(full.miR.list)
  full.miR.list = unique(full.miR.list)
  subset.sig.df = t(subset(t(i.sig.df), select = full.miR.list))
  cond.num = kappa(subset.sig.df, exact = TRUE)
  cond.num.results = rbind(cond.num.results, cond.num)
}

rownames(cond.num.results) = 1:nrow(cond.num.results)
iter.vector = c(10:100)
cond.num.results = cbind(cond.num.results,iter.vector)
colnames(cond.num.results) = c("cond.num","iteration")

df.cond.num.results = as.data.frame(cond.num.results)
min.cond.num = df.cond.num.results$iteration[df.cond.num.results$cond.num==min(df.cond.num.results[,1])]
min.cond.num

ordered.cond.num.results = df.cond.num.results[ order(df.cond.num.results[,1]), ]

write.table(ordered.cond.num.results, "OV_ordered.cond.num.results.txt", sep="\t", row.names = TRUE, col.names = TRUE)

# selected matrix with minimum condition number
full.miR.list = lapply(ordered.pList, function(x) rownames(head(x,min.cond.num)))
full.miR.list = unlist(full.miR.list)
full.miR.list = unique(full.miR.list)
final.subset.sig.df = t(subset(t(i.sig.df), select = full.miR.list))
cond.num = kappa(final.subset.sig.df, exact = TRUE)
cond.num

write.table(final.subset.sig.df, "OV_final.sig.matrix.txt", sep="\t", row.names = TRUE, col.names = TRUE)

### miR CIBERSORT ###

# import signature matrix
mir_sig_matrix = final.subset.sig.df
names(mir_sig_matrix) = gsub("A.matrix","",names(mir_sig_matrix))

# import miRNA mixture file
mir_mixture_file = mir_subset_df_t 

# CIBERSORT R script v1.03 (last updated 07-10-2015)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#
#       Options:
#       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii) QN = Quantile normalization of input mixture (default = TRUE)
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


#dependencies
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm
CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#do permutations
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  
  #read in data
  X <- sig_matrix
  Y <- mixture_file
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  #print(nulldist)
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  write.table(rbind(header,output), file="OV_miRNA_CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

results_OV_miRNA <- CIBERSORT(mir_sig_matrix, mir_mixture_file, perm=0, QN=FALSE)

# COMBINE mRNA and miRNA LEUKOCYTE SIGNATURES #

sig.matrix.intersect.cols = intersect(names(LM22), names(mir_sig_matrix))
combined.sig.matrix = rbind(LM22[,sig.matrix.intersect.cols], mir_sig_matrix[,sig.matrix.intersect.cols])

# import mRNA and miRNA mixtures

mRNA.mixture.file = mixture_subset_df_t
miRNA.mixture.file = mir_subset_df_t

mRNA.mixture.file[is.na(mRNA.mixture.file)] = 0
miRNA.mixture.file[is.na(miRNA.mixture.file)] = 0

# combine mixture files
intersect.cols = intersect(names(mRNA.mixture.file), names(miRNA.mixture.file))
merged.mixture.file = rbind(mRNA.mixture.file[,intersect.cols], miRNA.mixture.file[,intersect.cols])

# CIBERSORT R script v1.03 (last updated 07-10-2015)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#
#       Options:
#       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii) QN = Quantile normalization of input mixture (default = TRUE)
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


#dependencies
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm
CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#do permutations
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  
  #read in data
  X <- sig_matrix
  Y <- mixture_file
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  #print(nulldist)
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  write.table(rbind(header,output), file="OV_combined_CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

combined_OV.results <- CIBERSORT(combined.sig.matrix, merged.mixture.file, perm=0, QN=FALSE)

### mRNA SURVIVAL ANALYSIS ###

# subset LM22 based on column names of mir_sig_matrix
subset_LM22 = LM22[ ,colnames(LM22) %in% colnames(mir_sig_matrix)]

# transpose mRNA mixture for CIBERSORT
mixture_subset_df_t = as.data.frame(t(mixture_subset_df))

# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the ouput
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt

#dependencies
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm
CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#do permutations
doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE, absolute=FALSE, abs_method='sig.score'){
  
  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  #read in data
  X <- sig_matrix
  Y <- mixture_file
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  #if(max(Y) < 50) 
  Y <- 2^Y
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) out <- c(out, sum(w))
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  write.table(rbind(header,output), file="OV_survival_mRNA_CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
  else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
  obj
}

results_survival <- CIBERSORT(subset_LM22, mixture_subset_df_t, perm=0, QN=FALSE, absolute=FALSE)

# mRNA UNIVARIATE SURVIVAL ANALYSIS #

library(broom)
library(survival)
library(stringdist)

# import survival data
survival_df = read.delim("xena_OV_survival.txt", header=T)

# import deconvolution results
cibersort_df = read.delim("OV_survival_mRNA_CIBERSORT-Results.txt", header=T)

# convert time and status columns to numeric
survival_df$time = as.numeric(survival_df$time)
survival_df$status = as.numeric(survival_df$status)

# convert ID column to strings
colnames(cibersort_df)[1] = c("ID")
cibersort_df$ID = as.character(cibersort_df$ID)
survival_df$ID = as.character(survival_df$ID)

# remove punctuation in ID columns
survival_df$ID = gsub("[-]","",survival_df$ID)
cibersort_df$ID = gsub("[.]","",cibersort_df$ID)

# partial match function to match labels
# returns position of matched label in stringVector
ClosestMatch = function(string, stringVector){
  amatch(string, stringVector, method = "osa", weight = c(d=1, i=1, s=1, t=1), maxDist = 0.0001)
}  

# create new "matched" column of matched labels in survival_df
survival_df$ID.in.cibersort_df <- sapply(survival_df$ID, function(x) cibersort_df[ClosestMatch(x, cibersort_df$ID), 'ID'])

# merge two dataframes by matched labels
merged_df = merge(survival_df, cibersort_df, by.x = "ID.in.cibersort_df", by.y = "ID", all.x = TRUE, sort = TRUE)
intersect.names = intersect(survival_df$ID, cibersort_df$ID)

# cut unmatched rows for final dataset
merged_df = merged_df[1:length(intersect.names),]
merged_df = merged_df[,2:length(merged_df)]

# initialize results table
coxph_results_table = NULL

# loop through immune cell type for univariate coxph regression and add each result to results table
for (i in names(merged_df)[4:13]) {
  coxph_result = tidy(coxph(Surv(time, status) ~ merged_df[[i]], data = merged_df))
  coxph_results_table = rbind(coxph_results_table, coxph_result)
}

# add names and rearrange columns
coxph_row.names = as.data.frame(names(merged_df)[4:13])
coxph_results_table_names = cbind(coxph_results_table, coxph_row.names)
coxph_results_table_names = coxph_results_table_names[-c(1)]
coxph_results_table_names = coxph_results_table_names[c(7,1:6)]

write.table(coxph_results_table_names, "OV_mRNA_univariate_coxph_results.txt", sep="\t", row.names = FALSE, col.names = TRUE)

# miRNA UNIVARIATE SURVIVAL ANALYSIS #

# import deconvolution results
cibersort_df = read.delim("OV_miRNA_CIBERSORT-Results.txt", header=T)
colnames(cibersort_df)[1] = c("ID")
cibersort_df$ID = as.character(cibersort_df$ID)
cibersort_df$ID = gsub("[.]","",cibersort_df$ID)

# partial match function to match labels
# returns position of matched label in stringVector
ClosestMatch = function(string, stringVector){
  amatch(string, stringVector, method = "osa", weight = c(d=1, i=1, s=1, t=1), maxDist = 0.0001)
}  

# create new "matched" column of matched labels in survival_df
survival_df$ID.in.cibersort_df <- sapply(survival_df$ID, function(x) cibersort_df[ClosestMatch(x, cibersort_df$ID), 'ID'])

# merge two dataframes by matched labels
merged_df = merge(survival_df, cibersort_df, by.x = "ID.in.cibersort_df", by.y = "ID", all.x = TRUE, sort = TRUE)
intersect.names = intersect(survival_df$ID, cibersort_df$ID)

# cut unmatched rows for final dataset
merged_df = merged_df[1:length(intersect.names),]
merged_df = merged_df[,2:length(merged_df)]

# initialize results table
coxph_results_table = NULL

# loop through immune cell type for univariate coxph regression and add each result to results table
for (i in names(merged_df)[4:13]) {
  coxph_result = tidy(coxph(Surv(time, status) ~ merged_df[[i]], data = merged_df))
  coxph_results_table = rbind(coxph_results_table, coxph_result)
}

# add names and rearrange columns
coxph_row.names = as.data.frame(names(merged_df)[4:13])
coxph_results_table_names = cbind(coxph_results_table, coxph_row.names)
coxph_results_table_names = coxph_results_table_names[-c(1)]
coxph_results_table_names = coxph_results_table_names[c(7,1:6)]

write.table(coxph_results_table_names, "OV_miRNA_univariate_coxph_results.txt", sep="\t", row.names = FALSE, col.names = TRUE)

# mRNA-miRNA UNIVARIATESURVIVAL ANALYSIS #

# import deconvolution results
cibersort_df = read.delim("OV_combined_CIBERSORT-Results.txt", header=T)
colnames(cibersort_df)[1] = c("ID")
cibersort_df$ID = as.character(cibersort_df$ID)
cibersort_df$ID = gsub("[.]","",cibersort_df$ID)

# partial match function to match labels
# returns position of matched label in stringVector
ClosestMatch = function(string, stringVector){
  amatch(string, stringVector, method = "osa", weight = c(d=1, i=1, s=1, t=1), maxDist = 0.0001)
}  

# create new "matched" column of matched labels in survival_df
survival_df$ID.in.cibersort_df <- sapply(survival_df$ID, function(x) cibersort_df[ClosestMatch(x, cibersort_df$ID), 'ID'])

# merge two dataframes by matched labels
merged_df = merge(survival_df, cibersort_df, by.x = "ID.in.cibersort_df", by.y = "ID", all.x = TRUE, sort = TRUE)
intersect.names = intersect(survival_df$ID, cibersort_df$ID)

# cut unmatched rows for final dataset
merged_df = merged_df[1:length(intersect.names),]
merged_df = merged_df[,2:length(merged_df)]

# initialize results table
coxph_results_table = NULL

# loop through immune cell type for univariate coxph regression and add each result to results table
for (i in names(merged_df)[4:13]) {
  coxph_result = tidy(coxph(Surv(time, status) ~ merged_df[[i]], data = merged_df))
  coxph_results_table = rbind(coxph_results_table, coxph_result)
}

# add names and rearrange columns
coxph_row.names = as.data.frame(names(merged_df)[4:13])
coxph_results_table_names = cbind(coxph_results_table, coxph_row.names)
coxph_results_table_names = coxph_results_table_names[-c(1)]
coxph_results_table_names = coxph_results_table_names[c(7,1:6)]

write.table(coxph_results_table_names, "OV_combined_univariate_coxph_results.txt", sep="\t", row.names = FALSE, col.names = TRUE)

