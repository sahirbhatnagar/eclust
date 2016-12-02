#' Not in helper function
#'
#' @description Helper function which does the opposite of "in"

"%ni%" <- Negate("%in%")

#' Simulate Covariates With Exposure Dependent Correlations
#'
#' @description This is a wrapper of the \code{\link[WGCNA]{simulateDatExpr}}
#'   function which simulates data in a modular structure (i.e. in blocks). This
#'   function simulates data in 5 blocks referred to as Turquoise, Blue, Red,
#'   Green and Yellow, separately for exposed (E=1) and unexposed (E=0)
#'   observations.
#'
#' @param n number of observations
#' @param p total number of predictors to simulate
#' @param exposed binary numeric vector of length \code{n} with 0 for unexposed
#'   and 1 for exposed
#' @param rho numeric value representing the expected correlation between green
#'   module and red module
#' @param ... arguments passed to the \code{\link[WGCNA]{simulateDatExpr}} function
#' @return \code{n x p} matrix of simulated data
#' @examples
#' d0 <- simModule(n = 100, p = 1000, rho = 0, exposed = FALSE,
#'                 modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
#'                 minCor = 0.01,
#'                 maxCor = 1,
#'                 corPower = 1,
#'                 propNegativeCor = 0.3,
#'                 backgroundNoise = 0.5,
#'                 signed = FALSE,
#'                 leaveOut = 1:4)
#'
#' d1 <- simModule(n = 100, p = 1000, rho = 0.90, exposed = TRUE,
#'                 modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
#'                 minCor = 0.4,
#'                 maxCor = 1,
#'                 corPower = 0.3,
#'                 propNegativeCor = 0.3,
#'                 backgroundNoise = 0.5,
#'                 signed = FALSE)
#'
#' X <- rbind(d0$datExpr, d1$datExpr) %>%
#'  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
#'  magrittr::set_rownames(paste0("Subject",1:n))
#' dim(X)

simulate_modules <- function(n, p, rho, exposed, ...) {

  if (exposed) {
    #Step 1: simulate the seed module eigengenes
    sMEturquoise <- rnorm(n)

    #expected cor(sMEblue,sMEturquoise) = 0.60
    sMEblue <- 0.60 * sMEturquoise + sqrt(1 - 0.60 ^ 2) * rnorm(n)

    sMEyellow <- rnorm(n)

    sMEgreen <- rnorm(n)

    #expected cor(e.continuous,seed.ME)=0.95
    temp0 <- rho[1] * sMEgreen + sqrt(1 - rho[1] ^ 2) * rnorm(n)

    #expected cor(y.continuous,seed.ME) <- -0.95
    sMEred <- rho[1] * temp0 + sqrt(1 - rho[1] ^ 2) * rnorm(n)

    datsME <- data.frame(sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow)

    dat1 <- WGCNA::simulateDatExpr(eigengenes = datsME, nGenes = p, ...)
  } else {

    #Step 1: simulate the seed module eigengenes
    sMEturquoise <- rnorm(n)

    #expected cor(sMEblue,sMEturquoise) = 0.60
    sMEblue <- 0.60 * sMEturquoise + sqrt(1 - 0.60 ^ 2) * rnorm(n)

    sMEyellow <- rnorm(n)

    sMEgreen <- rnorm(n)

    #expected cor(e.continuous,seed.ME)=0.95
    temp0 <- rho[1] * sMEgreen + sqrt(1 - rho[1] ^ 2) * rnorm(n)

    #expected cor(y.continuous,seed.ME) <- -0.95
    sMEred <- rho[1] * temp0 + sqrt(1 - rho[1] ^ 2) * rnorm(n)

    datsME <- data.frame(sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow)

    dat1 <- WGCNA::simulateDatExpr(eigengenes = datsME, nGenes = p, ...)

  }

  return(dat1)
}



fisherTransform <- function (n1, r1, n2, r2) {
  num1a <- which(r1 >= 0.99)
  num2a <- which(r2 >= 0.99)
  r1[num1a] <- 0.99
  r2[num2a] <- 0.99
  num1b <- which(r1 <= -0.99)
  num2b <- which(r2 <= -0.99)
  r1[num1b] <- -0.99
  r2[num2b] <- -0.99
  # atanh (inverse hyperbolic tangent) simplifies to
  # 0.5 * log(1+r)/log(1-r) , for r < 1
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  dz <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
  pv <- 2 * (1 - pnorm(abs(dz)))
  return(list(diff = dz, pval = pv))
}



#' Calculate Fisher's Z test for correlations
fisherZ <- function(n0, cor0, n1, cor1) {

  # n0 = 50
  # n1 = 50
  # cor0 = corrX0
  # cor1 = corrX1

  # by default this doesnt include the diagonal
  # this collapses the correlation matrix by columns
  ccc0 <- as.vector(cor0[lower.tri(cor0)])
  ccc1 <- as.vector(cor1[lower.tri(cor1)])

  p <- nrow(cor1)

  # number of Z statistics to calculate (p choose 2)
  geneNames <- rownames(cor1)

  zstat <- fisherTransform(n0, ccc0, n1, ccc1)$diff

  # convert vector to symmetric matrix
  zMat <- diag(p)
  zMat[lower.tri(zMat)] <- zstat
  zMat <- zMat + t(zMat) - diag(diag(zMat))
  dimnames(zMat) <- list(geneNames,geneNames)
  class(zMat) <- c("similarity", class(zMat))
  return(zMat)
}


#' Cluster similarity matrix
#'
#' @description Return cluster membership of each predictor. This function is
#'   called internally by the \code{\link{generate_data}} and
#'   \code{\link{generate_data_mars}} function
#'
#' @param x similarity matrix. must have non-NULL dimnames i.e., the rows and
#'   columns should be labelled, e.g. "Gene1, Gene2, ..."
#' @param expr gene expression data (training set). rows are people, columns are
#'   genes
#' @param exprTest gene expression test set
#' @param distanceMethod  one of "euclidean","maximum","manhattan", "canberra",
#'   "binary","minkowski" to be passed to \code{\link[stats]{dist}} function. If
#'   missing, then this function will take 1-x as the dissimilarity measure.
#'   This functionality is for diffCorr,diffTOM, fisherScore matrices which need
#'   to be converted to a distance type matrix.
#' @param clustMethod Cluster the data using hierarchical clustering or
#'   prototype clustering. Defaults \code{clustMethod="hclust"}. Other option is
#'   \code{\link[protoclust]{protoclust}}, however this package must be
#'   installed before proceeding with this option
#' @param cutMethod what method to use to cut the dendrogram. \code{'dynamic'}
#'   refers to \code{\link[dynamicTreeCut]{}} library. \code{'gap'} is
#'   Tibshirani's gap statistic \code{\link[cluster]{clusGap}} using the
#'   \code{'Tibs2001SEmax'} rule. \code{'fixed'} is a fixed number specified by the
#'   \code{nClusters} argument
#' @param nClusters number of clusters. Only used if \code{cutMethod = fixed}
#' @param K.max the maximum number of clusters to consider, must be at least
#'   two. Only used if \code{cutMethod='gap'}
#' @param B integer, number of Monte Carlo (“bootstrap”) samples. Only used if \code{cutMethod='gap'}
#' @param method the agglomeration method to be used. This should be (an
#'   unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'   "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
#'   or "centroid" (= UPGMC).
#' @export
cluster_similarity <- function(x,
                              expr,
                              exprTest,
                              distanceMethod,
                              clustMethod = c("hclust", "protoclust"),
                              cutMethod = c("dynamic","gap", "fixed"),
                              nClusters,
                              method = c("complete", "average", "ward.D2",
                                         "single", "ward.D", "mcquitty",
                                         "median", "centroid"),
                              K.max = 10, B = 50, nPC) {

  # x = corrX ; expr = X
  # exprTest = X[sample(seq_len(nrow(X)),nrow(X), replace = TRUE ),]
  # dim(X) ; dim(expr) ; dim(exprTest)
  # clustMethod = c("hclust")
  # cutMethod = c("dynamic")
  # nClusters = 6
  # method = c("complete")
  # summary = c("pca")
  # K.max = 10; B = 50
  # distance = as.dist(1 - x)

  geneNames <- dimnames(x)[[1]]
  p <- nrow(x)
  method <- match.arg(method)
  cutMethod <- match.arg(cutMethod)
  clustMethod <- match.arg(clustMethod)


  if (clustMethod=="protoclust") {
    if (!requireNamespace("protoclust", quietly = TRUE)) {
      stop("protoclust package needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }

  if (cutMethod=="gap") {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("cluster package needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }


  distance <- if (missing(distanceMethod)) {
    as.dist(1 - x)
  } else dist(x = x, method = distanceMethod)

  hc <- switch(clustMethod,
               hclust = {
                 hclust(distance, method = method)
               },
               protoclust = {
                 protoclust(distance)
               }
  )

  #plot(hc)
  # create cluster function used if Gap statistic is requested
  # its as.dist(x) here because I am passing the
  # 1-x matrix to the cluster::clusGap function
  if (cutMethod == "gap") {

    FUNcluster <- if (missing(distanceMethod)) {
      switch(clustMethod,
             hclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(
                   cutree(
                     hclust(as.dist(xMat), method = method), k = k
                   )
                 )
               })
             },
             protoclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(protoclust::protocut(protoclust(as.dist(xMat)),
                                                 k = k)$cl)})
             }
      )
    } else {
      switch(clustMethod,
             hclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(cutree(hclust(dist(xMat, method = distanceMethod),
                                          method = method), k = k))})
             },
             protoclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(protoclust::protocut(
                   protoclust(dist(xMat, method = distanceMethod)),
                   k = k)$cl)})
             }
      )
    }
    #return(FUNcluster)
  }


  clustAssignment <- switch(cutMethod,
                            dynamic = {
                              if (clustMethod == "hclust") {
                                dynamicTreeCut::cutreeDynamic(
                                  hc,
                                  method = "hybrid",
                                  distM = as.matrix(distance),
                                  #cutHeight = 0.995,
                                  deepSplit = 1,
                                  pamRespectsDendro = T,
                                  minClusterSize = 50)
                              } else {
                                hcMod <- hc
                                class(hcMod) <- "hclust"
                                dynamicTreeCut::cutreeDynamic(
                                  hcMod,
                                  distM = as.matrix(distance),
                                  #cutHeight = 0.995,
                                  deepSplit = 1,
                                  method = "hybrid",
                                  pamRespectsDendro = T,
                                  minClusterSize = 50)
                              }
                            },
                            gap = {
                              if (clustMethod == "hclust") {
                                gapResult <- cluster::clusGap(1 - x,
                                                              FUNcluster = FUNcluster,
                                                              K.max = K.max,
                                                              B = B)
                                nClustGap <- cluster::maxSE(f = gapResult$Tab[, "gap"],
                                                   SE.f = gapResult$Tab[, "SE.sim"],
                                                   method = "Tibs2001SEmax",
                                                   SE.factor = 1)
                                cutree(hc, nClustGap)

                              } else {
                                gapResult <- cluster::clusGap(1 - x,
                                                              FUNcluster = FUNcluster,
                                                              K.max = K.max,
                                                              B = B)
                                nClustGap <- cluster::maxSE(f = gapResult$Tab[, "gap"],
                                                   SE.f = gapResult$Tab[, "SE.sim"],
                                                   method = "Tibs2001SEmax",
                                                   SE.factor = 1)
                                protocut(hc, k = nClustGap)[["cl"]]
                              }
                            },
                            fixed = {
                              if (clustMethod == "hclust") {
                                cutree(hc, nClusters)
                              } else protoclust::protocut(hc, k = nClusters)[["cl"]]
                            }
  )

  # check if all cluster groups are 0 which means no cluster
  # assignment and everyone is in their own group
  # plot(clustAssignment)
  clusters <- data.table(gene = geneNames,
                         cluster = if (all(clustAssignment == 0))
                           1:p else clustAssignment)
  #setkey(clusters, "cluster")

  # convert cluster numbers to colors which define modules
  clusters[, module := WGCNA::labels2colors(cluster)]
  clusters[, table(cluster,module)]


  # note that the align argument acts as follows if equal to "along average"
  # which is the default: it take the correlation between the average expression
  # in a module and the 1st eigenvector in a module and checks if its less
  # than 0, if its less than 0, then the moduleEigengenes function multiplies
  # the 1st eigenvector by -1, else it returns the unmodified 1st eigenvector
  # note that moduleEigengenes function returns the 1st eigenvector which is
  # equivalent to the rotation returned by prcomp, and what is used in
  # predict.prcomp to calculate the actual PCs.
  # to calculate PC's the following are all equivalent:
  # all.equal((expr %*% prcomp.object$rotation)[,1],
  # predict(prcomp.object)[,1],prcomp.object$x[,1])
  #
  # these are equivalent
  # p <- WGCNA::moduleEigengenes(expr = expr[, clusters$gene],
  #                              colors = clusters$module,
  #                              align = "",
  #                              scale = FALSE)
  # l <- prcomp(t(expr[, which(clusters$module %in% "blue")]), scale. = FALSE,
  # center = FALSE)
  #
  # plot(l$rotation[,1,drop=F],p$eigengenes[,"MEblue"])

  # this plots the eigenvector against the average expression
  # to show the effect of the "along average" argument
  # cbind(pp$PC,pp$averageExpr) %>%
  #   mutate(id = 1:n) %>%
  #   gather(type, value, -id) %>%
  #   separate(type, c("type","module")) %>%
  #   spread(type,value) %>%
  #   magrittr::set_colnames(c("id","module","average", "PC")) %>%
  #   ggplot(.,aes(x = average, y = PC)) + geom_point() + facet_grid(~module) +
  #   theme_bw()

  pp <- extractPC(x_train = expr[, clusters$gene],
                  x_test = exprTest[, clusters$gene],
                  colors = clusters$module,
                  scale = TRUE, nPC = nPC)

  # clusters
  # pp %>% names
  # pp$PCTest
  #
  # pp$varExplained
  # pp$averageExpr
  # pp$eigengenes
  # pp$PC


  return(list(clusters = clusters, pcInfo = pp))

}


#' @rdname simulated_data
#' @export
sim_data <- function(n , n0 , p , genes,
                     E, signal_to_noise_ratio = 1,
                     include_interaction = F,
                     beta = NULL) {

  # number of subjects with E=1
  #n1 = n - n0

  # not used for now
  # The coefficients are constructed to have alternating signs and to be exponentially
  # decreasing. (Friedman et al 2010, Journal of Statistical Software)
  # beta <- (-1)^{1:p} * exp(-(2*1:p-1)/100)
  # beta <- (-1)^{1:p} * exp(-(2*1:p - 1 )/600)*sin(1:p)

  #   genes = X
  #   signal_to_noise_ratio = 4
  #   n0 = n1 = 100
  #   E = c(rep(0,n0),rep(1, n1))
  #   beta = c(rnorm(1000),0, rep(0,1000));include_interaction = T
  #   beta = c(rnorm(1000));include_interaction = F

  if (include_interaction) {
    DT <- cbind(genes,E) %>% as.data.table()
    alloc.col(DT,2*p + 1) %>% invisible()
    indx <- grep('Gene', colnames(DT))

    for (j in indx){
      set(DT, i = NULL, j = paste0("Gene",j,":E"), value = DT[[j]]*DT[['E']])
    }
  } else {
    DT <- cbind(genes,E) %>% as.data.table()
  }

  y.star <- {DT %>% as.matrix()} %*% beta
  error <- rnorm(n)
  k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error)))

  y <- y.star + k*error

  result <- if (include_interaction) as.matrix(cbind(y,DT)) else as.matrix(cbind(y,DT))
  colnames(result)[1] <- "Y"
  class(result) <- append(class(result), "expression")

  return(result)
}

#' Generate data and test and training sets for simulation study
#'
#' @description create a function that takes as input, the number of genes, the
#'   true beta vector, the gene expression matrix created from the
#'   generate_blocks function and returns a list of data matrix, as well as
#'   correlation matrices, TOM matrices, cluster information, training and test
#'   data
#' @note this function calls the \code{sim_data} to generate phenotype as a
#'   function of the gene expression data. This function also returns other
#'   information derived from the simulated data including the test and training
#'   sets, the correlation and TOM matrices and the clusters.
#' @note the PCs and averages need to be calculated in the fitting functions,
#'   because these will change based on the CV fold
#' @return list of (in the following order) \describe{ \item{beta_truth}{}
#'   \item{distance}{} \item{DT}{data.table of simulated data from the
#'   \code{sim_data} function} \item{Y}{} \item{X0}{} \item{X1}{}
#'   \item{X_train}{} \item{X_test}{} \item{Y_train}{} \item{Y_test}{}
#'   \item{DT_train}{} \item{DT_test}{} \item{S0}{} \item{n_clusters}{}
#'   \item{clustered_genes_train}{} \item{clustered_genes_test}{}
#'   \item{clusters}{} \item{tom_train_all}{} \item{tom_train_diff}{}
#'   \item{tom_train_e1}{} \item{tom_train_e0}{} \item{corr_train_all}{}
#'   \item{corr_train_diff}{} \item{corr_train_e1}{} \item{corr_train_e0}{}
#'   \item{mse_null}{} }
#'
#' @param p number of genes in design matrix
#' @param X gene expression matrix of size n x p using the
#'   \code{generate_blocks} function
#' @param beta true beta coefficient vector
#' @param n total number of subjects
#' @param n0 total number of subjects with E=0
#' @param signal_to_noise_ratio signal to noise ratio, default is 4
#' @inheritParams cluster_similarity
#' @param cluster_distance character representing which matrix from the training
#'   set that you want to use to cluster the genes. Must be one of the following
#'   \itemize{ \item corr, corr0, corr1, tom, tom0, tom1, diffcorr, difftom,
#'   corScor, tomScor, fisherScore }
#' @param EclustDistance character representing which matrix from the training
#'   set that you want to use to cluster the genes based on the environment.
#'   See \code{cluster_distance} for
#'   avaialble options. Should be different from \code{cluster_distance}. For
#'   example, if \code{cluster_distance=corr} and
#'   \code{EclustDistance=fisherScore}. That is, one should be based on
#'   correlations ignoring the environment, and the other should be based on
#'   correlations accounting for the environment. This function will always
#'   return this add on
#'
#' @examples
#' \dontrun{
#' p = 1000
#' n=200;n0=100
#' beta_genes <- c(runif(50,0.9,1.1),
#'                 runif(50, -1.1,-0.9),
#'                 rep(0,900))
#' # gene expression matrix used in sim_data function
#' X <- mapply(generate_blocks,
#'             rho_E0 = c(-0.70, runif(8, 0.01,0.05), 0.70),
#'             rho_E1 = c(0.70, runif(8, 0.01, 0.05), 0.70),
#'             MoreArgs = list(block_size = 100, n = n, n0 = n0), SIMPLIFY = F) %>%
#'   do.call(cbind, . ) %>%
#'   magrittr::set_colnames(paste0("Gene", 1:1000)) %>%
#'   magrittr::set_rownames(paste0("Subject",1:200))
#'
#' cluster_distance <- "corr"
#' generate_data(p = p, n = n, n0 = n0, X = X, beta_genes = beta_genes, cluster_distance = "corr")
#' }
#' @rdname simulated_data
#' @export

generate_data <- function(p, X, beta,
                          cluster_distance = c("corr", "corr0", "corr1", "tom",
                                               "tom0", "tom1", "diffcorr",
                                               "difftom","corScor", "tomScor",
                                               "fisherScore"),
                          n, n0, include_interaction = F,
                          signal_to_noise_ratio = 1,
                          eclust_distance = c("fisherScore", "corScor", "diffcorr",
                                              "difftom"),
                          cluster_method = c("hclust", "protoclust"),
                          cut_method = c("dynamic","gap", "fixed"),
                          distance_method = c("euclidean","maximum", "manhattan",
                                              "canberra", "binary", "minkowski"),
                          n_clusters,
                          agglomeration_method = c("complete", "average", "ward.D2",
                                                   "single", "ward.D", "mcquitty",
                                                   "median", "centroid"),
                          nPC = 1,
                          K.max = 10, B = 10) {

  # p = p; X = X ; beta = beta
  # n = n; n0 = n0
  # cluster_distance = "corr"
  # include_interaction = F
  # signal_to_noise_ratio = 0.5
  # cluster_method = "hclust" ; cut_method = "dynamic";agglomeration_method="complete";
  # distance_method = "euclidean"
  # eclust_distance = "diffcorr"; nPC = 1


  agglomeration_method <- match.arg(agglomeration_method)
  cut_method <- match.arg(cut_method)
  cluster_method <- match.arg(cluster_method)
  distance_method <- match.arg(distance_method)
  cluster_distance <- match.arg(cluster_distance)
  eclust_distance <- match.arg(eclust_distance)


  names(beta) <- if (include_interaction) {
    c(paste0("Gene",1:p),"E", paste0("Gene",1:p,":E"))
  } else c(paste0("Gene",1:p),"E")

  # total true beta vector: this includes all the betas for the genes, then the
  # environment beta, then their interactions if interaction is true.
  # This is used to calculate the model error. This is the same as beta,
  # but in matrix form
  beta_truth <- as.matrix(beta)

  # Gene names belonging to the active set
  S0 <- names(beta)[which(beta != 0)]

  n1 <- n - n0

  message("Creating data and simulating response")

  DT <- as.data.frame(sim_data(n = n, n0 = n0, p = p, genes = X,
                               include_interaction = include_interaction,
                               E = c(rep(0,n0), rep(1, n1)),
                               beta = beta,
                               signal_to_noise_ratio = signal_to_noise_ratio))
  dim(DT)

  Y <- as.matrix(DT[,"Y"])

  #remove response from X0 and X1
  X0 <- as.matrix(DT[which(DT$E == 0),-1])
  X1 <- as.matrix(DT[which(DT$E == 1),-1])

  # partition-data
  trainIndex <- caret::createDataPartition(DT$E, p = .5, list = FALSE, times = 1)
  DT_train <- DT[trainIndex,]
  DT_test <- DT[-trainIndex,]

  # X_train and X_test contain the environment variable
  X_train <- DT_train[,-1] %>% as.matrix
  Y_train <- DT_train[, 1]
  X_test <- DT_test[,-1] %>% as.matrix
  Y_test <- DT_test[, 1]

  mse_null <- crossprod(mean(Y_test) - Y_test)/length(Y_test)

  # gene expression data
  genes_e0 <- DT_train[which(DT_train$E == 0),paste0("Gene",1:p)] %>% as.matrix
  genes_e1 <- DT_train[which(DT_train$E == 1),paste0("Gene",1:p)] %>% as.matrix
  genes_all <- rbind(genes_e0,genes_e1)

  message("Calculating similarity matrices")

  # gene expression data
  genes_all_test <- DT_test[,paste0("Gene",1:p)] %>% as.matrix

  corr_train_e0 <- WGCNA::cor(genes_e0)
  corr_train_e1 <- WGCNA::cor(genes_e1)
  corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
  corr_train_all <- WGCNA::cor(genes_all)

  tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0)
  dimnames(tom_train_e0)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_e0)[[2]] <- dimnames(corr_train_all)[[2]]

  tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1)
  dimnames(tom_train_e1)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_e1)[[2]] <- dimnames(corr_train_all)[[2]]

  tom_train_diff <- abs(tom_train_e1 - tom_train_e0)
  dimnames(tom_train_diff)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_diff)[[2]] <- dimnames(corr_train_all)[[2]]

  tom_train_all <- WGCNA::TOMsimilarityFromExpr(genes_all)
  dimnames(tom_train_all)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_all)[[2]] <- dimnames(corr_train_all)[[2]]




  # corScor and Fisher Score matrices
  alpha <- 2
  Scorr <- abs(corr_train_e0 + corr_train_e1 - alpha * corr_train_all)
  class(Scorr) <- c("similarity", class(Scorr))

  # Stom <- abs(tom_train_e1 + tom_train_e0 - alpha * tom_train_all)
  # class(Stom) <- c("similarity", class(Stom))

  fisherScore <- fisherZ(n0 = n0, cor0 = corr_train_e0,
                         n1 = n1, cor1 = corr_train_e1)

  # class(tom_train_all) <- append(class(tom_train_all), "similarity")
  # class(tom_train_diff) <- append(class(tom_train_diff), "similarity")
  # class(tom_train_e1) <- append(class(tom_train_e1), "similarity")
  # class(tom_train_e0) <- append(class(tom_train_e0), "similarity")
  class(corr_train_all) <- append(class(corr_train_all), "similarity")
  class(corr_train_diff) <- append(class(corr_train_diff), "similarity")
  class(corr_train_e1) <- append(class(corr_train_e1), "similarity")
  class(corr_train_e0) <- append(class(corr_train_e0), "similarity")

  message("Creating CV folds from training data")

  # Folds for Cross validation
  folds_train <- caret::createFolds(Y_train, k = 10, list = T)
  DT_train_folds <- lapply(folds_train, function(i) DT_train[-i,])
  X_train_folds <- lapply(DT_train_folds, function(i) i[,-grep("Y",colnames(i))])
  Y_train_folds <- lapply(DT_train_folds, function(i) i[,grep("Y",colnames(i))])

  message(sprintf("Calculating number of clusters based on %s using %s with %s
                  linkage and the %s to determine the number of clusters",
                  cluster_distance, cluster_method, agglomeration_method, cut_method))

  # clusters based on cluster_distance argument
  similarity <- switch(cluster_distance,
                       corr = corr_train_all,
                       corr0 = corr_train_e0,
                       corr1 = corr_train_e1,
                       diffcorr = corr_train_diff,
                       difftom = tom_train_diff,
                       tom0 = tom_train_e0,
                       tom1 = tom_train_e1,
                       tom = tom_train_all,
                       corScor = Scorr,
                       tomScor = Stom,
                       fisherScore = fisherScore)

  # results for clustering, PCs and averages for each block
  # the only difference here is the distance_method arg
  res <- if (cluster_distance %in% c("diffcorr","difftom",
                                     "corScor", "tomScor","fisherScore")) {
    cluster_similarity(x = similarity,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      distanceMethod = distance_method,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  } else {
    cluster_similarity(x = similarity,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  }

  message(paste("Calculating number of environment clusters based on ",
                eclust_distance))

  # clusters based on eclust_distance
  similarityEclust <- switch(eclust_distance,
                             corr = corr_train_all,
                             corr0 = corr_train_e0,
                             corr1 = corr_train_e1,
                             diffcorr = corr_train_diff,
                             difftom = tom_train_diff,
                             tom0 = tom_train_e0,
                             tom1 = tom_train_e1,
                             tom = tom_train_all,
                             corScor = Scorr,
                             tomScor = Stom,
                             fisherScore = fisherScore)


  resEclust <- if (eclust_distance %in% c("diffcorr","difftom",
                                          "corScor", "tomScor","fisherScore")) {
    cluster_similarity(x = similarityEclust,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      distanceMethod = distance_method,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  } else {
    cluster_similarity(x = similarityEclust,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  }


  # we need to combine the cluster information here
  # this is based on cluster_distance only
  clustersAll <- copy(res$clusters)
  n_clusters_All <- res$pcInfo$nclusters

  message(sprintf("There are %d clusters derived from the %s similarity matrix",
                  n_clusters_All, cluster_distance))

  # this is based on eclust_distance only
  n_clusters_Eclust <- resEclust$pcInfo$nclusters
  clustersEclust <- copy(resEclust$clusters)

  message(sprintf("There are %d clusters derived from the %s environment similarity matrix",
                  n_clusters_Eclust, eclust_distance))

  # this is based on both
  n_clusters_Addon <- n_clusters_All + n_clusters_Eclust

  message(sprintf("There are a total of %d clusters derived from the %s
                  similarity matrix and the %s environment similarity matrix",
                  n_clusters_Addon,cluster_distance,eclust_distance))

  # check if any of the cluster numbers in clustersEclust are 0
  # if there are, then add n_clusters+1 to each module number in
  # clustersEclust, else just add n_clusters. this is to control for the
  # situation where there are some clusters numbers of 0 which would cause
  # identical cluster numbers in the clusters and clustersEclust data
  if (clustersEclust[,any(cluster==0)]) {
    clustersEclust[,cluster := cluster + n_clusters_All + 1 ]
  } else {
    clustersEclust[,cluster := cluster + n_clusters_All ]
  }

  # this contains the clusters from the cluster_distance (e.g. corr matrix)
  # and the clusters from the eclust_distance (e.g. fisherScore)
  clustersAddon <- rbindlist(list(clustersAll, clustersEclust))

  # need to calculate penalty factors for group lasso
  # I put all main effects and interactions of a given module in the same group
  # and the size of the penalty factor is sqrt(size of module), where the
  # size of the module includes both main and interaction effects
  # environment should get penalized, in the original simulation 1
  # it was not being penalized which is maybe why it was performing well
  if (include_interaction) {

    gene_groups = copy(clustersAll)
    gene_groups[, gene := paste0(gene,":E")]
    gene_groups <- rbind(clustersAll,gene_groups) %>% setkey(cluster)

    pf_temp <- gene_groups[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)

    gene_groups_inter <- rbind(pf_temp[gene_groups],
                               data.table(cluster = n_clusters_All, N = 1,
                                          pf = 1, gene = "E", module = "empty"))
    # gglasso needs groups number consecutively 1, 2,3 ...
    gene_groups_inter[, cluster:=cluster+1]
    setkey(gene_groups_inter, cluster)

    gene_groups_Addon = copy(clustersAddon)
    gene_groups_Addon[, gene := paste0(gene,":E")]
    gene_groups_Addon <- rbind(clustersAddon, gene_groups_Addon) %>% setkey(cluster)

    pf_temp_Addon <- gene_groups_Addon[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)

    gene_groups_inter_Addon <- rbind(pf_temp_Addon[gene_groups_Addon],
                                     data.table(cluster = n_clusters_Addon, N = 1,
                                                pf = 1, gene = "E", module = "empty"))
    # gglasso needs groups number consecutively 1, 2,3 ...
    gene_groups_inter_Addon[, cluster:=cluster+1]
    setkey(gene_groups_inter_Addon, cluster)
  }

  DT <- DT %>% as.matrix
  class(DT) <- append(class(DT),"eset")

  result <- list(beta_truth = beta_truth,
                 similarity = similarity,
                 similarityEclust = similarityEclust,
                 DT = DT,
                 Y = Y, X0 = X0, X1 = X1, X_train = X_train, X_test = X_test,
                 Y_train = Y_train, Y_test = Y_test, DT_train = DT_train,
                 DT_test = DT_test, S0 = S0,
                 n_clusters_All = n_clusters_All,
                 n_clusters_Eclust = n_clusters_Eclust,
                 n_clusters_Addon = n_clusters_Addon,
                 clustersAll = clustersAll,
                 clustersAddon = clustersAddon,
                 clustersEclust = clustersEclust,
                 gene_groups_inter = if (include_interaction) gene_groups_inter else NULL,
                 gene_groups_inter_Addon = if (include_interaction) gene_groups_inter_Addon else NULL,
                 tom_train_all = tom_train_all, tom_train_diff = tom_train_diff,
                 tom_train_e1 = tom_train_e1,tom_train_e0 = tom_train_e0,
                 corr_train_all = corr_train_all,
                 corr_train_diff = corr_train_diff,
                 corr_train_e1 = corr_train_e1, corr_train_e0 = corr_train_e0,
                 fisherScore = fisherScore,
                 corScor = Scorr,
                 # corTom = Stom,
                 mse_null = mse_null, DT_train_folds = DT_train_folds,
                 X_train_folds = X_train_folds, Y_train_folds = Y_train_folds)
  return(result)
}



prepare_data <- function(data, response = "Y", exposure = "E", probe_names) {

  # data = cbind(pcTrain, Y = Y[trainIndex], E = E[trainIndex])


  # ===========================================================

  # Check for sensible dataset
  ## Make sure you have response, exposure.
  if (!(response %in% colnames(data))) stop(sprintf("response argument specified as %s but this column not found in 'data' data.frame", response))
  if (!(exposure %in% colnames(data))) stop(sprintf("exposure argument specified as %s but this column not found in 'data' data.frame", exposure))
  if (!missing(probe_names)) {
    if (!(probe_names %in% colnames(data))) stop(sprintf("probe_names argument specified as %s but this column not found in 'data' data.frame", probe_names))
  }

  # if missing main_effect_names, assume everything except response and exposure
  # are main effects
  if (missing(probe_names)) {
    probe_names <- setdiff(colnames(data), c(response, exposure))
  }

  # rename response to be Y and exposure to be E
  colnames(data)[which(colnames(data) == response)] <- "Y"
  colnames(data)[which(colnames(data) == exposure)] <- "E"

  x_mat <- model.matrix(as.formula(paste0("~(", paste0(probe_names, collapse="+"), ")*E - 1")), data = data)


  # reformulate(paste0("~(", paste0(colnames(pcTrain)[1:5], collapse="+"), ")*E"))


  interaction_names <- grep(":", colnames(x_mat), value = T)
  main_effect_names <- setdiff(colnames(x_mat), interaction_names)

  return(list(X = x_mat, Y = data[["Y"]], E = data[["E"]],
              main_effect_names = main_effect_names,
              interaction_names = interaction_names))
  # x_mat


}



#' This is a modified version of firstPC which was actually giving the first 2 PCs
#' with no other option. This function is more flexible and the nPC argument is used.
#' currently only nPC = 1 and nPC = 2 are supported
extractPC <- function(x_train, colors, x_test,
                      y_train, y_test,
                      impute = TRUE, nPC,
                      excludeGrey = FALSE,
                      grey = if (is.numeric(colors)) 0 else "grey",
                      subHubs = TRUE, trapErrors = FALSE,
                      returnValidOnly = trapErrors, softPower = 6, scale = TRUE,
                      verbose = 0, indent = 0) {


  # x_train = result[["X_train"]] ; x_test = result[["X_test"]];
  # x_train_mod <- x_train %>% as.data.frame
  # x_test_mod = x_test %>% as.data.frame
  # gene_groups = result[["clustersAll"]]
  # x_train = x_train_mod[,gene_groups$gene];
  # colors = gene_groups$cluster;
  # x_test = x_test_mod[,gene_groups$gene]
  # impute = TRUE; nPC = 2; align = "along average";
  # excludeGrey = FALSE; grey = if (is.numeric(colors)) 0 else "grey";
  # subHubs = TRUE; trapErrors = FALSE; returnValidOnly = trapErrors;
  # softPower = 6; scale = TRUE; verbose = 0; indent = 0;

  if (is.null(x_train)) {
    stop("moduleEigengenes: Error: x_train is NULL. ")
  }
  if (is.null(colors)) {
    stop("moduleEigengenes: Error: colors is NULL. ")
  }
  if (is.null(dim(x_train)) || length(dim(x_train)) != 2)
    stop("moduleEigengenes: Error: x_train must be two-dimensional.")
  if (dim(x_train)[2] != length(colors))
    stop("moduleEigengenes: Error: ncol(x_train) and length(colors) must be equal (one color per gene).")
  if (is.factor(colors)) {
    nl = nlevels(colors)
    nlDrop = nlevels(colors[, drop = TRUE])
    if (nl > nlDrop)
      stop(paste("Argument 'colors' contains unused levels (empty modules). ",
                 "Use colors[, drop=TRUE] to get rid of them."))
  }

  # maxVarExplained = 10
  # if (nPC > maxVarExplained)
  #   warning(paste("Given nPC is too large. Will use value",
  #                 maxVarExplained))
  #
  # nVarExplained = min(nPC, maxVarExplained)

  modlevels = levels(factor(colors))

  if (excludeGrey)
    if (sum(as.character(modlevels) != as.character(grey)) >
        0) {
      modlevels = modlevels[as.character(modlevels) !=
                              as.character(grey)]
    } else {
      stop(paste("Color levels are empty. Possible reason: the only color is grey",
                 "and grey module is excluded from the calculation."))
    }

  # these are the loadings aka the first and second eigenvector for each module
  # length of these vectors will vary depending on the size of the module
  eigenVectors <- vector("list", nPC*length(modlevels))

  # these are the actual PC's aka the data %*% eigenvector
  # each column will be a n-dimensional vector.. i.e. a value for each person
  #  this will contain the first 2 PCs for each module
  PC <- data.frame(matrix(NA, nrow = dim(x_train)[[1]],
                          ncol = nPC*length(modlevels)))

  PCTest <- data.frame(matrix(NA, nrow = dim(x_test)[[1]],
                              ncol = nPC*length(modlevels)))

  #   PLS <- data.frame(matrix(NA, nrow = dim(x_train)[[1]],
  #                           ncol = nPC*length(modlevels)))
  #
  #   PLSTest <- data.frame(matrix(NA, nrow = dim(x_test)[[1]],
  #                               ncol = nPC*length(modlevels)))

  # list to store prcomp objects
  prcompObj <- vector("list", length(modlevels))

  # list to store pls objects
  # plsObj <- vector("list", length(modlevels))


  # this is the average expression in a module for each subject
  # so this is a n x length(modlevels) matrix
  averExpr <- data.frame(matrix(NA, nrow = dim(x_train)[[1]],
                                ncol = length(modlevels)))

  averExprTest <- data.frame(matrix(NA, nrow = dim(x_test)[[1]],
                                    ncol = length(modlevels)))

  varExpl <- vector("double", nPC*length(modlevels))

  # validMEs = rep(TRUE, length(modlevels))
  # validAEs = rep(FALSE, length(modlevels))

  # these are the means and sds used for subsequent predictions
  means = vector("list", length(modlevels))
  sds = vector("list", length(modlevels))

  # isPC = rep(TRUE, length(modlevels))
  # isHub = rep(FALSE, length(modlevels))
  validColors = colors

  # names(eigenVectors) = paste(moduleColor.getMEprefix(), modlevels,
  #                          sep = "")
  names(PC) = paste(rep(paste0("pc",seq_len(nPC)), length(modlevels)),
                    rep(modlevels, each = nPC), sep = "_")
  names(averExpr) = paste("avg", modlevels, sep = "")
  #   names(PCTest) = paste(rep(paste0("pc",seq_len(nPC)), length(modlevels)),
  #                         rep(modlevels, each = nPC), sep = "_")
  names(averExprTest) = paste("avg", modlevels, sep = "")

  for (i in seq_len(length(modlevels))) {
    # i=2
    if (verbose > 1)
      printFlush(paste("moduleEigengenes : Working on ME for module",
                       modlevels[i]))
    modulename = modlevels[i]
    restrict1 = as.character(colors) == as.character(modulename)
    if (verbose > 2)
      printFlush(paste(spaces, " ...", sum(restrict1),
                       "genes"))

    datModule <- as.matrix(x_train[, restrict1])
    datModuleTest <- as.matrix(x_test[, restrict1])

    # xy_train <- data.frame(Y = as.matrix(y_train), x_train[, restrict1])
    # xy_test <- data.frame(Y = as.matrix(y_test), x_test[, restrict1])

    # dim(datModule)
    # dim(t(datModule))
    # dim(x_train)

    # using prcomp first (need to use untransposed data!)
    prcompObj[[i]] <- prcomp(datModule, center = scale, scale. = scale)

    # plsObj[[i]] <- pls::plsr(Y ~ ., ncomp = nPC, data = xy_train, validation = "CV")

    # plot(prcompObj[[i]])
    # View(stats:::prcomp.default)
    # prcompObj[[i]]$x %>% dim
    # prcompObj[[i]] %>% names
    # prcompObj[[i]]$rotation %>% dim

    if (nPC == 1) {

      eigenVectors[[i]] <- prcompObj[[i]]$rotation[,1, drop = F]
      averExpr[,i] <- rowMeans(datModule, na.rm = TRUE)
      averExprTest[,i] <- rowMeans(datModuleTest, na.rm = TRUE)

      varExpl[[i]] <- factoextra::get_eigenvalue(prcompObj[[i]])[1,"variance.percent"]

      # corAve = cor(averExpr[,i], prcompObj[[i]]$rotation[,1],
      #              use = "p")
      # if (!is.finite(corAve)) corAve = 0
      # if (corAve < 0) prcompObj[[i]]$rotation[,1] = -prcompObj[[i]]$rotation[,1]

      PC[, i] <- predict(prcompObj[[i]])[,1]
      PCTest[, i] <- predict(prcompObj[[i]], newdata = datModuleTest)[,1]

      # PLS[, i] <- predict(plsObj[[i]], ncomp = nPC, type = "scores")
      # PLSTest[, i] <- predict(plsObj[[i]], ncomp = nPC, type = "scores", newdata = xy_test)[,1]

    } else if (nPC == 2) {
      eigenVectors[[2*i-1]] <- prcompObj[[i]]$rotation[,1, drop = F]
      eigenVectors[[2*i]] <- prcompObj[[i]]$rotation[,2, drop = F]
      averExpr[,i] <- rowMeans(datModule, na.rm = TRUE)
      averExprTest[,i] <- rowMeans(datModuleTest, na.rm = TRUE)

      varExpl[[2*i-1]] <- factoextra::get_eigenvalue(prcompObj[[i]])[1,"variance.percent"]
      varExpl[[2*i]] <- factoextra::get_eigenvalue(prcompObj[[i]])[2,"variance.percent"]
      # corAve = cor(averExpr[,i], prcompObj[[i]]$rotation[,1],
      #              use = "p")
      # if (!is.finite(corAve)) corAve = 0
      # if (corAve < 0) prcompObj[[i]]$rotation[,1] = -prcompObj[[i]]$rotation[,1]

      PC[, 2*i-1] <- predict(prcompObj[[i]])[,1]
      PC[, 2*i] <- predict(prcompObj[[i]])[,2]
      # PCTest[, 2*i-1] <- predict(prcompObj[[i]], newdata = datModuleTest)[,1]
      # PCTest[, 2*i] <- predict(prcompObj[[i]], newdata = datModuleTest)[,2]

      # plot(PC[, i], prcompObj[[i]]$x[,1])
      #means[i] <- prcompObj[[i]]$center
      #sds[i] <- prcompObj[[i]]$scale
    }

  }

  list(eigengenes = eigenVectors, averageExpr = averExpr,
       averageExprTest = averExprTest,
       varExplained = varExpl, validColors = validColors,
       PC = PC, PCTest = PCTest, prcompObj = prcompObj,
       # PLS = PLS, PLSTest = PLSTest,
       nclusters = length(modlevels))
}


# obj is an earth object
get.used.pred.names <- function(obj) {
  any1 <- function(x) any(x != 0) # like any but no warning if x is double
  names(which(apply(obj$dirs[obj$selected.terms, , drop=FALSE], 2, any1)))
}


#' @rdname simulated_data
#' @export
sim_data_mars <- function(n , n0 , p , genes,
                          E, signal_to_noise_ratio = 1,
                          truemodule, nActive) {

  # nActive
  # truemodule = truemodule1
  # genes = X
  # E = c(rep(0,n0),rep(1, n1))

  DT <- cbind(genes,E) %>% as.data.table()

  x1 <- genes[,which(truemodule %in% 3)[1:(nActive/2)]]
  u1 <- svd(x1)$u[,1]

  x2 <- genes[,which(truemodule %in% 4)[1:(nActive/2)]]
  u2 <- svd(x2)$u[,1]

  y.star <- 0.1*(u1 + u2 + E) + 4 * pmax(u1-0.01, 0) * pmax(u2-0.05, 0) * E
  error <- rnorm(n)
  k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error)))

  y <- y.star + k*error

  result <- as.matrix(cbind(y,DT))
  colnames(result)[1] <- "Y"
  class(result) <- append(class(result), "expression")

  return(result)
}


generate_data_mars <- function(p, X, beta,
                               truemodule,
                               nActive,
                               cluster_distance = c("corr", "corr0", "corr1", "tom",
                                                    "tom0", "tom1", "diffcorr",
                                                    "difftom","corScor", "tomScor",
                                                    "fisherScore"),
                               n, n0, include_interaction = F,
                               signal_to_noise_ratio = 1,
                               eclust_distance = c("fisherScore", "corScor", "diffcorr",
                                                   "difftom"),
                               cluster_method = c("hclust", "protoclust"),
                               cut_method = c("dynamic","gap", "fixed"),
                               distance_method = c("euclidean","maximum", "manhattan",
                                                   "canberra", "binary", "minkowski"),
                               n_clusters,
                               agglomeration_method = c("complete", "average", "ward.D2",
                                                        "single", "ward.D", "mcquitty",
                                                        "median", "centroid"),
                               nPC = 1,
                               K.max = 10, B = 10) {

  # p = p; X = X ; beta = beta
  # n = n; n0 = n0
  # cluster_distance = "corr"
  # include_interaction = F
  # signal_to_noise_ratio = 0.5
  # cluster_method = "hclust" ; cut_method = "dynamic";agglomeration_method="complete";
  # distance_method = "euclidean"
  # eclust_distance = "diffcorr"; nPC = 1


  agglomeration_method <- match.arg(agglomeration_method)
  cut_method <- match.arg(cut_method)
  cluster_method <- match.arg(cluster_method)
  distance_method <- match.arg(distance_method)
  cluster_distance <- match.arg(cluster_distance)
  eclust_distance <- match.arg(eclust_distance)


  names(beta) <- if (include_interaction) {
    c(paste0("Gene",1:p),"E", paste0("Gene",1:p,":E"))
  } else c(paste0("Gene",1:p),"E")

  # total true beta vector: this includes all the betas for the genes, then the
  # environment beta, then their interactions if interaction is true.
  # This is used to calculate the model error. This is the same as beta,
  # but in matrix form
  beta_truth <- as.matrix(beta)

  # Gene names belonging to the active set
  S0 <- names(beta)[which(beta != 0)]

  n1 <- n - n0

  message("Creating data and simulating response for MARS model")

  DT <- as.data.frame(sim_data_mars(n = n, n0 = n0, p = p, genes = X,
                                    truemodule = truemodule,
                                    nActive = nActive,
                                    E = c(rep(0,n0), rep(1, n1)),
                                    signal_to_noise_ratio = signal_to_noise_ratio))
  dim(DT)

  Y <- as.matrix(DT[,"Y"])

  #remove response from X0 and X1
  X0 <- as.matrix(DT[which(DT$E == 0),-1])
  X1 <- as.matrix(DT[which(DT$E == 1),-1])

  # partition-data
  trainIndex <- caret::createDataPartition(DT$E, p = .5, list = FALSE, times = 1)
  DT_train <- DT[trainIndex,]
  DT_test <- DT[-trainIndex,]

  # X_train and X_test contain the environment variable
  X_train <- DT_train[,-1] %>% as.matrix
  Y_train <- DT_train[, 1]
  X_test <- DT_test[,-1] %>% as.matrix
  Y_test <- DT_test[, 1]

  mse_null <- crossprod(mean(Y_test) - Y_test)/length(Y_test)

  # gene expression data
  genes_e0 <- DT_train[which(DT_train$E == 0),paste0("Gene",1:p)] %>% as.matrix
  genes_e1 <- DT_train[which(DT_train$E == 1),paste0("Gene",1:p)] %>% as.matrix
  genes_all <- rbind(genes_e0,genes_e1)

  message("Calculating similarity matrices")

  # gene expression data
  genes_all_test <- DT_test[,paste0("Gene",1:p)] %>% as.matrix

  corr_train_e0 <- WGCNA::cor(genes_e0)
  corr_train_e1 <- WGCNA::cor(genes_e1)
  corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
  corr_train_all <- WGCNA::cor(genes_all)

  tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0)
  dimnames(tom_train_e0)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_e0)[[2]] <- dimnames(corr_train_all)[[2]]

  tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1)
  dimnames(tom_train_e1)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_e1)[[2]] <- dimnames(corr_train_all)[[2]]

  tom_train_diff <- abs(tom_train_e1 - tom_train_e0)
  dimnames(tom_train_diff)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_diff)[[2]] <- dimnames(corr_train_all)[[2]]

  tom_train_all <- WGCNA::TOMsimilarityFromExpr(genes_all)
  dimnames(tom_train_all)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_all)[[2]] <- dimnames(corr_train_all)[[2]]




  # corScor and Fisher Score matrices
  alpha <- 2
  Scorr <- abs(corr_train_e0 + corr_train_e1 - alpha * corr_train_all)
  class(Scorr) <- c("similarity", class(Scorr))

  # Stom <- abs(tom_train_e1 + tom_train_e0 - alpha * tom_train_all)
  # class(Stom) <- c("similarity", class(Stom))

  fisherScore <- fisherZ(n0 = n0, cor0 = corr_train_e0,
                         n1 = n1, cor1 = corr_train_e1)

  # class(tom_train_all) <- append(class(tom_train_all), "similarity")
  # class(tom_train_diff) <- append(class(tom_train_diff), "similarity")
  # class(tom_train_e1) <- append(class(tom_train_e1), "similarity")
  # class(tom_train_e0) <- append(class(tom_train_e0), "similarity")
  class(corr_train_all) <- append(class(corr_train_all), "similarity")
  class(corr_train_diff) <- append(class(corr_train_diff), "similarity")
  class(corr_train_e1) <- append(class(corr_train_e1), "similarity")
  class(corr_train_e0) <- append(class(corr_train_e0), "similarity")

  message("Creating CV folds from training data")

  # Folds for Cross validation
  folds_train <- caret::createFolds(Y_train, k = 10, list = T)
  DT_train_folds <- lapply(folds_train, function(i) DT_train[-i,])
  X_train_folds <- lapply(DT_train_folds, function(i) i[,-grep("Y",colnames(i))])
  Y_train_folds <- lapply(DT_train_folds, function(i) i[,grep("Y",colnames(i))])

  message(sprintf("Calculating number of clusters based on %s using %s with %s
                  linkage and the %s to determine the number of clusters",
                  cluster_distance, cluster_method, agglomeration_method, cut_method))

  # clusters based on cluster_distance argument
  similarity <- switch(cluster_distance,
                       corr = corr_train_all,
                       corr0 = corr_train_e0,
                       corr1 = corr_train_e1,
                       diffcorr = corr_train_diff,
                       difftom = tom_train_diff,
                       tom0 = tom_train_e0,
                       tom1 = tom_train_e1,
                       tom = tom_train_all,
                       corScor = Scorr,
                       tomScor = Stom,
                       fisherScore = fisherScore)

  # results for clustering, PCs and averages for each block
  # the only difference here is the distance_method arg
  res <- if (cluster_distance %in% c("diffcorr","difftom",
                                     "corScor", "tomScor","fisherScore")) {
    cluster_similarity(x = similarity,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      distanceMethod = distance_method,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  } else {
    cluster_similarity(x = similarity,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  }

  message(paste("Calculating number of environment clusters based on ",
                eclust_distance))

  # clusters based on eclust_distance
  similarityEclust <- switch(eclust_distance,
                             corr = corr_train_all,
                             corr0 = corr_train_e0,
                             corr1 = corr_train_e1,
                             diffcorr = corr_train_diff,
                             difftom = tom_train_diff,
                             tom0 = tom_train_e0,
                             tom1 = tom_train_e1,
                             tom = tom_train_all,
                             corScor = Scorr,
                             tomScor = Stom,
                             fisherScore = fisherScore)


  resEclust <- if (eclust_distance %in% c("diffcorr","difftom",
                                          "corScor", "tomScor","fisherScore")) {
    cluster_similarity(x = similarityEclust,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      distanceMethod = distance_method,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  } else {
    cluster_similarity(x = similarityEclust,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  }


  # we need to combine the cluster information here
  # this is based on cluster_distance only
  clustersAll <- copy(res$clusters)
  n_clusters_All <- res$pcInfo$nclusters

  message(sprintf("There are %d clusters derived from the %s similarity matrix",
                  n_clusters_All, cluster_distance))

  # this is based on eclust_distance only
  n_clusters_Eclust <- resEclust$pcInfo$nclusters
  clustersEclust <- copy(resEclust$clusters)

  message(sprintf("There are %d clusters derived from the %s environment similarity matrix",
                  n_clusters_Eclust, eclust_distance))

  # this is based on both
  n_clusters_Addon <- n_clusters_All + n_clusters_Eclust

  message(sprintf("There are a total of %d clusters derived from the %s
                  similarity matrix and the %s environment similarity matrix",
                  n_clusters_Addon,cluster_distance,eclust_distance))

  # check if any of the cluster numbers in clustersEclust are 0
  # if there are, then add n_clusters+1 to each module number in
  # clustersEclust, else just add n_clusters. this is to control for the
  # situation where there are some clusters numbers of 0 which would cause
  # identical cluster numbers in the clusters and clustersEclust data
  if (clustersEclust[,any(cluster==0)]) {
    clustersEclust[,cluster := cluster + n_clusters_All + 1 ]
  } else {
    clustersEclust[,cluster := cluster + n_clusters_All ]
  }

  # this contains the clusters from the cluster_distance (e.g. corr matrix)
  # and the clusters from the eclust_distance (e.g. fisherScore)
  clustersAddon <- rbindlist(list(clustersAll, clustersEclust))

  # need to calculate penalty factors for group lasso
  # I put all main effects and interactions of a given module in the same group
  # and the size of the penalty factor is sqrt(size of module), where the
  # size of the module includes both main and interaction effects
  # environment should get penalized, in the original simulation 1
  # it was not being penalized which is maybe why it was performing well
  if (include_interaction) {

    gene_groups = copy(clustersAll)
    gene_groups[, gene := paste0(gene,":E")]
    gene_groups <- rbind(clustersAll,gene_groups) %>% setkey(cluster)

    pf_temp <- gene_groups[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)

    gene_groups_inter <- rbind(pf_temp[gene_groups],
                               data.table(cluster = n_clusters_All, N = 1,
                                          pf = 1, gene = "E", module = "empty"))
    # gglasso needs groups number consecutively 1, 2,3 ...
    gene_groups_inter[, cluster:=cluster+1]
    setkey(gene_groups_inter, cluster)

    gene_groups_Addon = copy(clustersAddon)
    gene_groups_Addon[, gene := paste0(gene,":E")]
    gene_groups_Addon <- rbind(clustersAddon, gene_groups_Addon) %>% setkey(cluster)

    pf_temp_Addon <- gene_groups_Addon[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)

    gene_groups_inter_Addon <- rbind(pf_temp_Addon[gene_groups_Addon],
                                     data.table(cluster = n_clusters_Addon, N = 1,
                                                pf = 1, gene = "E", module = "empty"))
    # gglasso needs groups number consecutively 1, 2,3 ...
    gene_groups_inter_Addon[, cluster:=cluster+1]
    setkey(gene_groups_inter_Addon, cluster)
  }

  DT <- DT %>% as.matrix
  class(DT) <- append(class(DT),"eset")

  result <- list(beta_truth = beta_truth,
                 similarity = similarity,
                 similarityEclust = similarityEclust,
                 DT = DT,
                 Y = Y, X0 = X0, X1 = X1, X_train = X_train, X_test = X_test,
                 Y_train = Y_train, Y_test = Y_test, DT_train = DT_train,
                 DT_test = DT_test, S0 = S0,
                 n_clusters_All = n_clusters_All,
                 n_clusters_Eclust = n_clusters_Eclust,
                 n_clusters_Addon = n_clusters_Addon,
                 clustersAll = clustersAll,
                 clustersAddon = clustersAddon,
                 clustersEclust = clustersEclust,
                 gene_groups_inter = if (include_interaction) gene_groups_inter else NULL,
                 gene_groups_inter_Addon = if (include_interaction) gene_groups_inter_Addon else NULL,
                 tom_train_all = tom_train_all, tom_train_diff = tom_train_diff,
                 tom_train_e1 = tom_train_e1,tom_train_e0 = tom_train_e0,
                 corr_train_all = corr_train_all,
                 corr_train_diff = corr_train_diff,
                 corr_train_e1 = corr_train_e1, corr_train_e0 = corr_train_e0,
                 fisherScore = fisherScore,
                 corScor = Scorr,
                 # corTom = Stom,
                 mse_null = mse_null, DT_train_folds = DT_train_folds,
                 X_train_folds = X_train_folds, Y_train_folds = Y_train_folds)
  return(result)
}



#' @param ... arguments passed to the \code{cluster_similarity} function
#' @param distance_method argument passed to the dist function for calculating the distance
#' for the clusters based on the corr,corr1,corr0, tom, tom0, tom1 matrices
cluster_kmeans <- function(data,
                           response,
                           exposure,
                           train_index,
                           test_index,
                           min_cluster_size = 50,
                           cluster_distance = c("corr", "corr0", "corr1", "tom",
                                                "tom0", "tom1", "diffcorr",
                                                "difftom", "fisherScore"),
                           eclust_add = TRUE,
                           eclust_distance = c("fisherScore", "corScor", "diffcorr",
                                               "difftom"),
                           distance_method = c("euclidean","maximum","manhattan", "canberra",
                                               "binary","minkowski"), ...) {

  # data <- clusters[[1]]
  #   data <- t(placentaALL[filterd_probes[1:500],])
  #   min_cluster_size <- 50
  #   exposure <- E
  #   train_index <- trainIndex
  #   test_index <- testIndex
  #   nPC <- 2
  #   cluster_distance <- "tom"
  #   cluster_method <- "hclust"
  #   cut_method <- "dynamic"
  #   agglomeration_method <- "average"
  #   distance_method <- "euclidean"
  #   eclust_add <- TRUE
  #   eclust_distance <- "difftom"
  #   cut_method <- "dynamic" #,"gap", "fixed")
  #   # summary <- match.arg(summary)
  #   response <- Y

  # ==================================================================
  args <- list(...)
  xtrain <- data[train_index,]
  xtest <- data[test_index,]
  etrain <- exposure[train_index]

  ytrain <- response[train_index]
  ytest <- response[test_index]

  xtrainE0 <- xtrain[which(etrain == 0), ]
  xtrainE1 <- xtrain[which(etrain == 1), ]

  corrTrain <- WGCNA::cor(xtrain)

  #   message(sprintf("Calculating number of clusters based on %s using %s with %s
  #                   linkage and the %s method to determine the number of clusters",
  #                   cluster_distance, cluster_method, agglomeration_method, cut_method))

  #############################################################################
  #               CORRELATION/TOM CLUSTERS                                    #
  #############################################################################

  # this cannot be based on difftom, diffcorr, fisherScore
  if (cluster_distance %in% c("difftom", "diffcorr", "fisherScore")) stop(message("cluster_distance must be one of corr, corr0, corr1, tom, tom0, tom1"))

  # clusters based on cluster_distance argument
  similarity <- switch(cluster_distance,
                       corr = corrTrain,
                       corr0 = WGCNA::cor(xtrainE0),
                       corr1 = WGCNA::cor(xtrainE1),
                       tom0 = {
                         tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                         dimnames(tomTrainE0)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainE0)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainE0
                       },
                       tom1 = {
                         tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                         dimnames(tomTrainE1)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainE1)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainE1
                       },
                       tom = {
                         tomTrainAll <- WGCNA::TOMsimilarityFromExpr(xtrain)
                         dimnames(tomTrainAll)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainAll)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainAll
                       })

  # results for clustering
  # note that the cluster_similarity returns the PCs but you dont need this info
  # because it is calculated in the clust_fun fitting function
  # we just need to provide the clust_fun function the group membership for all the data
  res <- cluster_similarity(x = similarity,
                            x_train = xtrain,
                            x_test = xtest,
                            y_train = ytrain,
                            y_test = ytest,
                            distance = dist(x = similarity, method = distance_method), ...)

  #############################################################################
  #               ECLUST CLUSTERS                                             #
  #############################################################################

  message(paste("Calculating number of environment clusters based on",eclust_distance))

  # clusters based on eclust_distance
  similarityEclust <- switch(eclust_distance,
                             corr = corrTrain,
                             corr0 = WGCNA::cor(xtrainE0),
                             corr1 = WGCNA::cor(xtrainE1),
                             diffcorr = {
                               corr0 <- WGCNA::cor(xtrainE0)
                               corr1 <- WGCNA::cor(xtrainE1)
                               abs(corr1 - corr0)
                             },
                             difftom = {
                               tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                               tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                               tomTrainDiff <- abs(tomTrainE1 - tomTrainE0)
                               dimnames(tomTrainDiff)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainDiff)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainDiff
                             },
                             tom0 = {
                               tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                               dimnames(tomTrainE0)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainE0)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainE0
                             },
                             tom1 = {
                               tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                               dimnames(tomTrainE1)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainE1)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainE1
                             },
                             tom = {
                               tomTrainAll <- WGCNA::TOMsimilarityFromExpr(xtrain)
                               dimnames(tomTrainAll)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainAll)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainAll
                             },
                             fisherScore = {
                               n0 <- nrow(xtrainE0)
                               n1 <- nrow(xtrainE1)
                               corr0 <- WGCNA::cor(xtrainE0)
                               corr1 <- WGCNA::cor(xtrainE1)
                               fisherZ(n0 = n0, cor0 = corr0,
                                       n1 = n1, cor1 = corr1)
                             })

  resEclust <- if (eclust_distance %in% c("diffcorr","difftom","fisherScore")) {
    cluster_similarity(x = similarityEclust,
                       x_train = xtrain,
                       x_test = xtest,
                       y_train = ytrain,
                       y_test = ytest, ...)
  } else {
    cluster_similarity(x = similarityEclust,
                       x_train = xtrain,
                       x_test = xtest,
                       y_train = ytrain,
                       y_test = ytest,
                       distance = dist(x = similarity, method = distance_method), ...)
  }

  # we need to combine the cluster information here
  # this is based on cluster_distance only
  clustersAll <- copy(res$clusters)
  n_clusters_All <- res$pcInfo$nclusters

  message(sprintf("There are %d clusters derived from the %s similarity matrix",
                  n_clusters_All, cluster_distance))

  # this is based on eclust_distance only
  n_clusters_Eclust <- resEclust$pcInfo$nclusters
  clustersEclust <- copy(resEclust$clusters)

  message(sprintf("There are %d clusters derived from the %s environment similarity matrix",
                  n_clusters_Eclust, eclust_distance))

  # this is based on both
  n_clusters_Addon <- n_clusters_All + n_clusters_Eclust

  message(sprintf("There are a total of %d clusters derived from the %s
                  similarity matrix and the %s environment similarity matrix",
                  n_clusters_Addon,cluster_distance,eclust_distance))

  # check if any of the cluster numbers in clustersEclust are 0
  # if there are, then add n_clusters+1 to each module number in
  # clustersEclust, else just add n_clusters. this is to control for the
  # situation where there are some clusters numbers of 0 which would cause
  # identical cluster numbers in the clusters and clustersEclust data
  if (clustersEclust[,any(cluster==0)]) {
    clustersEclust[,cluster := cluster + n_clusters_All + 1 ]
  } else {
    clustersEclust[,cluster := cluster + n_clusters_All ]
  }

  # this contains the clusters from the cluster_distance (e.g. TOM matrix)
  # and the clusters from the eclust_distance (e.g. diffTOM)
  clustersAddon <- rbindlist(list(clustersAll, clustersEclust))
  # clustersAddon[, table(cluster, module)]

  # these are only derived on the main effects genes.. E is only included in the model
  # this is the clusters based on tom and difftom
  PC_and_avg_Addon <- extractPC(x_train = xtrain[,clustersAddon$gene],
                                colors = clustersAddon$cluster,
                                x_test = xtest[,clustersAddon$gene],
                                y_train = ytrain,
                                y_test = ytest, nPC = args$nPC)

  # this is the clusters based on tom only
  PC_and_avg_All <- extractPC(x_train = xtrain[,clustersAll$gene],
                              colors = clustersAll$cluster,
                              x_test = xtest[,clustersAll$gene],
                              y_train = ytrain,
                              y_test = ytest, nPC = args$nPC)

  # n.clusters <- PC_and_avg_Addon$nclusters

  # this contains either the averages or PCs for each module in a data.frame
  #   clust_data_Addon <- switch(summary,
  #                        avg = PC_and_avg_Addon$averageExpr,
  #                        pc = PC_and_avg_Addon$PC)
  #
  #   clust_data_All <- switch(summary,
  #                            avg = PC_and_avg_All$averageExpr,
  #                            pc = PC_and_avg_All$PC)


  return(list(clustersAddon = PC_and_avg_Addon,
              clustersAll = PC_and_avg_All,
              etrain = etrain))
}

