#' Fit Penalized Regression Models on Simulated Cluster Summaries
#'
#' @description This function creates summaries of the given clusters (e.g. 1st
#'   PC or average), and then fits a penalized regression model on those
#'   summaries. To be used with simulated data where the 'truth' is known i.e.,
#'   you know which features are associated with the response. This function was
#'   used to produce the simulation results in Bhatnagar et al. 2016. Can run
#'   lasso, elasticnet, SCAD or MCP models
#' @param x_train \code{ntrain x p} matrix of simulated training set where
#'   \code{ntrain} is the number of training observations  and \code{p} is total
#'   number of predictors. This matrix needs to have named columns representing
#'   the feature names or the gene names
#' @param x_test \code{ntest x p} matrix of simulated training set where
#'   \code{ntest} is the number of training observations  and \code{p} is total
#'   number of predictors. This matrix needs to have named columns representing
#'   the feature names or the gene names
#' @param y_train numeric vector of length \code{ntrain} representing the
#'   responses for the training subjects. If continuous then you must set
#'   \code{exp_family = "gaussion"}. For \code{exp_family="binomial"} should be
#'   either a factor with two levels, or a two-column matrix of counts or
#'   proportions (the second column is treated as the target class; for a
#'   factor, the last level in alphabetical order is the target class)
#' @param y_test numeric vector of length \code{ntest} representing the
#'   responses for the test subjects. If continuous then you must set
#'   \code{exp_family = "gaussion"}. For \code{exp_family="binomial"} should be
#'   either a factor with two levels, or a two-column matrix of counts or
#'   proportions (the second column is treated as the target class; for a
#'   factor, the last level in alphabetical order is the target class).
#' @param s0 chracter vector of the active feature names, i.e., the features in
#'   \code{x_train} that are truly associated with the response.
#' @param summary the summary of each cluster. Can be the principal component or
#'   average. Default is \code{summary = "pc"} which takes the first
#'   \code{number_pc} principal components. Currently a maximum of 2 principal
#'   components can be chosen.
#' @param model Regression model to be fit on cluster summaries. Default is
#'   \code{model="lasso"} which corresponds to glmnet mixing parameter
#'   \code{alpha=1}. \code{model="elasticnet"} corresponds to glmnet mixing
#'   parameter \code{alpha=0.5}, \code{model="mcp"} and \code{model="scad"} are
#'   the non-convex models from the \code{\link[ncvreg]{}} package
#' @param exp_family Response type. See details for \code{y_train} argument
#'   above.
#' @param gene_groups data.frame that contains the group membership for each
#'   feature. The first column is called 'gene' and the second column should be
#'   called 'cluster'. The 'gene' column identifies the features and must be the
#'   same identifiers in the \code{x_train,x_test} matrices. The 'cluster'
#'   column is a numeric integer indicating the cluster group membership.  A
#'   cluster group membership of 0 implies the feature did not cluster into any
#'   group.
#' @param topgenes List of features to keep if \code{filter=TRUE}. Default is
#'   \code{topgenes = NULL} which means all features are kept for the analysis
#' @param stability Should stability measures be calculated. Default is
#'   \code{stability=FALSE}. See details
#' @param filter Should analysis be run on a subset of features. Default is
#'   \code{filter = FALSE}
#' @param include_E Should the environment variable be included in the
#'   regression analysis. Default is \code{include_E = TRUE}
#' @param include_interaction Should interaction effects between the features in
#'   \code{x_train} and the environment variable be fit. Default is
#'   \code{include_interaction=TRUE}
#' @param clust_type Method used to cluster the features. This is used for
#'   naming the output only and has no consequence for the results.
#'   \code{clust_type = "CLUST"} is the default which means that the environment
#'   varible was not used in the clustering step. \code{clust_type = "ECLUST"}
#'   means that the environment variable was used in the clustering aspect.
#' @param number_pc Number of principal components if \code{summary = "pc"}.
#'   Default is \code{number_pc = 1}. Can be either 1 or 2.
#' @details The stability of feature importance is defined as the variability of
#'   feature weights under perturbations of the training set, i.e., small
#'   modifications in the training set should not lead to considerable changes
#'   in the set of important covariates (Toloşi, L., & Lengauer, T. (2011)). A
#'   feature selection algorithm produces a weight, a ranking, and a subset of
#'   features. In the CLUST and ECLUST methods, we defined a predictor to be
#'   non-zero if its corresponding cluster representative weight was non-zero.
#'   Using 10-fold cross validation (CV), we evaluated the similarity between
#'   two features and their rankings using Pearson and Spearman correlation,
#'   respectively. For each CV fold we re-ran the models and took the average
#'   Pearson/Spearman correlation of the 10 choose 2 combinations of estimated
#'   coefficients vectors. To measure the similarity between two subsets of
#'   features we took the average of the Jaccard distance in each fold. A
#'   Jaccard distance of 1 indicates perfect agreement between two sets while no
#'   agreement will result in a distance of 0.
#' @references Toloşi, L., & Lengauer, T. (2011). \emph{Classification with
#'   correlated features: unreliability of feature ranking and solutions.
#'   Bioinformatics, 27(14), 1986-1994.}
#' @references Bhatnagar, SR., Yang, Y., Blanchette, M., Bouchard, L.,
#'   Khundrakpam, B., Evans, A., Greenwood, CMT. (2016+). \emph{An analytic
#'   approach for interpretable predictive models in high dimensional data, in
#'   the presence of interactions with exposures
#'   \href{http://sahirbhatnagar.com/slides/manuscript1_SB_v4.pdf}{Preprint}}
#' @references Langfelder, P., Zhang, B., & Horvath, S. (2008). \emph{Defining
#'   clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for
#'   R. Bioinformatics, 24(5), 719-720.}
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#'   \emph{Regularization Paths for Generalized Linear Models via Coordinate
#'   Descent, \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}}
#' @references Breheny, P. and Huang, J. (2011) \emph{Coordinate descent
#'   algorithms for nonconvex penalized regression, with applications to
#'   biological feature selection. Ann. Appl. Statist., 5: 232-253.}
#' @note \code{number_pc=2} will not work if there is only one feature in an
#'   estimated cluster
#' @importFrom pacman p_load
#' @import data.table
#' @import magrittr
#' @return This function has two different outputs depending on whether
#'   \code{stability = TRUE} or \code{stability = FALSE}
#'
#'   If \code{stability = TRUE} then this function returns a \code{p x 2}
#'   data.frame or data.table of regression coefficients without the intercept.
#'   The output of this is used for subsequent calculations of stability.
#'
#'   If \code{stability = FALSE} then returns a vector with the following
#'   elements (See Table 3: Measures of Performance in Bhatnagar et al (2016+)
#'   for definitions of each measure of performance): \item{mse or AUC}{Test set
#'   mean squared error if \code{exp_family = "gaussion"} or test set Area under
#'   the curve if \code{exp_family = "binomial"} calculated using the
#'   \code{\link[pROC]{roc}} function} \item{RMSE}{Square root of the mse. Only
#'   applicable if \code{exp_family = "gaussion"}} \item{Shat}{Number of
#'   non-zero estimated regression coefficients. The non-zero estimated
#'   regression coefficients are referred to as being selected by the model}
#'   \item{TPR}{true positive rate} \item{FPR}{false positive rate}
#'   \item{Correct Sparsity}{Correct true positives + correct true negative
#'   coefficients divided by the total number of features}
#'   \item{CorrectZeroMain}{Proportion of correct true negative main effects}
#'   \item{CorrectZeroInter}{Proportion of correct true negative interactions}
#'   \item{IncorrectZeroMain}{Proportion of incorrect true negative main
#'   effects} \item{IncorrectZeroInter}{Proportion of incorrect true negative
#'   interaction effects} \item{nclusters}{number of estimated clusters by the
#'   \code{\link[dynamicTreeCut]{cutreeDynamic}} function}
#' @export
#'
#' @examples
#'
s_pen_clust <- function(x_train,
                        x_test,
                        y_train,
                        y_test,
                        s0,
                        gene_groups,
                        summary = c("pc","avg"),
                        model = c("lasso", "elasticnet", "scad", "mcp"),
                        exp_family = c("gaussian","binomial"),
                        filter = F,
                        topgenes = NULL,
                        stability = F,
                        include_E = T,
                        include_interaction = T,
                        clust_type = c("CLUST","ECLUST"),
                        number_pc = 1) {

  # result[["clustersAddon"]] %>% print(nrows=Inf)
  # result[["clustersAddon"]][, table(cluster, module)]
  # result %>% names
  # stability = F; gene_groups = result[["clustersAddon"]];
  # x_train = result[["X_train"]] ; x_test = result[["X_test"]];
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # dim(x_train)
  # filter = F; filter_var = F; include_E = T; include_interaction = F;
  # s0 = result[["S0"]]; p = p ;true_beta = result[["beta_truth"]]
  # model = "lasso"; summary = "pc"; topgenes = NULL; clust_type="clust"; number_pc = 1
  # exp_family="binomial"

  clust_type <- match.arg(clust_type)
  summary <- match.arg(summary)
  model <- match.arg(model)
  exp_family <- match.arg(exp_family)

  message(sprintf("Summary measure: %s, Model: %s, Cluster Type: %s",
                  summary, model, clust_type))

  if (include_E == F & include_interaction == T) stop("include_E needs to be
                                                      TRUE if you want to include
                                                      interactions")

  if (is.null(topgenes) & filter == T) stop("Argument topgenes is missing but
                                            filter is TRUE. You need to provide
                                            a filtered list of genes if filter
                                            is TRUE")

  # train data which includes the relevant (filtered or not filtered genes
  # and E or not E)
  x_train_mod <- if (filter & !include_E) {
    x_train[, topgenes] %>% as.data.frame
  } else if (!filter & include_E) {
    x_train %>% as.data.frame
  } else if (!filter & !include_E) {
    x_train[,which(colnames(x_train) %ni% "E")] %>% as.data.frame
  } else if (filter & include_E) {
    x_train[, c(topgenes,"E")] %>% as.data.frame
  }

  # test data
  x_test_mod = if (filter & !include_E) {
    x_test[, topgenes] %>% as.data.frame
  } else if (!filter & include_E) {
    x_test %>% as.data.frame
  } else if (!filter & !include_E) {
    x_test[,which(colnames(x_test) %ni% "E")] %>% as.data.frame
  } else if (filter & include_E) {
    x_test[, c(topgenes,"E")] %>% as.data.frame
  }

  # these are only derived on the main effects genes.. E is only included in the model
  PC_and_avg <- extractPC(x_train = x_train_mod[,gene_groups$gene],
                          colors = gene_groups$cluster,
                          x_test = x_test_mod[,gene_groups$gene],
                          number_pc = number_pc)

  n.clusters <- PC_and_avg$nclusters

  # this contains either the averages or PCs for each module in a data.frame
  clust_data <- switch(summary,
                       avg = PC_and_avg$averageExpr,
                       pc = PC_and_avg$PC)

  ml.formula <- if (include_interaction & include_E) {
    paste0("y_train ~","(",paste0(colnames(clust_data), collapse = "+"),")*E") %>% as.formula
  } else if (!include_interaction & include_E) {
    paste0("y_train ~",paste0(colnames(clust_data), collapse = "+"),"+E") %>% as.formula
  } else if (!include_interaction & !include_E) {
    paste0("y_train ~",paste0(colnames(clust_data), collapse = "+")) %>% as.formula
  }

  # this is the same as ml.formula, except without the response.. this is used for
  # functions that have the x = and y = input instead of a formula input
  model.formula <- if (include_interaction & include_E) {
    paste0("~ 0+(",paste0(colnames(clust_data), collapse = "+"),")*E") %>% as.formula
  } else if (!include_interaction & include_E) {
    paste0("~0+",paste0(colnames(clust_data), collapse = "+"),"+E") %>% as.formula
  } else if (!include_interaction & !include_E) {
    paste0("~0+",paste0(colnames(clust_data), collapse = "+")) %>% as.formula
  }

  # this is the design matrix based on model.formula
  X.model.formula <- model.matrix(model.formula, data = if (include_E) {
    cbind(clust_data,x_train_mod[,"E", drop = F])
  } else clust_data %>% as.data.frame)

  df <- X.model.formula %>% cbind(y_train) %>% as.data.frame()

  clust_train_model <- switch(model,
                              lasso = {if (n.clusters != 1) {
                                pacman::p_load(char = "glmnet")
                                glmnet::cv.glmnet(x = X.model.formula, y = y_train, alpha = 1, family = exp_family)
                              } else NA },
                              elasticnet = {if (n.clusters != 1) {
                                pacman::p_load(char = "glmnet")
                                glmnet::cv.glmnet(x = X.model.formula, y = y_train, alpha = 0.5, family = exp_family)
                              } else NA },
                              scad = {
                                pacman::p_load(char = "ncvreg")
                                ncvreg::cv.ncvreg(X = X.model.formula, y = y_train,
                                          family = "gaussian", penalty = "SCAD")
                              },
                              mcp = {
                                pacman::p_load(char = "ncvreg")
                                ncvreg::cv.ncvreg(X = X.model.formula, y = y_train,
                                          family = "gaussian", penalty = "MCP")
                              }
                              # ,
                              # shim = {
                              #   cv.shim(x = X.model.formula, y = y_train,
                              #           main.effect.names = c(colnames(clust_data), if (include_E) "E"),
                              #           interaction.names = setdiff(colnames(X.model.formula),c(colnames(clust_data),"E")),
                              #           max.iter = 120, initialization.type = "ridge",
                              #           verbose = FALSE, parallel = TRUE, nfolds = 10)
                              # }
                              )
  # plot(clust_train_model)

  # here we give the coefficient stability on the clusters and not the individual genes
  coefs <- switch(model,
                  lasso = {
                    # need to return all 0's if there is only 1 cluster since lasso
                    # wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {
                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        data.table::as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  elasticnet = {
                    # need to return all 0's if there is only 1 cluster since lasso
                    # wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {
                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        data.table::as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  scad = {
                    # need to return all 0's if there is only 1 cluster since lasso
                    # wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {

                      coef(clust_train_model, lambda = clust_train_model$lambda.min) %>%
                        as.matrix %>%
                        data.table::as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  mcp = {
                    # need to return all 0's if there is only 1 cluster since lasso
                    # wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {

                      coef(clust_train_model, lambda = clust_train_model$lambda.min) %>%
                        as.matrix %>%
                        data.table::as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  }
                  # ,
                  # shim = {
                  #   dat <- data.table::data.table(Gene = colnames(X.model.formula),
                  #                                 coef.est = rep(0, ncol(X.model.formula)))
                  #   if (n.clusters != 1) {
                  #
                  #     coef(clust_train_model, s = "lambda.min") %>%
                  #       as.matrix %>%
                  #       as.data.table(keep.rownames = TRUE) %>%
                  #       magrittr::set_colnames(c("Gene","coef.est"))
                  #   } else dat
                  # }
                  )

  if (stability) {
    # remove intercept for stability measures
    return(coefs %>% magrittr::extract(-1, , drop = F))
  } else {

    non_zero_clusters <- coefs[-1, , ][coef.est != 0] %>%
      magrittr::use_series("Gene")

    # need to determine which of non_zero_cluters are main effects and which
    # are interactions
    non_zero_clusters_interactions <- grep(":",non_zero_clusters, value = T)

    # this checks if the environment is non-zero
    non_zero_environment <- grep("^E", non_zero_clusters, value = T,
                                 ignore.case = TRUE)
    non_zero_clusters_main_effects <- setdiff(non_zero_clusters,
                                              c(non_zero_clusters_interactions,
                                                non_zero_environment))

    # this includes the environment if the environment is non-zero
    n.non_zero_clusters <- coefs[-1, , ][coef.est != 0] %>%
      magrittr::use_series("Gene") %>%
      length

    # need to get the genes corresponding to the non-zero clusters
    # NOTE: this also includes non-zero cluster:Environment interactions

    # genes corresponding to non-zero main effect clusters
    # this list might not be unique if clust_type="Addon" because the same gene
    # can be in different clusters
    clust.S.hat.main <- gene_groups[cluster %in%
                                      as.numeric(
                                        unlist(
                                          stringr::str_extract_all(
                                            non_zero_clusters_main_effects, "(\\d+)$")
                                        )
                                      ),gene]

    # identical(gene_groups[cluster %in% c(3,12),gene],
    #           clust.S.hat.main)
    # identical(unique(clust.S.hat.main), clust.S.hat.main)
    # table(clust.S.hat.main)

    # this is the same as gene_groups, but the gene names contain E
    # so that we can extract the interactions corresponding to the chose clusters
    gene_groups_E <- copy(gene_groups)
    gene_groups_E[,gene:=paste0(gene,":E")]

    clust.S.hat.interaction <- gene_groups_E[cluster %in%
                                               as.numeric(
                                                 unlist(
                                                   stringr::str_extract_all(
                                                     stringr::str_extract_all(non_zero_clusters_interactions,"^.*?(?=:)"),
                                                     "(\\d+)$")
                                                 )
                                               ),gene]

    # this represents all the genes corresponding to the non-zero PC or avg
    # this list might not be unique if clust_type="Addon"
    # identical(unique(clust.S.hat), clust.S.hat)
    # I will double count if a model takes a gene more than once. ie.
    # if the same gene gets selected twice, then this will contribute 2 to the
    # number of non-zero estimated coefficients
    clust.S.hat <- c(clust.S.hat.main, non_zero_environment,
                     clust.S.hat.interaction)


    clust_data_test <- switch(summary,
                              avg = PC_and_avg$averageExprTest,
                              pc = PC_and_avg$PCTest)

    # need intercept for prediction
    model.formula_test <- if (include_interaction & include_E) {
      paste0("~ 1+(",paste0(colnames(clust_data_test), collapse = "+"),")*E") %>% as.formula
    } else if (!include_interaction & include_E) {
      paste0("~1+",paste0(colnames(clust_data_test), collapse = "+"),"+E") %>% as.formula
    } else if (!include_interaction & !include_E) {
      paste0("~1+",paste0(colnames(clust_data_test), collapse = "+")) %>% as.formula
    }


    # this includes the intercept!
    X.model.formula_test <- model.matrix(model.formula_test,
                                         data = if (include_E) {
                                           cbind(clust_data_test,x_test_mod[,"E", drop = F])
                                         } else clust_data_test %>% as.data.frame)

    # True Positive Rate
    clust.TPR <- length(intersect(clust.S.hat, s0))/length(s0)

    # True negatives
    trueNegs <- setdiff(colnames(x_train_mod), s0)
    # identical(setdiff(colnames(x_train_mod), s0), setdiff(colnames(x_train), s0))

    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(x_train_mod),clust.S.hat)

    # how many of the terms identified by the model as zero, were actually zero
    # use to calculate correct sparsity as defined by Witten et al in the
    # Cluster Elastic Net paper Technometrics 2013
    C1 <- sum(modelIdentifyZero %in% trueNegs)
    C2 <- length(intersect(clust.S.hat, s0))
    clust.correct_sparsity <- (C1 + C2)/(ncol(x_train_mod))

    # this is from Interaction Screening for Ultrahigh Dimensional Data by ning hao and hao helen zhang
    true.interaction_names <- grep(":", s0, value = T)
    true.main_effect_names <- setdiff(s0, true.interaction_names)

    all.interaction_names <- grep(":", colnames(x_train_mod), value = T)
    all.main_effect_names <- setdiff(colnames(x_train_mod), all.interaction_names)

    true.negative_main_effects <- setdiff(all.main_effect_names, true.main_effect_names)
    true.negative_interaction_effects <- setdiff(all.interaction_names, true.interaction_names)

    (clust.correct_zeros_main_effects <- sum(setdiff(all.main_effect_names, c(clust.S.hat.main, non_zero_environment)) %in% true.negative_main_effects)/ length(true.negative_main_effects))
    (clust.correct_zeros_interaction_effects <- sum(setdiff(all.interaction_names, clust.S.hat.interaction) %in% true.negative_interaction_effects)/ length(true.negative_interaction_effects))

    (clust.incorrect_zeros_main_effects <- sum(setdiff(all.main_effect_names, c(clust.S.hat.main, non_zero_environment)) %in% true.main_effect_names)/ length(true.main_effect_names))
    (clust.incorrect_zeros_interaction_effects <- sum(setdiff(all.interaction_names, clust.S.hat.interaction) %in% true.interaction_names)/ length(true.interaction_names))

    # False Positive Rate = FP/(FP + TN) = FP / True number of 0 coefficients
    (clust.FPR <- sum(clust.S.hat %ni% s0)/(sum(clust.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs)))

    # Mean Squared Error
    (clust.mse <- crossprod(X.model.formula_test %*% coefs$coef.est - y_test)/length(y_test))


    if (exp_family == "binomial") {
      pred_response <- predict(clust_train_model, newx = X.model.formula_test[,-which(colnames(X.model.formula_test)=="(Intercept)")],
                               type = "response", s = "lambda.min")

      clust.AUC <- pROC::roc(y_test,as.numeric(pred_response))$auc %>% as.numeric()

    }

    # Root Mean Squared Error
    (clust.RMSE <- sqrt(crossprod(X.model.formula_test %*% coefs$coef.est - y_test)/length(y_test)))

    # remove intercept for prediction error formula given by ||X\beta - X\hat{\beta}||_2
    # given in Witten 2013 Cluster ENET paper in Technometrics
    # (clust.test_set_pred_error <- sqrt(crossprod(as.matrix(x_test_mod) %*% as.numeric(true_beta) - X.model.formula_test[,-1] %*% coefs$coef.est[-1])))

    # mse.null
    (mse_null <- crossprod(mean(y_test) - y_test)/length(y_test))

    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    # clust.r2 <- (mse_null - clust.mse)/mse_null

    # clust.adj.r2 <- 1 - (1 - clust.r2)*(nrow(x_test) - 1)/(nrow(x_test) - n.non_zero_clusters - 1)


    ls <- if (exp_family=="binomial") {
      list(clust.mse = clust.mse,
           clust.RMSE,
           clust.AUC = clust.AUC,
           clust.S.hat = length(clust.S.hat),
           clust.TPR = clust.TPR,
           clust.FPR = clust.FPR,
           clust.correct_sparsity = clust.correct_sparsity,
           clust.correct_zeros_main_effects,
           clust.correct_zeros_interaction_effects,
           clust.incorrect_zeros_main_effects,
           clust.incorrect_zeros_interaction_effects,
           n.clusters) } else if (exp_family=="gaussian") {
             list(clust.mse = clust.mse,
                  clust.RMSE,
                  clust.S.hat = length(clust.S.hat),
                  clust.TPR = clust.TPR,
                  clust.FPR = clust.FPR,
                  clust.correct_sparsity = clust.correct_sparsity,
                  clust.correct_zeros_main_effects,
                  clust.correct_zeros_interaction_effects,
                  clust.incorrect_zeros_main_effects,
                  clust.incorrect_zeros_interaction_effects,
                  n.clusters)
           }

    if (exp_family=="binomial") {

      names(ls) <- c(paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_AUC"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_nclusters")
      )

    } else if (exp_family=="gaussian") {
      names(ls) <- c(paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"),
                     paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_nclusters")
      )
    }

    return(ls)

  }
}




#' Fit Penalized Regression Models on Simulated Data
#'
#' @description This function can run penalized regression models on the
#'   untransformed design matrix. To be used with simulated data where the
#'   'truth' is known i.e., you know which features are associated with the
#'   response. This function was used to produce the simulation results in
#'   Bhatnagar et al. 2016. Can run lasso, elasticnet, SCAD or MCP models
#'
#' @inheritParams s_pen_clust
#' @return This function has two different outputs depending on whether
#'   \code{stability = TRUE} or \code{stability = FALSE}
#'
#'   If \code{stability = TRUE} then this function returns a \code{p x 2}
#'   data.frame or data.table of regression coefficients without the intercept.
#'   The output of this is used for subsequent calculations of stability.
#'
#'   If \code{stability = FALSE} then returns a vector with the following
#'   elements (See Table 3: Measures of Performance in Bhatnagar et al (2016+)
#'   for definitions of each measure of performance): \item{mse or AUC}{Test set
#'   mean squared error if \code{exp_family = "gaussion"} or test set Area under
#'   the curve if \code{exp_family = "binomial"} calculated using the
#'   \code{\link[pROC]{roc}} function} \item{RMSE}{Square root of the mse. Only
#'   applicable if \code{exp_family = "gaussion"}} \item{Shat}{Number of
#'   non-zero estimated regression coefficients. The non-zero estimated
#'   regression coefficients are referred to as being selected by the model}
#'   \item{TPR}{true positive rate} \item{FPR}{false positive rate}
#'   \item{Correct Sparsity}{Correct true positives + correct true negative
#'   coefficients divided by the total number of features}
#'   \item{CorrectZeroMain}{Proportion of correct true negative main effects}
#'   \item{CorrectZeroInter}{Proportion of correct true negative interactions}
#'   \item{IncorrectZeroMain}{Proportion of incorrect true negative main
#'   effects} \item{IncorrectZeroInter}{Proportion of incorrect true negative
#'   interaction effects}
#' @details The stability of feature importance is defined as the variability of
#'   feature weights under perturbations of the training set, i.e., small
#'   modifications in the training set should not lead to considerable changes
#'   in the set of important covariates (Toloşi, L., & Lengauer, T. (2011)). A
#'   feature selection algorithm produces a weight, a ranking, and a subset of
#'   features. In the CLUST and ECLUST methods, we defined a predictor to be
#'   non-zero if its corresponding cluster representative weight was non-zero.
#'   Using 10-fold cross validation (CV), we evaluated the similarity between
#'   two features and their rankings using Pearson and Spearman correlation,
#'   respectively. For each CV fold we re-ran the models and took the average
#'   Pearson/Spearman correlation of the 10 choose 2 combinations of estimated
#'   coefficients vectors. To measure the similarity between two subsets of
#'   features we took the average of the Jaccard distance in each fold. A
#'   Jaccard distance of 1 indicates perfect agreement between two sets while no
#'   agreement will result in a distance of 0.
#' @references Toloşi, L., & Lengauer, T. (2011). \emph{Classification with
#'   correlated features: unreliability of feature ranking and solutions.
#'   Bioinformatics, 27(14), 1986-1994.}
#' @references Bhatnagar, SR., Yang, Y., Blanchette, M., Bouchard, L.,
#'   Khundrakpam, B., Evans, A., Greenwood, CMT. (2016+). \emph{An analytic
#'   approach for interpretable predictive models in high dimensional data, in
#'   the presence of interactions with exposures
#'   \href{http://sahirbhatnagar.com/slides/manuscript1_SB_v4.pdf}{Preprint}}
#' @references Langfelder, P., Zhang, B., & Horvath, S. (2008). \emph{Defining
#'   clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for
#'   R. Bioinformatics, 24(5), 719-720.}
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#'   \emph{Regularization Paths for Generalized Linear Models via Coordinate
#'   Descent, \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}}
#' @references Breheny, P. and Huang, J. (2011) \emph{Coordinate descent
#'   algorithms for nonconvex penalized regression, with applications to
#'   biological feature selection. Ann. Appl. Statist., 5: 232-253.}
#' @importFrom pacman p_load
#' @import data.table
#' @import magrittr
#' @export
#'
#' @examples
#' \dontrun{hello}
s_pen_separate <- function(x_train,
                      x_test,
                      y_train,
                      y_test,
                      s0,
                      expFamily = c("gaussian","binomial"),
                      model = c("lasso", "elasticnet", "scad", "mcp"),
                      topgenes = NULL,
                      stability = F,
                      filter = F,
                      include_E = T,
                      include_interaction = T){

  # stability = F; x_train = result[["X_train"]] ; x_test = result[["X_test"]] ;
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = F;
  # s0 = result[["S0"]]; p = p ;
  # model = "lasso"; topgenes = NULL; true_beta = result[["beta_truth"]]

  # stability = F; x_train = result_interaction[["X_train"]] ; x_test = result_interaction[["X_test"]] ;
  # y_train = result_interaction[["Y_train"]] ; y_test = result_interaction[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result_interaction[["S0"]]; p = 1000 ;
  # model = "scad"; topgenes = NULL; true_beta = result_interaction[["beta_truth"]]

  # result[["clustersAddon"]] %>% print(nrows=Inf)
  # result[["clustersAddon"]][, table(cluster, module)]
  # result %>% names
  # stability = F; gene_groups = result[["clustersAll"]];
  # x_train = result[["X_train"]] ; x_test = result[["X_test"]];
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result[["S0"]]; p = p ;true_beta = result[["beta_truth"]]
  # model = "lasso"; summary = "pc"; topgenes = NULL; clust_type="clust"; nPC = 1

  # model: "scad", "mcp", "lasso", "elasticnet", "ridge"
  # filter: T or F based on univariate filter
  expFamily <- match.arg(expFamily)

  print(paste(model,"filter = ",
              filter, "filter_var = ",
              include_E, "include_interaction = ",
              include_interaction, sep = " "))

  if (include_E == F & include_interaction == T) stop("include_E needs to be
                                                      TRUE if you want to include
                                                      interactions")
  #   if (filter == F & include_interaction == T) stop("Interaction can only be run
  #                                                      if filter is TRUE.
  #                                                      This is to avoid exceedingly
  #                                                      large models")
  if (is.null(topgenes) & filter == T) stop("Argument topgenes is missing but
                                            filter is TRUE. You need to provide
                                            a filtered list of genes if filter
                                            is TRUE")

  #gene.names <- colnames(x_train)[which(colnames(x_train) %ni% "E")]

  # penalization model
  pen_model <- switch(model,
                      {
                        pacman::p_load(char = "glmnet")
                        lasso = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 1, family = expFamily)
                        },
                      elasticnet = {
                        pacman::p_load(char = "glmnet")
                        glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 0.5, family = expFamily)
                        },
                      ridge = {
                        pacman::p_load(char = "glmnet")
                        glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 0, family = expFamily)
                        },
                      scad = {
                        pacman::p_load(char = "ncvreg")
                        ncvreg::cv.ncvreg(X = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train,
                        family = "gaussian", penalty = "SCAD")},
                      mcp = {
                        pacman::p_load(char = "ncvreg")
                        ncvreg::cv.ncvreg(X = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train,
                        family = "gaussian", penalty = "MCP")}
  )

  # here we give the coefficient stability on the individual genes
  coefs <- coef(pen_model, s = "lambda.min") %>%
    as.matrix %>%
    data.table::as.data.table(keep.rownames = TRUE) %>%
    magrittr::set_colnames(c("Gene","coef.est")) %>%
    magrittr::extract(-1,)


  if (stability) {
    # remove intercept for stability measures
    return(coefs)
  } else {

    pen.S.hat <- coefs[coef.est != 0] %>% magrittr::use_series("Gene")
    pen.S.hat.interaction <- grep(":", pen.S.hat, value = T)
    pen.S.hat.main <- setdiff(pen.S.hat, pen.S.hat.interaction)

    pen.pred <- if (model %in% c("lasso","elasticnet","ridge")) {
      predict(pen_model, newx =  if (!include_E) as.matrix(x_test[,-grep("E", colnames(x_test))]) else
        as.matrix(x_test), s = "lambda.min") } else if (model %in% c("scad","mcp")) {
          predict(pen_model, X =  if (!include_E) as.matrix(x_test[,-grep("E", colnames(x_test))]) else
            as.matrix(x_test),
            lambda = pen_model$lambda.min)
        }

    # True Positive Rate
    pen.TPR <- length(intersect(pen.S.hat, s0))/length(s0)

    # True negatives
    trueNegs <- setdiff(colnames(x_train), s0)

    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(x_train),pen.S.hat)

    # how many of the terms identified by the model as zero, were actually zero
    # use to calculate correct sparsity as defined by Witten et al in the
    # Cluster Elastic Net paper Technometrics 2013
    C1 <- sum(modelIdentifyZero %in% trueNegs)
    C2 <- length(intersect(pen.S.hat, s0))
    correct_sparsity <- (C1 + C2)/(ncol(x_train))

    # this is from Interaction Screening for Ultrahigh Dimensional Data by ning hao and hao helen zhang
    true.interaction_names <- grep(":", s0, value = T)
    true.main_effect_names <- setdiff(s0, true.interaction_names)

    all.interaction_names <- grep(":", colnames(x_train), value = T)
    all.main_effect_names <- setdiff(colnames(x_train), all.interaction_names)

    true.negative_main_effects <- setdiff(all.main_effect_names, true.main_effect_names)
    true.negative_interaction_effects <- setdiff(all.interaction_names, true.interaction_names)

    (pen.correct_zeros_main_effects <- sum(setdiff(all.main_effect_names, pen.S.hat.main) %in% true.negative_main_effects)/ length(true.negative_main_effects))
    (pen.correct_zeros_interaction_effects <- sum(setdiff(all.interaction_names, pen.S.hat.interaction) %in% true.negative_interaction_effects)/ length(true.negative_interaction_effects))

    (pen.incorrect_zeros_main_effects <- sum(setdiff(all.main_effect_names, pen.S.hat) %in% true.main_effect_names)/ length(true.main_effect_names))
    (pen.incorrect_zeros_interaction_effects <- sum(setdiff(all.interaction_names, pen.S.hat.interaction) %in% true.interaction_names)/ length(true.interaction_names))

    # False Positive Rate = FP/(FP + TN) = FP / True number of 0 coefficients
    (pen.FPR <- sum(pen.S.hat %ni% s0)/(sum(pen.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs)))


    # Mean Squared Error
    (pen.mse <- crossprod(pen.pred - y_test)/length(y_test))

    # Root Mean Squared Error
    (pen.RMSE <- sqrt(crossprod(pen.pred - y_test)/length(y_test)))

    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)

    if (expFamily == "binomial") {
      pred_response <- predict(pen_model, newx =  if (!include_E) as.matrix(x_test[,-grep("E", colnames(x_test))]) else
        as.matrix(x_test), s = "lambda.min", type = "response")

      pen.AUC <- pROC::roc(y_test,as.numeric(pred_response))$auc %>% as.numeric()

    }

    ls <- if (expFamily == "binomial") {
      list(pen.mse = as.numeric(pen.mse),
           pen.RMSE = as.numeric(pen.RMSE),
           pen.AUC = pen.AUC,
           pen.S.hat = length(pen.S.hat),
           pen.TPR = pen.TPR,
           pen.FPR = pen.FPR,
           correct_sparsity,
           pen.correct_zeros_main_effects,
           pen.correct_zeros_interaction_effects,
           pen.incorrect_zeros_main_effects,
           pen.incorrect_zeros_interaction_effects
      ) } else if (expFamily=="gaussian") {

        list(pen.mse = as.numeric(pen.mse),
             pen.RMSE = as.numeric(pen.RMSE),
             pen.S.hat = length(pen.S.hat),
             pen.TPR = pen.TPR,
             pen.FPR = pen.FPR,
             correct_sparsity,
             pen.correct_zeros_main_effects,
             pen.correct_zeros_interaction_effects,
             pen.incorrect_zeros_main_effects,
             pen.incorrect_zeros_interaction_effects
        )
      }

    names(ls) <- if (expFamily == "binomial") {
      c(paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_AUC"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
    } else if (expFamily=="gaussian") {
      c(paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
        paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
    }

    return(ls)
  }

}






#' Fit Multivariate Adaptive Regression Splines on Simulated Data
#'
#' @description This function can run Friedman's MARS models on the
#'   untransformed design matrix. To be used with simulated data where the
#'   'truth' is known i.e., you know which features are associated with the
#'   response. This function was used to produce the simulation results in
#'   Bhatnagar et al. 2016. Uses caret functions to tune the degree and the
#'   nprune parameters
#'
#' @param model Type of non-linear model to be fit. Currently only Friedman's
#'   MARS is supported.
#' @param ... other parameters passed to \code{\link[caret]{trainControl}}
#'   function
#' @inheritParams s_pen_clust
#' @return This function has two different outputs depending on whether
#'   \code{stability = TRUE} or \code{stability = FALSE}
#'
#'   If \code{stability = TRUE} then this function returns a \code{p x 2}
#'   data.frame or data.table of regression coefficients without the intercept.
#'   The output of this is used for subsequent calculations of stability.
#'
#'   If \code{stability = FALSE} then returns a vector with the following
#'   elements (See Table 3: Measures of Performance in Bhatnagar et al (2016+)
#'   for definitions of each measure of performance): \item{mse or AUC}{Test set
#'   mean squared error if \code{exp_family = "gaussion"} or test set Area under
#'   the curve if \code{exp_family = "binomial"} calculated using the
#'   \code{\link[pROC]{roc}} function} \item{RMSE}{Square root of the mse. Only
#'   applicable if \code{exp_family = "gaussion"}} \item{Shat}{Number of
#'   non-zero estimated regression coefficients. The non-zero estimated
#'   regression coefficients are referred to as being selected by the model}
#'   \item{TPR}{true positive rate} \item{FPR}{false positive rate}
#'   \item{Correct Sparsity}{Correct true positives + correct true negative
#'   coefficients divided by the total number of features}
#'   \item{CorrectZeroMain}{Proportion of correct true negative main effects}
#'   \item{CorrectZeroInter}{Proportion of correct true negative interactions}
#'   \item{IncorrectZeroMain}{Proportion of incorrect true negative main
#'   effects} \item{IncorrectZeroInter}{Proportion of incorrect true negative
#'   interaction effects}
#' @details This function first does 10 fold cross-validation to tune the degree
#'   (either 1 or 2) using the \code{\link[caret]{train}} function with
#'   \code{method="earth"} and nprune is fixed at 1000. Then the
#'   \code{\link[earth]{earth}} function is used, with \code{nk = 1000} and
#'   \code{pmethod = "backward"} to fit the MARS model using the optimal degree
#'   from the 10-fold CV.
#' @export
#'
#' @examples
#' \dontrun{hello}
s_mars_separate <- function(x_train,
                       x_test,
                       y_train,
                       y_test,
                       s0,
                       model = c("MARS"),
                       expFamily = c("gaussian", "binomial"),
                       topgenes = NULL,
                       stability = F,
                       filter = F,
                       include_E = T,
                       include_interaction = T, ...){

  # stability = F; x_train = result[["X_train"]] ; x_test = result[["X_test"]] ;
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result[["S0"]]; p = p ;
  # model = "MARS"; topgenes = NULL; true_beta = result[["beta_truth"]]

  # stability = F; x_train = result_interaction[["X_train"]] ; x_test = result_interaction[["X_test"]] ;
  # y_train = result_interaction[["Y_train"]] ; y_test = result_interaction[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result_interaction[["S0"]]; p = 1000 ;
  # model = "scad"; topgenes = NULL; true_beta = result_interaction[["beta_truth"]]

  # result[["clustersAddon"]] %>% print(nrows=Inf)
  # result[["clustersAddon"]][, table(cluster, module)]
  # result %>% names
  # stability = F; gene_groups = result[["clustersAll"]];
  # x_train = result[["X_train"]] ; x_test = result[["X_test"]];
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result[["S0"]]; p = p ;true_beta = result[["beta_truth"]]
  # model = "lasso"; summary = "pc"; topgenes = NULL; clust_type="clust"; nPC = 1

  # model: "scad", "mcp", "lasso", "elasticnet", "ridge"
  # filter: T or F based on univariate filter

  expFamily <- match.arg(expFamily)

  print(paste(model,"filter = ", filter, "include_E = ",
              include_E, "include_interaction = ", include_interaction, sep = " "))

  if (include_E == F & include_interaction == T) stop("include_E needs to be
                                                      TRUE if you want to include
                                                      interactions")
  #   if (filter == F & include_interaction == T) stop("Interaction can only be run
  #                                                      if filter is TRUE.
  #                                                      This is to avoid exceedingly
  #                                                      large models")
  if (is.null(topgenes) & filter == T) stop("Argument topgenes is missing but
                                            filter is TRUE. You need to provide
                                            a filtered list of genes if filter
                                            is TRUE")

  #gene.names <- colnames(x_train)[which(colnames(x_train) %ni% "E")]

  # mars model
  mars_model <- switch(model,
                       MARS = {
                         pacman::p_load(char = "caret")
                         pacman::p_load(char = "earth")

                         fitControl <-  caret::trainControl(method = "cv",
                                                            verboseIter = FALSE, ...)

                         marsGrid <- expand.grid(.degree = 1:2, .nprune = 1000)

                         switch(expFamily,
                                gaussian = {
                                  mars_tuned <- caret::train(as.matrix(x_train),
                                                             y_train,
                                                             method = "earth",
                                                             trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                                                             tuneGrid = marsGrid,
                                                             trControl = fitControl)

                                  earth::earth(x = as.matrix(x_train),
                                               y = y_train,
                                               keepxy = TRUE,
                                               pmethod = "backward",
                                               nk = 1000,
                                               degree = mars_tuned$bestTune$degree,
                                               trace = 1, nfold = 10) },
                                binomial = {

                                  mars_tuned <- caret::train(as.matrix(x_train),
                                                             as.factor(y_train),
                                                             method = "earth",
                                                             trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                                                             glm=list(family=binomial),
                                                             tuneGrid = marsGrid,
                                                             trControl = fitControl)

                                  earth::earth(x = as.matrix(x_train),
                                               y = as.factor(y_train),
                                               keepxy = TRUE,
                                               pmethod = "backward",
                                               nk = 1000,
                                               glm=list(family=binomial),
                                               degree = mars_tuned$bestTune$degree,
                                               trace = 1, nfold = 10)
                                })
                       }
  )


  # selected genes
  # coef(mars_model)
  # get.used.pred.names(mars_model)
  #
  # plot(mars_model, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)
  # plot.earth.models(mars_model$cv.list, which=1)
  # plot(mars_model)
  # plot(mars_model, which=1,
  #      col.mean.infold.rsq="blue", col.infold.rsq="lightblue",
  #      col.grsq=0, col.rsq=0, col.vline=0, col.oof.vline=0)


  # ONLY Jaccard index can be calculated for MARS
  # since get.used.pred.names returns only the non-zero coefficients,
  # we give a coefficient of 1 here so that the stability calculation works
  # because it takes non-zero coef estimates

  coefs <- data.frame(get.used.pred.names(mars_model), rep(1, length(get.used.pred.names(mars_model))), stringsAsFactors = FALSE) %>%
    data.table::as.data.table(keep.rownames = FALSE) %>%
    magrittr::set_colnames(c("Gene","coef.est"))

  if (stability) {
    # remove intercept for stability measures
    return(coefs)
  } else {

    mars.S.hat <- get.used.pred.names(mars_model)
    mars.S.hat.interaction <- grep(":", mars.S.hat, value = T)
    mars.S.hat.main <- setdiff(mars.S.hat, mars.S.hat.interaction)

    mars.pred <- predict(mars_model, newdata = x_test, trace = 4)

    if (expFamily == "binomial") {
      pred_response <- predict(mars_model, newdata = x_test, trace = 4,
                               type = "response")

      mars.AUC <- pROC::roc(y_test,as.numeric(pred_response))$auc %>% as.numeric()

    }

    # True Positive Rate
    mars.TPR <- length(intersect(mars.S.hat, s0))/length(s0)

    # True negatives
    trueNegs <- setdiff(colnames(x_train), s0)

    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(x_train),mars.S.hat)

    # how many of the terms identified by the model as zero, were actually zero
    # use to calculate correct sparsity as defined by Witten et al in the
    # Cluster Elastic Net paper Technometrics 2013
    C1 <- sum(modelIdentifyZero %in% trueNegs)
    C2 <- length(intersect(mars.S.hat, s0))
    correct_sparsity <- (C1 + C2)/(ncol(x_train))

    # this is from Interaction Screening for Ultrahigh Dimensional Data by ning hao and hao helen zhang
    true.interaction_names <- grep(":", s0, value = T)
    true.main_effect_names <- setdiff(s0, true.interaction_names)

    all.interaction_names <- grep(":", colnames(x_train), value = T)
    all.main_effect_names <- setdiff(colnames(x_train), all.interaction_names)

    true.negative_main_effects <- setdiff(all.main_effect_names, true.main_effect_names)
    true.negative_interaction_effects <- setdiff(all.interaction_names, true.interaction_names)

    (mars.correct_zeros_main_effects <- sum(setdiff(all.main_effect_names, mars.S.hat.main) %in% true.negative_main_effects)/ length(true.negative_main_effects))
    (mars.correct_zeros_interaction_effects <- sum(setdiff(all.interaction_names, mars.S.hat.interaction) %in% true.negative_interaction_effects)/ length(true.negative_interaction_effects))

    (mars.incorrect_zeros_main_effects <- sum(setdiff(all.main_effect_names, mars.S.hat) %in% true.main_effect_names)/ length(true.main_effect_names))
    (mars.incorrect_zeros_interaction_effects <- sum(setdiff(all.interaction_names, mars.S.hat.interaction) %in% true.interaction_names)/ length(true.interaction_names))

    # False Positive Rate = FP/(FP + TN) = FP / True number of 0 coefficients
    (mars.FPR <- sum(mars.S.hat %ni% s0)/(sum(mars.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs)))

    # # False Positive Rate
    # mars.FPR <- sum(mars.S.hat %ni% s0)/(p - length(s0))

    # Mean Squared Error
    (mars.mse <- crossprod(mars.pred - y_test)/length(y_test))

    # Root Mean Squared Error
    (mars.RMSE <- sqrt(crossprod(mars.pred - y_test)/length(y_test)))

    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)

    ls <- switch(expFamily,
                 gaussian = {list(mars.mse = as.numeric(mars.mse),
                                  mars.RMSE = as.numeric(mars.RMSE),
                                  mars.S.hat = length(mars.S.hat),
                                  mars.TPR = mars.TPR,
                                  mars.FPR = mars.FPR,
                                  correct_sparsity,
                                  mars.correct_zeros_main_effects,
                                  mars.correct_zeros_interaction_effects,
                                  mars.incorrect_zeros_main_effects,
                                  mars.incorrect_zeros_interaction_effects)},
                 binomial = {list(mars.mse = as.numeric(mars.mse),
                                  mars.RMSE = as.numeric(mars.RMSE),
                                  mars.AUC = mars.AUC,
                                  mars.S.hat = length(mars.S.hat),
                                  mars.TPR = mars.TPR,
                                  mars.FPR = mars.FPR,
                                  correct_sparsity,
                                  mars.correct_zeros_main_effects,
                                  mars.correct_zeros_interaction_effects,
                                  mars.incorrect_zeros_main_effects,
                                  mars.incorrect_zeros_interaction_effects)

                 })

    names(ls) <- switch(expFamily,
                        gaussian = {
                          c(paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))},
                        binomial = {
                          c(paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_AUC"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
                        })
    return(ls)
  }
}




#' Fit MARS Models on Simulated Cluster Summaries
#'
#' @description This function creates summaries of the given clusters (e.g. 1st
#'   PC or average), and then runs Friedman's MARS on those
#'   summaries. To be used with simulated data where the 'truth' is known i.e.,
#'   you know which features are associated with the response. This function was
#'   used to produce the simulation results in Bhatnagar et al. 2016.
#'
#' @param model Type of non-linear model to be fit. Currently only Friedman's
#'   MARS is supported.
#' @inheritParams s_pen_clust
#' @details This function first does 10 fold cross-validation to tune the degree
#'   (either 1 or 2) using the \code{\link[caret]{train}} function with
#'   \code{method="earth"} and nprune is fixed at 1000. Then the
#'   \code{\link[earth]{earth}} function is used, with \code{nk = 1000} and
#'   \code{pmethod = "backward"} to fit the MARS model using the optimal degree
#'   from the 10-fold CV.
#' @return This function has two different outputs depending on whether
#'   \code{stability = TRUE} or \code{stability = FALSE}
#'
#'   If \code{stability = TRUE} then this function returns a \code{p x 2}
#'   data.frame or data.table of regression coefficients without the intercept.
#'   The output of this is used for subsequent calculations of stability.
#'
#'   If \code{stability = FALSE} then returns a vector with the following
#'   elements (See Table 3: Measures of Performance in Bhatnagar et al (2016+)
#'   for definitions of each measure of performance): \item{mse or AUC}{Test set
#'   mean squared error if \code{exp_family = "gaussion"} or test set Area under
#'   the curve if \code{exp_family = "binomial"} calculated using the
#'   \code{\link[pROC]{roc}} function} \item{RMSE}{Square root of the mse. Only
#'   applicable if \code{exp_family = "gaussion"}} \item{Shat}{Number of
#'   non-zero estimated regression coefficients. The non-zero estimated
#'   regression coefficients are referred to as being selected by the model}
#'   \item{TPR}{true positive rate} \item{FPR}{false positive rate}
#'   \item{Correct Sparsity}{Correct true positives + correct true negative
#'   coefficients divided by the total number of features}
#'   \item{CorrectZeroMain}{Proportion of correct true negative main effects}
#'   \item{CorrectZeroInter}{Proportion of correct true negative interactions}
#'   \item{IncorrectZeroMain}{Proportion of incorrect true negative main
#'   effects} \item{IncorrectZeroInter}{Proportion of incorrect true negative
#'   interaction effects}
#' @export
#'
#' @examples
#' \dontrun{hello}
s_mars_clust <- function(x_train,
                         x_test,
                         y_train,
                         y_test,
                         s0,
                         summary = c("pc","avg"),
                         model = c("MARS"),
                         expFamily = c("gaussian","binomial"),
                         gene_groups,
                         true_beta = NULL,
                         topgenes = NULL,
                         stability = F,
                         filter = F,
                         include_E = F,
                         include_interaction = F,
                         clust_type = c("clust","Eclust","Addon"),
                         nPC = 1) {

  # result[["clustersAddon"]] %>% print(nrows=Inf)
  # result[["clustersAddon"]][, table(cluster, module)]
  # result %>% names
  # stability = F; gene_groups = result[["clustersAddon"]];
  # x_train = result[["X_train"]] ; x_test = result[["X_test"]];
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # dim(x_train)
  # filter = F; filter_var = F; include_E = T; include_interaction = F;
  # s0 = result[["S0"]]; p = p ;true_beta = result[["beta_truth"]]
  # model = "MARS"; summary = "pc"; topgenes = NULL; clust_type="Eclust"; nPC = 1
  # expFamily = "binomial"

  clust_type <- match.arg(clust_type)
  summary <- match.arg(summary)
  model <- match.arg(model)
  expFamily <- match.arg(expFamily)

  message(sprintf("Summary measure: %s, Model: %s, Cluster Type: %s",
                  summary, model, clust_type))

  if (include_E == F & include_interaction == T) stop("include_E needs to be
                                                      TRUE if you want to include
                                                      interactions")

  if (is.null(topgenes) & filter == T) stop("Argument topgenes is missing but
                                            filter is TRUE. You need to provide
                                            a filtered list of genes if filter
                                            is TRUE")

  # train data which includes the relevant (filtered or not filtered genes
  # and E or not E)
  x_train_mod <- if (filter & !include_E) {
    x_train[, topgenes] %>% as.data.frame
  } else if (!filter & include_E) {
    x_train %>% as.data.frame
  } else if (!filter & !include_E) {
    x_train[,which(colnames(x_train) %ni% "E")] %>% as.data.frame
  } else if (filter & include_E) {
    x_train[, c(topgenes,"E")] %>% as.data.frame
  }

  # test data
  x_test_mod = if (filter & !include_E) {
    x_test[, topgenes] %>% as.data.frame
  } else if (!filter & include_E) {
    x_test %>% as.data.frame
  } else if (!filter & !include_E) {
    x_test[,which(colnames(x_test) %ni% "E")] %>% as.data.frame
  } else if (filter & include_E) {
    x_test[, c(topgenes,"E")] %>% as.data.frame
  }

  # these are only derived on the main effects genes.. E is only included in the model
  PC_and_avg <- extractPC(x_train = x_train_mod[,gene_groups$gene],
                          colors = gene_groups$cluster,
                          x_test = x_test_mod[,gene_groups$gene],
                          nPC = nPC)

  n.clusters <- PC_and_avg$nclusters

  # this contains either the averages or PCs for each module in a data.frame
  clust_data <- switch(summary,
                       avg = PC_and_avg$averageExpr,
                       pc = PC_and_avg$PC)

  ml.formula <- if (include_interaction & include_E) {
    paste0("y_train ~","(",paste0(colnames(clust_data), collapse = "+"),")*E") %>% as.formula
  } else if (!include_interaction & include_E) {
    paste0("y_train ~",paste0(colnames(clust_data), collapse = "+"),"+E") %>% as.formula
  } else if (!include_interaction & !include_E) {
    paste0("y_train ~",paste0(colnames(clust_data), collapse = "+")) %>% as.formula
  }

  # this is the same as ml.formula, except without the response.. this is used for
  # functions that have the x = and y = input instead of a formula input
  model.formula <- if (include_interaction & include_E) {
    paste0("~ 0+(",paste0(colnames(clust_data), collapse = "+"),")*E") %>% as.formula
  } else if (!include_interaction & include_E) {
    paste0("~0+",paste0(colnames(clust_data), collapse = "+"),"+E") %>% as.formula
  } else if (!include_interaction & !include_E) {
    paste0("~0+",paste0(colnames(clust_data), collapse = "+")) %>% as.formula
  }

  # this is the design matrix based on model.formula
  X.model.formula <- model.matrix(model.formula, data = if (include_E) {
    cbind(clust_data,x_train_mod[,"E", drop = F])
  } else clust_data %>% as.data.frame)

  df <- X.model.formula %>% cbind(y_train) %>% as.data.frame()

  clust_train_model <- switch(model,
                              MARS = {

                                fitControl <-  trainControl(method = "cv",
                                                            # number = 25,
                                                            # repeats = 3,
                                                            verboseIter = FALSE)

                                marsGrid <- expand.grid(.degree = 1:2, .nprune = 1000)

                                switch(expFamily,
                                       gaussian = {
                                         mars_tuned <- train(X.model.formula,
                                                             y_train,
                                                             method = "earth",
                                                             trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                                                             tuneGrid = marsGrid,
                                                             trControl = fitControl)

                                         earth::earth(x = X.model.formula,
                                                      y = y_train,
                                                      keepxy = TRUE,
                                                      pmethod = "backward",
                                                      nk = 1000,
                                                      degree = mars_tuned$bestTune$degree,
                                                      trace = 4, nfold = 10) },
                                       binomial = {

                                         mars_tuned <- train(X.model.formula,
                                                             as.factor(y_train),
                                                             method = "earth",
                                                             trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                                                             glm=list(family=binomial),
                                                             tuneGrid = marsGrid,
                                                             trControl = fitControl)

                                         earth::earth(x = X.model.formula,
                                                      y = as.factor(y_train),
                                                      keepxy = TRUE,
                                                      pmethod = "backward",
                                                      nk = 1000,
                                                      glm=list(family=binomial),
                                                      degree = mars_tuned$bestTune$degree,
                                                      trace = 4, nfold = 10)
                                       })
                              }
  )


  # summary(clust_train_model)
  # ONLY Jaccard index can be calculated for MARS
  # since get.used.pred.names returns only the non-zero coefficients,
  # we give a coefficient of 1 here so that the stability calculation works
  # because it takes non-zero coef estimates

  coefs <- data.frame(get.used.pred.names(clust_train_model), rep(1, length(get.used.pred.names(clust_train_model))), stringsAsFactors = F) %>%
    as.data.table(keep.rownames = FALSE) %>%
    magrittr::set_colnames(c("Gene","coef.est"))

  if (stability) {
    # remove intercept for stability measures
    return(coefs)
  } else {

    non_zero_clusters <- coefs[coef.est != 0] %>%
      magrittr::use_series("Gene")

    # need to determine which of non_zero_cluters are main effects and which
    # are interactions
    non_zero_clusters_interactions <- grep(":",non_zero_clusters, value = T)

    # this checks if the environment is non-zero
    non_zero_environment <- grep("^E", non_zero_clusters, value = T,
                                 ignore.case = TRUE)
    non_zero_clusters_main_effects <- setdiff(non_zero_clusters,
                                              c(non_zero_clusters_interactions,
                                                non_zero_environment))

    # this includes the environment if the environment is non-zero
    n.non_zero_clusters <- coefs[coef.est != 0] %>%
      magrittr::use_series("Gene") %>%
      length

    # need to get the genes corresponding to the non-zero clusters
    # NOTE: this also includes non-zero cluster:Environment interactions

    # genes corresponding to non-zero main effect clusters
    # this list might not be unique if clust_type="Addon" because the same gene
    # can be in different clusters
    clust.S.hat.main <- gene_groups[cluster %in%
                                      as.numeric(
                                        unlist(
                                          stringr::str_extract_all(
                                            non_zero_clusters_main_effects, "(\\d+)$")
                                        )
                                      ),gene]

    # identical(gene_groups[cluster %in% c(3,12),gene],
    #           clust.S.hat.main)
    # identical(unique(clust.S.hat.main), clust.S.hat.main)
    # table(clust.S.hat.main)

    # this is the same as gene_groups, but the gene names contain E
    # so that we can extract the interactions corresponding to the chose clusters
    gene_groups_E <- copy(gene_groups)
    gene_groups_E[,gene:=paste0(gene,":E")]

    clust.S.hat.interaction <- gene_groups_E[cluster %in%
                                               as.numeric(
                                                 unlist(
                                                   stringr::str_extract_all(
                                                     stringr::str_extract_all(non_zero_clusters_interactions,"^.*?(?=:)"),
                                                     "(\\d+)$")
                                                 )
                                               ),gene]

    # this represents all the genes corresponding to the non-zero PC or avg
    # this list might not be unique if clust_type="Addon"
    # identical(unique(clust.S.hat), clust.S.hat)
    # I will double count if a model takes a gene more than once. ie.
    # if the same gene gets selected twice, then this will contribute 2 to the
    # number of non-zero estimated coefficients
    clust.S.hat <- c(clust.S.hat.main, non_zero_environment,
                     clust.S.hat.interaction)


    clust_data_test <- switch(summary,
                              avg = PC_and_avg$averageExprTest,
                              pc = PC_and_avg$PCTest)

    if (summary=="pc") { colnames(clust_data_test) <- colnames(clust_data) }

    # need intercept for prediction
    model.formula_test <- if (include_interaction & include_E) {
      paste0("~ 1+(",paste0(colnames(clust_data_test), collapse = "+"),")*E") %>% as.formula
    } else if (!include_interaction & include_E) {
      paste0("~1+",paste0(colnames(clust_data_test), collapse = "+"),"+E") %>% as.formula
    } else if (!include_interaction & !include_E) {
      paste0("~1+",paste0(colnames(clust_data_test), collapse = "+")) %>% as.formula
    }


    # this includes the intercept!
    X.model.formula_test <- model.matrix(model.formula_test,
                                         data = if (include_E) {
                                           cbind(clust_data_test,x_test_mod[,"E", drop = F])
                                         } else clust_data_test %>% as.data.frame)

    # True Positive Rate
    clust.TPR <- length(intersect(clust.S.hat, s0))/length(s0)

    # True negatives
    trueNegs <- setdiff(colnames(x_train_mod), s0)
    # identical(setdiff(colnames(x_train_mod), s0), setdiff(colnames(x_train), s0))

    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(x_train_mod),clust.S.hat)

    # how many of the terms identified by the model as zero, were actually zero
    # use to calculate correct sparsity as defined by Witten et al in the
    # Cluster Elastic Net paper Technometrics 2013
    C1 <- sum(modelIdentifyZero %in% trueNegs)
    C2 <- length(intersect(clust.S.hat, s0))
    clust.correct_sparsity <- (C1 + C2)/(ncol(x_train_mod))

    # this is from Interaction Screening for Ultrahigh Dimensional Data by ning hao and hao helen zhang
    true.interaction_names <- grep(":", s0, value = T)
    true.main_effect_names <- setdiff(s0, true.interaction_names)

    all.interaction_names <- grep(":", colnames(x_train_mod), value = T)
    all.main_effect_names <- setdiff(colnames(x_train_mod), all.interaction_names)

    true.negative_main_effects <- setdiff(all.main_effect_names, true.main_effect_names)
    true.negative_interaction_effects <- setdiff(all.interaction_names, true.interaction_names)

    (clust.correct_zeros_main_effects <- sum(setdiff(all.main_effect_names, c(clust.S.hat.main, non_zero_environment)) %in% true.negative_main_effects)/ length(true.negative_main_effects))
    (clust.correct_zeros_interaction_effects <- sum(setdiff(all.interaction_names, clust.S.hat.interaction) %in% true.negative_interaction_effects)/ length(true.negative_interaction_effects))

    (clust.incorrect_zeros_main_effects <- sum(setdiff(all.main_effect_names, c(clust.S.hat.main, non_zero_environment)) %in% true.main_effect_names)/ length(true.main_effect_names))
    (clust.incorrect_zeros_interaction_effects <- sum(setdiff(all.interaction_names, clust.S.hat.interaction) %in% true.interaction_names)/ length(true.interaction_names))

    # False Positive Rate = FP/(FP + TN) = FP / True number of 0 coefficients
    (clust.FPR <- sum(clust.S.hat %ni% s0)/(sum(clust.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs)))

    # Mean Squared Error

    # mars predict automatically adds the intercept if its not already in the dataset
    # the code below proves this
    # mars.pred1 <- predict(clust_train_model, newdata = X.model.formula_test[,-1], trace = 4)
    # mars.pred2 <- predict(clust_train_model, newdata = X.model.formula_test, trace = 4)
    # plot(mars.pred1, mars.pred2) ; identical(mars.pred1,mars.pred2)
    mars.pred <- predict(clust_train_model, newdata = X.model.formula_test, trace = 4)

    # Mean Squared Error
    (clust.mse <- crossprod(mars.pred - y_test)/length(y_test))

    # Root Mean Squared Error
    (clust.RMSE <- sqrt(crossprod(mars.pred - y_test)/length(y_test)))

    # remove intercept for prediction error formula given by ||X\beta - X\hat{\beta}||_2
    # given in Witten 2013 Cluster ENET paper in Technometrics
    # (clust.test_set_pred_error <- sqrt(crossprod(as.matrix(x_test_mod) %*% as.numeric(true_beta) - X.model.formula_test[,-1] %*% coefs$coef.est[-1])))

    # mse.null
    (mse_null <- crossprod(mean(y_test) - y_test)/length(y_test))

    if (expFamily == "binomial") {
      pred_response <- predict(clust_train_model, newdata = X.model.formula_test, trace = 4,
                               type = "response")

      clust.AUC <- pROC::roc(y_test,as.numeric(pred_response))$auc %>% as.numeric()

    }


    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    # clust.r2 <- (mse_null - clust.mse)/mse_null

    # clust.adj.r2 <- 1 - (1 - clust.r2)*(nrow(x_test) - 1)/(nrow(x_test) - n.non_zero_clusters - 1)


    ls <- switch(expFamily,
                 gaussian = {list(clust.mse = as.numeric(clust.mse),
                                  clust.RMSE = as.numeric(clust.RMSE),
                                  clust.S.hat = length(clust.S.hat),
                                  clust.TPR = clust.TPR,
                                  clust.FPR = clust.FPR,
                                  clust.correct_sparsity,
                                  clust.correct_zeros_main_effects,
                                  clust.correct_zeros_interaction_effects,
                                  clust.incorrect_zeros_main_effects,
                                  clust.incorrect_zeros_interaction_effects)},
                 binomial = {list(clust.mse = as.numeric(clust.mse),
                                  clust.RMSE = as.numeric(clust.RMSE),
                                  clust.AUC = clust.AUC,
                                  clust.S.hat = length(clust.S.hat),
                                  clust.TPR = clust.TPR,
                                  clust.FPR = clust.FPR,
                                  clust.correct_sparsity,
                                  clust.correct_zeros_main_effects,
                                  clust.correct_zeros_interaction_effects,
                                  clust.incorrect_zeros_main_effects,
                                  clust.incorrect_zeros_interaction_effects)

                 })

    names(ls) <- switch(expFamily,
                        gaussian = {
                          c(paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))},
                        binomial = {
                          c(paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_AUC"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                            paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
                        })
    return(ls)

  }
}
