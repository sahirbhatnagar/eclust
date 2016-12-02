# number_pc=2 will not work if there is only one gene in a module!
# number_pc=2 will not work if there is only one gene in a module!

#' Fit Penalized Regression Models on Simulated Cluster Summaries
#'
#' @description This function creates summaries of the given clusters (e.g. 1st
#'   PC or average), and then fits a penalized regression model on those
#'   summaries. To be used with simulated data where the 'truth' is known i.e.,
#'   you know which features are associated with the response. This function was
#'   used to produce the simulation results in Bhatnagar et al. 2016
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
#' @param s0 chracter vector of
#' @param summary the summary of each cluster. Can be the principal component or
#'   average. Default is \code{summary = "pc"} which takes the first
#'   \code{number_pc} principal components. Currently a maximum of 2 principal
#'   components can be chosen.
#' @param model Regression model to be fit on cluster summaries. Default is
#'   \code{model="lasso"}. Can also be
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
#'
#'
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
#'   If \code{stability = FALSE} then returns a vector with the following elements:
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#'  \item{}{}
#' @export
#'
#' @examples
#' dsfds
s_clust_reg <- function(x_train,
                        x_test,
                        y_train,
                        y_test,
                        s0,
                        summary = c("pc","avg"),
                        model = c("lasso", "scad", "mcp", "elasticnet"),
                        exp_family = c("gaussian","binomial"),
                        gene_groups,
                        topgenes = NULL,
                        stability = F,
                        filter = F,
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
