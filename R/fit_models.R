# nPC=2 will not work if there is only one gene in a module!
clust_fun <- function(x_train,
                      x_test,
                      y_train,
                      y_test,
                      s0,
                      summary = c("pc","avg"),
                      model = c("lm", "lasso", "scad", "mcp", "elasticnet", "shim"),
                      gene_groups,
                      true_beta = NULL,
                      topgenes = NULL,
                      stability = F,
                      filter = F,
                      include_E = F,
                      include_interaction = F,
                      p = 1000,
                      filter_var = F,
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
  # model = "lasso"; summary = "pc"; topgenes = NULL; clust_type="clust"; nPC = 1

  clust_type <- match.arg(clust_type)
  summary <- match.arg(summary)
  model <- match.arg(model)

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
                              lm = lm(ml.formula , data = df),
                              lasso = {if (n.clusters != 1) {
                                cv.glmnet(x = X.model.formula, y = y_train, alpha = 1)
                              } else NA },
                              elasticnet = {if (n.clusters != 1) {
                                cv.glmnet(x = X.model.formula, y = y_train, alpha = 0.5)
                              } else NA },
                              scad = {
                                cv.ncvreg(X = X.model.formula, y = y_train,
                                          family = "gaussian", penalty = "SCAD")
                              },
                              mcp = {
                                cv.ncvreg(X = X.model.formula, y = y_train,
                                          family = "gaussian", penalty = "MCP")
                              },
                              shim = {
                                require(doMC)
                                registerDoMC(cores = 4)
                                cv.shim(x = X.model.formula, y = y_train,
                                        main.effect.names = c(colnames(clust_data), if (include_E) "E"),
                                        interaction.names = setdiff(colnames(X.model.formula),c(colnames(clust_data),"E")),
                                        max.iter = 120, initialization.type = "ridge",
                                        verbose = FALSE, parallel = TRUE, nfolds = 10)
                              })
  # plot(clust_train_model)

  # here we give the coefficient stability on the clusters and not the individual genes
  coefs <- switch(model,
                  lm = data.table::data.table(Gene = names(clust_train_model$coefficients),
                                              coef.est = coef(clust_train_model)),
                  lasso = {
                    # need to return all 0's if there is only 1 cluster since lasso
                    # wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {
                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
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
                        as.data.table(keep.rownames = TRUE) %>%
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
                        as.data.table(keep.rownames = TRUE) %>%
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
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  shim = {
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {

                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  })

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


    ls <- list(clust.mse = clust.mse,
               clust.RMSE,
               # clust.r2 = clust.r2,
               # clust.adj.r2 = clust.adj.r2,
               clust.S.hat = length(clust.S.hat),
               clust.TPR = clust.TPR,
               clust.FPR = clust.FPR,
               clust.correct_sparsity = clust.correct_sparsity,
               clust.correct_zeros_main_effects,
               clust.correct_zeros_interaction_effects,
               clust.incorrect_zeros_main_effects,
               clust.incorrect_zeros_interaction_effects,
               n.clusters)


    names(ls) <- c(paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                   # paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   # paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
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
    return(ls)

  }
}


pen_fun <- function(x_train,
                    x_test,
                    y_train,
                    y_test,
                    s0,
                    model,
                    true_beta,
                    topgenes = NULL,
                    stability = F,
                    filter = F,
                    include_E = F,
                    include_interaction = F,
                    p = 1000,
                    filter_var = F){

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

  print(paste(model,"filter = ", filter, "filter_var = ",filter_var, "include_E = ", include_E, "include_interaction = ", include_interaction, sep = " "))

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
                      lasso = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 1),
                      elasticnet = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 0.5),
                      ridge = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 0),
                      scad = ncvreg::cv.ncvreg(X = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train,
                        family = "gaussian", penalty = "SCAD"),
                      mcp = ncvreg::cv.ncvreg(X = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train,
                        family = "gaussian", penalty = "MCP")
  )

  # plot(pen_model)
  # oracle penalization model
  # pen_model_oracle <- switch(model,
  #                            lasso = glmnet::cv.glmnet(x = as.matrix(x_train[,s0]), y = y_train, alpha = 1),
  #                            elasticnet = glmnet::cv.glmnet(x = as.matrix(x_train[,s0]), y = y_train, alpha = 0.5),
  #                            ridge = glmnet::cv.glmnet(x = as.matrix(x_train[,s0]), y = y_train, alpha = 0),
  #                            scad = ncvreg::cv.ncvreg(X = x_train[,s0], y = y_train,
  #                                                     family = "gaussian", penalty = "SCAD"),
  #                            mcp = ncvreg::cv.ncvreg(X = x_train[,s0], y = y_train,
  #                                                    family = "gaussian", penalty = "MCP")
  # )

  # here we give the coefficient stability on the individual genes
  coefs <- coef(pen_model, s = "lambda.min") %>%
    as.matrix %>%
    as.data.table(keep.rownames = TRUE) %>%
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

    # pen.pred.oracle <- if (model %in% c("lasso","elasticnet","ridge")) {
    #   predict(pen_model_oracle, newx = x_test[,s0], s = "lambda.min") } else if (model %in% c("scad","mcp")) {
    #     predict(pen_model_oracle, X = x_test[,s0],
    #             lambda = pen_model_oracle$lambda.min)
    #   }

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

    # # False Positive Rate
    # pen.FPR <- sum(pen.S.hat %ni% s0)/(p - length(s0))

    # Mean Squared Error
    (pen.mse <- crossprod(pen.pred - y_test)/length(y_test))

    # Root Mean Squared Error
    (pen.RMSE <- sqrt(crossprod(pen.pred - y_test)/length(y_test)))

    # Mean Squared Error Oracle
    # pen.mse.oracle <- crossprod(pen.pred.oracle - y_test)/length(y_test)

    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)

    # sqrt(mse_null)

    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    # pen.r2 <- (mse_null - pen.mse)/mse_null

    # pen.adj.r2 <- 1 - (1 - pen.r2)*(length(y_test) - 1)/(length(y_test) - length(pen.S.hat) - 1)

    # model error
    # identical(true_beta %>% rownames(),coefs[["Gene"]])
    # pen.model.error <- {(true_beta - coefs[["coef.est"]]) %>% t} %*% WGCNA::cor(x_test[,coefs[["Gene"]]]) %*% (true_beta - coefs[["coef.est"]])

    ls <- list(pen.mse = as.numeric(pen.mse),
               pen.RMSE = as.numeric(pen.RMSE),
               # pen.r2 = as.numeric(pen.r2),
               # pen.adj.r2 = as.numeric(pen.adj.r2),
               pen.S.hat = length(pen.S.hat),
               pen.TPR = pen.TPR,
               pen.FPR = pen.FPR,
               # pen.relative.mse = pen.mse/pen.mse.oracle,
               # pen.model.error = pen.model.error,
               correct_sparsity,
               pen.correct_zeros_main_effects,
               pen.correct_zeros_interaction_effects,
               pen.incorrect_zeros_main_effects,
               pen.incorrect_zeros_interaction_effects
    )
    names(ls) <- c(paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                   # paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   # paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                   paste0("pen_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
    return(ls)
  }

}



mars_fun <- function(x_train,
                     x_test,
                     y_train,
                     y_test,
                     s0,
                     model = c("MARS"),
                     true_beta,
                     topgenes = NULL,
                     stability = F,
                     filter = F,
                     include_E = F,
                     include_interaction = F,
                     p = 1000,
                     filter_var = F, ...){

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

  print(paste(model,"filter = ", filter, "filter_var = ",filter_var, "include_E = ", include_E, "include_interaction = ", include_interaction, sep = " "))

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
                       MARS = earth::earth(x = as.matrix(x_train),
                                           y = y_train,
                                           keepxy = TRUE,
                                           pmethod = "backward",
                                           nk = 1000,
                                           degree = 3,
                                           trace = 4))
  # ,
  # nfold = 5))

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
    as.data.table(keep.rownames = FALSE) %>%
    magrittr::set_colnames(c("Gene","coef.est"))

  if (stability) {
    # remove intercept for stability measures
    return(coefs)
  } else {

    mars.S.hat <- get.used.pred.names(mars_model)
    mars.S.hat.interaction <- grep(":", mars.S.hat, value = T)
    mars.S.hat.main <- setdiff(mars.S.hat, mars.S.hat.interaction)

    mars.pred <- predict(mars_model, newdata = x_test, trace = 4)

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

    # Mean Squared Error Oracle
    # mars.mse.oracle <- crossprod(mars.pred.oracle - y_test)/length(y_test)

    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)

    # sqrt(mse_null)

    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    # mars.r2 <- (mse_null - mars.mse)/mse_null

    # mars.adj.r2 <- 1 - (1 - mars.r2)*(length(y_test) - 1)/(length(y_test) - length(mars.S.hat) - 1)

    # model error
    # identical(true_beta %>% rownames(),coefs[["Gene"]])
    # mars.model.error <- {(true_beta - coefs[["coef.est"]]) %>% t} %*% WGCNA::cor(x_test[,coefs[["Gene"]]]) %*% (true_beta - coefs[["coef.est"]])

    ls <- list(mars.mse = as.numeric(mars.mse),
               mars.RMSE = as.numeric(mars.RMSE),
               # mars.r2 = as.numeric(mars.r2),
               # mars.adj.r2 = as.numeric(mars.adj.r2),
               mars.S.hat = length(mars.S.hat),
               mars.TPR = mars.TPR,
               mars.FPR = mars.FPR,
               # mars.relative.mse = mars.mse/mars.mse.oracle,
               # mars.model.error = mars.model.error,
               correct_sparsity,
               mars.correct_zeros_main_effects,
               mars.correct_zeros_interaction_effects,
               mars.incorrect_zeros_main_effects,
               mars.incorrect_zeros_interaction_effects
    )
    names(ls) <- c(paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                   # paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   # paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
    return(ls)
  }
}


mars_clust_fun <- function(x_train,
                           x_test,
                           y_train,
                           y_test,
                           s0,
                           summary = c("pc","avg"),
                           model = c("MARS"),
                           gene_groups,
                           true_beta = NULL,
                           topgenes = NULL,
                           stability = F,
                           filter = F,
                           include_E = F,
                           include_interaction = F,
                           p = 1000,
                           filter_var = F,
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

  clust_type <- match.arg(clust_type)
  summary <- match.arg(summary)
  model <- match.arg(model)

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
                              MARS = earth::earth(x = X.model.formula,
                                                  y = y_train,
                                                  keepxy = TRUE,
                                                  pmethod = "backward",
                                                  nk = 1000,
                                                  degree = 3,
                                                  trace = 4))
  # ,
  # nfold = 5))

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

    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    # clust.r2 <- (mse_null - clust.mse)/mse_null

    # clust.adj.r2 <- 1 - (1 - clust.r2)*(nrow(x_test) - 1)/(nrow(x_test) - n.non_zero_clusters - 1)


    ls <- list(clust.mse = clust.mse,
               clust.RMSE,
               # clust.r2 = clust.r2,
               # clust.adj.r2 = clust.adj.r2,
               clust.S.hat = length(clust.S.hat),
               clust.TPR = clust.TPR,
               clust.FPR = clust.FPR,
               clust.correct_sparsity = clust.correct_sparsity,
               clust.correct_zeros_main_effects,
               clust.correct_zeros_interaction_effects,
               clust.incorrect_zeros_main_effects,
               clust.incorrect_zeros_interaction_effects,
               n.clusters)


    names(ls) <- c(paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                   # paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   # paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
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
    return(ls)

  }
}
