library(eclust)
library(magrittr)

parametersDf <- expand.grid(rho = c(0.90),
                            p = c(1000),
                            # SNR = c(0.2,1,2),
                            SNR = 1,
                            n = c(200), # this is the total train + test sample size
                            # nActive = c(300), # must be even because its being split among two modules
                            #n0 = 200,
                            cluster_distance = c("tom","corr"),
                            Ecluster_distance = c("difftom", "diffcorr"),
                            rhoOther = 0.6,
                            betaMean = c(2),
                            alphaMean = c(1),
                            betaE = 3,
                            includeInteraction = TRUE,
                            includeStability = TRUE,
                            distanceMethod = "euclidean",
                            clustMethod = "hclust",
                            #cutMethod = "gap",
                            cutMethod = "dynamic",
                            agglomerationMethod = "average",
                            # agglomerationMethod = "ward.D",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.10)
parametersDf <- parametersDf[which(parametersDf$cluster_distance=="tom" & parametersDf$Ecluster_distance=="difftom" |
                                     parametersDf$cluster_distance=="corr" & parametersDf$Ecluster_distance=="diffcorr"),]
parameterIndex = 2
simulationParameters <- parametersDf[parameterIndex,, drop = F]

print(simulationParameters)

## ---- generate-data ----
message("generating data")

p <- simulationParameters[,"p"];
n <- simulationParameters[,"n"];
n0 <- simulationParameters[,"n0"];
SNR <- simulationParameters[,"SNR"]
n1 <- n - n0
cluster_distance <- simulationParameters[,"cluster_distance"]
Ecluster_distance <- simulationParameters[,"Ecluster_distance"]
rhoOther <- simulationParameters[,"rhoOther"];
betaMean <- simulationParameters[,"betaMean"];
betaE <- simulationParameters[,"betaE"]
alphaMean <- simulationParameters[,"alphaMean"];
rho <- simulationParameters[,"rho"];
nActive <- simulationParameters[,"nActive"];
includeInteraction <- simulationParameters[,"includeInteraction"]
includeStability <- simulationParameters[,"includeStability"]
distanceMethod <- simulationParameters[,"distanceMethod"]
clustMethod <- simulationParameters[,"clustMethod"]
cutMethod <- simulationParameters[,"cutMethod"]
agglomerationMethod <- simulationParameters[,"agglomerationMethod"]
K.max <- simulationParameters[,"K.max"]
B <- simulationParameters[,"B"]

# in this simulation its blocks 3 and 4 that are important
# leaveOut:  optional specification of modules that should be left out
# of the simulation, that is their genes will be simulated as unrelated
# ("grey"). This can be useful when simulating several sets, in some which a module
# is present while in others it is absent.
d0 <- s_modules(n = n0, p = p, rho = c(0,0), exposed = FALSE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.01,
                maxCor = 1,
                corPower = 1,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = FALSE,
                leaveOut = 1:4)

d1 <- s_modules(n = n1, p = p, rho = c(rho, rho), exposed = TRUE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.4,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = FALSE)

# these should be the same. if they arent, its because I removed the red and
# green modules from the E=0 group
truemodule0 <- d0$setLabels
t0 <- table(truemodule0)
truemodule1 <- d1$setLabels
t1 <- table(truemodule1)
table(truemodule0,truemodule1)

# Convert labels to colors for plotting
moduleColors <- WGCNA::labels2colors(truemodule1)
table(moduleColors, truemodule1)

# does not contain E, needs to bind unexposed first
X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))

dim(X)


# pheatmap::pheatmap(cor(X))
# pheatmap::pheatmap(cor(d1$datExpr))
# pheatmap::pheatmap(cor(d0$datExpr))
# pheatmap(cor(d1$datExpr)-cor(d0$datExpr))
# pheatmap(WGCNA::TOMsimilarityFromExpr(X))
# pheatmap(WGCNA::TOMsimilarityFromExpr(d1$datExpr))
# pheatmap(WGCNA::TOMsimilarityFromExpr(d1$datExpr)-WGCNA::TOMsimilarityFromExpr(d0$datExpr))

# betaMainEffect <- vector("double", length = p)
# betaMainInteractions <- vector("double", length = p)
# # first assign random uniform to every gene in cluster 3 and 4,
# # then randomly remove so that thers only nActive left
# betaMainEffect[which(truemodule1 %in% 3:4)] <- runif(sum(truemodule1 %in% 3:4),
#                                                      betaMean - 0.1, betaMean + 0.1)
#
# # randomly set some coefficients to 0 so that there are only nActive non zero
# betaMainEffect <- replace(betaMainEffect,
#                           sample(which(truemodule1 %in% 3:4), sum(truemodule1 %in% 3:4) - nActive,
#                                  replace = FALSE), 0)
#
# betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)
#
# beta <- c(betaMainEffect,
#           betaE,
#           betaMainInteractions)
# plot(beta)

betaMainEffect <- vector("double", length = p)

betaMainInteractions <- vector("double", length = p)

# the first nActive/2 in the 3rd block are active
betaMainEffect[which(truemodule1 %in% 3)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

# the first nActive/2 in the 4th block are active
betaMainEffect[which(truemodule1 %in% 4)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean+2 - 0.1, betaMean+2 + 0.1)

betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)

# must be in this order!!!! main effects, E, and then interactions... this order is being used
# by the generate_data function
beta <- c(betaMainEffect,
          betaE,betaMainInteractions)

plot(beta)

simulated_data <- s_generate_data(p = p,
                          X = X,
                          beta = beta,
                          include_interaction = includeInteraction,
                          cluster_distance = cluster_distance,
                          n = n,
                          n0 = n0,
                          eclust_distance = Ecluster_distance,
                          signal_to_noise_ratio = SNR,
                          distance_method = distanceMethod,
                          cluster_method = clustMethod,
                          cut_method = cutMethod,
                          agglomeration_method = agglomerationMethod,
                          K.max = K.max, B = B, nPC = 1)

# devtools::use_data(simulated_data)


simulated_data$similarity %>% dim
simulated_data$DT %>% dim
simulated_data$DT %>% str
pryr::object_size(simulated_data$DT)

pryr::object_size(simulated_data$DT[,1:1002])

simdata <- simulated_data$DT[,c(1,1002,2:1001)]
simdata[1:5,1:5]
simdata %>% dim
devtools::use_data(simdata, overwrite = TRUE)
devtools::check()
