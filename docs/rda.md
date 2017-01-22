Real data analysis functions (`r_`)
===================================

We will use the `data(tcgaov)` dataset included in this package, which contains a subset of the TCGA mRNA Ovarian serous cystadenocarcinoma data generated using Affymetrix HTHGU133a arrays. See `?tcgaov` for details about the data.

In the example below we use the `r_cluster_data` to create the environment based clusters, and their summaries. We then use the `r_prepare_data` function to get it into proper form for regression routines such as `earth::earth`, `glmnet::cv.glmnet`, and `ncvreg::ncvreg`.

Extract the relevant data
-------------------------

``` r
# load the data
data("tcgaov")
tcgaov[1:5,1:6, with = FALSE]
```

    ##              rn subtype E status   OS    ABCA8
    ## 1: TCGA-04-1331       4 1      1 1336 3.684824
    ## 2: TCGA-04-1332       1 0      1 1247 7.892982
    ## 3: TCGA-04-1335       3 1      1   55 5.193188
    ## 4: TCGA-04-1336       3 1      0 1495 3.055437
    ## 5: TCGA-04-1337       1 0      1   61 3.149427

``` r
# use log survival as the response
Y <- log(tcgaov[["OS"]])

# specify the environment variable
E <- tcgaov[["E"]]

# specify the matrix of genes only
genes <- as.matrix(tcgaov[,-c("OS","rn","subtype","E","status"),with = FALSE])

# for this example the training set will be all subjects.
# change `p` argument to create a train and test set.
trainIndex <- drop(caret::createDataPartition(Y, p = 1, list = FALSE, times = 1))
testIndex <- trainIndex
```

Cluster the data and calculate cluster representations
------------------------------------------------------

We cluster the genes using the correlation matrix (specified by `cluster_distance = "corr"`) and the difference of the exposure dependent correlation matrices (specified by `eclust_distance = "diffcorr"`)

``` r
cluster_res <- r_cluster_data(data = genes,
                              response = Y,
                              exposure = E,
                              train_index = trainIndex,
                              test_index = testIndex,
                              cluster_distance = "corr",
                              eclust_distance = "diffcorr",
                              measure_distance = "euclidean",
                              clustMethod = "hclust",
                              cutMethod = "dynamic",
                              method = "average",
                              nPC = 1,
                              minimum_cluster_size = 30)
```

    ## 

    ##  ..cutHeight not given, setting it to 10.9  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

    ## Calculating number of environment clusters based on diffcorr

    ##  ..cutHeight not given, setting it to 0.923  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

    ## There are 7 clusters derived from the corr similarity matrix

    ## There are 4 clusters derived from the diffcorr environment similarity matrix

    ## There are a total of 11 clusters derived from the corr
    ##                   similarity matrix and the diffcorr environment similarity matrix

``` r
# the number of clusters determined by the similarity matrices specified
# in the cluster_distance and eclust_distance arguments. This will always be larger
# than cluster_res$clustersAll$nclusters which is based on the similarity matrix
# specified in the cluster_distance argument
cluster_res$clustersAddon$nclusters
```

    ## [1] 11

``` r
# the number of clusters determined by the similarity matrices specified
# in the cluster_distance argument only
cluster_res$clustersAll$nclusters
```

    ## [1] 7

``` r
# what's in the cluster_res object
names(cluster_res)
```

    ## [1] "clustersAddon"            "clustersAll"             
    ## [3] "etrain"                   "clustersAddonMembership" 
    ## [5] "clustersAllMembership"    "clustersEclustMembership"

Prepare data for input in any regression routine
------------------------------------------------

Now we use the `r_prepare_data` function, where we are using the average expression from each cluster as feaand their interaction with E as features in the regression model:

``` r
# prepare data for use with earth function
avg_eclust_interaction <- r_prepare_data(
  data = cbind(cluster_res$clustersAddon$averageExpr, 
               Y = Y[trainIndex],
               E = E[trainIndex]),
  response = "Y", exposure = "E")

head(avg_eclust_interaction[["X"]])
```

    ##       avg1     avg2     avg3     avg4     avg5     avg6     avg7     avg8
    ## 1 5.567323 4.603389 5.645931 5.463546 6.330674 5.017668 5.723670 4.537860
    ## 2 6.416015 5.135543 5.915361 5.412958 6.539799 5.514176 6.330553 4.949290
    ## 3 4.751661 4.428877 6.064181 4.553495 5.647454 5.434861 4.854128 4.923794
    ## 4 4.488082 5.360132 6.184347 4.440893 7.747494 5.475239 5.509369 4.919032
    ## 5 6.661257 5.316210 5.873579 4.951189 5.769990 5.254726 6.868405 4.618736
    ## 6 6.206893 5.783308 6.017647 4.862209 6.001902 5.370070 6.680607 4.444538
    ##       avg9    avg10    avg11 E   avg1:E   avg2:E   avg3:E   avg4:E
    ## 1 5.504560 5.769372 4.604983 1 5.567323 4.603389 5.645931 5.463546
    ## 2 6.195038 6.119165 4.787305 0 0.000000 0.000000 0.000000 0.000000
    ## 3 4.829277 5.675971 4.227867 1 4.751661 4.428877 6.064181 4.553495
    ## 4 4.855203 6.117120 4.754397 1 4.488082 5.360132 6.184347 4.440893
    ## 5 6.282899 5.814544 5.048062 0 0.000000 0.000000 0.000000 0.000000
    ## 6 6.000493 6.044401 5.190256 0 0.000000 0.000000 0.000000 0.000000
    ##     avg5:E   avg6:E   avg7:E   avg8:E   avg9:E  avg10:E  avg11:E
    ## 1 6.330674 5.017668 5.723670 4.537860 5.504560 5.769372 4.604983
    ## 2 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
    ## 3 5.647454 5.434861 4.854128 4.923794 4.829277 5.675971 4.227867
    ## 4 7.747494 5.475239 5.509369 4.919032 4.855203 6.117120 4.754397
    ## 5 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
    ## 6 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000

Fit a regression model
----------------------

At this stage, you can decide which regression model to use. Here we choose the MARS model from the `earth` package, but you may choose regression models from any number of packages (e.g. see the [extensive list of models](https://topepo.github.io/caret/available-models.html) of models available in the `caret` package).

``` r
fit_earth <- earth::earth(x = avg_eclust_interaction[["X"]], 
                          y = avg_eclust_interaction[["Y"]], 
                          pmethod = "backward", 
                          keepxy = TRUE, 
                          degree = 2, 
                          trace = 1, 
                          nk = 1000)
```

    ## x[511,23] with colnames avg1 avg2 avg3 avg4 avg5 avg6 avg7 avg8 avg9 avg10 avg11 ...
    ## y[511,1] with colname avg_eclust_interaction[["Y"]]
    ## Forward pass term 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 
    ##      30, 32, 34, 36, 38
    ## RSq changed by less than 0.001 at 37 terms, 35 terms used (DeltaRSq 0)
    ## After forward pass GRSq -0.176 RSq 0.193
    ## Prune method "backward" penalty 3 nprune null: selected 6 of 35 terms, and 4 of 23 preds
    ## After pruning pass GRSq 0.01 RSq 0.0579

``` r
coef(fit_earth)
```

    ##                         (Intercept) h(avg9-5.58701) * h(6.54361-avg5:E) 
    ##                            6.718508                            1.070107 
    ##   h(avg1-5.99657) * h(avg7-6.84195)   h(avg1-5.99657) * h(avg7-6.70675) 
    ##                           67.850654                          -39.786576 
    ##   h(avg1-5.99657) * h(avg7-6.98039)                     h(avg9-5.50456) 
    ##                          -28.720223                           -5.293762

You can also install the `plotmo` package to visualise the relationships between the hinge functions and the response using `plotmo::plotmo(fit_earth)`.

Determine the features that have been selected
----------------------------------------------

The `u_extract_selected_earth` is a utility function in this package to extract the selected predictors from the MARS model:

``` r
u_extract_selected_earth(fit_earth)
```

    ## [1] "avg1"   "avg7"   "avg9"   "avg5:E"

We that genes in clusters 1, 7 and 9 were selected. We also see that the interaction between the genes in cluster 5 and the environment was selected and has the highest variable importance. We can see the genes involved using the `cluster_res$clustersAddonMembership` object:

``` r
# Genes in cluster 5
cluster_res$clustersAddonMembership[cluster %in% 5]
```

    ##         gene cluster module
    ##  1: APOBEC3G       5  green
    ##  2:    APOL6       5  green
    ##  3:   CXCL10       5  green
    ##  4:   CXCL11       5  green
    ##  5:     GBP1       5  green
    ##  6:     HCP5       5  green
    ##  7:     IL15       5  green
    ##  8:    PSMB9       5  green
    ##  9:  RARRES3       5  green
    ## 10:     TAP1       5  green
    ## 11:     BST2       5  green
    ## 12:   BTN3A2       5  green
    ## 13:       C2       5  green
    ## 14:     CBR3       5  green
    ## 15: CCDC109B       5  green
    ## 16:     HPSE       5  green
    ## 17:  HTATIP2       5  green
    ## 18:    IFI16       5  green
    ## 19:    IFI27       5  green
    ## 20:    IFI35       5  green
    ## 21:    IFI44       5  green
    ## 22:   IFI44L       5  green
    ## 23:    IFIH1       5  green
    ## 24:    IFIT1       5  green
    ## 25:    IFIT2       5  green
    ## 26:    IFIT3       5  green
    ## 27:   IFITM1       5  green
    ## 28:   IL15RA       5  green
    ## 29:     IRF1       5  green
    ## 30:     IRF7       5  green
    ## 31:    ISG15       5  green
    ## 32:    ISG20       5  green
    ## 33:    LAMP3       5  green
    ## 34:     LAP3       5  green
    ## 35:   LGALS9       5  green
    ## 36:      MX1       5  green
    ## 37:     OAS1       5  green
    ## 38:     OAS2       5  green
    ## 39:     OAS3       5  green
    ## 40:   PLSCR1       5  green
    ## 41:   PSMB10       5  green
    ## 42:    PSMB8       5  green
    ## 43:    PSME2       5  green
    ## 44:    RSAD2       5  green
    ## 45:     RTP4       5  green
    ## 46:    SAMD9       5  green
    ## 47:  SLC15A3       5  green
    ## 48:    SP100       5  green
    ## 49:     TAP2       5  green
    ## 50:    TAPBP       5  green
    ## 51:     TLR3       5  green
    ## 52:  TMEM140       5  green
    ## 53:  TNFAIP8       5  green
    ## 54:  TNFSF10       5  green
    ## 55:   TRIM14       5  green
    ## 56:   UBE2L6       5  green
    ## 57:     XAF1       5  green
    ##         gene cluster module

``` r
# variable importance
earth::evimp(fit_earth)
```

    ##              nsubsets   gcv    rss
    ## avg5:E              4  84.7  100.0
    ## avg11-unused        3 -49.6   73.9
    ## avg1                2 100.0>  81.8>
    ## avg7                2 100.0   81.8
    ## avg9                2 100.0   81.8
