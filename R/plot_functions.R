#' Function to generate heatmap
#'
#' @description Plots a heatmap of a similarity matrix such as a correlation
#'   matrix or a TOM matrix. This function is a plotting method for an object of
#'   class similarity. These objects are returned by the
#'   \code{\link{s_generate_data}} and \code{\link{s_generate_data_mars}}
#'   functions
#' @param x an object of class similarity. This is a p x p symmetric matix such
#'   as a correlation matrix or a TOM matrix, where p is the number of genes
#' @param color colors for the heatmap. By default it uses the \code{viridis}
#'   color scheme. The \code{viridis} package needs to be installed.
#' @param truemodule a numeric vector of length p where p is the number of
#'   genes, giving the module membership. By default, 0 = Grey, 1 = Turquoise, 2
#'   = Blue, 3 = Red, 4 = Green, and 5 = Yellow. This information is used for
#'   annotating the heatmap
#' @param active a binary vector of length p (where p is the number of genes)
#'   where 0 means that gene is not related to the response, and 1 means that
#'   the gene is associated to the response.
#' @param ... other arguments passed to the pheatmap function
#'
#' @note this function is only meant to be used with output from the
#'   \code{\link{s_generate_data}} and \code{\link{s_generate_data_mars}}
#'   functions, since it assumes a fixed number of modules.
#' @return a heatmap of a similarity matrix
#' @examples
#' \dontrun{
#' corrX <- cor(simdata[,c(-1,-2)])
#' class(corrX) <- append(class(corrX), "similarity")
#' plot(corrX, truemodule = c(rep(1:5, each=150), rep(0, 250)))
#' }
#' @export

plot.similarity <- function(x,
                            color = viridis::viridis(100),
                            truemodule,
                            active, ...){

  if (missing(active)) {
    annotation_col <- data.frame(
      module = factor(truemodule,
                      labels = c("Grey","Turquoise","Blue","Red",
                                 "Green","Yellow")))

    rownames(annotation_col) <- dimnames(x)[[2]]
    ann_colors <- list(
      module = c(Turquoise = "turquoise",
                 Blue = "blue",
                 Red = "red",
                 Green = "green",
                 Yellow = "yellow",
                 Grey = "grey90")
    )
  } else {
    annotation_col <- data.frame(
      module = factor(truemodule,
                      labels = c("Grey","Turquoise","Blue","Red",
                                 "Green","Yellow")),
      active = factor(active,
                      levels = c(0,1),
                      labels = c("no","yes")))

    rownames(annotation_col) <- dimnames(x)[[2]]
    ann_colors <- list(
      module = c(Turquoise = "turquoise",
                 Blue = "blue",
                 Red = "red",
                 Green = "green",
                 Yellow = "yellow",
                 Grey = "grey90"),
      active = c(no = "black",
                 yes = "orange")
    )
  }

  pheatmap::pheatmap(x,
           show_rownames = F, show_colnames = F,
           color = color,
           annotation_col = annotation_col,
           annotation_row = annotation_col,
           annotation_colors = ann_colors,
           annotation_names_row = FALSE,
           annotation_names_col = TRUE,
           drop_levels = FALSE,
           annotation_legend = TRUE, ...)

}
