# #' Function to generate heatmap
# #'
# #' @description x matrix of true correlation (P x P matrix where P is the number
# #'   of genes). Must be object of class similarity
#
# plot.similarity <- function(x,
#                             color = viridis(100),
#                             truemodule,
#                             active, ...){
#
#   if (missing(active)) {
#     annotation_col <- data.frame(
#       module = factor(truemodule,
#                       labels = c("Grey","Turquoise","Blue","Red",
#                                  "Green","Yellow")))
#
#     rownames(annotation_col) <- dimnames(x)[[2]]
#     ann_colors <- list(
#       module = c(Turquoise = "turquoise",
#                  Blue = "blue",
#                  Red = "red",
#                  Green = "green",
#                  Yellow = "yellow",
#                  Grey = "grey90")
#     )
#   } else {
#     annotation_col <- data.frame(
#       module = factor(truemodule,
#                       labels = c("Grey","Turquoise","Blue","Red",
#                                  "Green","Yellow")),
#       active = factor(active,
#                       levels = c(0,1),
#                       labels = c("no","yes")))
#
#     rownames(annotation_col) <- dimnames(x)[[2]]
#     ann_colors <- list(
#       module = c(Turquoise = "turquoise",
#                  Blue = "blue",
#                  Red = "red",
#                  Green = "green",
#                  Yellow = "yellow",
#                  Grey = "grey90"),
#       active = c(no = "black",
#                  yes = "orange")
#     )
#   }
#
#   pheatmap(x,
#            show_rownames = F, show_colnames = F,
#            color = color,
#            annotation_col = annotation_col,
#            annotation_row = annotation_col,
#            annotation_colors = ann_colors,
#            annotation_names_row = FALSE,
#            annotation_names_col = TRUE,
#            drop_levels = FALSE,
#            annotation_legend = TRUE, ...)
#
#
#   # breaks = seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 101) ,
#   # legend_breaks = round(seq(min(min_max_heat$V2), max(min_max_heat$V1),
#   # length.out = 12),1),
#   # legend_labels = round(seq(min(min_max_heat$V2), max(min_max_heat$V1),
#   # length.out = 12),1),
#   # drop_levels = FALSE, ...)
# }
