#' @title fast PCA
#' @description perform PCA from xcms object
#' @author Yonghui Dong
#' @param dat sample ion intensity matrix, row sample, column feature.
#' @param Group sample group information
#' @param centering centering, default = TRUE
#' @param scaling scaling method, default is scaling = "none". You can choose "auto" or "pareto"
#' @param x PCA X axis, default is PC1
#' @param y PCA Y axis, defult is PC2
#' @param exclude exclude some classes of samples
#' @param size dot size
#' @param interactive should interactive figure be plotted? default = TRUE. If you want to save the result
#' in high resolution, use non interative plot.
#' @param scale_group select groups needs to be scaled.
#' @param scale_factor the scale factor, default = 1.
#' @param ... other parameters
#' @importFrom ggplot2 theme_bw autoplot
#' @importFrom graphics abline legend par plot title
#' @importFrom plotly ggplotly
#' @importFrom stats prcomp sd
#' @import ggfortify
#' @export
#'@examples
#' dat <- matrix(runif(2*300), ncol = 2, nrow = 300)
#' Group <- rep_len(LETTERS[1:3], 300)
#' out <- viewPCA(dat, Group = Group)

viewPCA <- function(dat, Group = NULL, centering = T, scaling = "none", x = 1, y = 2,
                 size = 1.5, exclude = NULL, scale_group = NULL, scale_factor = 1,
                 interactive = T, ...) {

  #(1) check input
  Group <- as.factor(Group)
  if(is.null(Group)){stop("Please include group information")}
  if(length(levels(Group)) <= 1){stop("At least two sample groups should be included")}
  if(length(Group) != nrow(dat)){stop("Missing group informaiton detected")}
  if(!(toupper(scaling) %in% c("PARETO", "AUTO", "NONE")) == TRUE) {stop("wrong scaling method")}
  dat <- cbind.data.frame(Group = Group, dat)

  ##(2) exclude groups
  if (is.null(exclude) == FALSE){
    dat <- dat[which(!(dat$Group %in% (exclude))), ]
  }

  ##(3) add scaling factor. for instance, some samples were diluted 10 times, so the intensity should be multiplied by 10
  if (is.null(scale_group) == FALSE){
    dat[which(dat$Group %in% (scale_group)), -1] <- dat[which(dat$Group %in% (scale_group)), -1] * scale_factor
  }

  #(3) apply scaling
  dat2 <- dat[, -1] ## temporally remove Group
  if(toupper(scaling) == "PARETO") {
    dat3 <- apply(dat2, 2, function(x) x/sqrt(sd(x)))
  } else if(toupper(scaling) == "AUTO") {
    dat3 <- scale(dat2, center = FALSE, scale = TRUE)
  } else if(toupper(scaling) == "NONE") {
    dat3 = dat2
  }

  ##(4) plot
  dat3 <- cbind.data.frame(Group = dat$Group, dat3)
  p <- suppressMessages(autoplot(prcomp(dat3[, -1], center = centering), data = dat3, colour = "Group",
                size = size, x = x, y = y, ...) +
    theme_bw())
  p2 <- ggplotly(p)

  ##(5) choose to use interactive or non-interactive plot
  if(isTRUE(interactive) == TRUE) {
    return(p2)
  } else {
    return(p)
  }
}
