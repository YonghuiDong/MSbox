#' @title view TIC
#' @description view variations of TIC among samples
#' @param dat sample ion intensity matrix
#' @param Batch sample batch information. If missing, all the samples will be considered from the same batch
#' @param Seq sample sequence with each batch. If missing, the Seq will be automatically assigned according to sample order
#' @param Group sample group information
#' @importFrom reshape2 melt
#' @importFrom ggplot2
#' @return a box plot
#' @export
#' @examples
#' dat <- matrix(runif(100*9), ncol = 100, nrow = 10)
#' myGroup <- rep_len(LETTERS[1:5], 10)
#' mySeq <- c(1:10)
#' out <- viewTIC(dat, Group = myGroup)

viewTIC <- function(x, Group = NULL, Batch = NULL, Seq = NULL, Trans = NULL){

   #(1) check input
  Group <- as.factor(Group)
  if(is.null(Group)){stop("Please include group information")}
  if(is.null(Batch)){Batch = 1}
  if(is.null(Seq)){Seq = 1:nrow(x)}
  if(length(Group) != nrow(x)){stop("Missing group informaiton detected")}

  ##(2) transform
  x[x == 0] <- 1
  if(toupper(Trans) == "LOG2"){x = log2(x)}
  if(toupper(Trans) == "LOG10"){x = log10(x)}

  #(3) view TIC
  ## order the data
  datOrder <- paste("B", Batch, "S", Seq, "G", Group, sep = "")
  dat2 <- cbind.data.frame(Order = datOrder, x)
  dat2 <- dat2[order(dat2$Order), ]
  ## reshape data
  dat3 <- suppressMessages({reshape2::melt(dat2)})
  dat3 <- dat3[, -2]
  colnames(dat3) <- c("Sample", "Intensity")
  ## to plot
  p <- ggplot2::ggplot(dat3, aes(x = Sample, y = Intensity, color = Sample)) +
    geom_boxplot() +
    theme_bw() +
    xlab("Sample (in Batch-Seq-Group order)") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1))

  return(p)
}


