#' @title View TIC
#' @description View variations of TIC among samples
#' @param x sample ion intensity matrix,row sample, column feature.
#' @param Batch sample batch information. If missing, all the samples will be considered from the same batch
#' @param Seq sample sequence with each batch. If missing, the Seq will be automatically assigned according to sample order
#' @param Group sample group information
#' @param resultBy show the result by Batch or by Group (default).
#' @param Trans How should data be transformed, "LOG2", "LOG10", or NULL transformation?
#' @importFrom reshape2 melt
#' @import ggplot2
#' @return a box plot
#' @export
#' @examples
#' dat <- matrix(runif(100*9), ncol = 100, nrow = 27)
#' myGroup <- rep_len(LETTERS[1:3], 27)
#' myBatch <- rep(1:3, each = 9, times = 1)
#' mySeq <- c(1:27)
#' out <- viewTIC(dat, Group = myGroup, Batch = myBatch, resultBy = "Batch")

viewTIC <- function(x, Seq = NULL, Batch = NULL, Group = NULL, Trans = "none", resultBy = "Group"){

   #(1) check input
  value = NULL
  Group <- as.factor(Group)
  if(is.null(Group)){stop("Please include group information")}
  if(is.null(Seq)){Seq = 1:nrow(x)}
  if(length(Group) != nrow(x)){stop("Missing group informaiton detected")}
  if(!(resultBy %in% c("Group", "Batch"))){stop("Please choose if you want to show the result by Batch or by Group?")}

  ##(2) transform
  x[x == 0] <- 1
  if(toupper(Trans) == "LOG2"){x = log2(x)}
  if(toupper(Trans) == "LOG10"){x = log10(x)}

  #(3) view TIC
  ## order the data
  if(is.null(Batch)) {
    dat2 <- cbind.data.frame(Seq = Seq, Group = Group, x)
    dat2 <- dat2[order(dat2$Seq, dat2$Group), ]
    dat3 <- melt(dat2, id = c("Seq", "Group"))
  } else {
    dat2 <- cbind.data.frame(Seq = Seq, Batch = Batch, Group = Group, x)
    dat2 <- dat2[order(dat2$Seq, dat2$Batch, dat2$Group), ]
    dat3 <- melt(dat2, id = c("Seq", "Batch", "Group"))
  }

  ## to plot
  dat3$Group <- as.factor(dat3$Group)
  dat3$Seq <- as.factor(dat3$Seq)
  if(resultBy == "Group"){
    p <- ggplot2::ggplot(dat3, aes(x = Seq, y = value, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
      theme_bw() +
      xlab("Sample (injection sequence)") +
      ylab("Intensity")
  } else if (resultBy == "Batch"){
    dat3$Batch <- as.factor(dat3$Batch)
    p <- ggplot2::ggplot(dat3, aes(x = Seq, y = value, fill = Batch)) +
      geom_boxplot(alpha = 0.5) +
      facet_grid(. ~ Batch, scales = "free_x", space = "free_x") +
      theme_bw() +
      xlab("Sample (injection sequence)") +
      ylab("Intensity")
    }

  return(p)
}


