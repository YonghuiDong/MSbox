#' @title Get the sample name which has the max ion intensity
#' @description get the sample name which has the max ion intensity
#' @param x sample ion intensity matrix, row sample, column feature.
#' @return a data frame
#' @export
#' @examples
#' dat <- cbind.data.frame(mz = c(100, 101, 300), mz2 = c(0, 0 , 1), mz3 = c(1, 9, 1))
#' rownames(dat) <- c("A", "B", "C")
#' out <- getSname(dat)

getSname <- function(x){
  cat("\n- Geting sample name with max intensity...\n")
  ## get the sample name which has the max intensity of each m/z.
  ## It is useful for further manual identification
  tmp <- t(x)
  if(is.null(colnames(tmp))) stop("No sample names found")
  sample_index <- as.matrix(apply(tmp, 1, which.max))
  sample_index <- as.vector(sample_index, mode = "numeric")
  sample_name <- colnames(tmp)[sample_index]
  return(sample_name)
  cat("done");
}


