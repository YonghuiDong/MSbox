#' @title Get the sample name which has the max ion intensity
#' @description get the sample name which has the max ion intensity
#' @param x sample ion intensity matrix, row sample, column feature.
#' @return a data frame
#' @export
#' @examples
#' dat <- cbind.data.frame(mz = c(100, 101, 300), mz2 = c(0, 0 , 1), mz3 = c(1, 9, 1))
#' rownames(dat) <- c("A", "B", "C")
#' out <- getMax(dat)

getMax <- function(x){
  cat("\n- Geting sample name with max intensity...\n")
  ## get the sample name which has the max intensity of each m/z
  ## get the max peak intensity/area
  ## It is useful for further manual identification
  tmp <- t(x)
  if(is.null(colnames(tmp))) stop("No sample names found")
  sample_index <- as.matrix(apply(tmp, 1, which.max))
  sample_index <- as.vector(sample_index, mode = "numeric")
  sample_name <- colnames(tmp)[sample_index]
  max_Int <- apply(tmp, 1, max)
  max_Result <- cbind.data.frame(maxSample = sample_name, maxInt = max_Int)
  return(max_Result)
  cat("done");
}


