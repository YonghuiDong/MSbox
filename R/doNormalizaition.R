#' @title normalize the data
#' @description perform normalization
#' @param x sample ion intensity matrix
#' @param Method normalization method, "TIC", "median",
#' @importFrom utils combn
#' @return a dataframe
#' @export
#' @examples
#' dat <- matrix(runif(2*300), ncol = 2, nrow = 300)
#' ret <- doNormalization(dat, method = )

doNormalization <- function(x){
  cat ("\n- Performing normalization...\n")
  #(1) check input
  if(is.null(Group)){stop("Please include group information")}
  if(length(levels(Group)) <= 1){stop("At least two sample groups should be included")}
  if(length(myGroup) != nrow(dat)){stop("Missing group informaiton detected")}

  # calculate FC
  i <- split(1:nrow(x), Group)
  mean_int <- sapply(i, function(i){colMeans(x[i, ])})
  x <- t(mean_int)
  j <- combn(levels(Group), 2)
  f_change1 <- x[j[1,],] / x[j[2,],]
  f_change2 <- x[j[2,],] / x[j[1,],]
  ## remove NaN in f_change Matrix
  f_change <- rbind(f_change1, f_change2)
  f_change[is.nan(f_change)] <- 0
  rownames(f_change) <- c(paste("Fold_", j[1,], "_vs_", j[2,], sep = ''),
                          paste("Fold_", j[2,], "_vs_", j[1,], sep = ''))
  ret <- as.data.frame(t(f_change))
  return(ret)
}
