#' @title perform normalization
#' @description perform normalization
#' @param x sample ion intensity matrix
#' @param method normalization method: (1) LBME: linear baseline normalization based on mean values;
#' (2) LBMD: linear baseline normalization based on median values; (3) PQN: probabilistic quotient normalization;
#' (4) QT: quantile normalization; (5) TIC: total ion current normalization.
#' @return normalized data matrix
#' @importFrom stats median
#' @export
#' @examples
#' dat <- matrix(runif(100*10), ncol = 100, nrow = 10)
#' out <- doNormalization(dat, method = "PQN" )

doNormalization <-function(x, method = NULL){
  #(1) check input
  if(is.null(method)) {stop("Please select a normalization method")}
  if (!(toupper(method) %in% c("LBME", "LBMD", "PQN", "QT", "TIC")))
  {stop("Invalid normalization method")}
  x[x == 0] <- 1 # make sure that data do not contain any zeros
  x = t(x) # feature in row and sample in column

  #(2)perform normalization
  ##(2.1) linear baseline normalization based on mean values
  if(toupper(method) =="LBME"){
    linear.baseline <- apply(x, 1, function(x) median(x, na.rm = T)) #compute baseline
    baseline.mean <- mean(linear.baseline)
    sample.means <- apply(x, 2, mean)
    linear.scaling <- baseline.mean/sample.means
    norm.metabo.data <- t(t(x)*linear.scaling)
  }

  ##(2.2) linear baseline normalization based on median values
  if (toupper(method) =="LBMD"){
    linear.baseline <- apply(x, 1, function(x) median(x, na.rm = T)) #compute baseline
    baseline.median<-median(linear.baseline)
    sample.medians<-apply(x, 2, function(x) median(x, na.rm = T))
    linear.scaling<-baseline.median/sample.medians
    norm.metabo.data <- t(t(x)*linear.scaling)
  }

  #(2.3) PQN
  if(toupper(method) == "PQN"){
    reference <- apply(x, 1, function(x) median(x, na.rm = T))
    quotient <- x/reference
    quotient.median <-apply(quotient, 2, function(x) median(x, na.rm = T))
    norm.metabo.data <-t(t(x)/quotient.median)
  }

  #(2.4) QT
  if(toupper(method) == "QT"){
    df_rank <- apply(x, 2, rank, ties.method = "min")
    df_sorted <- data.frame(apply(x, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    index_to_mean <- function(my_index, my_mean){
      return(my_mean[my_index])
    }
    norm.metabo.data <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
    rownames(norm.metabo.data) <- rownames(x)
  }

  #(2.5) TIC
  if(toupper(method) == "TIC"){
    norm.metabo.data = t(t(x)/colSums(x))
  }

  return(t(norm.metabo.data))
}
