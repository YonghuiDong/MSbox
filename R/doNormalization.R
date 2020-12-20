#' @title perform normalization
#' @description perform normalization
#' @param x sample ion intensity matrix
#' @param method normalization method: (1) LBME: linear baseline normalization based on mean values;
#' (2) LBMD: linear baseline normalization based on median values; (3) PQN: probabilistic quotient normalization;
#' (4): MEDIAN: median normalization; (5) TIC: TIC normalization
#' @return normalized data matrix
#' @export
#' @examples
#' dat <- matrix(runif(100*9), ncol = 100, nrow = 10)
#' out <- doNormalization(dat, method = "PQN" )

doNormalization <-function(x, method = NULL){
  #(1) check input
  if(is.null(method)) {stop("Please select a normalization method")}
  if (!(toupper(method) %in% c("LBME", "LBMD", "PQN", "MEDIAN", "TIC")))
  {stop("Invalid normalization method")}

  x = t(x)
  #(2)perform normalization
  ##(2.1) linear baseline normalization based on mean values
  if(toupper(method) =="LBME"){
    linear.baseline <- apply(x, 1, median) #compute baseline
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

  ##(2.4) Median
  if(toupper(method) == "MEDIAN"){
    norm_vector <- apply(x, 1, median)
    norm.metabo.data <- x / norm_vector
  }

  ##(2.5) TIC
  if(toupper(method) == "TIC"){
    norm_vector <- rowSums(x)
    norm.metabo.data <- x / norm_vector
  }

  return(t(norm.metabo.data))
}
