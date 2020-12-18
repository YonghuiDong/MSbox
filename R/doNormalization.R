

doNormalization <- function(x, method = "Median"){
  #(1) check input
  method <- toupper(method)
  if (!is.element(method, c("MEDIAN", "MEAN", "TIC"))) {stop("Invalid normalization method")}

  #(2) get norm_vector
  if (method == "MEDIAN") {
    norm_vector <- apply(x, 1, median, na.rm = TRUE)
  } else if (method == "MEAN") {
    norm_vector <- rowMeans(x, na.rm = TRUE)
  } else if (method == "TIC") {
    norm_vector <- rowSums(x, na.rm = TRUE)
  }

  #(3) normalize
  datNorm <- x / norm_vector

  return(datNorm)

}


