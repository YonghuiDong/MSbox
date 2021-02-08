#' @title Performing statistics
#' @description performing statistics, including calculating fold change, p-values and VIP values
#' @param x sample ion intensity matrix, row sample, column feature.
#' @param Group sample group information
#' @return a dataframe with statistical information
#' @export
#' @examples
#' dat <- matrix(runif(2*300), ncol = 2, nrow = 300)
#' rownames(dat) <- 1:dim(dat)[1]
#' myGroup <- rep_len(LETTERS[1:3], 300)
#' ret <- doStat(dat, Group = myGroup)

doStat <- function(x, Group = NULL){
  cat("\nPerforming statistics: \n")
  maxSample <- getMax(x)
  myCV <- getCV(x, Group = Group)
  myFC <- getFC(x, Group = Group)
  myP <- getP(x, Group = Group)
  myOPLSDA <- getOPLSDA(x, Group = Group)
  myStat <- cbind(maxSample, myCV, myFC, myP, myOPLSDA)
  cat("\nStatistical analysis done. \n")
  rownames(myStat) <- colnames(x)
  return(myStat)
}
