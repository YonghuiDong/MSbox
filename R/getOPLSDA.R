#' @title calculate VIP values
#' @description calculate VIP values using among different samples groups using OPLSDA.
#' @param x sample ion intensity matrix, row sample, column feature.
#' @param Group sample group information
#' @importFrom utils combn
#' @importFrom ropls opls
#' @return a dataframe with vip values
#' @export
#' @examples
#' dat <- matrix(runif(2*300), ncol = 2, nrow = 300)
#' myGroup <- rep_len(LETTERS[1:3], 300)
#' out <- getOPLSDA(dat, Group = myGroup)

getOPLSDA <- function(x, Group = NULL){
  cat("\n- Calculating VIP values...\n")
  #(1) check input
  Group <- as.factor(Group)
  if(is.null(Group)){stop("Please include group information")}
  if(length(levels(Group)) <= 1){stop("At least two sample groups should be included")}
  if(length(Group) != nrow(x)){stop("Missing group informaiton detected")}

  #(2) perform OPLS
  n <- levels(Group)
  dat2 <- cbind.data.frame(x, Group = Group)
  mylist <- combn(n, 2, FUN = function(x) subset(dat2, Group %in% x), simplify = FALSE)
  ## change mylist names
  listnames <- combn(n, 2, simplify = FALSE)
  names(mylist) <- lapply(listnames, function(x) paste("OPLSDA", paste(x, collapse="_vs_"), sep = "_"))
  ## perform OPLS-DA using opls package
  oplsdafun <- function(i) {
    myVIP <- opls(x = subset(i, select = -Group),
                  y = as.character(i$Group),
                  predI = 1,
                  orthoI = 1,
                  permI = 10,
                  crossvalI = min(nrow(i), 7),
                  fig.pdfC = "none",
                  info.txtC = "none")
    round(myVIP@vipVn, 2)
  }
  list_VIP <- lapply(mylist, oplsdafun)
  ## convert result to dataframe and keep the name
  return(data.frame(list_VIP))
}
