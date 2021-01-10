#' @title molecular mass
#' @description calculate accurate molecular mass
#' @author Yonghui Dong
#' @param F chemical formula, case insensitive
#' @param caseSensitive if case sensitive is `FALSE` (default), the elements are seperated by numbers.
#' for instance, Carbon dioxyde can be written as 'c1o2' or any combination of the two elements in lower or
#' upper cases. However, the number of elements should be clearly stated in the chemical formula. if case
#' sensitive is `TRUE`, the elements are seperated by upper case letters. For instance, Carbon dioxyde must
#' be written as 'C1O2' or `CO2`. You don't meed to write the number of the element if it is 1.
#' @importFrom stats aggregate
#' @importFrom crayon cyan bgMagenta
#' @export
#' @examples
#'  mass('C7h7o1')
#'  mass('C7H7O', caseSensitive = TRUE)
#'  mass(c('C7H7O4', 'C'), caseSensitive = TRUE) # vector input
#'  mass(c('c7h7O4', 'c1'))

mass <- function(F, caseSensitive = FALSE) {

  #(1) read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  element.max$Class <- toupper(element.max$Class)

  #(2) allow the function to a vector of input
  if (length(F) > 1) {
    F = as.list(F)
    acc_mass <- sapply(F, mass, caseSensitive = caseSensitive)
    names(acc_mass) <- F
    return(acc_mass)
  }

  #(3) main function
  ## If caseSensitive == T, split the mass formula based on upper case letters and add missing 1.
  if(isTRUE(caseSensitive)) {
    F <- gsub("([A-Z][a-z]?)(?!\\d)","\\11", F, perl = TRUE)
  }
  ## split the mass formula
  v1 <- strsplit(F, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  ## convert the first letter of atom to upper case. case insensitive.
  atom <- toupper(atom)
  ## check the input
  ###(1) check element
  atom_match <- match(atom, element.max$Class)
  if(any(is.na(atom_match))) {stop("Warning: certain element not found")}
  ## special elements that we need to pay attention
  ## for example, we want to know mass of C(carbon)O(oxygen), if your forget to write as C1O1,
  ## it will be calculated for element Co (cobalt)
  Exlist <- c("Co", "Cs", "Cu", "Ho", "In", "Ni", "Os", "Sn", "Si")
  exmatch <- Exlist[which(toupper(Exlist) %in% atom)]
  if(length(exmatch) > 0)
  {cat(cyan("Attention: are following elements indeed in your formula: "), bgMagenta(exmatch), "\n")}
  ###(2) check element number
  num <- as.numeric(v1[c(FALSE, TRUE)])
  if (length(atom) == length(num)) {
    atom_mass <- element.max$Mass[match(atom, element.max$Class)]
    accurate_mass <- sum(atom_mass*num)
    ## set 6 decimal
    accurate_mass = formatC(accurate_mass, digits = 6, format = "f")
    accurate_mass <- as.numeric(accurate_mass)
    return(accurate_mass)
  }
  else
    message('Wrong chemical formula. Numbers of some elements missing?')
}
