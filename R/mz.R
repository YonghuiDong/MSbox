#' @title Calculate accurate mass-to-charge ratio
#' @description Calculate accurate mass-to-charge ratio (m/z)
#' @author Yonghui Dong
#' @param m chemical formula of an ion, case insensitive
#' @param z charge
#' @param caseSensitive if case sensitive is `FALSE` (default), the elements are seperated by numbers.
#' for instance, Carbon dioxyde can be written as 'c1o2' or any combination of the two elements in lower or
#' upper cases. However, the number of elements should be clearly stated in the chemical formula. if case
#' sensitive is `TRUE`, the elements are seperated by upper case letters. For instance, Carbon dioxyde must
#' be written as 'C1O2' or `CO2`. You don't meed to write the number of the element if it is 1.
#' @importFrom crayon cyan bgMagenta
#' @export
#' @examples
#'  mz('C7h7o1', z = 1)
#'  mz('C7H7O', z = 1, caseSensitive = TRUE)
#'  mz(c('C7H7O4', 'C'), z = -1, caseSensitive = TRUE) # vector input
#'  mz(c('c7h7O4', 'c1'), z = -1)

mz <- function(m, z, caseSensitive = FALSE) {
  options(digits = 12)

  #(1) check input
  if(z == 0) {stop("Warning: charge z = 0 ?")}
  if(z%%1 != 0) {stop("Warning: charge z must be integer")}
  if(is.numeric(z) == FALSE) {stop("Warning: charge z shoule be numeric!")}

  #(2) read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  element.max$Class <- toupper(element.max$Class)
  e <- 0.000548597 # mass of an electron

  #(3) allow the function to a vector of input
  if (length(m) > 1) {
    m = as.list(m)
    acc_mz <- sapply(m, mz, z = z, caseSensitive = caseSensitive)
    names(acc_mz) <- m
    return(acc_mz)
  }

  #(4) main function
  ## If caseSensitive == T, split the mass formula based on upper case letters and add missing 1.
  if(isTRUE(caseSensitive)) {
    m <- gsub("(?<=[A-Z])(?=[A-Z])|(?<=[a-z])(?=[A-Z])|(?<=[A-Za-z])$", "1", m, perl = TRUE)
  }
  ## format input
  ## split the mass formula
  v1 <- strsplit(m, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  # convert the first letter of atom to upper case. case insensitive.
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
    accurate_mz <- (sum(atom_mass*num)-z*e)/abs(z)
    accurate_mz = formatC(accurate_mz, digits = 6, format = "f")
    accurate_mz <- as.numeric(accurate_mz)
    return(accurate_mz)
  }
  else
    message('Wrong chemical formula. Are numbers of some elements missing?')

  #(3) allow the function to a vector of input
  if (length(m) > 1) {
    m = as.list(m)
    acc_mz <- sapply(m, mz, z = z)
    names(acc_mz) <- m
    return(acc_mz)
  }
}
