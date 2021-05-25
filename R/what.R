#' @title search for m/z in from the idiom metabolomics database
#' @description tentative metabolite identification based on m/z value search
#' @author Yonghui Dong
#' @param myMZ  m/z values
#' @param ppm mass tolerance, default value = 10
#' @param mode ionization mode, either positive '+' or negative '-'
#' @param useDB which database to use, HMDB or KEGG? default is HMDB
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' @examples
#' a = what(133.014, mode = '-', ppm = 10)

what <- function (myMZ, mode = NULL, ppm = 5, useDB = "HMDB") {

  ##(1) input check
  if(is.numeric(myMZ) == FALSE) {stop("warning: mass to charge ratio mz shoule be numeric!")}
  if(mode != "+" & mode !="-") {stop("warning: ion mode invalid. Choose '+' or '-'.")}
  if(!(toupper(useDB) %in% c("HMDB", "KEGG"))) {stop("warning: selected database does not exist")}
  cat("Searching started...\n")

  ##(2) search in database
  Result <- vector("list", length(mz))
  ## select DB
  if(mode == '-' & toupper(useDB) == "HMDB") {DB <- as.data.frame(sysdata$HMDB_neg)}
  if(mode == '+' & toupper(useDB) == "HMDB") {DB <- as.data.frame(sysdata$HMDB_pos)}
  if(mode == '-' & toupper(useDB) == "KEGG") {DB <- as.data.frame(sysdata$KEGG_neg)}
  if(mode == '+' & toupper(useDB) == "KEGG") {DB <- as.data.frame(sysdata$KEGG_pos)}

  pb <- txtProgressBar(min = 0, max = length(myMZ), style = 3)
  for (i in 1:length(myMZ)) {
    # set progress bar
    setTxtProgressBar(pb, i)
    Result[[i]] <- subset(DB, abs(mz-myMZ[i]) * 10^6 /mz <= ppm)
    if(dim(Result[[i]])[1] > 0) {
      ppms <- (Result[[i]]$mz - myMZ[i]) * 10^6/Result[[i]]$mz
      Result[[i]] <- cbind(QuerryID = i, Search = myMZ[i], Result[[i]], ppm = round(ppms, 2))
      row.names(Result[[i]]) <- NULL
    } else{
      Result[[i]] <- NULL
    }
  }
  search_result <- do.call(rbind.data.frame, Result)
  ## rm columns with all NAs, this is specially for KEGG database
  search_result <- search_result[, !apply(is.na(search_result), 2, all)]

  ##(3) check if Result is empty, and return result
  if(nrow(search_result) == 0) {
    message('Not Found, Unknown')
  } else {
    return(search_result)
  }
}
