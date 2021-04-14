#' @title Search in customized database
#' @description search in customized database based on accurate m/z and RT
#' @author Yonghui Dong
#' @param DF input file, should contain at least a column named mz
#' @param DB database, should contain at least a column named mz
#' @param ppm mass tolerance, default 5ppm
#' @param RT retention time tolerance, default 0.2min
#' @param useRT should RT be considered during database search?
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' @examples
#' \dontrun{
#' searchDB(DF, DB)
#' }

searchDB <- function(DF, DB, ppm = 5, RT = 0.2, useRT = FALSE){

  #(1) input check
  ## rename DF and DB
  names(DF)[names(DF) == "mz"] <- "My.mz"
  names(DF)[names(DF) == "rt"] <- "My.RT"
  names(DB)[names(DB) == "mz"] <- "DB.mz"
  names(DB)[names(DB) == "rt"] <- "DB.RT"
  if(!'My.mz' %in% colnames(DF)) {stop("Column 'mz' is not found in your data")}
  if(!'DB.mz' %in% colnames(DB)) {stop("Column 'mz' is not found in your database")}
  cat("Searching started...\n")

  #(2) define function and search
  Result <- vector(mode = "list", length = length(DF$My.mz))
  pb <- txtProgressBar(min = 0, max = length(DF$My.mz), style = 3)

  for (i in 1:length(DF$My.mz)) {
    ## progress bar
    setTxtProgressBar(pb, i)

    if(useRT == TRUE) {
      if(!'My.RT' %in% colnames(DF)) {stop("Column 'rt' is not found in your data")}
      if(!'DB.RT' %in% colnames(DB)) {stop("Column 'rt' is not found in your database")}
      ## filter according to ppm and RT
      Result[[i]] <- subset(DB, abs(DB$DB.mz - DF$My.mz[i]) * 10^6 / DB$DB.mz <= ppm & abs(DB$DB.RT - DF$My.RT[i]) < RT)
      if(dim(Result[[i]])[1] > 0) {
        ppms <- (Result[[i]]$DB.mz - DF$My.mz[i]) * 10^6 / Result[[i]]$DB.mz
        rts <- Result[[i]]$DB.RT - DF$My.RT[i]
        Result[[i]] <- cbind(QuerryID = i, DF[i, ], Result[[i]], ppm = round(ppms, 2), RT_dif = round(rts, 2))
        row.names(Result[[i]]) <- NULL
      } else{
        Result[[i]] <- NULL
      }

    } else {
      Result[[i]] <- subset(DB, abs(DB$DB.mz - DF$My.mz[i]) * 10^6 / DB$DB.mz <= ppm)
      if(dim(Result[[i]])[1] > 0) {
        ppms <- (Result[[i]]$DB.mz - DF$My.mz[i]) * 10^6 / Result[[i]]$DB.mz
        Result[[i]] <- cbind(QuerryID = i, DF[i, ], Result[[i]], ppm = round(ppms, 2))
        row.names(Result[[i]]) <- NULL
      } else{
        Result[[i]] <- NULL
      }
    }
  }

  close(pb)

  #(3) format and return result
  myResult <- do.call(rbind.data.frame, Result)
  return(myResult)
}

