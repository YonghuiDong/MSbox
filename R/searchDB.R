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
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
  Result <- vector(mode = "list", length = length(DF$My.mz))
  pb <- txtProgressBar(min = 0, max = length(DF$My.mz), style = 3)

  for (i in 1:length(DF$My.mz)) {
    ## progress bar
    setTxtProgressBar(pb, i)

    DB.list <- expand.grid.df(i, DF[i, ], DB)
    colnames(DB.list)[1] <- "QueryID"
    if(useRT == TRUE) {
      if(!'My.RT' %in% colnames(DF)) {stop("Column 'rt' is not found in your data")}
      if(!'DB.RT' %in% colnames(DB)) {stop("Column 'rt' is not found in your database")}
      ## filter according to ppm and RT
      cal_ppm <- with(DB.list, (DB.mz - My.mz) * 10^6 / DB.mz)
      cal_ppm <-  round(cal_ppm, digits = 2)
      RT_diff <- with(DB.list, abs(My.RT - DB.RT))
      RT_diff <- round(RT_diff, digits = 2)
      DB.list <- cbind(DB.list, Cal.ppm = cal_ppm, RT.dif = RT_diff)
      Result[[i]] = DB.list[(abs(DB.list$Cal.ppm) <= ppm) & DB.list$RT.dif <= RT, ]
      row.names(Result[[i]]) <- NULL
    } else {
      cal_ppm <- with(DB.list, (DB.mz - My.mz) * 10^6 / DB.mz)
      cal_ppm <-  round(cal_ppm, digits = 2)
      DB.list <- cbind(DB.list, Cal.ppm = cal_ppm)
      Result[[i]] = DB.list[(abs(DB.list$Cal.ppm) <= ppm), ]
      row.names(Result[[i]]) <- NULL
    }
  }

  close(pb)

  #(3) format and return result
  myResult <- do.call(rbind.data.frame, Result)
  return(myResult)
}

