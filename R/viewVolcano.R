#' @title View volcano plot
#' @description View volcano plot.
#' @author Yonghui Dong
#' @param result result from doStat() function
#' @param compare_group which groups you want to compare, i.e. c("WT", "Treat1"), only two groups are allowed
#' @param FC select fold change values, default = 2
#' @param p select p value, default = 0.05
#' @importFrom plotly plot_ly layout add_segments %>%
#' @export
#'@examples
#'\dontrun{
#'viewVolcano(result, compare_group = c("WT", "JA"))
#'}

viewVolcano <- function(result, compare_group, FC = 2, p = 0.05) {

  #(1)check input
  if(!is.numeric(FC)){stop("invalid calss of fold change: not numeric")}
  if(!is.numeric(p)){stop("invalid calss of p value: not numeric")}
  if(p > 1){stop("p value should between 0 and 1")}
  if(FC < 0){stop("FC value cannot be negative value ")}
  ## suppress "no visible binding for global variable" notes in ggplot2
  Compoud_QueryID <- Fold_change <- aes <- my_cor <- p.value <- NULL

  #(2) get index
  F_id <- paste("Fold_", compare_group[1], "_vs_", compare_group[2], sep = "")
  P_id1 <- paste("pValue_", compare_group[1], "_vs_", compare_group[2], sep = "")
  P_id2 <- paste("pValue_", compare_group[2], "_vs_", compare_group[1], sep = "")

  #(3) plot
  ## select which p value colname exist.
  if ((P_id1 %in% colnames(result)) == TRUE) {
      P_id = P_id1
  } else if ((P_id2 %in% colnames(result)) == TRUE) {
      P_id = P_id2
  } else {stop("The selected groups do not exist")}

  ## prepare data frame for plot
  F_iden <- result[, F_id]
  P_iden <- result[, P_id]
  new_iden <- cbind.data.frame(Fold_change = log2(F_iden), p.value = -log10(P_iden))

  ## sometimes there are inf, NAN values in FC and p-value, so i replace these values with max(FC) and 0 for FC and p, respectively
  if(length(new_iden[!is.finite(new_iden$Fold_change), ]$Fold_change) > 0)
    {new_iden[!is.finite(new_iden$Fold_change), ]$Fold_change <- max(abs(new_iden[is.finite(new_iden$Fold_change),]$Fold_change))}
  if(length(new_iden[!is.finite(new_iden$p.value), ]$p.value) > 0)
    {new_iden[!is.finite(new_iden$p.value), ]$p.value <- 0}

    ## add color information
    new_iden$my_cor = "NS"
    con1 <- new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value < -log10(p)),]
    if(dim(con1)[1] > 0)
    {new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value < -log10(p)),]$my_cor =
      paste("|FC|", ">=", FC, "&", "p", ">", p, sep = " ")
    }

    con2 <- new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value >= -log10(p)),]
    if(dim(con2)[1] > 0)
    {
      new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value >= -log10(p)),]$my_cor =
        paste("|FC|", ">=", FC, "&", "p", "<=", p, sep = " ")
    }

    con3 <- new_iden[(abs(new_iden$Fold_change) < log2(FC) & new_iden$p.value >= -log10(p)),]
    if(dim(con3)[1] > 0)
    {
      new_iden[(abs(new_iden$Fold_change) < log2(FC) & new_iden$p.value >= -log10(p)),]$my_cor =
        paste("|FC|", "<", FC, "&", "p", "<=", p, sep = " ")
    }

    ##plot
    new_iden$names <- rownames(result)
    pal <- c("steelblue1", "Tomato", "springgreen3", "Grey")
    myplot <- plot_ly(new_iden[, -4], x = ~ Fold_change, y = ~ p.value, type = 'scatter', mode = 'markers',
                      color = ~ my_cor, colors = pal,
                      hoverinfo = 'text',
                      text = ~ paste('</br> Feature: ', new_iden$names,
                                     '</br> FC: ', formatC(2^Fold_change, format = "e", digits = 2),
                                     '</br> p.value: ', formatC(10^(-p.value), format = "e", digits = 2))
    )

    final = myplot %>%
      layout(xaxis = list(zeroline = F,  title = "Log2(FC)"),
             yaxis = list(zeroline = F,  title = "-Log10(p.value)")) %>%
      add_segments(x = log2(FC), xend = log2(FC),
                   y = min(new_iden$p.value), yend = max(new_iden$p.value),
                   line = list(dash = "dash", color = "black"), showlegend = FALSE) %>%
      add_segments(x = -log2(FC), xend = -log2(FC),
                   y = min(new_iden$p.value), yend = max(new_iden$p.value),
                   line = list(dash = "dash", color = "black"), showlegend = FALSE) %>%
      add_segments(x = min(new_iden$Fold_change), xend = max(new_iden$Fold_change),
                   y = -log10(p), yend = -log10(p),
                   line = list(dash = "dash", color = "black"), showlegend = FALSE)

  #(4) final plots
  final
}
