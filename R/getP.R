#' @title get p-values
#' @description get p-values from Post Hoc analysis
#' @param x sample ion intensity matrix, row sample, column feature.
#' @param Group sample group information
#' @importFrom stats as.formula TukeyHSD formula terms terms.formula t.test
#' @return a data frame
#' @export
#' @examples
#' dat <- matrix(runif(2*300), ncol = 2, nrow = 300)
#' myGroup <- rep_len(LETTERS[1:3], 300)
#' out <- getP(dat, Group = myGroup)

getP <- function(x, Group = NULL){
  cat("\n- Calculating p-values...\n")
  #(1) check input
  x <- as.data.frame(x)
  Group <- as.factor(Group)
  if(is.null(Group)){stop("Please include group information")}
  if(length(levels(Group)) <= 1){stop("At least two sample groups should be included")}
  if(length(Group) != nrow(x)){stop("Missing group informaiton detected")}

  #(2) statistical test
  ##(2.1) T-test
  if(length(levels(Group)) == 2){
    f <- as.factor(Group)
    j <- combn(levels(f), 2)
    out <- as.data.frame(sapply(x, function(x) t.test(x ~ Group, alternative = "two.sided")$p.value))
    colnames(out) <- paste0("AdjPvalue_", j[1,], "_vs_", j[2,])
    return(out)
  } else {
  ##(2.2) ANOVA
    response_names <- names(x)
    form <- as.formula(sprintf("cbind(%s) ~ Group", toString(response_names)))
    fit <- do.call("aov", list(formula = form, data = quote(x)))
    aov_hack <- fit
    aov_hack[c("coefficients", "fitted.values")] <- NULL
    aov_hack[c("contrasts", "xlevels")] <- NULL
    attr(aov_hack$model, "terms") <- NULL
    class(aov_hack) <- c("aov", "lm")
    ## post hoc
    N <- length(response_names)
    result <- vector("list", N)
    for (i in 1:N) {
      ## change response variable in the formula
      aov_hack$call[[2]][[2]] <- as.name(response_names[i])
      ## change residuals
      aov_hack$residuals <- as.data.frame(fit$residuals)[, i]
      ## change effects
      aov_hack$effects <- as.data.frame(fit$effects)[, i]
      ## change "terms" object and attribute
      old_tm <- terms(fit)  ## old "terms" object in the model
      old_tm[[2]] <- as.name(response_names[i])  ## change response name in terms
      new_tm <- terms.formula(formula(old_tm))  ## new "terms" object
      aov_hack$terms <- new_tm  ## replace `aov_hack$terms`
      ## replace data in the model frame
      aov_hack$model[1] <- data.frame(fit$model[[1]][, i])
      names(aov_hack$model)[1] <- response_names[i]
      ## run `TukeyHSD` on `aov_hack`
      result[[i]] <- TukeyHSD(aov_hack)$Group[, 4]
    }
    ## output
    col_num <- as.numeric(summary(result)[1])
    #list_names <- names(result[[1]])
    output <- data.frame(matrix(unlist(result), ncol = col_num, byrow = TRUE))
    #colnames(output) <- list_names
    f <- as.factor(Group)
    j <- combn(levels(f), 2)
    colnames(output) <- paste0("AdjPvalue_", j[1,], "_vs_", j[2,])
    return(output)
    }
  }


