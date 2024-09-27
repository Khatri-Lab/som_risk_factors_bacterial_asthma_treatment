##################################
# Code from
# https://rdrr.io/bioc/SAGx/src/R/JT.test.r
##################################

# The Lehmann Nonparametrics: Statistical Methods Based on Ranks p. 233 #
# Assumes that groups are coded in increasing numerical order #
# The suggestions posted on R list by Christopher Andrews (SUNY Buffalo, Department of Biostatistics) are gratefully acknowledged

# new idea strsplit(as.character(tet(data = A)), split = "=")[[2]]

JT.test <- function(data,class, labs = NULL, alternative = c("two-sided","decreasing","increasing"), ties = FALSE){
  # data = disc_run_plot_data$adj.level[disc_run_plot_data$phenotype !="Unspecified"]
  # class = disc_run_plot_data$phenotype[disc_run_plot_data$phenotype !="Unspecified"]
  # labs = severity_levels_adj 
  # alternative = "increasing"
  # data <- Response;class <- Treatment;alternative = "two-sided"
  # decreasing means that the null hypothesis states that the trend is decreasing in higher class 
  alternative <- match.arg(alternative)
  calls <- strsplit(as.character(match.call()), split = "=")
  if(!is.ordered(class)) {
    class <- as.ordered(class)
    cat("class was not an ordered factor.  Redefined to be one.\n")
  }
  if(is(data, "exprSet")) stop("Pls update data to ExpressionSet, see Help")
  if(is(data, "ExpressionSet")){ pDataX <- pData(data);class <- pDataX[,paste(class)];data <- exprs(data)}
  if(is.factor(data)) {data <- as.numeric(data);data <- as.matrix(data);factor.ind <- TRUE} else {factor.ind <- FALSE;data <- as.matrix(data)}
  if(!is.null(dim(data))) n.obs <- ncol(data) else n.obs <- length(data)
  if(!(n.obs == length(class))) {data <- t(data);n.obs <- ncol(data)}
  class.tab <- unique(class);class.tab <- class.tab[order(class.tab)]
  if(min(dim(data))==1){
    sums <- 0
    upper <- length(class.tab)-1
    for(i in 1:upper){
      for(j in seq(i+1,upper+1)){
        x <- t(as.matrix(data[,class == class.tab[i]]));y <- t(as.matrix(data[,class == class.tab[j]]))
        n.x <- ncol(x)
        ranked.x <- rank(c(x,y))
        r2= sum(ranked.x[1:n.x])
        sums <- sums + r2 - n.x * (n.x + 1) / 2  
      }
    }
  } else {
    sums <- 0
    upper <- length(class.tab)-1
    for(i in 1:upper){
      for(j in seq(i+1,upper+1)){
        x <- as.matrix(data[,class == class.tab[i]]);y <- as.matrix(data[,class == class.tab[j]])
        n.x <- ncol(x)
        ranked.x <- t(as.matrix(apply(cbind(x,y),1,rank)))
        r2= rowSums(as.matrix(ranked.x[,1:n.x]))
        sums <- sums + r2 - n.x * (n.x + 1) / 2  
      }
    }
  }
  ni <- table(class)
  EH <- (n.obs^2 - sum(ni^2))/4
  if(ties == FALSE) STDH <-sqrt((n.obs^2*(2*n.obs + 3)-sum(ni^2*(2*ni + 3)))/72 ) else{
    dj <- list(apply(data, 1, table));term1 <- sapply(dj, function(x) sum(x*(x - 1)*(2*x + 5)));
    term2 <- sapply(dj, function(x) sum(x*(x - 1)*(x - 2)));term3 <- sapply(dj, function(x) sum(x*(x - 1))) 
    STDH <-sqrt(
      (n.obs*(n.obs - 1)*(2*n.obs + 5) - sum(ni*(ni - 1)*(2*ni + 5)) - term1)/72 +
        sum(ni*(ni - 1)*(ni - 2)*term2)/(36*n.obs*(n.obs - 1)*(n.obs - 2)) +
        sum(ni*(ni - 1))*term3/(8*n.obs*(n.obs - 1)))
  }
  
  # continuity correction remains
  
  ps <- pnorm((EH - sums)/STDH)
  S1 <- (EH - sums)/EH # 09Feb2008 added, see reference Flandre and O'Quigley
  pvalues <- switch(alternative,
                    "two-sided" = 2*pmin(ps,1-ps), "decreasing" = ps, "increasing" = 1-ps)
  
  # pvalues <- 2*pmin(pnorm((sums-EH)/STDH),1-pnorm((sums-EH)/STDH))
  
  # utres <- t(apply(data, 1, function(x) c(2*min(pnorm((sum.stat(x,class)-EH)/STDH),1-pnorm((sum.stat(x,class)-EH)/STDH)),tapply(x,class,median),cor(rank(x),rank(class) ) )))
  
  medians <- t(apply(data, 1, function(x) c(tapply(x,class,median),cor(rank(x),rank(class) ) )))
  if(factor.ind) medians <- class.tab[medians,-ncol(medians)]
  # utres <- data.frame(pvalues, medians)
  if(is.null(labs)) level <- levels(class) else level <- labs
  alternative <- paste(alternative, paste(level,collapse=switch(alternative,
                                                                two.sided = ", ", decreasing= " > ", increasing = " < ")), sep=": ")
  data.name <- paste(calls[[2]], "by", calls[[3]])
  ifelse(is.null(labs),colnames(medians)[-ncol(medians)] <- class.tab, colnames(medians)[-ncol(medians)] <- labs)
  colnames(medians)[ncol(medians)] <- "rank correlation"
  rownames(medians) <- rownames(data)
  utres <- list(statistic = sums, parameter = NULL, p.value = pvalues, method = "Jonckheere-Terpstra", 
                null.value = NULL, alternative = alternative, medians = medians, S1 = S1, data.name = data.name)
  class(utres) <- c("JT-test","htest")
  names(utres$statistic) <- "J"
  return(utres)
}