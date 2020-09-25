#
### Function to assess the relevance of the 4th corner approach for toxicity and habitat modelling.
#
sequentialAICcblm <- function(y,ZU,k) {
  included <- numeric(0)
  candidates <- 1L:ncol(ZU)
  X <- as.data.frame(ZU)
  cnX <- colnames(X)
  while(TRUE) {
    AICc2 <- rep(NA,ncol(X))
    lm1 <- lm(as.formula(paste("y~",if(length(included)) paste(cnX[included],collapse="+") else "1",sep="")), data=X)
    k1 <- length(lm1$coef) ; AICc1 <- AIC(lm1,k=k) + (2*k1*(k1+1)/(length(y)-k1-1))
    for(i in candidates) {
      lm2 <- lm(as.formula(paste("y~",if(length(included)) paste(paste(cnX[included],collapse="+"),"+",sep="") else "",cnX[i],sep="")), data=X)
      k2 <- length(lm2$coef) ; AICc2[i] <- AIC(lm2,k=k) + (2*k2*(k2+1)/(length(y)-k2-1))
    }
    if(min(AICc2,na.rm=TRUE) < AICc1) {
      included <- c(included,candidates[candidates==which.min(AICc2)])
      candidates <- candidates[candidates!=which.min(AICc2)]
    } else {
      lm1$AICc <- AICc1
      return(lm1)
    }  
  }
}
#
Psquare <- function(observed,predicted) {
  n <- length(observed)
  if(n!=length(predicted)) stop("'observed' and 'predicted' must be of the same length")
  return(1-(sum((observed-predicted)^2)/n)/var(observed))
}
#

