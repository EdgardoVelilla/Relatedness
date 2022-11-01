quadratic <- function(p1,
                      a=NULL, 
					            d=NULL,				  
					            factor=NULL) {
library(data.table)
if(is.null(a)) a <- 1
if(is.null(d)) d <- 1
if(is.null(factor)) factor <- 1
  p <- p1
  q <- 1 - p 
  alpha <- a + d*(q - p) 
  # quadratic components
  sigma2.ar <- 2*p*q*alpha^2*factor
  sigma2.dr <- (2*p*q*d)^2*factor
  sigma2.di <- (4*p*q*(p^3 + q^3)*d^2 -(2*p*q*a)^2)*factor
  cov.adi <- 2*p*q*(p-q)*alpha*d*factor
  ID <- -2*p*q*factor
  ID2 <- ID^2*factor 
  comp=c(sigma2.ar,sigma2.dr,sigma2.di,cov.adi,ID)
  quad=c("ar","dr","di","covADI","ID")
  compon <- data.table(comp, quad)
  # breeding values
  a1 <- q*alpha
  a2 <- -p*alpha
  BV <- c(a1,a2)
  # dominance deviation
  d11 <- -2*q^2*d
  d12 <- 2*p*q*d
  d22 <- -2*p^2*d
  DV <- c(d11,d12,d22)
  return(list(BV= BV, 
		      DV= DV, 
		      comp= compon)) 
}