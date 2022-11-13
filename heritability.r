herit <- function(pedigree, 
                  p1, 
		  a=NULL, 
		  d=NULL, 
		  sigma2.e) {
				  
library(data.table)
library(pedigree)
source("quadratic.r")
				  
if(is.null(a)) a <- 1
if(is.null(d)) d <- 1				  
quadr <- quadratic(p1, a, d)[['comp']]
  sigma2.ar <- quadr[quad=='ar', comp]
  sigma2.dr <- quadr[quad=='dr', comp]
  sigma2.di <- quadr[quad=='di', comp]
  cov.adi <- quadr[quad=='covADI', comp]
  id <- quadr[quad=='ID', comp]

f <- as.vector(pedigree::calcInbreeding(
	     pedigree[, c(1L:3L)]))	  
mu.f <- mean(f) 
var.f <- mu.f*(1 - mu.f)

h2 <- ((1 + mu.f)*sigma2.ar + 2*mu.f*cov.adi)/
     (((1 + mu.f)*sigma2.ar + 2*mu.f*cov.adi + sigma2.e))

H2 <- ((1 + mu.f)*sigma2.ar + 2*mu.f*cov.adi 
      + (1 - mu.f)*sigma2.dr +mu.f*sigma2.di+var.f*id^2)/
      ((1 + mu.f)*sigma2.ar + 2*mu.f*cov.adi 
      + (1 - mu.f)*sigma2.dr +mu.f*sigma2.di+var.f*id^2 
      + sigma2.e)
 return(list(h2= h2, 
	     H2= H2)) 

}
