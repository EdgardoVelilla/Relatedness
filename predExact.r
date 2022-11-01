

predExact <- function(pedigree,
                      data,
		      h2=NULL,
		      p1, 
		      a=NULL, 
		      d=NULL, 
		      factor=NULL) {
library(Matrix)
library(data.table)
library(pedigree)
source("makeRelatedness.r")
source("makeKinship.nonadj.r")
source("Miscellaneous.r")
source("quadratic.r")

quadr <- quadratic(p1, a, d)[['comp']]
  sigma2.ar <- quadr[quad=='ar', comp]
  sigma2.dr <- quadr[quad=='dr', comp]
  sigma2.di <- quadr[quad=='di', comp]
  cov.adi <- quadr[quad=='covADI', comp]
  id <- quadr[quad=='ID', comp]

f <- as.vector(pedigree::calcInbreeding(
	     pedigree[, c(1L:3L)]))	  
mu.f <- mean(f) 
var.f <- var(f)	
if(is.null(h2)) h2 <- 0.35					  	
sigma2.e <- ((1 + mu.f)*sigma2.ar + 2*mu.f*cov.adi
         - h2*(1 + mu.f )*sigma2.ar - h2*2*mu.f*cov.adi)/h2  

y <- as.vector(data[, pheno]) 
n <- dim(data)[1] 
m <- dim(pedigree)[1]
X <- ones.matrix(n.rows= n, n.cols= 1)
p <- dim(X)[2] 
R <- W <- Z <- Matrix::.symDiagonal(n)
Rinv <- solve(sigma2.e*R) 
Wf <- as(W%*%f, "sparseMatrix") 

XpRinvX <- crossprod(X, Rinv)%*%X 
XpRinvWf <- crossprod(X, Rinv)%*%Wf 
XpRinvZ <- crossprod(X, Rinv)%*%Z 
XpRinvW <- crossprod(X, Rinv)%*%W 

fpWpRinvX <- t(XpRinvWf)
ZpRinvX <- t(XpRinvZ)
WpRinvX <- t(XpRinvW)

fpWpRinvWf <- crossprod(Wf, Rinv)%*%Wf 
ZpRinvWf <- crossprod(Z, Rinv)%*%Wf 
WpRinvWf <- crossprod(W, Rinv)%*%Wf 
fpWpRinvZ <- crossprod(Wf, Rinv)%*%Z  
ZpRinvZ <- crossprod(Z, Rinv)%*%Z 

fpWpRinvW <- t(WpRinvWf)
fpWpRinvW <- t(WpRinvWf) 
WpRinvW <- crossprod(W, Rinv)%*%W 

XpRinvy <- crossprod(X, Rinv)%*%y
ZpRinvy <- crossprod(Z, Rinv)%*%y
WpRinvy <- crossprod(W, Rinv)%*%y
fpWpRinvy <- crossprod(Wf, Rinv)%*%y

relatedness <- makeRelatedness(pedigree)
  A <- relatedness[['A']]
  Dr <- relatedness[['Dr']]
  DI <- relatedness[['Di']]
  ADI <- relatedness[['ADi']]
  U <- relatedness[['ID']]

G11 <- sigma2.ar*A
G12 <- cov.adi*ADI
G21 <- t(G12)
G22 <- sigma2.dr*Dr + sigma2.di*DI + id^2*U

G11inv <- solve.m(G11)
Q <- as(G22 - G21%*%G11inv%*%G12, "sparseMatrix")
Qinv <- solve(Q)
T <- G11inv%*%G12

P11 <- G11inv + T%*%Qinv%*%t(T)
P12 <- -T%*%Qinv
P21 <- t(P12)
P22 <- Qinv 

# the Mixed Model 
# y ~ 1mu + Zfid + Za + Zd + e   

# The mixed model equations (MME) to solved are:

# | XpRinvX   (pxp)    XpRinvWf     (px1)   XpRinvZ       (pxm)  XpRinvW        (pxm) | XpRinvy   (px1) | #  p equations
# | fpWpRinvX (1xp)    fpWpRinvWf   (1x1)   fpWpRinvZ     (1xm)  fpWpRinvW      (1xm) | fpWpRinvy (1x1) | #  1 equations
# | ZpRinvX   (mxp)    ZpRinvWf     (mx1)   ZpRinvZ + P11 (mxm)  P12            (mxm) | ZpRinvy   (mx1) | #  m equation
# | WpRinvX   (mxp)    WpRinvWf     (mx1)   P21           (mxm)  WpRinvW + P22  (mxm) | WpRinvy   (mx1) | #  m equations

c1 <- cbind(XpRinvX, XpRinvWf, XpRinvZ, XpRinvW) 
c2 <- cbind(fpWpRinvX, fpWpRinvWf, fpWpRinvZ, fpWpRinvW) 
c3 <- cbind(ZpRinvX, ZpRinvWf, ZpRinvZ + P11, P12) 
c4 <- cbind(WpRinvX, WpRinvWf, P21, WpRinvW + P22) 

# the coefficient matrix
C <- rbind(c1, c2, c3, c4) 
RHS <- as(rbind(XpRinvy, fpWpRinvy, ZpRinvy, WpRinvy), "sparseMatrix")  	
BLUP <- as.data.table(as.matrix(solve(C)%*%RHS))
mu <- BLUP[1]
id <- unlist(BLUP[ 2, 1])
dom.id <- as.data.table(
          as.matrix(Wf*id))
setnames(dom.id,'id')  
BLUP_dom <- BLUP[seq(m+p+2, 2*m+p+1)]
setnames(BLUP_dom,'BLUP_dom')
dom <- cbind(dom.id, BLUP_dom)
dom[, dom.adj:= id + BLUP_dom]
BLUP_add <- BLUP[seq(p+2, m+p+1)]
setnames(BLUP_add,'BLUP_add') 
TreeID <-  pedigree[[1]]
BLUPs <- data.table(cbind(TreeID, BLUP_add, dom), 
           stringsAsFactors = FALSE)	
BLUPs[, mu:= as.double(rep(mu, times=m))]
BLUPs[, y_hat:= mu + BLUP_add + dom.adj] 
BLUPs <- data[BLUPs, on=.(TreeID= TreeID)] 
setcolorder(BLUPs, c(1:3,5:9, 4,10))

##############################################################################################################
#                                                                                                            #
#      Random mating is assumed. Hence, dominance relatedness matrix is obtained from Cockerham (1954)       #
#                                                                                                            #
##############################################################################################################

Drinv <- makeKinship.nonadj(pedigree, Dinv=TRUE)[['Dinv']]
Drinv <- 1/sigma2.dr*Drinv

# The mixed model equations (MME) solved are:

# | XpRinvX   (pxp)    XpRinvWf     (px1)   XpRinvZ          (pxm)  XpRinvW           (pxm)     | XpRinvy   (px1) | #  p equations
# | fpWpRinvX (1xp)    fpWpRinvWf   (1x1)   fpWpRinvZ        (1xm)  fpWpRinvW         (1xm)     | fpWpRinvy (1x1) | #  1 equations
# | ZpRinvX   (mxp)    ZpRinvWf     (mx1)   ZpRinvZ + Ainv   (mxm)  WpRinvW           (mxm)     | ZpRinvy   (mx1) | #  m equation
# | WpRinvX   (mxp)    WpRinvWf     (mx1)   ZpRinvZ          (mxm)  WpRinvW + sigma.Dinv  (mxm) | WpRinvy   (mx1) | #  m equations

e1 <- cbind(XpRinvX, XpRinvWf, XpRinvZ, XpRinvW) 
e2 <- cbind(fpWpRinvX, fpWpRinvWf, fpWpRinvZ, fpWpRinvW) 
e3 <- cbind(ZpRinvX, ZpRinvWf, ZpRinvZ + G11inv, WpRinvW ) 
e4 <- cbind(WpRinvX, WpRinvWf, ZpRinvZ, WpRinvW + Drinv) 

# the coefficient matrix
CC <- rbind(e1, e2, e3, e4)   	
BLUP2 <- as.data.table(as.matrix(solve(CC)%*%RHS))
mu2 <- BLUP2[1]
id2 <- unlist(BLUP2[2, 1])
dom.id2 <- as.data.table(
          as.matrix(Wf*id2))
setnames(dom.id2,'id2')  
BLUP_dom2 <- BLUP2[seq(m+p+2, 2*m+p+1)]
setnames(BLUP_dom2,'BLUP_dom2')
dom2 <- cbind(dom.id2, BLUP_dom2)
dom2[, dom2.adj:= id2 + BLUP_dom2]
BLUP_add2 <- BLUP2[seq(p+2, m+p+1)]
setnames(BLUP_add2,'BLUP_add2') 
TreeID <-  pedigree[[1]]
BLUP2s <- data.table(cbind(TreeID, BLUP_add2, dom2), 
           stringsAsFactors = FALSE)	
BLUP2s[, mu2:= as.double(rep(mu2, times=m))]
BLUP2s[, y_hat2:= mu2 + BLUP_add2 + dom2.adj] 
BLUP.both <- BLUPs[BLUP2s, on=.(TreeID= TreeID)] 
BLUP.both[, f:= f]
setcolorder(BLUP.both, c(1:3,17,4,11,6,13,5,12,7,14,8,15,9,10,16))
idep <- data.table(id,id2) 
  return(list(BLUP= BLUP.both, 
              id= idep,
	      quadr= quadr))  
}
