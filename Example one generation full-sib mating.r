rm(list=ls(all=TRUE))
gc()

# installing missing packages 
list.of.packages <- c("Matrix", "data.table", "MatrixModels", "pedigree", "DiagrammeR", "bdsmatrix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

####      Table 2      ####
library(data.table)
source("quadratic.r")
Table.2 <- quadratic(p1=0.3, a=1, d=1)
# breeding values
bv <- Table.2[['BV']]
# dominance deviation
dv <- Table.2[['DV']]
# quadratic components
compon <- Table.2[['comp']]
  sigma2.ar <- compon[quad=='ar', comp]
  sigma2.dr <- compon[quad=='dr', comp]
  sigma2.di <- compon[quad=='di', comp]
  cov.adi <- compon[quad=='covADI', comp]
  id <- compon[quad=='ID', comp]
  
####     Figure 4      ####  
  
source("flowgene.r")  
# Hypothetical pedigree that consider one generation of full sib mating
mum <- c(0,0,2,2,4,4)
dad <- c(0,0,1,1,3,3)
# generate the pedigree for Figure 4
ped_Figure4 <- data.table(TreeID=1:6, mum, dad)
flowgene(ped_Figure4, filename="fsmating", ratio=0.7)
library(DiagrammeR)
# plot Figure 4
grViz("fsmating.dot")

####     Table 3      ####  

source("Idcoefs.r")
Table.3 <- idcoefs(ped_Figure4)   
# make full Table 3                 
Table.3[, ks:= J1 + 0.5*(J3 + J5 + J7) + 0.25*J8]

####     Table 4      #### 

# function to create several relatedness matrices based on coancestry coefficient
source("makeKinship.nonadj.r")
# building the dominance relatedness matrix based on Equation 30 (Cockerham 1954)
D.Cockerham <- makeKinship.nonadj(ped_Figure4)[['D']]
# function to create several relatedness matrices based condensed coefficients of 
# identity (Jacquard 1974)

source("makeRelatedness.r")
relatedness <- makeRelatedness(ped_Figure4)
# dominance relatedness matrix in the outbred population based on 
# condensed coefficients of identity (Jacquard 1974)
DR <- relatedness[['Dr']]

####     Table 5      #### 

# the Numerator Relatedness Matrix (NRM) 
AR <- relatedness[['A']]
# the dominance relatedness matrix in the outbred population 
DR 

####     Table 6      #### 

# the dominance relatedness matrix in the completely inbred population (DI) 
DI <- relatedness[['Di']]
# the covariance relatedness matrix between additive and dominance in the homozygous 
# population  
ADI <- relatedness[['ADi']]

####     Table 7      #### 

# the relatedness matrix due to sum of squared inbreeding depression in the homozygous 
# population 
U <- relatedness[['ID']]
# full dominance relatedness matrix: DR + DI + ID
D.full <- DR + DI + U

####     Table 8      #### 

library(pedigree)
# calculation of the inbreeding coefficient
f <- as.vector(pedigree::calcInbreeding(ped_Figure4))
# calculation inbreeding coefficient's mean	 
mu.f <- mean(f) 
# calculation inbreeding coefficient's variance	
var.f <- mu.f*(1- mu.f)	

source("heritability.r")	
herit <- herit(ped_Figure4, 
              p1=0.3, 
	      a=1, 
	      d=1, 
	      sigma2.e=1.5)
			
h2 <- herit[['h2']]
H2 <- herit[['H2']]  

####     Table 9      #### 
              
# phenotypic records
y <- c(17.0,21.0,18.1,20.3,16.1,15.6)
# data (pedigree plus phenotypic records)
data <- copy(ped_Figure4)
data[, pheno:= y]

Table.9 <- data[, f:= f][, c(1:3,5,4)]

####     Table 10      #### 

# details of sub-matrices are specified in Equation 32 
G11 <- sigma2.ar*AR
G12 <- cov.adi*ADI
G21 <- t(G12)
G22 <- sigma2.dr*DR + sigma2.di*DI + id^2*U 
G.e <- rbind(cbind(G11, G12),cbind(G21, G22))

####     Table 11      #### 

source("Miscellaneous.r")
m <- dim(G11)[1]
Zeros <- zeros.matrix(m, m)
G.a <- rbind(cbind(sigma2.ar*AR, Zeros),
             cbind(Zeros, sigma2.dr*D.Cockerham))
			 
####     Table 12      #### 

source("predExact.r")
 
# predictions for both approach
pred.fsm <- predExact(ped_Figure4, data, h2=h2, p1=0.3, a=1, d=1)
# inbreeding depression
id <- pred.fsm[['id']]

# data & pedigree, inbreeding coefficient, BLUE and BLUP
BLUPs <- pred.fsm[['BLUP']]

# subscript=2 refers to approximate approach (Equation 34)
Table.12 <- BLUPs[, c(1:10)]
Table.13 <- BLUPs[, c(1:3,11:17)]

####     Equation 41      #### 

# generate the pedigree for Figure 6
mum <- c(rep(0,4),2,2,4,6,7,9)
dad <- c(rep(0,4),1,1,3,5,6,8)
ped_Figure6 <- data.table(TreeID=1:10, mum, dad)

# make the numerator relatedness matrix (NRM or A)						  
A <- makeKinship.nonadj(ped_Figure6, 
                        NRM= TRUE)[['A']]

# Generalized decomposition of A matrix
# using Henderson notation A = TDT'
library(bdsmatrix)
tmp <- gchol(as.matrix(A)) 
D <- as(diag(diag(tmp)), "sparseMatrix") # Mendelian relatedness matrix 
T <- as(as.matrix(tmp), "dtCMatrix") # matrix T 

# checking decomposition A=TDT’
T%*%D%*%t(T)



