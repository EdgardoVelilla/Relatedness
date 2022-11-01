#################################################################################################################################################################################################
#'                                                                                                                                                                                             '#
#'                                                                                                                                                                                             '#
#' makeRelatedness function                                                                                                                                                                    '#
#'                                                                                                                                                                                             '#
#' Author: Edgardo Velilla P.                                                                                                                                                                  '#
#' email{edgardo.velilla@cmpc.cl}                                                                                                                                                              '#
#' Created: 31-Ago-2022                                                                                                                                                                        '#
#'                                                                                                                                                                                             '#
#' General description                                                                                                                                                                         '#
#'                                                                                                                                                                                             '#
#' This function generates the five relatedness matrices defined under general inbreeding conditions (Harris 1964; Jacquard 1974). That is, the additive genetic relatedness matrix (A),       '#
#' the dominance relatedness matrix (Dr) in the outbred and homozygous population (Di), the covariance relatedness matrix between additive and dominance (ADi) in the corresponding homozygous '#
#' population, and the relatedness matrix due to sum of squared inbreeding depression (U) based on Jacquard’s nine condensed coefficients of identity for all pairwise combinations between   '#
#' individuals i and j (Jacquard 1974).                                                                                                                                                        '#
#'                                                                                                                                                                                             '#
#' Argument                                                                                                                                                                                    '#
#'                                                                                                                                                                                             '#                                                                                                                                                                                                                                                                                                                                                                                   '#
#' @param  pedigree                                                                                                                                                                            '#
#'                                                                                                                                                                                             '#
#' A data.table/dataframe/matrix where the first three columns correspond to the identifiers for the individual, mother and father, respectively. The row giving the pedigree of an individual '#
#' must appears before any row where that individual appears as parent. Founders use 0 (zero) in the parental columns.                                                                         '#
#'                                                                                                                                                                                             '#
#' @return                                                                                                                                                                                     '#
#'                                                                                                                                                                                             '#
#' A list with the following matrices in sparse matrix format from a pedigree frame.: A, Dr, Di, ADi, and ID                                                                                   '#
#'                                                                                                                                                                                             '#
#' References:                                                                                                                                                                                 '#
#'                                                                                                                                                                                             '#
#' Harris D.L. (1964). Genotypic covariances between inbred relatives. Genetics 50:1319–1348.                                                                                                  '#
#'                                                                                                                                                                                             '#
#' Jacquard A. (1974) The Genetic Structure of Populations. Springer- Verlag, New York.                                                                                                        '#
#'                                                                                                                                                                                             '#
#################################################################################################################################################################################################


makeRelatedness <- function(pedigree) {
						 
library("Matrix")
source("Idcoefs.r")

 if (is.data.frame(pedigree)) {
      pedigree <- as.data.table(pedigree)
 } else if (is.matrix(pedigree)) { 
	  pedigree <- as.data.table(pedigree)
 } else if (!is.data.table(pedigree)) {
	 stop("nothing to do...")
 }

pedigree <- pedigree[, c(1L:3L)]	
setnames(pedigree, c('TreeID','mum', 'dad'))  
cols <- c('TreeID','mum', 'dad') 
pedigree [, (cols):= lapply(.SD, 
  as.integer), .SDcols= cols]	

output <- idcoefs(pedigree)	
  
m <- dim(pedigree)[1]
###  vector indexes (ux, uy) in the arrangement for each position of individual x/y in output data.table ###
ux <- as.vector(unique(output[, x]))
idx <- seq(1:length(ux))
index.x <- data.table(x=ux, idx=idx)
output <- index.x[output, on=.(x=x)]
uy <- as.vector(unique(output[, y]))
idy <- seq(1:length(uy))
index.y <- data.table(y=uy, idy=idy)
output <- index.y[output, on=.(y=y)]
cols <- c('idx', 'idy') 
output[, (cols):= lapply(.SD, as.double), 
  .SDcols= cols]

# bulding the numerator relationship matrix from 
# condensed coefficients of identity (Jacquard 1974)
ksp <- 2*(with(output,                                       
		  Matrix::sparseMatrix(
		  i=idx,                
		  j=idy, x=J1 + 0.5*(J3 + J5 + J7)
		    + 0.25*J8,
		  dims=c(m, m),                              
		  dimnames = list(pedigree[[1]],             
		  pedigree[[1]]),                            
		  triangular = FALSE,                        
		  check = TRUE)))                            
As <- forceSymmetric(new(
        "dgCMatrix",                        
		i=ksp@i,                                     
		p=ksp@p,                                     
		Dim=ksp@Dim,                                 
		x=ksp@x,                                     
		Dimnames=list(ksp@Dimnames[[1]],             
		ksp@Dimnames[[2]]))) 
As <- drop0(As, tol = 1e-15,                   
		    is.Csparse = NA) 		
                                                             
# bulding the dominance relationship matrix in the outbred population (DR)  
# J7 represents the expected fraternity coefficient in terms of condensed 
# identity coefficients between two non-inbred individuals                                                           
Dr <- (with(output,                                          
	   Matrix::sparseMatrix(
	   i=idx,                       
	   j=idy, x=J7,                                      
	   dims=c(m, m),                                     
	   dimnames = list(pedigree[[1]],                    
	   pedigree[[1]]),                                   
	   triangular = FALSE,                               
	   check = TRUE)))                                   
Drs <- forceSymmetric(new(
       "dgCMatrix",                       
	   i=Dr@i, p=Dr@p,                                   
	   Dim=Dr@Dim,                                       
	   x=Dr@x,                                           
	   Dimnames=list(Dr@Dimnames[[1]],                   
	   Dr@Dimnames[[2]])))                           
Drs <- drop0(Drs, tol = 1e-15,                   
		    is.Csparse = NA)          
                                                             
# bulding the dominance relationship matrix 
# in the completely inbred population (DI)                                                             
Di <- (with(output,                                          
	   Matrix::sparseMatrix(
	   i=idx,                       
	   j=idy,                                            
	   x=J1,                                             
	   dims=c(m, m),                                     
           dimnames = list(pedigree[[1]],                            
	   pedigree[[1]]),                                   
	   triangular = FALSE,                               
	   check = TRUE)))                                   
Dis <- forceSymmetric(
        new("dgCMatrix",                       
		i=Di@i,                                      
		p=Di@p,                                      
		Dim=Di@Dim,                                  
		x=Di@x,                                      
		Dimnames=list(Di@Dimnames[[1]],              
		Di@Dimnames[[2]])))                      
Dis <- drop0(Dis, tol = 1e-15,                   
		    is.Csparse = NA) 		            
                                                             
# bulding the covariance relationship matrix
# between additive and dominance effects in 
# the homozygous population (ADI)                                                            
ADi <- (with(output,                                         
		Matrix::sparseMatrix(
		i=idx,                  
		j=idy,                                       
		x= 4*J1 + J3 + J5,                           
		dims=c(m, m),                                
                dimnames = list(pedigree[[1]],                           
		pedigree[[1]]),                              
		triangular = FALSE,                          
		check = TRUE)))                              
ADis <- forceSymmetric(
        new("dgCMatrix",                      
		i=ADi@i,                                     
		p=ADi@p,                                     
		Dim=ADi@Dim,                                 
		x=ADi@x,                                     
		Dimnames=list(ADi@Dimnames[[1]],             
		ADi@Dimnames[[2]])))                     
ADis <- drop0(ADis, tol = 1e-15,                   
		    is.Csparse = NA) 

# calculate the inbreeding coefficients 
# for each individual (fx, fy)				
output[, fx := J1 + J2 + J3 + J4]                             
output[, fy := J1 + J2 + J5 + J6]                                                                                                                                                      

# bulding the relationship matrix 
# for the square of the inbreeding depression (ID)                                                              
ID <- (with(output,                                          
	   Matrix::sparseMatrix(
	   i=idx,                       
	   j=idy,                                            
	   x= J1 + J2 - fx*fy,                               
	   dims=c(m, m),                                     
           dimnames = list(pedigree[[1]],                            
	   pedigree[[1]]),                                   
	   triangular = FALSE,                               
	   check = TRUE)))                                   
IDs <- forceSymmetric(
       new("dgCMatrix",                       
	   i=ID@i,                                           
	   p=ID@p,                                           
	   Dim=ID@Dim,                                       
	   x=ID@x,                                           
           Dimnames=list(ID@Dimnames[[1]],                           
	   ID@Dimnames[[2]])))                           
IDs <- drop0(IDs, tol = 1e-15,                   
		    is.Csparse = NA) 		           
		  
  return(list(A=As, 
              Dr=Drs, 
	      Di=Dis, 
	      ADi=ADis,
	      ID=IDs)) 
}

    
