#################################################################################################################################################################################################
#'                                                                                                                                                                                             '#
#'                                                                                                                                                                                             '#
#' forASReml function                                                                                                                                                                          '#
#'                                                                                                                                                                                             '#
#' Author: Edgardo Velilla P.                                                                                                                                                                  '#
#' email{edgardo.velilla@cmpc.cl}                                                                                                                                                              '#
#' Created: 12-Mar-2021                                                                                                                                                                        '#
#' Modified: 26-Mar-2022                                                                                                                                                                       '#
#' References:                                                                                                                                                                                 '#
#' Butler, D. G., Cullis, B.R., A. R. Gilmour, Gogel, B.G. and Thompson, R. 2017. ASReml-R Reference Manual Version 4. VSN International Ltd, Hemel Hempstead, HP1 1ES, UK.                    '#
#'                                                                                                                                                                                             '#
#################################################################################################################################################################################################


#################################################################################################################################################################################################
#'                                                                                                                                                                                             '#
#' General description:                                                                                                                                                                        '#
#' Convert a sparse matrix held in three-column coordinate form to sparseMatrix to asreml’s giv format (ASReml-R version 4.1.0.126)                                                            '#
#'                                                                                                                                                                                             '#
#' ASReml-R requires a very SPECIFIC format to read a known relatedness or inverse relatedness matrix from an external source. This symmetric inverse matrix (or just a symmetric matrix)      '#
#' can be included as (Butler at al. 2017); a) a sparse inverse covariance matrix held in three column co-ordinate form in row major order. This triplet matrix must have class ginv or have   '#
#' attribute INVERSE set to TRUE. b) a sparse relatedness matrix held in three column coordinate form in row major order. If the attribute INVERSE is not set then FALSE is assumed; a         '#
#' attribute must be set. c) a matrix (or Matrix object) with a dimnames attribute giving the levels of the model term being defined. This may be a relatedness matrix or its inverse; if an   '#
#' inverse, it must have an attribute INVERSE set to TRUE.                                                                                                                                     '#
#' In all cases the matrix object must have an attribute rowNames, a character vector that uniquely identifies each row of the relatedness matrix.                                             '#
#'                                                                                                                                                                                             '#
#' Arguments                                                                                                                                                                                   '#
#'                                                                                                                                                                                             '#
#' @param Ginv                                                                                                                                                                                 '#
#'                                                                                                                                                                                             '#
#' A matrix containing an inverse relatedness matrix (or just a relatedness matrix) in three column coordinate form with a rowNames attribute.                                                 '#
#'                                                                                                                                                                                             '#
#' @param ginv                                                                                                                                                                                 '#
#'                                                                                                                                                                                             '#
#' Logical indicating if the relatedness matrix correspond to an inverse; the default is @ginv=TRUE.                                                                                           '#
#'                                                                                                                                                                                             '#
#' @param rowNames                                                                                                                                                                             '#
#'                                                                                                                                                                                             '#
#' A character vector set as the "rowNames" attribute of the input sparse matrix. The default is the dimnames(Gsinv)[[1] attribute of the leading dimension of Gsinv.                          '#
#'                                                                                                                                                                                             '#
#' @param colnames                                                                                                                                                                             '#
#'                                                                                                                                                                                             '#
#' A character vector that defines column's names for the Ginv matrix. The default vector is c("row", "column", "Gsinv")                                                                       '#
#'                                                                                                                                                                                             '#
#' @param identity                                                                                                                                                                             '#
#'                                                                                                                                                                                             '#
#' Logical indicating if the relatedness matrix correspond to an identity matrix; the default is @identity=FALSE.                                                                              '#
#'                                                                                                                                                                                             '#
#' @return                                                                                                                                                                                     '#
#'                                                                                                                                                                                             '#
#' A three column sparse coordinate form matrix in row major order with attribute "INVERSE" set to TRUE and attribute rowNames.                                                                '#
#'                                                                                                                                                                                             '#
#'                                                                                                                                                                                             '#
#' @example                                                                                                                                                                                    '#
#' Example 12.1 (pag 206) from "Linear Models for the Prediction of Animal Breeding Values" by Raphael A. Mrode, 3rd Edition (2014)                                                            '#
#'                                                                                                                                                                                             '#
#' # pedigre information                                                                                                                                                                       '#
#' sire <- c(0,0,0,0,1,3,6,0,3,3,6,6)                                                                                                                                                          '#
#' dam <- c(0,0,0,0,2,4,5,5,8,8,8,8)                                                                                                                                                           '#
#' ID <- c(seq(1:12)) # individuals                                                                                                                                                            '#
#' pedm9 <- data.table(ID=ID, mum=dam, dad=sire)                                                                                                                                               '#
#' Dinv <- makeKinship.nonadj(pedm9, Dinv=TRUE)[['Dinv']]                                                                                                                                      '#
#' library("Matrix")                                                                                                                                                                           '#
#' Dinv.T <- as(Dinv, "dgTMatrix") # convert sparse matrix to class of sparse matrix stored as "Triplets Form"                                                                                 '#
#' Dinv.as <- forASReml(Dinv.T)                                                                                                                                                                '#
#'                                                                                                                                                                                             '#
#################################################################################################################################################################################################

forASReml <- function(Ginv, 
                      ginv=TRUE, 
                      rowNames=NULL, 
		      colnames=c("row", "column", "Gsinv"),
		      identity=FALSE) {    

library("Matrix")
library(data.table)
source("Miscellaneous.r")
     
 if(ginv) {  
    if(inherits(Ginv, "dgCMatrix")) {
	    Ginv <- as(Ginv, "dgTMatrix") 
	} else if(inherits(Ginv, "dsCMatrix")) {
	    Ginv <- as(Ginv, "dgTMatrix") 
        } else if(inherits(Ginv, "dgRMatrix")) {
            Ginv <- as(Ginv, "dgTMatrix") 
        } else if(inherits(Ginv, "matrix")) {
	    Ginv <- as(as(Ginv, "sparseMatrix"), "dgTMatrix")   
        } else if(inherits(Ginv, "TsparseMatrix")) { 
            Ginv <- Ginv 
        } else {
            stop(substitute(Ginv), " must be a Sparse Matrix in Triplet Form")
	}	
    ginv <- data.table(
	        row=Ginv@i + 1, # at sparse matrices, zero is based...!
		col=Ginv@j + 1, 
		Ginv=Ginv@x) 
    ginv <- ginv[which(ginv[, row] >= ginv[, col]), ] # create a lower triangular matrix
    ginv <- as.matrix(ginv[order(ginv[, row], ginv[, col]), ])
    attr(ginv, "INVERSE") <- TRUE 
    class(ginv) <- c(class(ginv), "ginv")
    if(is.null(rowNames)) {
      rowNames <- as.character(Ginv@Dimnames[[1]])  
      attr(ginv, "rowNames") <- rowNames
    } 
    else {
      attr(ginv, "rowNames") <- as.character(rowNames)
      names(ginv) <- colnames
      if(length(dimnames(Ginv))== 0) {
            warning(substitute(Ginv), "  has no dimnames attribute, ASReml´s call 
	      require rowNames attribute")
            dimnames(Ginv) <- list(seq(1, nrow(Ginv)), 
	                               seq(1, ncol(Ginv)))	
      }	
    } 
	return(ginv)
 } else {
    if(inherits(Ginv, "dsCMatrix")) {
	  if(!identity) Ginv <- check.PD(Ginv)
	    ginv <- solve.perm(Ginv)
	    ginv <- as(ginv, "dgTMatrix")
	    ginv <- data.table(
	      row=ginv@i + 1,
	      col=ginv@j + 1, 
	      Ginv=ginv@x) 
	      ginv <- ginv[which(ginv[, row] >= ginv[, col]), ]   
	    ginv <- as.matrix(ginv[order(ginv[, row], ginv[, col]), ])
	    attr(ginv, "INVERSE") <- TRUE 
            class(ginv) <- c(class(ginv), "ginv")
              if(is.null(rowNames)) {
                rowNames <- as.character(Ginv@Dimnames[[1]])  
                attr(ginv, "rowNames") <- rowNames
	      } 
	      else {
		attr(ginv, "rowNames") <- as.character(rowNames)
	        names(ginv) <- colnames
              }
	      if(length(dimnames(Ginv))== 0) {
                 warning(substitute(Ginv), "  has no dimnames attribute, ASReml´s call 
	        require rowNames attribute")
                dimnames(Ginv) <- list(seq(1, nrow(Ginv)), 
	                               seq(1, ncol(Ginv)))	
              }	
	    return(ginv)
        } else if(inherits(Ginv, "dgTMatrix")) {
	    Ginv <- as(Ginv, "dgCMatrix") 
	} else if(inherits(Ginv, "dsCMatrix")) {
	    Ginv <- as(Ginv, "dgCMatrix") 
	} else if(inherits(Ginv, "dgRMatrix")) {
            Ginv <- as(Ginv, "dgCMatrix") 
        } else if(inherits(Ginv, "matrix")) {
	    Ginv <- as(as(Ginv, "sparseMatrix"), "dgCMatrix")   
        } else if(inherits(Ginv, "dgCMatrix")) { 
            Ginv <- Ginv 
        } else {
           stop(substitute(Ginv), "  must be a compressed column 
	      format from class dgCMatrix Sparse Matrix")
	}		  
	if(!identity) Ginv <- check.PD(Ginv)
	ginv <- solve.perm(Ginv)
	ginv <- as(ginv, "dgTMatrix")
	ginv <- data.table(
	    row=ginv@i + 1,
	    col=ginv@j + 1, 
	    Ginv=ginv@x) 
        ginv <- ginv[which(ginv[, row] >= ginv[, col]), ] 
	ginv <- as.matrix(ginv[order(ginv[, row], ginv[, col]), ])
	attr(ginv, "INVERSE") <- TRUE 
        class(ginv) <- c(class(ginv), "ginv")
        if(is.null(rowNames)) {
          rowNames <- as.character(Ginv@Dimnames[[1]])  
          attr(ginv, "rowNames") <- rowNames
	} 
	else {
	     attr(ginv, "rowNames") <- as.character(rowNames)
	     names(ginv) <- colnames
	  if(length(dimnames(Ginv))== 0) {
             warning(substitute(Ginv), "  has no dimnames attribute, ASReml´s call 
	      require rowNames attribute")
             dimnames(Ginv) <- list(seq(1, nrow(Ginv)), 
	                            seq(1, ncol(Ginv)))	
      }
    }   
   return(ginv)
  }   
 
}









