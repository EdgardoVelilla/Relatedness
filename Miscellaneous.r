###############################################################################################################################################
#'                                                                                                                                           '#
#' miscelaneous function                                                                                                                     '#
#'                                                                                                                                           '#
#' Author: Edgardo Velilla P.                                                                                                                '#
#' email{edgardo.velilla@cmpc.cl}                                                                                                            '#
#' Created: 05-Mar-2022                                                                                                                      '#
#' Modified: 10-Jul-2022                                                                                                                     '#
#'                                                                                                                                           '#
###############################################################################################################################################


###############################################################################################################################################
#'                                                                                                                                           '#
#'                                              Solving sparse positive definite linear equations                                            '#
#'                                                                                                                                           '#
#'                                                                                                                                           '#
#'                                                        Sparse Cholesky factorization                                                      '#
#'                                                                                                                                           '#
#' If A is sparse and positive definite, it is usually factored as                                                                           '#
#'                                                                                                                                           '#
#'                                                                A = PLL'P'                                                                 '#
#' where                                                                                                                                     '#
#' L is a lower triangular with positive diagonal elements. Thus, L is called the Cholesky factor of A, P a permutation matrix, L' and P'    '#
#' are transpose matrices of L and P.                                                                                                        '#
#'                                                                                                                                           '#
#' Interpretation: we permute the rows and columns of A and Cholesky factor so that                                                          '#
#'                                                                                                                                           '#
#'                                                               P'AP = LL'                                                                  '#
#'                                                                                                                                           '#
#' * choice of P greatly affects the sparsity L, which in turn can greatly affect the efficiency of solving Ax = b                           '#
#'                                                                                                                                           '#
#' * many heristic methods exist for selecting good permutation matrices. For instance, Cholesky() function implemented in the Matrix        '#
#'   package uses the approximate minimum degree algorithm to produce large blocks of zeros in the matrix (Davis 1996, 2006).                '#
#'                                                                                                                                           '#
#'                                                                                                                                           '#
#'                                                Solving sparse positive definite equations                                                 '#
#'                                                                                                                                           '#
#' Solve Ax = b with A sparse positive definite matrix via factorization A= PLL'P' (Boyd and Vandenberghe 2018; Vandenberghe 2021)           '#
#'                                                                                                                                           '#
#'  Algorithm                                                                                                                                '#
#'  1. compute sparse Cholesky factorization A= PLL'P'                                                                                       '#
#'  2. Permute right-hand side (RHS): c:= P'b                                                                                                '#
#'  3. solve Ly= c by forward substitution                                                                                                   '#
#'  4. solve L'z= y by backward substitution                                                                                                 '#
#'  5. permute solution: x:= Pz                                                                                                              '#
#'                                                                                                                                           '#
#'                                                                                                                                           '#
#' References:                                                                                                                               '#
#'                                                                                                                                           '#
#' Timothy A. Davis (2006). Direct Methods for Sparse Linear Systems, SIAM Series “Fundamentals of Algorithms”. p. 21-22.                    '#
#'                                                                                                                                           '#
#' Tim Davis (1996). An approximate minimal degree ordering algorithm, SIAM J. Matrix Analysis and Applications, 17, 4, 886–905.             '#
#'                                                                                                                                           '#
#' Lecture notes from ECE133A (2021) - Applied Numerical Computing, Prof. L. Vandenberghe, UCLA.                                             '#
#' https://www.seas.ucla.edu/~vandenbe/ee133a.html. 7. Linear Equations; 13. Cholesky factorization                                          '#
#'                                                                                                                                           '#
#' S. Boyd and L. Vandenberghe (2018). Introduction to Applied Linear Algebra – Vectors, Matrices, and Least Squares.                        '#
#' p. 207 - 210. https://web.stanford.edu/~boyd/vmls/                                                                                        '#
#'                                                                                                                                           '#
#'                                                                                                                                           '#
###############################################################################################################################################

###############################################################################################################################################
#'                                                                                                                                           '#
#'                                     Some lessons learned about sparse/dense covariance matrices                                           '#
#'                                                                                                                                           '#
#' * The use of sparse matrix methods provided by the "Matrix" package (Bates and Maechler 2019) take advantage of the sparseness of the     '#
#'   covariance matrix. The Matrix package provides S4 classes and methods to objects from these classes provide efficient acces to LAPACK   '#
#'   (Anderson et al. 1999) and accelerated BLAS (Dongarra et al., 1988; Blackford et al., 2002) libraries for dense matrices. The sparse    '#
#'   matrix  methods use CHOLMOD (Davis 2005a), CSparse (Davis, 2005b) and other parts (AMD, COLAMD) of Tim Davis’ “SuiteSparse” collection  '#
#'   of sparse matrix libraries, many of which also use BLAS.                                                                                '#
#'                                                                                                                                           '#
#' * It is almost always a bad idea to compute an inverse; it is more work than needed if we only want to solve one system of equations. It  '#
#'   is usually better to compute and use an appropriate decomposition (eg. functions solve.MME2(), solve.MME(), solve.CHMperm, solve.perm() '#
#'   and solve.m()).                                                                                                                         '#
#'                                                                                                                                           '#
#' * The most natural decomposition to use for covariance matrices is the Cholesky decomposition. Once we have the Cholesky factorization    '#
#'   we do not need the inverse.                                                                                                             '#
#'                                                                                                                                           '#
#' * R functions for computing the Cholesky decomposition of a matrix A are chol() and Cholesky() (from "Matrix" package) and produce the    '#
#'   upper and lower triangular matrices R and L, respectively, such that A = R'R= LL'.                                                      '#
#'                                                                                                                                           '#
#' * Using "forwardsolve" with a sparse matrix "silently" converts our "sparse matrix" to a "dense matrix".                                  '#
#'                                                                                                                                           '#
#' * The Cholesky factorization of the sparse covariance matrix is also sparse, and the sparse matrix method for the "solve" generic         '#
#'   function will take advantage of this and compute the solution z to R'z = y efficiently (e.g., functions solve.MME2(), solve.CHMperm(),  '#
#'   and solve.perm()).                                                                                                                      '#
#'                                                                                                                                           '#
#' References:                                                                                                                               '#
#'                                                                                                                                           '#
#' Anderson, E., Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra, J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney, and D.      '#
#' Sorensen (1999). LAPACK users’ guide, 3rd ed. Society for Industrial and Applied Mathematics, Philadelphia, PA.                           '#
#'                                                                                                                                           '#
#' Douglas Bates and Martin Maechler (2019). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.2-18.                  '#
#' https://CRAN.R-project.org/package=Matrix                                                                                                 '#
#'                                                                                                                                           '#
#' Blackford, L., J. Demmel, J. Dongarra, I. Duff, S. Hammarling, G. Henry, M. Heroux, L. Kaufman, A. Limsdaine, A. Petitet, R. Pozo, K.     '#
#' Remington, and R. C. Whaley. 2002. An updated set of Basic Linear Algebra Subprograms (BLAS). ACM Trans. Math. Softw. 28:135–151.         '#
#'                                                                                                                                           '#
#' Tim Davis. CHOLMOD: sparse supernodal Cholesky factorization and update/downdate. http://www.cise.ufl.edu/research/sparse/cholmod, 2005a. '#
#'                                                                                                                                           '#
#' Tim Davis. CSparse: a concise sparse matrix package. http://www.cise.ufl.edu/research/sparse/CSparse, 2005b.                              '#
#'                                                                                                                                           '#
#' Dongarra, J. J., J. D. Croz, S. Hammarling, and R. J. Hanson. 1988. An extended set of FORTRAN Basic Linear Algebra Subprograms. ACM      '#
#' Trans. Math. Softw. 14:1–17.                                                                                                              '#
#'                                                                                                                                           '#
###############################################################################################################################################

###############################################################################################################################################

# function suggested for solving "sparse" linear equations system of the form CX = b, where C held in the class "dgCMatrix", i.e., general numeric 
# sparse matrix in the compressed sparse column format.
 	 
solve.MME2 <- function(C, b, CheckPD = FALSE){
library(Matrix)
library(data.table)
        
  C <- forceSymmetric(new("dgCMatrix",
		 i=C@i,
		 p=C@p,
		 Dim=C@Dim,
		 x=C@x,
		 Dimnames=list(C@Dimnames[[1L]],
		 C@Dimnames[[2L]])))
  C <- drop0(C, tol = 1e-15,
	    is.Csparse = NA)
  if(CheckPD) C <- check.PD(C)
  R <- Matrix::chol(C, pivot= TRUE)
  Rp <- t(R)
  P <- as(order(attr(R, 'pivot')), 'pMatrix')
  Pp <- t(P)
  c <- Pp%*%b
  y <- solve(Rp, c) # class 'dgeMatrix'
  z0 <- solve(R, y) # class 'dgCMatrix'
  z1 <- as.data.table(as.matrix(P%*%z0))
 z1[]
}

# function suggested for solving "dense" linear equations system of the form CX = b, where C held in the class "dgCMatrix", i.e., general numeric  
# sparse matrix in the compressed sparse column format 	
	
solve.MME <- function(C, b, CheckPD= FALSE){
library(Matrix)
 
  C <- forceSymmetric(new("dgCMatrix",
		 i=C@i,
		 p=C@p,
		 Dim=C@Dim,
		 x=C@x,
		 Dimnames=list(C@Dimnames[[1L]],
		 C@Dimnames[[2L]])))
  C <- drop0(C, tol = 1e-15,
	    is.Csparse = NA)
  if(CheckPD) C <- check.PD(C)
  R <- chol(C, pivot=TRUE)
  P <- as(order(attr(R, 'pivot')),
	        'pMatrix')
  z <- as.vector(
	     backsolve(R, forwardsolve(R,
	     t(P)%*%b,
	     upper.tri = TRUE,
	     transpose = TRUE)))
  sol <- as.vector(P%*%z)
 sol[]
}

# function suggested for solving sparse "symmetric" linear equations system of the form C%*%C^(-1)= I, where C held in the class "dsCMatrix", i.e., 
# symmetric real sparse matrix in the compressed sparse column format. That is, only the upper or the lower triangle of C is stored, and I is 
# an identity matrix. 

solve.CHMperm <- function(C, CheckPD= FALSE){
library(Matrix)
      
	  if(is.null(C@Dimnames[[1L]]))
	  dimnames(C) <- list(seq(1L, nrow(C)),
	                      seq(1L, ncol(C)))
	  m <- dim(C)[1L]
	  b <- Diagonal(m)
      if(CheckPD) C <- check.PD(C)
      CHMp <- Cholesky(C, perm = TRUE)
	  L <- as(CHMp, "Matrix") # L is the Cholesky factor such that C=t(P)%*%L%*%t(L)%*%P <- P'LL'P
	  P <- as(CHMp, "pMatrix") # where P is a permutation matrix
	  z <- solve(t(L), solve(L, t(P)%*%b, # Cinv= P’R’RP=(RP)’(RP) where (RP)' <- t(PR) with  R <- t(L)=L'
		   system="L"), system="Lt")
	  ginv <- P%*%z
      ginv@Dimnames[[1L]] <- C@Dimnames[[1L]]
      ginv@Dimnames[[2L]] <- C@Dimnames[[2L]]
  ginv[]
}

# function suggested for solving "sparse" linear equations system  of the form C%*%C^(-1)= I, where C held in the class "dgCMatrix", i.e., general
# numeric sparse matrix in the compressed sparse column format and I is an identity matrix.

solve.perm <- function(C, CheckPD= FALSE){
library(Matrix)

      if(is.null(C@Dimnames[[1L]]))
	  dimnames(C) <- list(seq(1L, nrow(C)),
	                      seq(1L, ncol(C)))
	  C <- forceSymmetric(new(
	        "dgCMatrix",
		    i=C@i,
		    p=C@p,
		    Dim=C@Dim,
		    x=C@x,
		    Dimnames=list(C@Dimnames[[1L]],
		    C@Dimnames[[2L]])))
	   C <- drop0(C, tol = 1e-15,
			 is.Csparse = NA)
	   if(CheckPD) C <- check.PD(C)
              m <- dim(C)[1L]
	      b <- Diagonal(m)
         CHMp <- Cholesky(C, perm = TRUE)
		 L <- as(CHMp, "Matrix") # L is the Cholesky factor such that C=t(P)%*%L%*%t(L)%*%P <- P'LL'P
		 P <- as(CHMp, "pMatrix") # P is a permuation matrix
		 z <- solve(t(L), solve(L, t(P)%*%b,
		      system="L"), system="Lt") # Cinv= P’R’RP=(RP)’(RP) where (RP)' <- t(PR) with  R <- t(L)=L'
		 Cinv <- P%*%z
         Cinv@Dimnames[[1L]] <- C@Dimnames[[1L]]
         Cinv@Dimnames[[2L]] <- C@Dimnames[[2L]]
  Cinv[]
}

# function suggested for solving "symmetric dense" linear equations system  of the form C%*%C^(-1)= I, where C held in the class "dgCMatrix", i.e., real
# matrix in general storage mode and I is an identity matrix.

solve.m <- function(C, CheckPD= FALSE){
library(Matrix)
      
     if(is.null(C@Dimnames[[1L]]))
	  dimnames(C) <- list(seq(1L, nrow(C)),
	                      seq(1L, ncol(C)))
     m <- dim(C)[1L]
	 b <- Diagonal(m)
	 if(CheckPD) C <- check.PD(C)
	 R <- chol(C, pivot=TRUE)
	 P <- as(order(attr(R, 'pivot')),
	        'pMatrix')
     c <- t(P)%*%b
     z <- backsolve(R, forwardsolve(R, c,
	        upper.tri = TRUE,
	        transpose = TRUE))
	 sol <- P%*%z
	 colnames(sol) <- rownames(sol) <- C@Dimnames[[1L]]
   sol[]
}

###############################################################################################################################################


###############################################################################################################################################
#'                                                                                                                                           '#
#'                                                             General functions                                                             '#
#'                                                                                                                                           '#
###############################################################################################################################################

# function to add a "generation" field to the pedigree
gen.add <- function(pedigree) {
library(data.table)
  pedigree <- pedigree[, c(1L:3L)]
  setnames(pedigree, old = colnames(pedigree),
    new = c('TreeID','mum', 'dad'))
  ped <- copy(pedigree) # get a copy of pedigree to avoid side effect of :=
  ped[, gen := 0L]
        for (i in 0:nrow(ped)) {
          ped[mum %in% TreeID[gen == i] | dad %in% TreeID[gen == i],
		  gen := gen + 1L]
          if (.Last.updated == 0L) break
        }
  ped[]
}

# function to add a "cross" field to the pedigree
makeFam <- function(pedigree){
     ped.cross <- copy(pedigree[, c(1L:3L)])
	 setnames(ped.cross,
                 c('TreeID','mum', 'dad'))
     ped.cross[,
	  cross:= ifelse((mum==0L | dad==0L),
	  NA,
          paste(pmin(mum, dad),
	  pmax(mum, dad),
	  sep = "x"))]
     ped.cross[, cross:= factor(cross,
        levels = unique(cross))]
	ped.cross[]
}

# function to make a matrix of ones in the compressed sparse column format ("dgCMatrix")
ones.matrix <- function(n.rows, n.cols) {
library(Matrix)
    m <- Matrix::Matrix(
	  nrow= n.rows,
	  ncol= n.cols,
          data= 1L,
          sparse= TRUE)
  return(m)
 }

# function to make a matrix of zeros in the compressed sparse column format ("dgCMatrix")
zeros.matrix <- function(n.rows, n.cols) {
library(Matrix)
    m <- Matrix::Matrix(
	  nrow= n.rows,
	  ncol= n.cols,
          data= 0L,
          sparse= TRUE)
  return(m)
 }  

# function to create the relatedness matrix corresponding to non-additive Mendelian Sampling term "non-adjuted" by inbreeding in
# the class "dsCMatrix", i.e., symmetric real sparse matrix in the compressed sparse column format.
Mend.m <- function(pedigree) {
library(Matrix)

if (is.data.frame(pedigree)) {
    pedigree <- as.data.table(pedigree)
 } else if (is.matrix(pedigree)) {
	   pedigree <- as.data.table(pedigree)
 } else if (is.data.table(pedigree)) {
	   pedigree <- as.data.table(pedigree)
 } else {
	   stop("nothing to do...")
 }
  pedi <- pedigree[, c(1L: 3L)]
  setnames(pedi,
    c('TreeID','mum', 'dad'))
  m0 <- length(pedi[mum==0L | dad==0L, TreeID])
  s <- dim(pedi)[1L] - m0
  Im <- .symDiagonal(m0 + s)
  Im[cbind(seq(m0 + 1L, m0 + s),
	    seq(m0 + 1L, m0 + s))] <- 0.75
  Im@Dimnames[[1L]] <- Im@Dimnames[[2L]] <- as.character(pedi[[1]])
 Im[]
}

# function to create the relatedness matrix corresponding to non-additive Mendelian Sampling term "adjuted" by inbreeding ("f") in
# the class "dsCMatrix", i.e., symmetric real sparse matrix in the compressed sparse column format.

Mend.adj <- function(pedigree, f= NULL) {
library(Matrix)
library(pedigree)

if (is.data.frame(pedigree)) {
    pedigree <- as.data.table(pedigree)
 } else if (is.matrix(pedigree)) {
	   pedigree <- as.data.table(pedigree)
 } else if (is.data.table(pedigree)) {
	   pedigree <- as.data.table(pedigree)
 } else {
	   stop("nothing to do...")
 }
  if(is.null(f)) {
   f <- as.vector(calcInbreeding(
	     pedigree[, c(1L:3L)]))
  }
  pedi <- pedigree[, c(1L: 3L)]
  setnames(pedi,
    c('TreeID','mum', 'dad'))
  m0 <- length(pedi[mum==0L | dad==0L, TreeID])
  s <- dim(pedigree)[1L] - m0
  M <- .symDiagonal(m0 + s) # relatedness matrix (Im) for non-additive Mendelian Sampling term
  M[cbind(seq(m0 + 1L, m0 + s),
	  seq(m0 + 1L, m0+s))] <- 0.75*(1-f[seq(m0 + 1L, m0 + s)])
  M@Dimnames[[1L]] <- M@Dimnames[[2L]] <- as.character(pedi[[1L]])
 M[]
}

# A simple (but useful) function to check if sigma is positive-definite, if not, fixed it...!
check.PD <- function(C){
library(Matrix)
    eigen <- eigen(C, symmetric = TRUE)$values
    PD <- ifelse(any(eigen < 0L), 'TRUE', 'FALSE')
    if(PD) {
      F.pd <- as(nearPD(C, corr=FALSE,
	          keepDiag = FALSE,
	          conv.tol = 1e-5,
		  maxit = 300)$mat,
		  'sparseMatrix')
           } else {
               F.pd <- C
           }
  F.pd[]
}

# function to build the family design matrix (Z.fam) "dsparseModelMatrix"

Z.fam <-function(trial){
library(Matrix)
library(MatrixModels)
library(data.table)
 
    if (is.data.frame(trial)) {
       pedigree <- as.data.table(trial)
    } else if (is.matrix(trial)) {
	     pedigree <- as.data.table(trial)
    } else if (is.data.table(trial)) {
	    trial <- trial
    } else {
	  stop("nothing to do...")
    }
    trial[, cross:= factor(cross,
         levels = unique(cross))]
    form <- formula(~ cross -1)
    termsf <- terms(form,
                keep.order = TRUE)
    mf <- model.frame(
	   termsf, data=trial,
	   na.action= na.pass)
    Zfam <- MatrixModels::model.Matrix(
	      form, mf, sparse=TRUE)
    rownames(Zfam) <- c(paste0("TreeID",
	                trial[[1L]]))
  return(Zfam)
}


# function to build the design matrix (Z)
# assumption: all individuals have records

Z.mat <- function(pedigree){
library(Matrix)
library(MatrixModels)
library(data.table)
 if (is.data.frame(pedigree)) {
      pedigree <- as.data.table(pedigree)
 } else if (is.matrix(pedigree)) {
	  pedigree <- as.data.table(pedigree)
 } else if (is.data.table(pedigree)) {
	  pedigree <- pedigree
 } else {
	  stop("nothing to do...")
 }
if(!any(names(pedigree)== 'gen')) pedigree <- gen.add(pedigree)
ped.trial <- pedigree[gen > 0L,]
ped.trial <- ped.trial[, c(1L:3L)]
setnames(ped.trial, c('TreeID','mum', 'dad'))
pedigree <- pedigree[, c(1L:3L)]
ped.trial[, TreeID:= factor(TreeID,
      levels = unique(TreeID))]
form <- formula(~ TreeID -1)
termsf <- terms(form, keep.order = TRUE)
mf <- model.frame(
	    termsf, data=ped.trial,
	    na.action= na.pass)
Zdata <- MatrixModels::model.Matrix(
	       form, mf, sparse=TRUE)
m <- dim(pedigree)[1L]
n <- dim(ped.trial)[1L]
Zbase <- zeros.matrix(
	       n.rows=n,
	       n.cols=m-n)
Z <- cbind(Zbase, Zdata)
Z@Dimnames[[2L]] <- c(paste0("TreeID",
	  pedigree[[1L]]))
Z@Dimnames[[1L]] <- c(paste0("TreeID",
	  ped.trial[[1L]]))
  return(Z)
}


X.mat <- function(trial, effect, DP=TRUE){
library(Matrix)
library(MatrixModels
library(data.table)
  
if (is.data.frame(trial)) {
      trial <- as.data.table(trial)
 } else if (is.matrix(trial)) {
	  trial <- as.data.table(trial)
 } else if (is.data.table(trial)) {
	  trial <- trial
 } else {
	  stop("nothing to do...")
 }
trial.tmp <- copy(trial)
setnames(trial.tmp, 1, "ID")
if(any(names(trial.tmp)== effect)) {
  j <- ifelse(effect== "gen.block", j <- 20L,
	   ifelse(effect== "gen", j <- 5L,
	   ifelse(effect== "Xtype",  j <- 7L,
		 stop("\n effect must be named block or gen or XType \n"))))
}
col <- names(trial.tmp)[j]
trial.tmp[, factor:= trial.tmp[, .SD, .SDcols = col]]
trial.tmp[, factor:= factor(factor,
  levels = unique(factor))]
form <- formula(~ factor -1)
termsf <- terms(form,
  keep.order = TRUE)
mf <- model.frame(termsf,
  data= trial.tmp,
  na.action= na.pass)
X <- MatrixModels::model.Matrix(
  form, mf, sparse=TRUE)
X@Dimnames[[1L]] <- as.character(trial.tmp[, ID])
X@Dimnames[[2L]] <- as.vector(as.matrix(trial.tmp[,
  list(lev=levels(factor)), ]))
  if(DP) {
    r <- dim(X)[2]
    X <- X[, seq(2, r), drop=FALSE]
    return(X)
  } else {
  	return(X)
  }
}
	
X.mat0 <- function(trial){
trial.tmp <- copy(trial)
trial.tmp[, gen:= as.factor(gen)]
form <- formula(~ gen -1)
termsf <- terms(form, keep.order = TRUE)
mf <- model.frame(termsf,
  data= trial.tmp,
  na.action= na.pass)
Zg <- MatrixModels::model.Matrix(
  form, mf, sparse=TRUE)
Zg[is.na(Zg)] <- 0
rownames(Zg) <- c(paste0("gen", trial[[5]]))
r <- dim(Zg)[2]
Zg <- Zg[, seq(2, r), drop=FALSE]  # remove first level to get n-1 level of the factor gen. Otherwise, C matrix is singular...!
}


# function to build the incidence matrix (Z.dom) relating an individual in the pedigree to its family (if it proceeds)
Z.dom <-function(pedigree, init.ghost=NULL){
library(Matrix)
library(MatrixModels)
library(data.table)
 
   if (is.data.frame(pedigree)) {
      pedigree <- as.data.table(pedigree)
   } else if (is.matrix(pedigree)) {
	  pedigree <- as.data.table(pedigree)
   } else if (is.data.table(pedigree)) {
	  pedigree <- as.data.table(pedigree)
   } else {
	  stop("nothing to do...")
   }
    pedi <- pedigree[, c(1L: 3L)]
    setnames(pedi, c('TreeID','mum', 'dad'))
	# copy() to avoid change on child data.table (datafam)
	# affecting by parental data.table (pedigree)
	datafam <- copy(pedi)
    #init.ghost <- datafam[mum==0L & dad==0L & TreeID >= 1e5,
	#  .SD[which.min(TreeID)]][, TreeID]
	if(is.null(init.ghost)) init.ghost <- 1e6
	datafam[dad < init.ghost, cross:= ifelse((mum==0L | dad==0L),
	  NA,
      paste(pmin(mum, dad),
	  pmax(mum, dad),
	  sep = "x"))]
    datafam[, cross:= factor(cross,
      levels = unique(cross))]
    form <- formula(~ cross -1)
    termsf <- terms(form,
	  keep.order = TRUE)
    mf <- model.frame(
	  termsf,
	  data=datafam,
	  na.action= na.pass)
    Zdom <- MatrixModels::model.Matrix(
	  form, mf, sparse=TRUE)
	Zdom <- as(Zdom, "sparseMatrix")
    Zdom@Dimnames[[1L]] <- c(paste0("TreeID",
	  pedigree[[1L]]))
	fam <- as.vector(unique(datafam[, cross]))
	Zdom@Dimnames[[2L]] <- fam[!is.na(fam)]
  return(list(Zdom=Zdom,
              datafam=datafam))
}  

# function to built the inverse of the numerator relationship matrix (Ainv) according to the algorithm given by Quaas (1976) from pedigree package
Ainverse <- function(pedigree)  {
library(data.table)
library(Matrix)
library(pedigree)

    makeAinv(pedigree)
	output <- fread("Ainv.txt")
    setnames(output,
	  old = colnames(output),
    new = c('x','y', 'ai'))
	m <- dim(pedigree)[1L]
	ux <- as.vector(unique(output[, x]))
	idx <- base::order(ux)
	index.x <- data.table(x=ux, idx=idx)
	output <- index.x[output, on=.(x=x)]
	uy <- as.vector(unique(output[, y]))
	idy <- base::order(uy)
	index.y <- data.table(y=uy, idy=idy)
	output <- index.y[output, on=.(y=y)]
	cols <- c('idx', 'idy')
    output[, (cols):= lapply(.SD, as.integer),
	  .SDcols= cols]
	Ainv <- (with(output,
			  Matrix::sparseMatrix(i=idx,
			  j=idy,
			  x=ai,
			  dims=c(m, m),
			  dimnames = list(pedigree[[1L]],
			  pedigree[[1L]]),
			  triangular = FALSE,
			  check = TRUE)))
	Ainv.s <- forceSymmetric(Ainv, uplo="L")
	Ainv.s <- drop0(Ainv.s, tol = 1e-15,
			     is.Csparse = NA)
 return(Ainv.s)
}







