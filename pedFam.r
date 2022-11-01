
# Convert a full pedigree into a family pedigree keeping key links between families across generations

pedFam <- function(pedigree) {
library("data.table")

 if (is.data.frame(pedigree)) {
      pedigree <- as.data.table(pedigree)
 } else if (is.matrix(pedigree)) { 
	  pedigree <- as.data.table(pedigree)
 } else if (is.data.table(pedigree)) {
	  pedigree <- as.data.table(pedigree)
 } else {
	  stop("nothing to do...")  
 }	
 pedigree <- pedigree[, c(1:3)]	
 setnames(pedigree, 	
   c('TreeID','mum', 'dad'))    
 ped <- makeFam(pedigree)[!is.na(cross), ]	  
 setkey(ped, cross)
 ped <- ped[, .SD[1L], by = key(ped)] 
 keep <- ped[, TreeID]
 # function to prune the pedigree but retaining the ancestors of specified subset of trees
   cutPed <- function(pedigree, keep){	
     ped.tmp <- copy(pedigree)
	 id.keep <- keep
     n.id <- length(id.keep) + 1L
       while(length(id.keep)!= n.id){
         n.id <- length(id.keep)
         id.keep <- union(na.omit(c(unlist(ped.tmp[, 
		   .(mum, dad)][match(id.keep, 
		   ped.tmp[, TreeID]), ]))), id.keep)
       }
     ped.tmp <- ped.tmp[sort(match(id.keep, 
	   ped.tmp[, TreeID])), ]
	 ped.tmp <- gen.add(ped.tmp)
	ped.tmp[]
   }
 pedfam0 <- cutPed(pedigree, keep)
 pedfam0[order(gen)]
 pedfam0[, gen:= NULL]
 
 # second round...! (to purge those trees that were retained as a 2nd option, but will not be parents in the next generation) #
 pedfam.sel <- makeFam(pedfam0)  
 pedfam.sel <- pedfam.sel[, .SD[1L:5L], 
   by = cross][!is.na(TreeID) & !is.na(cross), ]	
 mum.check <- match(pedfam.sel[, TreeID], 
   pedfam.sel[, mum], nomatch = NA)
 pedfam.sel[, mum2:= mum.check]
 dad.check <- match(pedfam.sel[, TreeID], 
   pedfam.sel[, dad], nomatch = NA)
 pedfam.sel[, dad2:= dad.check]
 pedfam.sel[, 
   score:= ifelse((!is.na(mum2) | !is.na(dad2)),
   1, 0)]
 setorder(pedfam.sel, cross, -score)
 setkey(pedfam.sel, cross) 
 ped.parents <- copy(pedfam.sel)  
 
 # random tree's selection to represent its own family (1st option) #
 pedfam.sel <- pedfam.sel[, .SD[1L], 
   by = key(pedfam.sel)] 
 setcolorder(pedfam.sel, c("TreeID", "mum", "dad"))  
 
 # additional check if 2nd, 3rd,..., 10th next trees are also parent in next generation #                               
 ped.parents <- ped.parents[score==1, .SD[2L:10L], 
   by = cross][!is.na(TreeID) & !is.na(cross), ]
 setcolorder(ped.parents, c("TreeID", "mum", "dad"))   
 
 pedfam <- rbind(pedfam0[mum==0L & dad==0L, c(1L:3L)], 
                 pedfam.sel[, c(1L:3L)],
				 ped.parents[, c(1L:3L)])			 
 pedfam <- gen.add(pedfam)
 setorder(pedfam, gen)
pedfam[]
  } 
