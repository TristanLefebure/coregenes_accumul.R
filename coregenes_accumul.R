#8/14/2009, TL
require(vegan)

#Count the number of genes in common (need a matrix as input)
coregenes.old <- function(mat) {
	mat <- as.matrix(mat)
	nvect <- ncol(mat)
	if(nvect == 1) {
		length(which(mat[,1] == 1))
	}
	else {
		incommon <- which(mat[,1] == 1)
		for (i in 2:nvect) {
			incommon <- intersect(which(mat[,i] == 1),incommon)
		}
		length(incommon)
	}
}

coregenes.revemp <- function(mat) {
	mat <- as.matrix(t(mat))
	nvect <- ncol(mat)
	if(nvect == 1) {
		length(which(mat[,1] == 1))
	}
	else {
		incommon <- which(mat[,1] == 1)
		for (i in 2:nvect) {
			incommon <- intersect(which(mat[,i] == 1),incommon)
		}
		length(incommon)
	}
}


#Count the number of genes in common (need a matrix as input).
#The matrix does not have to be a binary matrix anymore
#CAREFULL: taxa in rows, genes in cols
#new option to allow some genomes to miss some genes
coregenes <- function(mat, nmiss=0) {
	mat <- as.matrix(mat)
	#if not binary matrix, make it binary
	if(any(mat > 1)) {
	    mat[mat >1] <- 1
	}
	#if a single col (i.e. there was a single genome samples
	if(ncol(mat) == 1) {
	    sum(mat)
	}
	else {
		
	    sum(colSums(mat) >= (nrow(mat) - nmiss) )
# 	    sum(apply(mat, 2, sum) == nrow(mat))
	}
}


#Function to obtain a matrix of permutation of gene in commom. Build using the function specaccum of vegan
#CAREFULL: taxa in rows, genes in cols
coreaccum <- function (comm, method="random",permutations = 100, nmiss=0, ...) {
	x <- comm
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	perm <- array(dim = c(n, permutations))

	#the function running coregenes
	generator <- function(x,n,ind,nmiss=nmiss) {
		result <- vector(length=n)
# 		rmat <- t(x[ind,])
		rmat <- x[ind,]
		for (j in 1:n) {
# 			message(j)
			result[j] <- coregenes(rmat[1:j,], nmiss=nmiss)
# 			result[j] <- coregenes.revemp(rmat[1:j,])
# 			message(result[j])
			
		}
		result
	}

	#run genetor with the different permutations
	for (i in 1:permutations) {
		message(i)
		perm[, i] <- generator(x,n,sample(n),nmiss=nmiss)
	}
	sites <- 1:n
	specaccu <- apply(perm, 1, mean)
	sdaccu <- apply(perm, 1, sd)

	out <- list(call = match.call(), method = method, sites = sites,
	richness = specaccu, sd = sdaccu, perm = perm)
	class(out) <- "specaccum"
	out
}


#function to calculate the number of new genes found when analysing a new genome
newaccum <-  function (comm, method="random",permutations = 100, ...) {
	x <- comm
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	perm <- array(dim = c(n-1, permutations))

	#if not binary matrix, make it binary
	if(any(x > 1)) {
	    x[x >1] <- 1
	}

	#run genetor with the different permutations
	for (i in 1:permutations) {
		message(i)
		pan.accumul <- rowSums(apply(x[sample(n),], 2, cumsum) > 0)
		for (j in 1:n-1) {
		    perm[j, i] <-  pan.accumul[j+1] -  pan.accumul[j]
		}
	}
	sites <- 1:(n-1)
	specaccu <- apply(perm, 1, mean)
	sdaccu <- apply(perm, 1, sd)

	out <- list(call = match.call(), method = method, sites = sites,
	richness = specaccu, sd = sdaccu, perm = perm)
	class(out) <- "specaccum"
	out
}


#a function to return the bxp stat of an acccumulation object and the mean
bxp.stats <- function(x) {
    x.st <- matrix(ncol=6, nrow=length(x$sites))
    for (i in 1:nrow(x$perm)) {
	bx.st <- boxplot.stats(x$perm[i,])
	ave <- mean(x$perm[i,])
	x.st[i,] <- c(bx.st$stats, ave)
    }
    x.st
}

#a function to add accumul object stats to a plot (i.e. the  one returned by bxp.stats)
add.range <- function(x, eps=0, col=1, segwd=1, lwd=1, lty=1, mean=FALSE, log=FALSE, quantile=c(1,5), type='b', elty=2, xrange=NULL,...) {
    n <- nrow(x)
    if(is.null(xrange)) { xrange <- 1:n }
    
    min = quantile[1]
    max = quantile[2]
#     eps <- 0.1
    m <- x
#     color <- col
    if(log==TRUE) { #transform 0 to 0.5
	m[m == 0] <- 0.5
    }
    if(mean==TRUE) {
	lines(xrange +eps, m[,6], col=col, lwd=lwd, lty=lty)
    } else { lines(xrange +eps, m[,3], col=col, lwd=lwd, lty=lty) }
    
    if(type=='b') {
	segments(xrange +eps, m[,min], 1:n +eps, m[,max], col=col, lwd=segwd, lty=elty)
	segments(xrange -0.2 +eps, m[,min], 1:n +0.2 +eps, m[,min], col=col, lwd=segwd)
	segments(xrange -0.2 +eps, m[,max], 1:n +0.2 +eps, m[,max], col=col, lwd=segwd)
    }
    if(type=='l') {
   	lines(xrange +eps, m[,min], col=col, lwd=segwd, lty=elty)
	lines(xrange +eps, m[,max], col=col, lwd=segwd, lty=elty)
    }
#     if(type=='n') {
#     }
}
