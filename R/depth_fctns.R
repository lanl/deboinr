#' Calculates the normalizing term of a beta distribution.
#'
#' @param n Integer.
#' @param p Integer.
#'
#' @noRd
#' 
binom_norm_term=function(n,p){
  if (n<p){
    quot_of_factorials = 0
  } else{
    quot_of_factorials=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))
  }
  return(quot_of_factorials)
}

#' Calculates classical functional band depth for an ensemble of PDFs.
#' Assumes that each row is an observation and each column is a point on the
#' PDF curve.
#'
#' @param densities_matrix Matrix.
#'
#' @noRd
#' 
fBD2 = function(densities_matrix){
	p              = nrow(densities_matrix)
	n              = ncol(densities_matrix)
	rank_of_matrix = apply(densities_matrix,1,rank)
	min_value_row  = apply(rank_of_matrix,1,min)-1
	max_value_row  = n-apply(rank_of_matrix,1,max)
	return((max_value_row*min_value_row+n-1)/binom_norm_term(n,2))
}

#' Calculates modified functional band depth for an ensemble of PDFs.
#' Assumes that each row is an observation and each column is a point on the
#' PDF curve.
#'
#' @param densities_matrix Matrix.
#'
#' @noRd
#' 
fMBD = function(densities_matrix){
	p              = nrow(densities_matrix)
	n              = ncol(densities_matrix)
	rank_of_matrix = apply(densities_matrix,1,rank)
	min_value_row  = rank_of_matrix-1
	max_value_row  = n-rank_of_matrix
	return((rowSums(max_value_row*min_value_row)/p+n-1)/binom_norm_term(n,2))
}

#' Calculates the depth of every curve in the ensemble.  Returns the corresponding
#' depth, the order of the curves, and the indices of any outliers.
#'
#' @param densities_matrix Matrix.
#' @param x_grid Vector.
#' @param method 'MBD' or 'BD2'
#' @param factor Positive Scalar.  Scaling factor for IQR rule.
#'
#' @noRd
#' 
get_depth_of_curves=function(densities_matrix,
                             x_grid,
                             method,
				                     factor=1.5 ){
	tp   = nrow(densities_matrix)
	n    = ncol(densities_matrix)

	if (method=='BD2'){depth=fBD2(densities_matrix)}
	else if (method=='MBD'){depth=fMBD(densities_matrix)}

	order_of_depths = order(depth,decreasing=TRUE)
	index_of_median = which(depth==max(depth))
	m               = ceiling(n*0.5)
		
	# Murph added 11/1/23:
	# The LQD space seems to really squish up at the ends.  I'll
	# cut off these parts for the outlier detection stuff.
	dens_temp = densities_matrix[15:(nrow(densities_matrix)-15),]
	center    = dens_temp[,order_of_depths[1:m]]
	out       = dens_temp[,order_of_depths[(m+1):n]]
	inf       = apply(center,1,min)
	sup       = apply(center,1,max)
	
	dist      = factor*(sup-inf)
	upper     = sup+dist
	lower     = inf-dist
	outly     = (dens_temp<=lower)+(dens_temp>=upper)
	outcol    = colSums(outly)
	remove    = (outcol>0)
	column    = 1:n
	outpoint  = column[remove==1]
	
	return(list(depth    = depth,
	            outpoint = outpoint,
	            medcurve = index_of_median))
}



