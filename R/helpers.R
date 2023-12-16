# Helper functions for the DeBoinR main method.
## Author: Alexander C. Murph and Justin D. Strait
## Date: October 2023
library(parallel)
library(KernSmooth)
library(pracma)

#' Calculate the pairwise hellinger distances for all matrices.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_hellinger_dist_mat = function(x_grid, 
                                  densities_matrix,
                                  num_cores, 
                                  center_PDFs){
  n           = nrow(densities_matrix)
  p           = ncol(densities_matrix)
  pdfs_matrix = matrix(0,n,p)
  dist_matrix = matrix(0,n,n)
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Make the hellinger distance calculation into a function for parallelization.
  f_hellinger = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          pdf_1       = densities_matrix[i,]
          pdf_2       = densities_matrix[j,]
          diff_vector = sqrt(pdf_1) - sqrt(pdf_2)
          dist_value  = sqrt(sum(diff_vector**2)) / sqrt(2)
          flag        = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }
  
  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_hellinger, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_hellinger)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }
  
  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}

#' Calculate the median and distances from the median for every function according to
#' hellinger distance.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_hellinger_dist_from_median = function(x_grid, 
                                  densities_matrix,
                                  num_cores,
                                  center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,1,n)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix  = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Get the cross-sectional median:
  cross_median = apply(densities_matrix, 2, median)
  for(par_idx in 1:n){
    pdf_1                  = densities_matrix[par_idx,]
    pdf_2                  = cross_median
    dist_value             = sqrt(sum(( pdf_1 - pdf_2)^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  # The overall median will be the one that minimizes the L2 distance
  # from the cross functional median:
  cross_median = densities_matrix[min(which.min(dist_matrix[1,])),]
  dist_matrix  = matrix(0,1,n)
  for(par_idx in 1:n){
    pdf_1       = densities_matrix[par_idx,]
    pdf_2       = cross_median
    diff_vector = sqrt(pdf_1) - sqrt(pdf_2)
    dist_value  = sqrt(sum(diff_vector**2)) / sqrt(2)
    dist_matrix[1,par_idx] = dist_value
  }

  return(list(dist_matrix = dist_matrix, median_curve = orig_densities_matrix[min(which.min(dist_matrix[1,])),], median_in_set = TRUE))
}

#' Calculate the pairwise L2 distance on the normalized LQD space for all pdfs.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_nLQD_dist_mat = function(x_grid, 
                             densities_matrix,
                             num_cores,
                             center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  nlqds_matrix          = matrix(0,n,length(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Start by transforming the densities_matrix into the nLQD space:
  f_par = function(x){pdf_to_nlqd(alpha_mix(densities_matrix[x,], 0.1), x_grid)}
  if(num_cores>1){
    nlqds_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    nlqds_list = lapply(1:n, f_par)
  }
  for(i_idx in 1:n){
    nlqds_matrix[i_idx,] = nlqds_list[[i_idx]]
  }
  
  # Make the nLQD distance calculation into a function for parallelization.
  f_nLQD = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          nLQD_1     = nlqds_matrix[i,]
          nLQD_2     = nlqds_matrix[j,]
          dist_value = sqrt(sum((nLQD_1 - nLQD_2)^2))
          flag       = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }
  
  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_nLQD, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_nLQD)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }
  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}

#' Calculate the median and distances from the median for every function according to
#' L2 distance on the LQD space.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_nLQD_dist_from_median = function(x_grid, 
                                     densities_matrix,
                                     num_cores,
                                     center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,1,n)
  nlqds_matrix          = matrix(0,n,length(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  ######################
  #### So this is a weird situation in the case where we center the PDFs.
  #### I am going to recalculate the median twice in this case: once for
  #### the normalized PDFs and once for the non-normalized.
  
  # Start by transforming the densities_matrix into the nLQD space:
  f_par = function(x){pdf_to_lqd(alpha_mix(densities_matrix[x,], 0.1), x_grid)}
  if(num_cores>1){
    nlqds_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    nlqds_list = lapply(1:n, f_par)
  }
  for(i_idx in 1:n){
    nlqds_matrix[i_idx,] = nlqds_list[[i_idx]]
  }
  
  # Get the cross-sectional median:
  cross_median = apply(nlqds_matrix, 2, median)
  for(par_idx in 1:n){
    nLQD_1                 = nlqds_matrix[par_idx,]
    nLQD_2                 = cross_median
    dist_value             = sqrt(sum((nLQD_1 - nLQD_2)^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  # The cross-section median will not map back properly to the space of PDFs,
  # since it was calculated on a transformed version of the centered PDF space.
  # I am doing the distance calculations on this sort of doubly-transformed space,
  # but the final cross-section median I give will be on the non-normalized PDFs.
  if(center_PDFs){
    f_par = function(x){pdf_to_lqd(alpha_mix(orig_densities_matrix[x,], 0.1), x_grid)}
    if(num_cores>1){
      nlqds_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      nlqds_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      nlqds_matrix[i_idx,] = nlqds_list[[i_idx]]
    }
    
    # Get the cross-sectional median that will actually map back to the
    # space of unnormalized PDFs:
    cross_median = apply(nlqds_matrix, 2, median)
  }

  cross_median = recover_pdf(cdf_to_pdf(lqd_to_cdf(cross_median, x_grid), x_grid), 0.1, x_grid)
  
  return(list(dist_matrix = dist_matrix, median_curve = cross_median, median_in_set = FALSE))
}


#' Calculate the pairwise fisher rao distance for all pdfs.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_fisher_rao_dist_mat = function(x_grid,
                                  densities_matrix,
                                  num_cores,
                                  center_PDFs){
  #g1, g2 is just arccos( \int_0^1 sqrt( g1(t) ) * sqrt( g2(t) ) dt )
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  sq_dens_matrix        = sqrt(densities_matrix)
  d_sq_pdf              = mean(diff(x_grid))
  dist_matrix           = matrix(0,n,n)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }

  # Make the hellinger distance calculation into a function for parallelization.
  f_fish_rao_dist = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          sqpdf_1     = sq_dens_matrix[i,]
          sqpdf_2     = sq_dens_matrix[j,]
          dist_value  = acos(d_sq_pdf * sum( sqpdf_1 * sqpdf_2 ))
          flag        = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }
  
  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_fish_rao_dist, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_fish_rao_dist)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }
  
  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}

#' Calculate the median and distances from the median for every function according to
#' fisher-rao distance
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_fisher_rao_dist_from_median = function(x_grid,
                                   densities_matrix,
                                   num_cores,
                                   center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  sq_dens_matrix        = sqrt(densities_matrix)
  d_sq_pdf              = mean(diff(x_grid))
  dist_matrix           = matrix(0,1,n)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Get the cross-sectional median:
  cross_median = apply(densities_matrix, 2, median)
  for(par_idx in 1:n){
    pdf_1                  = densities_matrix[par_idx,]
    pdf_2                  = cross_median
    dist_value             = sqrt(sum(( pdf_1 - pdf_2)^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  # The overall median will be the one that minimizes the L2 distance
  # from the cross functional median:
  cross_median = sq_dens_matrix[min(which.min(dist_matrix[1,])),]
  dist_matrix  = matrix(0,1,n)
  for(par_idx in 1:n){
    sqpdf_1     = sq_dens_matrix[par_idx,]
    sqpdf_2     = cross_median
    dist_value  = acos(d_sq_pdf * sum( sqpdf_1 * sqpdf_2 ))
    
    dist_matrix[1,par_idx] = dist_value
  }
  
  return(list(dist_matrix = dist_matrix, median_curve = orig_densities_matrix[min(which.min(dist_matrix[1,])),], median_in_set = TRUE))
}



#' Calculate the median and distances from the median for every function according to
#' wasserstein distance.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_wasserstein_dist_mat = function(x_grid, 
                                    densities_matrix,
                                    num_cores,
                                    center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  cdfs_matrix           = matrix(0,n,length(x_grid))
  d_cdf                 = mean(diff(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Start by transforming the densities_matrix into the cdf space:
  f_par = function(x){pdf_to_cdf(densities_matrix[x,], x_grid)}
  if(num_cores>1){
    cdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    cdfs_list = lapply(1:n, f_par)
  }
  
  for(i_idx in 1:n){
    cdfs_matrix[i_idx,] = cdfs_list[[i_idx]]
  }


  f_wass_dist = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          cdf_1     = cdfs_matrix[i,]
          cdf_2     = cdfs_matrix[j,]
          dist_value  = d_cdf * sum(abs(cdf_1 - cdf_2))
          flag        = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }
  
  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_wass_dist, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_wass_dist)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }
  
  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}

#' Calculate pairwise wasserstein distance between all pdfs.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_wasserstein_dist_from_median = function(x_grid, 
                                    densities_matrix,
                                    num_cores,
                                    center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,1,n)
  cdfs_matrix           = matrix(0,n,length(x_grid))
  d_cdf                 = mean(diff(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Start by transforming the densities_matrix into the cdf space:
  f_par = function(x){pdf_to_cdf(densities_matrix[x,], x_grid)}
  if(num_cores>1){
    cdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    cdfs_list = lapply(1:n, f_par)
  }
  
  for(i_idx in 1:n){
    cdfs_matrix[i_idx,] = cdfs_list[[i_idx]]
  }
  
  # Get the cross-sectional median:
  cross_median = apply(densities_matrix, 2, median)
  for(par_idx in 1:n){
    pdf_1                  = densities_matrix[par_idx,]
    pdf_2                  = cross_median
    dist_value             = sqrt(sum(( pdf_1 - pdf_2)^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  # The overall median will be the one that minimizes the L2 distance
  # from the cross functional median:
  cross_median = cdfs_matrix[min(which.min(dist_matrix[1,])),]
  dist_matrix  = matrix(0,1,n)
  for(par_idx in 1:n){
    cdf_1       = cdfs_matrix[par_idx,]
    cdf_2       = cross_median
    dist_value  = d_cdf * sum(abs(cdf_1 - cdf_2))
    
    dist_matrix[1,par_idx] = dist_value
  }
  
  return(list(dist_matrix = dist_matrix, median_curve = orig_densities_matrix[min(which.min(dist_matrix[1,])),], median_in_set = TRUE))
}


#' Calculate the pairwise total variation distance between all PDFs.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_TV_dist_dist_mat = function(x_grid, 
                                densities_matrix,
                                num_cores,
                                center_PDFs){
  n                     = nrow(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  p                     = ncol(densities_matrix)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Make the hellinger distance calculation into a function for parallelization.
  f_tv_dist = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          PDF_1       = densities_matrix[i,]
          PDF_2       = densities_matrix[j,]
          dist_value  = 0.5 * sum(abs(PDF_1 - PDF_2))
          flag        = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }
  
  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_tv_dist, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_tv_dist)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }

  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}

#' Calculate the median and distances from the median for every function according to
#' total variation distance
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_TV_dist_dist_from_median = function(x_grid, 
                                densities_matrix,
                                num_cores,
                                center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  
  # Get the cross-sectional median:
  cross_median = apply(densities_matrix, 2, median)
  for(par_idx in 1:n){
    pdf_1                  = densities_matrix[par_idx,]
    pdf_2                  = cross_median
    dist_value             = sqrt(sum(( pdf_1 - pdf_2)^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  # The overall median will be the one that minimizes the L2 distance
  # from the cross functional median:
  cross_median = densities_matrix[min(which.min(dist_matrix[1,])),]
  dist_matrix  = matrix(0,1,n)
  for(par_idx in 1:n){
    PDF_1                  = densities_matrix[par_idx,]
    PDF_2                  = cross_median
    dist_value             = 0.5 * sum(abs(PDF_1 - PDF_2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  return(list(dist_matrix = dist_matrix, median_curve = orig_densities_matrix[min(which.min(dist_matrix[1,])),], median_in_set = TRUE))
}

#' Calculate the Bayes distance between all centered log ratio-ed pdfs.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_CLR_dist_mat = function(x_grid, 
                            densities_matrix,
                            num_cores,
                            center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  clr_matrix            = matrix(0,n,length(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Start by transforming the densities_matrix into the nLQD space:
  f_par = function(x){pdf_to_clr(alpha_mix(densities_matrix[x,], 0.1), x_grid)}

  if(num_cores>1){
    clr_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    clr_list = lapply(1:n, f_par)
  }
  for(i_idx in 1:n){
    clr_matrix[i_idx,] = clr_list[[i_idx]]
  }
  
  # Make the nLQD distance calculation into a function for parallelization.
  f_clr = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          # We will just do L2 distance from here:
          clr_1      = clr_matrix[i,]
          clr_2      = clr_matrix[j,]
          dist_value = sqrt(sum( (clr_1 - clr_2)^2 ) )
          flag       = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }

  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_clr, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_clr)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }
  
  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}

#' Determines the median and calculates distance from the median for all functions
#' according to Bayes distance.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_CLR_dist_from_median = function(x_grid, 
                            densities_matrix,
                            num_cores,
                            center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,1,n)
  clr_matrix            = matrix(0,n,length(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Start by transforming the densities_matrix into the CLR space:
  f_par = function(x){pdf_to_clr(alpha_mix(densities_matrix[x,], 0.1), x_grid)}
  
  if(num_cores>1){
    clr_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    clr_list = lapply(1:n, f_par)
  }
  for(i_idx in 1:n){
    clr_matrix[i_idx,] = clr_list[[i_idx]]
  }
  
  # Get the cross-sectional median:
  cross_median = apply(clr_matrix, 2, median)
  for(par_idx in 1:n){
    clr_1                  = clr_matrix[par_idx,]
    clr_2                  = cross_median
    dist_value             = sqrt(sum(( clr_1 - clr_2 )^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  # The overall median will be the one that minimizes the L2 distance
  # from the cross functional median:
  cross_median = clr_matrix[min(which.min(dist_matrix[1,])),]
  dist_matrix  = matrix(0,1,n)
  for(par_idx in 1:n){
    clr_1                  = clr_matrix[par_idx,]
    clr_2                  = cross_median
    dist_value             = sqrt(sum(( clr_1 - clr_2 )^2))
    dist_matrix[1,par_idx] = dist_value
  }
  
  return(list(dist_matrix = dist_matrix, median_curve = orig_densities_matrix[min(which.min(dist_matrix[1,])),], median_in_set = TRUE))
}

#' Transform PDFs into nLQDs, then feed these into the get_depth_of_curves function from the
#' fda package (using BD2).
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param k Float.  The factor by which to expand the IQR when calculating outliers.
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_BD_fboxplot = function(x_grid, 
                           densities_matrix,
                           k,
                           num_cores,
                           center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  nlqds_matrix          = matrix(0,n,length(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }

  # Start by transforming the densities_matrix into the nLQD space:
  f_par = function(x){pdf_to_lqd(alpha_mix(densities_matrix[x,], 0.1), x_grid)}
  if(num_cores>1){
    nlqds_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    nlqds_list = lapply(1:n, f_par)
  }
  for(i_idx in 1:n){
    nlqds_matrix[i_idx,] = nlqds_list[[i_idx]]
  }
  t_nlqds_matrix = t(nlqds_matrix)
  depth_output = get_depth_of_curves(t_nlqds_matrix, x_grid = x_grid, method = "BD2",
                        factor = k)
  return( list(depth_output=depth_output) )
}

#' Transform PDFs into nLQDs, then feed these into the get_depth_of_curves function from the
#' fda package (using MBD).
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param k Float.  The factor by which to expand the IQR when calculating outliers.
#' @param num_cores Integer. Number of cores used if parallelizing.
#'
#' @noRd
#' 
get_MBD_fboxplot = function(x_grid, 
                            densities_matrix,
                            k,
                            num_cores,
                            center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  nlqds_matrix          = matrix(0,n,length(x_grid))
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Start by transforming the densities_matrix into the nLQD space:
  f_par = function(x){pdf_to_lqd(alpha_mix(densities_matrix[x,], 0.1), x_grid)}
  if(num_cores>1){
    nlqds_list = mclapply(1:n, f_par, mc.cores = num_cores)
  } else {
    nlqds_list = lapply(1:n, f_par)
  }
  for(i_idx in 1:n){
    nlqds_matrix[i_idx,] = nlqds_list[[i_idx]]
  }
  
  t_nlqds_matrix = t(nlqds_matrix)
  depth_output   = get_depth_of_curves(t_nlqds_matrix, x_grid = x_grid, method = "MBD",
                                       factor = k)
  
  return( list(depth_output=depth_output) )
}

#' Calculate pairwise distances for pdfs with a user-defined function.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param user_dist_function R Function.  User-defined.
#' @param num_cores Integer.  Number of cores to use if parallelizing.
#'
#' @noRd
#' 
get_user_defined_dist_mat = function(x_grid, 
                                    densities_matrix,
                                    user_dist_function,
                                    num_cores,
                                    center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,n,n)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Make the user distance calculation into a function for parallelization.
  f_user = function(x){
    count = 0
    flag  = FALSE
    for(i in 1:n){
      for(j in (i+1):n){
        count = count + 1
        if(count == x){
          pdf_1      = densities_matrix[i,]
          pdf_2      = densities_matrix[j,]
          dist_value = user_dist_function(pdf_1, pdf_2)
          flag       = TRUE
          break
        }
      }
      if(flag){
        break
      }
    }
    return(list(i = i, j = j, value = dist_value))
  }
  
  if(num_cores > 1){
    # Parallelize calculating each element of the distance matrix.
    num_of_unique_elements = n*(n-1)/2
    par_vector             = mclapply(1:num_of_unique_elements, f_user, mc.cores = num_cores)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
    
  } else {
    # Do not parallelize
    num_of_unique_elements = n*(n-1)/2
    par_vector             = lapply(1:num_of_unique_elements, f_user)
    for(par_idx in 1:num_of_unique_elements){
      dist_matrix[par_vector[[par_idx]]$i,par_vector[[par_idx]]$j] = par_vector[[par_idx]]$value
    }
  }
  
  dist_matrix = dist_matrix + t(dist_matrix)
  return(dist_matrix)
}


#' Calculates median and distance from the median using user-defined distance metric.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. rows are densities and columns are values at x_grid values
#' @param user_dist_function R Function.  User-defined.
#' @param num_cores Integer.  Number of cores to use if parallelizing.
#'
#' @noRd
#' 
get_user_defined_dist_from_median = function(x_grid, 
                                             densities_matrix,
                                             user_dist_function,
                                             num_cores,
                                             center_PDFs){
  n                     = nrow(densities_matrix)
  p                     = ncol(densities_matrix)
  dist_matrix           = matrix(0,1,n)
  orig_densities_matrix = densities_matrix
  
  pdfs_matrix = matrix(0,n,p)
  if(center_PDFs){
    # In this case, I'll align the mode of every PDF.
    # Start by transforming the densities_matrix into the cdf space:
    f_par = function(x){center_pdf(densities_matrix[x,], x_grid)}
    if(num_cores>1){
      pdfs_list = mclapply(1:n, f_par, mc.cores = num_cores)
    } else {
      pdfs_list = lapply(1:n, f_par)
    }
    for(i_idx in 1:n){
      pdfs_matrix[i_idx,] = pdfs_list[[i_idx]]
    }
    densities_matrix = pdfs_matrix
  }
  
  # Get the cross-sectional median:
  cross_median = apply(densities_matrix, 2, median)
  for(par_idx in 1:n){
    pdf_1      = densities_matrix[par_idx,]
    pdf_2      = cross_median
    dist_value = user_dist_function(pdf_1, pdf_2)
    dist_matrix[1,par_idx] = dist_value
  }
  
  # Again, I'll let the median on the non-centered space be the true cross median.
  if(center_PDFs){
    # Get the cross-sectional median:
    cross_median = apply(orig_densities_matrix, 2, median)
  }
  
  return(list(dist_matrix = dist_matrix, median_curve = cross_median, median_in_set = FALSE))
}

center_pdf = function(pdf, x_grid){
  # Pretend i have a single density for now.
  shifted_dens = rep(0, times = length(x_grid))
  idx_of_mode  = which.max(pdf)
  shift_num    = floor(length(x_grid)/2) - idx_of_mode
  
  if(shift_num < 0){
    # Then we shifted all values to the left.
    shifted_dens[1:(length(pdf)+shift_num)] = pdf[(-shift_num+1):length(pdf)]
    
  } else if (shift_num > 0) {
    # Then we shifted all values to the right.
    shifted_dens[(shift_num+1):length(pdf)] = pdf[1:(length(pdf)-shift_num)]
  }
  
  # Okay i really need to debug the above.
  
  # So I could linearly interpolate values to make this shift more concious
  # of cutting off values at the edges, but I really think that this will be
  # unnecessary.  Maybe I'll just make a check that catches if the majority of
  # the mass is as the edge and tell the researcher to "not center, or center
  # yourself."
  
  # Next, renormalize the density.
  normalized_dens = normalize_density(shifted_dens, x_grid)
  
  return(normalized_dens)
}

#' Normalize a pdf
#'
#' @param pdf Vector.  Single pdf with same length as x_grid.
#' @param x_grid Vector. X values of the PDF
#' @param norm Logical. Whether or not to normalize the CDF.
#'
#' @noRd
#' 
normalize_density = function(pdf, x_grid){
  cdf = cumtrapz(x_grid, pdf)
  return(pdf / cdf[length(cdf)])
}

#' Convert pdf to cdf
#'
#' @param pdf Vector.  Single pdf with same length as x_grid.
#' @param x_grid Vector. X values of the PDF
#' @param norm Logical. Whether or not to normalize the CDF.
#'
#' @noRd
#' 
pdf_to_cdf = function(pdf, x_grid, norm=TRUE){
  cdf = cumtrapz(x_grid, pdf)
  if(norm){ 
    cdf = cdf/cdf[length(cdf)] 
  }
  return(cdf)
}

#' Convert cdf to pdf
#'
#' @param cdf Vector.  Single cdf with same length as x_grid. 
#' @param x_grid Vector. X values of the PDF 
#'
#' @noRd
#' 
cdf_to_pdf = function(cdf, x_grid){
  pdf = gradient(as.numeric(cdf), x_grid)
  return(pdf)
}

#
#'  Convert cdf to quantile function. Uses linear interpolation
#'
#' @param cdf Vector.  Single cdf with same length as x_grid.  
#' @param x_grid Vector. X values of the PDF  
#'
#' @noRd
#' 
cdf_to_quant = function(cdf, x_grid){
  interp_fn = stats::approxfun(cdf, x_grid)
  quant = interp_fn(x_grid)
  return(quant)
}


#' Convert log quantile density to cdf.  Uses linear interpolation
#'
#' @param Vector. Single lqd with same length as x_grid.  
#' @param Vector. X values of the PDF  
#' @param lb Float. Linear bias to add to transformation.
#'
#' @noRd
#' 
lqd_to_cdf = function(lqd, x_grid, lb=0){
  p_grid = seq(0, 1, length.out=length(lqd))
  F_inv = lb + as.numeric(cumtrapz(p_grid, exp(lqd))/trapz(p_grid, exp(lqd)))
  cdf = stats::approx(F_inv, p_grid, xout=x_grid, yleft=0, yright=1)
  
  return(cdf$y)
}


#' Convert pdf to normalized LQD. Uses linear interpolation
#'
#' @param pdf Vector.  Single pdf with same length as x_grid.  
#' @param x_grid x_grid Vector. X values of the PDF  
#'
#' @noRd
#' 
pdf_to_lqd = function(pdf, x_grid){
  interp_fn = approxfun(x_grid, pdf)
  cdf       = pdf_to_cdf(pdf, x_grid, norm=TRUE)
  quant     = cdf_to_quant(cdf, x_grid)
  
  lqd  = -log(interp_fn(quant))
  nlqd = lqd
  
  return(nlqd)
}


#' Convert pdf to normalized LQD. Uses linear interpolation
#'
#' @param pdf Vector.  Single pdf with same length as x_grid.  
#' @param quant  Vector.  Single quantile function with same length as x_grid.  
#' @param x_grid x_grid Vector. X values of the PDF  
#'
#' @noRd
#' 
pdf_to_nlqd = function(pdf, x_grid){
  interp_fn = approxfun(x_grid, pdf)
  cdf       = pdf_to_cdf(pdf, x_grid, norm=TRUE)
  quant     = cdf_to_quant(cdf, x_grid)

  lqd  = -log(interp_fn(quant))
  nlqd = lqd / sum(lqd)
  
  return(nlqd)
}


#' Convert pdf to centered log ratio.
#'
#' @param pdf Vector.  Single pdf with same length as x_grid.  
#' @param x_grid x_grid Vector. X values of the PDF  
#'
#' @noRd
#' 
pdf_to_clr = function(pdf, x_grid){
  # Note that due to the log function here, you will likely want to mix the pdf
  # a uniform pdf before using this method.
  d_pdf = mean(diff(x_grid))
  average_log = d_pdf * sum(log(pdf)) / length(x_grid)
  clr         = log(pdf) - average_log
  return(clr)
}

#' alpha-mixture pdf with uniform
#'
#' @param pdf Vector.  Single pdf with same length as x_grid.  
#' @param alpha Float in (0,1).  Mixture parameter for pdf and uniform pdf.
#'
#' @noRd
#' 
alpha_mix = function(pdf, alpha){
  pdf_star = (1-alpha)*pdf + alpha*1
  return(pdf_star)
}

#' Recover original pdf from alpha-mixed pdf
#'
#' @param pdf_star Vector.  Single alpha-mixed pdf with same length as x_grid.  
#' @param alpha Float in (0,1).  Mixture parameter for pdf and uniform pdf.
#' @param x_grid x_grid Vector. X values of the PDF  
#'
#' @noRd
#' 
recover_pdf = function(pdf_star, alpha, x_grid){
  W = trapz(x_grid, abs((pdf_star-alpha)/(1-alpha)))
  pdf = abs((pdf_star-alpha)/(1-alpha))/W
  return(pdf)
}

#' Create the ggplot graph object that graphs the different regions via shading.
#'
#' @param curve_data Matrix. Data created by summarize_densities.R.
#'
#' @noRd
#' 
graph_polygons = function(curve_data){
  x       = NULL
  y       = NULL
  id      = NULL
  density = NULL
  
  positions3 = curve_data[which(curve_data$type == "outlier"),]
  if(!(nrow(positions3)==0)){
    positions3$y = positions3$density
    positions3$density = NULL
    p = ggplot(positions3, aes(x = x, y = y)) + geom_line(aes(group = id), color = "#D55E00", linetype = "dashed")
  }
  
  
  # The outside center 50%:
  lower_50_vertical_bounds = curve_data[which(curve_data$type == "Lower of outside center 50%"),]
  upper_50_vertical_bounds = curve_data[which(curve_data$type == "Upper of outside center 50%"),]
  x_vals = unique(lower_50_vertical_bounds$x)
  positions = NULL
  for( idx in 1:(length(x_vals)-1) ){
    x1 = x_vals[idx]
    x2 = x_vals[idx+1]
    y1 = lower_50_vertical_bounds[which(lower_50_vertical_bounds$x == x1),]$density
    y2 = upper_50_vertical_bounds[which(upper_50_vertical_bounds$x == x1),]$density
    y3 = lower_50_vertical_bounds[which(lower_50_vertical_bounds$x == x2),]$density
    y4 = upper_50_vertical_bounds[which(upper_50_vertical_bounds$x == x2),]$density
    
    position = data.frame( id = rep(idx, times = 4),
                           x = c(x1,x1,x2,x2),
                           y = c(y1,y2,y4,y3) )
    positions = rbind(positions, position)
  }
  if(nrow(positions3)==0){
    p = ggplot(positions, aes(x = x, y = y)) + geom_polygon(data = positions, aes(group = id), fill = 'blue', alpha = 0.4)
  } else {
    p = p + geom_polygon(data = positions, aes(group = id), fill = 'blue', alpha = 0.3)
  }
  
  
  # The center 50%:
  lower_50_vertical_bounds = curve_data[which(curve_data$type == "Lower of center 50%"),]
  upper_50_vertical_bounds = curve_data[which(curve_data$type == "Upper of center 50%"),]
  x_vals = unique(lower_50_vertical_bounds$x)
  positions2 = NULL
  for( idx in 1:(length(x_vals)-1) ){
    x1 = x_vals[idx]
    x2 = x_vals[idx+1]
    y1 = lower_50_vertical_bounds[which(lower_50_vertical_bounds$x == x1),]$density
    y2 = upper_50_vertical_bounds[which(upper_50_vertical_bounds$x == x1),]$density
    y3 = lower_50_vertical_bounds[which(lower_50_vertical_bounds$x == x2),]$density
    y4 = upper_50_vertical_bounds[which(upper_50_vertical_bounds$x == x2),]$density
    
    position = data.frame( id = rep(idx, times = 4),
                            x = c(x1,x1,x2,x2),
                            y = c(y1,y2,y4,y3) )
    positions2 = rbind(positions2, position)
  }
  p = p + geom_polygon(data = positions2, aes(group = id), fill = 'yellow', alpha = 0.6)
  
  # Now draw the median curve:
  positions4 = curve_data[which(curve_data$type == "Median"), ]
  p = p + geom_line(data = positions4, aes(x = x, y = density), color = 'black', linewidth=2) + 
      xlab("x_grid") + 
      ylab("density") +
      ggtitle(paste("Regions according to order determined by chosen distance")) + 
    theme_bw()
  

  # And return the full graph:
  return(p)
}

