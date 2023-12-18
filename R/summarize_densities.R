# Orders functions, creates boxplots, calculates outliers for PDFs.
## Author: Alexander C. Murph
## Date: October 2023
require(ggplot2)
library(gridExtra)
library(dplyr)

#' Orders a data-set consisting of probability density functions on the same 
#' x-grid.  
#' Visualizes a boxplot of these functions based on the notion of
#' distance determined by the user.  
#' Reports outliers based on the distance chosen and k value.
#'
#' @param x_grid Vector. X values of the PDF
#' @param densities_matrix Matrix. A n x p matrix where rows are individual PDFs and p matches the length of x_grid.
#' @param distance Character. The distance metric to use for the pairwise distances, or one of the two band depth options.
#' @param median_type Character. Whether the cross-median or the geometric median should be used.
#' @param center_PDFs Logical. Whether or not the modes of all the PDFs should be aligned prior to performing any calculations.
#' @param user_dist R Function. User-defined function that takes in two PDFs as vectors and returns a non-negative float corresponding to a distance between them.
#' @param k Float.  The factor by which to expand the IQR when calculating outliers.
#' @param num_cores Integer.  The number of cores to use if parallelizing the distance matrix calculations.
#'
#' @return An deboinr object containing the following:
#' \itemize{
#'   \item density_order. Vector of indices corresponding to rows of densities_matrix that sort from closest to furthest from the median PDF.
#'   \item outliers. Vector of indices corresponding to rows of densities_matrix that are determined to be outliers.
#'   \item box_plot. ggplot object of graphic output by calling this method.
#' }
#' @export
#'
#' @examples
#'
#' example_data = DeBoinR::pdf_data[1:100,]
#' xx = deboinr(DeBoinR::x_grid,
#'              as.matrix(example_data),
#'              distance = "hellinger",
#'              median_type = 'cross',
#'              center_PDFs = TRUE,
#'              num_cores = 1
#' )
#' 
#' print("about to print DeBoinR object...")
#' print(xx)
deboinr = function(x_grid, 
                   densities_matrix, 
                   distance = c("hellinger",
                                "nLQD",
                                "fisher_rao",
                                "TV_dist",
                                "CLR",
                                "wasserstein",
                                "BD_fboxplot",
                                "MBD_fboxplot",
                                "user_defined"),
                   median_type = c("cross",
                                   "geometric"),
                   center_PDFs = FALSE,
                   user_dist = NULL,
                   k = 1.5,
                   num_cores = 1
                   ){
  if( !("matrix"%in%class(densities_matrix)) ) stop("densities_matrix must be a matrix.")
  if(length(x_grid)!=ncol(densities_matrix)) stop("The length of x_grid must equal the number of columns in densities_matrix.")
  if(distance=="user_defined"&is.null(user_dist)) stop("Please provide an R function for user_dist.")
  if(length(distance)>1) distance = "hellinger"
  if(length(median_type)>1) median_type = "cross"
  
  type    = NULL
  x       = NULL
  density = NULL
  num     = NULL

  # For use in the geometric calculation:
  pairwise_distances    = NULL
  # For use in the cross-sectional calculation
  distances_from_median = NULL
  
  # Get pairwise distances between densities, or the return values from fboxplot package.
  if(distance == "hellinger"){
    if(median_type=='geometric'){
      pairwise_distances    = get_hellinger_dist_mat(x_grid, densities_matrix, num_cores, center_PDFs)
    } else {
      distances_from_median = get_hellinger_dist_from_median(x_grid, densities_matrix, num_cores, center_PDFs)
    }
  } else if (distance == "nLQD"){
    if(median_type=='geometric') {
      pairwise_distances    = get_nLQD_dist_mat(x_grid, densities_matrix, num_cores, center_PDFs)
    } else {
      distances_from_median = get_nLQD_dist_from_median(x_grid, densities_matrix, num_cores, center_PDFs)
    }
  } else if (distance == "fisher_rao"){
    if(median_type=='geometric') {
      pairwise_distances    = get_fisher_rao_dist_mat(x_grid, densities_matrix, num_cores, center_PDFs)
    } else {
      distances_from_median = get_fisher_rao_dist_from_median(x_grid, densities_matrix, num_cores, center_PDFs)
    }
  } else if (distance == "TV_dist"){
    if(median_type=='geometric') {
      pairwise_distances    = get_TV_dist_dist_mat(x_grid, densities_matrix, num_cores, center_PDFs)
    } else {
      distances_from_median = get_TV_dist_dist_from_median(x_grid, densities_matrix, num_cores, center_PDFs)
    }
  } else if (distance == "CLR"){
    if(median_type=='geometric') {
      pairwise_distances    = get_CLR_dist_mat(x_grid, densities_matrix, num_cores, center_PDFs)
    } else {
      distances_from_median = get_CLR_dist_from_median(x_grid, densities_matrix, num_cores, center_PDFs)
    }
  } else if (distance == "wasserstein"){
    if(median_type=='geometric') {
      pairwise_distances    = get_wasserstein_dist_mat(x_grid, densities_matrix, num_cores, center_PDFs)
    } else {
      distances_from_median = get_wasserstein_dist_from_median(x_grid, densities_matrix, num_cores, center_PDFs)
    }
  } else if (distance == "BD_fboxplot"){
    BD_return_values = get_BD_fboxplot(x_grid, densities_matrix, k, num_cores, center_PDFs)
  } else if (distance == "MBD_fboxplot"){
    BD_return_values   = get_MBD_fboxplot(x_grid, densities_matrix, k, num_cores, center_PDFs)
  } else if (distance == "user_defined") {
    # Check to see if user put in a valid distance function:
    check_value = user_dist(densities_matrix[1,], densities_matrix[2,])
    if(!is.numeric(check_value)){stop("user_dist(densities_matrix[1,], densities_matrix[2,]) did not return a single numeric.")}
    if(check_value < 0){stop("user_dist(densities_matrix[1,], densities_matrix[2,]) did not return a non-negative numeric.")}
    if(median_type=='cross'){
      warning("Make sure the cross-sectional median belongs to the metric space defined by your distance function!  If not, use 'geometric' option.")
      distances_from_median = get_user_defined_dist_from_median(x_grid, densities_matrix, user_dist, num_cores, center_PDFs)
    } else {
      pairwise_distances    = get_user_defined_dist_mat(x_grid, densities_matrix, user_dist, num_cores, center_PDFs)
    }
  } else {
    stop("Invalid value for distance.")
  }
  my_object = list()
  epsilon   = 1e-10
    
  # Determine which type of median calculation we are doing.
  if( !is.null(distances_from_median) ){
    ## In this scenario, we calculate the cross-sectional median
    if(distances_from_median$median_in_set){
      dist_from_median   = distances_from_median$dist_matrix
      median_curve_index = min(which.min(dist_from_median))
    } else {
      densities_matrix   = rbind(densities_matrix, distances_from_median$median_curve)
      median_curve_index = nrow(densities_matrix)
      dist_from_median   = matrix(c(distances_from_median$dist_matrix[1,], 0), nrow = 1)
    }
    order_of_curves    = order(dist_from_median)
    center_50_percent  = order_of_curves[1:floor(length(order_of_curves)/2)]
    
    # Find outliers.
    IQR                = max(dist_from_median[center_50_percent])
    outlier_cutoff     = k*IQR + IQR
    idx_of_outliers    = which(dist_from_median > outlier_cutoff)
    
  } else if (!is.null(pairwise_distances)) {
    # With the pairwise distances, determine the median as the geometric.
    summed_differences = colSums(pairwise_distances)
    
    # Use the median to order the functions.
    median_curve_index = min(which.min(summed_differences))
    dist_from_median   = pairwise_distances[median_curve_index,]
    
    order_of_curves    = order(dist_from_median)
    center_50_percent  = order_of_curves[1:floor(length(order_of_curves)/2)]
    
    # Find outliers.
    IQR                = max(dist_from_median[center_50_percent])
    outlier_cutoff     = k*IQR + IQR
    idx_of_outliers    = which(dist_from_median > outlier_cutoff)
  } else {
    # Use the output from the BD depth function to get the ordering of the curves.
    functional_depth   = BD_return_values$depth_output$depth
    order_of_curves    = order(functional_depth, decreasing = TRUE)
    center_50_percent  = order_of_curves[1:floor(length(order_of_curves)/2)]
    median_curve_index = BD_return_values$depth_output$medcurve
    # Find outliers.
    # IQR                = min(functional_depth[center_50_percent])
    
    idx_of_outliers    = BD_return_values$depth_output$outpoint
    # outlier_cutoff     = IQR - k*IQR
    
  }

  # Create dataframes for each of the graph objects
  ### Data to print all covers, but color them by ordering.
  curves_by_colors = NULL
  n                = nrow(densities_matrix)
  n_grid           = length(x_grid)

  for(idx in 1:n){
    if(idx%in%idx_of_outliers) {
      temp_row = data.frame(density = densities_matrix[idx,], x = x_grid, 
                            type = rep("outlier", times = n_grid),
                            num = rep(idx, times = n_grid))
    } else if(idx%in%center_50_percent){
      temp_row = data.frame(density = densities_matrix[idx,], x = x_grid, 
                            type = rep("center 50%", times = n_grid),
                            num = rep(idx, times = n_grid))
      if(idx == min(median_curve_index)){
        temp_row$type = "median"
      }
    } else {
      temp_row = data.frame(density = densities_matrix[idx,], x = x_grid, 
                            type = rep("outside 50%, non-outlier", times = n_grid),
                            num = rep(idx, times = n_grid))
    }
    curves_by_colors = rbind(curves_by_colors, temp_row)
  }
  # I want the median to be at the end...
  idx_of_med       = which(curves_by_colors$type == "median")
  med_row          = curves_by_colors[idx_of_med, ]
  curves_by_colors = curves_by_colors[-idx_of_med, ]
  curves_by_colors = rbind(curves_by_colors, med_row)
  
  ### Data to create VERTICAL bands:
  # for center 50%
  center_50_data  = curves_by_colors[which(curves_by_colors$type == "center 50%"),]
  outlier_data    = curves_by_colors[which(curves_by_colors$type == "outlier"),]
  outside_50_data = curves_by_colors[which(curves_by_colors$type == "outside 50%, non-outlier"),]
  
  vertical_bounds = NULL

  for(x_val in x_grid){
    # Get upper and lower bounds on center 50%:
    temp_data       = center_50_data[which(center_50_data$x == x_val),]
    temp_row        = data.frame(density = min(temp_data$density), x = x_val, type = "Lower of center 50%", id = 1)
    vertical_bounds = rbind(vertical_bounds, temp_row)
    temp_row        = data.frame(density = max(temp_data$density), x = x_val, type = "Upper of center 50%", id = 1)
    vertical_bounds = rbind(vertical_bounds, temp_row)
    
    # Get upper and lower bounds on outside 50% that aren't outliers:
    temp_data       = outside_50_data[which(outside_50_data$x == x_val),]
    temp_row        = data.frame(density = min(temp_data$density), x = x_val, type = "Lower of outside center 50%", id = 1)
    vertical_bounds = rbind(vertical_bounds, temp_row)
    temp_row        = data.frame(density = max(temp_data$density), x = x_val, type = "Upper of outside center 50%", id = 1)
    vertical_bounds = rbind(vertical_bounds, temp_row)
  }
  
  if(length(which(curves_by_colors$type == "outlier"))>0){
    outlier_data     = curves_by_colors[which(curves_by_colors$type == "outlier"),]
    outlier_data$id  = outlier_data$num
    outlier_data$num = NULL
    vertical_bounds  = rbind(vertical_bounds, outlier_data)
  }
  
  median_curve_data  = curves_by_colors[which(curves_by_colors$type == "median"),]
  median_curve       = data.frame(density = median_curve_data$density, x = x_grid, type = "Median", id = 1)
  vertical_bounds    = rbind(vertical_bounds, median_curve)

  # Create the graph object.
  if(nrow(outlier_data)==0){
    p_curves = ggplot(curves_by_colors %>% filter(type != "median")) + 
      geom_line(aes(x = x, y = density, color = type, 
                    alpha = type, linewidth = type, group = num)) +
      geom_line(data = curves_by_colors %>% filter(type == "median"), aes(x = x, y = density, color = type, 
                                                                          alpha = type, linewidth = type, group = num)) + 
      ggtitle("Classified Raw Curves") + 
      theme_bw() + 
      xlab("x_grid") + 
      ylab("density") +
      scale_color_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier"),
                         values=c("black", "yellow", "blue")) + 
      scale_alpha_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier"),
                         values=c(1, 1, 0.3)) + 
      scale_linewidth_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier"),
                             values=c(2, 0.5, 0.5)) +
      guides(alpha = "none") +
      guides(linewidth = "none")  
  } else {
    p_curves = ggplot(curves_by_colors %>% filter(type != "median") ) + 
      geom_line(aes(x = x, y = density, color = type, 
                    alpha = type, linetype = type, linewidth = type, group = num)) +
      geom_line(data = curves_by_colors %>% filter(type == "median"), aes(x = x, y = density, color = type, 
                    alpha = type, linetype = type, linewidth = type, group = num)) +
      ggtitle("Classified Raw Curves") + 
      theme_bw() + 
      xlab("x_grid") + 
      ylab("density") +
      scale_color_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier", "outlier"),
                         values=c("black", "yellow", "blue", "#D55E00") ) + 
      scale_alpha_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier", "outlier"),
                         values=c(1, 1, 0.3, 1)) + 
      scale_linewidth_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier", "outlier"),
                         values=c(2, 0.5, 0.5, 0.5)) + 
      scale_linetype_manual(breaks = c("median", "center 50%", "outside 50%, non-outlier", "outlier"),
                             values=c("solid", "solid", "solid", "dashed")) +
      guides(alpha = "none") +
      guides(linewidth = "none") +
      guides(linetype = "none") 
  }
  
  p_areas     = graph_polygons(vertical_bounds)
  lst_p       = list(p_curves, p_areas)
  
  # Compile things into a deboinr object.
  my_object$box_plot       = lst_p
  my_object$outliers       = idx_of_outliers
  my_object$density_order  = order_of_curves
  class(my_object)         = "DeBoinR"
  
  # Return the object:
  return(my_object)
  
}


#' Print function for a DeBoinR object.  Prints ggplot graphs and other output values.
#'
#'
#' @param x deboinr object. Fit from DeBoinR main method.
#' @param ... Additional plotting arguments.
#'
#' @exportS3Method
print.DeBoinR = function(x, ...)
{
  x_grid           = NULL
  k                = NULL
  densities_matrix = NULL
  
  lst_p = x$box_plot
  
  gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))
  
  cat("\n     DeBoinR object\n")
  cat("----------------------------------\n")
  cat("--> Order of Densities: \n")
  print(x$density_order)
  cat("----------------------------------\n")
  cat("--> Indices of Outliers: \n")
  print(x$outliers)
}





