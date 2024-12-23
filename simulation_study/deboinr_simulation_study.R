# DeBoinR Simulation Study
## Author: Alexander Murph
## October 2023
setwd("~/DeBoinR/R")
source("summarize_densities.R")
source("helpers.R")
library(dgumbel)
library(reshape) 
library(gridExtra)
library(parallel)

# We'll perform this as a simulation study:
# sim_num = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# print(paste("Running simulation number:", sim_num))

sim_num = commandArgs()
print(sim_num)
print(sim_num[6])
sim_num = as.integer(sim_num[6])

if(is.na(sim_num)){
  sim_num = 3
}


set.seed(sim_num)

# For now, I'll set this.  In the future I should figure out how to detect it
# on the Chimera system.
num_cores = 1

# Parameters for all simulations.
n_of_densities        = 200
num_of_outliers       = 10

# Parameters for all simulations.
n_of_densities        = 200
num_of_outliers       = 10
num_of_central_curves = n_of_densities - num_of_outliers
x_grid                = seq(from=0, to=1, length.out = 200)
distances             = c("hellinger", "nLQD", "fisher_rao", "TV_dist", "CLR",
                          "wasserstein", "BD_fboxplot", "MBD_fboxplot")
median_types          = c("cross", "geometric")
normalizations        = c(TRUE, FALSE)

###################################
##### Create the two datasets #####
###################################

### Shift Outlier
# Create varying pdf data where one curve is just shifted.
magnitude_of_shift = 0.07
pdf_data_shift = NULL
for(idx in 1:num_of_central_curves){
  noise_around_location = rnorm(1,0,sd=0.01)
  noise_around_scale    = rnorm(1,0,sd=0.01)
  temp_pdf              = dgumbel( x_grid, 
                                   location = (0.25 + noise_around_location), 
                                   scale = (0.05 + noise_around_scale) )
  pdf_data_shift        = rbind(pdf_data_shift, temp_pdf)
}

for(idx in 1:num_of_outliers){
  noise_around_location = rnorm(1,0,sd=0.01)
  noise_around_scale    = rnorm(1,0,sd=0.01)
  temp_pdf              = dgumbel( x_grid, 
                                   location = (0.25 + magnitude_of_shift + noise_around_location), 
                                   scale = (0.05 + noise_around_scale) )
  pdf_data_shift        = rbind(pdf_data_shift, temp_pdf)
}
# melt_pdf         = melt(pdf_data_shift)
# melt_pdf$x_grid  = rep(x_grid, each = n_of_densities)
# melt_pdf$pdf_num = rep(1:n_of_densities, times = length(x_grid))
# melt_pdf$outlier = rep(c(rep(0,times=num_of_central_curves), rep(1,times=num_of_outliers)), times = length(x_grid))
# g1 = ggplot(melt_pdf, aes(x = x_grid, y = value, color = as.factor(outlier), group = as.factor(pdf_num))) +
#                       geom_line()  + theme_bw() + theme(legend.position = 'none', axis.title = element_text(size = 15)) +
#                       scale_colour_manual(values = c("0" = "grey", "1" = "black")) + xlab("") + ylab("Horizontal-Shift Outliers")

### Shape Outlier
# Create varying pdf data where one curve is just shifted.
magnitude_of_shift = 0.05
pdf_data_shape     = NULL
for(idx in 1:num_of_central_curves){
  noise_around_location = rnorm(1,0,sd=0.01)
  noise_around_scale    = rnorm(1,0,sd=0.01)
  temp_pdf              = dgumbel( x_grid, 
                                   location = (0.25 + noise_around_location), 
                                   scale = (0.05 + noise_around_scale) )
  pdf_data_shape        = rbind(pdf_data_shape, temp_pdf)
}

for(idx in 1:num_of_outliers){
  noise_around_location = rnorm(1,0,sd=0.01)
  noise_around_scale    = rnorm(1,0,sd=0.01)
  alpha                 = rbeta(1,3,3)
  temp_pdf              = alpha * dgumbel( x_grid, 
                                           location = (0.25 + magnitude_of_shift + noise_around_location), 
                                           scale = (0.05 + noise_around_scale) ) +
                          (1-alpha) * dgumbel( x_grid, 
                                               location = (0.25 - magnitude_of_shift + noise_around_location), 
                                               scale = (0.05 + noise_around_scale) ) 
  pdf_data_shape        = rbind(pdf_data_shape, temp_pdf)
}
# melt_pdf         = melt(pdf_data_shape)
# melt_pdf$x_grid  = rep(x_grid, each = n_of_densities)
# melt_pdf$pdf_num = rep(1:n_of_densities, times = length(x_grid))
# melt_pdf$outlier = rep(c(rep(0,times=num_of_central_curves), rep(1,times=num_of_outliers)), times = length(x_grid))
# g2 = ggplot(melt_pdf, aes(x = x_grid, y = value, color = as.factor(outlier), group = as.factor(pdf_num))) +
#   geom_line()  + theme_bw() + theme(legend.position = 'none', axis.title = element_text(size = 15) ) +
#   scale_colour_manual(values = c("0" = "grey", "1" = "black")) + xlab("") + ylab("Shape Outliers")
# 
# lst_p = list(g1,g2)
# gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))


#######################################################
##### See which methods can detect these outliers #####
#######################################################

f_perform_detection = function(x, data, grid, sim_number, outlier_type){
  dist_num = (x-1)%%8 + 1
  med_num  = as.numeric(x > 16) + 1
  center   = c(rep(1,times=8),rep(2,times=8),rep(1,times=8),rep(2,times=8))
  
  xx = deboinr(x_grid,
               as.matrix(data),
               distance = distances[dist_num],
               median_type = median_types[med_num],
               center_PDFs = normalizations[center[x]],
               num_cores = 1
  )
  
  outliers      = xx$outliers
  density_order = xx$density_order
  
  write.csv(outliers, file=paste("outliers_data/sim_num", sim_number, "dist", 
                                 distances[dist_num], "median", median_types[med_num], 
                                 "outlier", outlier_type, "centered", center[x], ".csv", sep="_"))
  
  write.csv(density_order, file=paste("density_orders/sim_num", sim_number, "dist", 
                                 distances[dist_num], "median", median_types[med_num], 
                                 "outlier", outlier_type, "centered", center[x], ".csv", sep="_"))
}

xx = deboinr(x_grid,
             as.matrix(pdf_data_shift),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = FALSE,
             num_cores = 1
)
g_wass_geo_notcentered_shift = xx$box_plot[[2]] + xlab("") +
  ggtitle("Horizontal-Shift Outliers") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + ylab("")

xx = deboinr(x_grid,
             as.matrix(pdf_data_shift),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = TRUE,
             num_cores = 1
)
g_wass_geo_centered_shift = xx$box_plot[[2]] + xlab("") + ggtitle("") + ylab("")

xx = deboinr(x_grid,
             as.matrix(pdf_data_shift),
             distance = "MBD_fboxplot",
             median_type = 'cross',
             center_PDFs = FALSE,
             num_cores = 1
)
g_mbd_cross_notcentered_shift = xx$box_plot[[2]] + ggtitle("") + ylab("") 

xx = deboinr(x_grid,
             as.matrix(pdf_data_shape),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = FALSE,
             num_cores = 1
)
g_wass_geo_notcentered_shape = xx$box_plot[[2]] + ylab("") + xlab("") +
  ggtitle("Shape Outliers") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

xx = deboinr(x_grid,
             as.matrix(pdf_data_shape),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = TRUE,
             num_cores = 1
)
g_wass_geo_centered_shape = xx$box_plot[[2]] + ylab("") + xlab("") + ggtitle("")

xx = deboinr(x_grid,
             as.matrix(pdf_data_shape),
             distance = "MBD_fboxplot",
             median_type = 'cross',
             center_PDFs = FALSE,
             num_cores = 1
)
g_mbd_cross_notcentered_shape = xx$box_plot[[2]] + ylab("") + ggtitle("")


lst_p = list(g_wass_geo_notcentered_shift, g_wass_geo_centered_shift, g_mbd_cross_notcentered_shift,
             g_wass_geo_notcentered_shape, g_wass_geo_centered_shape, g_mbd_cross_notcentered_shape)

gridExtra::grid.arrange(lst_p[[1]],lst_p[[2]], lst_p[[3]], lst_p[[4]], lst_p[[5]], lst_p[[6]], 
                        layout_matrix = matrix(c(1:6), byrow = FALSE, ncol = 2))

grid.arrange(arrangeGrob(lst_p[[1]], lst_p[[4]], ncol = 2, 
                         left = 'Wasserstein Distance,\n Geometric Median,\n without Centering'), 
             arrangeGrob(lst_p[[2]], lst_p[[5]], ncol = 2, 
                         left = 'Wasserstein Distance,\n Geometric Median,\n with Centering'),
             arrangeGrob(lst_p[[3]], lst_p[[6]], ncol = 2, 
                         left = 'Modified Band Depth,\n Cross-Section Median,\n without Centering'), 
             ncol=1)




setwd("~/DeBoinR_study_files")
## All simulations for horizontal shift
f_sim_horizontal = function(x){ f_perform_detection(x, pdf_data_shift, x_grid, sim_num, "horizontalShift") }
par_vector       = lapply(1:32, f_sim_horizontal)

## All simulations for shape
f_sim_shape      = function(x){ f_perform_detection(x, pdf_data_shape, x_grid, sim_num, "shape") }
par_vector       = lapply(1:32, f_sim_shape)






