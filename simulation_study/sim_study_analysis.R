# Analysis of the Simulation Study Results
## Author: Alexander C. Murph
## November 2023
library(ggplot2)
library(ggpubr)
library(ggpattern)
library(reshape)

# Parameters for all simulations.
n_of_densities        = 200
all_densities         = 1:200
num_of_outliers       = 10
sim_number            = 3
n_of_simulations      = 100
num_of_central_curves = n_of_densities - num_of_outliers
true_classifiers      = c( rep(0, times = num_of_central_curves), 
                           rep(1, times = num_of_outliers) )
x_grid                = seq(from=0, to=1, length.out = 200)
renamed_distances     = c("Hellinger", "L2 on LQD", "Fisher-Rao", "TV", "Bayes",
                          "Wasserstein", "Band Depth", "Mod. Band Depth")
distances             = c("Hellinger", "nLQD", "fisher_rao", "TV_dist", "CLR",
                          "wasserstein", "BD_fboxplot", "MBD_fboxplot")
ordered_distances     = c("L2 on LQD","Bayes","Wasserstein","TV","Hellinger","Fisher-Rao", 
                           "Band Depth", "Mod. Band Depth")
median_types          = c("cross", "geometric")
centered              = c('centered', 'not centered')

################################################
## Let's look at the horizontal shift data first
horizontal_shift_data          = NULL
horizontal_shift_data_combined = NULL
horizontal_shift_data_metrics  = NULL

outlier_type          = 'horizontalShift'
for(sim_num in 1:n_of_simulations){
  for(exp_num in 1:8){
    for(cent_num in 1:2){
      dist_num                 = (exp_num-1)%%8 + 1
      outliers_temp_cross      = read.csv(file=paste("outliers_data/sim_num", sim_num, "dist", 
                                           distances[dist_num], "median", "cross", 
                                           "outlier", outlier_type, "centered", cent_num, ".csv", sep="_"))
      temp_row_cross           = data.frame(distance = renamed_distances[dist_num],
                                            median_type = 'cross',
                                            Centering = centered[cent_num],
                                            sim_number = sim_num)
      temp_row_cross           = cbind(temp_row_cross, t(as.data.frame(as.numeric(all_densities%in%outliers_temp_cross$x))))
      rownames(temp_row_cross) = NULL
      horizontal_shift_data    = rbind(horizontal_shift_data, temp_row_cross)
      
      outliers_temp_geometric      = read.csv(file=paste("outliers_data/sim_num", sim_num, "dist", 
                                                         distances[dist_num], "median", "geometric", 
                                                         "outlier", outlier_type, "centered", cent_num, ".csv", sep="_"))
      temp_row_geometric           = data.frame(distance = renamed_distances[dist_num],
                                                median_type = 'geometric',
                                                Centering = centered[cent_num],
                                                sim_number = sim_num)
      temp_row_geometric           = cbind(temp_row_geometric, t(as.data.frame(as.numeric(all_densities%in%outliers_temp_geometric$x))))
      rownames(temp_row_geometric) = NULL
      horizontal_shift_data        = rbind(horizontal_shift_data, temp_row_geometric)
      
      # Combine both metrics in a way to be used for the tile plots:
      row_1_info   = as.vector(unlist(temp_row_cross[c(as.character(1:length(x_grid)))]))
      row_2_info   = as.vector(unlist(temp_row_geometric[c(as.character(1:length(x_grid)))]))
      combined_row = rep('neither', times = length(x_grid))
      combined_row = as.vector(ifelse( ((row_1_info+row_2_info)==2), 'both', 'neither' ))
      
      combined_row[row_1_info==1] = ifelse( (combined_row[row_1_info==1]=='both'), 'both', 'cross')
      combined_row[row_2_info==1] = ifelse( (combined_row[row_2_info==1]=='both'), 'both', 'geometric')
      
      demographics            = data.frame(sim_number = sim_num, distance = renamed_distances[dist_num], Centering = centered[cent_num])
      combined_row            = cbind(demographics, t(as.data.frame(combined_row)))
      rownames(combined_row)  = NULL
      
      combined_row$sim_number        = as.integer(sim_num)
      combined_row$distance          = renamed_distances[dist_num]
      combined_row$Centering         = centered[cent_num]
      horizontal_shift_data_combined = rbind(horizontal_shift_data_combined, combined_row)
  
      # Collect the metrics data:
      accuracy_cross     = sum(row_1_info == true_classifiers) / n_of_densities
      tpr_cross          = sum(row_1_info[(num_of_central_curves+1):n_of_densities] == 1) / num_of_outliers
      fpr_cross          = sum(row_1_info[1:num_of_central_curves] == 1) / num_of_central_curves
      
      accuracy_geometric = sum(row_2_info == true_classifiers) / n_of_densities
      tpr_geometric      = sum(row_2_info[(num_of_central_curves+1):n_of_densities] == 1) / num_of_outliers
      fpr_geometric      = sum(row_2_info[1:num_of_central_curves] == 1) / num_of_central_curves
      temp_row1          = data.frame(sim_number = sim_num, distance = renamed_distances[dist_num], 
                                      Centering = centered[cent_num],
                                       median_type = 'cross', accuracy = accuracy_cross,
                                       tpr = tpr_cross, fpr = fpr_cross)
      temp_row2          = data.frame(sim_number = sim_num, distance = renamed_distances[dist_num], 
                                      Centering = centered[cent_num],
                                      median_type = 'geometric', accuracy = accuracy_geometric,
                                      tpr = tpr_geometric, fpr = fpr_geometric)
      
      horizontal_shift_data_metrics = rbind(horizontal_shift_data_metrics, temp_row1, temp_row2)
    }
  }
}

# # Let's try this grid plot idea I've had:
# temp_horizontal = horizontal_shift_data[which( (horizontal_shift_data$Centering=='centered')&
#                                                  (horizontal_shift_data$median_type=='cross')&
#                                                  (horizontal_shift_data$sim_number==sim_number) ),]
# temp_horizontal$median_type = NULL
# temp_horizontal$sim_number  = NULL
# temp_horizontal$Centering   = NULL
# temp_melt = melt(temp_horizontal)
# 
# variable = unique(temp_melt$variable)
# g1 = ggplot(temp_melt, aes(x=variable, y=distance)) +
#   geom_tile(aes(fill=value)) +
#   scale_fill_gradient(low="white", high="blue")+theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")+ 
#   scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) 
# 
# temp_horizontal = horizontal_shift_data[which( (horizontal_shift_data$median_type=='centered')&
#                                                  (horizontal_shift_data$median_type=='geometric')&
#                                                  (horizontal_shift_data$sim_number==sim_number) ),]
# temp_horizontal$median_type = NULL
# temp_horizontal$sim_number  = NULL
# temp_horizontal$Centering   = NULL
# temp_melt = melt(temp_horizontal)
# 
# variable = unique(temp_melt$variable)
# g2 = ggplot(temp_melt, aes(x=variable, y=distance)) +
#   geom_tile(aes(fill=value)) +
#   scale_fill_gradient(low="white", high="blue")+theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")+ 
#   scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) 
# 
# lst_p = list(g1,g2)
# gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))

## Let's make a tile plot with all of this information at once.
temp_melt = horizontal_shift_data_combined[which(horizontal_shift_data_combined$sim_number==sim_number),]
temp_melt$sim_number = NULL
temp_melt = melt(temp_melt, id = c(1,2))
variable = unique(temp_melt$variable)
ggplot(temp_melt, aes(x=variable, y=distance)) +
  geom_tile(aes(fill=as.factor(value))) + 
  scale_fill_manual(values = c('neither'='white', 'cross'='grey', 'geometric'='brown', 'both'='black'))+theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #, legend.position="none"
  scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) + guides(fill=guide_legend(title="Median Type")) +
  xlab("Density Observation") + ylab("Method")

temp_melt_horizontal_shift = temp_melt
variable_horizontal_shift  = variable


################################################
## Let's look at the shape data 
shape_data          = NULL
shape_data_combined = NULL
shape_data_metrics  = NULL
outlier_type        = 'shape'
for(sim_num in 1:n_of_simulations){
  for(exp_num in 1:8){
    for(cent_num in 1:2){
      dist_num                 = (exp_num-1)%%8 + 1
      
      # Stuff for cross median
      outliers_temp_cross      = read.csv(file=paste("outliers_data/sim_num", sim_num, "dist", 
                                                     distances[dist_num], "median", "cross", 
                                                     "outlier", outlier_type, "centered", cent_num, ".csv", sep="_"))
      temp_row_cross           = data.frame(distance = renamed_distances[dist_num], Centering = centered[cent_num],
                                            median_type = 'cross',
                                            sim_number = sim_num)
      temp_row_cross           = cbind(temp_row_cross, t(as.data.frame(as.numeric(all_densities%in%outliers_temp_cross$x))))
      rownames(temp_row_cross) = NULL
      shape_data               = rbind(shape_data, temp_row_cross)
      
      # Stuff for shape median
      outliers_temp_geometric      = read.csv(file=paste("outliers_data/sim_num", sim_num, "dist", 
                                                         distances[dist_num], "median", "geometric", 
                                                         "outlier", outlier_type, "centered", cent_num, ".csv", sep="_"))
      temp_row_geometric           = data.frame(distance = renamed_distances[dist_num], Centering = centered[cent_num],
                                                median_type = 'geometric',
                                                sim_number = sim_num)
      temp_row_geometric           = cbind(temp_row_geometric, t(as.data.frame(as.numeric(all_densities%in%outliers_temp_geometric$x))))
      rownames(temp_row_geometric) = NULL
      shape_data                   = rbind(shape_data, temp_row_geometric)
      
      row_1_info   = as.vector(unlist(temp_row_cross[c(as.character(1:length(x_grid)))]))
      row_2_info   = as.vector(unlist(temp_row_geometric[c(as.character(1:length(x_grid)))]))
      combined_row = rep('neither', times = length(x_grid))
      combined_row = as.vector(ifelse( ((row_1_info+row_2_info)==2), 'both', 'neither' ))
      
      combined_row[row_1_info==1] = ifelse( (combined_row[row_1_info==1]=='both'), 'both', 'cross')
      combined_row[row_2_info==1] = ifelse( (combined_row[row_2_info==1]=='both'), 'both', 'geometric')
      
      demographics            = data.frame(sim_number = sim_num, distance = renamed_distances[dist_num], Centering = centered[cent_num])
      combined_row            = cbind(demographics, t(as.data.frame(combined_row)))
      rownames(combined_row)  = NULL
      
      combined_row$sim_number = as.integer(sim_num)
      combined_row$distance   = renamed_distances[dist_num]
      combined_row$Centering  = centered[cent_num]
      shape_data_combined     = rbind(shape_data_combined, combined_row)
      
      # Collect the metrics data:
      accuracy_cross     = sum(row_1_info == true_classifiers) / n_of_densities
      tpr_cross          = sum(row_1_info[(num_of_central_curves+1):n_of_densities] == 1) / num_of_outliers
      fpr_cross          = sum(row_1_info[1:num_of_central_curves] == 1) / num_of_central_curves
      
      accuracy_geometric = sum(row_2_info == true_classifiers) / n_of_densities
      tpr_geometric      = sum(row_2_info[(num_of_central_curves+1):n_of_densities] == 1) / num_of_outliers
      fpr_geometric      = sum(row_2_info[1:num_of_central_curves] == 1) / num_of_central_curves
      temp_row1          = data.frame(sim_number = sim_num, distance = renamed_distances[dist_num], Centering = centered[cent_num], 
                                      median_type = 'cross', accuracy = accuracy_cross,
                                      tpr = tpr_cross, fpr = fpr_cross)
      temp_row2          = data.frame(sim_number = sim_num, distance = renamed_distances[dist_num], Centering = centered[cent_num], 
                                      median_type = 'geometric', accuracy = accuracy_geometric,
                                      tpr = tpr_geometric, fpr = fpr_geometric)
      shape_data_metrics = rbind(shape_data_metrics, temp_row1, temp_row2)
    }
  }
}

# # Let's try this grid plot idea I've had:
# temp_horizontal             = shape_data[which( (shape_data$median_type=='cross')&(shape_data$sim_number==sim_number) ),]
# temp_horizontal$median_type = NULL
# temp_horizontal$sim_number  = NULL
# temp_melt                   = melt(temp_horizontal)
# variable                    = unique(temp_melt$variable)
# 
# g1 = ggplot(temp_melt, aes(x=variable, y=distance)) +
#   geom_tile(aes(fill=value)) +
#   scale_fill_gradient(low="white", high="black")+theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")+ 
#   scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) 
# 
# temp_horizontal = shape_data[which( (shape_data$median_type=='geometric')&(shape_data$sim_number==sim_number) ),]
# temp_horizontal$median_type = NULL
# temp_horizontal$sim_number = NULL
# temp_melt = melt(temp_horizontal)
# 
# variable = unique(temp_melt$variable)
# g2 = ggplot(temp_melt, aes(x=variable, y=distance)) +
#   geom_tile(aes(fill=value)) +
#   scale_fill_gradient(low="white", high="black")+theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")+ 
#   scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) 
# 
# lst_p = list(g1,g2)
# gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))

temp_melt = shape_data_combined[which(shape_data_combined$sim_number==sim_number),]
temp_melt$sim_number = NULL
temp_melt = melt(temp_melt, id = c(1,2))
variable = unique(temp_melt$variable)
ggplot(temp_melt, aes(x=variable, y=distance)) +
  geom_tile(aes(fill=as.factor(value))) + scale_fill_manual(values = c('neither'='white', 'cross'='grey', 'geometric'='brown', 'both'='black'))+theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #, legend.position="none"
  scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) + guides(fill=guide_legend(title="Median Type")) +
  xlab("Density Observation") + ylab("Method")

temp_melt_shape = temp_melt
variable_shape  = variable

melt_data_total                   = rbind(temp_melt_horizontal_shift, temp_melt_shape)
# melt_data_total$outlier_type = c(rep('Horizontal Shift Outliers', times = nrow(temp_melt_horizontal_shift)),
#                                       rep('Shape Outliers', times = nrow(temp_melt_shape)) )
melt_data_total$outlier_type = c(ifelse(temp_melt_horizontal_shift$Centering=='centered', "Shift, Centered", "Shift, Not Centered"),
                                 ifelse(temp_melt_shape$Centering=='not centered', "Shape, Centered", "Shape, Not Centered"))

## This will likely be the final plot:
melt_data_total$distance = factor(melt_data_total$distance)
melt_data_total$distance = factor(melt_data_total$distance, levels=ordered_distances)
ggplot(melt_data_total, aes(x=variable, y=distance)) +
  geom_tile(aes(fill=as.factor(value))) +  scale_fill_manual(values = c('neither'='white', 
                                                                        'cross'='grey', 'geometric'='brown', 
                                                                        'both'='black'))+theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #, legend.position="none"
  scale_x_discrete(breaks = variable[c(F,F,F,F,T)]) + guides(fill=guide_legend(title="Median Type")) +
  xlab("Density Observation") + ylab("Method") + facet_grid(rows = vars(outlier_type)) + 
  geom_vline(xintercept = 190, color = 'red', linewidth = 1, linetype = 'dashed') +
  scale_y_discrete(limits = rev(levels(melt_data_total$distance)))


################################################
################################################
# I should now calculate Accuracy, TPRs, FPRs:
################################################

# Let's put them all together!
metrics_data              = rbind(horizontal_shift_data_metrics, shape_data_metrics)
# melt_data_total$outlier_type = c(ifelse(horizontal_shift_data_metrics$Centering=='centered', "Horizontal Shift, Centered", "Horizontal Shift, Not Centered"),
#                                  ifelse(shape_data_metrics$Centering=='not centered', "Shape, Centered", "Shape, Not Centered"))

metrics_data$outlier_type = c(rep('Horizontal Shift Outliers', times = nrow(horizontal_shift_data_metrics)),
                              rep('Shape Outliers', times = nrow(shape_data_metrics)) )
metrics_data$distance     = factor(metrics_data$distance)
metrics_data$distance     = factor(metrics_data$distance, levels=ordered_distances)

metrics_data_horizontal_shift = metrics_data[which(metrics_data$outlier_type == "Horizontal Shift Outliers"),]

# g1 = ggplot(metrics_data_horizontal_shift, aes(x = distance, y = accuracy, fill = median_type)) + 
#   theme_bw() + geom_boxplot() + 
#   scale_fill_manual(values = c('cross'='white', 'geometric'='grey'))+
#   ylab("Accuracy") + xlab("Method")
# 
# g2 = ggplot(metrics_data_horizontal_shift, aes(x = distance, y = tpr, fill = median_type)) + 
#   theme_bw() + geom_boxplot() + 
#   scale_fill_manual(values = c('cross'='white', 'geometric'='grey'))+
#   ylab("True Positive Rate") + xlab("Method")
# 
# g3 = ggplot(metrics_data_horizontal_shift, aes(x = distance, y = fpr, fill = median_type)) + 
#   theme_bw() + geom_boxplot() + 
#   scale_fill_manual(values = c('cross'='white', 'geometric'='grey'))+
#   ylab("False Positive Rate") + xlab("Method")
# 

g11 =
  ggplot(metrics_data_horizontal_shift, aes(x = distance, y = accuracy, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white",
                       pattern_angle = 45, pattern_density = 0.15, pattern_spacing = 0.015, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(pattern_fill = "white", fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
    ylab("Accuracy") + xlab("Method") + theme_bw() 

leg = as_ggplot(get_legend(g11))

g1 =
  ggplot(metrics_data_horizontal_shift, aes(x = distance, y = accuracy, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white",
                       pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("Accuracy") + xlab("Method") + theme_bw() + theme(legend.position = "none")

# temp_metrics_data_horizontal_shift           = metrics_data_horizontal_shift
# temp_metrics_data_horizontal_shift$Centering = ifelse(temp_metrics_data_horizontal_shift$Centering=="centered", "        ", "            ")
# 
# g1 =
#   ggplot(temp_metrics_data_horizontal_shift, aes(x = distance, y = accuracy, pattern = Centering, fill = median_type)) +
#   geom_boxplot()+
#   scale_fill_manual(name = "", values = c("white", "grey"), labels = c("","")) +
#   geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", 
#                        pattern_angle = 45, pattern_density = 0.15, pattern_spacing = 0.015, pattern_key_scale_factor = 0.6) +
#   guides(name = "", pattern = guide_legend(name = "", override.aes = list(pattern = "none", fill = "white", color = "white")), 
#          fill = guide_legend(name = "", override.aes = list(pattern = "none", fill = "white", color = "white"))) +
#   ylab("Accuracy") + xlab("Method") + theme_bw()

g2 =
  ggplot(metrics_data_horizontal_shift, aes(x = distance, y = tpr, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", 
                       pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("True Positive Rate") + xlab("Method") + theme_bw() + theme(legend.position = "none")

g3 =
  ggplot(metrics_data_horizontal_shift, aes(x = distance, y = fpr, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white",
                       pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("False Positive Rate") + xlab("Method") + theme_bw() + theme(legend.position = "none")

lst_p = list(g1,g2,g3,leg)
gridExtra::grid.arrange(lst_p[[1]], NULL, lst_p[[2]], lst_p[[4]], lst_p[[3]], NULL, 
                        layout_matrix = matrix(c(1:6), byrow = TRUE, ncol = 2),
                        widths = c(5,0.75))

metrics_data_shape = metrics_data[which(metrics_data$outlier_type == "Shape Outliers"),]
# g1 = ggplot(metrics_data_shape, aes(x = distance, y = accuracy, fill = median_type)) + 
#   theme_bw() + geom_boxplot() + 
#   scale_fill_manual(values = c('cross'='white', 'geometric'='grey'))+
#   ylab("Accuracy") + xlab("Method")
# 
# g2 = ggplot(metrics_data_shape, aes(x = distance, y = tpr, fill = median_type)) + 
#   theme_bw() + geom_boxplot() + 
#   scale_fill_manual(values = c('cross'='white', 'geometric'='grey'))+
#   ylab("True Positive Rate") + xlab("Method")
# 
# g3 = ggplot(metrics_data_shape, aes(x = distance, y = fpr, fill = median_type)) + 
#   theme_bw() + geom_boxplot() + 
#   scale_fill_manual(values = c('cross'='white', 'geometric'='grey'))+
#   ylab("False Positive Rate") + xlab("Method")

g11 =
  ggplot(metrics_data_horizontal_shift, aes(x = distance, y = accuracy, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white",
                       pattern_angle = 45, pattern_density = 0.15, pattern_spacing = 0.015, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(pattern_fill = "white",fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("Accuracy") + xlab("Method") + theme_bw() 

leg = as_ggplot(get_legend(g11))

g1 =
  ggplot(metrics_data_shape, aes(x = distance, y = accuracy, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", 
                       pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("Accuracy") + xlab("Method") + theme_bw() + theme(legend.position = "none") + theme(legend.position = "none")

g2 =
  ggplot(metrics_data_shape, aes(x = distance, y = tpr, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", 
                       pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("True Positive Rate") + xlab("Method") + theme_bw() + theme(legend.position = "none")

g3 =
  ggplot(metrics_data_shape, aes(x = distance, y = fpr, pattern = Centering, fill = median_type)) +
  geom_boxplot()+
  scale_fill_manual(name = "Median Type", values = c("white", "grey")) +
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", pattern_fill = "white", 
                       pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylab("False Positive Rate") + xlab("Method") + theme_bw() + theme(legend.position = "none")

lst_p = list(g1,g2,g3,leg)
gridExtra::grid.arrange(lst_p[[1]], NULL, lst_p[[2]], lst_p[[4]], lst_p[[3]], NULL, 
                        layout_matrix = matrix(c(1:6), byrow = TRUE, ncol = 2),
                        widths = c(5,1))





