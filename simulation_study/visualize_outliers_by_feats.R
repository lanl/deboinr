## Script to visualize graph features of the DFNworks curves and color
## these features by whether or not our different methods identify them
## as outliers.
#### Author: Alexander C. Murph
library(GGally)
library(ggplot2)
# Example for deboinr using Justin's breakthrough curves.
setwd("~/DeBoinR/R")
source("summarize_densities.R")
source("helpers.R")
library(reshape)
load(file = "../data/pdf_data.RData")
load(file = "../data/x_grid.RData")

# First grab the feature data:
setwd("~/Documents/UNEs/code_from_justin/loading_graph_stuff")
feat_data = read.csv("features_of_graphs.csv")
setwd("~/DeBoinR_study_files")
feat_data$X = NULL

# Let's start with doing this analysis using the wasserstein metric w the geometric mean.
xx = deboinr(x_grid,
             as.matrix(pdf_data),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = TRUE,
             num_cores = 1
)
print("about to print DeBoinR object...")
print(xx)

# Add the outlier distinctions to the feature data.
feat_data$outlier = as.factor(as.numeric(feat_data$graph_num%in%xx$outliers))

# This plot gave me some good ideas.  Now let's extract the stuff I actually want.
# ggpairs(feat_data, aes(colour = outlier, alpha = 0.4))
gglist = list()
for(col_idx in 2:(ncol(feat_data)-1)){
  temp_graph_data = feat_data[,c(col_idx,ncol(feat_data))]
  
  gglist[[col_idx-1]] = ggplot(temp_graph_data, aes(x = !!sym(names(temp_graph_data)[1] ), fill = outlier ) ) + geom_density(alpha = 0.5)
  
}
pdf("feats_by_outliers/wasserstein_geo_centered_dens.pdf", width = 12, height = 8) # Open a new pdf file
grid.arrange(grobs = gglist, ncol=3, top = "Distribution of Graph Features for Outliers Determined by Wasserstein Dist with Geometric Mean, PDFs Centered")
dev.off()

gglist = list()
for(col_idx in 2:(ncol(feat_data)-1)){
  temp_graph_data = feat_data[,c(col_idx,ncol(feat_data))]
  
  gglist[[col_idx-1]] = ggplot(temp_graph_data, aes(x = !!sym(names(temp_graph_data)[1] ), fill = outlier ) ) + 
    geom_boxplot(alpha = 0.5) + coord_flip()
  
}
pdf("feats_by_outliers/wasserstein_geo_centered_box.pdf", width = 12, height = 8) # Open a new pdf file
grid.arrange(grobs = gglist, ncol=3, top = "Distribution of Graph Features for Outliers Determined by Wasserstein Dist with Geometric Mean, PDFs Centered")
dev.off()



# Let's start with doing this analysis using the wasserstein metric w the geometric mean.
xx = deboinr(x_grid,
             as.matrix(pdf_data),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = FALSE,
             num_cores = 1
)
print("about to print DeBoinR object...")
print(xx)

# Add the outlier distinctions to the feature data.
feat_data$outlier = as.factor(as.numeric(feat_data$graph_num%in%xx$outliers))

# This plot gave me some good ideas.  Now let's extract the stuff I actually want.
# ggpairs(feat_data, aes(colour = outlier, alpha = 0.4))
gglist = list()
for(col_idx in 2:(ncol(feat_data)-1)){
  temp_graph_data = feat_data[,c(col_idx,ncol(feat_data))]
  
  gglist[[col_idx-1]] = ggplot(temp_graph_data, aes(x = !!sym(names(temp_graph_data)[1] ), fill = outlier ) ) + geom_density(alpha = 0.5)
  
}
pdf("feats_by_outliers/wasserstein_geo_not_centered_dens.pdf", width = 12, height = 8) # Open a new pdf file
grid.arrange(grobs = gglist, ncol=3, top = "Distribution of Graph Features for Outliers Determined by Wasserstein Dist with Geometric Mean, PDFs Not Centered")
dev.off()

# This plot gave me some good ideas.  Now let's extract the stuff I actually want.
# ggpairs(feat_data, aes(colour = outlier, alpha = 0.4))
gglist = list()
for(col_idx in 2:(ncol(feat_data)-1)){
  temp_graph_data = feat_data[,c(col_idx,ncol(feat_data))]
  
  gglist[[col_idx-1]] = ggplot(temp_graph_data, aes(x = !!sym(names(temp_graph_data)[1] ), fill = outlier ) ) + 
    geom_boxplot(alpha = 0.5) + coord_flip()
  
}
pdf("feats_by_outliers/wasserstein_geo_not_centered_box.pdf", width = 12, height = 8) # Open a new pdf file
grid.arrange(grobs = gglist, ncol=3, top = "Distribution of Graph Features for Outliers Determined by Wasserstein Dist with Geometric Mean, PDFs Not Centered")
dev.off()




# Let's start with doing this analysis using the wasserstein metric w the geometric mean.
xx = deboinr(x_grid,
             as.matrix(pdf_data),
             distance = "MBD_fboxplot",
             median_type = 'cross',
             center_PDFs = FALSE,
             num_cores = 1
)
print("about to print DeBoinR object...")
print(xx)

# Add the outlier distinctions to the feature data.
feat_data$outlier = as.factor(as.numeric(feat_data$graph_num%in%xx$outliers))

# This plot gave me some good ideas.  Now let's extract the stuff I actually want.
# ggpairs(feat_data, aes(colour = outlier, alpha = 0.4))
gglist = list()
for(col_idx in 2:(ncol(feat_data)-1)){
  temp_graph_data = feat_data[,c(col_idx,ncol(feat_data))]
  
  gglist[[col_idx-1]] = ggplot(temp_graph_data, aes(x = !!sym(names(temp_graph_data)[1] ), fill = outlier ) ) + geom_density(alpha = 0.5)
  
}
pdf("feats_by_outliers/MBD_cross_not_centered_dens.pdf", width = 12, height = 8) # Open a new pdf file
grid.arrange(grobs = gglist, ncol=3, top = "Distribution of Graph Features for Outliers Determined by Modified Band Width with Cross Mean, PDFs Not Centered")
dev.off()

gglist = list()
for(col_idx in 2:(ncol(feat_data)-1)){
  temp_graph_data = feat_data[,c(col_idx,ncol(feat_data))]
  
  gglist[[col_idx-1]] = ggplot(temp_graph_data, aes(x = !!sym(names(temp_graph_data)[1] ), fill = outlier ) ) + 
    geom_boxplot(alpha = 0.5) + coord_flip()
  
}
pdf("feats_by_outliers/MBD_cross_not_centered_box.pdf", width = 12, height = 8) # Open a new pdf file
grid.arrange(grobs = gglist, ncol=3, top = "Distribution of Graph Features for Outliers Determined by Modified Band Width with Cross Mean, PDFs Not Centered")
dev.off()


