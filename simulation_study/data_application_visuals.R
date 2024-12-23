# Example for deboinr using Justin's breakthrough curves.
setwd("~/DeBoinR/R")
source("summarize_densities.R")
source("helpers.R")
library(reshape)

# pdf_data        = read.csv("../../pdf_example_data.csv")
# pdf_data$X      = NULL
# x_grid          = unlist(pdf_data[1,])
# pdf_data        = pdf_data[2:nrow(pdf_data),]
# save(pdf_data, file = "../data/pdf_data.RData")
# save(x_grid, file = "../data/x_grid.RData")

load(file = "../data/pdf_data.RData")
load(file = "../data/x_grid.RData")

# xx = deboinr(x_grid,
#              as.matrix(pdf_data),
#              distance = "hellinger",
#              median_type = 'geometric',
#              center_PDFs = TRUE,
#              num_cores = 1
# )

xx = deboinr(x_grid,
       as.matrix(pdf_data),
       distance = "wasserstein",
       median_type = 'geometric',
       center_PDFs = TRUE,
       num_cores = 1
)
print("about to print DeBoinR object...")
print(xx)
g1 = xx$box_plot[[2]]

idx_list                    = 1:nrow(pdf_data)
g1_outliers_data            = data.frame(t(pdf_data))
g1_outliers_data            = g1_outliers_data[,c(which(!(idx_list%in%xx$outliers)),which(idx_list%in%xx$outliers))]
g1_outliers_data            = melt(g1_outliers_data)
g1_outliers_data$is_outlier = rep( c(rep(0, times = sum(!(idx_list%in%xx$outliers))), rep(1, times = sum(idx_list%in%xx$outliers))), each = length(x_grid) )
g1_outliers_data$time       = rep(x_grid, times = nrow(pdf_data))

g1_outliers = ggplot(g1_outliers_data, aes(x = time, y = value, group = as.factor(variable), color = as.character(is_outlier))) + 
                      geom_line() + 
                      scale_color_manual(values = c('0'='grey', '1'='black'))+theme_bw()+ theme(legend.position = 'none')

xx = deboinr(x_grid,
             as.matrix(pdf_data),
             distance = "wasserstein",
             median_type = 'geometric',
             center_PDFs = FALSE,
             num_cores = 1
)
print("about to print DeBoinR object...")
print(xx)
g2 = xx$box_plot[[2]]

idx_list                    = 1:nrow(pdf_data)
g2_outliers_data            = data.frame(t(pdf_data))
g2_outliers_data            = g2_outliers_data[,c(which(!(idx_list%in%xx$outliers)),which(idx_list%in%xx$outliers))]
g2_outliers_data            = melt(g2_outliers_data)
g2_outliers_data$is_outlier = rep( c(rep(0, times = sum(!(idx_list%in%xx$outliers))), rep(1, times = sum(idx_list%in%xx$outliers))), each = length(x_grid) )
g2_outliers_data$time       = rep(x_grid, times = nrow(pdf_data))

g2_outliers = ggplot(g2_outliers_data, aes(x = time, y = value, group = as.factor(variable), color = as.character(is_outlier))) + 
  geom_line() + 
  scale_color_manual(values = c('0'='grey', '1'='black'))+theme_bw()+ theme(legend.position = 'none')



xx = deboinr(x_grid,
             as.matrix(pdf_data),
             distance = "MBD_fboxplot",
             median_type = 'cross',
             center_PDFs = FALSE,
             num_cores = 1
)
print("about to print DeBoinR object...")
print(xx)
g3 = xx$box_plot[[2]]

idx_list                    = 1:nrow(pdf_data)
g3_outliers_data            = data.frame(t(pdf_data))
g3_outliers_data            = g3_outliers_data[,c(which(!(idx_list%in%xx$outliers)),which(idx_list%in%xx$outliers))]
g3_outliers_data            = melt(g3_outliers_data)
g3_outliers_data$is_outlier = rep( c(rep(0, times = sum(!(idx_list%in%xx$outliers))), rep(1, times = sum(idx_list%in%xx$outliers))), each = length(x_grid) )
g3_outliers_data$time       = rep(x_grid, times = nrow(pdf_data))

g3_outliers = ggplot(g3_outliers_data, aes(x = time, y = value, group = as.factor(variable), color = as.character(is_outlier))) + 
  geom_line() + 
  scale_color_manual(values = c('0'='grey', '1'='black'))+theme_bw()+ theme(legend.position = 'none')




g1          = g1 + ggtitle("") + xlab("log breakthrough time")
g2          = g2 + ggtitle("") + xlab("log breakthrough time")
g3          = g3 + ggtitle("") + xlab("log breakthrough time")
g1_outliers = g1_outliers + ggtitle("") + xlab("log breakthrough time") + ylab("density")
g2_outliers = g2_outliers + ggtitle("") + xlab("log breakthrough time") + ylab("density")
g3_outliers = g3_outliers + ggtitle("") + xlab("log breakthrough time") + ylab("density")
lst_p = list(g2,g2_outliers, g1,g1_outliers ,g3,g3_outliers)
# gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], lst_p[[3]], lst_p[[4]], lst_p[[5]], lst_p[[6]], 
#                         layout_matrix = matrix(c(1,2,3,4,5,6), byrow = TRUE, ncol = 2))

grid.arrange(arrangeGrob(lst_p[[1]], lst_p[[2]], ncol = 2, top = 'Wasserstein Distance on the PDF Space using the Geometric Median without Centering'), 
             arrangeGrob(lst_p[[3]], lst_p[[4]], ncol = 2, top = 'Wasserstein Distance on the PDF Space using the Geometric Median with Centering'),
             arrangeGrob(lst_p[[5]], lst_p[[6]], ncol = 2, top = 'Modified Band Depth on the PDF Space using the Cross-Section Median without Centering'), 
             ncol=1)

