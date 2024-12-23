# Example for deboinr using Justin's breakthrough curves.
# setwd("~/DeBoinR/R")
# source("summarize_densities.R")
# source("helpers.R")
# source("depth_fctns.R")

# pdf_data        = read.csv("../../pdf_example_data.csv")
# pdf_data$X      = NULL
# x_grid          = unlist(pdf_data[1,])
# pdf_data        = pdf_data[2:nrow(pdf_data),]
# save(pdf_data, file = "../data/pdf_data.RData")
# save(x_grid, file = "../data/x_grid.RData")

load(file = "../data/pdf_data.RData")
load(file = "../data/x_grid.RData")

library(DeBoinR)
xx = deboinr(DeBoinR::x_grid,
             as.matrix(DeBoinR::pdf_data),
             distance = "hellinger",
             median_type = 'geometric',
             center_PDFs = TRUE,
             num_cores = 1
             )

print(xx)

xx = deboinr(x_grid,
             as.matrix(pdf_data),
             distance = "hellinger",
             median_type = 'geometric',
             center_PDFs = TRUE,
             num_cores = 1
)

g1 = xx$box_plot[1]

g2 = xx$box_plot[2]

xx = deboinr(x_grid,
       as.matrix(pdf_data),
       distance = "nLQD",
       median_type = 'geometric',
       num_cores = 1
)
# 
# print("about to print DeBoinR object...")
# print(xx)
# g1 = xx$box_plot[[2]]
# 
# 
# xx = deboinr(x_grid,
#              as.matrix(pdf_data),
#              distance = "wasserstein",
#              median_type = 'cross',
#              num_cores = 1
# )

print("about to print DeBoinR object...")
print(xx)
g2 = xx$box_plot[[2]]

g1 = g1 + ggtitle("L2 Distance on the LQD Space using the Geometric Mean")
g2 = g2 + ggtitle("Wasserstein Distance on the PDF Space using the Cross-Section Mean")
lst_p = list(g1,g2)
gridExtra::grid.arrange(lst_p[[1]], lst_p[[2]], layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))



