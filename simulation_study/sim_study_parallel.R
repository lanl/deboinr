# File to run 100 instances of the simulation study in parallel.

sim_study - function(x){
  system(paste("Rscript deboinr_simulation_study.R", x))
}