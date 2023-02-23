load("results/CurrentEnv.RData")
print("Data loaded")
library(IlesShowUSomething)
# source("cmulti_fit_joint_cpp.R")
# bootsamps <- read_rds("results/bootstrapsample_14.rds")
# bootsamps <- readr::read_rds("results/bootstrapsample_1.rds")

# for (b in 1:bootreps){

run_bootstrap <- function(b){
  print(b)
  bootsamps <- sample(1:nsurvey,nsurvey,replace=TRUE)
  Yboot <- Yarray[bootsamps,,]
  rboot <- rarray[bootsamps,]
  tboot <- tarray[bootsamps,]
  X1boot <- X1[bootsamps,]
  X2boot <- X2[bootsamps,]
  start <- Sys.time()
  fit <- cmulti_fit_joint(Yboot,
                               rboot,
                               tboot,
                               X1 = X1boot, # Design matrix for tau
                               X2 = X2boot  # Design matrix for phi
  )
  end <- Sys.time()
  print(end-start)
  return(fit)
  
  
}

fits <- purrr::map(1:bootreps, run_bootstrap)

readr::write_rds(fits, "output/fits_bootstraps2.rds")
