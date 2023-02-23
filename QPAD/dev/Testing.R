# list2env(readr::read_rds("dev/testingdat.rds"), envir = environment())
load("dev/BreakingBad.RData")
Yarray_x = aperm(Yarray, c(2,3,1))
nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
# Rcpp::sourceCpp("joing_funs_dh_testing.cpp")
library(IlesShowUSomething)

# runs <- matrix(runif(1e4, -1,1),ncol = 5)
# readr::write_rds(runs, 'dev/runs.rds')
runs <- readr::read_rds("dev/runs.rds")

# t1 <- nll_fun(inits, X1, X2, tau_params, nsurvey, Yarray_x, tarray, rarray, nrint, ntint, max_r, Ysum, nlimit )
# t2 <- nll.fun(inits)
# t1
# t2
# # print(t1,t2)
# t1$true/sum(t1$true)
# t2[[1]]-t1$CDF_binned
# t2[[2]]-t1$tmp1
# t2[[3]]-t1$p_matrix
# 
# for(k in 1:nsurvey){
#   Y <- Yarray[k,,]
#   nrow(Y) == nrint[k]
# }

# Yarraylist <- vector("list", nsurvey)
# for(k in 1:nsurvey){
#   Yarraylist[[k]] <- as.matrix(Yarray[k, 1:nrint[k],1:ntint[k] ])
# }  

vec <- c(Yarray) 
# vecno_NA <- vec[!is.na(vec)]
# arr2 <- replace(Yarray_x, TRUE, vec + 100)

# vec[addr(0,1,0, nsurvey, max(nrint) )]
# Yarray[1,2,1]
# out <- TRUE
# for(k in 1:(nsurvey)){
#   for(i in 1:(nrint[k])){
#   for(j in 1:(ntint[k])){
#    cat(vec[addr(k-1,i-1,j-1, nsurvey, max(nrint) )])
#     print(Yarray[k,i,j])
#   }
#   }
# }
# vec[addr(3, 0,0, nsurvey, max(nrint) )+1]
# Yarray[4,1,1]
nll_fun(runs[1,], X1, X2, tau_params, nsurvey, 
        Yarray_x, tarray, rarray, nrint,
        ntint, max_r, Ysum, nlimit )

out <- purrr::map_dbl(1:nrow(runs), 
               ~{print(.x);print(nll_fun(runs[.x,], X1, X2, tau_params, nsurvey, 
                                   Yarray_x, tarray, rarray, nrint,
                        ntint, max_r, Ysum, nlimit ))}
                 )


# 
runs <- readr::read_rds("dev/runs.rds")
runs[345,]
k <- 6
Yarray_x[,,k]
Ysum[k]
tau_params
ntint[k]
tarray[k,]
nrint[k]
rarray[k,]
max_r[k]
tau_est <- X1[k,] %*% runs[32,1:2]
phi_est <- X2[k,] %*% runs[32,3:5]
# # 
# # for(i in 1:100){
# nll_fun(runs[33,], X1, X2, tau_params, nsurvey,
#                   Yarray_x, tarray, rarray, nrint,
#                   ntint, max_r, Ysum, nlimit )
# }


