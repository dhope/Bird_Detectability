list2env(readr::read_rds("testingdat.rds"), envir = environment())
Yarray_x = aperm(Yarray, c(2,3,1))
nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
Rcpp::sourceCpp("joing_funs_dh_testing.cpp")
t1 <- nll_fun(inits, X1, X2, tau_params, nsurvey, Yarray_x, tarray, rarray, nrint, ntint, max_r, Ysum, nlimit )
t2 <- nll.fun(inits)
t1
t2
# print(t1,t2)
t1$true/sum(t1$true)
t2[[1]]-t1$CDF_binned
t2[[2]]-t1$tmp1
t2[[3]]-t1$p_matrix

