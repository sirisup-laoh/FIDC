source("FIDC_strm_clst.R")
rm(list = setdiff(ls(), lsf.str()))
#KDD 
#load("strm_kdd.Rdata")
#rg=0.5; wnd_len = 10000; lambda = 0.2e-3
#NSL_kdd 
#load("strm_nslkdd.Rdata")
#rg=0.5; wnd_len = 10000; lambda = 0.2e-3
#covr type
#load("strm_covtype.Rdata")
# rg = 0.4; wnd_len = 5000; lambda= 0.4e-3
strm_clst<-FIDC_strm_clst(rg, data_pnt,act_classa, wind_len, lambda)
clst_pnt<-strm_clst[[1]]; 
pf_mat<-strm_clst[[2]]