library(R.matlab)
rm(list = setdiff(ls(), lsf.str()))
####load dataset##########
load("data_multiple.Rdata")
rg<- 2.5;
#load("data_imgsgmnt.Rdata")
#rg<- 0.16;
#load("data_satimage.Rdata")
#rg<-0.26
#load("data_pendigit.Rdata")
#rg<-0.34
#load("data_waveform.Rdata")
#rg=0.46
#load("data_waveform_noise.mat")
#rg=0.82
#load('spambased.RData')
#rg<-0.2


source("FIDC_trd_clst.R")
clst_inf <- FIDC_trd_clst(rg, data_pnt=data_pnt)
mcs_cent<-clst_inf[[1]]; mcs_wgh<-clst_inf[[2]]; Adj_lst<-clst_inf[[3]]; x_c<-clst_inf[[4]];


source("MVSA.R")
clst_res <- MVSA(np=NaN,x_c,mcs_wgh,Adj_lst)
cylclst<-clst_res[[1]]; clst_pnt<-clst_res[[2]]; cluster_label<-clst_res[[3]];
source("clst_idx.R")
pfmn_idx<- clst_idx(act_classa,clst_pnt)
nmi<-pfmn_idx[[1]]; ar<-pfmn_idx[[2]]; ri<-pfmn_idx[[3]]
c(nmi,ar,ri)
sum(clst_pnt==0)

 