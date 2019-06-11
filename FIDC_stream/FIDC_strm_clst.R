FIDC_strm_clst <- function(rg, data_pnt,act_classa, wind_len, lmbd){
source("mcs_assgn.R")
source("remove_sngl_st.R")
source("mrg_hll_V2.R")
source("clst_idx.R")
x_feat<-dim(data_pnt)[2]
x_ttl<-dim(data_pnt)[1]
pf_mat<-matrix(nrow=0,ncol=4)
Adj_mat=matrix(nrow=0,ncol=0)
x_c<-rep(0,x_ttl)
mcs_cent<- matrix(nrow=0,ncol=x_feat)
mcs_wgh <- vector(mode="numeric")
tm_stmp <- vector(mode="numeric")
mntnc_len <- wnd_len/2
data_wgh<-rep(1,dim(data_pnt)[1])
clst_pnt_all=vector()
#create first mcs
x<-data_pnt[1,]
x_wgh<-data_wgh[1]
mcs_wgh<-c(mcs_wgh,x_wgh)
mcs_cent<-rbind(mcs_cent,x)
Adj_mat<-matrix(0,nrow=1,ncol=1)
tm_stmp <- c(tm_stmp,1)
x_c[1]<-1
for (j in 2:x_ttl)
{
  x<-data_pnt[j,]
  x_wgh<-data_wgh[j]
  #   mcs_updt<-mcs_assgn(mcs_cent, mcs_wgh, Adj_mat, x, x_c, x_wgh, tm_stmp, j,rg)
  #   mcs_cent<-mcs_updt[[1]]; Adj_mat<-mcs_updt[[2]]; x_c<-mcs_updt[[3]]; tm_stmp<-mcs_updt[[4]];
  #   mcs_wgh<-mcs_updt[[5]];
  m=1.5
  mcs_no<-dim(mcs_cent)[1]
  x_rep<-matrix(rep(x,mcs_no),byrow=TRUE,nrow=mcs_no)
  dist2x<-sqrt(rowSums((x_rep-mcs_cent)^2))+0.0001
  min_idx<-which.min(dist2x)
  #mcs_in<-dist2x<rg
  mcs_near<-dist2x<2*rg
  mcs_far<-dist2x>5*rg
  
  if (any(mcs_near))
  {
    #distribute load to the relative mcs
    #print("distribute")
    mcs_nridx<-mcs_near
    rel_dist <- dist2x[mcs_nridx]
    mem_den <- 1/sum(rel_dist^(-1/m))
    mem_fnc <-mem_den/(rel_dist^(1/m))
    mcs_wgh_near <- mcs_wgh[mcs_nridx]
    mcs_near_ln <-length(mcs_wgh_near)
    nmrt <- diag(mcs_wgh_near,nrow=mcs_near_ln,ncol=mcs_near_ln)%*%mcs_cent[mcs_nridx,]+
      as.matrix(mem_fnc*x_wgh)%*%t(as.matrix(x))
    dnmrt <- mem_fnc*x_wgh+mcs_wgh_near
    
    mcs_wgh[mcs_near] <- dnmrt
    dnmrt_ln<-length(dnmrt)
    mcs_cent[mcs_near,]<- diag(dnmrt^-1,nrow=dnmrt_ln,ncol=dnmrt_ln)%*%nmrt
    if (any(is.nan(mcs_cent))){
      print(j)
      break
    }
    Adj_mat[mcs_near,mcs_near]<-Adj_mat[mcs_near,mcs_near]+1
    Adj_mat[mcs_near,mcs_far]<-0
    Adj_mat[mcs_far,mcs_near]<-0
    x_c[j]<-min_idx #assign to the closest mcs
    tm_stmp[mcs_near]<-j
  } else
  {
    #generate a new mcs for point x
    mcs_wgh=c(mcs_wgh,x_wgh)
    tm_stmp<-c(tm_stmp,j)
    mcs_cent=rbind(mcs_cent,x)
    if (any(is.nan(mcs_cent))){
      print(j)
      break
    }
    mcs_no<-dim(mcs_cent)[1]
    Adj_mat<-cbind(Adj_mat,rep(0,mcs_no-1))
    Adj_mat<-rbind(Adj_mat,rep(0,mcs_no))
    x_c[j]<-mcs_no
  }
  Main_flg = FALSE
  Clst_flg = FALSE
  #maintainance
  if ((j %% mntnc_len)==0){
    #perform maintenance
    dcrs_fctr <- 2^(-lmbd*(j-tm_stmp))
    mcs_wgh <- mcs_wgh*dcrs_fctr
    #remove mcs with weight less than 1.5
    crnt_wn_idx = (1+ ((j %/% wnd_len)*wnd_len)):j
    Main_flg = TRUE
    #     clst_inf = remove_sngl_st(mcs_cent, mcs_wgh, Adj_mat, tm_stmp, x_c[crnt_wn_idx])
    #     mcs_cent<-clst_inf[[1]]; mcs_wgh<-clst_inf[[2]]; Adj_mat<-clst_inf[[3]]; 
    #     tm_stmp<-clst_inf[[4]]; x_c[crnt_wn_idx]<-clst_inf[[5]];  
    clst_inf = remove_sngl_st(mcs_cent, mcs_wgh, Adj_mat, tm_stmp, x_c)
    mcs_cent<-clst_inf[[1]]; mcs_wgh<-clst_inf[[2]]; Adj_mat<-clst_inf[[3]]; 
    tm_stmp<-clst_inf[[4]]; x_c<-clst_inf[[5]];    
  }#if (j %% mntnc_len)
  #measure cluster quality
  
  
  if ((j %% wnd_len)==0){
    #perform final clustering
    #convert to adj_list
    Adj_lst = list()
    for (i in 1:dim(Adj_mat)[1]){
      Adj_lst[[i]]=which(Adj_mat[i,]>0)
    }
    crnt_wn_idx = (1+((j %/% wnd_len)-1)*wnd_len):j
    Clst_flg = TRUE
    clst_res <- mrg_hll_V2(np=NaN,x_c[crnt_wn_idx],mcs_wgh,Adj_lst)
    cylclst <- clst_res[[1]]; clst_pnt<-clst_res[[2]]; cluster_label<-clst_res[[3]];
    #measure clustering performance
    clst_pfm <- clst_idx(clst_pnt,act_classa[crnt_wn_idx])
    nmi<-clst_pfm[[1]]; ar<-clst_pfm[[2]]; ri<-clst_pfm[[3]];
    pf_mat<-rbind(pf_mat,c(j,nmi,ar,ri))
    clst_pnt_all=c(clst_pnt_all,clst_pnt)
  }
}#for (j in 2:x_ttl)
FIDC_strm_clst<-list(clst_pnt_all,pf_mat)
return(FIDC_strm_clst)
}