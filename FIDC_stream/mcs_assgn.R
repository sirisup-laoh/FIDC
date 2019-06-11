mcs_assgn <-function(mcs_cent, mcs_wgh, Adj_mat, x, x_c, x_wgh, tm_stmp, j, rg){
  #created 01/08/2019
  m=1.5
  mcs_no<-dim(mcs_cent)[1]
  x_rep<-matrix(rep(x,mcs_no),byrow=TRUE,nrow=mcs_no)
  dist2x<-sqrt(rowSums((x_rep-mcs_cent)^2))#+0.001
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
    mcs_no<-dim(mcs_cent)[1]
    Adj_mat<-cbind(Adj_mat,rep(0,mcs_no-1))
    Adj_mat<-rbind(Adj_mat,rep(0,mcs_no))
    x_c[j]<-mcs_no
  }
  mcs_assgn<-list(mcs_cent, Adj_mat,x_c,tm_stmp, mcs_wgh)
#} 
}