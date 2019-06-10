FIDC_trd_clst <- function(rg, data_pnt){
  #created 01/08/2019
  m <-1.5
  rg<-rg/2
  x_feat<-dim(data_pnt)[2]
  x_ttl<-dim(data_pnt)[1]
  data_wgh=rep(1,dim(data_pnt)[1])
  Adj_mat=matrix(nrow=0,ncol=0)
  x_c<-rep(0,x_ttl)
  mcs_cent<-matrix(nrow=0,ncol=x_feat)
  mcs_wgh<-vector(mode="numeric")
  for(j in 1:x_ttl)
  {
    x<-data_pnt[j,]
    x_wgh<-data_wgh[j]
    mcs_no<-dim(mcs_cent)[1]
    if (mcs_no>0)
    {  
      x_rep<-matrix(rep(x,mcs_no),byrow=TRUE,nrow=mcs_no)
      dist2x<-sqrt(rowSums((x_rep-mcs_cent)^2))+0.001
      min_idx<-which.min(dist2x)
      mcs_in<-dist2x<rg
      mcs_near<-dist2x<2*rg
      mcs_far<-dist2x>5*rg
      # if (any(mcs_in))
      # {
      #   #place inside a micro-cluster
      #   min_cent<-mcs_cent[min_idx,]
      #   min_cent<-(x+min_cent*mcs_wgh[min_idx])/(mcs_wgh[min_idx]+x_wgh)
      #   mcs_wgh[min_idx]<-mcs_wgh[min_idx]+x_wgh
      #   mcs_cent[min_idx,]<-min_cent
      #   Adj_mat[mcs_near,mcs_near]<-Adj_mat[mcs_near,mcs_near]+1
      #   x_c[j]<-min_idx #assign to the closest mcs
      # } else 
      
      if (any(mcs_near))
      {
        #distribute load to the relative mcs
        #print("distribute")
        rel_dist <- dist2x[mcs_near]
        mem_den <- 1/sum(rel_dist^(-1/m))
        mem_fnc <-mem_den/(rel_dist^(1/m))
        mcs_wgh_near <- mcs_wgh[mcs_near]
        mcs_near_ln <-length(mcs_wgh_near)
        nmrt <- diag(mcs_wgh_near,nrow=mcs_near_ln,ncol=mcs_near_ln)%*%mcs_cent[mcs_near,]+
          as.matrix(mem_fnc*x_wgh)%*%t(as.matrix(x))
        dnmrt <- mem_fnc*x_wgh+mcs_wgh_near
        mcs_wgh[mcs_near] <- dnmrt
        dnmrt_ln<-length(dnmrt)
        mcs_cent[mcs_near,]<- diag(dnmrt^-1,nrow=dnmrt_ln,ncol=dnmrt_ln)%*%nmrt
        Adj_mat[mcs_near,mcs_near]<-Adj_mat[mcs_near,mcs_near]+1
        Adj_mat[mcs_near,mcs_far]<-0
        Adj_mat[mcs_far,mcs_near]<-0
        x_c[j]<-min_idx #assign to the closest mcs
      } else
      {
        #generate a new mcs for point x
        mcs_wgh=c(mcs_wgh,x_wgh)
        mcs_cent=rbind(mcs_cent,x)
        mcs_no<-dim(mcs_cent)[1]
        Adj_mat<-cbind(Adj_mat,rep(0,mcs_no-1))
        Adj_mat<-rbind(Adj_mat,rep(0,mcs_no))
        x_c[j]<-mcs_no
        #print("genearate")
      }
    } 
    else
    {
      #generate the 1st mcs for point x
      mcs_wgh=c(mcs_wgh,x_wgh)
      mcs_cent=rbind(mcs_cent,x)
      Adj_mat<-matrix(0,nrow=1,ncol=1)
      x_c[j]<-dim(mcs_cent)[1]
    }
  }#for(j in 1:x_ttl)
  #fuzzy_mcsfd<-list(mcs_cent,mcs_wgh,Adj_mat,x_c)
  #Remove_sngl
  rmv_mcs<-which(mcs_wgh<=1.5)
  rmn_mcs<-setdiff(1:length(mcs_wgh),rmv_mcs)
  p_clst<-rmn_mcs
  Adj_op<-Adj_mat[p_clst,p_clst]
  mcs_wghop <- mcs_wgh[p_clst]                 
  mcs_centop <- mcs_cent[p_clst,]  
  x_cop <- x_c
  x_cop[which(x_c %in% rmv_mcs)] <- 0
  #rearrange x_c in the remaining mcs
  for (i in 1:length(p_clst)){
    if (p_clst[i]!=i){
      x_cop[x_c==p_clst[i]]<-i
    }
  }
  mcs_cent<-mcs_centop; mcs_wgh<-mcs_wghop; Adj_mat<-Adj_op; x_c<-x_cop;
  #convert to adj_list
  Adj_lst = list()
  for (i in 1:dim(Adj_mat)[1]){
    Adj_lst[[i]]=which(Adj_mat[i,]>0)
  }
  FIDC_trd_clst<-list(mcs_cent,mcs_wgh,Adj_lst,x_c)
  }