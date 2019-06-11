remove_sngl_st<-function(mcs_cent, mcs_wgh, Adj_mat, tm_stmp, x_c){
  #for streaming data we need tm_stmp as well
  rmv_mcs<-which(mcs_wgh<=1.5)
  rmn_mcs<-setdiff(1:length(mcs_wgh),rmv_mcs)
  p_clst<-rmn_mcs
  Adj_op<-Adj_mat[p_clst,p_clst]
  mcs_wghop <- mcs_wgh[p_clst]                 
  mcs_centop <- mcs_cent[p_clst,]
  tm_stmpop <- tm_stmp[p_clst]
  x_cop <- x_c
  x_cop[which(x_c %in% rmv_mcs)] <- 0
  #rearrange x_c in the remaining mcs
  for (i in 1:length(p_clst)){
    if (p_clst[i]!=i){
      x_cop[x_c==p_clst[i]]<-i
    }
  }
  remove_sngl_st<-list(mcs_centop, mcs_wghop, Adj_op,tm_stmpop, x_cop)
  return(remove_sngl_st)
}