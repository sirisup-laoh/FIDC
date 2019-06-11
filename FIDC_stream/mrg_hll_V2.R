mrg_hll_V2<-function(np=NaN,x_c,mcs_wgh,Adj_lst){
  #generate cluster
  if (!is.nan(np))
  {  
    p_clst <- which(mcs_wgh>np)
    non_active<-setdiff(1:length(mcs_wgh),p_clst)
    micro_len = length(p_clst)
    Adj_op = Adj_lst
    idx_decend<-order(mcs_wgh[p_clst],decreasing = TRUE)
    clstId = 1; cylclst=list();
    cluster_label = rep(0,micro_len);
    mark_clst = rep(0,length(mcs_wgh));
    for (i in 1:micro_len)
    {
      if (mark_clst[p_clst[idx_decend[i]]]==0)
      {
        seeds = p_clst[idx_decend[i]]; k=1;
        while (k<=length(seeds))
        {
          seedp <- setdiff(Adj_op[[seeds[k]]],non_active)
          seedp = seedp[mark_clst[seedp]==0]           
          seedp <- setdiff(seedp,seeds);
          seeds <- c(seeds, seedp);
          #cannot use union here because it would lose order
          k=k+1;
        }
        
        
        if (length(seeds)>=1)
        {
          cluster_label[seeds] = clstId;
          mark_clst[seeds] <- 1;
          cylclst <- c(cylclst,list(seeds));
          clstId = clstId+1;
        } else
        {
          #skip small cluster
          mark_clst[seeds]<- 1;
        }
      } 
    }
  
  } else
  {
    #No np
    # remove small micro_cluster
    
    p_clst<-1:length(mcs_wgh)
    Adj_op<-Adj_lst[p_clst]
    
    
    #sorting to find the largest weight
    idx_decend<-order(mcs_wgh,decreasing = TRUE)
    clstId = 1; cylclst=list();
    cluster_label = rep(0,length(mcs_wgh));
    micro_len=length(mcs_wgh)
    mark_clst = rep(0,length(mcs_wgh));
    for (i in 1:length(p_clst))
    {
      if (mark_clst[idx_decend[i]]==0)
      {
        seeds <- p_clst[idx_decend[i]]
        seedp <- Adj_op[[seeds]]
        if (length(seedp)>=0) #for removing singleton mcs
        {
          thresh2= 0.13*mcs_wgh[idx_decend[i]]
          #0.13 = 95%, 0.26= 90%
          k=1
          while (k<=length(seeds))
          {
            seedp <- Adj_op[[seeds[k]]]
            seedp = seedp[mark_clst[seedp]==0]           
            seedp = seedp[mcs_wgh[seedp]>thresh2]            
            seedp <- setdiff(seedp,seeds)
            seeds <- c(seeds, seedp)
            #cannot use union here because it would lose order
            k=k+1;
          }#endwhile
        }#endif
        if (length(seeds)>1)
        {
          
          #cluster_label[seeds] = clstId;
          mark_clst[seeds] <- 1;
          cylclst <- c(cylclst,list(seeds));
          clstId = clstId+1;
        } else
        {
          #skip small cluster
          mark_clst[seeds]<- 1;
        } #endelse
      }#endif
      
    }#endfor
    #calculate the size of each tentative cluster and order descendingly
    cylclst_size<-unlist(lapply(cylclst,function(x) sum(mcs_wgh[unlist(x)])))
    idx_decend<-order(cylclst_size,decreasing = TRUE)
    merge_cyl=rep(0,length(cylclst))
    
    #run from large mcs to small mcs
    for (k in length(cylclst):2){
      for (j in (k-1):1){
        #check if two hill connect
        #k is a small hill, j is a larger hill
        #hill1 = larger hill, hill2=smaller hill
        hill1<-cylclst[[idx_decend[j]]]
        hill2<-cylclst[[idx_decend[k]]]
        avg_hill1<-sum(mcs_wgh[hill1])/length(hill1)
        avg_hill2<-sum(mcs_wgh[hill2])/length(hill2)
        lst_cnnct<- c((hill2 %in% unlist(Adj_lst[hill1])),(hill2 %in% unlist(Adj_lst[hill1])))
        #allow merging when two clusters are densly compatible
        flag_cnsst<-(abs(avg_hill1-avg_hill2)<= (0.5*avg_hill1))#0.3 
        if (any(lst_cnnct)){  
          if (flag_cnsst){
            #merge hill2 to hill1 then break
            cylclst[[idx_decend[j]]]<-union(hill2,hill1)
            merge_cyl[idx_decend[k]]<-j
            
            break()
          } else {
            #mark hill2 as border area and remove hill2
            merge_cyl[idx_decend[k]]<- -1
            break()
          }#endif
        } #endif (flag_cnnct)
      }#endfor j
    }#endfor k
    
    #remove merged mcs
    rmv_cyl<-which(merge_cyl!=0)
    rmn_cyl<-setdiff(1:length(cylclst),rmv_cyl)
    cylclst<-cylclst[rmn_cyl]
    
    #remove small clusters (outliers)
    
    ttl_cyl<-length(cylclst); avg_mcs=vector(mode = "numeric", length=ttl_cyl); ttl_mcs=avg_mcs;
    for (i in 1:ttl_cyl)
    {
      ttl_mcs[i]<-sum(mcs_wgh[cylclst[[i]]])
      avg_mcs[i]<-ttl_mcs[i]/length(cylclst[[i]])
      
    }
    dns_value = sqrt(ttl_mcs^2+avg_mcs^2)
    rmn_clst = which(dns_value>0.01*max(dns_value))
    cylclst = cylclst[rmn_clst] 

    #assign cluster label to each mcs
    for (i in 1:length(cylclst))
    {
      cluster_label[cylclst[[i]]]<-i
    }
    
    
    
  }#endelse
  #Reassign points 
  clst_pnt<-rep(0,length(x_c))
  clst_unq<-unique(cluster_label)
  for (i in 1:length(clst_unq)){
    ind <- which(cluster_label==clst_unq[i])  
    for (j in 1:length(ind)){
      clst_pnt[x_c==ind[j]]<-clst_unq[i]
    }
  }
  mrg_hll_V2<-list(cylclst,clst_pnt,cluster_label)
return(mrg_hll_V2)
}