# clst_clss = clst_pnt[clst_pnt!=0];
# act_clss = act_clssa[clst_pnt!=0];
# act_unq = unique(act_clss); act_len = length(act_unq);
# clst_unq = unique(clst_clss); clst_len = length(clst_unq);
# C_mat = matrix(0,nrow=act_len,ncol=clst_len);
# for (i in 1:act_len){
#   for (j in 1:clst_len) {
#     C_mat[i,j] = sum(act_clss==act_unq[i]&clst_clss==clst_unq[j]);
#   }
# }
# N = length(act_clss);
# 
# 
# Cdotj = rowSums(C_mat); logCdotj = log10(Cdotj/N);
# Cidot = colSums(C_mat); logCidot = log10(Cidot/N);
# denon = sum(Cdotj*logCdotj)+sum(Cidot*logCidot);
# Cdenon = Cdotj%*%t(Cidot);#Cidot*Cdotj;
# idx <- (C_mat!=0);
# num = -2*sum(C_mat[idx]*log10(C_mat[idx]*N/Cdenon[idx]));
# 
# nmi = num/denon; nmi

clst_idx<-function(clst_pnt,act_clssa){
  
idx_rem<-which(clst_pnt!=0)
act_clss<-act_clssa[idx_rem]; clst_clss<-clst_pnt[idx_rem];
rw_mx<-max(c(act_clss,clst_clss))
c_mat<-matrix(0,ncol=rw_mx,nrow=rw_mx)
sq_len<-length(act_clss)
for (i in 1:sq_len){
  c_mat[act_clss[i],clst_clss[i]]<-c_mat[act_clss[i],clst_clss[i]]+1
}
N<-sq_len
Cidot<-rowSums(c_mat); logCidot = log10(Cidot/N); 
idx1<-!is.infinite(logCidot);
Cdotj<-colSums(c_mat); logCdotj = log10(Cdotj/N);
idx2<-!is.infinite(logCdotj);
denon <- sum(Cdotj[idx2]*logCdotj[idx2])+sum(Cidot[idx1]*logCidot[idx1]);
Cdenon <- Cidot%*%t(Cdotj);
idx <- c_mat!=0;
num <- -2*sum(sum(c_mat[idx]*log10(c_mat[idx]*N/Cdenon[idx])));
nmi <- num/denon; nmi
nis=sum(rowSums(c_mat)^2);
njs=sum(colSums(c_mat)^2);
t1=choose(N,2);
t2=sum(sum(c_mat^2));
t3=.5*(nis+njs);
nc=(N*(N^2+1)-(N+1)*nis-(N+1)*njs+2*(nis*njs)/N)/(2*(N-1));
A=t1+t2-t3;
D=  -t2+t3;
if (t1==nc){
ar=0;  		
} else {
  ar=(A-nc)/(t1-nc);	
}
ri=A/t1
clst_idx<-list(nmi,ar,ri)
}