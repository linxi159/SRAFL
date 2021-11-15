
test=read.table("symbol_KICH_counts.txt",sep="\t",header=F,check.names=F)
rt=as.matrix(test)
rownames(rt)=rt[,1]

dimnames=list(rownames(rt))
#data=matrix(nrow=nrow(rt),dimnames=dimnames)
data1=matrix(nrow=nrow(rt),ncol=ncol(rt),dimnames=dimnames)
data2=matrix(nrow=nrow(rt),ncol=ncol(rt),dimnames=dimnames)
#data[,1]=rt[,2]
#a = rt[1,2]
#'0' == substring(a, 14, 14)

cnt1=0
cnt2=0
for (i in 2:ncol(rt)){
  a = rt[1,i]
  if ( '0' == substring(a, 14, 14)){#tumor sample
    cnt1=cnt1+1
    data1[,cnt1]=rt[,i]}
  else{#normal sample
    cnt2=cnt2+1
    data2[,cnt2]=rt[,i]
  }
}

exp1=data1[2:nrow(data1),1:cnt1]
dimnames1=list(rownames(exp1),data1[1,1:cnt1])
data11=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)

exp2=data2[2:nrow(data2),1:cnt2]
dimnames2=list(rownames(exp2),data2[1,1:cnt2])
data22=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)

write.table(data11,file="symbol_KICH_counts-T65.txt",sep="\t",quote=F,col.names=T)
write.table(data22,file="symbol_KICH_counts-N24.txt",sep="\t",quote=F,col.names=T)

