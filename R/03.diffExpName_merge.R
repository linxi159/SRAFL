
setwd("~/r_workplace/kidney_cancer/data/diffExpName_merge")

rt1=read.table("diffExpName_KICH.txt",sep="\t",header=T,check.names=F) 
rt2=read.table("diffExpName_KIRC.txt",sep="\t",header=T,check.names=F) 
rt3=read.table("diffExpName_KIRP.txt",sep="\t",header=T,check.names=F) 
a=rt1[,1]
b=rt2[,1]
c=rt3[,1]

re1=intersect(x=a,y=b)
re2=intersect(x=re1,y=c)
write.table(re2,file="diffExpName_intersect.txt",sep="\t",quote=F,col.names=F)

re3=union(x=a,y=b)
re4=union(x=re3,y=c)
write.table(re4,file="diffExpName_uniont.txt",sep="\t",quote=F,col.names=F)



