
setwd("~/r_workplace/kidney_cancer/data/exp_merge_selection_FPKM")

test1=read.table("merge_RCRPCH_FPKM-T535T289T65.txt",sep="\t",header=T,check.names=F)

rt=as.matrix(test1)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
exp1=t(exp)
dimnames=list(rownames(exp1),colnames(exp1))
data=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)

write.table(data,file="expression.txt",sep="\t",quote=F,col.names=T)

