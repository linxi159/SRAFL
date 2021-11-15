
library('org.Hs.eg.db')
library('clusterProfiler')
library(biomaRt)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
gene<-read.table("sampleExp_id.txt",header=T,sep=",")[,1]
gene
gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_gene_id",values = gene, mart = mart)
outTab=data.frame(gene_name)
write.table(file="gene_name.txt",outTab,sep="\t",quote=F,row.names=F)

#merge两组 修改合并的ID使其一致。
test1=read.table("gene_name.txt",sep="\t",header=T)
test2=read.table("sampleExp.txt",sep="\t",header=T)
test<-merge (test1,test2,by="id")
outTab=data.frame(test)
write.table(file="merge.txt",outTab,sep="\t",quote=F,row.names=F) 

#normalize
library("limma")
rt=read.table("merge.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,3]
exp=rt[,4:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

write.table(data,file="symbol.txt",sep="\t",quote=F,col.names=T)


