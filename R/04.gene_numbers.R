###########################################
#bar plot
##########################################
#读入数据
n_genes=read.table("gene_numbers.txt",header=T,row.names=1,sep="\t")
#加载ggplot2
#BiocManager::install("ggplot2")
library(ggplot2)

n_genes$name=rownames(n_genes)
ggplot(n_genes, aes(x=name, y=n_genes)) + 
  geom_bar(stat="identity")

#改变颜色
ggplot(n_genes, aes(x=name, y=n_genes)) + 
  geom_bar(stat="identity", fill="lightblue", colour="red",width=0.3)

###########################################
#add text
###########################################
ggplot(n_genes, ylabel='a', aes(x=name, y=n_genes,fill=groups)) + 
  geom_bar(stat="identity",width=0.5,aes(fill=groups))+
  geom_text(aes(label=n_genes),vjust=-0.8,colour="black")+
  ylim(0,950)+
  xlab("Land cover classes")+
  ylab("SOC (g C/m2/yr)")



