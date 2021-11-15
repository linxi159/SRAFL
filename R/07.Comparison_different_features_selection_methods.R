#######################################################
#折线图
#######################################################
library("reshape2")
library("ggplot2")

#读入数据
data=read.table("Ave_Accuracy_chi2_RFE_f1-20_svm_rf_100_times.csv",header=T,sep=",")

ggplot(data,aes(x=Gene_number,y=Accuracy,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Gene number") +
  ylab("Accuracy") 

#保存图片
ggsave(file="Comparision_dif_fea_sel_methods.pdf",path="./",width=7,height=5)
