#######################################################
#折线图
#######################################################
library("reshape2")
library("ggplot2")

#读入数据
data=read.table("Ave_accuracy_data1_frequency_1-65_100_times_svm_rf.csv",header=T,sep=",")

ggplot(data,aes(x=Sampling_frequency_threshold,y=Accuracy,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Sampling frequency threshold") +
  ylab("Accuracy") 

#保存图片
ggsave(file="Comparision_different_frequency_and_methods.pdf",path="./",width=7,height=5)
