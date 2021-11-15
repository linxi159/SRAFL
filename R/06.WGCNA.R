library("caret")
library("WGCNA")
library("stringr")

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
#enableWGCNAThreads()
disableWGCNAThreads()


#source("./Rcode_function/function.R")
#Remove abnormal samples for WGCNA
get_abnormal_samples<-function(expro.upper){
  datExpr=as.data.frame(t(expro.upper));
  gsg = goodSamplesGenes(datExpr, verbose = 3);
  gsg$allOK
  sampleTree = hclust(dist(datExpr), method = "average")
  #plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  clust = cutreeStatic(sampleTree, cutHeight =20000, minSize = 10)
  table(clust)
  keepSamples = (clust==1)
  datExpr = datExpr[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  return(datExpr)
}

#Initial setting
set.seed(123)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
corType = "pearson"
type = "unsigned"
robustY = ifelse(corType=="pearson",T,F)


#Read data
trait=read.csv('Clinical.txt',sep = '\t',row.names = 1)
expro.upper=read.csv('merge_RCRPCH_FPKM-T535T289T65.txt',sep = '\t',row.names = 1)
#Data preprocessing
datExpr<-get_abnormal_samples(expro.upper)
rm(expro.upper)

## WGCNA CODE
#[1] source("Rcode_WGCNA\\pickSoftThreshold.R")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#Set the network construction parameter selection range, 
#calculate the scaleless distribution topology matrix

# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresft sholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="scale-free fit index",type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")


#[2] source("Rcode_WGCNA\\cluster.R")
net = blockwiseModules(datExpr,power = sft$powerEstimate, maxBlockSize = 7000,
                       TOMType = "unsigned", deepSplit =2 ,minModuleSize =20,
                       reassignThreshold = 0, mergeCutHeight =0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FPKM-TOM",
                       verbose = 3)
table(net$colors)
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors", 
                                    "GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#[3]source("Rcode_WGCNA\\Eigengene cluster.R")
# module eigengene
MEs = net$MEs 
MEs_col = MEs 
colnames(MEs_col) = paste0("ME", labels2colors( 
  as.numeric(str_replace_all(colnames(MEs),"ME","")))) 
MEs_col = orderMEs(MEs_col) 

#sizeGrWindow(5, 9)
plotEigengeneNetworks(MEs_col, NULL, marDendro = c(0.5,3,3,3), marHeatmap = c(2,2,0.5,2), 
                      plotDendrograms = T, xLabelsAngle = 90)

#[4]source("Rcode_WGCNA\\relation.R")
trait=read.csv('Clinical.txt',sep = '\t',row.names = 1)
# Read in phenotype data
traitData <- read.table(file="Clinical.txt", sep='\t', header=T, row.names=1, 
                        check.names=FALSE, comment='',quote="") 
sampleName = rownames(datExpr) 

library(do)
rownames(traitData) = Replace(data=rownames(traitData),from='-',to='.')

traitData = traitData[match(sampleName, rownames(traitData)), ] 

### Module association with phenotype data
if (corType=="pearsoon") { 
  modTraitCor = cor(MEs_col, traitData, use = "p") 
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
}else { 
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY) 
  modTraitCor = modTraitCorP$bicor 
  modTraitP = modTraitCorP$p}


## Pearson correlation was used for individual columns with zero (or missing) 
## MAD. 
textMatrix = paste(signif(modTraitCor, 2)) 
dim(textMatrix) = dim(modTraitCor)
par(mar=c(2,8,2,2))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), cex.lab = 0.6,ySymbols = colnames(MEs_col), 
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, 
               setStdMargins = FALSE, cex.text =0.6, zlim = c(-1,1))


#[5]source("Rcode_WGCNA\\TOM.R")
TOM = TOMsimilarityFromExpr(datExpr,power = sft$powerEstimate,corType=corType,networkType=type)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
rm(dissTOM)
diag(plotTOM) = 0
TOMplot(plotTOM, net$dendrograms, mergedColors,col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
#TOMplot(plotTOM, net$dendrograms, mergedColors,main = "Network heatmap plot, all genes",
#        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
rm(plotTOM)

#[6]提取指定模块的基因名
# red pink yellow // brown magenta
# Select module
module1 = "red";
module2 = "pink";
module3 = "yellow";
module4 = "brown";
module5 = "magenta";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名

red_inModule = (mergedColors==module1);
red_modProbes = probes[red_inModule]

pink_inModule = (mergedColors==module2);
pink_modProbes = probes[pink_inModule];

yellow_inModule = (mergedColors==module3);
yellow_modProbes = probes[yellow_inModule];

brown_inModule = (mergedColors==module4);
brown_modProbes = probes[brown_inModule];

magenta_inModule = (mergedColors==module5);
magenta_modProbes = probes[magenta_inModule]

#取用模块 red 82 + pink 48，130个基因。
genes = c(red_modProbes,pink_modProbes)
data_ = datExpr[,genes]

Tdata <- read.table(file="Clinical.txt", sep='\t', header=T, row.names=1, 
                    check.names=FALSE, comment='',quote="") 
rownames(data_) = rownames(Tdata)

rt=data_
rt=as.matrix(rt)
dimnames=list(rownames(rt),colnames(rt))
data1=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)

write.table(data1,file="results_RCRPCH_FPKM-T535T289T65-Genes130.txt",sep="\t",quote=F,col.names=T)

genes_red = c(red_modProbes)
genes_pink = c(pink_modProbes)

rt_=as.matrix(genes_red)
rt__=as.matrix(genes_pink)

write.table(rt_,file="Genes_red.txt",sep="\t",quote=F,col.names=T)
write.table(rt__,file="Genes_pink.txt",sep="\t",quote=F,col.names=T)



