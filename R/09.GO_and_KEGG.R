library("GOplot")
library("DOSE")
library("clusterProfiler")
library("ggplot2")

#gene_GO2<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','ENPP7P8','OR2T10','SLC47A1','UQCRB','NTS','SLC6A3','CYP4A22','C5orf46')
#gene_GO1<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','MID2','MRPL14','ENPP7P8','ANG',
#           'OR2T10','SLC47A1','UQCRB','MIR9-3HG','NTS','SLC6A3','TRAPPC9','LINC01428','CYP4A22','C5orf46')
#gene_GO3<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','MID2','MRPL14','ENPP7P8','ANG',
#           'OR2T10','SLC47A1','UQCRB','MIR9-3HG','NTS','SLC6A3')
#gene_GO4<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','ENPP7P8',
#           'OR2T10','NTS','SLC6A3')
#gene_GO5<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','ENPP7P8',
#           'OR2T10','NTS','SLC6A3','CYP4A22','C5orf46')
#gene_GO6<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','ENPP7P8','OR2T10','SLC47A1','UQCRB','NTS','SLC6A3','TRAPPC9','LINC01428','CYP4A22','C5orf46')

#gene_GO7<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','MID2','MRPL14','ENPP7P8','ANG',
#           'OR2T10','SLC47A1','UQCRB','MIR9-3HG','NTS','SLC6A3','CYP4A22','C5orf46')

gene_GO<-c("LINC00887",'TTC21B-AS1','SLC47A1P1','SLC10A2','ENPP7P8','OR2T10','SLC47A1','UQCRB','NTS','SLC6A3')

eg = bitr(gene_GO, fromType="SYMBOL", toType="ENTREZID",OrgDb="org.Hs.eg.db", drop = TRUE)

target_gene_id = as.character(eg[,2])
display_number = c(1, 13, 21)
ego_MF <- enrichGO(gene = target_gene_id,
                OrgDb = org.Hs.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
#ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1],]
ego_result_MF <- as.data.frame(ego_MF)[,]

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
#ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_CC <- as.data.frame(ego_CC)[,]

ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
#ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
ego_result_BP <- na.omit(as.data.frame(ego_BP)[,])

category=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
                  rep("molecular function", display_number[3])), 
                levels=c("molecular function", "cellular component", "biological process"))


#KEGG Enrichment
kk<- enrichKEGG(gene= target_gene_id,
                 organism = "hsa",keyType = "kegg", pvalueCutoff = 0.15,
                 pAdjustMethod = "BH",qvalueCutoff = 0.20)

#display_number2 = c(8)
#ego_result_kegg <- na.omit(as.data.frame(kk )[1:display_number2[1], ])
ego_result_kegg <- na.omit(as.data.frame(kk )[,])
    
write.table(ego_result_MF, file="ego_result_MF.csv",sep=",",quote=F,col.names=T,row.names=F)
write.table(ego_result_CC, file="ego_result_CC.csv",sep=",",quote=F,col.names=T,row.names=F)
write.table(ego_result_BP, file="ego_result_BP.csv",sep=",",quote=F,col.names=T,row.names=F)
write.table(ego_result_kegg, file="ego_result_kegg.csv",sep=",",quote=F,col.names=T,row.names=F)



  
