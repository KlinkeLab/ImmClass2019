##############################################Info#########################################
############# Title:          Sample Selection
############# Description:    This Script is Used to Select Samples from Different Studied
############# Author:         Arezo Torang
############# Date:           13-Nov-2018
###########################################################################################-
####Samples for Immune Cell Types      (B Cell, CD4, CD8, Mono, Neu, NK, DC, M1, M2, MQ)-----
rm(list = ls())
library(dplyr)
library(annotables)
gt = grch38             #human genes table
##Linsley:  B Cell, Mono, Neu, NK, CD4, CD8 :    114 Samples:   45099 Genes---------------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60424
Linsley <- read.table("GSE60424_GEOSubmit_FC1to11_normalized_counts.txt", 
                      header = TRUE, sep = "", dec = ".")
#Selceting samples
Main_Linsley <- data.frame(Linsley$lib311, Linsley$lib297, Linsley$lib290, Linsley$lib228,
                           Linsley$lib310, Linsley$lib296, Linsley$lib289, Linsley$lib227,
                           Linsley$lib309, Linsley$lib295, Linsley$lib288, Linsley$lib226,
                           Linsley$lib314, Linsley$lib300, Linsley$lib293, Linsley$lib231,
                           Linsley$lib312, Linsley$lib298, Linsley$lib291, Linsley$lib229,
                           Linsley$lib313, Linsley$lib299, Linsley$lib292, Linsley$lib230,
                           Linsley$lib232, Linsley$lib233, Linsley$lib234, Linsley$lib235,
                           Linsley$lib236, Linsley$lib237, Linsley$lib238, Linsley$lib239,
                           Linsley$lib240, Linsley$lib241, Linsley$lib242, Linsley$lib243,
                           Linsley$lib244, Linsley$lib245, Linsley$lib246, Linsley$lib247,
                           Linsley$lib248, Linsley$lib249, Linsley$lib250, Linsley$lib251,
                           Linsley$lib252, Linsley$lib253, Linsley$lib254, Linsley$lib255,
                           Linsley$lib256, Linsley$lib257, Linsley$lib258, Linsley$lib260,
                           Linsley$lib261, Linsley$lib262, Linsley$lib263, Linsley$lib264,
                           Linsley$lib266, Linsley$lib267, Linsley$lib268, Linsley$lib269,
                           Linsley$lib270, Linsley$lib271, Linsley$lib273, Linsley$lib274,
                           Linsley$lib275, Linsley$lib276, Linsley$lib277, Linsley$lib278,
                           Linsley$lib280, Linsley$lib281, Linsley$lib282, Linsley$lib283,
                           Linsley$lib284, Linsley$lib285, Linsley$lib302, Linsley$lib303,
                           Linsley$lib304, Linsley$lib305, Linsley$lib306, Linsley$lib307,
                           Linsley$lib316, Linsley$lib317, Linsley$lib318, Linsley$lib319,
                           Linsley$lib320, Linsley$lib321, Linsley$lib323, Linsley$lib324,
                           Linsley$lib325, Linsley$lib326, Linsley$lib327, Linsley$lib329,
                           Linsley$lib330, Linsley$lib331, Linsley$lib332, Linsley$lib333,
                           Linsley$lib334, Linsley$lib336, Linsley$lib337, Linsley$lib338,
                           Linsley$lib339, Linsley$lib340, Linsley$lib342, Linsley$lib343,
                           Linsley$lib344, Linsley$lib345, Linsley$lib346, Linsley$lib347,
                           Linsley$lib349, Linsley$lib350, Linsley$lib351, Linsley$lib352,
                           Linsley$lib353, Linsley$lib354)
#Name Samples
colnames(Main_Linsley) = c("B1_L","B2_L","B3_L","B4_L","Mo1_L","Mo2_L","Mo3_L","Mo4_L","Ne1_L",
                           "Ne2_L","Ne3_L","Ne4_L","Nk1_L","Nk2_L","Nk3_L","Nk4_L","CD4_1_L",
                           "CD4_2_L","CD4_3_L","CD4_4_L","CD8_1_L","CD8_2_L","CD8_3_L","CD8_4_L",
                           "Ne5_L","Mo5_L","B5_L","CD4_5_L","CD8_5_L","Ne6_L","Mo6_L","B6_L",
                           "CD4_6_L","CD8_6_L","Ne7_L","Mo7_L","B7_L","CD4_7_L","CD8_7_L",
                           "Nk5_L","Ne8_L","Mo8_L","B8_L","CD4_8_L","CD8_8_L","Ne9_L","Mo9_L",
                           "B9_L","CD4_9_L","CD8_9_L","Nk6_L","Ne10_L","Mo10_L","B10_L",
                           "CD4_10_L","CD8_10_L","Ne11_L","Mo11_L","B11_L","CD4_11_L","CD8_11_L",
                           "Nk7_L","Ne12_L","Mo12_L","B12_L","CD4_12_L","CD8_12_L","Nk8_L",
                           "Ne13_L","Mo13_L","B13_L","CD4_13_L","CD8_13_L","Nk9_L","Ne14_L",
                           "Mo14_L","B14_L","CD4_14_L","CD8_14_L","Nk10_L","Ne15_L","Mo15_L",
                           "B15_L","CD4_15_L","CD8_15_L","Nk11_L","Ne16_L","Mo16_L","B16_L",
                           "CD4_16_L","CD8_16_L","Ne17_L","Mo17_L","B17_L","CD4_17_L","CD8_17_L",
                           "Nk12_L","Ne18_L","Mo18_L","B18_L","CD4_18_L","CD8_18_L","Ne19_L",
                           "Mo19_L","B19_L","CD4_19_L","CD8_19_L","Nk13_L","Ne20_L","Mo20_L",
                           "B20_L","CD4_20_L","CD8_20_L","Nk14_L")

Main_Linsley<-cbind(Linsley[,1], Main_Linsley[ , order(colnames(Main_Linsley))])
Names<- colnames(Main_Linsley)[-1]

#Add Gene Name with Ensname 
Main_Linsley_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Main_Linsley[,1]),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Linsley_ID)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Linsley_ID)[1],ncol=116))
for (i in 1:dim(Main_Linsley_ID)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,1] %in% Main_Linsley_ID[i,1]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Linsley_ID[i,]
    Main[j,] <-data.frame(Main_Linsley_ID[i,],
                          Main_Linsley[which(Main_Linsley[,1] %in% Main_Linsley_ID[i,1]),2:115])
  } else {print("Repetition")}
}
Main_Linsley <- Main[1:j,]
colnames(Main_Linsley) <-c("Genes","ID",Names)

#write the File
write.table(Main_Linsley, file = "Linsley.txt", sep = "\t", quote = FALSE)



##Hoek:     B Cell, Mono, Neu, NK, mDC, T Cell:  18  Samples:   31042 Genes-------------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64655
Hoek <- read.table("GSE64655_Normalized_transcript_expression_in_human_immune_cells.txt", 
                   header = TRUE, sep = "\t", dec = ".")

#Selceting samples
Main_Hoek <- data.frame(as.numeric(gsub(",", "", Hoek$X.32)),as.numeric(gsub(",", "", Hoek$X.36)),
                        as.numeric(gsub(",", "", Hoek$X.8)),as.numeric(gsub(",", "", Hoek$X.12)),
                        as.numeric(gsub(",", "", Hoek$X.9)),as.numeric(gsub(",", "", Hoek$X.10)),
                        as.numeric(gsub(",", "", Hoek$X.11)),as.numeric(gsub(",", "", Hoek$X.13)),
                        as.numeric(gsub(",", "", Hoek$X.14)), as.numeric(gsub(",", "", Hoek$X.15)),
                        as.numeric(gsub(",", "", Hoek$X.16)),as.numeric(gsub(",", "", Hoek$X.20)),
                        as.numeric(gsub(",", "", Hoek$X.24)),as.numeric(gsub(",", "", Hoek$X.28)),
                        as.numeric(gsub(",", "", Hoek$X.48)),as.numeric(gsub(",", "", Hoek$X.52)),
                        as.numeric(gsub(",", "", Hoek$X.40)),as.numeric(gsub(",", "", Hoek$X.44)))

#Name Samples
colnames(Main_Hoek) = c("B1_H","B2_H","D1_H","D2_H","D3_H","D4_H","D5_H","D6_H","D7_H","D8_H",
                        "Mo1_H","Mo2_H","Ne1_H","Ne2_H","Nk1_H","Nk2_H","T1_H","T2_H")
Main_Hoek<-cbind(Hoek[5:dim(Main_Hoek)[1],c(1,58)], 
                 Main_Hoek[5:dim(Main_Hoek)[1], order(colnames(Main_Hoek))])
colnames(Main_Hoek)[1:2] <- c("Genes","ID")
Names<- colnames(Main_Hoek)

#Add Gene Name with Ensname 
Main_Hoek_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Main_Hoek[,1]),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Hoek_ID)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Hoek_ID)[1],ncol=20))
for (i in 1:dim(Main_Hoek_ID)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,1] %in% Main_Hoek_ID[i,1]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Hoek_ID[i,]
    Main[j,] <-data.frame(Main_Hoek_ID[i,],
                          Main_Hoek[which(Main_Hoek[,1] %in% Main_Hoek_ID[i,1]),3:20])
  } else {print("Repetition")}
}
Main_Hoek <- Main[1:j,]
colnames(Main_Hoek) <-Names

#write the File
write.table(Main_Hoek, file = "Hoek.txt", sep = "\t", quote = FALSE)



##Beyer:    M1, M2:                              6   Samples:   18876 Genes------------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36952
Beyer <- read.table("GSE36952_RPKM.txt", header = TRUE, sep = "\t", dec = ".")

#Selceting samples
Main_Beyer <- Beyer[,1:7]

#Add Gene Name with Ensname 
Main_Beyer_Genes <-data.frame(gt[which(as.matrix(gt[,3]) %in% Main_Beyer$ID),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Beyer_Genes)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Beyer_Genes)[1],ncol=8))
for (i in 1:dim(Main_Beyer_Genes)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,2] %in% Main_Beyer_Genes[i,2]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Beyer_Genes[i,]
    Main[j,] <-data.frame(Main_Beyer_Genes[i,],
                          Main_Beyer[which(Main_Beyer[,1] %in% Main_Beyer_Genes[i,2]),2:7])
  }else{print("Repetition")}
}
Main_Beyer <- Main[1:18876,]
colnames(Main_Beyer) = c("Genes","ID","M1_1_B","M1_2_B","M1_3_B","M2_1_B","M2_2_B","M2_3_B")
write.table(Main_Beyer, file = "Beyer.txt", sep = "\t", quote = FALSE)



##Zhao:     MQ, NK, CD8, CD4:                    14  Samples:   18308 Genes----------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84697
Zhao <- read.table("GSE84697_run_RNASTAR_pipe.txt", header = TRUE, sep = "\t", dec = ".")

#Selceting samples
Main_Zhao <- data.frame(Zhao$GENE,Zhao$s_RCC290_1,Zhao$s_RCC344hg_1,Zhao$s_RCC351_1,
                        Zhao$s_RCC290_2,Zhao$s_RCC343_2,Zhao$s_RCC344cc_1,Zhao$s_RCC344hg_2,
                        Zhao$s_RCC351_2,Zhao$s_RCC343_4,Zhao$s_RCC344cc_3,Zhao$s_RCC344hg_4,
                        Zhao$s_RCC356_1,Zhao$s_RCC344cc_2,Zhao$s_RCC351_4)

#Add Gene Name with Ensname 
Main_Zhao_Genes <-data.frame(gt[which(as.matrix(gt[,3]) %in% Main_Zhao$ID),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Zhao_Genes)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Zhao_Genes)[1],ncol=16))
for (i in 1:dim(Main_Zhao_Genes)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,2] %in% Main_Zhao_Genes[i,2]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Zhao_Genes[i,]
    Main[j,] <-data.frame(Main_Zhao_Genes[i,],
                          Main_Zhao[which(Main_Zhao[,1] %in% Main_Zhao_Genes[i,2]),2:15])
  } else{print("Repetition")}
}
Main_Zhao <- Main[1:18308,]

colnames(Main_Zhao) = c("Genes","ID","CD4_1_Z","CD4_2_Z","CD4_3_Z","CD8_1_Z",
                        "CD8_2_Z","CD8_3_Z","CD8_4_Z","CD8_5_Z","MQ1_Z","MQ2_Z",
                        "MQ3_Z","MQ4_Z","NK1_Z","NK2_Z")
write.table(Main_Zhao, file = "Zhao.txt", sep = "\t", quote = FALSE)



##Daniel:   B Cell, CD4, CD8, NK, Mono:          20  Samples:   22978 Genes----------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246
Daniel <- read.table("GSE74246_RNAseq_All_Counts.txt", header = TRUE, sep = "\t", dec = ".")

Main_Daniel <- data.frame(Daniel$X_TranscriptID,Daniel[,c(41:44,29:36,25:28,37:40)])

Main_Daniel_Genes <-data.frame(gt[which(as.matrix(gt[,3]) %in% Main_Daniel[,1]),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Daniel_Genes)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Daniel_Genes)[1],ncol=22))
for (i in 1:dim(Main_Daniel_Genes)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,2] %in% Main_Daniel_Genes[i,2]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Daniel_Genes[i,]
    Main[j,] <-data.frame(Main_Daniel_Genes[i,],
                          Main_Daniel[which(Main_Daniel[,1] %in% Main_Daniel_Genes[i,2]),2:21])
  }
  else{
    print("Repetition")
  }
}
Main_Daniel <- Main[1:j,]
colnames(Main_Daniel) = c("Genes","ID","B1_D","B2_D","B3_D","B4_D","CD4_1_D","CD4_2_D",
                          "CD4_3_D","CD4_4_D","CD8_1_D","CD8_2_D","CD8_3_D","CD8_4_D",
                          "Mo1_D","Mo2_D","Mo3_D","Mo4_D","Nk1_D","Nk2_D","Nk3_D","Nk4_D")
write.table(Main_Daniel, file = "Daniel.txt", sep = "\t", quote = FALSE)



##Powell:   mDC, Mono:                           6   Samples:   21436 Genes-------------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70106
Powell <- read.table("GSE70106_counts.txt", header = TRUE, sep = "\t", dec = ".")
Main_Powell <-data.frame(Powell$Feature,Powell$hgnc_symbol,Powell$SM_97.CD1c_RP_1_1_CD1c,
                         Powell$SM_28.CD1c_RP_1_1_CD1c,Powell$SM_D36.COK_RP_1_1_CD1c,
                         Powell$SM_97.CD14.Mo_RP_1_1_CD14,Powell$SM_28.CD14.Mo_RP_1_1_CD14,
                         Powell$SM_D36.14.MP_RP_1_1_CD14)
Main_Powell<-cbind(Powell[5:dim(Main_Powell)[1],c(1,58)], 
                   Main_Powell[5:dim(Main_Powell)[1], order(colnames(Main_Powell))])

Main_Powell_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Main_Powell[,1]),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Powell_ID)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Powell_ID)[1],ncol=8))
for (i in 1:dim(Main_Powell_ID)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,1] %in% Main_Powell_ID[i,1]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Powell_ID[i,]
    Main[j,] <-data.frame(Main_Powell_ID[i,],
                          Main_Powell[which(Main_Powell[,1] %in% Main_Powell_ID[i,1]),3:8])
  }
  else{
    print("Repetition")
  }
}
Main_Powell <- Main[1:j,]
colnames(Main_Powell) <- c("Genes","ID","D1_P","D2_P","D3_P","Mo1_P","Mo2_P","Mo3_P")
write.table(Main_Powell, file = "Powell.txt", sep = "\t", quote = FALSE)



##Xue:      MQ, M1, M2:                          9   samples:   11806 Genes----------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55536
Xue <- read.table("GSE55536_HMDM_IPS_IPSDM.norm.expr.txt", header = TRUE, sep = "\t", dec = ".")
Main_Xue <-data.frame(Xue$gene,Xue$HMDM_M1_rep1,Xue$HDMM_M1_rep2,
                      Xue$HMDM_M1_rep3,Xue$HMDM_M2_rep1,Xue$HMDM_M2_rep2,
                      Xue$HMDM_M2_rep3,Xue$HMDM_MAC_rep1,Xue$HMDM_MAC_rep2,Xue$HMDM_MAC_rep3)

Main_Xue_Genes <-data.frame(gt[which(as.matrix(gt[,3]) %in% Main_Xue[,1]),c(1,3)])
j=0
Repeat_free=data.frame(matrix(nrow=dim(Main_Xue_Genes)[1],ncol=2))
Main=data.frame(matrix(nrow=dim(Main_Xue_Genes)[1],ncol=11))
for (i in 1:dim(Main_Xue_Genes)[1]) {
  if( dim(Repeat_free[which(Repeat_free[,2] %in% Main_Xue_Genes[i,2]),])[1]==0){
    j=j+1
    Repeat_free[j,] <-Main_Xue_Genes[i,]
    Main[j,] <-data.frame(Main_Xue_Genes[i,],
                          Main_Xue[which(Main_Xue[,1] %in% Main_Xue_Genes[i,2]),2:10])
  }
  else{
    print("Repetition")
  }
}
Main_Xue <- Main[1:j,]

colnames(Main_Xue) = c("Genes","ID","M1_1_X","M1_2_X","M1_3_X","M2_1_X",
                       "M2_2_X","M2_3_X","MQ1_X","MQ2_X","MQ3_X")
write.table(Main_Xue, file = "Xue.txt", sep = "\t", quote = FALSE)



####Merging Immune Cell Samples------------------
Total_Genes <- rbind(as.matrix(Linsley[,1]),as.matrix(Hoek[,1]),
                     as.matrix(Beyer[,1]),as.matrix(Zhao[,1]),
                     as.matrix(Daniel[,1]),as.matrix(Powell[,1]),
                     as.matrix(Xue[,1]))
#make unique list of genes
Total_Genes <- data.frame(Total_Genes[!duplicated(Total_Genes)]) 
Main_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Total_Genes[,1]),c(1,3)])
Repeat_Free_Index<-seq_along(Main_ID[,1])[!duplicated(Main_ID[,1])]
Total_Genes<-Main_ID[Repeat_Free_Index,]
colnames(Total_Genes) <-c("Genes","ID")
Total <- merge(Total_Genes,Linsley[,-2], all = TRUE)
Total <- merge(Total,Hoek[,-2], all = TRUE)
Total <- merge(Total,Beyer[,-2], all = TRUE)
Total <- merge(Total,Zhao[,-2], all = TRUE)
Total <- merge(Total,Daniel[,-2], all = TRUE)
Total <- merge(Total,Powell[,-2], all = TRUE)
Total <- merge(Total,Xue[,-2], all = TRUE)
write.table(Total, file = "Total.txt", sep = "\t", quote = FALSE)


####Samples for T helper Cells         (Th1, Th2, Th17, Th0, iTreg)-----
##Kanduri:  Th1,Th2:                             6   Samples:   18880 Genes-------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71645
Kanduri <- read.table("GSE71645_Expression_Normalized_DESeq.txt", 
                      header = TRUE, sep = "\t", dec = ".")
Main_Kanduri<-Kanduri[,c(1,5:10)]
colnames(Main_Kanduri)[1]<-"Entrez"
library("org.Hs.eg.db")
Genes_Entrez <- as.list(org.Hs.egENSEMBL2EG)
Genes_Entrez<-data.frame(unlist(Genes_Entrez))
Genes_Entrez <- data.frame(row.names(Genes_Entrez),Genes_Entrez)
colnames(Genes_Entrez)<-c("Genes","Entrez")
Repeat_Free_Index<-seq_along(Genes_Entrez[,2])[!duplicated(Genes_Entrez[,2])]
Genes_Entrez<-Genes_Entrez[Repeat_Free_Index,]
Main_Kanduri_Entrez<-merge(Genes_Entrez,Main_Kanduri)
Main_Kanduri_Entrez<-Main_Kanduri_Entrez[,c(2,1,3:8)]

Main_Kanduri_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Main_Kanduri_Entrez[,1]),c(1,3)])
Repeat_Free_Index<-seq_along(Main_Kanduri_ID[,2])[!duplicated(Main_Kanduri_ID[,2])]
Main_Kanduri_ID<-Main_Kanduri_ID[Repeat_Free_Index,]
colnames(Main_Kanduri_ID)<-c("Genes","ID")
Main_Kanduri<-merge(Main_Kanduri_ID,Main_Kanduri_Entrez)
Main_Kanduri<-Main_Kanduri[,-3]
colnames(Main_Kanduri)<-c("Genes","ID","Th1_1_K","Th1_2_K","Th1_3_K","Th2_1_K","Th2_2_K","Th2_3_K")
write.table(Main_Kanduri, file = "Kanduri.txt", sep = "\t", quote = FALSE)



##Spurlock: Th1,Th2,Th17:                        12  Samples:   46609 Genes------------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66261
Spurlock <- read.table("GSE66261_read_count_table.txt", header = TRUE, sep = "\t", dec = ".")
Main_Spurlock<-Spurlock[,c(1:3,6,9,12,4,7,10,13,5,8,11,14)]
colnames(Main_Spurlock)[1]<-c("Genes")
Main_Spurlock_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Main_Spurlock[,1]),c(1,3)])
Repeat_Free_Index<-seq_along(Main_Spurlock_ID[,2])[!duplicated(Main_Spurlock_ID[,2])]
Main_Spurlock_ID<-Main_Spurlock_ID[Repeat_Free_Index,]
colnames(Main_Spurlock_ID)<-c("Genes","ID")
Main_Spurlock<-merge(Main_Spurlock_ID,Main_Spurlock)
Main_Spurlock<-Main_Spurlock[,-3]
colnames(Main_Spurlock)<-c("Genes","ID","Th1_1_S","Th1_2_S","Th1_3_S","Th1_4_S","Th2_1_S","Th2_2_S",
                           "Th2_3_S","Th2_4_S","Th17_1_S","Th17_2_S","Th17_3_S","Th17_4_S")
write.table(Main_Spurlock, file = "Spurlock.txt", sep = "\t", quote = FALSE)



##Rautio:   Th0,iTreg:                           30  Samples:   35143 Genes------------
####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96538
Rautio <- read.table("GSE96538_rpkm.txt", header = TRUE, sep = "\t", dec = ".")
Main_Rautio<-data.frame(row.names(Rautio),Rautio[,c(7:11,28:32,49:53,17:21,38:42,59:63)])
ENS = data.frame("data"=Main_Rautio[,1])
ENS$data = substr(ENS$data,1,15)
Main_Rautio <- data.frame(ENS,Main_Rautio[2:length(Main_Rautio)])
colnames(Main_Rautio)[1]<-c("Genes")
Main_Rautio_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Main_Rautio[,1]),c(1,3)])
Repeat_Free_Index<-seq_along(Main_Rautio_ID[,2])[!duplicated(Main_Rautio_ID[,2])]
Main_Rautio_ID<-Main_Rautio_ID[Repeat_Free_Index,]
colnames(Main_Rautio_ID)<-c("Genes","ID")
Main_Rautio<-merge(Main_Rautio_ID,Main_Rautio)
colnames(Main_Rautio)<-c("Genes","ID","Th0_1_R","Th0_2_R","Th0_3_R", "Th0_4_R",
                         "Th0_5_R", "Th0_6_R", "Th0_7_R", "Th0_8_R", "Th0_9_R",
                         "Th0_10_R","Th0_11_R","Th0_12_R","Th0_13_R","Th0_14_R",
                         "Th0_15_R","iTreg1_R","iTreg2_R","iTreg3_R","iTreg4_R",
                         "iTreg5_R","iTreg6_R","iTreg7_R","iTreg8_R","iTreg9_R",
                         "iTreg10_R","iTreg11_R","iTreg12_R","iTreg13_R",
                         "iTreg14_R","iTreg15_R")
write.table(Main_Rautio, file = "Rautio.txt", sep = "\t", quote = FALSE)



# ##Revu:     Th0,Th17:                            5   Samples:   26020 Genes-------
# ####https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110097
# Revu <- read.table("GSE110097_FPKM_values_all_Revu.txt", header = TRUE, sep = "\t", dec = ".")
# Main_Revu<-data.frame(Revu[,c(1:3,10,6,7,13)])
# colnames(Main_Revu)<-c("Genes","ID","Th0_1_Re","Th0_2_Re",
#                          "Th17_1_Re","Th17_2_Re","Th17_3_Re")
# Genes=data.frame(gt[which(as.matrix(gt[,3]) %in% Main_Revu[,2]),c(1,3)])
# colnames(Genes)<-c("Genes","ID")
# Main_Revu<-merge(Genes,Main_Revu[,-1])
# Repeat_Free_Index<-seq_along(Main_Revu[,1])[!duplicated(Main_Revu[,1])]
# Main_Revu<-Main_Revu[Repeat_Free_Index,]
# Main_Revu<-Main_Revu[,c(2,1,3:7)]
# write.table(Main_Revu, file = "Revu.txt", sep = "\t", quote = FALSE)



####Merging T helper Cell Samples-------------
Total_Genes <- rbind(as.matrix(Kanduri[,1]),as.matrix(Spurlock[,1]),
                     as.matrix(Rautio[,1]))
Total_Genes <- data.frame(Total_Genes[!duplicated(Total_Genes)])
Main_ID <-data.frame(gt[which(as.matrix(gt[,1]) %in% Total_Genes[,1]),c(1,3)])
Repeat_Free_Index<-seq_along(Main_ID[,1])[!duplicated(Main_ID[,1])]
Total_Genes<-Main_ID[Repeat_Free_Index,]
colnames(Total_Genes) <-c("Genes","ID")
Total_Th <- merge(Total_Genes,Kanduri[,-2], all = TRUE)
Total_Th <- merge(Total_Th,Spurlock[,-2], all = TRUE)
Total_Th <- merge(Total_Th,Rautio[,-2], all = TRUE)

write.table(Total_Th, file = "Total_Th.txt", sep = "\t", quote = FALSE)
