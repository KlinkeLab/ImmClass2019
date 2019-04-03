##############################################Info#########################################
############# Title:          Data Preparation and Normalization
############# Description:    This Script is Used to Prepare Data and Normalize Data
############# Author:         Arezo Torang
############# Date:           13-Nov-2018
###########################################################################################-




rm(list = ls())

###Importing Data-----------------------------------------------------------




Total <- read.table("Total.txt", sep = "\t", quote = "",
                    stringsAsFactors=FALSE)

Total_Th <- read.table("Total_Th.txt", sep = "\t", quote = "",
                       stringsAsFactors=FALSE)




####Keeping Genes That Have >5 Values---------------------------------------
####for at Least 5 Samples in Immune Cells




ZeroGenes=data.frame()
Main=data.frame()

Total_0<-Total
Total_0[is.na(Total_0)] <- 0


for (i in 1:dim(Total)[1]) {
  
 if (sum(sum(Total_0[i,-c(1,2)]>=5))>4) {
   Main <-rbind(Main,i)
   
 } else {
   ZeroGenes <- rbind(ZeroGenes,i)
 }
}


MainData <-Total[unlist(Main),]
write.table(MainData, file = "MainData.txt", 
            sep = "\t", quote = FALSE)




####Keep Genes That Have >5 Values------------------------------------------
####for at Least 4 Samples in T helper Cells




ZeroGenes=data.frame()
Main=data.frame()

Total_0<-Total_Th
Total_0[is.na(Total_0)] <- 0


for (i in 1:dim(Total_Th)[1]) {
  
  if (sum(sum(Total_0[i,-c(1,2)]>=5))>3) {
    Main <-rbind(Main,i)
    
  } else {
    ZeroGenes <- rbind(ZeroGenes,i)
  }
}


MainData_Th <-Total_Th[unlist(Main),]
write.table(MainData_Th, file = "MainData_Th.txt", 
            sep = "\t", quote = FALSE)




####Separation of Testing & Training Samples--------------------------------
##(Train:108,Test:77)(Train_Th:31,Test_Th:17)




##Immune Cells
TrainData<-MainData[,c(1:12,23:32,43:52,63:72,83:92,103:112,
                       117, 119:124, 127, 129, 131, 135:136, 
                       138:139, 141:142,144:146,149:151,153, 
                       155 : 157, 159:161, 163:165, 167:169,
                       171 : 173, 175:176, 178:179, 181:182,
                       184 : 185, 187, 188)]


#Training Data For Immune Cells
write.table(TrainData,"TrainData.txt") 

TestData<-MainData[,-c(3:12,23:32,43:52,63:72,83:92,103:112,
                       117, 119:124, 127, 129, 131, 135:136, 
                       138:139, 141:142, 144:146, 149 : 151,
                       153, 155 : 157, 159 : 161, 163 : 165, 
                       167 : 169, 171:173, 175:176, 178:179,
                       181:182, 184 : 185, 187, 188)]




##T helper Cells
TrainData_Th<-MainData_Th[,c(1:4,6:7,9:11,13:15,17:19,23:25,
                             28:30,33:35,38:40,43:45,48:50)]


#Training Data For T helper Cells
write.table(TrainData_Th,"TrainData_Th.txt") 

TestData_Th<-MainData_Th[,-c(3:4,6:7,9:11,13:15,17:19,23:25,
                             28:30,33:35,38:40,43:45,48:50)]




####Normalization----------------------------------------------------------




####Transcript Length & GC content for Immune Cell Samples 
####(Using GC Content for Normalization)
####https://git.embl.de/rogon/predoc_2017/blob/master/Practical/mart_export.txt




Gene_Lenght <- read.table("mart_export.txt", header = TRUE, 
                          sep = "\t", dec = ".")

G1 <- Gene_Lenght[which(Gene_Lenght[,1] %in% MainData[,1]),]
Repeat_Free <- data.frame(matrix(0,nrow=dim(G1)[1],ncol=4))
R <- data.frame(matrix(0,nrow=dim(G1)[1],ncol=1))
j=0


for (i in 1:dim(G1)[1]) {
  
 if (sum(which(R[,1] %in% as.character(G1[i,1])))==0){
   j <- j+1; print(i)
   R[j,1] <- as.character(G1[i,1])
   m <- which(G1[,1] %in% G1[i,1])
   
   if (length(m) > 1) {
       a <- i-1+which.max(G1[m, 3])
       Repeat_Free[j,] <- G1[a,]
       
   } else {Repeat_Free[j,] <- G1[i,]}
 }
}


Repeat_Free <- data.frame(R[1:j,1],Repeat_Free[1:j,3:4])
colnames(Repeat_Free) <- c("Genes", "Transcript_Length", "GC_Content")
Repeat_Free <- Repeat_Free[order(Repeat_Free[,1]),]

write.table(Repeat_Free, file = "Length_GC.txt", 
            sep = "\t", quote = FALSE)


Length_GC <- read.table("Length_GC.txt", sep = "\t", quote = "")
Data_Length_GC <- merge(MainData,Length_GC)
TrainData_Length_GC <- merge(TrainData,Length_GC)
TestData_Length_GC <- merge(TestData,Length_GC)




####Transcript Length & GC content For T helper Cell Samples



Gene_Lenght <- read.table("mart_export.txt", header = TRUE,
                          sep = "\t", dec = ".")

G1 <- Gene_Lenght[which(Gene_Lenght[,1] %in% MainData_Th[,1]),]
Repeat_Free <- data.frame(matrix(0,nrow=dim(G1)[1],ncol=4))
R <- data.frame(matrix(0,nrow=dim(G1)[1],ncol=1))
j=0


for (i in 1:dim(G1)[1]) {
  
  if (sum(which(R[,1] %in% as.character(G1[i,1])))==0){
    j <- j+1
    print(i)
    R[j,1] <- as.character(G1[i,1])
    m <- which(G1[,1] %in% G1[i,1])
    
    if (length(m) > 1) {
      a <- i-1+which.max(G1[m, 3])
      Repeat_Free[j,] <- G1[a,]
      
    } else {Repeat_Free[j,] <- G1[i,]}
  }
}


Repeat_Free <- data.frame(R[1:j,1],Repeat_Free[1:j,3:4])
colnames(Repeat_Free) <- c("Genes", "Transcript_Length", "GC_Content")
Repeat_Free <- Repeat_Free[order(Repeat_Free[,1]),]

write.table(Repeat_Free, file = "Length_GC_Th.txt", 
            sep = "\t", quote = FALSE)


Length_GC_Th <- read.table("Length_GC_Th.txt", sep = "\t", quote = "")
Data_Th_Length_GC <- merge(MainData_Th,Length_GC_Th)
TrainData_Th_Length_GC <- merge(TrainData_Th,Length_GC_Th)
TestData_Th_Length_GC <- merge(TestData_Th,Length_GC_Th)




####Normalization Using EDASeq

library(EDASeq)



###Training Samples for Immune Cells



boxplot(as.matrix(log2(TrainData_Length_GC[,3:110]+1)))

dataWithin <- withinLaneNormalization(x=as.matrix(TrainData_Length_GC[,3:110]),
                                    y=TrainData_Length_GC[,112], offset=FALSE,
                                    round=TRUE, which = "full")

boxplot(log2(dataWithin+1))

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

boxplot(log2(dataNorm+1))

Normalized_TrainData <- data.frame(TrainData_Length_GC[,1:2],dataNorm)




####Deleting Zeros after Normalization


ZeroGenes=data.frame()
Main=data.frame()

Total_0<-Normalized_TrainData
Total_0[is.na(Total_0)] <- 0

for (i in 1:dim(Total_0)[1]) {
  
 if (sum(sum(Total_0[i,3:dim(Total_0)[2]]>=5))>4) {
   Main <-rbind(Main,i)
   print(i)
   
 } else {ZeroGenes <- rbind(ZeroGenes,i)
   print("Zero")}
}

Normalized_TrainData <-Normalized_TrainData[unlist(Main),]

write.table(Normalized_TrainData, file = "Normalized_TrainData.txt", 
            sep = "\t", quote = FALSE)



####Testing Samples for Immune Cells



boxplot(as.matrix(log2(TestData_Length_GC[,3:81]+1)))

dataWithin <- withinLaneNormalization(x=as.matrix(TestData_Length_GC[,3:81]),
                                      y=TestData_Length_GC[,83], offset=FALSE,
                                      round=TRUE, which = "full")

boxplot(log2(dataWithin+1))

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

boxplot(log2(dataNorm+1))

Normalized_TestData <- data.frame(TestData_Length_GC[,1:2],dataNorm)
Normalized_TestData <-Normalized_TestData[unlist(Main),]

write.table(Normalized_TestData, file = "Normalized_TestData.txt", 
            sep = "\t", quote = FALSE)




####Training Samples for T helper Cells


n=ncol(TrainData_Th_Length_GC)
boxplot(as.matrix(log2(TrainData_Th_Length_GC[,3:n-2]+1)))

dataWithin <- withinLaneNormalization(x=as.matrix(TrainData_Th_Length_GC[,3:(n-2)]),
                                    y=TrainData_Th_Length_GC[,n], offset=FALSE,
                                    round=TRUE, which = "full")

boxplot(log2(dataWithin+1))

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

boxplot(log2(dataNorm+1))

Normalized_TrainData_Th <- data.frame(TrainData_Th_Length_GC[,1:2],dataNorm)




####Deleting Zeros after Normalization



ZeroGenes=data.frame()
Main=data.frame()

Total_0<-Normalized_TrainData_Th
Total_0[is.na(Total_0)] <- 0

for (i in 1:dim(Total_0)[1]) {
  
  if (sum(sum(Total_0[i,3:dim(Total_0)[2]]>=5))>3) {
    Main <-rbind(Main,i)
    print(i)
    
  } else {ZeroGenes <- rbind(ZeroGenes,i)
  print("Zero")}
}


Normalized_TrainData_Th <-Normalized_TrainData_Th[unlist(Main),]
write.table(Normalized_TrainData_Th, file = "Normalized_TrainData_Th.txt",
            sep = "\t", quote = FALSE)



####Testing Samples For T helper Cells



boxplot(as.matrix(log2(TestData_Th_Length_GC[,3:19]+1)))

dataWithin <- withinLaneNormalization(x=as.matrix(TestData_Th_Length_GC[,3:19]),
                                      y=TestData_Th_Length_GC[,21], offset=FALSE,
                                      round=TRUE, which = "full")

boxplot(log2(dataWithin+1))

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

boxplot(log2(dataNorm+1))

Normalized_TestData_Th <- data.frame(TestData_Th_Length_GC[,1:2],dataNorm)
Normalized_TestData_Th <-Normalized_TestData_Th[unlist(Main),]

write.table(Normalized_TestData_Th, file = "Normalized_TestData_Th.txt", 
            sep = "\t", quote = FALSE)
