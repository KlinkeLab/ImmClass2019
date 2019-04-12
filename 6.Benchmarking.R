#Benchmarking of gene lists against CIBERSORT gene lists

library(xlsx)
library(matlab)
library(MESS)
library(ggplot2)

# Importing Data ----------------------------------------------------------


Full_Test_Data <- read.table("Normalized_TestData.txt")
Full_Test_Data_Th <- read.table("Normalized_TestData_Th.txt")


# Assigning Gene Lists ----------------------------------------------------

#Assigning appropriate genes to each cell subset (Cibersort)
Full_Ciber <- read.xlsx('CiberList.xlsx', sheetIndex = 1, stringsAsFactors = FALSE)
Full_Ciber <- Full_Ciber[,-24]
B <- Full_Ciber[union(which(Full_Ciber$B.cells.naive == 1), which(Full_Ciber$B.cells.memory == 1)),1]
CD4 <- Full_Ciber[union(union(which(Full_Ciber$T.cells.CD4.naive == 1), which(Full_Ciber$T.cells.CD4.memory.resting == 1)), which(Full_Ciber$T.cells.CD4.memory.activated == 1)),1]
CD8 <- Full_Ciber[which(Full_Ciber$T.cells.CD8== 1),1]
DC <- Full_Ciber[union(which(Full_Ciber$Dendritic.cells.resting == 1), which(Full_Ciber$Dendritic.cells.activated == 1)),1]
M1 <- Full_Ciber[which(Full_Ciber$Macrophages.M1== 1),1]
M2 <- Full_Ciber[which(Full_Ciber$Macrophages.M2== 1),1]
M0 <- Full_Ciber[which(Full_Ciber$Macrophages.M0== 1),1]
Neu <- Full_Ciber[which(Full_Ciber$Neutrophils== 1),1]
NK <- Full_Ciber[union(which(Full_Ciber$NK.cells.resting == 1), which(Full_Ciber$NK.cells.activated == 1)),1]
Mono <- Full_Ciber[which(Full_Ciber$Monocytes == 1),1]

#creating vectors of equal size
MAX <- max(length(B), length(CD4), length(CD8), length(DC), length(M1), length(M2), length(M0), length(Neu), length(NK), length(Mono))
length(B) <- MAX
length(CD4) <- MAX
length(CD8) <- MAX
length(DC) <- MAX 
length(M1) <- MAX
length(M2) <- MAX
length(M0) <- MAX
length(Neu) <- MAX
length(NK) <- MAX 
length(Mono) <- MAX

# Th0 
# Th1
# Th2
# Th17
# Treg
# No comparable Cibersort categories for above

temp <- cbind(B,CD4,CD8,DC,M1,M2,M0,Neu,NK,Mono)
Ciber <- as.data.frame(temp)
rm(temp,B,CD4,CD8,DC,M1,M2,M0,Neu,NK,Mono)

Full_LogR <- read.xlsx('Separated_Genes_Modified.xlsx', sheetIndex = 1, stringsAsFactors = FALSE)
temp <- Full_LogR
colnames(temp) <- c(colnames(Ciber),'Th0','Th1','Th2','Th17','Treg')
LogR <- temp
rm(temp)

#Similar to above, repeated for single cell signatures
temp1<-cell.sig$B.cell
temp2<-cell.sig$T.CD4
temp3<-cell.sig$T.CD8
temp4<-cell.sig$Macrophage
temp5<-cell.sig$NK
MAX <- max(length(temp1), length(temp2),length(temp3),length(temp4),length(temp5))
length(temp1) <- MAX
length(temp2) <- MAX
length(temp3) <- MAX
length(temp4) <- MAX
length(temp5) <- MAX

SC <- t(rbind(temp1,temp2,temp3,temp4,temp5))
rm(temp1,temp2,temp3,temp4,temp5)
colnames(SC) <- colnames(Ciber)[c(1,2,3,7,9)]

# write.xlsx(as.table(SC), file = "C:/Users/Paraag/Desktop/Human Gene Selection/sccellsig.xlsx")

#reading in old heuristic signatures
Old <- read.xlsx('Oldlist.xlsx', sheetIndex = 1, stringsAsFactors = FALSE)


# Thresholding of Data (Convert to Binary Data) ---------------------------

# BEGINNING OF LOOP FOR 2D DENSITY PLOT
LogR_ROC_mat_all <- matrix(data = 0, nrow = 2, ncol = 4100)
rownames(LogR_ROC_mat_all) <- c('TPR','FPR')

Ciber_ROC_mat_all <- matrix(data = 0, nrow = 2, ncol = 4100)
rownames(Ciber_ROC_mat_all) <- c('TPR','FPR')

LogR2_ROC_mat_all <- matrix(data = 0, nrow = 2, ncol = 4100)
rownames(LogR2_ROC_mat_all) <- c('TPR','FPR')

SC_ROC_mat_all <- matrix(data = 0, nrow = 2, ncol = 4100)
rownames(SC_ROC_mat_all) <- c('TPR','FPR')

LogR3_ROC_mat_all <- matrix(data = 0, nrow = 2, ncol = 4100)
rownames(LogR3_ROC_mat_all) <- c('TPR','FPR')

Old_ROC_mat_all <- matrix(data = 0, nrow = 2, ncol = 4100)
rownames(Old_ROC_mat_all) <- c('TPR','FPR')

for (plot_counter in 0:99){

# Threshold <- 0.81
Threshold <- 2.48


Test_Data <- Full_Test_Data[,-1]
Test_Data <-log2(as.matrix(data.frame((Test_Data[,c(2:80)])))+1)
Test_Data[is.na(Test_Data)]<- -1
TN_Test_Data<-as.matrix(Test_Data[sample(nrow(Test_Data)),]); TN_Test_Data<-as.matrix(TN_Test_Data[,sample(ncol(Test_Data))])

Histogram_Data <- Test_Data

Test_Data <- Test_Data > Threshold
TN_Test_Data <- TN_Test_Data > Threshold

Test_Data <- as.data.frame(Test_Data)
Test_Data <- Test_Data[order(colnames(Test_Data))]
Test_Data <- cbind(Full_Test_Data[,2], Test_Data)
Test_Data <- Test_Data[-which(is.na(Test_Data[,1])),]

TN_Test_Data <- as.data.frame(TN_Test_Data)
TN_Test_Data <- TN_Test_Data[order(colnames(TN_Test_Data))]
TN_Test_Data <- cbind(Full_Test_Data[,2], TN_Test_Data)

# Fisher Exact Testing -----------------------------------------
# (Gene List# of Tables per Sample per benchmark)

#Logistic Regression Test Data
P_val_LogR_Test <- matrix(data = NA, nrow = (ncol(Test_Data) - 1), ncol = ncol(Ciber))
colnames(P_val_LogR_Test) <- colnames(Ciber)
rownames(P_val_LogR_Test) <- colnames(Test_Data)[2:ncol(Test_Data)]

for(j in 1:ncol(P_val_LogR_Test)){
  IN <- Test_Data[which(Test_Data[,1] %in% LogR[,j]),]
  OUT <- Test_Data[-which(Test_Data[,1] %in% LogR[,j]),]
  for(i in 1:nrow(P_val_LogR_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_LogR_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

#Logistic Regression Scrambled Data
P_val_LogR_TN_Test <- matrix(data = NA, nrow = (ncol(TN_Test_Data) - 1), ncol = ncol(Ciber))
colnames(P_val_LogR_TN_Test) <- colnames(Ciber)
rownames(P_val_LogR_TN_Test) <- colnames(TN_Test_Data)[2:ncol(TN_Test_Data)]

for(j in 1:ncol(P_val_LogR_TN_Test)){
  IN <- TN_Test_Data[which(TN_Test_Data[,1] %in% LogR[,j]),]
  OUT <- TN_Test_Data[-which(TN_Test_Data[,1] %in% LogR[,j]),]
  for(i in 1:nrow(P_val_LogR_TN_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_LogR_TN_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

#Cibersort Test Data
P_val_Ciber_Test <- matrix(data = NA, nrow = (ncol(Test_Data) - 1), ncol = ncol(Ciber))
colnames(P_val_Ciber_Test) <- colnames(Ciber)
rownames(P_val_Ciber_Test) <- colnames(Test_Data)[2:ncol(Test_Data)]

for(j in 1:ncol(P_val_Ciber_Test)){
  IN <- Test_Data[which(Test_Data[,1] %in% Ciber[,j]),]
  OUT <- Test_Data[-which(Test_Data[,1] %in% Ciber[,j]),]
  for(i in 1:nrow(P_val_Ciber_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_Ciber_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

#Cibersort Scrambled Data
P_val_Ciber_TN_Test <- matrix(data = NA, nrow = (ncol(TN_Test_Data) - 1), ncol = ncol(Ciber))
colnames(P_val_Ciber_TN_Test) <- colnames(Ciber)
rownames(P_val_Ciber_TN_Test) <- colnames(TN_Test_Data)[2:ncol(TN_Test_Data)]

for(j in 1:ncol(P_val_Ciber_TN_Test)){
  IN <- TN_Test_Data[which(TN_Test_Data[,1] %in% Ciber[,j]),]
  OUT <- TN_Test_Data[-which(TN_Test_Data[,1] %in% Ciber[,j]),]
  for(i in 1:nrow(P_val_Ciber_TN_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_Ciber_TN_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

#SC Test Data
P_val_SC_Test <- matrix(data = NA, nrow = (ncol(Test_Data) - 1), ncol = ncol(SC))
colnames(P_val_SC_Test) <- colnames(SC)
rownames(P_val_SC_Test) <- colnames(Test_Data)[2:ncol(Test_Data)]

for(j in 1:ncol(P_val_SC_Test)){
  IN <- Test_Data[which(Test_Data[,1] %in% SC[,j]),]
  OUT <- Test_Data[-which(Test_Data[,1] %in% SC[,j]),]
  for(i in 1:nrow(P_val_SC_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_SC_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

#SC Scrambled Data
P_val_SC_TN_Test <- matrix(data = NA, nrow = (ncol(TN_Test_Data) - 1), ncol = ncol(SC))
colnames(P_val_SC_TN_Test) <- colnames(SC)
rownames(P_val_SC_TN_Test) <- colnames(TN_Test_Data)[2:ncol(TN_Test_Data)]

for(j in 1:ncol(P_val_SC_TN_Test)){
  IN <- TN_Test_Data[which(TN_Test_Data[,1] %in% SC[,j]),]
  OUT <- TN_Test_Data[-which(TN_Test_Data[,1] %in% SC[,j]),]
  for(i in 1:nrow(P_val_SC_TN_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_SC_TN_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

P_val_SC_Test <- P_val_SC_Test[c(1:37,58:59,71:77),]
P_val_SC_TN_Test <- P_val_SC_TN_Test[c(1:37,58:59,71:77),]

P_val_LogR2_Test <- P_val_LogR_Test[c(1:37,58:59,71:77),]
P_val_LogR2_TN_Test <- P_val_LogR_TN_Test[c(1:37,58:59,71:77),]

#Old Test Data
P_val_Old_Test <- matrix(data = NA, nrow = (ncol(Test_Data) - 1), ncol = ncol(Old))
colnames(P_val_Old_Test) <- colnames(Old)
rownames(P_val_Old_Test) <- colnames(Test_Data)[2:ncol(Test_Data)]

for(j in 1:ncol(P_val_Old_Test)){
  IN <- Test_Data[which(Test_Data[,1] %in% Old[,j]),]
  OUT <- Test_Data[-which(Test_Data[,1] %in% Old[,j]),]
  for(i in 1:nrow(P_val_Old_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_Old_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

#Old Scrambled Data
P_val_Old_TN_Test <- matrix(data = NA, nrow = (ncol(TN_Test_Data) - 1), ncol = ncol(Old))
colnames(P_val_Old_TN_Test) <- colnames(Old)
rownames(P_val_Old_TN_Test) <- colnames(TN_Test_Data)[2:ncol(TN_Test_Data)]

for(j in 1:ncol(P_val_Old_TN_Test)){
  IN <- TN_Test_Data[which(TN_Test_Data[,1] %in% Old[,j]),]
  OUT <- TN_Test_Data[-which(TN_Test_Data[,1] %in% Old[,j]),]
  for(i in 1:nrow(P_val_Old_TN_Test)){
    mat <- matrix(data = NA, nrow = 2, ncol = 2)
    
    mat[1,1] <- sum(IN[,i+1] == TRUE)
    mat[1,2] <- sum(IN[,i+1] == FALSE)
    mat[2,1] <- sum(OUT[,i+1] == TRUE)
    mat[2,2] <- sum(OUT[,i+1] == FALSE)
    
    P_val_Old_TN_Test[i,j] <- fisher.test(x = mat, alternative = "greater",conf.int = FALSE)[[1]]
  }
} 

P_val_Old_Test <- P_val_Old_Test[c(41:44,58:59,71:77),]
P_val_Old_TN_Test <- P_val_Old_TN_Test[c(41:44,58:59,71:77),]

P_val_LogR3_Test <- P_val_LogR_Test[c(41:44,58:59,71:77),]
P_val_LogR3_TN_Test <- P_val_LogR_TN_Test[c(41:44,58:59,71:77),]


# Contingency Tables for ROC Curve ----------------------------------------
LogR_ROC_mat <- matrix(data = 0, nrow = 2, ncol = 41)
# colnames(LogR_ROC_mat) <- c('0','0.000001','0.00001','0.0001','0.001','0.01','0.1','0.2','0.3','0.4','0.5','1')
rownames(LogR_ROC_mat) <- c('TPR','FPR')
LogR_ROC_mat[1:2,1] <- c(0,0)
LogR_ROC_mat[1:2,41] <- c(1,1)

Ciber_ROC_mat <- matrix(data = 0, nrow = 2, ncol = 41)
# colnames(Ciber_ROC_mat) <- c('0','0.000001','0.00001','0.0001','0.001','0.01','0.1','0.2','0.3','0.4','0.5','1')
rownames(Ciber_ROC_mat) <- c('TPR','FPR')
Ciber_ROC_mat[1:2,1] <- c(0,0)
Ciber_ROC_mat[1:2,41] <- c(1,1)

LogR2_ROC_mat <- matrix(data = 0, nrow = 2, ncol = 41)
# colnames(LogR2_ROC_mat) <- c('0','0.000001','0.00001','0.0001','0.001','0.01','0.1','0.2','0.3','0.4','0.5','1')
rownames(LogR2_ROC_mat) <- c('TPR','FPR')
LogR2_ROC_mat[1:2,1] <- c(0,0)
LogR2_ROC_mat[1:2,41] <- c(1,1)

SC_ROC_mat <- matrix(data = 0, nrow = 2, ncol = 41)
# colnames(SC_ROC_mat) <- c('0','0.000001','0.00001','0.0001','0.001','0.01','0.1','0.2','0.3','0.4','0.5','1')
rownames(SC_ROC_mat) <- c('TPR','FPR')
SC_ROC_mat[1:2,1] <- c(0,0)
SC_ROC_mat[1:2,41] <- c(1,1)

Old_ROC_mat <- matrix(data = 0, nrow = 2, ncol = 41)
# colnames(Old_ROC_mat) <- c('0','0.000001','0.00001','0.0001','0.001','0.01','0.1','0.2','0.3','0.4','0.5','1')
rownames(Old_ROC_mat) <- c('TPR','FPR')
Old_ROC_mat[1:2,1] <- c(0,0)
Old_ROC_mat[1:2,41] <- c(1,1)

LogR3_ROC_mat <- matrix(data = 0, nrow = 2, ncol = 41)
# colnames(LogR3_ROC_mat) <- c('0','0.000001','0.00001','0.0001','0.001','0.01','0.1','0.2','0.3','0.4','0.5','1')
rownames(LogR3_ROC_mat) <- c('TPR','FPR')
LogR3_ROC_mat[1:2,1] <- c(0,0)
LogR3_ROC_mat[1:2,41] <- c(1,1)

counter = 1

for(ROC_Threshold in c(0.000001,0.00001,0.0001,
                       0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
                       0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
                       0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)) {
  counter <- counter + 1
  
  #Ciber Comparison LogR Table (1) --------------------------------------------------
  LogR_Tab <- matrix(data = 0,nrow = 2, ncol = 2)
  colnames(LogR_Tab) <- c('Condition A', 'Not A')
  rownames(LogR_Tab) <- c('Test says A', 'Test says Not A')
  
  for(i in 1: (nrow(P_val_LogR_Test)-2) ){
    if (i < 13){  #B
      if(P_val_LogR_Test[i,1] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,1] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 25){ #CD4
      if(P_val_LogR_Test[i,2] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,2] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 38){ #CD8
      if(P_val_LogR_Test[i,3] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,3] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 41){ #DC
      if(P_val_LogR_Test[i,4] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,4] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 43){ #M1
      if(P_val_LogR_Test[i,5] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,5] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 45){ #M2
      if(P_val_LogR_Test[i,6] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,6] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 58){ #Mono
      if(P_val_LogR_Test[i,10] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,10] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 60){ #M0
      if(P_val_LogR_Test[i,7] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,7] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 71){ #Neu
      if(P_val_LogR_Test[i,8] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,8] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    } else if (i < 78){ #NK
      if(P_val_LogR_Test[i,9] < ROC_Threshold){
        LogR_Tab[1,1] <- LogR_Tab[1,1] + 1
      } else {
        LogR_Tab[2,1] <- LogR_Tab[2,1] + 1
      }
      
      if(P_val_LogR_TN_Test[i,9] < ROC_Threshold){
        LogR_Tab[1,2] <- LogR_Tab[1,2] + 1
      } else {
        LogR_Tab[2,2] <- LogR_Tab[2,2] + 1
      }
    }
  }
  
  LogR_ROC_mat[1,counter] <- LogR_Tab[1,1]/(LogR_Tab[1,1] + LogR_Tab[2,1])
  LogR_ROC_mat[2,counter] <- LogR_Tab[1,2]/(LogR_Tab[1,2] + LogR_Tab[2,2])
  
  
  #Ciber Table ---------------------------
  Ciber_Tab <- matrix(data = 0,nrow = 2, ncol = 2)
  colnames(Ciber_Tab) <- c('Condition A', 'Not A')
  rownames(Ciber_Tab) <- c('Test says A', 'Test says Not A')
  
  for(i in 1: (nrow(P_val_Ciber_Test)-2) ){
    if (i < 13){  #B
      if(P_val_Ciber_Test[i,1] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,1] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 25){ #CD4
      if(P_val_Ciber_Test[i,2] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,2] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 38){ #CD8
      if(P_val_Ciber_Test[i,3] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,3] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 41){ #DC
      if(P_val_Ciber_Test[i,4] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,4] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 43){ #M1
      if(P_val_Ciber_Test[i,5] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,5] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 45){ #M2
      if(P_val_Ciber_Test[i,6] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,6] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 58){ #Mono
      if(P_val_Ciber_Test[i,10] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,10] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 60){ #M0
      if(P_val_Ciber_Test[i,7] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,7] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 71){ #Neu
      if(P_val_Ciber_Test[i,8] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,8] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    } else if (i < 78){ #NK
      if(P_val_Ciber_Test[i,9] < ROC_Threshold){
        Ciber_Tab[1,1] <- Ciber_Tab[1,1] + 1
      } else {
        Ciber_Tab[2,1] <- Ciber_Tab[2,1] + 1
      }
      
      if(P_val_Ciber_TN_Test[i,9] < ROC_Threshold){
        Ciber_Tab[1,2] <- Ciber_Tab[1,2] + 1
      } else {
        Ciber_Tab[2,2] <- Ciber_Tab[2,2] + 1
      }
    }
  }
  
  Ciber_ROC_mat[1,counter] <- Ciber_Tab[1,1]/(Ciber_Tab[1,1] + Ciber_Tab[2,1])
  Ciber_ROC_mat[2,counter] <- Ciber_Tab[1,2]/(Ciber_Tab[1,2] + Ciber_Tab[2,2])
  
  #Single Cell Comparison LogR Table (2) ---------------------------
  LogR2_Tab <- matrix(data = 0,nrow = 2, ncol = 2)
  colnames(LogR2_Tab) <- c('Condition A', 'Not A')
  rownames(LogR2_Tab) <- c('Test says A', 'Test says Not A')
  
  for(i in 1: (nrow(P_val_LogR2_Test))){
    if (i < 13){  #B
      if(P_val_LogR2_Test[i,1] < ROC_Threshold){
        LogR2_Tab[1,1] <- LogR2_Tab[1,1] + 1
      } else {
        LogR2_Tab[2,1] <- LogR2_Tab[2,1] + 1
      }
      
      if(P_val_LogR2_TN_Test[i,1] < ROC_Threshold){
        LogR2_Tab[1,2] <- LogR2_Tab[1,2] + 1
      } else {
        LogR2_Tab[2,2] <- LogR2_Tab[2,2] + 1
      }
    } else if (i < 25){ #CD4
      if(P_val_LogR2_Test[i,2] < ROC_Threshold){
        LogR2_Tab[1,1] <- LogR2_Tab[1,1] + 1
      } else {
        LogR2_Tab[2,1] <- LogR2_Tab[2,1] + 1
      }
      
      if(P_val_LogR2_TN_Test[i,2] < ROC_Threshold){
        LogR2_Tab[1,2] <- LogR2_Tab[1,2] + 1
      } else {
        LogR2_Tab[2,2] <- LogR2_Tab[2,2] + 1
      }
    } else if (i < 38){ #CD8
      if(P_val_LogR2_Test[i,3] < ROC_Threshold){
        LogR2_Tab[1,1] <- LogR2_Tab[1,1] + 1
      } else {
        LogR2_Tab[2,1] <- LogR2_Tab[2,1] + 1
      }
      
      if(P_val_LogR2_TN_Test[i,3] < ROC_Threshold){
        LogR2_Tab[1,2] <- LogR2_Tab[1,2] + 1
      } else {
        LogR2_Tab[2,2] <- LogR2_Tab[2,2] + 1
      }
    } else if (i < 40){ #M0
      if(P_val_LogR2_Test[i,4] < ROC_Threshold){
        LogR2_Tab[1,1] <- LogR2_Tab[1,1] + 1
      } else {
        LogR2_Tab[2,1] <- LogR2_Tab[2,1] + 1
      }
      
      if(P_val_LogR2_TN_Test[i,4] < ROC_Threshold){
        LogR2_Tab[1,2] <- LogR2_Tab[1,2] + 1
      } else {
        LogR2_Tab[2,2] <- LogR2_Tab[2,2] + 1
      }
    } else if (i < 47){ #NK
      if(P_val_LogR2_Test[i,5] < ROC_Threshold){
        LogR2_Tab[1,1] <- LogR2_Tab[1,1] + 1
      } else {
        LogR2_Tab[2,1] <- LogR2_Tab[2,1] + 1
      }
      
      if(P_val_LogR2_TN_Test[i,5] < ROC_Threshold){
        LogR2_Tab[1,2] <- LogR2_Tab[1,2] + 1
      } else {
        LogR2_Tab[2,2] <- LogR2_Tab[2,2] + 1
      }
    }
  }
  
  LogR2_ROC_mat[1,counter] <- LogR2_Tab[1,1]/(LogR2_Tab[1,1] + LogR2_Tab[2,1])
  LogR2_ROC_mat[2,counter] <- LogR2_Tab[1,2]/(LogR2_Tab[1,2] + LogR2_Tab[2,2])
  
  #Single Cell Table ------------------------------
  SC_Tab <- matrix(data = 0,nrow = 2, ncol = 2)
  colnames(SC_Tab) <- c('Condition A', 'Not A')
  rownames(SC_Tab) <- c('Test says A', 'Test says Not A')
 
  for(i in 1: (nrow(P_val_SC_Test))){
    if (i < 13){  #B
      if(P_val_SC_Test[i,1] < ROC_Threshold){
        SC_Tab[1,1] <- SC_Tab[1,1] + 1
      } else {
        SC_Tab[2,1] <- SC_Tab[2,1] + 1
      }
      
      if(P_val_SC_TN_Test[i,1] < ROC_Threshold){
        SC_Tab[1,2] <- SC_Tab[1,2] + 1
      } else {
        SC_Tab[2,2] <- SC_Tab[2,2] + 1
      }
    } else if (i < 25){ #CD4
      if(P_val_SC_Test[i,2] < ROC_Threshold){
        SC_Tab[1,1] <- SC_Tab[1,1] + 1
      } else {
        SC_Tab[2,1] <- SC_Tab[2,1] + 1
      }
      
      if(P_val_SC_TN_Test[i,2] < ROC_Threshold){
        SC_Tab[1,2] <- SC_Tab[1,2] + 1
      } else {
        SC_Tab[2,2] <- SC_Tab[2,2] + 1
      }
    } else if (i < 38){ #CD8
      if(P_val_SC_Test[i,3] < ROC_Threshold){
        SC_Tab[1,1] <- SC_Tab[1,1] + 1
      } else {
        SC_Tab[2,1] <- SC_Tab[2,1] + 1
      }
      
      if(P_val_SC_TN_Test[i,3] < ROC_Threshold){
        SC_Tab[1,2] <- SC_Tab[1,2] + 1
      } else {
        SC_Tab[2,2] <- SC_Tab[2,2] + 1
      }
    } else if (i < 40){ #M0
      if(P_val_SC_Test[i,4] < ROC_Threshold){
        SC_Tab[1,1] <- SC_Tab[1,1] + 1
      } else {
        SC_Tab[2,1] <- SC_Tab[2,1] + 1
      }
      
      if(P_val_SC_TN_Test[i,4] < ROC_Threshold){
        SC_Tab[1,2] <- SC_Tab[1,2] + 1
      } else {
        SC_Tab[2,2] <- SC_Tab[2,2] + 1
      }
    } else if (i < 47){ #NK
      if(P_val_SC_Test[i,5] < ROC_Threshold){
        SC_Tab[1,1] <- SC_Tab[1,1] + 1
      } else {
        SC_Tab[2,1] <- SC_Tab[2,1] + 1
      }
      
      if(P_val_SC_TN_Test[i,5] < ROC_Threshold){
        SC_Tab[1,2] <- SC_Tab[1,2] + 1
      } else {
        SC_Tab[2,2] <- SC_Tab[2,2] + 1
      }
    }
  }
  
  SC_ROC_mat[1,counter] <- SC_Tab[1,1]/(SC_Tab[1,1] + SC_Tab[2,1])
  SC_ROC_mat[2,counter] <- SC_Tab[1,2]/(SC_Tab[1,2] + SC_Tab[2,2])
  
  #Oldlist Comparison LogR Table (3) -------------------------
  LogR3_Tab <- matrix(data = 0,nrow = 2, ncol = 2)
  colnames(LogR3_Tab) <- c('Condition A', 'Not A')
  rownames(LogR3_Tab) <- c('Test says A', 'Test says Not A')
  
  for(i in 1: (nrow(P_val_LogR3_Test)) ){
    if (i < 3){  #M1
      if(P_val_LogR3_Test[i,5] < ROC_Threshold){
        LogR3_Tab[1,1] <- LogR3_Tab[1,1] + 1
      } else {
        LogR3_Tab[2,1] <- LogR3_Tab[2,1] + 1
      }
      
      if(P_val_LogR3_TN_Test[i,5] < ROC_Threshold){
        LogR3_Tab[1,2] <- LogR3_Tab[1,2] + 1
      } else {
        LogR3_Tab[2,2] <- LogR3_Tab[2,2] + 1
      }
    } else if (i < 5){ #M2
      if(P_val_LogR3_Test[i,6] < ROC_Threshold){
        LogR3_Tab[1,1] <- LogR3_Tab[1,1] + 1
      } else {
        LogR3_Tab[2,1] <- LogR3_Tab[2,1] + 1
      }
      
      if(P_val_LogR3_TN_Test[i,6] < ROC_Threshold){
        LogR3_Tab[1,2] <- LogR3_Tab[1,2] + 1
      } else {
        LogR3_Tab[2,2] <- LogR3_Tab[2,2] + 1
      }
    } else if (i < 7){ #M0
      if(P_val_LogR3_Test[i,7] < ROC_Threshold){
        LogR3_Tab[1,1] <- LogR3_Tab[1,1] + 1
      } else {
        LogR3_Tab[2,1] <- LogR3_Tab[2,1] + 1
      }
      
      if(P_val_LogR3_TN_Test[i,7] < ROC_Threshold){
        LogR3_Tab[1,2] <- LogR3_Tab[1,2] + 1
      } else {
        LogR3_Tab[2,2] <- LogR3_Tab[2,2] + 1
      }
    } else if (i < 14){ #NK
      if(P_val_LogR3_Test[i,9] < ROC_Threshold){
        LogR3_Tab[1,1] <- LogR3_Tab[1,1] + 1
      } else {
        LogR3_Tab[2,1] <- LogR3_Tab[2,1] + 1
      }
      
      if(P_val_LogR3_TN_Test[i,9] < ROC_Threshold){
        LogR3_Tab[1,2] <- LogR3_Tab[1,2] + 1
      } else {
        LogR3_Tab[2,2] <- LogR3_Tab[2,2] + 1
      }
    }
  }
  LogR3_ROC_mat[1,counter] <- LogR3_Tab[1,1]/(LogR3_Tab[1,1] + LogR3_Tab[2,1])
  LogR3_ROC_mat[2,counter] <- LogR3_Tab[1,2]/(LogR3_Tab[1,2] + LogR3_Tab[2,2])
  
  #Oldlist Table --------------------------
  Old_Tab <- matrix(data = 0,nrow = 2, ncol = 2)
  colnames(Old_Tab) <- c('Condition A', 'Not A')
  rownames(Old_Tab) <- c('Test says A', 'Test says Not A')
  
  for(i in 1: (nrow(P_val_Old_Test)) ){
    if (i < 3){  #M1
      if(P_val_Old_Test[i,1] < ROC_Threshold){
        Old_Tab[1,1] <- Old_Tab[1,1] + 1
      } else {
        Old_Tab[2,1] <- Old_Tab[2,1] + 1
      }
      
      if(P_val_Old_TN_Test[i,1] < ROC_Threshold){
        Old_Tab[1,2] <- Old_Tab[1,2] + 1
      } else {
        Old_Tab[2,2] <- Old_Tab[2,2] + 1
      }
    } else if (i < 5){ #M2
      if(P_val_Old_Test[i,2] < ROC_Threshold){
        Old_Tab[1,1] <- Old_Tab[1,1] + 1
      } else {
        Old_Tab[2,1] <- Old_Tab[2,1] + 1
      }
      
      if(P_val_Old_TN_Test[i,2] < ROC_Threshold){
        Old_Tab[1,2] <- Old_Tab[1,2] + 1
      } else {
        Old_Tab[2,2] <- Old_Tab[2,2] + 1
      }
    } else if (i < 7){ #M0
      if(P_val_Old_Test[i,3] < ROC_Threshold){
        Old_Tab[1,1] <- Old_Tab[1,1] + 1
      } else {
        Old_Tab[2,1] <- Old_Tab[2,1] + 1
      }
      
      if(P_val_Old_TN_Test[i,3] < ROC_Threshold){
        Old_Tab[1,2] <- Old_Tab[1,2] + 1
      } else {
        Old_Tab[2,2] <- Old_Tab[2,2] + 1
      }
    } else if (i < 14){ #NK
      if(P_val_Old_Test[i,4] < ROC_Threshold){
        Old_Tab[1,1] <- Old_Tab[1,1] + 1
      } else {
        Old_Tab[2,1] <- Old_Tab[2,1] + 1
      }
      
      if(P_val_Old_TN_Test[i,4] < ROC_Threshold){
        Old_Tab[1,2] <- Old_Tab[1,2] + 1
      } else {
        Old_Tab[2,2] <- Old_Tab[2,2] + 1
      }
    }
  }
  Old_ROC_mat[1,counter] <- Old_Tab[1,1]/(Old_Tab[1,1] + Old_Tab[2,1])
  Old_ROC_mat[2,counter] <- Old_Tab[1,2]/(Old_Tab[1,2] + Old_Tab[2,2])
  
}


LogR_ROC_mat_all[,c((plot_counter*41 + 1):(plot_counter*41 + 41))] <- LogR_ROC_mat
Ciber_ROC_mat_all[,c((plot_counter*41 + 1):(plot_counter*41 + 41))] <- Ciber_ROC_mat
LogR2_ROC_mat_all[,c((plot_counter*41 + 1):(plot_counter*41 + 41))] <- LogR2_ROC_mat
SC_ROC_mat_all[,c((plot_counter*41 + 1):(plot_counter*41 + 41))] <- SC_ROC_mat
LogR3_ROC_mat_all[,c((plot_counter*41 + 1):(plot_counter*41 + 41))] <- LogR3_ROC_mat
Old_ROC_mat_all[,c((plot_counter*41 + 1):(plot_counter*41 + 41))] <- Old_ROC_mat
}
#END OF LOOP FOR 2D DENSITY PLOT

X_Curve <- t(rbind(LogR_ROC_mat_all[2,],Ciber_ROC_mat_all[2,]))
Y_Curve <- t(rbind(LogR_ROC_mat_all[1,],Ciber_ROC_mat_all[1,]))
AUC_LogR <- auc(x=LogR_ROC_mat[2,], y=LogR_ROC_mat[1,])
AUC_Ciber <- auc(x=Ciber_ROC_mat[2,], y=Ciber_ROC_mat[1,])

Color=c("red","blue")
matplot(X_Curve,Y_Curve, type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")
legend("bottomright", legend = c("Logistic Regression","CIBERSORT"),col=Color, lty=c(1,2))

# write.xlsx(P_val_Ciber_Test, file = "C:/Users/Paraag/Desktop/Human Gene Selection/CiberTest.xlsx")
# write.xlsx(P_val_Ciber_TN_Test, file = "C:/Users/Paraag/Desktop/Human Gene Selection/CiberTNTest.xlsx")
# write.xlsx(P_val_LogR_Test, file = "C:/Users/Paraag/Desktop/Human Gene Selection/LogRTest.xlsx")
# write.xlsx(P_val_LogR_TN_Test, file = "C:/Users/Paraag/Desktop/Human Gene Selection/LogRTNTest.xlsx")

#CONTOUR PLOT
# LogR_Density_Data <- t(LogR_ROC_mat_all)
# LogR_Density_Data <- as.data.frame(LogR_Density_Data[,c(2,1)])
# colnames(LogR_Density_Data) <- c("x","y")
# 
# ggplot(LogR_Density_Data, aes(x=x, y=y) ) + geom_density_2d()
# 
# Ciber_Density_Data <- t(Ciber_ROC_mat_all)
# Ciber_Density_Data <- as.data.frame(Ciber_Density_Data[,c(2,1)])
# colnames(Ciber_Density_Data) <- c("x","y")
# 
# ggplot(Ciber_Density_Data, aes(x=x, y=y) ) + geom_density_2d()

LogR_CI <- matrix(data = 0, nrow = 3, ncol = 41)
Ciber_CI <- matrix(data = 0, nrow = 3, ncol = 41)
LogR2_CI <- matrix(data = 0, nrow = 3, ncol = 41)
SC_CI <- matrix(data = 0, nrow = 3, ncol = 41)
LogR3_CI <- matrix(data = 0, nrow = 3, ncol = 41)
Old_CI <- matrix(data = 0, nrow = 3, ncol = 41)

LogR_temp <- matrix(data = 0, nrow = 2, ncol = 100)
Ciber_temp <- matrix(data = 0, nrow = 2, ncol = 100)
LogR2_temp <- matrix(data = 0, nrow = 2, ncol = 100)
SC_temp <- matrix(data = 0, nrow = 2, ncol = 100)
LogR3_temp <- matrix(data = 0, nrow = 2, ncol = 100)
Old_temp <- matrix(data = 0, nrow = 2, ncol = 100)

for (offset in 0:40){
  counter <- 1
  for (i in seq(from = 1, to = 4100, by = 41)){
  pos <- offset + i
  
  LogR_temp[,counter] <- LogR_ROC_mat_all[,pos]
  Ciber_temp[,counter] <- Ciber_ROC_mat_all[,pos]
  LogR2_temp[,counter] <- LogR2_ROC_mat_all[,pos]
  SC_temp[,counter] <- SC_ROC_mat_all[,pos]
  LogR3_temp[,counter] <- LogR3_ROC_mat_all[,pos]
  Old_temp[,counter] <- Old_ROC_mat_all[,pos]
  
  counter <- counter + 1
  }
  LogR_CI[,(offset + 1)] <- t(c(t.test(LogR_temp[2,])$conf.int[1:2], LogR_temp[1,1]))
  Ciber_CI[,(offset + 1)] <- t(c(t.test(Ciber_temp[2,])$conf.int[1:2], Ciber_temp[1,1]))
  LogR2_CI[,(offset + 1)] <- t(c(t.test(LogR2_temp[2,])$conf.int[1:2], LogR2_temp[1,1]))
  SC_CI[,(offset + 1)] <- t(c(t.test(SC_temp[2,])$conf.int[1:2], SC_temp[1,1]))
  LogR3_CI[,(offset + 1)] <- t(c(t.test(LogR3_temp[2,])$conf.int[1:2], LogR3_temp[1,1]))
  Old_CI[,(offset + 1)] <- t(c(t.test(Old_temp[2,])$conf.int[1:2], Old_temp[1,1]))
}

matplot(c(LogR_CI[1,-41],LogR_CI[2,-41]),c(LogR_CI[3,-41],LogR_CI[3,-41]), type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")

matplot(c(Ciber_CI[1,-41],Ciber_CI[2,-41]),c(Ciber_CI[3,-41],Ciber_CI[3,-41]), type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")

matplot(c(LogR2_CI[1,-41],LogR2_CI[2,-41]),c(LogR2_CI[3,-41],LogR2_CI[3,-41]), type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")

matplot(c(SC_CI[1,-41],SC_CI[2,-41]),c(SC_CI[3,-41],SC_CI[3,-41]), type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")

matplot(c(LogR3_CI[1,-41],LogR3_CI[2,-41]),c(LogR3_CI[3,-41],LogR3_CI[3,-41]), type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")

matplot(c(Old_CI[1,-41],Old_CI[2,-41]),c(Old_CI[3,-41],Old_CI[3,-41]), type="b", col=Color, pch = 1, lwd = 1.5, xlab = "1-Specificity", ylab = "Sensitivity")

hist(Histogram_Data)

plot(LogR_CI[1,-41],LogR_CI[3,-41],type = 'n', ylim = c(0,1),xlim = c(0,1),
     ylab = 'Sensitivity', xlab = '1-Specificity', main = 'LogR vs. Cibersort')
lines(c(0,LogR_CI[1,-41]), c(0,LogR_CI[3,-41]), col = 'red')
lines(c(0,LogR_CI[2,-41]), c(0,LogR_CI[3,-41]), col = 'red')
polygon(c(LogR_CI[1,-41], rev(LogR_CI[2,-41])), c(LogR_CI[3,-41], rev(LogR_CI[3,-41])),
        col = "red", border = NA)

lines(c(0,Ciber_CI[1,-41]), c(0,Ciber_CI[3,-41]), col = 'blue')
lines(c(0.01,Ciber_CI[2,-41]), c(0.01,Ciber_CI[3,-41]), col = 'blue')
polygon(c(Ciber_CI[1,-41], rev(Ciber_CI[2,-41])), c(Ciber_CI[3,-41], rev(Ciber_CI[3,-41])),
        col = "blue", border = NA)

plot(LogR2_CI[1,-41],LogR2_CI[3,-41],type = 'n', ylim = c(0,1),xlim = c(0,1),
     ylab = 'Sensitivity', xlab = '1-Specificity', main = 'LogR vs Singe Cell')
lines(c(0,LogR2_CI[1,-41]), c(0,LogR2_CI[3,-41]), col = 'red')
lines(c(0,LogR2_CI[2,-41]), c(0,LogR2_CI[3,-41]), col = 'red')
polygon(c(LogR2_CI[1,-41], rev(LogR2_CI[2,-41])), c(LogR2_CI[3,-41], rev(LogR2_CI[3,-41])),
        col = "red", border = NA)

lines(c(0,SC_CI[1,-41]), c(0,SC_CI[3,-41]), col = 'green')
lines(c(0.01,SC_CI[2,-41]), c(0.01,SC_CI[3,-41]), col = 'green')
polygon(c(SC_CI[1,-41], rev(SC_CI[2,-41])), c(SC_CI[3,-41], rev(SC_CI[3,-41])),
        col = "green", border = NA)

plot(LogR3_CI[1,-41],LogR3_CI[3,-41],type = 'n', ylim = c(0,1),xlim = c(0,1),
     ylab = 'Sensitivity', xlab = '1-Specificity', main = 'LogR vs Heuristic')
lines(c(0,LogR3_CI[1,-41]), c(0,LogR3_CI[3,-41]), col = 'red')
lines(c(0,LogR3_CI[2,-41]), c(0,LogR3_CI[3,-41]), col = 'red')
polygon(c(LogR3_CI[1,-41], rev(LogR3_CI[2,-41])), c(LogR3_CI[3,-41], rev(LogR3_CI[3,-41])),
        col = "red", border = NA)

lines(c(0,Old_CI[1,-41]), c(0,Old_CI[3,-41]), col = 'yellow')
lines(c(0.01,Old_CI[2,-41]), c(0.01,Old_CI[3,-41]), col = 'yellow')
polygon(c(Old_CI[1,-41], rev(Old_CI[2,-41])), c(Old_CI[3,-41], rev(Old_CI[3,-41])),
        col = "yellow", border = NA)


