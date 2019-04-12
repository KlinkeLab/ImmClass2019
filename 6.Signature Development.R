##############################################Info#########################################
############# Title:          Gene Signatures Development
############# Description:    This Script develops Ten Gene Signatures for Immune Cells
#############                 and Five Gene Signatures for T helper Cells along with 
#############                 Related Heatmaps
############# Author:         Arezo Torang
############# Date:           8-Apr-2019
###########################################################################################-



rm(list = ls())

###Importing Data-----------------------------------------------------------




load("Gene Selection.RData")

Normalized_TrainData <- read.table("Normalized_TrainData.txt", sep = "\t",
                                   quote = "",stringsAsFactors=FALSE)

Normalized_TestData <- read.table("Normalized_TestData.txt", sep = "\t", 
                                  quote = "",stringsAsFactors=FALSE)

Normalized_TrainData_Th <- read.table("Normalized_TrainData_Th.txt", 
                                      sep = "\t", quote = "",
                                      stringsAsFactors=FALSE)

Normalized_TestData_Th <- read.table("Normalized_TestData_Th.txt", 
                                     sep = "\t", quote = "",
                                     stringsAsFactors=FALSE)






####Gene Signature Development for Immune Cells-----------


Main<-Normalized_TrainData[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Normalized_TrainData<-data.frame(Normalized_TrainData[,1:2],Main)

Main<-Normalized_TestData[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Normalized_TestData<-data.frame(Normalized_TestData[,1:2],Main)


## B Cell Gene Signature: 0:Other, 1:BCell------------


library(glmnet)
library(matlab)

Index <-coef(Lasso,s=0.0001)[["1"]]@i[-1]

x <-log2(as.matrix(data.frame(t(Normalized_TrainData[Index,3:110])))+1)
Test <-log2(as.matrix(data.frame(t(Normalized_TestData[Index,3:79])))+1)

y<-rbind(ones(14, ncol = 1),0*ones(94, ncol = 1))
y_Test<-rbind(ones(12, ncol = 1),0*ones(65, ncol = 1))

x[is.na(x)]<- -1
Test[is.na(Test)]<- -1



##Creating noisy data to make model neglect -1

X<-x
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y); print(i)
}

Lasso_BCell <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                         lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                         type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_BCell, newx=data.matrix(x),
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_BCell, newx=data.matrix(Test), 
                                      interval ="prediction",type="class", 
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy



##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0


for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_BCell, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_BCell, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_BCell, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


##ROC Curve for B Cell

library(DescTools)

pdf( "ROC_BCell.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", col=Color,
        lwd = 1.5, pch=15,xlab = "1-Specificity", 
        ylab = "Sensitivity",ylim = c(0,1))

Legend= c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
          "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
          "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
          "lambda=1e-4 (AUC=    )"))

legend("bottomright", legend = Legend, col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off() #End of Plot


library(pROC)
auc=c()

for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_BCell<-auc
auc_BCell


##Lambda=.01 has been selected for B Cell (28 Genes)

Index_BCell <-coef(Lasso_BCell,s=0.01)[["1"]]@i[-1]
Coef_BCell<-coef(Lasso_BCell,s=0.01)
Coef_BCell<-data.frame(Coef_BCell[["0"]]@x,Coef_BCell[["1"]]@x)
colnames(Coef_BCell)<-c("Other","B Cell")
rownames(Coef_BCell)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_BCell]
BCell_Genes<-as.character(Normalized_TrainData[Index,2])[Index_BCell]

##Heatmap

require(gplots)
require(RColorBrewer)

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_BCell.pdf", width = 5, height = 4 )
heatmap.2(as.matrix(Coef_BCell), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none")

dev.off()





## CD4 T Cell Gene Signature: 0:Other, 1:CD4-----------



y<-rbind(0*ones(14, ncol = 1),ones(15, ncol = 1),
         0*ones(79, ncol = 1))

y_Test<-rbind(0*ones(12, ncol = 1),ones(12, ncol = 1),
              0*ones(53, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_CD4 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                       lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_CD4, newx=data.matrix(x),
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)),arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_CD4, newx=data.matrix(Test),
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy



##Scrambling data to creat true negative samples and plot ROC


TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  prediction_Test <- data.frame(predict(Lasso_CD4, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_CD4, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_CD4, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_CD4.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", 
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", 
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_CD4<-auc
auc_CD4



##Lambda=.005 has been selected (30 Genes)


Index_CD4 <-coef(Lasso_CD4,s=0.005)[["1"]]@i[-1]
Coef_CD4<-coef(Lasso_CD4,s=0.005)
Coef_CD4<-data.frame(Coef_CD4[["0"]]@x,Coef_CD4[["1"]]@x)
colnames(Coef_CD4)<-c("Other","CD4")
rownames(Coef_CD4)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_CD4]
CD4_Genes<-as.character(Normalized_TrainData[Index,2])[Index_CD4]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_CD4.pdf", width = 5, height = 4 )
heatmap.2(as.matrix(Coef_CD4),breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55,
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()




## CD8 T Cell Gene Signature: 0:Other, 1:CD8-----------------



y<-rbind(0*ones(29, ncol = 1),ones(16, ncol = 1),
         0*ones(63, ncol = 1))

y_Test<-rbind(0*ones(24, ncol = 1),ones(13, ncol = 1),
              0*ones(40, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_CD8 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                       lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_CD8, newx=data.matrix(x),
                                 interval ="prediction",
                                 type="response", s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_CD8, newx=data.matrix(Test),
                                      interval ="prediction",
                                      type="class", s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy



##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  prediction_Test <- data.frame(predict(Lasso_CD8, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_CD8, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_CD8, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}

pdf( "ROC_CD8.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", col=Color,
        lwd = 1.5, pch=15,xlab = "1-Specificity",
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend ,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_CD8<-auc
auc_CD8


##Lambda=.01 has been selected (24 Genes)

Index_CD8 <-coef(Lasso_CD8,s=0.01)[["1"]]@i[-1]
Coef_CD8<-coef(Lasso_CD8,s=0.01)
Coef_CD8<-data.frame(Coef_CD8[["0"]]@x,Coef_CD8[["1"]]@x)
colnames(Coef_CD8)<-c("Other","CD8")
rownames(Coef_CD8)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_CD8]
CD8_Genes<-as.character(Normalized_TrainData[Index,2])[Index_CD8]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_CD8.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_CD8), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none")

dev.off()




## Monocyte Gene Signature: 0:Other, 1:Mono----------------



y<-rbind(0*ones(61, ncol = 1),ones(16, ncol = 1),
         0*ones(31, ncol = 1))

y_Test<-rbind(0*ones(44, ncol = 1),ones(13, ncol = 1),
              0*ones(20, ncol = 1))

##Creating noisy data to make model neglect -1

X<-x 
Y<-y
for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Mono <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                        lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Mono, newx=data.matrix(x),
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Mono, newx=data.matrix(Test),
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy



##Scrambling data to creat true negative samples and plot ROC


TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_Mono, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Mono, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Mono, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Mono.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_Mono<-auc
auc_Mono


##Lambda=.01 has been selected (25 Genes)

Index_Mono <-coef(Lasso_Mono,s=0.01)[["1"]]@i[-1]
Coef_Mono<-coef(Lasso_Mono,s=0.01)
Coef_Mono<-data.frame(Coef_Mono[["0"]]@x,Coef_Mono[["1"]]@x)
colnames(Coef_Mono)<-c("Other","Mono")
rownames(Coef_Mono)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_Mono]
Mono_Genes<-as.character(Normalized_TrainData[Index,2])[Index_Mono]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Mono.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Mono), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55,
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()



## Neutrophil Gene Signature: 0:Other, 1:Neu------------



y<-rbind(0*ones(82, ncol = 1),ones(11, ncol = 1),
         0*ones(15, ncol = 1))

y_Test<-rbind(0*ones(59, ncol = 1),ones(11, ncol = 1),
              0*ones(7, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Neu <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                       lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Neu, newx=data.matrix(x),
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Neu, newx=data.matrix(Test), 
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy



##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_Neu, newx=data.matrix(Test), 
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Neu, newx=data.matrix(Test), 
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Neu, newx=data.matrix(TN_Test), 
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0; FN=0; TN=0; FP=0; R=R+1
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Neu.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_Neu<-auc
auc_Neu


##Lambda=.01 has been selected (38 Genes)

Index_Neu <-coef(Lasso_Neu,s=0.01)[["1"]]@i[-1]
Coef_Neu<-coef(Lasso_Neu,s=0.01)
Coef_Neu<-data.frame(Coef_Neu[["0"]]@x,Coef_Neu[["1"]]@x)
colnames(Coef_Neu)<-c("Other","Neu")
rownames(Coef_Neu)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_Neu]
Neu_Genes<-as.character(Normalized_TrainData[Index,2])[Index_Neu]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Neu.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Neu), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.45, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()



## NK Cell Gene Signature: 0:Other, 1:Nk--------------



y<-rbind(0*ones(93, ncol = 1),ones(15, ncol = 1))
y_Test<-rbind(0*ones(70, ncol = 1),ones(7, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Nk <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Nk, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Nk, newx=data.matrix(Test), 
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


###Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_Nk, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Nk, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Nk, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Nk.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_NK<-auc
auc_NK


##Lambda=.01 has been selected (26 Genes)

Index_Nk <-coef(Lasso_Nk,s=0.01)[["1"]]@i[-1]
Coef_Nk<-coef(Lasso_Nk,s=0.01)
Coef_Nk<-data.frame(Coef_Nk[["0"]]@x,Coef_Nk[["1"]]@x)
colnames(Coef_Nk)<-c("Other","NK")
rownames(Coef_Nk)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_Nk]
Nk_Genes<-as.character(Normalized_TrainData[Index,2])[Index_Nk]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Nk.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Nk), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none")

dev.off()



## Dendritic Cell Gene Signature: 0:Other, 1:DC--------------



y<-rbind(0*ones(45, ncol = 1),ones(8, ncol = 1),
         0*ones(55, ncol = 1))

y_Test<-rbind(0*ones(37, ncol = 1),ones(3, ncol = 1),
              0*ones(37, ncol = 1))

##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_DC <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_DC, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_DC, newx=data.matrix(Test), 
                                      interval ="prediction",type="class", 
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1; R=0
  prediction_Test <- data.frame(predict(Lasso_DC, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_DC, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_DC, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}

pdf( "ROC_DC.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_DC<-auc
auc_DC



##Lambda=.05 has been selected (21 Genes)

Index_DC <-coef(Lasso_DC,s=0.01)[["1"]]@i[-1]
Coef_DC<-coef(Lasso_DC,s=0.01)
Coef_DC<-data.frame(Coef_DC[["0"]]@x,Coef_DC[["1"]]@x)
colnames(Coef_DC)<-c("Other","DC")
rownames(Coef_DC)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_DC]
DC_Genes<-as.character(Normalized_TrainData[Index,2])[Index_DC]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_DC.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_DC[-1,]), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55,
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none",dendrogram = "row") 

dev.off()




## M1 Gene Signature: 0:Other, 1:M1--------------



y<-rbind(0*ones(53, ncol = 1),ones(4, ncol = 1),
         0*ones(51, ncol = 1))

y_Test<-rbind(0*ones(40, ncol = 1),ones(2, ncol = 1),
              0*ones(35, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_M1 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_M1, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.05))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_M1, newx=data.matrix(Test),
                                      interval ="prediction",type="class",
                                      s=0.05))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_M1, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_M1, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_M1, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0; FN=0; TN=0; FP=0; R=R+1
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_M1.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity, Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", 
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_M1<-auc
auc_M1

##Lambda=.1 has been selected (5 Genes)

Index_M1 <-coef(Lasso_M1,s=0.1)[["1"]]@i[-1]
Coef_M1<-coef(Lasso_M1,s=0.1)
Coef_M1<-data.frame(Coef_M1[["0"]]@x,Coef_M1[["1"]]@x)
colnames(Coef_M1)<-c("Other","M1")
rownames(Coef_M1)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_M1]
M1_Genes<-as.character(Normalized_TrainData[Index,2])[Index_M1]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_M1.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_M1), breaks= seq(-.1, .1, length.out = 21), 
          scale="none",col = Color, margins=c(5,5), cexRow=0.55,
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()



## M2 Gene Signature: 0:Other, 1:M2--------------



y<-rbind(0*ones(57, ncol = 1),ones(4, ncol = 1),
         0*ones(47, ncol = 1))

y_Test<-rbind(0*ones(42, ncol = 1),ones(2, ncol = 1),
              0*ones(33, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_M2 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_M2, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.05))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_M2, newx=data.matrix(Test), 
                                      interval ="prediction",type="class",
                                      s=0.05))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_M2, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_M2, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_M2, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0; FN=0; TN=0; FP=0; R=R+1
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_M2.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", 
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", 
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_M2<-auc
auc_M2


##Lambda=.05 has been selected (7 Genes)

Index_M2 <-coef(Lasso_M2,s=0.05)[["1"]]@i[-1]
Coef_M2<-coef(Lasso_M2,s=0.05)
Coef_M2<-data.frame(Coef_M2[["0"]]@x,Coef_M2[["1"]]@x)
colnames(Coef_M2)<-c("Other","M2")
rownames(Coef_M2)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_M2]
M2_Genes<-as.character(Normalized_TrainData[Index,2])[Index_M2]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_M2.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_M2), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()




## Macrophage Gene Signature: 0:Other, 1:MΦ------------



y<-rbind(0*ones(77, ncol = 1),ones(5, ncol = 1),
         0*ones(26, ncol = 1))

y_Test<-rbind(0*ones(57, ncol = 1),ones(2, ncol = 1),
              0*ones(18, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_MQ <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_MQ, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_MQ, newx=data.matrix(Test), 
                                      interval ="prediction",type="class", 
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1; R=0
  prediction_Test <- data.frame(predict(Lasso_MQ, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_MQ, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_MQ, newx=data.matrix(TN_Test), 
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}

pdf( "ROC_MQ.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", 
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", 
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_MQ<-auc
auc_MQ


##Lambda=.05 has been selected (15 Genes)

Index_MQ <-coef(Lasso_MQ,s=0.05)[["1"]]@i[-1]
Coef_MQ<-coef(Lasso_MQ,s=0.05)
Coef_MQ<-data.frame(Coef_MQ[["0"]]@x,Coef_MQ[["1"]]@x)
colnames(Coef_MQ)<-c("Other","MQ")
rownames(Coef_MQ)[-1]<-as.character(Normalized_TrainData[Index,2])[Index_MQ]
MQ_Genes<-as.character(Normalized_TrainData[Index,2])[Index_MQ]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_MQ.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_MQ), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()




####Gene Signature Development for T Helper Cells-----------



Main<-Normalized_TrainData_Th[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Normalized_TrainData_Th<-data.frame(Normalized_TrainData_Th[,1:2],Main)
Main<-Normalized_TestData_Th[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Normalized_TestData_Th<-data.frame(Normalized_TestData_Th[,1:2],Main)



## Th0 Gene Signature: 0:Other, 1:Th0-----------



Index_Th <-coef(Lasso_Th,s=.05)[["1"]]@i[-1]

x <-log2(as.matrix(data.frame(t(Normalized_TrainData_Th[Index_Th,-c(1,2)])))+1)
Test <-log2(as.matrix(data.frame(t(Normalized_TestData_Th[Index_Th,-c(1,2)])))+1) 

y<-rbind(0*ones(9, ncol = 1),ones(9, ncol = 1),0*ones(13, ncol = 1))
y_Test<-rbind(0*ones(6, ncol = 1),ones(6, ncol = 1),0*ones(5, ncol = 1))

x[is.na(x)]<- -1
Test[is.na(Test)]<- -1


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Th0 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                       lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                   type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Th0, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Th0, newx=data.matrix(Test), 
                                      interval ="prediction",type="class",
                                      s=0.05))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy



##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(17,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1; R=0
  prediction_Test <- data.frame(predict(Lasso_Th0, newx=data.matrix(Test),
                                        interval ="prediction",type="response", 
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Th0, newx=data.matrix(Test),
                                              interval ="prediction",type="class", 
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Th0, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Th0.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity",
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_Th0<-auc 
auc_Th0


##Lambda=.001 has been selected (39 Genes)

Index_Th0 <-coef(Lasso_Th0,s=0.001)[["1"]]@i[-1]
Coef_Th0<-coef(Lasso_Th0,s=0.001)
Coef_Th0<-data.frame(Coef_Th0[["0"]]@x,Coef_Th0[["1"]]@x)
colnames(Coef_Th0)<-c("Other","Th0")
rownames(Coef_Th0)[-1]<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th0]
Th0_Genes<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th0]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Th0.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Th0), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.45, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()





## Th1 Gene Signature: 0:Other, 1:Th1-------



y<-rbind(0*ones(18, ncol = 1),ones(5, ncol = 1),
         0*ones(8, ncol = 1))

y_Test<-rbind(0*ones(12, ncol = 1),ones(2, ncol = 1),
              0*ones(3, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Th1 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                       lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Th1, newx=data.matrix(x),
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Th1, newx=data.matrix(Test),
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(17,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_Th1, newx=data.matrix(Test), 
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Th1, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Th1, newx=data.matrix(TN_Test), 
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0; FN=0; TN=0; FP=0; R=R+1
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Th1.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", 
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity",
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")
legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_Th1<-auc 
auc_Th1


##Lambda=.1 has been selected (12 Genes)

Index_Th1 <-coef(Lasso_Th1,s=0.1)[["1"]]@i[-1]
Coef_Th1<-coef(Lasso_Th1,s=0.1)
Coef_Th1<-data.frame(Coef_Th1[["0"]]@x,Coef_Th1[["1"]]@x)
colnames(Coef_Th1)<-c("Other","Th1")
rownames(Coef_Th1)[-1]<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th1]
Th1_Genes<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th1]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Th1.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Th1[-1,]), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, cexCol=1.0,
          key=TRUE, keysize=2,key.title = "Color Key",trace="none", 
          dendrogram = "row") 

dev.off()




## Th2 Gene Signature: 0:Other, 1:Th2-----------



y<-rbind(0*ones(26, ncol = 1),ones(5, ncol = 1))
y_Test<-rbind(0*ones(15, ncol = 1),ones(2, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Th2 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                       lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Th2, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Th2, newx=data.matrix(Test),
                                      interval ="prediction",type="class", 
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(17,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_Th2, newx=data.matrix(Test), 
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Th2, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Th2, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Th2.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity",
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_Th2<-auc 
auc_Th2


##Lambda=.01 has been selected (16 Genes)

Index_Th2 <-coef(Lasso_Th2,s=0.01)[["1"]]@i[-1]
Coef_Th2<-coef(Lasso_Th2,s=0.01)
Coef_Th2<-data.frame(Coef_Th2[["0"]]@x,Coef_Th2[["1"]]@x)
colnames(Coef_Th2)<-c("Other","Th2")
rownames(Coef_Th2)[-1]<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th2]
Th2_Genes<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th2]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Th2.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Th2[-1,]), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, cexCol=1.0, 
          key=TRUE, keysize=2,key.title = "Color Key",
          trace="none",dendrogram = "row")

dev.off()




## Th17 Gene Signature: 0:Other, 1:Th17--------------



y<-rbind(0*ones(23, ncol = 1),ones(3, ncol = 1),
         0*ones(5, ncol = 1))

y_Test<-rbind(0*ones(14, ncol = 1),ones(1, ncol = 1),
              0*ones(2, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_Th17 <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                        lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10) 

prediction <- data.frame(predict(Lasso_Th17, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_Th17, newx=data.matrix(Test),
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(17,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_Th17, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Th17, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Th17, newx=data.matrix(TN_Test),
                                           interval ="prediction",type="response", 
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0
    FN=0
    TN=0
    FP=0
    R=R+1
    
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_Th17.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l", col=Color,
        lwd = 1.5, pch=15, xlab = "1-Specificity",
        ylab = "Sensitivity",ylim = c(0,1))

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_Th17<-auc 
auc_Th17


##Lambda=.01 has been selected (14 Genes)

Index_Th17 <-coef(Lasso_Th17,s=0.01)[["1"]]@i[-1]
Coef_Th17<-coef(Lasso_Th17,s=0.01)
Coef_Th17<-data.frame(Coef_Th17[["0"]]@x,Coef_Th17[["1"]]@x)
colnames(Coef_Th17)<-c("Other","Th17")
rownames(Coef_Th17)[-1]<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th17]
Th17_Genes<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_Th17]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_Th17.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_Th17), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55,
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none")

dev.off()



## iTreg Gene Signature: 0:Other, 1:iTreg--------------



y<-rbind(ones(9, ncol = 1),0*ones(22, ncol = 1))
y_Test<-rbind(ones(6, ncol = 1),0*ones(11, ncol = 1))


##Creating noisy data to make model neglect -1

X<-x 
Y<-y

for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
  print(i)
}

Lasso_iTreg <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                         lambda = c(.1,.05,.01,.005,.001,.0005,.0001),
                       type.multinomial="grouped", nfolds=10)

prediction <- data.frame(predict(Lasso_iTreg, newx=data.matrix(x), 
                                 interval ="prediction",type="response",
                                 s=0.01))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso_iTreg, newx=data.matrix(Test), 
                                      interval ="prediction",type="class",
                                      s=0.01))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy


##Scrambling data to creat true negative samples and plot ROC

TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-2*ones(17,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0

for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  
  prediction_Test <- data.frame(predict(Lasso_iTreg, newx=data.matrix(Test),
                                        interval ="prediction",type="response",
                                        s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_iTreg, newx=data.matrix(Test),
                                              interval ="prediction",type="class",
                                              s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_iTreg, newx=data.matrix(TN_Test), 
                                           interval ="prediction",type="response",
                                           s=i))
  
  for (t in seq(0,1,length.out =100)) {
    TP=0; FN=0; TN=0; FP=0; R=R+1
    for (j in 1:length(TN_y)) {
      if (rowMaxs(as.matrix(prediction_Test[j,]))>t) {
        if (y_Test[j]==as.matrix(prediction_Test_Class[j,1])) {
          TP=TP+1
        } else {FN=FN+1}
      } else {FN=FN+1}
      if (rowMaxs(as.matrix(prediction_TN_Test[j,]))<t) {
        TN=TN+1
      } else {FP=FP+1}
    }
    
    Sensitivity[R,C]=TP/(TP+FN)
    Specificity[R,C]=TN/(TN+FP)
  }
}


pdf( "ROC_iTreg.pdf", width = 5, height =5 )

Color=c("green","purple","red","blue",
        "orange","black","yellow")

matplot(1-Specificity,Sensitivity, type="l",
        col=Color,lwd = 1.5, pch=15,
        xlab = "1-Specificity", 
        ylab = "Sensitivity")

Legend=c("lambda=1e-1 (AUC=    )","lambda=5e-2 (AUC=    )",
         "lambda=1e-2 (AUC=    )","lambda=5e-3 (AUC=    )",
         "lambda=1e-3 (AUC=    )","lambda=5e-4 (AUC=    )",
         "lambda=1e-4 (AUC=    )")

legend("bottomright", legend = Legend,col=Color,
       lty=1:7, lwd=.8, cex=0.6)

dev.off()


auc=c()
for (i in 1:7) {
  auc<-c(auc,AUC(1-Specificity[,i],Sensitivity[,i]))
}

auc_iTreg<-auc 
auc_iTreg


##Lambda=.005 has been selected (23 Genes)

Index_iTreg <-coef(Lasso_iTreg,s=0.005)[["1"]]@i[-1]
Coef_iTreg<-coef(Lasso_iTreg,s=0.005)
Coef_iTreg<-data.frame(Coef_iTreg[["0"]]@x,Coef_iTreg[["1"]]@x)
colnames(Coef_iTreg)<-c("Other","iTreg")
rownames(Coef_iTreg)[-1]<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_iTreg]
iTreg_Genes<-as.character(Normalized_TrainData_Th[Index_Th,2])[Index_iTreg]


##Heatmap

Color <- colorRampPalette(c("blue","white", "red"))(n = 20)

pdf( "heatmap_iTreg.pdf", width = 5, height = 4 )

heatmap.2(as.matrix(Coef_iTreg), breaks= seq(-.1, .1, length.out = 21),
          scale="none",col = Color,margins=c(5,5), cexRow=0.55, 
          cexCol=1.0, key=TRUE, keysize=2,key.title = "Color Key",
          trace="none") 

dev.off()




####Table of Gene Signatures-------



n <- max(length(BCell_Genes),length(CD4_Genes),length(CD8_Genes),
         length(DC_Genes),length(M1_Genes),length(M2_Genes),
         length(MQ_Genes),length(Neu_Genes),length(Nk_Genes),
         length(Mono_Genes),length(Th0_Genes),length(Th1_Genes),
         length(Th2_Genes),length(Th17_Genes),length(iTreg_Genes))

length(BCell_Genes)<-n
length(CD4_Genes)<-n
length(CD8_Genes)<-n
length(DC_Genes)<-n
length(M1_Genes)<-n
length(M2_Genes)<-n
length(MQ_Genes)<-n
length(Neu_Genes)<-n
length(Nk_Genes)<-n
length(Mono_Genes)<-n

length(Th0_Genes)<-n
length(Th1_Genes)<-n
length(Th2_Genes)<-n
length(Th17_Genes)<-n
length(iTreg_Genes)<-n


Separated_Genes<-data.frame(BCell_Genes,CD4_Genes,CD8_Genes,
                            DC_Genes,M1_Genes,M2_Genes,MQ_Genes,
                            Neu_Genes,Nk_Genes,Mono_Genes,
                            Th0_Genes,Th1_Genes,Th2_Genes,
                            Th17_Genes,iTreg_Genes)

write.xlsx(Separated_Genes,"Separated_Genes.xlsx")

a<-as.vector(as.matrix(Separated_Genes[,1:10]))
a<-data.frame(a[!is.na(a)])
a<-unique(a)
colnames(a)<-"ID"

Main_Selected_Genes<-merge(a,Selected_Genes)
rownames(Main_Selected_Genes)<-Main_Selected_Genes[,1]
Main<-Main_Selected_Genes[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Main_Selected_Genes<-data.frame(Main_Selected_Genes[,1:2],Main)
Main_Selected_Genes<-Main_Selected_Genes[,c(2,1,3:110)]

write.table(Main_Selected_Genes,"Main_Selected_Genes.txt")


a<-as.vector(as.matrix(Separated_Genes[,11:15]))
a<-data.frame(a[!is.na(a)])
a<-unique(a)
colnames(a)<-"ID"

Main_Selected_Genes_Th<-merge(a,Selected_Genes_Th)
rownames(Main_Selected_Genes_Th)<-Main_Selected_Genes_Th[,1]
Main<-Main_Selected_Genes_Th[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Main_Selected_Genes_Th<-data.frame(Main_Selected_Genes_Th[,1:2],Main)
Main_Selected_Genes_Th<-Main_Selected_Genes_Th[,c(2,1,3:33)]

write.table(Main_Selected_Genes_Th,"Main_Selected_Genes_Th.txt")



####Heatmap over Seprated Genes Expression-----------


## Expression for Immune Cells


a<-Main_Selected_Genes
a<-log2(a[,-c(1,2)]+1)
a[is.na(a)]<--1

Rowv <- rowMeans(as.matrix(a), na.rm = TRUE)
distr <- dist(as.matrix(a))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
ddr <- reorder(ddr, Rowv)
rowInd <- order.dendrogram(ddr)

plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")


Colv <- colMeans(as.matrix(a), na.rm = TRUE)
distc <- dist(as.matrix(t(a)))
hcc <- hclust(distc)
ddc <- as.dendrogram(hcc)
ddc <- reorder(ddc, Colv)
colInd <- order.dendrogram(ddc)

plot(ddc, axes = FALSE, xaxs = "i",leaflab = "none")



## Expression for T Helper Cells


b<-Main_Selected_Genes_Th
b<-log2(b[,-c(1,2)]+1)
b[is.na(b)]<--1

Rowv_Th <- rowMeans(as.matrix(b), na.rm = TRUE)
distr_Th <- dist(as.matrix(b))
hcr_Th <- hclust(distr_Th)
ddr_Th <- as.dendrogram(hcr_Th)
ddr_Th <- reorder(ddr_Th, Rowv_Th)
rowInd_Th <- order.dendrogram(ddr_Th)

plot(ddr_Th, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")


Colv_Th <- colMeans(as.matrix(b), na.rm = TRUE)
distc_Th <- dist(as.matrix(t(b)))
hcc_Th <- hclust(distc_Th)
ddc_Th <- as.dendrogram(hcc_Th)
ddc_Th <- reorder(ddc_Th, Colv_Th)
colInd_Th <- order.dendrogram(ddc_Th)

plot(ddc_Th, axes = FALSE, xaxs = "i",leaflab = "none")



##Heatmaps

pdf("Heatmap_Main_Genes.pdf", width = 17, height = 11)

nf <- layout(matrix(c(1,4,7,10,2,5,8,11,3,6,9,12),
                    3,4,byrow = TRUE), widths = c(2,7,2,7),
             heights = c(.5,.3,10), respect = TRUE)

##Immune Cells
par(mar = c(.5, 4, 0, 1))

image(as.matrix(seq(-.01, .01, length.out = 12)), axes = FALSE,
      col=c("lightblue",colorRampPalette(c("white","red"))(n = 11)))

mtext(c("NA","0","6+"), side = 3, las = 1,
      line = 0.5, cex = 1,
      at = c(-.07, .09, .98))

par(mar = c(.5, 0.5, 0, 1))
Color=c("deeppink","yellow","green","orange",
        "deepskyblue","Gold","gray","violet",
        "cyan","lightpink1")

image(as.matrix(x=1:10),col=Color, axes = FALSE)

x=c("DC","M1","BCell","Neu","NK","M2",
    expression(paste("M", phi)),"CD4",
    "CD8","Mono")

text(seq(0,1,1/(length(x)-1)), .1, x,cex = .8)

par(mar = c(.5, .5, 0, 1))
plot(ddr, horiz = TRUE, axes = FALSE,
     yaxs = "i", leaflab = "none")

par(mar = c(.5, .5, .5, 3.5))
plot(ddc, axes = FALSE, xaxs = "i",
     leaflab = "none")

par(mar = c(.5, 0.5, 0, 3.5))
table(substr(colnames(a),1,2))
x=c(3*ones(14,ncol=1),8*ones(15,ncol=1),
    9*ones(16,ncol=1),1*ones(8,ncol=1),
    2*ones(4,ncol=1),6*ones(4,ncol=1),
    10*ones(16,ncol=1),7*ones(5,ncol=1),
    4*ones(11,ncol=1),5*ones(15,ncol=1))

image(as.matrix(x[colInd]),col=Color, 
      axes = FALSE)

par(mar = c(.5, .5, 0, 3.5))
image(as.matrix(t(a[rowInd,colInd])),breaks= c(seq(-7, 6, length.out = 21),100),
      col=colorRampPalette(c("blue","white", "red"))(n = 21), axes = FALSE)

ylab <- rownames(a)[rowInd]
xlab <- colnames(a)[colInd]

mtext("Genes", side = 4, adj = 1,
      line = 1.5, cex = 2, at = .5)
# mtext(xlab, side = 3, las = 2, line = 0.5, cex = .4, at = seq(0, 1, 1/(length(xlab)-1)))



##T Helper Cells


par(mar = c(.5, 4, 0, 1))
image(as.matrix(seq(-.01, .01, length.out = 12)), axes = FALSE,
      col=c("lightblue",colorRampPalette(c("white","red"))(n = 11)))

mtext(c("NA","0","6+"), side = 3, 
      las = 1, line = 0.5, cex = 1,
      at = c(-.07, .09, .98))

par(mar = c(.5, 0.5, 0, 1))
Color=c("deeppink","yellow",
        "green","violet","cyan")

image(as.matrix(x=1:5),col=Color, axes = FALSE)
x=c("Th1","Th0","Th2","iTreg","Th17")
text(seq(0,1,1/(length(x)-1)), .1, x,cex = 1)

par(mar = c(.5, .5, 0, 1))
plot(ddr_Th, horiz = TRUE, axes = FALSE,
     yaxs = "i", leaflab = "none")

par(mar = c(.5, .5, .5, 3.5))
plot(ddc_Th, axes = FALSE, 
     xaxs = "i",leaflab = "none")

par(mar = c(.5, 0.5, 0, 3.5))
table(substr(colnames(b),1,4))
x=c(4*ones(9,ncol=1),2*ones(9,ncol=1),
    1*ones(5,ncol=1),5*ones(3,ncol=1),
    3*ones(5,ncol=1))

image(as.matrix(x[colInd_Th]),col=Color,
      axes = FALSE)

par(mar = c(.5, .5, 0, 3.5))
image(as.matrix(t(b[rowInd_Th,colInd_Th])), breaks= c(seq(-7, 6, length.out = 21),100),
      col=colorRampPalette(c("blue","white", "red"))(n = 21), axes = FALSE)

ylab <- rownames(b)[rowInd_Th]
xlab <- colnames(b)[colInd_Th]

mtext("Genes", side = 4, adj = 1, 
      line = 1.5, cex = 2, at = .5)

# mtext(xlab, side = 3, las = 2, line = 0.5, cex = .7, at = seq(0, 1, 1/(length(xlab)-1)))
dev.off()





save.image("Selected_Genes.RData")
