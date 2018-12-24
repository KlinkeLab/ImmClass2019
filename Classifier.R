##############################################Info#########################################
############# Title:          Classifier Development
############# Description:    This Script Generates One Classifier for Immune Cells and
#############                 One Classifier for T helper Cells
############# Author:         Arezo Torang
############# Date:           22-Dec-2018
###########################################################################################-



rm(list = ls())

###Importing Data-----------------------------------------------------------




Normalized_TrainData <- read.table("Normalized_TrainData.txt", 
                                   sep = "\t", quote = "",
                                   stringsAsFactors=FALSE)

Normalized_TestData <- read.table("Normalized_TestData.txt", 
                                  sep = "\t", quote = "",
                                  stringsAsFactors=FALSE)

Normalized_TrainData_Th <- read.table("Normalized_TrainData_Th.txt", 
                                      sep = "\t", quote = "",
                                      stringsAsFactors=FALSE)

Normalized_TestData_Th <- read.table("Normalized_TestData_Th.txt", 
                                     sep = "\t", quote = "",
                                     stringsAsFactors=FALSE)




####Classifier for Immune Cells---------------------------------------------
##Class Labels: 1:B Cell, 2:CD4, 3:CD8, 4:Mono, 
##5:Neu, 6:NK, 7:DC, 8:M1, 9:M2, 10:MQ




x <-log2(as.matrix(data.frame(t(Normalized_TrainData[,3:110])))+1)
Test <-log2(as.matrix(data.frame(t(Normalized_TestData[,c(3:62,65:81)])))+1) 

y<-rbind(ones(10, ncol = 1) , 2*ones(10, ncol = 1), 
         3*ones(10, ncol = 1),4*ones(10, ncol = 1),
         5*ones(10, ncol = 1),6*ones(10, ncol = 1),
         1,7*ones(6, ncol = 1),4 ,5 ,6 ,8 ,8, 9, 9,
         2 ,2 ,3 ,3,3,10,10,10,6,ones(3, ncol = 1),
         2*ones(3, ncol = 1) , 3*ones(3, ncol = 1),
         4*ones(3, ncol = 1) , 6*ones(3, ncol = 1),
         7,7,4,4,8,8,9,9,10,10)


y_Test<-rbind( ones(10, ncol = 1), 2*ones(10, ncol = 1),
              3*ones(10, ncol = 1),4*ones(10, ncol = 1),
              5*ones(10, ncol = 1), 6*ones(4, ncol = 1),
              1,7,7,4,5,6,8,9,2,3,3,10,6,1,2,3,4,6,7, 4,
              8,9,10)


##Change Missing Values to -1

x[is.na(x)]<- -1
Test[is.na(Test)]<- -1


##Create Noisy Data to Neglect -1 and Treat It as Missing Values 

X<-x ; Y<-y
for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),length(x)*.1), -1)))
  Y<-rbind(Y,y)
}


##Generation of Classifier (Elastic-Net Logistic Regression)

library(glmnet)
Lasso <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                   lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                   type.multinomial="grouped", nfolds=10)



##Accuracy of Classifier for Testing and Training Samples

prediction <- data.frame(predict(Lasso, newx=data.matrix(x), 
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
      arr.ind = TRUE)

prediction_Test <- data.frame(predict(Lasso, newx=data.matrix(Test), 
                                      interval ="prediction",
                                      type="class", s=0.0001))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy




####Scrambling Data to Creat True Negative Samples and Plot ROC-------------




TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-0*ones(77,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0


for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1
  R=0
  prediction_Test <- data.frame(predict(Lasso, newx=data.matrix(Test),
                                        interval ="prediction",
                                        type="response", s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso, newx=data.matrix(Test),
                                              interval ="prediction",
                                              type="class", s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso, newx=data.matrix(TN_Test),
                                           interval ="prediction",
                                           type="response", s=i))
  
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


##Plot ROC Curve

library(DescTools)

pdf( "ROC_New.pdf", width = 5, height =5 )

Color=c("red","orange","yellow","forestgreen","blue","black","violet")
matplot(1-Specificity,Sensitivity, type="l", col=Color,lwd = 1.5, 
        lty=c(1,2,1,4,5,6,2),xlab = "1-Specificity", 
        ylab = "Sensitivity", ylim = c(.5,1))

Legend=c("lambda=1e-1","lambda=5e-2","lambda=1e-2","lambda=5e-3",
         "lambda=1e-3 ","lambda=5e-4","lambda=1e-4")
legend("bottomright", legend = Legend , col = Color,
       lty = c(1,2,1,4,5,6,2), lwd = 1.5, cex = 0.6)

dev.off()




####Comparig ROC Curves, Lambda=0.0001 has been selected (452 Genes)--------




Index <-coef(Lasso,s=0.0001)[["1"]]@i[-1]
Selected_Genes<-Normalized_TrainData[Index, ]
row.names(Selected_Genes)<-Selected_Genes[,2]

write.table(Selected_Genes, file = "Selected_Genes.txt", 
            sep = "\t", quote = FALSE)


Coef<-coef(Lasso,s=0.0001)
Coef<-data.frame(Coef[["1"]]@x,Coef[["2"]]@x,Coef[["3"]]@x,
                 Coef[["4"]]@x,Coef[["5"]]@x,Coef[["6"]]@x,
                 Coef[["7"]]@x,Coef[["8"]]@x,Coef[["9"]]@x,
                 Coef[["10"]]@x)

colnames(Coef)<-c("B Cell", "CD4", "CD8", "Mono", "Neu",
                  "NK", "DC", "M1", "M2", "MQ")

rownames(Coef)[-1]<-as.character(Normalized_TrainData[Index,2])
write.table(Coef, file = "Coef.txt", sep = "\t", quote = FALSE)




####Classifier for T helper Cells-------------------------------------------
####Class Labels: 1:Th1, 2:Th2, 3:Th17, 4:Th0, 5:iTreg




x <-log2(as.matrix(data.frame(t(Normalized_TrainData_Th[,3:33])))+1)
Test <-log2(as.matrix(data.frame(t(Normalized_TestData_Th[,3:19])))+1) 

y<-rbind(1,1,2,2,1,1,1,2,2,2,3,3,3,4*ones(9, ncol = 1),5*ones(9, ncol = 1))
y_Test<-rbind(1,2,1,2,3,4*ones(6, ncol = 1),5*ones(6, ncol = 1))


##Change Missing Values to -1

x[is.na(x)]<- -1
Test[is.na(Test)]<- -1


##Create Noisy Data to Neglect -1 and Treat It as Missing Values 

X<-x ; Y<-y
for (i in 1:100) {
  X<-as.matrix(rbind(X,replace(x, sample(1:length(x),
                                        length(x)*.1), -1)))
  Y<-rbind(Y,y)
}


##Generation of Classifier (Elastic-Net Logistic Regression)

Lasso_Th <- cv.glmnet(X, Y, family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                      type.multinomial="grouped", nfolds=10) 



##Accuracy of Classifier for Testing and Training Samples

prediction <- data.frame(predict(Lasso_Th, newx=data.matrix(x), 
                                 interval ="prediction",
                                 type="response", s=0.05))

which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)),arr.ind = TRUE)
prediction_Test <- data.frame(predict(Lasso_Th, newx=data.matrix(Test), 
                                      interval ="prediction",
                                      type="class", s=0.05))

sum(y_Test==as.matrix(prediction_Test))/length(y_Test) #accuracy




#### Scrambling Data to Creat True Negative Samples and Plot ROC------------




TN_Test<-as.matrix(Test[sample(nrow(Test)),])
TN_Test<-as.matrix(TN_Test[,sample(ncol(Test))])
TN_y<-0*ones(17,ncol=1)

Sensitivity=matrix(nrow = 100,ncol = 7)
Specificity=matrix(nrow = 100,ncol = 7)
C=0


for (i in c(.1,.05,.01,.005,.001,.0005,.0001)) {
  C=C+1; R=0
  prediction_Test <- data.frame(predict(Lasso_Th, newx=data.matrix(Test), 
                                        interval ="prediction",
                                        type="response", s=i))
  
  prediction_Test_Class <- data.frame(predict(Lasso_Th, newx=data.matrix(Test), 
                                              interval ="prediction",
                                              type="class", s=i))
  
  prediction_TN_Test <- data.frame(predict(Lasso_Th, newx=data.matrix(TN_Test), 
                                           interval ="prediction",
                                           type="response", s=i))
  
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



##Plot ROC Curve

pdf( "ROC_Th_New.pdf", width = 5, height =5 )

Color=c("red","orange","yellow","forestgreen","blue","black","violet")
matplot(1-Specificity,Sensitivity, type="l", col=Color,lwd = 1.5,
        lty=c(1,2,1,4,5,6,2),xlab = "1-Specificity",
        ylab = "Sensitivity")

legend("bottomright", legend = Legend ,col=Color,
       lty = c(1,2,1,4,5,6,2), lwd=1.5, cex=0.6)


dev.off()




####Comparig ROC Curves, Lambda=0.05 has been selected (72 Genes)-----------




Index_Th <-coef(Lasso_Th,s=.05)[["1"]]@i[-1]
Selected_Genes_Th<-Normalized_TrainData_Th[Index_Th, ]
row.names(Selected_Genes_Th)<-Selected_Genes_Th[,2]

write.table(Selected_Genes_Th, file = "Selected_Genes_Th.txt", 
            sep = "\t", quote = FALSE)


Coef_Th<-coef(Lasso_Th,s=0.05)
Coef_Th<-data.frame(Coef_Th[["1"]]@x,Coef_Th[["2"]]@x,
                    Coef_Th[["3"]]@x,Coef_Th[["4"]]@x,
                    Coef_Th[["5"]]@x)

colnames(Coef_Th)<-c("Th1", "Th2", "Th17", "Th0", "iTreg")
rownames(Coef_Th)[-1]<-as.character(Normalized_TrainData_Th[Index_Th,2])

write.table(Coef_Th, file = "Coef_Th.txt", sep = "\t", quote = FALSE)
