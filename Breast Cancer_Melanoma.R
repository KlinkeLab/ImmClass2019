##############################################Info#########################################
############# Title:          Breast Cancer_Melanoma Analysis
############# Description:    This Script Investigate the Samples in Breast Cancer Dataset 
#############                 and Melanoma Dataset 
############# Author:         Arezo Torang
############# Date:           10-Apr-2019
###########################################################################################-



rm(list = ls())

###Importing Data-----------------------------------------------------------



load("Gene Selection.RData")

TrainData <- read.table("TrainData.txt")

Normalized_TrainData <- read.table("Normalized_TrainData.txt", 
                                   sep = "\t", quote = "",
                                   stringsAsFactors=FALSE)

Length_GC <- read.table("Length_GC.txt", sep = "\t", quote = "")

library("Matrix")
load("Melanoma SSRNAseq.RData")

####Breast Cancer-------------
#Data was obtained from: (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688)



BC_Total<-read.table("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",
                     sep = "\t", quote = "")

BC_Inf<-read.table("GSE75688_final_sample_information.txt",
                   sep = "\t", quote = "")

Samples<-t(BC_Total[1,])

BC_BCell<-BC_Inf[which(BC_Inf$V5=="Bcell"),1] #83 Sample
BC_BCell<-BC_Total[,which(Samples %in% BC_BCell)]

BC_TCell<-BC_Inf[which(BC_Inf$V5=="Tcell"),1] #54 Sample
BC_TCell<-BC_Total[,which(Samples %in% BC_TCell)]

BC_Myeloid<-BC_Inf[which(BC_Inf$V5=="Myeloid"),1] #38 Sample
BC_Myeloid<-BC_Total[,which(Samples %in% BC_Myeloid)]

BC_Bulk<-BC_Inf[which(BC_Inf$V5=="Immune"),1] #3 Sample
BC_Bulk<-BC_Total[,which(Samples %in% BC_Bulk)]

BC_Samples<-cbind(BC_Total[,1:2],BC_BCell,BC_TCell,
                  BC_Myeloid,BC_Bulk)

y=rbind(ones(83,ncol=1),2*ones(54,ncol=1),
        11*ones(38,ncol=1),12*ones(3,ncol=1))

a<-t(rbind(0,0,y))            #1:Bcell, 2:Tcell, 11:Myeloid, 12:Bulk
colnames(a)<-colnames(BC_Samples)
BC_Samples<-rbind(a,BC_Samples)

write.table(BC_Samples,"BC_Samples.txt")



##Normalization of Breast Cancer Samples------


library(EDASeq)

colnames(BC_Samples)[2]<-"ID"
Genes<-TrainData[,1:2]
a<-merge(Genes,BC_Samples,all.x = TRUE)
a<-a[,-3]

for (i in 3:ncol(a)) {
  a[,i]<-(as.numeric(a[,i]))
}

a<-a[,2:ncol(a)]
BC_Data<-merge(TrainData,a,all.x = TRUE)
Repeat_Free_Index<-seq_along(BC_Data[,1])[!duplicated(BC_Data[,1])]
BC_Data<-BC_Data[Repeat_Free_Index,]

Data_Length_GC <- merge(BC_Data,Length_GC)

dataWithin<-withinLaneNormalization(x=as.matrix(Data_Length_GC[,3:(ncol(Data_Length_GC)-2)]),
                                    y=Data_Length_GC[,ncol(Data_Length_GC)], offset=FALSE,
                                    round=TRUE,which = "full")

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

boxplot(log2(dataNorm+1))

Normalized_BC <- data.frame(Data_Length_GC[,1:2],
                            dataNorm[,109:ncol(dataNorm)])

Genes<-data.frame(Normalized_TrainData[,1])
colnames(Genes)<-"Genes"

Normalized_BC <- unique(merge(Genes,Normalized_BC,
                              all.x = TRUE))

write.table(Normalized_BC,file = "Normalized_BC.txt")


##Analysis over Breast Cancer Data------


y_Test=rbind(ones(83,ncol=1),2*ones(54,ncol=1),
             11*ones(38,ncol=1),12*ones(3,ncol=1))

Test <-log2(as.matrix(data.frame(t(Normalized_BC[,-c(1,2)])))+1)

Test[is.na(Test)]<- -1

library(glmnet)
prediction <- data.frame(predict(Lasso, newx=data.matrix(Test),
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

Class<-which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)),
             arr.ind = TRUE)

Class<-Class[order(Class[,1]),]

prediction_Class <-predict(Lasso, newx=data.matrix(Test),
                           interval ="prediction",type="class",
                           s=0.0001)


Correct=0
for (i in min(which(y_Test %in% 1)):max(which(y_Test %in% 1))){
  if(prediction_Class[i,1]==1){
    Correct=Correct+1
  }
}

for (i in min(which(y_Test %in% 2)):max(which(y_Test %in% 2))){
  if(prediction_Class[i,1]==2|prediction_Class[i,1]==3){
    Correct=Correct+1
  }
}

Total_Accuracy<-Correct/max(which(y_Test %in% 2))
Total_Accuracy #0.8759124

TCell=0
BCell=0
M=0
B=0
j=0
k=0

FN_TCell=c()
FN_BCell=c()
FN_M=c()
FN_B=c()

TP_TCell=c()
TP_BCell=c()
TP_M=c()
TP_B=c()

for (i in min(which(y_Test %in% 2)):max(which(y_Test %in% 2))){
  if(prediction_Class[i,1]==2|prediction_Class[i,1]==3){
    TCell=TCell+1
    k<-k+1
    TP_TCell[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_TCell[j]<-prediction_Class[i,1]}
}

j=0
k=0

for (i in min(which(y_Test %in% 1)):max(which(y_Test %in% 1))){
  if(prediction_Class[i,1]==1){
    BCell=BCell+1
    k<-k+1
    TP_BCell[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_BCell[j]<-prediction_Class[i,1]}
}

j=0
k=0

for (i in min(which(y_Test %in% 11)):max(which(y_Test %in% 11))){
  if(prediction_Class[i,1]==0){
    M=M+1
    k<-k+1
    TP_M[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_M[j]<-prediction_Class[i,1]}
}

j=0
k=0

for (i in min(which(y_Test %in% 12)):max(which(y_Test %in% 12))){
  if(prediction_Class[i,1]==0){
    B=B+1
    k<-k+1
    TP_B[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_B[j]<-prediction_Class[i,1]}
}


TCell/length(which(y_Test%in%2)) #0.9444444
BCell/length(which(y_Test%in%1)) #0.8313253
M/length(which(y_Test%in%11))    #0
B/length(which(y_Test%in%12))    #0


##Pie Chart 


library(tidyverse)
library(plotrix)

BCell<-table(prediction_Class[min(which(y_Test %in% 1)):max(which(y_Test %in% 1)),1])
TCell<-table(prediction_Class[min(which(y_Test %in% 2)):max(which(y_Test %in% 2)),1])
M<-table(prediction_Class[min(which(y_Test %in% 11)):max(which(y_Test %in% 11)),1])
B<-table(prediction_Class[min(which(y_Test %in% 12)):max(which(y_Test %in% 12)),1])

Col1<-c(rep("T Cell",dim(TCell)), rep("B Cell",dim(BCell)), 
        rep("Myeloid",dim(M)), rep("Bulk",dim(B)))

COl2<-c("CD8","CD4","DC","NK",
        "B Cell","CD8","CD4","DC",
        "CD8","MQ","M1","CD4","DC","M2","B Cell","Neu",
        "CD8","CD4")

COl3<-c(TCell[[2]],TCell[[1]],TCell[[4]],TCell[[3]],
        BCell[[1]],BCell[[3]],BCell[[2]],BCell[[4]],
        M[[4]],M[[2]],M[[7]],M[[6]],M[[3]],M[[8]],M[[1]],M[[5]],
        B[[2]],B[[1]])

BC_Result<-data.frame(Col1,COl2,COl3)

colnames(BC_Result)<-c("True_Label","Prediction","Number")

Labeled_Result<-data.frame(c("T Cell", "B Cell", "Myeloid","Bulk"),
                           c(length(which(y_Test%in%2)),
                             length(which(y_Test%in%1)),
                             length(which(y_Test%in%11)),
                             length(which(y_Test%in%12))))

colnames(Labeled_Result)<-c("Cell","Number")


pdf("pie_BC.pdf",7,7)

plot.new()

par(new = TRUE)
pie(BC_Result$Number, border = NA, radius = 1, labels = NA,
    col = c("firebrick2","firebrick3","plum2","yellow1",
            "royalblue1","firebrick2","firebrick3","plum2",
            "firebrick2","gray","black","firebrick3","plum2","green","royalblue1","gold",
            "firebrick2","firebrick3"), cex=1)

x=c()
x[1]<-BC_Result$Number[1]/2

for(i in 2:nrow(BC_Result)){
  x[i]=sum(BC_Result$Number[1:i-1])+
    BC_Result$Number[i]/2
}

lab<-c("CD8","CD4","DC"," ","B Cell","CD8","CD4",
       " ","CD8",expression(paste("M",phi)),"M1",
       " "," ","M2"," "," "," "," ")

pie.labels(x=0,y=0,(x/length(y_Test))*2*pi,
           labels= lab,radius=1.05,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=1.5)


par(new = TRUE)
pie(Labeled_Result$Number, border = NA,
    radius = .7, labels = NA,
    col = c("darkred","royalblue3",
            "darkorange2","yellow3"))

x=c()
x[1]<-Labeled_Result$Number[1]/2

for(i in 2:nrow(Labeled_Result)){
  x[i]=sum(Labeled_Result$Number[1:i-1])+
    Labeled_Result$Number[i]/2
}

lab<-c("T Cell", "B Cell", "Myeloid","Bulk")
pie.labels(x=0,y=0,(x/length(y_Test))*2*pi,
           labels= lab,radius=.35,bg=NA,
           border=NA, minangle=NA,boxed=FALSE,
           cex=1.2, col="beige")

dev.off()


plot(x=prediction$CD8[1:83],y=prediction$`B Cell`[1:83],
     xlim=c(0,1),ylim=c(0,1))

lines(c(0,1),c(0,1),col="blue", lwd=2, lty=2)









####Melanoma-------
#Data was obtained from: (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056)



##Normalization of Melanoma Samples------


#Annotated Melanoma Samples (Label = 1)

library(EDASeq)

Melanoma_Labeled<-data.frame(matrix(nrow=nrow(Melanoma_SC),ncol=1))

ind=c()
for (i in 2:ncol(Melanoma_SC)) {
  if(Melanoma_SC[3,i]!=0){
    Melanoma_Labeled<-data.frame(Melanoma_Labeled,Melanoma_SC[,i])
    ind<-c(ind,i)
    print(i)
  }
}

Melanoma_Labeled<-Melanoma_Labeled[,2:ncol(Melanoma_Labeled)]
colnames(Melanoma_Labeled)<-colnames(Melanoma_SC)[ind]
Melanoma_Labeled<-cbind(Melanoma_SC[,1],Melanoma_Labeled)
colnames(Melanoma_Labeled)[1]<-"ID"
Genes<-TrainData[,1:2]

a<-merge(Genes,Melanoma_Labeled,all.x = TRUE)
a<-a[,2:ncol(a)]
Melanoma_Data<-merge(TrainData,a,all.x = TRUE)

Repeat_Free_Index<-seq_along(Melanoma_Data[,1])[!duplicated(Melanoma_Data[,1])]
Melanoma_Data<-Melanoma_Data[Repeat_Free_Index,]

Length_GC <- read.table("Length_GC.txt", sep = "\t", quote = "")
Data_Length_GC <- merge(Melanoma_Data,Length_GC)

dataWithin <- withinLaneNormalization(x=as.matrix(Data_Length_GC[,3:(ncol(Data_Length_GC)-2)]),
                                      y=Data_Length_GC[,ncol(Data_Length_GC)], offset=FALSE,
                                      round=TRUE,which = "full")

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

boxplot(log2(dataNorm+1))

Normalized_Melanoma <- data.frame(Data_Length_GC[,1:2],
                                  dataNorm[,109:ncol(dataNorm)])

Genes<-data.frame(Normalized_TrainData[,1])
colnames(Genes)<-"Genes"

Normalized_Melanoma_Labeled <- unique(merge(Genes,Normalized_Melanoma,all.x = TRUE))


ind<-as.matrix(Melanoma_Labeled[3,-1])
colnames(ind)<-c()
Main<-Normalized_Melanoma_Labeled[,-c(1,2)]
a<-Melanoma_Labeled[,-1]
Main<-Main[ , order(ind)]; a<-a[ , order(ind)]
ind<-ind[order(ind)]

Main<-Main[,!ind %in% 4:5]
a<-a[,!ind %in% 4:5]
; ind<-ind[!ind %in% 4:5]

Normalized_Melanoma_Labeled<-data.frame(Normalized_Melanoma_Labeled[,1:2],Main)
Melanoma_Labeled<-data.frame(Melanoma_Labeled[,1],a)

write.table(Normalized_Melanoma_Labeled,file = "Normalized_Melanoma_Labeled.txt")
write.table(Melanoma_Labeled,"Melanoma_Labeled.txt")



#Unresolved Melanoma Samples (Label = 0)



Melanoma_Unlabeled<-data.frame(matrix(nrow=nrow(Melanoma_SC),ncol=1))
ind=c()

for (i in 2:ncol(Melanoma_SC)) {
  if(Melanoma_SC[3,i]==0){
    Melanoma_Unlabeled<-data.frame(Melanoma_Unlabeled,Melanoma_SC[,i])
    ind<-c(ind,i)
    print(i)
  }
}

Melanoma_Unlabeled<-Melanoma_Unlabeled[,2:ncol(Melanoma_Unlabeled)]
colnames(Melanoma_Unlabeled)<-colnames(Melanoma_SC)[ind]
Melanoma_Unlabeled<-cbind(Melanoma_SC[,1],Melanoma_Unlabeled)

colnames(Melanoma_Unlabeled)[1]<-"ID"
write.table(Melanoma_Unlabeled,file="Melanoma_Unlabeled.txt")

Genes<-TrainData[,1:2]
a<-merge(Genes,Melanoma_Unlabeled,all.x = TRUE)
a<-a[,2:ncol(a)]

Melanoma_Data<-merge(TrainData,a,all.x = TRUE)

Repeat_Free_Index<-seq_along(Melanoma_Data[,1])[!duplicated(Melanoma_Data[,1])]
Melanoma_Data<-Melanoma_Data[Repeat_Free_Index,]

Length_GC <- read.table("Length_GC.txt", sep = "\t", quote = "")
Data_Length_GC <- merge(Melanoma_Data,Length_GC)

dataWithin <- withinLaneNormalization(x=as.matrix(Data_Length_GC[,3:(ncol(Data_Length_GC)-2)]),
                                      y=Data_Length_GC[,ncol(Data_Length_GC)], 
                                      offset=FALSE,round=TRUE, which = "full")

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

Normalized_Melanoma_Unlabeled <- data.frame(Data_Length_GC[,1:2],
                                            dataNorm[,109:ncol(dataNorm)])

Genes<-data.frame(Normalized_TrainData[,1])
colnames(Genes)<-"Genes"

Normalized_Melanoma_Unlabeled <- unique(merge(Genes,Normalized_Melanoma_Unlabeled,
                                              all.x = TRUE))

write.table(Normalized_Melanoma_Unlabeled,file = "Normalized_Melanoma_Unlabeled.txt")


##Analysis over Labeled Melanoma Data----


y_Test<-as.matrix(Melanoma_Labeled[3,-1])
Test <-log2(as.matrix(data.frame(t(Normalized_Melanoma_Labeled[,-c(1,2)])))+1)
Test[is.na(Test)]<- -1

prediction <- data.frame(predict(Lasso, newx=data.matrix(Test),
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono",
                        "Neu","NK","DC","M1","M2","MQ")

Class<-which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
             arr.ind = TRUE)

prediction_Class <-predict(Lasso, newx=data.matrix(Test),
                           interval ="prediction",
                           type="class", s=0.0001)

Correct=0
for (i in min(which(y_Test %in% 1)):max(which(y_Test %in% 1))){
  if(prediction_Class[i,1]==2|prediction_Class[i,1]==3){
    Correct=Correct+1
  }
}

for (i in min(which(y_Test %in% 2)):max(which(y_Test %in% 2))){
  if(prediction_Class[i,1]==1){
    Correct=Correct+1
  }
}

for (i in min(which(y_Test %in% 3)):max(which(y_Test %in% 3))){
  if(prediction_Class[i,1]==10){
    Correct=Correct+1
  }
}

for (i in min(which(y_Test %in% 6)):max(which(y_Test %in% 6))){
  if(prediction_Class[i,1]==6){
    Correct=Correct+1
  }
}


Total_Accuracy<-Correct/length(y_Test)
Total_Accuracy #0.9619703

TCell=0
BCell=0
MQ=0
NK=0

j=0
k=0

FN_TCell=c()
FN_BCell=c()
FN_MQ=c()
FN_NK=c()

TP_TCell=c()
TP_BCell=c()
TP_MQ=c()
TP_NK=c()

for (i in min(which(y_Test %in% 1)):max(which(y_Test %in% 1))){
  if(prediction_Class[i,1]==2|prediction_Class[i,1]==3){
    TCell=TCell+1
    k<-k+1
    TP_TCell[k]<-prediction_Class[i,1]
  }else{j<-j+1;FN_TCell[j]<-prediction_Class[i,1]}
}

j=0
k=0
for (i in min(which(y_Test %in% 2)):max(which(y_Test %in% 2))){
  if(prediction_Class[i,1]==1){
    BCell=BCell+1
    k<-k+1
    TP_BCell[k]<-prediction_Class[i,1]
  }else{
    j<-j+1
    FN_BCell[j]<-prediction_Class[i,1]
  }
}

j=0
k=0
for (i in min(which(y_Test %in% 3)):max(which(y_Test %in% 3))){
  if(prediction_Class[i,1]==10){
    MQ=MQ+1
    k<-k+1
    TP_MQ[k]<-prediction_Class[i,1]
  }else{
    j<-j+1
    FN_MQ[j]<-prediction_Class[i,1]
  }
}

j=0
k=0
for (i in min(which(y_Test %in% 6)):max(which(y_Test %in% 6))){
  if(prediction_Class[i,1]==6){
    NK=NK+1
    k<-k+1
    TP_NK[k]<-prediction_Class[i,1]
  }else{
    j<-j+1
    FN_NK[j]<-prediction_Class[i,1]
  }
}

#Accuracy for Each Cell Type
TCell/length(which(y_Test%in%1)) #0.9758221
BCell/length(which(y_Test%in%2)) #0.9805825
MQ/length(which(y_Test%in%3))    #0.8412698
NK/length(which(y_Test%in%6))    #0.5192308



##Pie Chart for Labeled Cell Types (1=T,2=B,3=Macro,6=NK)

library(tidyverse)
library(plotrix)

Col1<-c(rep("T Cell",5), rep("B Cell",4),
        rep("MQ",4), rep("NK",2))

COl2<-c("CD8","CD4","B Cell","MQ","NK",
        "B Cell","CD8","MQ","DC",
        "MQ","B Cell","CD8","Mono",
        "NK","CD8")

COl3<-c(table(TP_TCell)[[2]],table(TP_TCell)[[1]],table(FN_TCell)[[1]],table(FN_TCell)[[2]],table(FN_TCell)[[3]],
        length(TP_BCell),table(FN_BCell)[[2]],table(FN_BCell)[[1]],table(FN_BCell)[[3]],
        length(TP_MQ),table(FN_MQ)[[1]],table(FN_MQ)[[2]],table(FN_MQ)[[3]],
        length(TP_NK),table(FN_NK)[[1]])

Melanoma_Result<-data.frame(Col1,COl2,COl3)

colnames(Melanoma_Result)<-c("True_Label","Prediction","Number")

Labeled_Result<-data.frame(c("T Cell", "B Cell", "MQ","NK"),
                           c(2068,515,126,52))

colnames(Labeled_Result)<-c("Cell","Number")


pdf("pie_Labeled.pdf",7,7)

plot.new()
par(new = TRUE)
pie(Melanoma_Result$Number, border = NA, radius = 1, labels = NA,
    col = c("firebrick2","firebrick3","royalblue1","darkorange","yellow1",
            "royalblue1","firebrick2","darkorange","grey",
            "darkorange","royalblue1","firebrick2","black",
            "yellow1","firebrick2"), cex=1.5)

x=c()
x[1]<-Melanoma_Result$Number[1]/2
for(i in 2:nrow(Melanoma_Result)){
  x[i]=sum(Melanoma_Result$Number[1:i-1])+
    Melanoma_Result$Number[i]/2
}

lab<-c("CD8","CD4"," "," "," ","B Cell"," "," ",
       " ",expression(paste("M",phi))," "," ",
       " ","NK"," ")

pie.labels(x=0,y=0,(x/length(y_Test))*2*pi,
           labels= lab,radius=1.05,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=1)


par(new = TRUE)
pie(Labeled_Result$Number, border = NA,
    radius = .7, labels = NA,
    col = c("darkred","royalblue3",
            "darkorange2","yellow3"))

x=c()
x[1]<-Labeled_Result$Number[1]/2

for(i in 2:nrow(Labeled_Result)){
  x[i]=sum(Labeled_Result$Number[1:i-1])+
    Labeled_Result$Number[i]/2
}

lab<-c("T Cell", "B Cell", expression(paste("M",phi)),"NK")

pie.labels(x=0.1,y=0,(x/length(y_Test))*2*pi,
           labels= lab,radius=.35,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=1.2, col="beige")

dev.off()



##Analysis over Unlabeled Melanoma Data---------
#(malignant(1=no,2=yes,0=unresolved))


y_Test<-as.matrix(Melanoma_Unlabeled[2,2:ncol(Melanoma_Unlabeled)])

colnames(y_Test)<-c()
Main<-Normalized_Melanoma_Unlabeled[,-c(1,2)]
Main<-Main[ , order(y_Test)]
y_Test<-y_Test[order(y_Test)]

Test <-log2(as.matrix(data.frame(t(Main)))+1)
Test[is.na(Test)]<- -1

prediction <- data.frame(predict(Lasso, newx=data.matrix(Test), 
                                 interval ="prediction",
                                 type="response", s=0.0001))

colnames(prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                        "NK","DC","M1","M2","MQ")

Class<-which(as.matrix(prediction)==rowMaxs(as.matrix(prediction)) ,
             arr.ind = TRUE)

Class<-Class[order(Class[,1]),]

Nonmalignant_Samples<-Main[,min(which(y_Test %in% 1)):max(which(y_Test %in% 1))]

Names<-data.frame(colnames(Nonmalignant_Samples))

Immune_Cells<-Nonmalignant_Samples[,-c(65,66,72,73,77,80,82,87,99,109,124,128,130,
                                       132:135,146:154,157:163,165,167,169,170,177,
                                       182,188,189,195,198,204,208,214,217:220,227,
                                       238:248,250,252,259,263,265,267,273,278,279,282,
                                       283,285,294,301:304,308:310,312,314,319,320,
                                       381:415)]

Test_Immune<-log2(as.matrix(data.frame(t(Immune_Cells)))+1)
Test_Immune[is.na(Test_Immune)]<- -1

Immune_Cells_Prediction<- data.frame(predict(Lasso, newx=data.matrix(Test_Immune), 
                                             interval ="prediction",type="response",
                                             s=0.0001))

colnames(Immune_Cells_Prediction)<-c("B Cell","CD4","CD8","Mono","Neu",
                                     "NK","DC","M1","M2","MQ")

Unresolved<-prediction[min(which(y_Test %in% 0)):max(which(y_Test %in% 0)),]
Nonmalignant<-prediction[min(which(y_Test %in% 1)):max(which(y_Test %in% 1)),]
Malignant<-prediction[min(which(y_Test %in% 2)):max(which(y_Test %in% 2)),]

Class_Immune<-which(as.matrix(Immune_Cells_Prediction)==
                      rowMaxs(as.matrix(Immune_Cells_Prediction)) ,
                    arr.ind = TRUE)

Class_Immune<-Class_Immune[order(Class_Immune[,1]),]
I<-data.frame(table(Class_Immune[,2]))

Class_Unresolved<-Class[min(which(y_Test %in% 0)):max(which(y_Test %in% 0)),]
Class_Unresolved<-Class_Unresolved[order(Class_Unresolved[,2]),]
U<-data.frame(table(Class_Unresolved[,2]))

Class_Nonmalignant<-Class[min(which(y_Test %in% 1)):max(which(y_Test %in% 1)),]
Class_Nonmalignant<-Class_Nonmalignant[order(Class_Nonmalignant[,2]),]
N<-data.frame(table(Class_Nonmalignant[,2]))

Class_Malignant<-Class[min(which(y_Test %in% 2)):max(which(y_Test %in% 2)),]
Class_Malignant<-Class_Malignant[order(Class_Malignant[,2]),]
M<-data.frame(table(Class_Malignant[,2]))



pdf("pie_Unlabeled_Immune_Cells.pdf",5,5)
pie(N[c(3,2,1,6,4,5),2], border = NA,radius = 1,cex=1.2,
    col = c("royalblue3","red","darkred","yellow",
            "green","darkorange2")[c(3,2,1,6,4,5)],
    labels=c("B Cell","CD4","CD8","NK","DC",
             expression(paste("M",phi)))[c(3,2,1,6,4,5)])

dev.off()


pdf("pie_Unlabeled.pdf",5,5)

par(mfrow=c(1,3))
pie(U[c(2,1,5,3,4),2] ,border = NA,radius = 1,cex=1.2,
    col = c("royalblue3","darkred","yellow",
            "green","darkorange2")[c(2,1,5,3,4)],
    labels=c("B Cell","CD8","NK","DC",
             expression(paste("M",phi)))[c(2,1,5,3,4)])

pie(N[c(3,2,1,6,4,5),2], border = NA,radius = 1,cex=1.2,
    col = c("royalblue3","red","darkred","yellow",
            "green","darkorange2")[c(3,2,1,6,4,5)],
    labels=c("B Cell","CD4","CD8","NK","DC",
             expression(paste("M",phi)))[c(3,2,1,6,4,5)])

pie(M[c(3,2,1,7,5,6,4),2], border = NA,radius = 1,cex=1.2,
    col = c("royalblue3","red","darkred","grey","yellow",
            "green","darkorange2")[c(3,2,1,7,5,6,4)],
    labels=c("B Cell","CD4","CD8","Mono","NK","DC",
             expression(paste("M",phi)))[c(3,2,1,7,5,6,4)])

dev.off()



##Analysis of NK Samples to See the Similarities with T Cells --------

X <-log2(as.matrix(data.frame(t(BigTrain[,-1])))+1)

PCA <- prcomp(t(X), retx = TRUE, scale. = FALSE)
PCA.data <- as.matrix(t(X)) %*% PCA$rotation

PCA1 <- prcomp(X, retx = TRUE, scale. = FALSE)
PCA1.data <- as.matrix(X) %*% PCA1$rotation

plot(PCA.data[,c(1,2)], type = "n")
text(PCA.data[,c(1,2)], labels =BigTrain[,1] )

lines(c(5, 5), c(-10, 10), col="blue", lwd=2, lty=2)
lines(c(-100, 100), c(-2, -2), col="blue", lwd=2, lty=2)
lines(c(-100, 100), c(1.5, 1.5), col="blue", lwd=2, lty=2)

KeepGene1 <- PCA.data[,1] >5
KeepGene2 <- PCA.data[,2] < -2
KeepGene3 <- PCA.data[,2] > 1.5

KeepGene<-c(which(KeepGene1),which(KeepGene2),
            which(KeepGene3))

KeepGene<-unique(KeepGene)
KeepGene<-KeepGene[order(KeepGene)]

plot(PCA.data[KeepGene,c(3,2)], type = "n")
text(PCA.data[KeepGene,c(3,2)], labels=BigTrain[KeepGene,1],
     cex = 0.8)

BigTrain[KeepGene,1]

#PCA Result Enrichment Link: http://amp.pharm.mssm.edu/Enrichr/enrich?dataset=6gipm



##Kde2d contour Plots for CD8 and NK


library(MASS)
KdCD8 <- kde2d(prediction$NK[min(which(y_Test==1)):max(which(y_Test==1))],
               prediction$CD8[min(which(y_Test==1)):max(which(y_Test==1))],
               lims = c(-.1,1,-.1,1), n = 300)

KdBCell <- kde2d(prediction$NK[min(which(y_Test==2)):max(which(y_Test==2))],
                 prediction$CD8[min(which(y_Test==2)):max(which(y_Test==2))],
                 lims = c(-.1,1,-.1,1), n = 300)

KdMQ <- kde2d(prediction$NK[min(which(y_Test==3)):max(which(y_Test==3))],
              prediction$CD8[min(which(y_Test==3)):max(which(y_Test==3))],
              lims = c(-.1,1,-.1,1), n = 300)

KdNK <- kde2d(prediction$NK[min(which(y_Test==6)):max(which(y_Test==6))], 
              prediction$CD8[min(which(y_Test==6)):max(which(y_Test==6))],
              n = 120)

Kd <- kde2d(prediction$NK[min(which(y_Test!=6)):max(which(y_Test!=6))],
            prediction$CD8[min(which(y_Test!=6)):max(which(y_Test!=6))],
            n = 120)


pdf("Contour.pdf",5,5)

contour(KdCD8,xlim=c(0,.9),ylim = c(0,.9),drawlabels = FALSE,
        col=colorRampPalette(c("#EBEBEB", "black"))(n = 14),
        lwd = 1.2,xlab="NK Probability",ylab="CD8 Probability")

par(new=TRUE)
contour(KdMQ, xlab = NA, ylab = NA,xlim=c(0,.9),ylim = c(0,.9),
        col=colorRampPalette(c("#D4D4FF", "blue"))(n = 7),
        drawlabels = FALSE,lwd = 1.2)

par(new=TRUE)
contour(KdBCell, xlab = NA, ylab = NA,xlim=c(0,.9),ylim = c(0,.9),
        col=colorRampPalette(c("#E5FFE5", "forestgreen"))(n = 11),
        drawlabels = FALSE,lwd = .7)

points(x=prediction$NK[min(which(y_Test==6)):max(which(y_Test==6))],
       y=prediction$CD8[min(which(y_Test==6)):max(which(y_Test==6))],
       type = "p",xlim=c(0,1),ylim=c(0,.9),col="red",pch="*")

legend("bottomright",  c('T Cell',expression(paste('M',phi)),
                         'B Cell','NK Cell'),
       col=c("#363636","blue","forestgreen", 'red'), 
       lty=c(1,1,1,NA),pch=c(NA,NA,NA,'*'), cex=0.8,
       lwd = 1.2)


dev.off()



##T Helper Normalization for CD4 Samples-------------


y_Test<-as.matrix(Melanoma_Labeled[3,-1])
Test_TCell <-log2(as.matrix(data.frame(t(Normalized_Melanoma_Labeled
                                         [,-c(1,2)][,which(y_Test %in% 1)]
                                         )))+1)

Test_TCell[is.na(Test_TCell)]<- -1

prediction_TCell <- data.frame(predict(Lasso, newx=data.matrix(Test_TCell),
                                       interval ="prediction",type="response",
                                       s=0.0001))

colnames(prediction_TCell)<-c("B Cell","CD4","CD8","Mono","Neu",
                              "NK","DC","M1","M2","MQ")

Class_TCell<-which(as.matrix(prediction_TCell)==rowMaxs(as.matrix(prediction_TCell)),
                   arr.ind = TRUE)

CD4_Names<-rownames(Class_TCell[which(as.vector(Class_TCell[,2]) %in% 2),])
CD4_Names<-CD4_Names[order(CD4_Names)]
write.table(CD4_Names,"CD4_Names.txt")


library(EDASeq)

CD4_Names<-(read.table("CD4_Names.txt"))
CD4_Names<-as.vector(CD4_Names[1:nrow(CD4_Names),1])

TrainData_Th <- read.table("TrainData_Th.txt")
Genes<-TrainData_Th[,1:2]
CD4<-Melanoma_Labeled[CD4_Names]
CD4<-data.frame(Melanoma_Labeled[,1],CD4)
colnames(CD4)[1]<-"ID"

a<-merge(Genes,CD4,all.x = TRUE)
a<-a[,2:ncol(a)]

CD4_Data<-merge(TrainData_Th,a,all.x = TRUE)
Repeat_Free_Index<-seq_along(CD4_Data[,1])[!duplicated(CD4_Data[,1])]

CD4_Data<-CD4_Data[Repeat_Free_Index,]

Length_GC_Th <- read.table("Length_GC_Th.txt", 
                           sep = "\t", quote = "")

Data_Length_GC <- merge(CD4_Data,Length_GC_Th)

dataWithin <- withinLaneNormalization(x=as.matrix(Data_Length_GC[,3:(ncol(Data_Length_GC)-2)]),
                                      y=Data_Length_GC[,ncol(Data_Length_GC)], 
                                      offset=FALSE,round=TRUE,which = "full")

dataNorm <- betweenLaneNormalization(dataWithin, which = "full")

Normalized_CD4_Th <- data.frame(Data_Length_GC[,1:2],
                                dataNorm[,32:ncol(dataNorm)])

Genes<-data.frame(Normalized_TrainData_Th[,1])
colnames(Genes)<-"Genes"

Normalized_CD4_Th <- unique(merge(Genes,Normalized_CD4_Th,all.x = TRUE))
write.table(Normalized_CD4_Th,file = "Normalized_CD4_Th.txt")



##Plot T helper Classification of CD4 For Each Patient--------------------

Test_CD4 <-log2(as.matrix(data.frame(t(Normalized_CD4_Th[,-c(1,2)])))+1)
Test_CD4[is.na(Test_CD4)]<- -1

prediction_CD4 <- data.frame(predict(Lasso_Th, newx=data.matrix(Test_CD4), 
                                     interval ="prediction",type="response",
                                     s=0.05))

colnames(prediction_CD4)<-c("Th1", "Th2", "Th17", "Th0", "iTreg")

Class_CD4<-which(as.matrix(prediction_CD4)==rowMaxs(as.matrix(prediction_CD4)) ,
                 arr.ind = TRUE)


table(substr(rownames(Class_CD4),3,4))

Mel53_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="53"),2])
Mel58_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="58"),2])
Mel59_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="59"),2])
Mel60_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="60"),2])
Mel65_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="65"),2])
Mel67_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="67"),2])
Mel71_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="71"),2])
Mel72_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="72"),2])
Mel74_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="74"),2])
Mel75_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="75"),2])
Mel78_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="78"),2])
Mel79_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="79"),2])
Mel80_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="80"),2])
Mel81_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="81"),2])
Mel82_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="82"),2])
Mel84_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="84"),2])
Mel88_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="88"),2])
Mel89_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="89"),2])
Mel94_CD4<-data.frame(Class_CD4[which(substr(rownames(Class_CD4),3,4)=="94"),2])

b<-rbind(as.matrix(table(Mel53_CD4)),as.matrix(table(Mel71_CD4)),
         as.matrix(table(Mel79_CD4)),as.matrix(table(Mel80_CD4)),
         as.matrix(table(Mel81_CD4)),as.matrix(table(Mel82_CD4)),
         as.matrix(table(Mel84_CD4)),as.matrix(table(Mel89_CD4)),
         
         as.matrix(table(Mel58_CD4)),as.matrix(table(Mel60_CD4)),
         as.matrix(table(Mel72_CD4)),as.matrix(table(Mel74_CD4)),
         as.matrix(table(Mel75_CD4)),as.matrix(table(Mel88_CD4)),
         as.matrix(table(Mel94_CD4)),
         
         as.matrix(table(Mel65_CD4)),as.matrix(table(Mel67_CD4)))


treatment_CD4<-data.frame(Treatment=c("Untreated","Resistant","NoInf"),
                          Number=c(sum(rbind(as.matrix(table(Mel53_CD4)),
                                             as.matrix(table(Mel71_CD4)),
                                             as.matrix(table(Mel79_CD4)),
                                             as.matrix(table(Mel80_CD4)),
                                             as.matrix(table(Mel81_CD4)),
                                             as.matrix(table(Mel82_CD4)),
                                             as.matrix(table(Mel84_CD4)),
                                             as.matrix(table(Mel89_CD4)))),
                                   sum(rbind(as.matrix(table(Mel58_CD4)),
                                             as.matrix(table(Mel60_CD4)),
                                             as.matrix(table(Mel72_CD4)),
                                             as.matrix(table(Mel74_CD4)),
                                             as.matrix(table(Mel75_CD4)),
                                             as.matrix(table(Mel88_CD4)),
                                             as.matrix(table(Mel94_CD4)))),
                                   sum(rbind(as.matrix(table(Mel65_CD4)),
                                             as.matrix(table(Mel67_CD4))))),
                          col=c("#3838FF","#FF3838","gray30"))

Untreated=colorRampPalette(c("#3838FF","white"))(15)
Resistant=colorRampPalette(c("#FF3838","white"))(15)
NoInf=colorRampPalette(c("gray30","white"))(15)

patients_CD4<-data.frame(Patients=c("Mel53","Mel71","Mel79","Mel80",
                                    "Mel81","Mel82","Mel84","Mel89",
                                    
                                    "Mel58","Mel60","Mel72","Mel74",
                                    "Mel75","Mel88","Mel94",
                                    
                                    "Mel65","Mel67"),
                         Number=c(nrow(Mel53_CD4),nrow(Mel71_CD4),
                                  nrow(Mel79_CD4),nrow(Mel80_CD4),
                                  nrow(Mel81_CD4),nrow(Mel82_CD4),
                                  nrow(Mel84_CD4),nrow(Mel89_CD4),
                                  
                                  nrow(Mel58_CD4),nrow(Mel60_CD4),
                                  nrow(Mel72_CD4),nrow(Mel74_CD4),
                                  nrow(Mel75_CD4),nrow(Mel88_CD4),
                                  nrow(Mel94_CD4),
                                  
                                  nrow(Mel65_CD4),nrow(Mel67_CD4)),
                         col=c(Untreated[3:10],Resistant[3:9],NoInf[3:4]))

Cell_Types_CD4<-data.frame(Cells=substr(rownames(b),1,3),
                           Number=rbind(as.matrix(table(Mel53_CD4)),
                                        as.matrix(table(Mel71_CD4)),
                                        as.matrix(table(Mel79_CD4)),
                                        as.matrix(table(Mel80_CD4)),
                                        as.matrix(table(Mel81_CD4)),
                                        as.matrix(table(Mel82_CD4)),
                                        as.matrix(table(Mel84_CD4)),
                                        as.matrix(table(Mel89_CD4)),
                                        
                                        as.matrix(table(Mel58_CD4)),
                                        as.matrix(table(Mel60_CD4)),
                                        as.matrix(table(Mel72_CD4)),
                                        as.matrix(table(Mel74_CD4)),
                                        as.matrix(table(Mel75_CD4)),
                                        as.matrix(table(Mel88_CD4)),
                                        as.matrix(table(Mel94_CD4)),
                                        
                                        as.matrix(table(Mel65_CD4)),
                                        as.matrix(table(Mel67_CD4))),
                           col=c(1:31))


library(plyr)
Cell_Types_CD4$col <- revalue(Cell_Types_CD4$Cells, c("1"="royalblue1",
                                                      "2"="firebrick3",
                                                      "3"="firebrick2",
                                                      "4"="black",
                                                      "5"="white"))

Cell_Types_CD4$Cells <- revalue(Cell_Types_CD4$Cells, c("1"="Th1",
                                                        "2"="Th2",
                                                        "3"="Th17",
                                                        "4"="Th0",
                                                        "5"="Treg"))



pdf("Melanoma_CD4_Patients.pdf",7,7)

pie(Cell_Types_CD4$Number, border = NA, 
    radius = 1, labels = NA,cex=1.5,
    col=as.vector(Cell_Types_CD4$col))

x=c()
x[1]<-Cell_Types_CD4$Number[1]/2

for(i in 2:nrow(Cell_Types_CD4)){
  x[i]=sum(Cell_Types_CD4$Number[1:i-1])+
    Cell_Types_CD4$Number[i]/2
}

lab<-Cell_Types_CD4$Cells

pie.labels(x=0,y=0,(x/sum(Cell_Types_CD4$Number))*2*pi,
           labels= lab,radius=1.05,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.6)



par(new = TRUE)
pie(patients_CD4$Number, border = "white",
    radius = 1, labels = NA,col = NA)

par(new = TRUE)
pie(patients_CD4$Number, border = NA,
    radius = .73, labels = NA,col = "white")

par(new = TRUE)
pie(patients_CD4$Number, border = NA, 
    radius = .69, labels = NA,
    col = as.vector(patients_CD4$col))

x=c()
x[1]<-patients_CD4$Number[1]/2
for(i in 2:nrow(patients_CD4)){
  x[i]=sum(patients_CD4$Number[1:i-1])+
    patients_CD4$Number[i]/2
}

lab<-patients_CD4$Patients

pie.labels(x=0,y=0,(x/sum(Cell_Types_CD4$Number))*2*pi,
           labels= lab,radius=.5,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.6, col="black")



par(new = TRUE)
pie(treatment_CD4$Number, border = NA,
    radius = .4, labels = NA,
    col = as.vector(treatment_CD4$col))

x=c()
x[1]<-treatment_CD4$Number[1]/2
for(i in 2:nrow(treatment_CD4)){
  x[i]=sum(treatment_CD4$Number[1:i-1])+
    treatment_CD4$Number[i]/2
}

lab<-treatment_CD4$Treatment
pie.labels(x=0,y=0,(x/sum(Cell_Types_CD4$Number))*2*pi,
           labels= lab,radius=.15,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.8, col="black")


dev.off()




##Plot T helper Classification of CD4 For Each Patient Based on Percentage--------------


Ptreatment_CD4<-data.frame(Treatment=c("Untreated","Resistant","NoInfo"),
                           Number=c(8,7,2),
                           col=c("#3838FF","#FF3838","gray30"))

Untreated=colorRampPalette(c("#3838FF","white"))(15)
Resistant=colorRampPalette(c("#FF3838","white"))(15)
NoInf=colorRampPalette(c("gray30","white"))(15)

Ppatients_CD4<-data.frame(Patients=c("Mel53","Mel71","Mel79","Mel80",
                                     "Mel81","Mel82","Mel84","Mel89",
                                     
                                     "Mel58","Mel60","Mel72","Mel74",
                                     "Mel75",        "Mel88","Mel94",   #Mel78
                                     
                                     "Mel65","Mel67"),                  #Mel59
                          Number=rep(1,17),
                          col=c(Untreated[3:10],Resistant[3:9],
                                NoInf[3:4]))

PCell_Types_CD4<-data.frame(Cells=substr(rownames(b),1,3),
                            Number=rbind(as.matrix(table(Mel53_CD4))/nrow(Mel53_CD4),
                                         as.matrix(table(Mel71_CD4))/nrow(Mel71_CD4),
                                         as.matrix(table(Mel79_CD4))/nrow(Mel79_CD4),
                                         as.matrix(table(Mel80_CD4))/nrow(Mel80_CD4),
                                         as.matrix(table(Mel81_CD4))/nrow(Mel81_CD4),
                                         as.matrix(table(Mel82_CD4))/nrow(Mel82_CD4),
                                         as.matrix(table(Mel84_CD4))/nrow(Mel84_CD4),
                                         as.matrix(table(Mel89_CD4))/nrow(Mel89_CD4),
                                         
                                         as.matrix(table(Mel58_CD4))/nrow(Mel58_CD4),
                                         as.matrix(table(Mel60_CD4))/nrow(Mel60_CD4),
                                         as.matrix(table(Mel72_CD4))/nrow(Mel72_CD4),
                                         as.matrix(table(Mel74_CD4))/nrow(Mel74_CD4),
                                         as.matrix(table(Mel75_CD4))/nrow(Mel75_CD4),
                                         as.matrix(table(Mel88_CD4))/nrow(Mel88_CD4),
                                         as.matrix(table(Mel94_CD4))/nrow(Mel94_CD4),
                                         
                                         as.matrix(table(Mel65_CD4))/nrow(Mel65_CD4),
                                         as.matrix(table(Mel67_CD4))/nrow(Mel67_CD4)),
                            col=c(1:31))


library(plyr)
PCell_Types_CD4$col <- revalue(PCell_Types_CD4$Cells, c("1"="royalblue1",
                                                        "2"="firebrick3",
                                                        "3"="firebrick2",
                                                        "4"="black",
                                                        "5"="white"))

PCell_Types_CD4$Cells <- revalue(PCell_Types_CD4$Cells, c("1"="Th1",
                                                          "2"="Th2",
                                                          "3"="Th17",
                                                          "4"="Th0",
                                                          "5"="Treg"))



pdf("Melanoma_Patients_Percentage_CD4.pdf",8,7)

pie(PCell_Types_CD4$Number, border = NA, 
    radius = 1, labels = NA, cex=1.5,
    col = as.vector(PCell_Types_CD4$col))

par(new = TRUE)
pie(Ppatients_CD4$Number, border = "white",
    radius = 1, labels = NA, lwd=10,
    col = NA)

par(new = TRUE)
pie(Ppatients_CD4$Number, border = NA, 
    radius = .73, labels = NA,
    col = "white")

par(new = TRUE)
pie(Ppatients_CD4$Number, border = NA,
    radius = .69, labels = NA,
    col = as.vector(Ppatients_CD4$col))

x=c()
x[1]<-Ppatients_CD4$Number[1]/2

for(i in 2:nrow(Ppatients_CD4)){
  x[i]=sum(Ppatients_CD4$Number[1:i-1])+
    Ppatients_CD4$Number[i]/2
}

lab<-Ppatients_CD4$Patients
pie.labels(x=0,y=0,(x/sum(PCell_Types_CD4$Number))*2*pi,
           labels= lab,radius=.5,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.9, col="black")


par(new = TRUE)
pie(Ptreatment_CD4$Number, border = NA,
    radius = .4, labels = NA,
    col = as.vector(Ptreatment_CD4$col))

x=c()
x[1]<-Ptreatment_CD4$Number[1]/2
for(i in 2:nrow(Ptreatment_CD4)){
  x[i]=sum(Ptreatment_CD4$Number[1:i-1])+
    Ptreatment_CD4$Number[i]/2
}

lab<-Ptreatment_CD4$Treatment
pie.labels(x=0,y=0,(x/sum(PCell_Types_CD4$Number))*2*pi,
           labels= lab,radius=.13,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=1.2, col="black")


colors=c("royalblue1","firebrick3","white","white","white")
Labels=c("Th1","Th2","Th17","Th0","Treg")
legend(1, .3, Labels, cex = 1, fill = colors,
       border = NA,box.col="white")

dev.off()




##Plot Immune Cell Classification of Melanoma Data For Each Patient--------



Mel<-data.frame(Normalized_Melanoma_Labeled,
                Normalized_Melanoma_Unlabeled[,-c(1,2)])

Test_Mel <-log2(as.matrix(data.frame(t(Mel[,-c(1,2)])))+1)
Test_Mel[is.na(Test_Mel)]<- -1

prediction_Mel <- data.frame(predict(Lasso, newx=data.matrix(Test_Mel),
                                     interval ="prediction",
                                     type="response", s=0.0001))

colnames(prediction_Mel)<-c("B Cell","CD4","CD8","Mono","Neu",
                            "NK","DC","M1","M2","MQ")

Class_Mel<-which(as.matrix(prediction_Mel)==rowMaxs(as.matrix(prediction_Mel)) ,
                 arr.ind = TRUE)

prediction_Class_Mel <-predict(Lasso, newx=data.matrix(Test_Mel),
                               interval ="prediction",
                               type="class", s=0.0001)

table(substr(rownames(Class_Mel),3,4))

Mel53<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="53"),2])
Mel58<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="58"),2])
Mel59<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="59"),2])
Mel60<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="60"),2])
Mel65<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="65"),2])
Mel67<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="67"),2])
Mel71<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="71"),2])
Mel72<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="72"),2])
Mel74<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="74"),2])
Mel75<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="75"),2])
Mel78<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="78"),2])
Mel79<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="79"),2])
Mel80<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="80"),2])
Mel81<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="81"),2])
Mel82<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="82"),2])
Mel84<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="84"),2])
Mel88<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="88"),2])
Mel89<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="89"),2])
Mel94<-data.frame(Class_Mel[which(substr(rownames(Class_Mel),3,4)=="94"),2])

a<-rbind(as.matrix(table(Mel53)),as.matrix(table(Mel71)),
         as.matrix(table(Mel79)),as.matrix(table(Mel80)),
         as.matrix(table(Mel81)),as.matrix(table(Mel82)),
         as.matrix(table(Mel84)),as.matrix(table(Mel89)),
         
         as.matrix(table(Mel58)),as.matrix(table(Mel60)),
         as.matrix(table(Mel72)),as.matrix(table(Mel74)),
         as.matrix(table(Mel75)),as.matrix(table(Mel78)),
         as.matrix(table(Mel88)),as.matrix(table(Mel94)),
         
         as.matrix(table(Mel59)),as.matrix(table(Mel65)),
         as.matrix(table(Mel67)))

treatment<-data.frame(Treatment=c("Untreated","Resistant","NoInf"),
                      Number=c(sum(rbind(as.matrix(table(Mel53)),
                                         as.matrix(table(Mel71)),
                                         as.matrix(table(Mel79)),
                                         as.matrix(table(Mel80)),
                                         as.matrix(table(Mel81)),
                                         as.matrix(table(Mel82)),
                                         as.matrix(table(Mel84)),
                                         as.matrix(table(Mel89)))),
                               sum(rbind(as.matrix(table(Mel58)),
                                         as.matrix(table(Mel60)),
                                         as.matrix(table(Mel72)),
                                         as.matrix(table(Mel74)),
                                         as.matrix(table(Mel75)),
                                         as.matrix(table(Mel78)),
                                         as.matrix(table(Mel88)),
                                         as.matrix(table(Mel94)))),
                               sum(rbind(as.matrix(table(Mel59)),
                                         as.matrix(table(Mel65)),
                                         as.matrix(table(Mel67))))),
                      col=c("#3838FF","#FF3838","gray30"))


Untreated=colorRampPalette(c("#3838FF","white"))(15)
Resistant=colorRampPalette(c("#FF3838","white"))(15)
NoInf=colorRampPalette(c("gray30","white"))(15)

patients<-data.frame(Patients=c("Mel53","Mel71","Mel79","Mel80",
                                "Mel81","Mel82","Mel84","Mel89",
                                
                                "Mel58","Mel60","Mel72","Mel74",
                                "Mel75","Mel78","Mel88","Mel94",
                                
                                "Mel59","Mel65","Mel67"),
                     Number=c(nrow(Mel53),nrow(Mel71),
                              nrow(Mel79),nrow(Mel80),
                              nrow(Mel81),nrow(Mel82),
                              nrow(Mel84),nrow(Mel89),
                              
                              nrow(Mel58),nrow(Mel60),
                              nrow(Mel72),nrow(Mel74),
                              nrow(Mel75),nrow(Mel78),
                              nrow(Mel88),nrow(Mel94),
                              
                              nrow(Mel59),nrow(Mel65),
                              nrow(Mel67)),
                     col=c(Untreated[3:10],Resistant[3:10],NoInf[3:5]))


Cell_Types<-data.frame(Cells=substr(rownames(a),1,3),
                       Number=rbind(as.matrix(table(Mel53)),
                                    as.matrix(table(Mel71)),
                                    as.matrix(table(Mel79)),
                                    as.matrix(table(Mel80)),
                                    as.matrix(table(Mel81)),
                                    as.matrix(table(Mel82)),
                                    as.matrix(table(Mel84)),
                                    as.matrix(table(Mel89)),
                                    
                                    as.matrix(table(Mel58)),
                                    as.matrix(table(Mel60)),
                                    as.matrix(table(Mel72)),
                                    as.matrix(table(Mel74)),
                                    as.matrix(table(Mel75)),
                                    as.matrix(table(Mel78)),
                                    as.matrix(table(Mel88)),
                                    as.matrix(table(Mel94)),
                                    
                                    as.matrix(table(Mel59)),
                                    as.matrix(table(Mel65)),
                                    as.matrix(table(Mel67))),
                       col=c(1:99))


library(plyr)
Cell_Types$col <- revalue(Cell_Types$Cells, c("1"="royalblue1",
                                              "2"="firebrick3",
                                              "3"="firebrick2",
                                              "4"="black",
                                              "5"="white",
                                              "6"="yellow1",
                                              "7"="grey",
                                              "8"="white",
                                              "9"="white",
                                              "10"="darkorange"))

Cell_Types$Cells <- revalue(Cell_Types$Cells, c("1"="BCell",
                                                "2"="CD4",
                                                "3"="CD8",
                                                "4"="Mono",
                                                "5"="Neu",
                                                "6"="NK",
                                                "7"="DC",
                                                "8"="M1",
                                                "9"="M2",
                                                "10"="MQ"))


pdf("Melanoma_Patients.pdf",7,7)

pie(Cell_Types$Number, border = NA,
    radius = 1, labels = NA, cex=1.5,
    col = as.vector(Cell_Types$col))

x=c()
x[1]<-Cell_Types$Number[1]/2

for(i in 2:nrow(Cell_Types)){
  x[i]=sum(Cell_Types$Number[1:i-1])+
    Cell_Types$Number[i]/2
}

lab<-c(rep(" ",8),as.vector(Cell_Types$Cells[9:12])," ",
       as.vector(Cell_Types$Cells[14:19]), rep(" ",18),
       as.vector(Cell_Types$Cells[38:42]), rep(" ",7),
       as.vector(Cell_Types$Cells[50:53]), rep(" ",12),
       as.vector(Cell_Types$Cells[66]), rep(" ",12),
       as.vector(Cell_Types$Cells[79:83]),
       rep(" ",nrow(Cell_Types)-84))

pie.labels(x=0,y=0,(x/sum(Cell_Types$Number))*2*pi,
           labels= lab,radius=1.05,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.6)


par(new = TRUE)
pie(patients$Number, border = NA, 
    radius = .73, labels = NA,
    col = "white")

par(new = TRUE)
pie(patients$Number, border = NA, 
    radius = .69, labels = NA,
    col = as.vector(patients$col))

x=c()
x[1]<-patients$Number[1]/2

for(i in 2:nrow(patients)){
  x[i]=sum(patients$Number[1:i-1])+
    patients$Number[i]/2
}

lab<-patients$Patients
pie.labels(x=0,y=0,(x/sum(Cell_Types$Number))*2*pi,
           labels= lab,radius=.5,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.6, col="black")


par(new = TRUE)
pie(treatment$Number, border = NA, 
    radius = .4, labels = NA,
    col = as.vector(treatment$col))

x=c()
x[1]<-treatment$Number[1]/2
for(i in 2:nrow(treatment)){
  x[i]=sum(treatment$Number[1:i-1])+
    treatment$Number[i]/2
}

lab<-treatment$Treatment
pie.labels(x=0,y=0,(x/sum(Cell_Types$Number))*2*pi,
           labels= lab,radius=.15,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.8, col="black")

dev.off()





##Plot Immune Cell Classification of Melanoma Data For Each Patient Based on Percentage--------------



Ptreatment<-data.frame(Treatment=c("Untreated","Resistant","NoInfo"),
                       Number=c(8,8,3),
                       col=c("#3838FF","#FF3838","gray30"))

Untreated=colorRampPalette(c("#3838FF","white"))(15)
Resistant=colorRampPalette(c("#FF3838","white"))(15)
NoInf=colorRampPalette(c("gray30","white"))(15)

Ppatients<-data.frame(Patients=c("Mel53","Mel71","Mel79","Mel80",
                                 "Mel81","Mel82","Mel84","Mel89",
                                 
                                 "Mel58","Mel60","Mel72","Mel74",
                                 "Mel75","Mel78","Mel88","Mel94",
                                 
                                 "Mel59","Mel65","Mel67"),
                      Number=rep(1,19),
                      col=c(Untreated[3:10],Resistant[3:10],
                            NoInf[3:5]))

PCell_Types<-data.frame(Cells=substr(rownames(a),1,3),
                        Number=rbind(as.matrix(table(Mel53))/nrow(Mel53),
                                     as.matrix(table(Mel71))/nrow(Mel71),
                                     as.matrix(table(Mel79))/nrow(Mel79),
                                     as.matrix(table(Mel80))/nrow(Mel80),
                                     as.matrix(table(Mel81))/nrow(Mel81),
                                     as.matrix(table(Mel82))/nrow(Mel82),
                                     as.matrix(table(Mel84))/nrow(Mel84),
                                     as.matrix(table(Mel89))/nrow(Mel89),
                                     
                                     as.matrix(table(Mel58))/nrow(Mel58),
                                     as.matrix(table(Mel60))/nrow(Mel60),
                                     as.matrix(table(Mel72))/nrow(Mel72),
                                     as.matrix(table(Mel74))/nrow(Mel74),
                                     as.matrix(table(Mel75))/nrow(Mel75),
                                     as.matrix(table(Mel78))/nrow(Mel78),
                                     as.matrix(table(Mel88))/nrow(Mel88),
                                     as.matrix(table(Mel94))/nrow(Mel94),
                                     
                                     as.matrix(table(Mel59))/nrow(Mel59),
                                     as.matrix(table(Mel65))/nrow(Mel65),
                                     as.matrix(table(Mel67))/nrow(Mel67)),
                        col=c(1:99))

library(plyr)
PCell_Types$col <- revalue(PCell_Types$Cells, c("1"="royalblue1",
                                                "2"="firebrick3",
                                                "3"="firebrick2",
                                                "4"="black",
                                                "5"="white",
                                                "6"="yellow1",
                                                "7"="grey",
                                                "8"="white",
                                                "9"="white",
                                                "10"="darkorange"))

PCell_Types$Cells <- revalue(PCell_Types$Cells, c("1"="BCell",
                                                  "2"="CD4",
                                                  "3"="CD8",
                                                  "4"="Mono",
                                                  "5"="Neu",
                                                  "6"="NK",
                                                  "7"="DC",
                                                  "8"="M1",
                                                  "9"="M2",
                                                  "10"="MQ"))


pdf("Melanoma_Patients_Percentage.pdf",8,7)

pie(PCell_Types$Number, border = NA, 
    radius = 1, labels = NA, cex=1.5,
    col = as.vector(PCell_Types$col))

par(new = TRUE)
pie(Ppatients$Number, border = NA, 
    radius = 1, labels = NA, lwd=10,
    col = NA)

par(new = TRUE)
pie(Ppatients$Number, border = NA,
    radius = .73, labels = NA,
    col = "white")

par(new = TRUE)
pie(Ppatients$Number, border = NA,
    radius = .69, labels = NA,
    col = as.vector(Ppatients$col))

x=c()
x[1]<-Ppatients$Number[1]/2

for(i in 2:nrow(Ppatients)){
  x[i]=sum(Ppatients$Number[1:i-1])+
    Ppatients$Number[i]/2
}

lab<-Ppatients$Patients
pie.labels(x=0,y=0,(x/sum(PCell_Types$Number))*2*pi,
           labels= lab,radius=.5,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=.8, col="black")


par(new = TRUE)
pie(Ptreatment$Number, border = NA, 
    radius = .4, labels = NA,
    col = as.vector(Ptreatment$col))

x=c()
x[1]<-Ptreatment$Number[1]/2

for(i in 2:nrow(Ptreatment)){
  x[i]=sum(Ptreatment$Number[1:i-1])+
    Ptreatment$Number[i]/2
}

lab<-Ptreatment$Treatment
pie.labels(x=0,y=0,(x/sum(PCell_Types$Number))*2*pi,
           labels= lab,radius=.12,bg=NA,border=NA,
           minangle=NA,boxed=FALSE,cex=1.2, col="black")

colors=c("royalblue1","firebrick3","firebrick2","black",
         "yellow1","grey","darkorange","white","white",
         "white")

Labels=c("BCell","CD4","CD8","Mono","NK","DC",
         expression(paste("M",phi)),"M1","M2",
         "Neu")

legend(1, .5, Labels, cex = 1, fill = colors,
       border = NA,box.col="white")


dev.off()

##ICIs Treatment Correlation by Immune Cells Composition---------



Ptreatment<-data.frame(Treatment=c("Untreated","Resistant","NoInfo"),
                       Number=c(8,8,3),
                       col=c("#3838FF","#FF3838","gray30"))

Untreated=colorRampPalette(c("#3838FF","white"))(15)
Resistant=colorRampPalette(c("#FF3838","white"))(15)
NoInf=colorRampPalette(c("gray30","white"))(15)

Ppatients<-data.frame(Patients=c("Mel53","Mel71","Mel79","Mel80",
                                 "Mel81","Mel82","Mel84","Mel89",
                                 
                                 "Mel58","Mel60","Mel72","Mel74",
                                 "Mel75","Mel78","Mel88","Mel94",
                                 
                                 "Mel59","Mel65","Mel67"),
                      Number=rep(1,19),
                      col=c(Untreated[3:10],Resistant[3:10],
                            NoInf[3:5]))

pMel53<-data.frame(Cells=rownames(as.matrix(table(Mel53))/nrow(Mel53)), 
                   Mel53=as.matrix(table(Mel53))/nrow(Mel53))
pMel71<-data.frame(Cells=rownames(as.matrix(table(Mel71))/nrow(Mel71)), 
                   Mel71=as.matrix(table(Mel71))/nrow(Mel71))
pMel79<-data.frame(Cells=rownames(as.matrix(table(Mel79))/nrow(Mel79)), 
                   Mel79=as.matrix(table(Mel79))/nrow(Mel79))
pMel80<-data.frame(Cells=rownames(as.matrix(table(Mel80))/nrow(Mel80)), 
                   Mel80=as.matrix(table(Mel80))/nrow(Mel80))
pMel81<-data.frame(Cells=rownames(as.matrix(table(Mel81))/nrow(Mel81)), 
                   Mel81=as.matrix(table(Mel81))/nrow(Mel81))
pMel82<-data.frame(Cells=rownames(as.matrix(table(Mel82))/nrow(Mel82)), 
                   Mel82=as.matrix(table(Mel82))/nrow(Mel82))
pMel84<-data.frame(Cells=rownames(as.matrix(table(Mel84))/nrow(Mel84)), 
                   Mel84=as.matrix(table(Mel84))/nrow(Mel84))
pMel89<-data.frame(Cells=rownames(as.matrix(table(Mel89))/nrow(Mel89)), 
                   Mel89=as.matrix(table(Mel89))/nrow(Mel89))

pMel58<-data.frame(Cells=rownames(as.matrix(table(Mel58))/nrow(Mel58)), 
                   Mel58=as.matrix(table(Mel58))/nrow(Mel58))
pMel60<-data.frame(Cells=rownames(as.matrix(table(Mel60))/nrow(Mel60)), 
                   Mel60=as.matrix(table(Mel60))/nrow(Mel60))
pMel72<-data.frame(Cells=rownames(as.matrix(table(Mel72))/nrow(Mel72)), 
                   Mel72=as.matrix(table(Mel72))/nrow(Mel72))
pMel74<-data.frame(Cells=rownames(as.matrix(table(Mel74))/nrow(Mel74)), 
                   Mel74=as.matrix(table(Mel74))/nrow(Mel74))
pMel75<-data.frame(Cells=rownames(as.matrix(table(Mel75))/nrow(Mel75)), 
                   Mel75=as.matrix(table(Mel75))/nrow(Mel75))
pMel78<-data.frame(Cells=rownames(as.matrix(table(Mel78))/nrow(Mel78)), 
                   Mel78=as.matrix(table(Mel78))/nrow(Mel78))
pMel88<-data.frame(Cells=rownames(as.matrix(table(Mel88))/nrow(Mel88)), 
                   Mel88=as.matrix(table(Mel88))/nrow(Mel88))
pMel94<-data.frame(Cells=rownames(as.matrix(table(Mel94))/nrow(Mel94)), 
                   Mel94=as.matrix(table(Mel94))/nrow(Mel94))

pMel59<-data.frame(Cells=rownames(as.matrix(table(Mel59))/nrow(Mel59)), 
                   Mel59=as.matrix(table(Mel59))/nrow(Mel59))
pMel65<-data.frame(Cells=rownames(as.matrix(table(Mel65))/nrow(Mel65)), 
                   Mel65=as.matrix(table(Mel65))/nrow(Mel65))
pMel67<-data.frame(Cells=rownames(as.matrix(table(Mel67))/nrow(Mel67)), 
                   Mel67=as.matrix(table(Mel67))/nrow(Mel67))


Patient_Table<-Reduce(function(x,y) merge(x = x, y = y, 
                                          by = "Cells",all=TRUE), 
                      list(pMel53,pMel71,pMel79,pMel80,
                           pMel81,pMel82,pMel84,pMel89,
                           
                           pMel58,pMel60,pMel72,pMel74,
                           pMel75,pMel78,pMel88,pMel94,
                           
                           pMel59,pMel65,pMel67))



Patient_Table$Cells<-revalue(Patient_Table$Cells,c("1"="BCell",
                                                   "2"="CD4",
                                                   "3"="CD8",
                                                   "4"="Mono",
                                                   "5"="Neu",
                                                   "6"="NK",
                                                   "7"="DC",
                                                   "8"="M1",
                                                   "9"="M2",
                                                   "10"="MQ"))

Patient_Table[is.na(Patient_Table)]<-0


#PCA

pc=prcomp(as.matrix(t(Patient_Table[,-1])))
pc.data <- as.matrix(t(Patient_Table[,-1]))%*%pc$rotation

par(mfrow=c(1,2))

plot(pc.data[,c(1,2)], type = "n")
text(pc.data[,c(1,2)], labels =colnames(Patient_Table[,-1]) ,
     col=c(rep("#3838FF",8),
           rep("#FF3838",8),
           rep("gray30",3)))

plot(pc.data[,c(2,3)], type = "n")
text(pc.data[,c(2,3)], labels =colnames(Patient_Table[,-1]) ,
     col=c(rep("#3838FF",8),
           rep("#FF3838",8),
           rep("gray30",3)))

par(mfrow=c(1,1))


#Clustring

dist1=dist(t(as.matrix(Patient_Table[,-1])))
hclust1=hclust(dist1)

library(rafalib)

myplclust(hclust1, labels = hclust1$labels, 
          lab.col = c(rep("#3838FF",8),
                      rep("#FF3838",8),
                      rep("gray30",3)), 
          hang = 0.1)



##ICIs Treatment Correlation by Immune Cells and T helper Cells Composition---------

real_Percentage_CD4<-rbind(as.matrix(table(Mel53_CD4))/nrow(Mel53),
                           as.matrix(table(Mel71_CD4))/nrow(Mel71),
                           as.matrix(table(Mel79_CD4))/nrow(Mel79),
                           as.matrix(table(Mel80_CD4))/nrow(Mel80),
                           as.matrix(table(Mel81_CD4))/nrow(Mel81),
                           as.matrix(table(Mel82_CD4))/nrow(Mel82),
                           as.matrix(table(Mel84_CD4))/nrow(Mel84),
                           as.matrix(table(Mel89_CD4))/nrow(Mel89),
                           
                           as.matrix(table(Mel58_CD4))/nrow(Mel58),
                           as.matrix(table(Mel60_CD4))/nrow(Mel60),
                           as.matrix(table(Mel72_CD4))/nrow(Mel72),
                           as.matrix(table(Mel74_CD4))/nrow(Mel74),
                           as.matrix(table(Mel75_CD4))/nrow(Mel75),
                           as.matrix(table(Mel88_CD4))/nrow(Mel88),
                           as.matrix(table(Mel94_CD4))/nrow(Mel94),
                           
                           as.matrix(table(Mel65_CD4))/nrow(Mel65),
                           as.matrix(table(Mel67_CD4))/nrow(Mel67))

Patient_Table_withTh<-rbind(Patient_Table,c(0,real_Percentage_CD4[c(1,3,5,7,9,10,11,13,15,
                                                                    17,19,21,23)],0,
                                            real_Percentage_CD4[c(25,26)],0,
                                            real_Percentage_CD4[c(28,30)]),
                            c(0,real_Percentage_CD4[c(2,4,6,8)],0,0,
                              real_Percentage_CD4[c(12,14,16,18,20,22,24)],0,0,
                              real_Percentage_CD4[27],0,
                              real_Percentage_CD4[c(29,31)]))

Zeros=data.frame(zeros(6,20))
colnames(Zeros)<-colnames(Patient_Table_withTh)
Table_All<-rbind(Patient_Table_withTh,Zeros)



#PCA

pc_withTh=prcomp(as.matrix(t(Table_All[,-1][-3,])))
pc_withTh.data <- as.matrix(t(Table_All[,-1][-3,]))%*%pc_withTh$rotation

par(mfrow=c(1,2))

plot(pc_withTh.data[,c(1,2)], type = "n")
text(pc_withTh.data[,c(1,2)], labels =colnames(Table_All[,-1]) ,
     col=c(rep("#3838FF",8),
           rep("#FF3838",8),
           rep("gray30",3)))
plot(pc_withTh.data[,c(2,3)], type = "n")
text(pc_withTh.data[,c(2,3)], labels =colnames(Table_All[,-1]) ,
     col=c(rep("#3838FF",8),
           rep("#FF3838",8),
           rep("gray30",3)))

par(mfrow=c(1,1))




#Scree plot

EigenVals <- pc_withTh$sdev^2
TotVar <- sum(EigenVals)
PrPC <- EigenVals/TotVar
#0.6713400590 0.2105285065 0.0832556624 0.0240638567 0.0083360867 
#0.0015505426 0.0007265323 0.0001987539

barplot(PrPC*100, names.arg = c(1:14), ylab = "Variance (%)", 
        xlab = "Principal Component", ylim = c(0,70))

TestNeg<-Table_All[,-1][-3,]
Neg<-as.matrix(TestNeg[sample(nrow(TestNeg)),])
Neg<-as.matrix(Neg[,sample(ncol(TestNeg))])

for (i in 1:2000) {
  Neg1<-as.matrix(TestNeg[sample(nrow(TestNeg)),])
  Neg1<-as.matrix(Neg1[,sample(ncol(TestNeg))])  
  Neg<-data.frame(Neg,Neg1)
}

NegPCA=prcomp(as.matrix(t(Neg)))
NegEigenVals <- NegPCA$sdev^2
NegTotVar <- sum(NegEigenVals)
NegPrPC <- NegEigenVals/NegTotVar
#0.0848862440 0.0833103930 0.0814907732 0.0803736057 0.0798195401 0.0792675918
#0.0784763528 0.0772660695 0.0767707571
#0.0736236839 0.0723748547 0.0704887107 0.0617471701 0.0001042534

barplot(NegPrPC*100, names.arg = c(1:14), ylab = "Variance (%)",
        xlab = "Principal Component", ylim = c(0,70))

lines(c(0,20), NegPrPC[1:2]*100, 
      col="red", lwd = 2, lty = 2)


#Clustring

dist1_withTh=dist(t(as.matrix(Patient_Table_withTh[,-1][-3,])))
hclust1_withTh=hclust(dist1_withTh)

library(rafalib)
myplclust(hclust1_withTh, labels = hclust1_withTh$labels, 
          lab.col = c(rep("#3838FF",8),
                      rep("#FF3838",8),
                      rep("gray30",3)), 
          hang = 0.1)



#Number of clusters (NbClust)

library(factoextra)
fviz_nbclust(t(as.matrix(Patient_Table_withTh[,-1][-3,])),
             FUNcluster = hcut, method = c("gap_stat"))

library(NbClust)
NbClust(t(as.matrix(Patient_Table_withTh[,-1][-3,])), 
        diss = NULL, distance = "euclidean", 
        min.nc = 2, max.nc =10, method = "complete",
        index =  "all", alphaBeale = 0.1) 


#pdf

pdf("PCA&Dendrogram_Melanoma.pdf",10,5)

par(mfrow=c(1,2))
plot(pc_withTh.data[,c(1,2)], type = "n")
text(pc_withTh.data[,c(1,2)], 
     labels =colnames(Table_All[,-1]) ,
     col=c(rep("#3838FF",8),rep("#FF3838",8),
           rep("gray30",3)))

colors=c("#3838FF","#FF3838","gray30")
Labels=c("Untreated","Resistant","NoInfo")
legend("topright", Labels, cex = .7, fill = colors,
       border = NA)

myplclust(hclust1_withTh, labels = hclust1_withTh$labels, 
          lab.col = c(rep("#3838FF",8),
                      rep("#FF3838",8),
                      rep("gray30",3)), 
          hang = 0.1, main=NA)
dev.off()
