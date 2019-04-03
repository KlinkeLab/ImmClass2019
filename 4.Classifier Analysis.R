##############################################Info#########################################
############# Title:          Classifiers Analysis
############# Description:    This Script Generates Heatmaps and Similarity plots 
#############                 Which are Presented in the Manuscript 
############# Author:         Arezo Torang
############# Date:           1-Apr-2019
###########################################################################################-



rm(list = ls())

####Importing Data-----------------------------------------------------------



Coef <- read.table("Coef.txt", sep = "\t", quote = "")

Coef_Th <- read.table("Coef_Th.txt", sep = "\t", quote = "")

Selected_Genes <- read.table("Selected_Genes.txt", sep = "\t", quote = "")

Selected_Genes_Th <- read.table("Selected_Genes_Th.txt", sep = "\t", quote = "")



####Immune Cell Classifier---------------------------

####Heatmaps of Classifier Coefficients and Expression matrix:
##Here we creat a heatmap for coefficients of the immune cell classifier
##and one for expression matrix of corresponding selected genes

##Coefficient's Dendrograms Generation



Rowv <- rowMeans(as.matrix(Coef[-1,]), na.rm = TRUE)
distr <- dist(as.matrix(Coef[-1,]))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
ddr <- reorder(ddr, Rowv)
rowInd <- order.dendrogram(ddr)
plot(ddr, horiz = TRUE, axes = FALSE, 
     yaxs = "i", leaflab = "none")



Colv <- colMeans(as.matrix(Coef[-1,]), na.rm = TRUE)
distc <- dist(as.matrix(t(Coef[-1,])))
hcc <- hclust(distc)
ddc <- as.dendrogram(hcc)
ddc <- reorder(ddc, Colv)
colInd <- order.dendrogram(ddc)
plot(ddc, axes = FALSE, xaxs = "i",leaflab = "none")




##Expression matrix's preparation



Main<-Selected_Genes[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Selected_Genes<-data.frame(Selected_Genes[,1:2],Main)
a<-Selected_Genes
a<-log2(a[,-c(1,2)]+1)
a[is.na(a)]<--1
row.names(a)<-Selected_Genes[,2]
a<-a[,c(46:53,54:57,1:14,83:93,94:108,
        58:61,78:82,15:29,30:45,62:77)] #Baced on Coef Dendrogram



##Creating Heatmap pdf for immune cell classifier



pdf("Heatmap.pdf", width = 17, height = 11)
nf <- layout(matrix(c(1,4,7,10,2,5,8,11,3,6,9,12),3,4,byrow = TRUE),
             widths = c(2,7,2,7),heights = c(.5,.3,10), respect = TRUE)



##Coeficient Heatmap----

#Color Key
par(mar = c(.5, 2, .5, 1))
image(as.matrix(seq(-.01, .01, length.out = 21)), axes = FALSE,
      col=colorRampPalette(c("blue","white", "red"))(n = 20))

mtext(c("-0.01","0","0.01+"), side = 3, 
      las = 1, line = 0.5, cex = 1,
      at = c(0.02, .525, .98))

par(mar = c(.5, 0.5, 0, 1))
image(as.matrix(x=1:2),col="white", axes = FALSE)

#Dendrogram
par(mar = c(.5, 0.5, 0, 1))
plot(ddr, horiz = TRUE, axes = FALSE,
     yaxs = "i", leaflab = "none")

mtext("Coefficients", side = 3, las = 1,
      line = .2, cex = 1.5,at =.25)

par(mar = c(.5, 0.5, .5, 3.5))
plot(ddc, axes = FALSE, xaxs = "i",leaflab = "none")

#Heatmap Body
par(mar = c(.5, 0.5, 0, 3.5))
Color=c("deeppink","yellow","green","orange",
        "deepskyblue","Gold","gray","violet",
        "cyan","lightpink1")

image(as.matrix(x=1:10),col=Color, axes = FALSE)

par(mar = c(.5, 0.5, 0, 3.5))
image(as.matrix(t(Coef[rowInd+1,colInd])), 
      breaks= c(-1,seq(-.01, .01, length.out = 19),1),
      col=colorRampPalette(c("blue","white", "red"))(n = 20), axes = FALSE)

ylab <- rownames(Coef)[-1]
xlab <- colnames(Coef)[colInd]
xlab[7]<-expression(paste("M", phi))

mtext("Genes", side = 4, adj= 1, line = 1.5, cex = 2, at = .5)
mtext(xlab, side = 3, las = 1, line = 0.5, cex = 1.2,
      at = seq(0, 1, 1/(length(xlab)-1)))




##Expression Heatmap------

#Color Key
par(mar = c(.5, 2, .5, 1))
image(as.matrix(seq(-.01, .01, length.out = 12)), axes = FALSE,
      col=c("lightblue",colorRampPalette(c("white","red"))(n = 11)))

mtext(c("NA","0","6+"), side = 3, 
      las = 1, line = 0.5, cex = 1,
      at = c(-.07, .09, .98))

par(mar = c(.5, 0.5, 0, 1))
image(as.matrix(x=1:2),col="white", axes = FALSE)

#Dendrograms Space
par(mar = c(0.5, 0.5, 0, 1))
image(as.matrix(1:10), axes = FALSE, col="white")

mtext("log2(Expression)", side = 3, 
      las = 1, line = 0.2, cex = 1.5,at = .5)

par(mar = c(0.5, 0.5, .5, 3.5))
image(as.matrix(1:10), axes = FALSE, col="white")

#Heatmap Body
par(mar = c(.5, 0.5, 0, 3.5))
Color=c("deeppink","yellow","green","orange",
        "deepskyblue","Gold","gray","violet",
        "cyan","lightpink1")

library(matlab)

x=c(ones(8, ncol = 1),2*ones(4, ncol = 1), 3*ones(14,ncol = 1),
    4*ones(11,ncol = 1),5*ones(15,ncol = 1),6*ones(4,ncol = 1),
    7*ones(5,ncol = 1),8*ones(15,ncol = 1),9*ones(16,ncol = 1),
    10*ones(16,ncol = 1))

image(as.matrix(x),col=Color, axes = FALSE,
      oldstyle=TRUE, useRaster=TRUE)

#Same order
par(mar = c(.5, 0.5, 0, 3.5))
image(as.matrix(t(a[rowInd,])), breaks= c(seq(-7, 6, length.out = 21),100),
      col=colorRampPalette(c("blue","white", "red"))(n = 21), axes = FALSE)

mtext("Genes", side = 4, adj = 1, 
      line = 1.5, cex = 2, at = .5)

mtext(xlab, side = 3, las = 1, line = 0.5, cex = 1.2, 
      at = c(4/length(x),10/length(x),19/length(x),
             32/length(x),44/length(x),54/length(x),
             58.5/length(x),69/length(x),84/length(x),
             100/length(x)))

dev.off() #End of Heatmaps



##Similarity Matrix for Selected Genes in Immune Cells-------------



Main<-Selected_Genes[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Selected_Genes<-data.frame(Selected_Genes[,1:2],Main)
a<-Selected_Genes
a<-log2(a[,-c(1,2)]+1)
a[is.na(a)]<--1
colnames(a)<-colnames(Selected_Genes)[3:ncol(Selected_Genes)]
dist1=dist(t(a))


##Similarity Matrix pdf

pdf("Similarity.pdf", width = 8, height = 8.5)
nf <- layout(matrix(c(1,6,2,3,4,5),3,2,byrow = FALSE), 
             widths = c(6,.5), heights = c(.5,6,.5),
             respect = TRUE)

x=c(ones(14,ncol = 1),2*ones(15,ncol = 1),
    3*ones(16,ncol = 1),4*ones(8,ncol = 1),
    5*ones(4,ncol = 1),6*ones(4,ncol = 1),
    7*ones(16,ncol = 1),8*ones(5,ncol = 1),
    9*ones(11,ncol = 1),10*ones(15,ncol = 1))

xlab <- c("BCell","CD4","CD8","DC","M1","M2","Mono",
          expression(paste("M", phi)),"Neu","NK")

Center=rep(0,10)
for(i in 1:10){
  Center[i]=sum(as.vector(table(x)[1:i-1]))+
    as.vector(table(x)[i]/2)
}


par(mar = c(.5, .5, 2, .5))
colramp=colorRampPalette(c(2,"white",3))(21)

image(as.matrix(seq(-.01, .01, length.out = 22)), 
      axes = FALSE, col=colramp)

mtext(c("Low","High"), side = 3, las = 1, line = 0.5,
      cex = 1, at = c(0.02, .98))

par(mar = c(.5, .5, .5, .5))
Color=c("green","violet","cyan", "deeppink",
        "yellow","Gold","lightpink1","gray",
        "orange","deepskyblue")

image(as.matrix(x),col=Color, axes = FALSE,
      oldstyle=TRUE, useRaster=TRUE)

par(mar = c(0, 0, 0, .5))
image(as.matrix(x=1:2),col="white", axes = FALSE)

par(mar = c(0, 0, 0, .5))
Color=c("green","violet","cyan", "deeppink",
        "yellow","Gold","lightpink1","gray",
        "orange","deepskyblue")

image(as.matrix(t(x)),col=Color, axes = FALSE,
      oldstyle=TRUE, useRaster=TRUE)

par(mar = c(.5, .5, .5, .5))
image(as.matrix(x=1:2),col="white", axes = FALSE)

par(mar = c(0, 0.5, 0, .5))
colramp=colorRampPalette(c(3,"white",2))(21)

image(as.matrix(dist1),col=colramp, axes = FALSE)

mtext(xlab, side = 1, las = 1, line = 2.1, 
      cex = 1.1, at = Center/length(x))

mtext(xlab, side = 4, las = 1, line = 1.2, 
      cex = 1.2, at = Center/length(x))


dev.off() #End of Plot





####T Helper Cell Classifier---------------------------

####Heatmaps of Classifier Coefficients and Expression matrix:
##Here we creat a heatmap for coefficients of the T helper classifier
##and one for expression matrix of corresponding selected genes

##Coefficient's Dendrograms Generation



Rowv <- rowMeans(as.matrix(Coef_Th[-1,]), na.rm = TRUE)
distr <- dist(as.matrix(Coef_Th[-1,]))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
ddr <- reorder(ddr, Rowv)
rowInd <- order.dendrogram(ddr)
plot(ddr, horiz = TRUE, axes = FALSE, 
     yaxs = "i", leaflab = "none")

Colv <- colMeans(as.matrix(Coef_Th[-1,]), na.rm = TRUE)
distc <- dist(as.matrix(t(Coef_Th[-1,])))
hcc <- hclust(distc)
ddc <- as.dendrogram(hcc)
ddc <- reorder(ddc, Colv)
colInd <- order.dendrogram(ddc)
plot(ddc, axes = FALSE, xaxs = "i",
     leaflab = "none")



##Expression matrix's preparation


Main<-Selected_Genes_Th[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Selected_Genes_Th<-data.frame(Selected_Genes_Th[,1:2],Main)
a<-Selected_Genes_Th
a<-log2(a[,-c(1,2)]+1)
a[is.na(a)]<--1
row.names(a)<-Selected_Genes_Th[,2]
a<-a[,c(19:23,10:18,27:31,1:9,24:26)] #Baced of Coef Dendrogram



##Creating Heatmap pdf for T Helper cell classifier



pdf("Heatmap_Th.pdf", width = 17, height = 6)
nf <- layout(matrix(c(1,4,7,10,2,5,8,11,3,6,9,12),3,4,byrow = TRUE),
             widths = c(2,7,2,7),heights = c(.5,.3,5), respect = TRUE)

##Coeficient Heatmaps--------------

#Color Key
par(mar = c(.5, 2, .5, 1))
image(as.matrix(seq(-.05, .05, length.out = 21)), axes = FALSE,
      col=colorRampPalette(c("blue","white", "red"))(n = 20))

mtext(c("-0.05","0","0.05+"), side = 3, 
      las = 1, line = 0.5, cex = 1,
      at = c(0.02, .525, .98))

par(mar = c(.5, 0.5, 0, 1))
image(as.matrix(x=1:2),col="white", axes = FALSE)

#Dendrogram
par(mar = c(.5, 0.5, 0, 1))
plot(ddr, horiz = TRUE, axes = FALSE, 
     yaxs = "i", leaflab = "none")

mtext("Coefficients", side = 3, las = 1,
      line = .2, cex = 1.5,at =.2)

par(mar = c(.5, 0.5, .5, 3.5))
plot(ddc, axes = FALSE, xaxs = "i",leaflab = "none")

#Heatmap Body
par(mar = c(.5, 0.5, 0, 3.5))
Color=c("deeppink","yellow","green","violet","cyan")

image(as.matrix(x=1:5),col=Color, axes = FALSE)

par(mar = c(.5, 0.5, 0, 3.5))
image(as.matrix(t(Coef_Th[rowInd+1,colInd])), 
      breaks= c(-1,seq(-.05, .05, length.out = 19),1),
      col=colorRampPalette(c("blue","white", "red"))(n = 20),
      axes = FALSE)

ylab <- rownames(Coef_Th)[-1]
xlab <- colnames(Coef_Th)[colInd]

mtext("Genes", side = 4, adj = 1, 
      line = 1.5, cex = 2, at = .55)

mtext(xlab, side = 3, las = 1, line = 0.5, cex = 1.2,
      at = seq(0, 1, 1/(length(xlab)-1)))


##Expression Heatmap--------------

#Color Key
par(mar = c(.5, 2, .5, 1))
image(as.matrix(seq(-.01, .01, length.out = 12)), axes = FALSE,
      col=c("lightblue",colorRampPalette(c("white","red"))(n = 11)))

mtext(c("NA","0","6+"), side = 3, las = 1, line = 0.5,
      cex = 1, at = c(-.07, .09, .98))

par(mar = c(.5, 0.5, 0, 1))
image(as.matrix(x=1:2),col="white", axes = FALSE)

#Dendrograms Space
par(mar = c(0.5, 0.5, 0, 1))
image(as.matrix(1:10), axes = FALSE, col="white")

mtext("log2(Expression)", side = 3, las = 1,
      line = 0.2, cex = 1.5,at = .5)

par(mar = c(0.5, 0.5, .5, 3.5))
image(as.matrix(1:10), axes = FALSE, col="white")

#Heatmap Body
par(mar = c(.5, 0.5, 0, 3.5))
Color=c("deeppink","yellow","green",
        "violet","cyan")

x=c( ones(5,ncol = 1), 2*ones(9,ncol = 1),
    3*ones(5,ncol = 1),4*ones(9,ncol = 1),
    5*ones(3,ncol = 1))

image(as.matrix(x),col=Color, axes = FALSE,
      oldstyle=TRUE, useRaster=TRUE)

par(mar = c(.5, 0.5, 0, 3.5))
image(as.matrix(t(a[rowInd,])), 
      breaks= c(seq(-7, 6, length.out = 21),100),
      col=colorRampPalette(c("blue","white", "red"))(n = 21), 
      axes = FALSE)

mtext("Genes", side = 4, adj = 1, line = 1.5,
      cex = 2, at = .55)

mtext(xlab, side = 3, las = 1, line = 0.5, cex = 1.2, 
      at = c(2/length(x),9.5/length(x),16.5/length(x),
             23.5/length(x),30/length(x)))

dev.off() #End of Heatmap



##Similarity Heatmaps for Selected Genes in T Helper Cells-----------



Main<-Selected_Genes_Th[,-c(1,2)]
Main<-Main[ , order(colnames(Main))]
Selected_Genes_Th<-data.frame(Selected_Genes_Th[,1:2],Main)
b<-Selected_Genes_Th
b<-log2(b[,-c(1,2)]+1)
b[is.na(b)]<--1
dist2=dist(t(b))


##Similarity Matrix pdf

pdf("Similarity_Th.pdf", width = 8, height = 8.5)
nf <- layout(matrix(c(1,6,2,3,4,5),3,2,byrow = FALSE),
             widths = c(6,.5),heights = c(.5,6,.5),
             respect = TRUE)

x=c(ones(9,ncol = 1),2*ones(9,ncol = 1),
    3*ones(5,ncol = 1),4*ones(3,ncol = 1),
    5*ones(5,ncol = 1))

xlab <- c("Treg","Th0","Th1","Th17","Th2")

Center=rep(0,5)
for(i in 1:5){
  Center[i]=sum(as.vector(table(x)[1:i-1]))+
    as.vector(table(x)[i]/2)
}

par(mar = c(.5, .5, 2, .5))
colramp=colorRampPalette(c(2,"white",3))(21)

image(as.matrix(seq(-.01, .01, length.out = 22)), 
      axes = FALSE, col=colramp)

mtext( c("Low", "High"), side = 3, las = 1, 
      line = 0.5, cex = 1,at = c(0.02, .98))

par(mar = c(.5, .5, .5, .5))
Color=c("violet","yellow","deeppink",
        "cyan","green")

image(as.matrix(x),col=Color, axes = FALSE,
      oldstyle=TRUE, useRaster=TRUE)

par(mar = c(0, 0, 0, .5))
image(as.matrix(x=1:2),col="white", axes = FALSE)

par(mar = c(0, 0, 0, .5))
image(as.matrix(t(x)),col=Color, axes = FALSE,
      oldstyle=TRUE, useRaster=TRUE)

par(mar = c(.5, .5, .5, .5))
image(as.matrix(x=1:2),col="white", axes = FALSE)

par(mar = c(0, 0.5, 0, .5))
colramp=colorRampPalette(c(3,"white",2))(21)

image(as.matrix(dist2),col=colramp, axes = FALSE)

mtext(xlab, side = 1, las = 1, line = 2.1, 
      cex = 1.1, at = Center/length(x))

mtext(xlab, side = 4, las = 1, line = 1.2, 
      cex = 1.2, at = Center/length(x))


dev.off() #End of plot


