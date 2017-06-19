# Spiec EASI
#See https://github.com/zdk123/SpiecEasi for more info


#Install all dependencies and packages to make it works 
# Can take up to 30 min if you have nothing installed...
#library (devtools)
#library(BIOM.utils)# On my machine "biom" is not properly install with this package, hence I had to install the 2 following ones
#library (RJSONIO)
#library (biom) #Downloaded from old CRAN repository (https://cran.r-project.org/src/contrib/Archive/biom/), not maintained but it works and it is needed. 

#Run the lines below only once
#source ("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#install_github("zdk123/SpiecEasi")

#library (SpiecEasi)
#library (phyloseq)
library(igraph)
library(Matrix)


#Read the data. NB non-normalized data!!
setwd ("/home/fonvielle/Documents/Mars/FT") # VM's path
list.files()

ct.I.minus.double <- t(ct.I.minus.double)

#The spiec-EASI computation, nothing to do actually, everything is done by just this line...
#See link,the git-hub webpage, at the beginning of the script for more info. 
Speak.E <- spiec.easi(ct.I.minus.double.1,method="mb", lambda.min.ratio=0.01,nlambda=30,
                      icov.select.params =list(rep.num=100,ncores=3))
Speak.E # network sparsitiy is given
Speak.E$sparsity 

#Weight of edges
Weigh.edges <- (symBeta(getOptBeta(Speak.E),mode="maxabs"))


#The network plot
For.Plot.A.1 <- adj2igraph(Speak.E$refit, vertex.attr=list("name", value=colnames(ct.I.minus.double.1))) # Convert adjency matrix to an igraph object
vsize.A.1 <- as.numeric(rowMeans(clr(ct.I.minus.double.1))+6) #Just a size estimation for proper ploting, nb there is the log centered transformation in it for proper size if nodes.  
am.coord.A.1 <- layout.fruchterman.reingold(For.Plot.A.1) # Function to set a force-directed drawing
For.Plot.A.1 <- set.vertex.attribute(For.Plot.A.1, "name", value=colnames(ct.I.minus.double.1))

par(mar=c(1,1,1,1))
plot(For.Plot.A.1,layout=am.coord.A.1,vertex.size=vsize.A.1,vertex.color=Color.HC.1, vertex.label="" )
#dev.off()

#write.graph(For.Plot.A.1,"total.graphml",format="graphml")

#Degree Statistics / Degree distribution
Degree.Stat <- degree_distribution(For.Plot)
plot(Degree.Stat)

Cluster <- read.csv("G_2 default node.csv")


Time.1 <- which(env.data.MARS$Time==3)
Time.2 <- which(env.data.MARS$Time==4)
Time.3 <- which(env.data.MARS$Time==5)
Time.4 <- which(env.data.MARS$Time==7)
Time.5 <- which(env.data.MARS$Time==8)

#Van Krev Diagram of Clusters (all)
Color.clust <- (Cluster[,1])
levels(Color.clust)
colfunc <- colorRampPalette(c("red","yellow","cyan", "green3","black"))
colfunc(5)
levels(Color.clust)<- colfunc(5)
Color.clust <- as.character(Color.clust)

element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))

par(mar=c(5,5,5,5))
#setEPS()
#postscript("Cluster-VK.eps")
plot.vankrevelen(element.matrix=element.matrix,color.criterion=Color.clust,plotting.order="juggle",
                 main.hl=paste("Cluster in VK space (",nrow(element.matrix),")",sep=""),cex=1)
dev.off()


#Each Cluster in an individual VK diagram 
Cluster.1 <- which (Cluster[,1]==1)

Cluster.1 <- match(mass.list.1[Cluster.1],paf[,1])
element.matrix<-paf[Cluster.1,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))



setEPS()
postscript("Cluster1-VK.eps")
plot.vankrevelen(element.matrix=element.matrix,plotting.order="juggle",
                 main.hl=paste("Cluster 1 (Red) in VK space (",nrow(element.matrix),")",sep=""),cex=1)
dev.off()
AIC1<- (1+element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-0.5*element.matrix$elH)/(element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-element.matrix$elN-element.matrix$elP)
Fck <- which(is.infinite(AIC1)==T)
mean(AIC1[-Fck])
sd (AIC1[-Fck])
DBE1 <- element.matrix$elC - (element.matrix$elH/2)+(element.matrix$elN/2)+1
mean(DBE1)
sd(DBE1)
C.N.1 <- (element.matrix$elC/element.matrix$elN)
C.N.1[is.infinite(C.N.1)]<- NA
mean(C.N.1,na.rm=T)
sd(C.N.1,na.rm=T)
C.P.1 <- element.matrix$elC/element.matrix$elP
C.P.1[is.infinite(C.P.1)]<- NA
mean(C.P.1,na.rm=T)
sd(C.P.1,na.rm=T)
length(which(is.na(C.N.1)==F))
length(which(is.na(C.P.1)==F))
mean(mass.list.1[Cluster.1])
sd(mass.list.1[Cluster.1])
C.N.Tot.1 <- sum(element.matrix$elC)/sum(element.matrix$elN)
C.P.Tot.1 <- sum(element.matrix$elC)/sum(element.matrix$elP)
############################################################
Cluster.2 <- which (Cluster[,1]==2)
Cluster.2 <- match(mass.list.1[Cluster.2],paf[,1])
element.matrix<-paf[Cluster.2,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))

setEPS()
postscript("Cluster2-VK.eps")
plot.vankrevelen(element.matrix=element.matrix,plotting.order="juggle",
                 main.hl=paste("Cluster 2 (Yellow) in VK space (",nrow(element.matrix),")",sep=""),cex=1)
dev.off()
AIC2<- (1+element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-0.5*element.matrix$elH)/(element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-element.matrix$elN-element.matrix$elP)
mean(AIC2)
sd (AIC2)
DBE2 <- element.matrix$elC - (element.matrix$elH/2)+(element.matrix$elN/2)+1
mean(DBE2)
sd(DBE2)
C.N.2 <- element.matrix$elC/element.matrix$elN
C.N.2[is.infinite(C.N.2)]<- NA
mean(C.N.2,na.rm=T)
sd(C.N.2,na.rm=T)
C.P.2 <- element.matrix$elC/element.matrix$elP
C.P.2[is.infinite(C.P.2)]<- NA
mean(C.P.2,na.rm=T)
sd(C.P.2,na.rm=T)
length(which(is.na(C.N.2)==F))
length(which(is.na(C.P.2)==F))
mean(mass.list.1[Cluster.2])
sd(mass.list.1[Cluster.2])
C.N.Tot.2 <- sum(element.matrix$elC)/sum(element.matrix$elN)
C.P.Tot.2 <- sum(element.matrix$elC)/sum(element.matrix$elP)
############################################################
Cluster.3 <- which (Cluster[,1]==3)
ct.I.minus.double.cl3 <- ct.I.minus.double[Cluster.3,]
Cluster.3 <- match(mass.list.1[Cluster.3],paf[,1])
element.matrix<-paf[Cluster.3,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))

setEPS()
postscript("Cluster3-VK.eps")
plot.vankrevelen(element.matrix=element.matrix,plotting.order="juggle",
                 main.hl=paste("Cluster 3 (Blue) in VK space (",nrow(element.matrix),")",sep=""),cex=1)
dev.off()

AIC3<- (1+element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-0.5*element.matrix$elH)/(element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-element.matrix$elN-element.matrix$elP)
mean(AIC3)
sd (AIC3)
DBE3 <- element.matrix$elC - (element.matrix$elH/2)+(element.matrix$elN/2)+1
mean(DBE3)
sd(DBE3)
C.N.3 <- element.matrix$elC/element.matrix$elN
C.N.3[is.infinite(C.N.3)]<- NA
mean(C.N.3,na.rm=T)
sd(C.N.3,na.rm=T)
C.P.3 <- element.matrix$elC/element.matrix$elP
C.P.3[is.infinite(C.P.3)]<- NA
mean(C.P.3,na.rm=T)
sd(C.P.3,na.rm=T)
length(which(is.na(C.N.3)==F))
length(which(is.na(C.P.3)==F))
mean(mass.list.1[Cluster.3])
sd(mass.list.1[Cluster.3])
C.N.Tot.3 <- sum(element.matrix$elC)/sum(element.matrix$elN)
C.P.Tot.3 <- sum(element.matrix$elC)/sum(element.matrix$elP)
############################################################
Cluster.4 <- which (Cluster[,1]==4)
Cluster.4 <- match(mass.list.1[Cluster.4],paf[,1])
element.matrix<-paf[Cluster.4,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))

setEPS()
postscript("Cluster4-VK.eps")
plot.vankrevelen(element.matrix=element.matrix,plotting.order="juggle",
                 main.hl=paste("Cluster 4 (Green) in VK space (",nrow(element.matrix),")",sep=""),cex=1)
dev.off()

AIC4<- (1+element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-0.5*element.matrix$elH)/(element.matrix$elC-0.5*element.matrix$elO-element.matrix$elS-element.matrix$elN-element.matrix$elP)
mean(AIC4)
sd (AIC4)
DBE4 <- element.matrix$elC - (element.matrix$elH/2)+(element.matrix$elN/2)+1
mean(DBE4)
sd(DBE4)
C.N.4 <- element.matrix$elC/element.matrix$elN
C.N.4[is.infinite(C.N.4)]<- NA
mean(C.N.4,na.rm=T)
sd(C.N.4,na.rm=T)
C.P.4 <- element.matrix$elC/element.matrix$elP
C.P.4[is.infinite(C.P.4)]<- NA
mean(C.P.4,na.rm=T)
sd(C.P.4,na.rm=T)
length(which(is.na(C.N.4)==F))
length(which(is.na(C.P.4)==F))
mean(mass.list.1[Cluster.4])
sd(mass.list.1[Cluster.4])
C.N.Tot.4 <- sum(element.matrix$elC)/sum(element.matrix$elN)
C.P.Tot.4 <- sum(element.matrix$elC)/sum(element.matrix$elP)
############################################################
Cluster.5 <- which (Cluster[,1]==5)
ct.I.minus.double.cl5 <- ct.I.minus.double[Cluster.5,]
Not.0.T1 <- which(rowSums(ct.I.minus.double.cl5[,Time.1])!=0)
Not.0.T2 <- which(rowSums(ct.I.minus.double.cl5[,Time.2])!=0)
Not.0.T3 <- which(rowSums(ct.I.minus.double.cl5[,Time.3])!=0)
Not.0.T4 <- which(rowSums(ct.I.minus.double.cl5[,Time.4])!=0)
Not.0.T5 <- which(rowSums(ct.I.minus.double.cl5[,Time.5])!=0)

Cluster.5 <- match(mass.list.1[Cluster.5],paf[,1])
element.matrix<-paf[Cluster.5,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))

setEPS()
postscript("Cluster5-VK.eps")
plot.vankrevelen(element.matrix=element.matrix,plotting.order="juggle",
                 main.hl=paste("Cluster 5 (Black) in VK space (",nrow(element.matrix),")",sep=""),cex=1)
dev.off()

AIC5<- (1+element.matrix$elC-(0.5*element.matrix$elO)-element.matrix$elS-(0.5*element.matrix$elH))/(element.matrix$elC-(0.5*element.matrix$elO)-element.matrix$elS-element.matrix$elN-element.matrix$elP)
mean(AIC5)
sd (AIC5)
DBE5 <- element.matrix$elC - (element.matrix$elH/2)+(element.matrix$elN/2)+1
mean(DBE5)
sd(DBE5)
C.N.5 <- element.matrix$elC/element.matrix$elN
C.N.5[is.infinite(C.N.5)]<- NA
mean(C.N.5,na.rm=T)
sd(C.N.5,na.rm=T)
C.P.5 <- element.matrix$elC/element.matrix$elP
C.P.5[is.infinite(C.P.5)]<- NA
mean(C.P.5,na.rm=T)
sd(C.P.5,na.rm=T)
length(which(is.na(C.N.5)==F))
length(which(is.na(C.P.5)==F))
mean(mass.list.1[Cluster.5])
sd(mass.list.1[Cluster.5])
C.N.Tot.5 <- sum(element.matrix$elC)/sum(element.matrix$elN)
C.P.Tot.5 <- sum(element.matrix$elC)/sum(element.matrix$elP)



###########Total############################
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))


AIC5<- (1+element.matrix$elC-(0.5*element.matrix$elO)-element.matrix$elS-(0.5*element.matrix$elH))/(element.matrix$elC-(0.5*element.matrix$elO)-element.matrix$elS-element.matrix$elN-element.matrix$elP)
Fck <- which(is.infinite(AIC5)==T)
mean(AIC5[-Fck])
sd (AIC5[-Fck])
DBE5 <- element.matrix$elC - (element.matrix$elH/2)+(element.matrix$elN/2)+1
mean(DBE5)
sd(DBE5)
C.N.5 <- element.matrix$elC/element.matrix$elN
C.N.5[is.infinite(C.N.5)]<- NA
mean(C.N.5,na.rm=T)
sd(C.N.5,na.rm=T)
C.P.5 <- element.matrix$elC/element.matrix$elP
C.P.5[is.infinite(C.P.5)]<- NA
mean(C.P.5,na.rm=T)
sd(C.P.5,na.rm=T)
length(which(is.na(C.N.5)==F))
length(which(is.na(C.P.5)==F))
mean(mass.list.1)
sd(mass.list.1)
C.N.Tot.5 <- sum(element.matrix$elC)/sum(element.matrix$elN)
C.P.Tot.5 <- sum(element.matrix$elC)/sum(element.matrix$elP)
boxplot(AIC1[-Fck],AIC2,AIC3,AIC4,AIC5)
t.test(AIC1[-Fck],AIC3)

