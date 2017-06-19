#FT-MS analysis

rm (list=ls(all=T))
library(gtools) 
library(vegan)
library(gam)
library(gamlss)
library(minpack.lm)
library(shape)
library(dataframes2xls)
library (MASS)

#round<-"_PreRecal"
#round<-"_RecalRd1"
round<-"_RecalRd2"

MagAbsMode<-"MagM"
#Get the data
#setwd(homepath<-"/Users/gabrielsinger/COLLABORATIONS/Jeremy Fonvielle/FTICRMS-MARS/")
setwd(homepath<- "/home/fonvielle/Documents/Mars/FT/")
source("Rcode/00.FTICRMS_R.functions.R")

#Get the original cross tab
load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.I999",".Rdata",sep=""))
#Get the "original" mass list
load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.mean999",".Rdata",sep=""))


nrow(ct.I)
BLK <- ct.I[,c(1:3)]
ct.I <- cbind(ct.I[,-c(1:3)],c(BLK))
colnames(ct.I)

#Remove double/triple analysed samples. 
#Choice between one or another is made based on the excel file. 
#If we did it twice, it was for a reason and one should be better than the other one
#If not, the second is selected in agreement with the student-learning effect. 
ct.I.minus.double <- ct.I[,-12]#D1B4
colnames(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-19]#D1C4
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-c(22,24)]#D2A1,D2A2
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-30]#D2B3
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-34]#D2C1
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-38]#D2C5
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-39]#D2C5
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-40]#D2C5
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-66]#D4A6
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-84]#DFA3
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-c(81,89)]#DFB1
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-96]#DFC2
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-99]#DFC5
names(ct.I.minus.double)

#Check if we have the right number
#Only 100 samples were processed
length(names(ct.I.minus.double[,-c(102:128)]))/21 #  Here, it should not be higher than 5 (should be 4.76)

#Remove photodegradation and mastermix
Tab.samples.Int<- ct.I.minus.double[,c(1:101)]
names(Tab.samples.Int)
head(Tab.samples.Int)

#Get the blanks (Why they are not with the other samples BTW?)
BLK.Samples <- ct.I.minus.double[,c(126:128)]
BLK.Samples <- t(BLK.Samples)
rownames(BLK.Samples)
BLK.Int <- colMeans(BLK.Samples,na.rm=T)

#Transpose dataframe
ct.I.minus.double <- t(Tab.samples.Int)


#Bind the m/z
#We do it here to be sure that we remove the correct mass afterward and avoid mistakes.
(mz.list <- as.vector(ct.mean$mz))
ct.I.minus.double1<- rbind(mz.list,ct.I.minus.double)


#Get the number of time a peak was found
nr.matches <- rep(NA, times=length(mz.list))

f <- 1458
  for (f in 1:length(mz.list)){
    Toto <- which(is.na(ct.I.minus.double1[,f])==F)
    nr.matches[f]<-length(Toto)
  }

#Bind the number of hit 
ct.I.minus.double1 <- rbind(nr.matches,ct.I.minus.double1)

#Removed peaks not appearing at least 4 times in the dataset.
#Those peaks are most likely not ecologically relevant peaks. 
hit.sup.3 <- which(nr.matches>3) # get them
ct.I.minus.double1 <- ct.I.minus.double1[,c(hit.sup.3)]#Remove them.

#Get back the mass list
mass.list <- ct.I.minus.double1[2,]

#Remove the masslist and the nr.matches. 
#We don't need nr.matches anymore and we got the mass list of interest above
ct.I.minus.double1<- ct.I.minus.double1[-c(1,2),]


#Remove peaks from the blanks
ct.I.minus.double1[is.na(ct.I.minus.double1)]<-0 #Replace all the NA by 0.
which(BLK.Int>0)
ct.I.minus.double1 <- ct.I.minus.double1[,-c(which(BLK.Int>0))]
mass.list <- mass.list[-c(which(BLK.Int>0))]

#Sum normalised the intensities
i=1
length(ct.I.minus.double1[1,])
tab.norm <- matrix(0,ncol=length(ct.I.minus.double1[29,]),nrow=101)#Here we create a matrix to collect the normalised intensities. We took the sample 29 just to have a length...

for (i in 1:101){ 
  sum.correct <- ct.I.minus.double1[i,]/sum((ct.I.minus.double1[i,]))
  tab.norm[i,] <- sum.correct
}


dev.off()
tomate <- hclust(vegdist(tab.norm, method="bray",na.rm=T))
tab.norm <- tab.norm[-c(97),]#Sample 97 is more than an outlier!! Even for a student artifact, it seems too much. 
tomate <- hclust(vegdist(tab.norm, method="bray",na.rm=T))
plot (tomate)


#############################################################################
#############################################################################
###################                         #################################
################### Dataset detection limit #################################
###################                         #################################
#############################################################################
#############################################################################

Sum.row <- rowSums(ct.I.minus.double1,na.rm=T)
length(Sum.row)
worst <- (which(Sum.row==min(Sum.row)))#Seems like the worst sample is the number 97, shame on you DFC3! 
Sum.row[worst]

Worst.sample <- as.vector(tab.norm[worst,])
#We just get and remove peaks not in the sample to be able to get the minimal intensity afterward
Not.in.sample <- which(Worst.sample==0)
Worst.sample <- Worst.sample[-Not.in.sample]
#Here we get the minimal intensity
lowest <- which (Worst.sample==min(Worst.sample))
Lowest.intensity <- Worst.sample[lowest]

tab.norm[tab.norm<=Lowest.intensity] <- 0
Sum.col <- colSums(tab.norm,na.rm=T)
Not.present.peaks <- which(Sum.col==0)
length(Not.present.peaks)
tab.norm <- tab.norm[,-c(Not.present.peaks)]
mass.list <- mass.list[-c(Not.present.peaks)]
mass.list[1]
ncol(tab.norm)
length(mass.list)#They have to match, otherwise it means there is something wrong up there. 
  

#############################################################################
#############################################################################

env.data.MARS <- read.table("env.mars.MAg.csv",header=T,sep=",")#We read the environmental variables
env.data.MARS <- env.data.MARS[-97,]#Outlier removed
#env.data.MARS <- env.data.MARS[c(1:100),]
#site.wo.b <- which(env.data.MARS$HuminFeed=="B")
#env.data.MARS.wo.b <- env.data.MARS[-c(site.wo.b),]
#tab.norm.wo.b <- tab.norm[-c(site.wo.b),]

NMDS <- metaMDS(tab.norm, distance="bray",k=3)# We performe the nMDS
  
nmds.species <- NMDS$species
nmds.sites <- NMDS$points

#write.csv (cbind(nmds.sites[-c(site.wo.b),],env.data.MARS$HuminFeed,env.data.MARS$Nutrients), "Sites-nmds.csv")

mass.lists <- as.factor (mass.list)
length(levels((mass.lists)))
colfunc <- colorRampPalette(c("plum2","pink3","ivory","lightsalmon"))
colfunc(length(levels((mass.lists))))
levels(mass.lists) <- c(colfunc(length(levels((mass.lists)))))
mass.lists <- as.character(mass.lists)


###################################################################################
###################################################################################
############### We create color vectors  ##########################################
###################################################################################
###################################################################################
Color.Time <- as.factor(env.data.MARS.wo.b$Time)
levels(Color.Time)
colfunc <- colorRampPalette(c("cyan","darkmagenta","darkred"))
colfunc(5)
levels(Color.Time) <- c(colfunc(5))
Color.Time <- as.character(Color.Time)

Color.Humin <- as.factor(env.data.MARS$HuminFeed)
levels(Color.Humin)
colfunc <- colorRampPalette(c("green3","black","orange"))
colfunc(3)
levels(Color.Humin) <- c(colfunc(3))
Color.Humin <- as.character(Color.Humin)


Color.Tn <- as.factor(env.data.MARS.wo.b$Nutrients)
levels(Color.Tn)
colfunc <- colorRampPalette(c("blue","lightblue","pink","red"))
colfunc(7)
levels(Color.Tn) <- c(colfunc(7))
Color.Tn <- as.character(Color.Tn)



Symbol.humics <- as.factor(env.data.MARS$HuminFeed)
levels(Symbol.humics) <- c("22","23","25")
Symbol.humics <- as.numeric(as.character(Symbol.humics))

###################################################################################
###################################################################################
############### Time is colored  ##################################################
###################################################################################
###################################################################################
dev.off()
par(mar=c(5,5,5,5))
plot(nmds.species[,2],nmds.species[,3],pch=21,bg=mass.lists,xlim=c(-1.75,1.75),ylim=c(-1.8,1.75),cex=2,cex.lab=2,cex.axis=2,xlab="Axis2",ylab="Axis3")
par(new=T)
plot(nmds.sites[,2]*10,nmds.sites[,3]*10, pch=Symbol.humics, 
     bg=Color.Time,xlim=c(-2.15,1.95),ylim=c(-1.9,1.9),
     cex=2.5,cex.lab=2,cex.axis=2,xlab="Axis2",ylab="Axis3")

names(env.data.MARS.wo.b)
env.scale <-as.data.frame(((env.data.MARS[,c(9:17,19,20,24,37:41)])))
env.scale$Tn<-as.numeric(env.scale$Tn)

env.scale <- as.data.frame(scale(env.scale,center=T, scale=T))
names(env.scale)
length(nmds.sites[,1])
length(env.scale$DOC)

Chemo.all <- diversity(tab.norm,index="shannon")
plot(Chemo.all)

env.scale <- cbind(env.scale,Chemo.all)

fit <- envfit(NMDS,env.scale,, permutations = 9999,na.rm = T,choices = c(1:3)) 

P_val_fit <- as.vector((fit$vectors$pvals))
Low.p.value <- which(P_val_fit<=0.05 & P_val_fit>=0.01)
Strong.p.value <- which(P_val_fit<=0.01 & P_val_fit>!0.01)


#ordiArrowMul (scores(fit, fill = 1, display = "vectors"))
arrows <- scores(fit,display = "vectors")* 1
Arrows(x0 =0,y0 =0,x1=arrows[c(Low.p.value),2],y1=arrows[c(Low.p.value),3],lwd=2,col="black",arr.length = 0.1,arr.width = 0.1)
text(x=arrows[c(Low.p.value),2]*1.2,y=arrows[c(Low.p.value),3]*1.2,label=c(names(env.scale[c(Low.p.value)])),cex=1)
Arrows(x0 =0,y0 =0,x1=arrows[c(Strong.p.value),2],y1=arrows[c(Strong.p.value),3],lwd=2,arr.length = 0.1,arr.width = 0.1,col="navy")
text(x=arrows[c(Strong.p.value),2]*1.2,y=arrows[c(Strong.p.value),3]*1.2,label=c(names(env.scale[c(Strong.p.value)])),col="navy",cex=1)





##############################################################################################################
############# Collect mole from A, B, C ######################################################################
##############################################################################################################

A <- which(as.numeric(env.data.MARS$HuminFeed)==1)
ct.I.A <- tab.norm[A,]
#Put 0 instead of NA
ct.I.A[is.na(ct.I.A)]<- 0
ct.I.A <- colSums(ct.I.A)
(not.present.in.A <- which(ct.I.A==0))#Must be different than 0 and not contain NA
length(not.present.in.A)
mass.list.1.A <- mass.list[-c(not.present.in.A)]


B <- which(as.numeric(env.data.MARS$HuminFeed)==2)
ct.I.B <- tab.norm[B,]
#Put 0 instead of NA
ct.I.B[is.na(ct.I.B)]<- 0
ct.I.B <- colSums(ct.I.B)
(not.present.in.B <- which(ct.I.B==0))#Must be different than 0 and not contain NA
mass.list.1.B <- mass.list [-c(not.present.in.B)]


C <- which(as.numeric(env.data.MARS$HuminFeed)==3)
ct.I.C <- tab.norm[C,]
#Put 0 instead of NA
ct.I.C[is.na(ct.I.C)]<- 0
ct.I.C <- colSums(ct.I.C)
(not.present.in.C <- which(ct.I.C==0))#Must be different than 0 and not contain NA
mass.list.1.C <- mass.list [-c(not.present.in.C)]


#########Only the C molecules#################################################################
C.a.Mole <- subset (mass.list.1.C, !(mass.list.1.C %in% mass.list.1.A))
C.Mole <- subset (C.a.Mole, !(C.a.Mole %in% mass.list.1.B))
C.plot.mole <- match(C.Mole,mass.list)

#########Only the A molecules###################################################################
A.c.Mole <- subset (mass.list.1.A, !(mass.list.1.A %in% mass.list.1.C))
A.Mole <- subset (A.c.Mole, !(A.c.Mole %in% mass.list.1.B))
A.plot.mole <- match(A.Mole,mass.list)

#########Only the A molecules###################################################################
B.c.Mole <- subset (mass.list.1.B, !(mass.list.1.B %in% mass.list.1.C))
B.Mole <- subset (B.c.Mole, !(B.c.Mole %in% mass.list.1.A))
B.plot.mole <- match(B.Mole,mass.list)


#setEPS(width=18, height=15)
#postscript("NMDS-MARS.eps")
par(mar=c(7,7,7,7))
plot(nmds.species[,2],nmds.species[,3],pch=21,
     bg="grey87",col="grey79",xlim=c(-1.75,1.75),ylim=c(-1.8,1.75),
     cex=2,cex.lab=3,cex.axis=3,xlab="Axis2",ylab="Axis3")
par(new=T)
plot(nmds.species[B.plot.mole,2],nmds.species[B.plot.mole,3],pch=21,
     bg="yellow",col="grey7",xlim=c(-1.75,1.75),ylim=c(-1.8,1.75),
     cex=2,cex.lab=3,cex.axis=3,xlab="Axis2",ylab="Axis3")
par(new=T)
plot(nmds.species[C.plot.mole,2],nmds.species[C.plot.mole,3],pch=21,
     bg="tomato2",col="grey7",xlim=c(-1.75,1.75),ylim=c(-1.8,1.75),
     cex=2,cex.lab=3,cex.axis=3,xlab="Axis2",ylab="Axis3")
par(new=T)
plot(nmds.species[A.plot.mole,2],nmds.species[A.plot.mole,3],pch=21,
     bg="cadetblue2",col="grey7",xlim=c(-1.75,1.75),ylim=c(-1.8,1.75),
     cex=2,cex.lab=3,cex.axis=3,xlab="Axis2",ylab="Axis3")

Arrows(x0 =0,y0 =0,x1=arrows[c(Low.p.value),2]*1.5,y1=arrows[c(Low.p.value),3]*1.5,lwd=2,col="black",arr.length = 0.5,arr.width = 0.3)
text(x=arrows[c(Low.p.value),2]*1.7,y=arrows[c(Low.p.value),3]*1.7,label=c(names(env.scale[c(Low.p.value)])),cex=1)
Arrows(x0 =0,y0 =0,x1=arrows[c(Strong.p.value),2]*1.6,y1=arrows[c(Strong.p.value),3]*1.6,lwd=2,arr.length = 0.5,arr.width = 0.3,col="navy")
text(x=arrows[c(Strong.p.value),2]*1.8,y=arrows[c(Strong.p.value),3]*1.8,label=c(names(env.scale[c(Strong.p.value)])),col="navy",cex=1)
legend(x=0.8,y=-1.35,legend = "Stress= 0.04",bty="n",cex=2.5)
#dev.off()
