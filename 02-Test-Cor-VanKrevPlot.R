#Correspondance nMDS and molForm

library(gtools) 
library(vegan)
library(gam)
library(gamlss)
library(minpack.lm)
library(shape)
library(dataframes2xls)
library (MASS)

#Put some names on object to make file reading easier

#round<-"_PreRecal"
#round<-"_RecalRd1"
round<-"_RecalRd2"
MagAbsMode<-"MagM"


#Get the data
#setwd(homepath<-"/Users/gabrielsinger/COLLABORATIONS/Jeremy Fonvielle/FTICRMS-MARS/")#This is Gabriel's path
setwd(homepath<- "/home/fonvielle/Documents/Mars/FT/")#This is Jeremy's IGB computer path
source("Rcode/00.FTICRMS_R.functions.R")

load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.mean.iso.FIN",round,".Rdata",sep=""))#CT.MEAN file
load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.I.iso.FIN",round,".Rdata",sep="")) # is same (!) #CT.I FILE
load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.mz.iso.FIN",round,".Rdata",sep="")) #CT.MZ FILE
#load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.SN.iso.FIN",round,".Rdata",sep=""))
#load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.ResPow.iso.FIN",round,".Rdata",sep=""))
#load(paste("Rdata/",MagAbsMode,".FormulaeAssignment_parms0b/ct.mySN.iso.FIN",round,".Rdata",sep=""))
load(paste("Rdata/MagM.FormulaeAssignment_parms0b/paf.05.RefTol.iso.duprem.blkrem.FIN_RecalRd2.Rdata",sep=""))#PAF FILE
paf

#We remove the samples analysed twice (same as in script 1 done before)
colnames(ct.I)
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
ct.I.minus.double <- ct.I.minus.double[,-c(88)]#DFB1
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-96]#DFC2
names(ct.I.minus.double)
ct.I.minus.double <- ct.I.minus.double[,-98]#DFC5
names(ct.I.minus.double)

ct.photo <- ct.I.minus.double[,c(104:123)]
ct.I.minus.double <- ct.I.minus.double[,c(1:100)]#Only the samples, no MasterMix no BLK, no Photodegradation


#We get the mass list
mass.list <- as.vector(ct.mean$MeanExpMass)
length(mass.list)


######### We match that to the paf file ########################################################
Not.interested <- which(is.na(match(mass.list,paf$MeanExpMass))==T)#The peaks that have nothing to do here, get out of the dataset =)
mass.list.1 <- mass.list[-c(Not.interested)]
which(is.na(match(mass.list.1,paf$MeanExpMass))==T)#If >0 Something is wrong
#Those mole are not in the dataset but are in the paf list...
Not.interested.2 <- which(is.na(match(paf$MeanExpMass,mass.list.1))==T)# Molecules present here and not in the PAF file are artifact from molecular formulas assignement
ct.I.minus.double <- ct.I.minus.double[-c(Not.interested),]
paf <- paf[-c(Not.interested.2),]
length(mass.list.1)
ct.I.minus.double[is.na(ct.I.minus.double)]<-0
write.csv(ct.I.minus.double,"CT-I-MARS.csv")

ct.I.minus.double[is.na(ct.I.minus.double)]<-0
Ct.scale <- t(ct.I.minus.double)
Ct.scale <-  scale(tab.norm.wo.b,center=T,scale=T)
head(tab.norm)

##############################################################################################################
############# Collect mole from A, B, C ######################################################################
############# Sing Mickael Jackson "ABC..." ##################################################################
############# Dumb repetitive coding as to be cleaned ########################################################
############## Next time do it without music #################################################################
##############################################################################################################

A <- which(as.numeric(env.data.MARS$HuminFeed)==1)
ct.I.A <- ct.I.minus.double[,A]
#Put 0 instead of NA
ct.I.A[is.na(ct.I.A)]<- 0
ct.I.A <- rowSums(ct.I.A)
(not.present.in.A <- which(ct.I.A==0))#Must be different than 0 and not contain NA
mass.list.1.A <- mass.list.1[-c(not.present.in.A)]


B <- which(as.numeric(env.data.MARS$HuminFeed)==2)
ct.I.B <- ct.I.minus.double[,B]
#Put 0 instead of NA
ct.I.B[is.na(ct.I.B)]<- 0
ct.I.B <- rowSums(ct.I.B)
(not.present.in.B <- which(ct.I.B==0))#Must be different than 0 and not contain NA
mass.list.1.B <- mass.list.1[-c(not.present.in.B)]


C <- which(as.numeric(env.data.MARS$HuminFeed)==3)
ct.I.C <- ct.I.minus.double[,C]
#Put 0 instead of NA
ct.I.C[is.na(ct.I.C)]<- 0
ct.I.C <- rowSums(ct.I.C)
(not.present.in.C <- which(ct.I.C==0))#Must be different than 0 and not contain NA
mass.list.1.C <- mass.list.1[-c(not.present.in.C)]


#########Only the C molecules#################################################################
C.a.Mole <- subset (mass.list.1.C, !(mass.list.1.C %in% mass.list.1.A))
C.Mole <- subset (C.a.Mole, !(C.a.Mole %in% mass.list.1.B))

#########Only the A molecules###################################################################
A.c.Mole <- subset (mass.list.1.A, !(mass.list.1.A %in% mass.list.1.C))
A.Mole <- subset (A.c.Mole, !(A.c.Mole %in% mass.list.1.B))

#########Only the A molecules###################################################################
A.c.Mole <- subset (mass.list.1.A, !(mass.list.1.A %in% mass.list.1.C))
A.Mole <- subset (A.c.Mole, !(A.c.Mole %in% mass.list.1.B))



### --> This is very disapointing, barely no molecule with a formula can be identify as belonging to only C or to only A
## =( 

###################################################
###################################################
### We will do the correlation and see ############
###################################################
###################################################
env.data.MARS <- read.table("env.mars.MAg.csv",header=T,sep=",")
env.data.MARS <- env.data.MARS[-97,]
#env.data.MARS <- env.data.MARS[-c(site.wo.b),]
which(env.data.MARS$HuminFeed=="A" & env.data.MARS$Nutrients==1 &env.data.MARS$Time)
#env.data.MARS <- env.data.MARS [-55,]

length(ct.I.minus.double[,1])
length(paf[,1])
length(C1.cor.Rho)

names(paf)
paf$DBE<-1+0.5*(2*paf$C-paf$H+paf$N+paf$P) # double-bond equivalents
paf$HC<-paf$H/paf$C
paf$OC<-paf$O/paf$C
paf$meas.err<-NULL # not needed

element.matrix<-paf[,4:10]
names(element.matrix)
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
#ct.I.minus.double <- ct.I.minus.double[,-c(site.wo.b)]
#ct.I.minus.double <- ct.I.minus.double [,-55]


ct.I.minus.double.t <- t(ct.I.minus.double)


Tata <- vegdist(ct.I.minus.double.t,method="bray")
plot(hclust(Tata))

i=1
length(ct.I.minus.double.t[1,])
tab.norm <- matrix(0,ncol=length(ct.I.minus.double.t[35,]),nrow=100)#Here we create a matrix to collect the normalised intensities. We took the sample 29 just to have a length...

for (i in 1:100){ 
  sum.correct <- ct.I.minus.double.t[i,]/(sum((ct.I.minus.double.t[i,])))
  tab.norm[i,] <- sum.correct
}

rownames(tab.norm)<-rownames(ct.I.minus.double.t) 

Div.A <- diversity (tab.norm[A,],index="shannon")
Div.B <- diversity (tab.norm[B,],index="shannon")
Div.C <- diversity (tab.norm[C,],index="shannon")

boxplot(Div.A,Div.B,Div.C,names = c("0","0.6","1.2"),
        xlab=expression("cDOM loaded (mg C " *L^-1 *")" ),
        ylab="Chemodiversity index",
        cex.lab=2,cex.axis=2,
        col=c("cadetblue2","yellow","tomato2"))

Tat0 <- vegdist(tab.norm,method="bray")
plot(hclust(Tat0))

MDS.PAF <- metaMDS(tab.norm,distance="bray",k=3)


MDS.PAF.Site <- MDS.PAF$points
MDS.PAF.species <- MDS.PAF$species

#Color.O/C
Color.OC <- factor(paf$OC)
length(levels(Color.OC))
Colfunt <- colorRampPalette(c("blue","cyan","pink","red"),bias=1,interpolate="linear")
Colfunt(208)
levels(Color.OC)<- c(Colfunt(208))
Color.OC <- as.character(Color.OC)



#Color.H/C
Color.HC <- factor(paf$HC)
length(levels(Color.HC))
Colfunt <- colorRampPalette(c("blue","cyan","pink","red"),bias=1,interpolate="linear")
Colfunt(290)
levels(Color.HC)<- c(Colfunt(290))
Color.HC <- as.character(Color.HC)

Color.Humin <- as.factor(env.data.MARS$HuminFeed)
levels(Color.Humin)
colfunc <- colorRampPalette(c("green3","white","orange"))
colfunc(3)
levels(Color.Humin) <- c(colfunc(3))
Color.Humin <- as.character(Color.Humin)

Symbol.humics <- as.factor(env.data.MARS$HuminFeed)
levels(Symbol.humics) <- c("22","23","25")
Symbol.humics <- as.numeric(as.character(Symbol.humics))

######### Correlation with C1
C1.cor <- rep (NA, times=length(tab.norm[1,]))
C1.cor.Rho <- rep (NA, times=length(tab.norm[1,]))
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,env.data.MARS$Fmax1norm,method = "spearman")
  C1.cor[i]<- x$p.value
  C1.cor <- C1.cor
  C1.cor.Rho[i] <- x$estimate
  C1.cor.Rho <- C1.cor.Rho
  
}

C1.cor.Rho[is.na(C1.cor.Rho)]<-0
Cor.C1.adjust <- p.adjust(C1.cor,method="bonferroni")
C1.Significant <- which(Cor.C1.adjust<0.05)
C1.cor.Rho.Sig <- C1.cor.Rho[c(C1.Significant)]
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
element.matrix <- element.matrix[c(C1.Significant),]

plot.vankrevelen(element.matrix=element.matrix,color.criterion=C1.cor.Rho.Sig,plotting.order="composition",
                 main.hl=paste("Correlation with C1 (",nrow(element.matrix),")",sep=""),cex=3)
dev.off()
######### Correlation with C3
C3.cor <- rep (NA, times=length(tab.norm[1,]))
C3.cor.Rho <- rep (NA, times=length(tab.norm[1,]))
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,env.data.MARS$Fmax3norm,method = "spearman")
  C3.cor[i]<- x$p.value
  C3.cor <- C3.cor
  C3.cor.Rho[i] <- x$estimate
  C3.cor.Rho <- C3.cor.Rho
  
}

C3.cor.Rho[is.na(C3.cor.Rho)]<-0
Cor.C3.adjust <- p.adjust(C3.cor,method="bonferroni")
C3.Significant <- which(Cor.C3.adjust<0.05)
C3.cor.Rho.Sig <- C3.cor.Rho[c(C3.Significant)]
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
element.matrix <- element.matrix[c(C3.Significant),]

plot.vankrevelen(element.matrix=element.matrix,color.criterion=C3.cor.Rho.Sig,plotting.order="juggle",
                 main.hl=paste("Correlation with C3 (",nrow(element.matrix),")",sep=""),cex=1)

############## Correlation with HIX ####################################
HIX.cor <- rep (NA, times=length(tab.norm[1,]))
HIX.cor.Rho <- rep (NA, times=length(tab.norm[1,]))
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,env.data.MARS$HIX2,method = "spearman")
  HIX.cor[i]<- x$p.value
  HIX.cor <- HIX.cor
  HIX.cor.Rho[i] <- x$estimate
  HIX.cor.Rho <- HIX.cor.Rho
  
}

HIX.cor.Rho[is.na(HIX.cor.Rho)]<-0
Cor.HIX.adjust <- p.adjust(HIX.cor,method="bonferroni")
Hix.Significant <- which(Cor.HIX.adjust<0.05)
HIX.cor.Rho.Sig <- HIX.cor.Rho[c(Hix.Significant)]
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=HIX.cor.Rho,plotting.order="juggle",
                 main.hl=paste("Correlation with HIX (",nrow(element.matrix),")",sep=""),cex=3,molgrpdef=TRUE)

################ Correlation with Chla #################################################################
Chla.cor <- rep (NA, times=length(tab.norm[1,]))
Chla.cor.Rho <- rep (NA, times=length(tab.norm[1,]))
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,env.data.MARS$Chla,method = "spearman")
  Chla.cor[i]<- x$p.value
  Chla.cor <- Chla.cor
  Chla.cor.Rho[i] <- x$estimate
  Chla.cor.Rho <- Chla.cor.Rho
  
}

Chla.cor.Rho[is.na(Chla.cor.Rho)]<-0
Cor.Chla.adjust <- p.adjust(Chla.cor,method="bonferroni")
Chla.Significant <- which(Cor.Chla.adjust<0.05)
Chla.cor.Rho.Sig <- Chla.cor.Rho[c(Chla.Significant)]
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=Chla.cor.Rho,plotting.order="juggle",
                 main.hl=paste("Correlation with Chla (",nrow(element.matrix),")",sep=""),cex=3)
element.matrix <- element.matrix[c(Chla.Significant),]

plot.vankrevelen(element.matrix=element.matrix,color.criterion=Chla.cor.Rho.Sig,plotting.order="juggle",
                 main.hl=paste("Correlation with Chla (",nrow(element.matrix),")",sep=""),cex=3)

##############################################################################################################
##############################################################################################################


Ox.cor <- rep (NA, times=length(tab.norm[1,]))
 Ox.cor.Rho <- rep (NA, times=length(tab.norm[1,]))
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,env.data.MARS[,15],method = "spearman")
   Ox.cor[i]<- x$p.value
   Ox.cor <-  Ox.cor
   Ox.cor.Rho[i] <- x$estimate
   Ox.cor.Rho <-  Ox.cor.Rho
  
}

 #par(mar=c(10,15,10,15))
 Ox.cor.Rho[is.na( Ox.cor.Rho)]<-0
Cor.Ox.adjust <- p.adjust( Ox.cor,method="bonferroni")
 Ox.Significant <- which(Cor.Ox.adjust<0.05)
 Ox.cor.Rho.Sig <-  Ox.cor.Rho[c( Ox.Significant)]
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
element.matrix <- element.matrix[c( Ox.Significant),]

plot.vankrevelen(element.matrix=element.matrix,color.criterion= Ox.cor.Rho.Sig,plotting.order="juggle",
                 main.hl=paste("Correlation with  Ox (",nrow(element.matrix),")",sep=""),cex=1.5,cex.lab=2,cex.axis=2)
text(x=0.1,y=2.49,label="n=429",cex=2)

##############################################################################################################
##############################################################################################################
Oxygen.saturation.cor <- rep (NA, times=length(tab.norm[1,]))
Oxygen.saturation.cor.Rho <- rep (NA, times=length(tab.norm[1,]))
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,env.data.MARS$Oxygen.saturation,method = "spearman")
  Oxygen.saturation.cor[i]<- x$p.value
  Oxygen.saturation.cor <- Oxygen.saturation.cor
  Oxygen.saturation.cor.Rho[i] <- x$estimate
  Oxygen.saturation.cor.Rho <- Oxygen.saturation.cor.Rho
  
}

Oxygen.saturation.cor.Rho[is.na(Oxygen.saturation.cor.Rho)]<-0
Cor.Oxygen.saturation.adjust <- p.adjust(Oxygen.saturation.cor,method="bonferroni")
Oxygen.saturation.Significant <- which(Cor.Oxygen.saturation.adjust<0.05)
Oxygen.saturation.cor.Rho.Sig <- Oxygen.saturation.cor.Rho[c(Oxygen.saturation.Significant)]
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=Oxygen.saturation.cor.Rho,plotting.order="composition",
                 main.hl=paste("Correlation with Oxygen.saturation (",nrow(element.matrix),")",sep=""),cex=3,molgrpdef=F)
element.matrix <- element.matrix[c(Oxygen.saturation.Significant),]
plot.vankrevelen(element.matrix=element.matrix,color.criterion="mol.grps",plotting.order="composition",
                 main.hl=paste("Correlation with Oxygen.saturation (",nrow(element.matrix),")",sep=""),cex=3,molgrpdef=F)
plot.vankrevelen(element.matrix=element.matrix,color.criterion=Oxygen.saturation.cor.Rho.Sig,plotting.order="composition",
                 main.hl=paste("Correlation with Oxygen.saturation (",nrow(element.matrix),")",sep=""),cex=3,molgrpdef=F)



dev.off()



