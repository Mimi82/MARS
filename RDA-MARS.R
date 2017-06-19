#MARS optics analyses
#The idea is to try to do something like in Graeber et al 2012 or to get inspire by it at least#


rm(list=ls(all=TRUE))

library(shape)
library (permute)
library(vegan)
setwd(homepath<- "/home/fonvielle/Documents/Mars/FT/")
list.files()

#get the data
data <- read.table("env.mars.MAg.csv",header=T,sep=",")
names(data)

# Set the treatment (Indices, components)
DOC_caracterisation<-data[,c(19,21,24,27,29,30,37:41)]

names(DOC_caracterisation)
str(DOC_caracterisation)

which(is.na(DOC_caracterisation),arr.ind=TRUE) # Check for NA but no NA here so we are fine

# check potential outliers and correlations
plot(DOC_caracterisation) 

## Second round of outlier removal ##

which(DOC_caracterisation$E2.to.E3==max(DOC_caracterisation$E2.to.E3))
DOC_caracterisation <- DOC_caracterisation[-25,]
data <- data [-25,]
plot(DOC_caracterisation) 

## Rename column name for plotting nice name... # 
colnames(DOC_caracterisation) <- c("SUVA", expression(S[275-295]), "E2 to E3", "FIX", "HIX", "Beta/Alpha", expression(C[1]), expression(C[2]), expression(C[3]), expression(C[4]), expression(C[5]))

################
#Start to normalise the data to compute the RDA
zDOC <- scale(DOC_caracterisation, center=T,scale=T) #z standardisation
zDOC


# Color (Time or nutrients) and symbol (Humics) for the Mars experiment. 
#Colors for time
colors.time<- as.factor(data$Time)
levels(colors.time)
colfunc <- colorRampPalette(c("cyan", "darkmagenta", "darkred"))
colfunc(9)
levels(colors.time)<-c(colfunc(9))
colors.time<-(as.character(colors.time))

#Color for nutrients
colors.Nutrients<- as.factor(data$Nutrients)
levels(colors.Nutrients)
colfunc <- colorRampPalette(c("cyan", "darkmagenta", "darkred"))
colfunc(7)
levels(colors.Nutrients)<-c(colfunc(7))
colors.Nutrients<-(as.character(colors.Nutrients))

#Color for Chla
colors.chloro<- as.factor(data$Chla)
levels(colors.chloro)
colfunc <- colorRampPalette(c("cyan", "darkmagenta", "darkred"))
colfunc(162)
levels(colors.chloro)<-c(colfunc(162))
colors.chloro<-(as.character(colors.chloro))


Color.Humin <- as.factor(data$HuminFeed)
levels(Color.Humin)
colfunc <- colorRampPalette(c("cadetblue2","yellow","tomato2"))
colfunc(3)
levels(Color.Humin) <- c(colfunc(3))
Color.Humin <- as.character(Color.Humin)


#Symbol for Humics
#A :pch=21, B: pch=23, C: pch=25, L: pch=22
Symbol.humic<- data$HuminFeed 
levels(Symbol.humic)<-c("21","23","25","22")
Symbol.humic<-as.numeric(as.character(Symbol.humic))



# compute the rda
HF <- as.numeric(data$HuminFeed)
Nutrients <- (data$Nutrients)
Time <- (data$Time)
treatment <- data.frame(HF,Nutrients, Time)
Partitionning <- varpart(zDOC,~HF,~Nutrients,~Time, data=treatment)
plot (Partitionning,digit=2)

rda<-rda(zDOC~HF+Nutrients+Time, data=treatment)

summary (rda)

#Test the rda
anova(rda)
anova(rda,first=TRUE)
anova(rda,by="axis",model="direct",data=data,perm.max=9999,step=1000) # The two first axis are significant 
anova(rda,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)
anova(rda,by="margin",model="direct",data=data,perm.max=9999,step=1000) # Still significant

#Even more test
permutest(rda)
anova.cca(rda)

#Plotting
plot(rda, scaling=1) #Tp is annoying for scaling
(species<-scores(rda,choices=c(1,2),display="sp",scaling=1))
text(species[,1:2],label=rownames(species),pos=4,cex=0.6)

#Plot the score with nice shape and arrows
#for the constraint the value is given by the summary of the rda (at the begining of the result)
(sites<-scores(rda,choices=c(1,2),display="sites",scaling=1,const=6.6048)) 
(constraints<-scores(rda,choices=c(1,2),display="bp",scaling=1))


#### cDOM is colored
plot(sites*2,pch=Symbol.humic,bg=Color.Humin,
     ylim=c(-2.5,2.5), xlim=c(-3,3),cex=2)
constraints <- constraints*2
Arrows(x0=0,y0=0,x1=constraints[,1],y1=constraints[,2],lwd=1.5,col="blue")
text(constraints[,1:2]*1.1,label=rownames(constraints),pos=4,cex=0.8,col="blue")

Arrows(x0=0,y0=0,x1=species[,1]/1.5,y1=species[,2]/1.5,lwd=1,arr.length=0) # Divided by 2.5 for scalling reason
text(species[,1:2]/1.5,label=rownames(species),pos=4,cex=0.6) # Divided by 1.5 for scalling

##### Nutrients are colored
plot(sites,pch=Symbol.humic,bg=colors.Nutrients,ylim=c(-2.5,2.5), xlim=c(-3,3))
Arrows(x0=0,y0=0,x1=constraints[,1],y1=constraints[,2],lwd=1.5,col="blue")
text(constraints[,1:2]*1.1,label=rownames(constraints),pos=4,cex=0.8,col="blue")

Arrows(x0=0,y0=0,x1=species[,1]/2,y1=species[,2]/2,lwd=1,arr.length=0) # Divided by 2.5 for scalling reason
text(species[,1:2]/2,label=rownames(species),pos=4,cex=0.6) # Divided by 2 for scalling

'
#Export in EPS#
setEPS()
postscript("RDA2-MARS1-2.eps")
plot(sites*2,pch=Symbol.humic,bg=Color.Humin,
     ylim=c(-2.5,2.5), xlim=c(-3,3),cex=2)
#constraints <- constraints*2
Arrows(x0=0,y0=0,x1=constraints[,1],y1=constraints[,2],lwd=1.5,col="blue")
text(constraints[,1:2]*1.1,label=rownames(constraints),pos=4,cex=0.8,col="blue")

Arrows(x0=0,y0=0,x1=species[,1]/1.5,y1=species[,2]/1.5,lwd=1,arr.length=0) # Divided by 2.5 for scalling reason
text(species[,1:2]/1.5,label=rownames(species),pos=4,cex=0.6) # Divided by 2 for scalling

dev.off()
'
