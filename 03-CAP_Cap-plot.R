#CAP Analysis

?capscale


##############################################
############ Humin Feed ######################
DOM.cap <- capscale(tab.norm~factor(as.numeric(env.data.MARS$HuminFeed)),dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000) 
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

CAP1.HF <- DOM.cap$CCA$u
CAP1.HF.1 <- DOM.cap$CCA$v

CAP1.HF.cor <- rep (NA, times=length(ct.I.minus.double[,1]))
CAP1.HF.cor.Rho <- rep (NA, times=length(ct.I.minus.double[,1]))
f<-1
for (i in 1:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,CAP1.HF[,1],method = "spearman")
  CAP1.HF.cor[i]<- x$p.value
  CAP1.HF.cor <- CAP1.HF.cor
  CAP1.HF.cor.Rho[i] <- x$estimate
  CAP1.HF.cor.Rho <- CAP1.HF.cor.Rho
  
}
dev.off()
par(mar=c(5,5,5,5))
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=CAP1.HF.cor.Rho,plotting.order="juggle",
                 main.hl=paste("Correlation with CAP1 (HF) (",nrow(element.matrix),")",sep=""),cex=1,molgrpdef=TRUE,cex.lab=3, cex.axis=2)



################################################
############ Nutrients #########################
################################################

DOM.cap <- capscale(tab.norm~env.data.MARS$Levels.Nut,dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000)
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

CAP1.N <- DOM.cap$CCA$u
CAP1.N.1 <- DOM.cap$CCA$v


CAP1.N.cor <- rep (NA, times=length(ct.I.minus.double[,1]))
CAP1.N.cor.Rho <- rep (NA, times=length(ct.I.minus.double[,1]))
i<-1
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,CAP1.N,method = "spearman")
  CAP1.N.cor[i]<- x$p.value
  CAP1.N.cor <- CAP1.N.cor
  CAP1.N.cor.Rho[i] <- x$estimate
  CAP1.N.cor.Rho <- CAP1.N.cor.Rho
  
}
CAP1.N.cor.Rho[is.na(CAP1.N.cor.Rho)==T]<-0
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=CAP1.N.cor.Rho,plotting.order="juggle",
                 main.hl=paste("Correlation with CAP1 (Nutrients) (",nrow(element.matrix),")",sep=""),cex=1,molgrpdef=TRUE)




################################################
############ Time ##############################
################################################

DOM.cap <- capscale(tab.norm~env.data.MARS$Time,dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000) 
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

CAP1.T <- DOM.cap$CCA$u
CAP1.T.1 <- DOM.cap$CCA$v


CAP1.T.cor <- rep (NA, times=length(ct.I.minus.double[,1]))
CAP1.T.cor.Rho <- rep (NA, times=length(ct.I.minus.double[,1]))
f<-1
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[,i])) ,CAP1.T,method = "spearman")
  CAP1.T.cor[i]<- x$p.value
  CAP1.T.cor <- CAP1.T.cor
  CAP1.T.cor.Rho[i] <- x$estimate
  CAP1.T.cor.Rho <- CAP1.T.cor.Rho
  
}

CAP1.T.cor.Rho[is.na(CAP1.T.cor.Rho)==T]<-0
element.matrix<-paf[,4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=CAP1.T.cor.Rho,plotting.order="juggle",
                 main.hl=paste("Correlation with CAP1 (",nrow(element.matrix),")",sep=""),cex=1,molgrpdef=TRUE)

plot(CAP1.T)



################################################
############ Time*HF ###########################

DOM.cap <- capscale(tab.norm~env.data.MARS$Time*(as.numeric(env.data.MARS$HuminFeed)),dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000) 
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

#CAP1.T*HF <- DOM.cap$CCA$u



################################################
############ HF*Nutrients ######################

DOM.cap <- capscale(tab.norm~env.data.MARS$Levels.Nut*env.data.MARS$Fmax1norm,dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000) 
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

CAP1.N.HF <- DOM.cap$CCA$u


################################################
############        All          ###############
################################################
DOM.cap <- capscale(tab.norm~as.factor(env.data.MARS$Nutrients)*as.factor(env.data.MARS$HuminFeed)+env.data.MARS$Time,dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000) 
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

CAP1.Alles <- DOM.cap$CCA$u


################################################
############ Nutrients only in A ###############
################################################

#Look at the nutrients effect only in the A series
A <- which(as.numeric(env.data.MARS$HuminFeed)==1)
Tan.A. <- env.data.MARS[A,]

DOM.cap <- capscale(tab.norm[A,]~Tan.A.$Levels.Nut,dist="bray",sqrt.dist= TRUE)
summary(DOM.cap)

anova.cca(DOM.cap)
anova(DOM.cap,by="axis",model="direct",data=data,perm.max=9999,step=1000) # The two first axis are significant 
anova(DOM.cap,by="terms",model="direct",data=data,perm.max=9999,step=1000)  # Same for the terms (carefull, here order matter)

CAP.A <- DOM.cap$CCA$wa
CAP.A.1 <- DOM.cap$CCA$v


CAP1.A.cor <- rep (NA, times=length(ct.I.minus.double[,1]))
CAP1.A.cor.Rho <- rep (NA, times=length(ct.I.minus.double[,1]))
f<-1
for (i in 2:length(tab.norm[1,])){
  x <- cor.test (as.vector(as.numeric(tab.norm[A,i])) ,CAP.A,method = "spearman")
  CAP1.A.cor[i]<- x$p.value
  CAP1.A.cor <- CAP1.A.cor
  CAP1.A.cor.Rho[i] <- x$estimate
  CAP1.A.cor.Rho <- CAP1.A.cor.Rho
  
}


CAP1.A.cor.Rho[is.na(CAP1.A.cor.Rho)==T]<-0
which(CAP1.A.cor<0.001)

element.matrix<-paf[which(CAP1.A.cor<0.001),4:10]
(names(element.matrix)<-paste("el",names(element.matrix),sep=""))
plot.vankrevelen(element.matrix=element.matrix,color.criterion=CAP1.A.cor.Rho[c(which(CAP1.A.cor<0.001))],plotting.order="juggle",
                 main.hl=paste("Correlation with Nutrients (Only A series) (",nrow(element.matrix),")",sep=""),cex=2,molgrpdef=TRUE)


