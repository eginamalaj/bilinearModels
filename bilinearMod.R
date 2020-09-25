#################################################################################################################
#
#
#     R code for the bilinear model calculations presented in:
#
#     Evolutionary patterns and physicochemical properties explain macroinvertebrate sensitivity to heavy metals
#
#     Egina Malaj, Guillaume Guenard, Ralf B. Schafer, Peter C. von der Ohe
#       
#     For questions contact: guillaume.guenard@gmail.com or eginamalaj@gmail.com
#     
#
#
##################################################################################################################
#
# rm(list = ls())
#
# Libraries Used
pkg<- c("doBy","reshape", "MPSEM","phytools") 
#
# Call packages
lapply(pkg, require, character.only=T)
#
#
#
#####
##### Load Data
#####
#
# Functions to calculate the AIC and prediction coefficient (here psquare).
source("Functions.R")
#
# Phylogenetic tree
load("data_tree.rda")
#
# Heavy metal data - Source: Malaj et al. (2012)
load("data_metals.RData")
dt_tox <- merge(dt_out_fin,Taxa[,c("Taxon","DNA_source")],by="Taxon")
#
# Add short names from the tree
dt_tox2 <- cbind(dt_tox,ShortNames=spnames[match(dt_tox[,"DNA_source"],spnames[,1L]),2L])
#
# Water chemistry
load("data_physicochem.RData")
#
# Attach W.Chemistry
wat_chem<-merge(physo2, dt_tox2, by.x="Test_Number", by.y="Test_Number", all.y=T)
nrow(wat_chem)==nrow(dt_tox2)
#
#
# Nr.Species for which chemistry is reported
#
wat_chem2<-wat_chem[!is.na(wat_chem$Temp)&!is.na(wat_chem$Hardness)&!is.na(wat_chem$pH),]
wat_chem3<-summaryBy(data=wat_chem2, Conc_mmgL~Taxon+Chemical)
wat_chem4<- reshape(wat_chem3, direction="wide",timevar="Chemical",idvar = "Taxon",v.names ="Conc_mmgL.mean")
wat_chem4 # 4 species for the 4 HM model.
nrow(wat_chem[!is.na(wat_chem$Temp),])/nrow(wat_chem) # 78 % have temp reported
nrow(wat_chem[!is.na(wat_chem$Hardness),])/nrow(wat_chem)# 45 % of the entries have hardness 
nrow(wat_chem[!is.na(wat_chem$pH),])/nrow(wat_chem)# 38 % have pH
#
# Restriction pH: 6.5 to 7.5
wat_chem5<-wat_chem[!is.na(wat_chem$pH) & wat_chem$pH<7.66 & wat_chem$pH>6.48,]
wat_chem6<-wat_chem[is.na(wat_chem$pH),]; wat_chem7<-rbind(wat_chem6,wat_chem5)
#
#
# Table for properties
#
# Add exposure duration which was always present in the dataset
load("data_exposure.RData")
wat_chem8<-merge(db_phylo, wat_chem7, by.x="Test_Number", by.y="Test_Number", all.y=T)
nrow(wat_chem8)==nrow(wat_chem7)
#
tab_wat_chem<-summaryBy(Temp+Hardness+pH+Duration~Taxon+Chemical, FUN=c(min,max), na.rm=TRUE, data=wat_chem8)
#
# Remove two species outside of the ph range
tab_wat_chem2<-tab_wat_chem[!(tab_wat_chem$Taxon=="Aphelenchus avenae") & !(tab_wat_chem$Taxon=="Chironomus tentans"),]
#
# Table for water chemistry and exposure time
#write.csv(tab_wat_chem2, file='Table_range-prop.csv') 
#
# 2 species are removed as a result of the pH restrictions. Remove the water chemistry columns
Chemical.data<-wat_chem7[!(wat_chem7$Taxon=="Aphelenchus avenae") & !(wat_chem7$Taxon=="Chironomus tentans"),-c(2:6)]
#
unique(Chemical.data$DNA_source)
#
#
# Table S3 paper
chem_mean<-Chemical.data
chem_mean$Conc_mmgL<-log10(chem_mean$Conc_mmgL)
chem_mean1<-summaryBy(data=chem_mean, Conc_mmgL~DNA_source+Chemical, FUN=c(mean,sd,length))
chem_mean2<-reshape(chem_mean1, idvar = "DNA_source", timevar = "Chemical", direction="wide")
chem_mean2[-1]<-round(chem_mean2[-1],3)
#write.csv(chem_mean2, file='TableS3.csv') 
#
# Chemical properties 
Chemical.properties<- read.csv("data_Properties.csv",sep=",", header=T) 
Chemical.properties[,"Chemical"] <- as.character(Chemical.properties[,"Chemical"])
#
#####
#####
#####
#
# Build the tree
#
ModelSpecies <- as.character(unique(Taxa[match(Chemical.data[,"Taxon"],Taxa[,"Taxon"]),"DNA_source"]))
treeChemical <- list()
treeChemical[[1L]] <- drop.tip(treeAll,spnames[-match(ModelSpecies,spnames[,1L]),2L])
#
#
# Aggregate & Format the Concentration.
Chemical.Table <- data.frame(ShortNames=treeChemical[[1L]]$tip.label,Maj_group=NA,Min_group=NA)
for(i in Chemical.properties[,"Chemical"])
  Chemical.Table[,i] <- NA
for(i in as.character(Chemical.Table[,"ShortNames"])) {
  # i <- as.character(Chemical.Table[,"ShortNames"])[7L]
  Chemical.Table[Chemical.Table[,"ShortNames"]==i,"Maj_group"] <- as.character(Chemical.data[Chemical.data[,"ShortNames"] == i,"Maj_group"][1L])
  Chemical.Table[Chemical.Table[,"ShortNames"]==i,"Min_group"] <- as.character(Chemical.data[Chemical.data[,"ShortNames"] == i,"Min_group"][1L])
  tmp <- Chemical.data[Chemical.data[,"ShortNames"] == i,]
  tmp[,"Conc_mmgL"] <- log10(tmp[,"Conc_mmgL"])
  tmp <- summaryBy(Conc_mmgL~Taxon+Chemical,keep.names=T,data=tmp)
  tmp <- summaryBy(Conc_mmgL~Chemical,data=tmp)
  Chemical.Table[Chemical.Table[,"ShortNames"]==i,match(as.character(tmp[,1L]),Chemical.properties[,"Chemical"])+3] <- tmp[,2L]
} ; rm(i,tmp)
#
# Make a list of Cases: Case 1 with 4 HM and Case2 with 6 HM.
treeChemical[[2L]] <- drop.tip(treeChemical[[1L]],as.character(Chemical.Table[!complete.cases(Chemical.Table),"ShortNames"]))
HMetal <- list(FourHM=c("Cd","Cu","Hg","Zn"),SixHM=c("Cd","Cu","Hg","Ni","Pb","Zn"))
names(treeChemical) <- names(HMetal)
#
#
# Figure Phylogenetic Tree w/ bootstraps.
### Plot the tree with bootstrap values at nodes (all nodes with >50% support were collapsed):
x11()
tmp <- treeChemical[[1L]] 
tmp$node.label<-as.numeric(tmp$node.label)
tmp$node.label<-round((tmp$node.label)/10,0); tmp$node.label[tmp$node.label=="100"] <- "*"
tmp$node.label <- paste(" ",tmp$node.label,sep="") ; tmp$tip.label <- paste("  ",tmp$tip.label,sep="")
plot(tmp,cex=1, no.margin = TRUE)
nodelabels(tmp$node.label, adj = c(1.2, -0.2),  cex = 0.8, frame = "n",font=2)
#write.tree(treeChemical[[1L]],file="treeChemical.treefile")
#
## Building model:
#
# 1. Get info from the Phylo tree as fig (package:ape) 2 directed graph. for both models
graphChemical <- list()
for(m in names(HMetal))
  graphChemical[[m]] <- Phylo2DirectedGraph(treeChemical[[m]])
rm(m)
#
# 2. Phylogenetic Eigenvector Maps from the directed graph.
BChemical <- list()
for(m in names(HMetal))
  BChemical[[m]] <- PEMInfluence(graphChemical[[m]])[graphChemical[[m]]$vertex$species,]
rm(m)
#
# 3. Get concentrations for 4 metals and 6 metal.
LC50Chemical <- list()
for(m in names(HMetal)) {
  LC50Chemical[[m]] <- as.matrix(Chemical.Table[,match(HMetal[[m]],colnames(Chemical.Table))])
  LC50Chemical[[m]] <- LC50Chemical[[m]][complete.cases(Chemical.Table[,HMetal[[m]]]),]
  rownames(LC50Chemical[[m]]) <- as.character(Chemical.Table[complete.cases(Chemical.Table[,HMetal[[m]]]),"ShortNames"])
} ; rm(m)
#
#
# 4. Calculate and manipulate PEM
PEMChemical <- list()
for(m in names(HMetal))
  PEMChemical[[m]] <- PEM.fitSimple(y=LC50Chemical[[m]],x=NULL,d="distance",sp="species",lower=0,upper=1,w=graphChemical[[m]])
rm(m)
#
for(m in names(HMetal)) {
  print(PEMChemical[[m]]$optim$par)
} ; rm(m)
#
## 5. Assuming neutral trait evolution
# For the paper Guenard et al.(2011)
if(FALSE) {
PEMChemical <- list()
for(m in names(HMetal))
  PEMChemical[[m]] <- PEM.forcedSimple(y=LC50Chemical[[m]],x=NULL,d="distance",sp="species",w=graphChemical[[m]],a=0,psi=1)
rm(m)
}
# 
### Should make a list of trials. Can add others: FirstHydro OxyRed LogFstab MOHSol XSQmr ZSQoverR AnoverIP SigmaSQ
descr <- c("FirstHydro","OxyRed","LogFstab","MOHSol","XSQmr","ZSQoverR","AnoverIP","SigmaSQ")
models <- list(FourHM=list(),SixHM=list())
for(i in 1L:8) {
  models[["FourHM"]][[descr[i]]] <- list(descr[i])
  models[["SixHM"]][[descr[i]]] <- list(descr[i])
}
for(i in 1L:7) {
  for(j in (i+1L):8L) {
    models[["SixHM"]][[paste(descr[c(i,j)],collapse="+")]] <- descr[c(i,j)]
  }
}
rm(descr,i,j)
#
# Generates all possible model combinations.
# Figure adapted from Guenard et al.(2011)
#
ResList <- list()
for(d in names(models)) {
  # d <- names(models)[1L]
  for(dd in names(models[[d]])) {
    # dd <- names(models[[d]])[1L]
    trial <- paste(d,dd,sep="-")
    ResList[[trial]] <- list()
    ResList[[trial]]$data <- d   
    ResList[[trial]]$descr <- unlist(models[[d]][[dd]])
    ResList[[trial]]$Z <- as.matrix(cbind(MainSpecies=1,Chemical.properties[match(HMetal[[ResList[[trial]]$data]],Chemical.properties[,"Chemical"]),ResList[[trial]]$descr,drop=FALSE]))
    rownames(ResList[[trial]]$Z) <- colnames(LC50Chemical[[ResList[[trial]]$data]])
    ResList[[trial]]$U <- as.matrix(cbind(MainPropert=1,PEMChemical[[ResList[[trial]]$data]]$u))
    ResList[[trial]]$ZkronU <- ResList[[trial]]$Z %x% ResList[[trial]]$U
    dimnames(ResList[[trial]]$ZkronU) <- list(paste(rep(rownames(ResList[[trial]]$U),nrow(ResList[[trial]]$Z)),rep(rownames(ResList[[trial]]$Z),each=nrow(ResList[[trial]]$U)),sep="_"),
                                              paste(rep(colnames(ResList[[trial]]$U),ncol(ResList[[trial]]$Z)),rep(colnames(ResList[[trial]]$Z),each=ncol(ResList[[trial]]$U)),sep="_"))
    ResList[[trial]]$Y <- LC50Chemical[[ResList[[trial]]$data]]
    dim(ResList[[trial]]$Y) <- c(prod(dim(LC50Chemical[[ResList[[trial]]$data]])),1L)
    dimnames(ResList[[trial]]$Y) <- list(paste(rep(rownames(LC50Chemical[[ResList[[trial]]$data]]),ncol(LC50Chemical[[ResList[[trial]]$data]])),
                                               rep(colnames(LC50Chemical[[ResList[[trial]]$data]]),each=nrow(LC50Chemical[[ResList[[trial]]$data]])),sep="X"),"y")
    ### Building bilinear models:
    ResList[[trial]]$BLMfit <- list()
    ResList[[trial]]$BLMfit$lm1 <- sequentialAICcblm(y=ResList[[trial]]$Y,ZU=ResList[[trial]]$ZkronU[,-1],k=2)
    ResList[[trial]]$BLMfit$summary <- summary(ResList[[trial]]$BLMfit$lm1)
    ResList[[trial]]$BLMfit$included <- coef(ResList[[trial]]$BLMfit$lm1)
    ResList[[trial]]$BLMfit$Cmat <- numeric(0)
    ResList[[trial]]$BLMfit$Cmat[c("(Intercept)",colnames(ResList[[trial]]$ZkronU)[-1])] <- 0
    ResList[[trial]]$BLMfit$Cmat[labels(ResList[[trial]]$BLMfit$included)] <- ResList[[trial]]$BLMfit$included
    dim(ResList[[trial]]$BLMfit$Cmat) <- c(ncol(ResList[[trial]]$U),ncol(ResList[[trial]]$Z))
    dimnames(ResList[[trial]]$BLMfit$Cmat) <- list(colnames(ResList[[trial]]$U),colnames(ResList[[trial]]$Z))
    ResList[[trial]]$BLMfit$CmatStrip <- ResList[[trial]]$BLMfit$Cmat[which(rowSums(ResList[[trial]]$BLMfit$Cmat!=0)!=0),]
    ResList[[trial]]$BLMfit$Yhat <- as.matrix(ResList[[trial]]$U%*%ResList[[trial]]$BLMfit$Cmat%*%t(ResList[[trial]]$Z))
    # Plotting:
    #quartz(height=5.0,width=7.25) # To have proportional squares run: 2L (width=9, height=2.6), 1L (width=7.25, height=5)
    tiff(file=paste("Fig2-",trial,".tiff",sep=""), width=9, height=2.6, 
         units = "in", compression="lzw", res = 600, type="cairo")
    layout(matrix(c(1,1,seq(2,length.out=ncol(LC50Chemical[[ResList[[trial]]$data]]))),1L,2L+ncol(LC50Chemical[[ResList[[trial]]$data]])))
    par(mar=c(4.25,2.1,2.2,0.1))
    tmp <- treeChemical[[ResList[[trial]]$data]]
    tmp$node.label<-as.numeric(tmp$node.label)
    tmp$node.label<-round((tmp$node.label)/10,0); tmp$node.label[tmp$node.label=="100"] <- "*"
    tmp$node.label <- paste(" ",tmp$node.label,sep="") ; tmp$tip.label <- paste("  ",tmp$tip.label,sep="")
    plot(tmp, cex=1.15)
    nodelabels(tmp$node.label, adj = c(1.1, -0.2),  cex = 0.9, frame = "n",font=2)
    #
    par(mar=c(4.1,0.6,2.1,0.8))
    for (i in HMetal[[ResList[[trial]]$data]]) {
      plot(NA, xlim=c(-0.5,6.25), xlab="", type = "n", main=i, yaxt="n", ylim=c(1,nrow(LC50Chemical[[ResList[[trial]]$data]])), xaxt="n")
      axis(1,at=c(0,2,4,6),labels=c("1",expression(10^2),expression(10^4),expression(10^6)))
      abline(h=1L:nrow(LC50Chemical[[ResList[[trial]]$data]]),col=gray(0.67))
      abline(v=-1:6,col=gray(0.67))
      points(y=1L:nrow(LC50Chemical[[ResList[[trial]]$data]]), x=LC50Chemical[[ResList[[trial]]$data]][,i], cex=2, pch=15)
      points(y=1L:nrow(LC50Chemical[[ResList[[trial]]$data]]), x=ResList[[trial]]$BLMfit$Yhat[,i] ,pch=1, cex=2)
    }
    mtext(expression(paste("                                          ",LC[50]," (",mu*g%.%L^{-1},")",sep="")),side=1,line=-1.5,outer=TRUE,cex=0.9,at=1/2)
    dev.off()
    n <- nrow(LC50Chemical[[ResList[[trial]]$data]]) ; m <- ncol(LC50Chemical[[ResList[[trial]]$data]]) ; p <- ncol(ResList[[trial]]$U) ; q <- ncol(ResList[[trial]]$Z)
    ResList[[trial]]$Ftests <- list()
    ResList[[trial]]$Ftests$Overall <- c(value=ResList[[trial]]$BLMfit$summary$fstatistic[1L],
                                         p=pf(ResList[[trial]]$BLMfit$summary$fstatistic[1L],
                                              ResList[[trial]]$BLMfit$summary$fstatistic[2L],
                                              ResList[[trial]]$BLMfit$summary$fstatistic[2L],lower.tail=FALSE))
    mask <- ResList[[trial]]$BLMfit$CmatStrip ; mask[] <- 1
    ResList[[trial]]$Ftests$RowTests <- matrix(NA,0L,4L)
    for (i in 2L:nrow(ResList[[trial]]$BLMfit$CmatStrip)) {
      maski <- mask ; maski[-i,] <- 0 ; maski[,1L] <- 0
      Yfiti <- ResList[[trial]]$U[,rownames(ResList[[trial]]$BLMfit$CmatStrip),drop=FALSE]%*%(ResList[[trial]]$BLMfit$CmatStrip*maski)%*%
                 t(ResList[[trial]]$Z[,colnames(ResList[[trial]]$BLMfit$CmatStrip),drop=FALSE])
      Fstat <- (n*m-sum(mask))*sum(Yfiti^2)/(sum(ResList[[trial]]$BLMfit$CmatStrip*maski!=0)*sum((LC50Chemical[[ResList[[trial]]$data]]-ResList[[trial]]$BLMfit$Yhat)^2))
      ResList[[trial]]$Ftests$RowTests <- rbind(ResList[[trial]]$Ftests$RowTests,
                                                round(c(Fstat,sum(ResList[[trial]]$BLMfit$CmatStrip*maski!=0),
                                                        n*m-sum(mask),pf(Fstat,sum(ResList[[trial]]$BLMfit$CmatStrip*maski!=0),n*m-sum(mask),lower.tail=FALSE)),3))
    }
    ResList[[trial]]$Ftests$ColTests <- matrix(NA,0L,4L)
    for (i in 2L:ncol(ResList[[trial]]$BLMfit$CmatStrip)) {
      maski <- mask ; maski[1L,] <- 0 ; maski[,-i] <- 0
      Yfiti <- ResList[[trial]]$U[,rownames(ResList[[trial]]$BLMfit$CmatStrip),drop=FALSE]%*%(ResList[[trial]]$BLMfit$CmatStrip*maski)%*%
                 t(ResList[[trial]]$Z[,colnames(ResList[[trial]]$BLMfit$CmatStrip),drop=FALSE])
      Fstat <- (n*m-sum(mask))*sum(Yfiti^2)/(sum(ResList[[trial]]$BLMfit$CmatStrip*maski!=0)*sum((LC50Chemical[[ResList[[trial]]$data]]-ResList[[trial]]$BLMfit$Yhat)^2))
      ResList[[trial]]$Ftests$ColTests <- rbind(ResList[[trial]]$Ftests$ColTests,
                                                round(c(Fstat,sum(ResList[[trial]]$BLMfit$CmatStrip*maski!=0),
                                                        n*m-sum(mask),pf(Fstat,sum(ResList[[trial]]$BLMfit$CmatStrip*maski!=0),n*m-sum(mask),lower.tail=FALSE)),3))
    }
  }
} ; rm(d,dd,i,n,m,p,q,mask,maski,trial,Fstat,Yfiti)
#
# All adj.r.squared
for(i in names(ResList))
  cat(i," : ",ResList[[i]]$BLMfit$summary$adj.r.squared,"\n")
#
round(cor(Chemical.properties[,-1L]),3)
#
# Results: 
#
# R2-adj - see text for rational
# FourHM-SigmaSQ:             0.66
# SixHM-SigmaSQ:              0.78 
#
# Marginal difference in 2 properties - 6HM discarded this option
#
# Export the results
#
Table1 <- data.frame(row.names=names(ResList),nsp=integer(length(ResList)),nhm=integer(length(ResList)),Chem1=character(length(ResList)),
                     Chem2=character(length(ResList)),MainPhylo=integer(length(ResList)),MainChem=integer(length(ResList)),
                     Inter=integer(length(ResList)),AdjRsq=numeric(length(ResList)),Focus=character(length(ResList)),stringsAsFactors=FALSE)
for(i in 1L:length(ResList)) {
  Table1[i,"nsp"] <- nrow(ResList[[i]]$U)
  Table1[i,"nhm"] <- nrow(ResList[[i]]$Z)
  Table1[i,"Chem1"] <- colnames(ResList[[i]]$Z)[2L]
  Table1[i,"Chem2"] <- if(!is.na(colnames(ResList[[i]]$Z)[3L])) colnames(ResList[[i]]$Z)[3L] else ""
  Table1[i,"MainPhylo"] <- sum(ResList[[i]]$BLMfit$Cmat[-1L,1L]!=0)
  Table1[i,"MainChem"] <- sum(ResList[[i]]$BLMfit$Cmat[1L,-1L]!=0)
  Table1[i,"Inter"] <- sum(ResList[[i]]$BLMfit$Cmat[-1L,-1L]!=0)
  Table1[i,"AdjRsq"] <- round(ResList[[i]]$BLMfit$summary$adj.r.squared,2)
} ; rm(i)
Table1[Table1[,"nhm"]==4L&Table1[,"AdjRsq"]==max(Table1[Table1[,"nhm"]==4L,"AdjRsq"]),"Focus"] <- "4HM-1Z"
Table1[Table1[,"nhm"]==6L&Table1[,"Chem2"]==""&Table1[,"AdjRsq"]==max(Table1[Table1[,"nhm"]==6L&Table1[,"Chem2"]=="","AdjRsq"]),"Focus"] <- "6HM-1Z"
Table1[Table1[,"nhm"]==6L&Table1[,"Chem2"]!=""&Table1[,"AdjRsq"]==max(Table1[Table1[,"nhm"]==6L&Table1[,"Chem2"]!="","AdjRsq"]),"Focus"] <- "6HM-2Z"
# write.csv(Table1, file='All_PhysicoChem.csv')  
#
# Chosen models include sigmasq (see paper for rational). All the calculations below are done 
# with sigmasq
#
# Add the sigmasq
Table1$Focus[Table1$Chem1=='SigmaSQ'&Table1$nhm==4]<-'4HM-1Z-1'
#
Table2 <- list()
for(f in Table1[Table1[,"Focus"]!="","Focus"]) {
  # f <- "4HM-1Z"
  m <- rownames(Table1)[Table1[,"Focus"]==f]
  tmp <- ResList[[m]]$BLMfit$CmatStrip ; tmp[tmp==0] <- NA
  tmp <- round(tmp,3) ; tmp[is.na(tmp)] <- "---"
  Table2[[m]] <- rbind(cbind(tmp,F="",p=""),F="",p="")
  tmp <- as.character(ResList[[m]]$Ftests$RowTests[,1L])
  tmp[tmp=="NaN"] <- ""
  Table2[[m]][1L:length(tmp),"F"] <- tmp
  tmp <- as.character(ResList[[m]]$Ftests$RowTests[,4L])
  tmp[as.numeric(tmp)<0.001] <- "<0.001" ; tmp[tmp=="NaN"] <- ""
  Table2[[m]][1L:length(tmp),"p"] <- tmp
  tmp <- as.character(ResList[[m]]$Ftests$ColTests[,1L])
  tmp[tmp=="NaN"] <- ""
  Table2[[m]]["F",1L:length(tmp)] <- tmp
  tmp <- as.character(ResList[[m]]$Ftests$ColTests[,4L])
  tmp[as.numeric(tmp)<0.001] <- "<0.001" ; tmp[tmp=="NaN"] <- ""
  Table2[[m]]["p",1L:length(tmp)] <- tmp
  Table2[[m]]["F","F"] <- round(ResList[[m]]$Ftests$Overall[1L],3)
  Table2[[m]]["p","p"] <- if(ResList[[m]]$Ftests$Overall[2L]>=0.001) round(ResList[[m]]$Ftests$Overall[2L],3) else "<0.001"
} ; rm(f,m,tmp)
for(i in 1L:length(Table2)) {
  # i <- 2L# relevant 2L - 4HM and 3L - 6HM
  write.csv(Table2[[i]],file=paste("Table 2",letters[i],".csv",sep=""))
} ; rm(i)
#
# 
### Foci: a list of selected scenarios. See Paper for details on selecting "SigmaSQ"
Foci <- list(Focus1=list(data="FourHM",descr="SigmaSQ"),
             Focus2=list(data="SixHM", descr="SigmaSQ"))
#
### Performing the leave-one-out jacknifes: Type 1A 
LC50_jknife <- list()
for(j in names(Foci)) {
  # j <- names(Foci)[1L]
  LC50_jknife[[j]]$Ypred <- LC50Chemical[[Foci[[j]]$data]] ; LC50_jknife[[j]]$Ypred[] <- NA
  LC50_jknife[[j]]$Z <- as.matrix(cbind(MainSpecies=1,Chemical.properties[match(HMetal[[Foci[[j]]$data]],
                                        Chemical.properties[,"Chemical"]),Foci[[j]]$descr,drop=FALSE]))
  rownames(LC50_jknife[[j]]$Z) <- colnames(LC50Chemical[[Foci[[j]]$data]])
  LC50_jknife[[j]]$JN <- list()
  for(i in 1L:nrow(LC50Chemical[[Foci[[j]]$data]])) {
    # i <- 1L
    LC50_jknife[[j]]$JN[[i]] <- list()
    LC50_jknife[[j]]$JN[[i]]$LC50 <- LC50Chemical[[Foci[[j]]$data]][-i,,drop=FALSE]
    LC50_jknife[[j]]$JN[[i]]$locs <- getGraphLocations(treeChemical[[Foci[[j]]$data]],rownames(LC50Chemical[[Foci[[j]]$data]])[i])
    LC50_jknife[[j]]$JN[[i]]$PEMJN <- PEM.forcedSimple(y=LC50_jknife[[j]]$JN[[i]]$LC50,x=NULL,w=LC50_jknife[[j]]$JN[[i]]$locs$x,d="distance",sp="species",a=0,psi=1)
    # LC50_jknife[[j]]$JN[[i]]$PEMJN <- PEM.fitSimple(y=LC50_jknife[[j]]$JN[[i]]$LC50,x=NULL,w=LC50_jknife[[j]]$JN[[i]]$locs$x,d="distance",sp="species",lower=0,upper=1)
    LC50_jknife[[j]]$JN[[i]]$U <- as.matrix(cbind(MainProp=1,LC50_jknife[[j]]$JN[[i]]$PEMJN$u))
    # cbind(attr(LC50_jknife[[j]]$JN[[i]]$locs$x,"vlabel")[LC50_jknife[[j]]$JN[[i]]$locs$x$vertex$species],rownames(LC50_jknife[[j]]$JN[[i]]$LC50))
    # cbind(rownames(LC50_jknife[[j]]$JN[[i]]$PEMJN$u),rownames(LC50_jknife[[j]]$JN[[i]]$LC50))
    n <- nrow(LC50_jknife[[j]]$JN[[i]]$LC50) ; m <- ncol(LC50_jknife[[j]]$JN[[i]]$LC50) ; p <- ncol(LC50_jknife[[j]]$JN[[i]]$U) ; q <- ncol(LC50_jknife[[j]]$Z)
    LC50_jknife[[j]]$JN[[i]]$ZU <- LC50_jknife[[j]]$Z%x%LC50_jknife[[j]]$JN[[i]]$U
    dimnames(LC50_jknife[[j]]$JN[[i]]$ZU) <- list(paste(rep(rownames(LC50_jknife[[j]]$JN[[i]]$U),nrow(LC50_jknife[[j]]$Z)),rep(rownames(LC50_jknife[[j]]$Z),each=nrow(LC50_jknife[[j]]$JN[[i]]$U)),sep="X"),
                                                  paste(rep(colnames(LC50_jknife[[j]]$JN[[i]]$U),ncol(LC50_jknife[[j]]$Z)),rep(colnames(LC50_jknife[[j]]$Z),each=ncol(LC50_jknife[[j]]$JN[[i]]$U)),sep="X"))    
    LC50_jknife[[j]]$JN[[i]]$y <- LC50_jknife[[j]]$JN[[i]]$LC50 ; dim(LC50_jknife[[j]]$JN[[i]]$y) <- c(prod(dim(LC50_jknife[[j]]$JN[[i]]$LC50)),1L)
    dimnames(LC50_jknife[[j]]$JN[[i]]$y) <- list(paste(rep(rownames(LC50_jknife[[j]]$JN[[i]]$LC50),ncol(LC50_jknife[[j]]$JN[[i]]$LC50)),
                                                       rep(colnames(LC50_jknife[[j]]$JN[[i]]$LC50),each=nrow(LC50_jknife[[j]]$JN[[i]]$LC50)),sep="X"),NULL)
    LC50_jknife[[j]]$JN[[i]]$lm1 <- sequentialAICcblm(y=LC50_jknife[[j]]$JN[[i]]$y,ZU=LC50_jknife[[j]]$JN[[i]]$ZU[,-1],k=2)
    LC50_jknife[[j]]$JN[[i]]$Cmat <- numeric(0) ; LC50_jknife[[j]]$JN[[i]]$Cmat[c("(Intercept)",colnames(LC50_jknife[[j]]$JN[[i]]$ZU)[-1])] <- 0
    LC50_jknife[[j]]$JN[[i]]$included <- coef(LC50_jknife[[j]]$JN[[i]]$lm1) ; LC50_jknife[[j]]$JN[[i]]$Cmat[labels(LC50_jknife[[j]]$JN[[i]]$included)] <- LC50_jknife[[j]]$JN[[i]]$included
    dim(LC50_jknife[[j]]$JN[[i]]$Cmat) <- c(ncol(LC50_jknife[[j]]$JN[[i]]$U),ncol(LC50_jknife[[j]]$Z))
    dimnames(LC50_jknife[[j]]$JN[[i]]$Cmat) <- list(colnames(LC50_jknife[[j]]$JN[[i]]$U),colnames(LC50_jknife[[j]]$Z))
    LC50_jknife[[j]]$JN[[i]]$CmatStrip <- LC50_jknife[[j]]$JN[[i]]$Cmat[which(rowSums(LC50_jknife[[j]]$JN[[i]]$Cmat!=0)!=0),]
    LC50_jknife[[j]]$JN[[i]]$scores <- Locations2PEMscores(LC50_jknife[[j]]$JN[[i]]$PEMJN,LC50_jknife[[j]]$JN[[i]]$locs)
    LC50_jknife[[j]]$JN[[i]]$scores$Uscores <- cbind(1,LC50_jknife[[j]]$JN[[i]]$scores$scores) ; colnames(LC50_jknife[[j]]$JN[[i]]$scores$Uscores) <- colnames(LC50_jknife[[j]]$JN[[i]]$U)
    # colnames(LC50_jknife[[j]]$JN[[i]]$CmatStrip)
    # as.matrix(LC50_jknife[[j]]$JN[[i]]$scores$Uscores[,rownames(LC50_jknife[[j]]$JN[[i]]$CmatStrip),drop=FALSE])%*%LC50_jknife[[j]]$JN[[i]]$CmatStrip%*%t(LC50_jknife[[j]]$Z[,colnames(LC50_jknife[[j]]$JN[[i]]$CmatStrip),drop=FALSE])
    LC50_jknife[[j]]$Ypred[i,] <- as.matrix(LC50_jknife[[j]]$JN[[i]]$scores$Uscores[,rownames(LC50_jknife[[j]]$JN[[i]]$CmatStrip),drop=FALSE])%*%
                                 LC50_jknife[[j]]$JN[[i]]$CmatStrip%*%t(LC50_jknife[[j]]$Z[,colnames(LC50_jknife[[j]]$JN[[i]]$CmatStrip),drop=FALSE])
  }
} ; rm(i,j,n,m,p,q)
#
#
for(i in names(Foci))
  LC50_jknife[[i]]$Psquare <- Psquare(as.numeric(LC50Chemical[[Foci[[i]]$data]]),as.numeric(LC50_jknife[[i]]$Ypred))
#
for(i in names(Foci))
  LC50_jknife[[i]]$lm<- lm(as.numeric(LC50_jknife[[i]]$Ypred)~as.numeric(LC50Chemical[[Foci[[i]]$data]])); rm(i)
#
for(i in names(Foci))
  cat(i,":",LC50_jknife[[i]]$Psquare,"\n")
#
# Results:
# P2 (same as q2 in the paper) for the models:
# FourHM-SigmaSQ:             0.51
# SixHM-SigmaSQ:              0.68
#
#
### Deviation factors (|d|)
devChemical <- list() ; ddisplChemical <- list()
for(i in names(Foci)) {
  # i <- names(Foci)[1L]
  devChemical[[i]] <- matrix(0,nrow(LC50_jknife[[i]]$Ypred),ncol(LC50_jknife[[i]]$Ypred),dimnames=dimnames(LC50_jknife[[i]]$Ypred))
  devChemical[[i]][LC50_jknife[[i]]$Ypred > LC50Chemical[[Foci[[i]]$data]]] <-
    ((10^(LC50_jknife[[i]]$Ypred-LC50Chemical[[Foci[[i]]$data]]))-1)[LC50_jknife[[i]]$Ypred > LC50Chemical[[Foci[[i]]$data]]]
  devChemical[[i]][LC50_jknife[[i]]$Ypred < LC50Chemical[[Foci[[i]]$data]]] <-
    -((10^(LC50Chemical[[Foci[[i]]$data]]-LC50_jknife[[i]]$Ypred))-1)[LC50_jknife[[i]]$Ypred < LC50Chemical[[Foci[[i]]$data]]]
  ddisplChemical[[i]] <- matrix(0,nrow(LC50_jknife[[i]]$Ypred),ncol(LC50_jknife[[i]]$Ypred),dimnames=dimnames(LC50_jknife[[i]]$Ypred))
  ddisplChemical[[i]][devChemical[[i]] > 0] <- log10(devChemical[[i]][devChemical[[i]] > 0]+1)
  ddisplChemical[[i]][devChemical[[i]] < 0] <- -log10(1-devChemical[[i]][devChemical[[i]] < 0])
}
#
#
# Deviation Factors Plotted
#
##### 4HM
dev_stat1<-data.frame(devChemical[[1]])
dev_stat1$Species<- rownames(dev_stat1)
#
# Deviating species
dev_stat1[abs(dev_stat1$Cd)>5,-(2:4)] # Cd
dev_stat1[abs(dev_stat1$Cu)>5,-c(1,3:4)]# Cu
dev_stat1[abs(dev_stat1$Hg)>5,-c(1:2,4)]# Hg
dev_stat1[abs(dev_stat1$Zn)>5,-c(1:3)]# Zn
#
# Overall - 73 %
1-(nrow(dev_stat1[abs(dev_stat1$Cd)>=5,])+
  nrow(dev_stat1[abs(dev_stat1$Cu)>=5,])+
  nrow(dev_stat1[abs(dev_stat1$Hg)>=5,])+
  nrow(dev_stat1[abs(dev_stat1$Zn)>=5,]))/(nrow(dev_stat1)*4)
#
# Individual - Change each heavymetal: 71% Cd; 65% Cu; 84% Hg; 74% Zn   
1-(nrow(dev_stat1[abs(dev_stat1$Cd)>=5,])/31)
#
tiff(file='Deviation Factor_4HM.tiff',width = 9, height = 7, units = "in", res = 250)  
#
layout(matrix(c(1,1,1,2,3,4,5),1,7))
par(mar=c(4.25,2.1,2.2,0.1))
plot(treeChemical[[1L]], cex=1.15)
par(mar=c(4.1,0.6,2.1,0.8))
for (i in HMetal[[1]]) {
  plot(NA, xlim=c(-75,75), xlab="", type = "n", main=i, yaxt="n", ylim=c(1,31), xaxt="n")
  axis(1,at=c(-75,-50,-25,0,25,50,75),labels=c(expression(-75),expression(-50),expression(-25),"0",expression(25),expression(50),expression(75)))
  #abline(h=1:33,col=gray(0.67),v=c(-50,-25,0,25,50))
  abline(v=c(-5,5),lty = 3)
  abline(v=0)
  points(y=1:31, x=devChemical[[1]][,i], cex=1.0, pch=17)
}
mtext(expression(paste("                                                                 Deviation factor",sep="")),side=1,line=-1.5,outer=TRUE,cex=0.9,at=1/2)
dev.off()
#
# 6HM
dev_stat2<-data.frame(devChemical[[2]])
dev_stat2$Species<- rownames(dev_stat2)
#
# Deviating species
dev_stat2[abs(dev_stat2$Cd)>5,-(2:6)] 
dev_stat2[abs(dev_stat2$Cu)>5,-c(1,3:6)]
dev_stat2[abs(dev_stat2$Hg)>5,-c(1:2,4:6)]
dev_stat2[abs(dev_stat2$Ni)>5,-c(1:3,5:6)]
dev_stat2[abs(dev_stat2$Pb)>5,-c(1:4,6)]
dev_stat2[abs(dev_stat2$Zn)>5,-c(1:5)]
#
# Overall - 74 % deviating
1-(nrow(dev_stat2[abs(dev_stat2$Cd)>=5,])+
   nrow(dev_stat2[abs(dev_stat2$Cu)>=5,])+
   nrow(dev_stat2[abs(dev_stat2$Hg)>=5,])+
   nrow(dev_stat2[abs(dev_stat2$Zn)>=5,])+
   nrow(dev_stat2[abs(dev_stat2$Ni)>=5,])+
   nrow(dev_stat2[abs(dev_stat2$Pb)>=5,]))/(nrow(dev_stat2)*6)
#
# Individual - change heavy metal. Cd: 79%, Cu:64%, Hg: 79%, Zn: 79%, Ni: 64%, Pb:79%
1-(nrow(dev_stat2[abs(dev_stat2$Cd)>=5,])/nrow(dev_stat2))
#
tiff(file='Deviation Factor_6HM.tiff', 
     width = 9, height = 3, units = "in", res = 250)    
#
layout(matrix(c(1,1,1,2,3,4,5,6,7),1,9))
par(mar=c(4.25,2.1,2.2,0.1))
plot(treeChemical[[2L]], cex=1.15, main="")
par(mar=c(4.1,0.6,2.1,0.8))
for (i in HMetal[[2]]) {
  plot(NA, xlim=c(-33,20), xlab="", type = "n", main=i, yaxt="n", ylim=c(1,14), xaxt="n")
  axis(1,at=c(-30,-15,0,15),labels=c(expression(-30),expression(-15),"0",expression(15)))
  #abline(h=1:33,col=gray(0.67),v=c(-50,-25,0,25,50))
  abline(v=c(-5,5),lty = 3)
  abline(v=0)
  points(y=1:14, x=devChemical[[2]][,i], cex=1.0, pch=17)
}
mtext(expression(paste("                                                                   Deviation factor",sep="")),side=1,line=-1.5,outer=TRUE,cex=0.9,at=1/2)
dev.off()
#
#
### Additive variance proportions for Foci (for the main models)
#
RsqList <- list()
for(trial in names(ResList)) {
  # trial <- names(ResList)[1L]
  RsqList[[trial]] <- list()
  RsqList[[trial]]$data <- ResList[[trial]]$data
  RsqList[[trial]]$descr <- ResList[[trial]]$descr
  #
  RsqList[[trial]]$Z <- as.matrix(cbind(MainSpecies=1,Chemical.properties[match(HMetal[[RsqList[[trial]]$data]],Chemical.properties[,"Chemical"]),RsqList[[trial]]$descr,drop=FALSE]))
  rownames(RsqList[[trial]]$Z) <- colnames(LC50Chemical[[RsqList[[trial]]$data]])
  RsqList[[trial]]$U <- as.matrix(cbind(MainPropert=1,PEMChemical[[RsqList[[trial]]$data]]$u))
  RsqList[[trial]]$ZkronU <- RsqList[[trial]]$Z %x% RsqList[[trial]]$U
  dimnames(RsqList[[trial]]$ZkronU) <- list(paste(rep(rownames(RsqList[[trial]]$U),nrow(RsqList[[trial]]$Z)),rep(rownames(RsqList[[trial]]$Z),each=nrow(RsqList[[trial]]$U)),sep="_"),
                                            paste(rep(colnames(RsqList[[trial]]$U),ncol(RsqList[[trial]]$Z)),rep(colnames(RsqList[[trial]]$Z),each=ncol(RsqList[[trial]]$U)),sep="_"))
  RsqList[[trial]]$ZkronU1 <- RsqList[[trial]]$Z[,1L,drop=FALSE] %x% RsqList[[trial]]$U
  dimnames(RsqList[[trial]]$ZkronU1) <- list(paste(rep(rownames(RsqList[[trial]]$U),nrow(RsqList[[trial]]$Z[,1L,drop=FALSE])),rep(rownames(RsqList[[trial]]$Z[,1L,drop=FALSE]),
                                                                                                                                  each=nrow(RsqList[[trial]]$U)),sep="_"),
                                             paste(rep(colnames(RsqList[[trial]]$U),ncol(RsqList[[trial]]$Z[,1L,drop=FALSE])),rep(colnames(RsqList[[trial]]$Z[,1L,drop=FALSE]),
                                                                                                                                  each=ncol(RsqList[[trial]]$U)),sep="_"))
  RsqList[[trial]]$ZkronU2 <- RsqList[[trial]]$Z %x% RsqList[[trial]]$U[,1L,drop=FALSE]
  dimnames(RsqList[[trial]]$ZkronU2) <- list(paste(rep(rownames(RsqList[[trial]]$U[,1L,drop=FALSE]),nrow(RsqList[[trial]]$Z)),rep(rownames(RsqList[[trial]]$Z),
                                                                                                                                  each=nrow(RsqList[[trial]]$U[,1L,drop=FALSE])),sep="_"),
                                             paste(rep(colnames(RsqList[[trial]]$U[,1L,drop=FALSE]),ncol(RsqList[[trial]]$Z)),rep(colnames(RsqList[[trial]]$Z),
                                                                                                                                  each=ncol(RsqList[[trial]]$U[,1L,drop=FALSE])),sep="_"))
  RsqList[[trial]]$Y <- ResList[[trial]]$Y
  RsqList[[trial]]$BLMfit <- list()
  RsqList[[trial]]$BLMfit$summary <- summary(sequentialAICcblm(y=RsqList[[trial]]$Y,ZU=RsqList[[trial]]$ZkronU[,-1],k=2))
  RsqList[[trial]]$BLMfit$summary1 <- summary(sequentialAICcblm(y=RsqList[[trial]]$Y,ZU=RsqList[[trial]]$ZkronU1[,-1],k=2))
  RsqList[[trial]]$BLMfit$summary2 <- summary(lm(RsqList[[trial]]$Y~RsqList[[trial]]$ZkronU2[,-1]))
} ; rm(trial)
#
for(trial in c("FourHM-SigmaSQ","SixHM-SigmaSQ")) {
  cat(trial,": Global :",RsqList[[trial]]$BLMfit$summary$adj.r.squared,"\n")
  cat(trial,": Phylogeny alone:",RsqList[[trial]]$BLMfit$summary1$adj.r.squared,"\n")
  cat(trial,": Properties alone:",RsqList[[trial]]$BLMfit$summary2$adj.r.squared,"\n\n\n")
} ; rm(trial)
#
#FourHM-SigmaSQ : Global : 0.66 
#FourHM-SigmaSQ : Phylogeny alone: 0.44 
#FourHM-SigmaSQ : Properties alone: 0.19
#
#SixHM-SigmaSQ : Global : 0.78 
#SixHM-SigmaSQ : Phylogeny alone: 0.43
#SixHM-SigmaSQ : Properties alone: 0.29
#
### Conclusion: In all case, phylogeny explains more than properties alone but using both with interactions gives the best models.
#
####### External-validation
#
# Type 1B - new species
# Predict the toxicity of Cd,Cu,Zn & Hg for species removed from the 6HM-Model
#
# Procedure is the same as for the leave_one_out cross validation. We don't leave 1 out but 
# the identified 14 species for which we need predictions for. Initially these species are removed from the
# complete 33-species model. A tree with the remaining 33-14 species is build. Then the 19-species-model 
# is used to predict.
#
Foci_n <- list(Focus1=list(data="FourHM",descr=c("SigmaSQ")))
#
spe_jknife <- list()
for(j in names(Foci_n)) {
  spe_jknife[[j]]$Ypred <- LC50Chemical[[Foci_n[[j]]$data]] ; spe_jknife[[j]]$Ypred[] <- NA
  spe_jknife[[j]]$Z <- as.matrix(cbind(MainSpecies=1,Chemical.properties[match(HMetal[[Foci_n[[j]]$data]],
                    Chemical.properties[,"Chemical"]),Foci_n[[j]]$descr,drop=FALSE]))
  rownames(spe_jknife[[j]]$Z) <- colnames(LC50Chemical[[Foci_n[[j]]$data]])
}
#  
spe_jknife[[j]]$LC50 <- LC50Chemical[[Foci_n[[j]]$data]][-c(20,11,29,5,4,31,24,7,6,30,23,2,10,9),,drop=FALSE]
spe_jknife[[j]]$locs <- getGraphLocations(treeChemical[[Foci_n[[j]]$data]],rownames(LC50Chemical[[Foci_n[[j]]$data]])[c(20,11,29,5,4,31,24,7,6,30,23,2,10,9)])
spe_jknife[[j]]$PEMJN <- PEM.forcedSimple(y=spe_jknife[[j]]$LC50,x=NULL,w=spe_jknife[[j]]$locs$x,d="distance",sp="species",a=0,psi=1)
spe_jknife[[j]]$U <- as.matrix(cbind(MainProp=1,spe_jknife[[j]]$PEMJN$u))
n <- nrow(spe_jknife[[j]]$LC50) ; m <- ncol(spe_jknife[[j]]$LC50) ; p <- ncol(spe_jknife[[j]]$U) ; q <- ncol(spe_jknife[[j]]$Z)
spe_jknife[[j]]$ZU <- spe_jknife[[j]]$Z%x%spe_jknife[[j]]$U
dimnames(spe_jknife[[j]]$ZU) <- list(paste(rep(rownames(spe_jknife[[j]]$U),nrow(spe_jknife[[j]]$Z)),rep(rownames(spe_jknife[[j]]$Z),each=nrow(spe_jknife[[j]]$U)),sep="X"),
                                              paste(rep(colnames(spe_jknife[[j]]$U),ncol(spe_jknife[[j]]$Z)),rep(colnames(spe_jknife[[j]]$Z),each=ncol(spe_jknife[[j]]$U)),sep="X"))    
spe_jknife[[j]]$y <- spe_jknife[[j]]$LC50 ; dim(spe_jknife[[j]]$y) <- c(prod(dim(spe_jknife[[j]]$LC50)),1L)
dimnames(spe_jknife[[j]]$y) <- list(paste(rep(rownames(spe_jknife[[j]]$LC50),ncol(spe_jknife[[j]]$LC50)),
                                                   rep(colnames(spe_jknife[[j]]$LC50),each=nrow(spe_jknife[[j]]$LC50)),sep="X"),NULL)
spe_jknife[[j]]$lm1 <- sequentialAICcblm(y=spe_jknife[[j]]$y,ZU=spe_jknife[[j]]$ZU[,-1],k=2)
spe_jknife[[j]]$Cmat <- numeric(0) ; spe_jknife[[j]]$Cmat[c("(Intercept)",colnames(spe_jknife[[j]]$ZU)[-1])] <- 0
spe_jknife[[j]]$included <- coef(spe_jknife[[j]]$lm1) ; spe_jknife[[j]]$Cmat[labels(spe_jknife[[j]]$included)] <- spe_jknife[[j]]$included
dim(spe_jknife[[j]]$Cmat) <- c(ncol(spe_jknife[[j]]$U),ncol(spe_jknife[[j]]$Z))
dimnames(spe_jknife[[j]]$Cmat) <- list(colnames(spe_jknife[[j]]$U),colnames(spe_jknife[[j]]$Z))
spe_jknife[[j]]$CmatStrip <- spe_jknife[[j]]$Cmat[which(rowSums(spe_jknife[[j]]$Cmat!=0)!=0),]
spe_jknife[[j]]$scores <- Locations2PEMscores(spe_jknife[[j]]$PEMJN,spe_jknife[[j]]$locs)
spe_jknife[[j]]$scores$Uscores <- cbind(1,spe_jknife[[j]]$scores$scores) ; colnames(spe_jknife[[j]]$scores$Uscores) <- colnames(spe_jknife[[j]]$U)
colnames(spe_jknife[[j]]$CmatStrip)
Ypred_s<-as.matrix(spe_jknife[[j]]$scores$Uscores[,rownames(spe_jknife[[j]]$CmatStrip),drop=FALSE])%*%spe_jknife[[j]]$CmatStrip%*%t(spe_jknife[[j]]$Z[,colnames(spe_jknife[[j]]$CmatStrip),drop=FALSE])
#
Ypred_s2<- data.frame(Ypred_s,Species=rownames(Ypred_s)) 
Ypred_s2_p<- melt(Ypred_s2, id=("Species"))
Yobs_s2<-melt(Chemical.Table[,-c(2:3)], id=("ShortNames"))
#
Ymer_s2<- merge(Yobs_s2,Ypred_s2_p, by.x=c("variable","ShortNames"), by.y=c("variable","Species"))
names(Ymer_s2)[c(1,3:4)]<-c("HM","Yobs","Ypred")
#
summary(lm(Ymer_s2$Ypred~Ymer_s2$Yobs))# 0.46
Psquare(Ymer_s2$Yobs,Ymer_s2$Ypred)# 0.47
#
# Deviation factor
over_HM<-transform(Ymer_s2[Ymer_s2$Ypred > Ymer_s2$Yobs,], d= (5^(Ypred-Yobs)-1))
under_HM<-transform(Ymer_s2[Ymer_s2$Ypred < Ymer_s2$Yobs,], d= 1-(5^(Yobs-Ypred)))
a_HM<-rbind(over_HM,under_HM); nrow(Ymer_s2)==nrow(a_HM)
a_HM[abs(a_HM$d)>5,] # deviating species
nrow(a_HM[abs(a_HM$d)<5,])/nrow(a_HM) # 0.88
#
#
# Type 2: New metals
# Use FourHM model to predict LC50 for Ni and Pb. 
#
# Ypred=Utarget*C*Z(transpozed)
Ypred <- as.data.frame(ResList[["FourHM-SigmaSQ"]]$U[,rownames(ResList[["FourHM-SigmaSQ"]]$BLMfit$CmatStrip)]%*%ResList[["FourHM-SigmaSQ"]]$BLMfit$CmatStrip%*%t(ResList[["SixHM-SigmaSQ"]]$Z[c("Ni","Pb"),,drop=FALSE]))
Yobs <- Chemical.Table[match(Chemical.Table[,"ShortNames"],rownames(Ypred)),c("Ni","Pb")]
#
Ypred_np<- data.frame(Ypred,Species=rownames(Ypred)) 
Ypred_np_p<- melt(Ypred_np, id=("Species"))
Yobs_np<-melt(Chemical.Table[,-c(2:3)], id=("ShortNames"))
#
Ymer_np<- merge(Yobs_np,Ypred_np_p, by.x=c("variable","ShortNames"), by.y=c("variable","Species"))
names(Ymer_np)[c(1,3:4)]<-c("HM","Yobs","Ypred")
Ymer_np2<- Ymer_np[!is.na(Ymer_np$Yobs),]
#
# Overall
summary(lm(Ymer_np2$Ypred~Ymer_np2$Yobs)) # 0.42
Psquare(Ymer_np2$Yobs,Ymer_np2$Ypred)# 0.37
#
#
# Deviation factor
over_spe<-transform(Ymer_np2[Ymer_np2$Ypred >= Ymer_np2$Yobs,], d= (5^(Ypred-Yobs))-1)
under_spe<-transform(Ymer_np2[Ymer_np2$Ypred < Ymer_np2$Yobs,], d= 1-(5^(Yobs-Ypred)))
a_spe<-rbind(over_spe,under_spe); nrow(Ymer_np2)==nrow(a_spe)
a_spe[abs(a_spe$d)>5,] # deviating species
nrow(a_spe[abs(a_spe$d)<5,])/nrow(a_spe) # 90
#
# Figure 3
#
tiff(file="Fig3 Model Validations-q2.tiff", 
     compression="lzw", width = 6, height = 6, units = "in", res = 1200, type="cairo")
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
        mar = c(4, 4, 1, 1)) # space for one row of text at ticks and to separate plots
         
### Type 1A: Jacknife
plot(NA,xlim=c(0,6),ylim=c(0,6),axes=T,frame = F, 
     xlab="",
     ylab=expression(paste("Predicted ",log[10],LC[50]," (",mu*g%.%L^{-1},")",sep="")),main="")
abline(0,1, lty = 2); abline(coef(LC50_jknife[[1L]]$lm))
for(j in 1L:ncol(LC50Chemical[[Foci[[1L]]$data]]))
  points(y=LC50_jknife[[1L]]$Ypred[,j],x=LC50Chemical[[Foci[[1L]]$data]][,j],pch=j,cex=1.0)
legend(x=4,y=2.5,pch=c(1:3,6),legend=c("Cd","Cu","Hg","Zn"),cex=0.9)
text(1,4.0, bquote(q^2 == .(format(LC50_jknife[["Focus1"]]$Psquare, digits = 2))),cex=1.0)
text(1,5.5, "A",cex=1.5)
#
#
plot(NA,xlim=c(0,6),ylim=c(0,6),axes=T,frame = F, 
     xlab="",
     ylab="",main="")
for(j in 1L:ncol(LC50Chemical[[Foci[[2L]]$data]]))
  points(y=LC50_jknife[[2L]]$Ypred[,j],x=LC50Chemical[[Foci[[2L]]$data]][,j],pch=j,cex=1.0)
abline(0,1, lty = 2); abline(coef(LC50_jknife[[2L]]$lm))
legend(x=4.5,y=3.2,pch=1L:6,legend=Chemical.properties[match(colnames(LC50Chemical[[Foci[[2L]]$data]]),Chemical.properties[,"Chemical"]),"Chemical"],cex=0.9)
text(1,4.0, bquote(q^2 == .(format(LC50_jknife[["Focus2"]]$Psquare, digits = 2))),cex=1.0)
text(1,5.5, "B",cex=1.5)
#
# Type 1B: New Species
#
plot(NA,xlim=c(0,6),ylim=c(0,6),axes=T,frame = F,
     xlab=expression(paste("Experimental ",log[10],LC[50]," (",mu*g%.%L^{-1},")",sep="")),
     ylab=expression(paste("Predicted ",log[10],LC[50]," (",mu*g%.%L^{-1},")",sep="")),main="")
#
# 
mod_sp<-lm(Ypred~Yobs, data=Ymer_s2[3:4])
abline(0,1, lty = 2);abline(coef(mod_sp))
points(y=Ymer_s2$Ypred[Ymer_s2$HM=="Cd"],x=Ymer_s2$Yobs[Ymer_s2$HM=="Cd"],cex=1.0,pch=1)
points(y=Ymer_s2$Ypred[Ymer_s2$HM=="Cu"],x=Ymer_s2$Yobs[Ymer_s2$HM=="Cu"],cex=1.0,pch=2)
points(y=Ymer_s2$Ypred[Ymer_s2$HM=="Hg"],x=Ymer_s2$Yobs[Ymer_s2$HM=="Hg"],cex=1.0,pch=3)
points(y=Ymer_s2$Ypred[Ymer_s2$HM=="Zn"],x=Ymer_s2$Yobs[Ymer_s2$HM=="Zn"],cex=1.0,pch=6)
legend(x=4,y=2.2,pch=c(1:3,6),legend=c("Cd","Cu","Hg","Zn"),cex=0.9)
text(1,4.0, bquote(q^2 == .(format(Psquare(Ymer_s2$Yobs,Ymer_s2$Ypred), digits = 2))),cex=1.0)
text(1,5.5, "C",cex=1.5)
#
#
# Type 2: New Metals
#
plot(NA,xlim=c(0,6),ylim=c(0,6),axes=T,frame = F,
     xlab=expression(paste("Experimental ",log[10],LC[50]," (",mu*g%.%L^{-1},")",sep="")),
     ylab="",main="")
#
mod_hm<-lm(Ypred~Yobs, data=Ymer_np2)
abline(0,1, lty = 2); abline(coef(mod_hm));  m <- "Ni";n <- "Pb"
points(y=Ypred[!is.na(Yobs[,m]),m],x=Yobs[!is.na(Yobs[,m]),m],pch=4,cex=1.0)
points(y=Ypred[!is.na(Yobs[,n]),n],x=Yobs[!is.na(Yobs[,n]),n],pch=5,cex=1.0)
legend(x=4,y=1.5,pch=c(4:5),legend=c("Ni","Pb"),cex=0.9)
text(1,4.0, bquote(q^2 == .(format(Psquare(Ymer_np2$Yobs,Ymer_np2$Ypred), digits = 2))),cex=1.0)
text(1,5.5, "D",cex=1.5)
#
dev.off()
#
#
#
#  End 
#
###################################################################################################################

 


