#Install and Load Necessary Packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages(minfi)
BiocManager::install("IlluminaHumanMethylationEPICmanifest")

library(minfi)
library(BiocManager)
library(IlluminaHumanMethylationEPICmanifest)
library(dplyr)

#Generate Beta Values
baseDir<-"_______"
RGset<-read.metharray.exp(base=baseDir,targets=NULL,recursive=TRUE)
annotation(RGset)["array"] = "IlluminaHumanMethylationEPIC"
annotation(RGset)["annotation"] = "iln12.hg19"
Mset<-preprocessRaw(RGset)
ratioset<-ratioConvert(Mset,what="both",keepCN=F)
b<-getBeta(ratioset)
b<-t(b)

#Add and Clean Merge File
mergecsv<-read.csv("_____",header=T)
B<-b
rownames(mergecsv) <- mergecsv$X
mergecsv$X<- NULL
mergedcsv <-merge(mergecsv,B, by= 0)
mergedcsv$Row.names <- NULL
rownames(mergedcsv) <- mergedcsv$Sample_Name
mergedcsv$Sample_Name <- NULL

#Save Merged Samples
saveRDS(mergedcsv, "___") 
mergedcsv <- readRDS("___")

#Heatmap
library(gplots)
library(RColorBrewer)
heatmap(mergedcsv, col = colorRampPalette(brewer.pal(9, ""))(256), scale="none", margins=c(10,10))

#load EPIC annotation
epic.ann <- readRDS("____Epic_Annotation.rds")
tmerged <- t(mergedcsv)
USEQcsv <- merge(epic.ann,mergedcsv, by= 0)
USEQcsv$Row.names <- NULL
USEQ <- na.omit(USEQcsv)

##Use USEQ Tools here (see evernote)
chr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

for(i in colnames(mydf)){
  setwd("~path to working directory")
  dir.create(i)
  setwd(i)
  x<-mydf[,c("CHR","MAPINFO",i)]
  for(i in chr){
    c<-subset(x,CHR == i)
    colnames(c)<-NULL
    write.table(c,i,sep="\t",row.names=FALSE)}
  filez<-list.files()
  for(i in filez){
    file.rename(i,paste(i,".sgr",sep=""))}
}

#Generate Sample Data Frames
pe <- USEQ[,c("CHR","MAPINFO","19483","17017","19467","17710","17911","18159")]
iugr <- USEQ[,c("CHR","MAPINFO","18481","16814","18230","17420","19247","16437")]
control <- USEQ[,c("CHR","MAPINFO","19235","18944","18733","18243","19510","19315")]

#Creating Density Plots for the [control] subset
tcontrol <- t(control)
for(i in colnames(tcontrol))
{
  list<-tcontrol[,i]
  Controldensity<-density(list)
  plot(Controldensity)
  png(filename=paste("/Users/timjenkinslab/Desktop/carterjosh/Densityplots",i,".png",sep=""))
  plot(Controldensity)
  dev.off()
}

#Creating Density Plots for the [pe] subset
tpe <- t(pe)
tpe <- na.omit(tpe)
for(i in colnames(tpe))
{
  list<-tpe[,i]
  Controldensity<-density(list)
  plot(Controldensity)
  png(filename=paste("/Users/timjenkinslab/Desktop/carterjosh/peplots",i,".png",sep=""))
  plot(Controldensity)
  dev.off()
}
#Creating Density Plots for the [sperm] subset
tsperm <- t(sperm)
tsperm <- na.omit(tsperm)
for(i in colnames(tsperm))
{
  list<-tsperm[,i]
  Controldensity<-density(list)
  plot(Controldensity)
  png(filename=paste("/Users/timjenkinslab/Desktop/carterjosh/spermplots",i,".png",sep=""))
  plot(Controldensity)
  dev.off()
}
#Creating Density Plots for the [iugr] subset
tiugr <- t(iugr)
tiugr <- na.omit(tiugr)
for(i in colnames(tiugr))
{
  list<-tiugr[,i]
  Controldensity<-density(list)
  plot(Controldensity)
  png(filename=paste("/Users/timjenkinslab/Desktop/carterjosh/iugrplots",i,".png",sep=""))
  plot(Controldensity)
  dev.off()
}
#Point Data Analysis of Control and [pe] as a loop to get all the cg values
pePointDataResults<-NULL
i <- NULL

for(i in rownames(control))
{
  ControlList<-control[i,]
  peTreatmentList<-pe[i,]
  r<-cbind(CG=i,PVal=(t.test(ControlList,peTreatmentList)$p.value))
  pePointDataResults<-rbind(pePointDataResults,r)
}
write.csv(pePointDataResults,"___.csv")

#Point Data Analysis of Control and [iugr] as a loop to get all the cg values
iugrPointDataResults<-NULL

for(i in rownames(control))
{
  ControlList<-control[i,]
  iugrTreatmentList<-iugr[i,]
  r<-cbind(CG=i,PVal=(t.test(ControlList,iugrTreatmentList)$p.value))
  iugrPointDataResults<-rbind(iugrPointDataResults,r)
}
write.csv(iugrPointDataResults,"____.csv")

#Generate Box Plots within USEQ regions
#Given Region
peimchr2<- subset(pe, CHR=="chr2")
peim2<- subset(peimchr2, MAPINFO >="20101139" & MAPINFO <= "20102465")
test<- peim2
test<- subset(test, select= -c(CHR,MAPINFO))
petestmean <- apply(test,1,mean)
controlimchr2<- subset(control, CHR=="chr2")
controlim1<- subset(controlimchr2, MAPINFO >="20101139" & MAPINFO <= "20102465")
testing<- controlim1
testing<- subset(testing, select= -c(CHR,MAPINFO))
controltestmean <- apply(testing,1,mean)
png("___.png")
boxplot(petestmean,controltestmean,main= "Methylation at 20101139-20102465", col=c("red","blue"))
legend("topleft", c("Preeclampsia","Control"), col=c("red","blue"), lwd=c("2","2"))
dev.off()
#end given region

#Generate Global Mean Methylation Values (for bar charts)
tcontrol<- t(control)
tpe<- t(pe)
tiugr<- t(iugr)
fiveURS<- readRDS("~/Desktop/placenta/5URT.rds")
exons<- readRDS("~/Desktop/placenta/Exons.rds")
islands<- readRDS("~/Desktop/placenta/Island.rds")
shores<- readRDS("~/Desktop/placenta/Shores.rds")

#clean annotation groups
exons <- unique(exons)
row.names(exons) <- exons$IlmnID
exons$IlmnID=NULL
head(exons)

islands <- unique(islands)
row.names(islands) <- islands$IlmnID
islands$IlmnID=NULL

shores <- unique(shores)
row.names(shores) <- shores$IlmnID
shores$IlmnID=NULL

fiveURS <- unique(fiveURS)
row.names(fiveURS) <- fiveURS$IlmnID
fiveURS$IlmnID=NULL

#annotating pe group
tpe <- t(pe)
na.tpe <- na.omit(tpe)
merged.pe.exons <- merge(na.tpe,exons,by=0)
merged.pe.exons$UCSC_RefGene_Group = NULL
rownames(merged.pe.exons) <- merged.pe.exons$Row.names
merged.pe.exons$Row.names <- NULL
mpe <- apply(merged.pe.exons,1,mean)
tmpe <- t(mpe)
mean(tmpe)

merged.pe <- merge(na.tpe,islands,by=0)
merged.pe$Relation_to_UCSC_CpG_Island = NULL
rownames(merged.pe) <- merged.pe$Row.names
merged.pe$Row.names <- NULL
ipe <- apply(merged.pe,1,mean)
t.ipe <- t(ipe)
mean(t.ipe)

merged.pe <- merge(na.tpe,shores,by=0)
merged.pe$Relation_to_UCSC_CpG_Island = NULL
rownames(merged.pe) <- merged.pe$Row.names
merged.pe$Row.names <- NULL
spe <- apply(merged.pe,1,mean)
t.spe <- t(spe)
mean(t.spe)

merged.pe <- merge(na.tpe,fiveURS,by=0)
merged.pe$UCSC_RefGene_Group = NULL
rownames(merged.pe) <- merged.pe$Row.names
merged.pe$Row.names <- NULL
fpe <- apply(merged.pe,1,mean)
t.fpe <- t(fpe)
mean(t.fpe)

#Regionalize Data
###create objects to hold cumulative results
results.mean<-NULL
regions<- read.csv("/Users/carternorton/Desktop/placenta/CpGIslands.csv")

###Parse relavent objects by chr
df <- USEQ[,c("CHR", "MAPINFO", "19483","17017","19467","17710","17911","18159","18481","16814","18230","17420","19247","16437","19235","18944","18733","18243","19510","19315")]
df <- rename(df, Location = MAPINFO)

chromosomes = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

for(i in chromosomes){
  subset_regions = subset(regions, CHR == i)
  subset_df = subset(df, CHR == i)

  for (j in rownames(subset_regions)){
    working = subset_regions[j,]
    Start = working$Start
    Stop = working$Stop
    r = subset(subset_df, Location >= Start & Location <= Stop)
    r$CHR = NULL
    r$Location = NULL
    l = length(rownames(r))
    r = apply(r,2,mean)
    r = t(r)
    r = cbind(CHR = i, Start = Start, Stop = Stop, CpG.Count = l, r)
    results.mean = rbind(results.mean,r)
  }
}

#Separate into groups
region.pe<- results.mean[, c("19483","17017","19467","17710","17911","18159")]
region.iugr<- results.mean[,c("18481","16814","18230","17420","19247","16437")]
region.control<- results.mean[,c("19235","18944","18733","18243","19510","19315")]

saveRDS(results.mean, "~___.rds")

new.region.control <- na.omit(region.control)
write.csv(new.region.control, "~/Desktop/placenta/new.region.control.csv")
new.region.control <- read.csv("~/Desktop/placenta/new.region.control.csv")
new.region.control$X <- NULL
avg.control <- apply(new.region.control,1, mean)

new.region.pe <- na.omit(region.pe)
write.csv(new.region.pe, "~/Desktop/placenta/new.region.pe.csv")
new.region.pe <- read.csv("~/Desktop/placenta/new.region.pe.csv")
new.region.pe$X <- NULL
avg.pe <- apply(new.region.pe,1, mean)

new.region.iugr<- na.omit(region.iugr)
write.csv(new.region.iugr, "~/Desktop/placenta/new.region.iugr.csv")
new.region.iugr <- read.csv("~/Desktop/placenta/new.region.iugr.csv")
new.region.iugr$X <- NULL
avg.iugr <- apply(new.region.iugr,1, mean)

#Specific Sample Means
region.pe1 <- results.mean[, c("19483")]
region.pe2 <- results.mean[, c("17017")]
region.pe3 <- results.mean[, c("19467")]
region.pe4 <- results.mean[, c("17710")]
region.pe5 <- results.mean[, c("17911")]
region.pe6 <- results.mean[, c("18159")]
region.control1 <- results.mean[, c("19235")]

# High Density Scatterplot with Color Transparency

##Not Averaged Plots
png("~/Desktop/PEvs.Control.png", width=600, height=600)
plot(region.pe,region.control, main="PE vs. Control", col=rgb(0,60,100,10,maxColorValue=255), pch=8)
dev.off()

png("~/Desktop/IUGRvs.Control.png", width=600, height=600)
plot(region.iugr,region.control, main="IUGR vs. Control", col=rgb(100,60,0,10,maxColorValue=255), pch=8)
dev.off()

png("~/Desktop/PEvs.IUGR.png", width=600, height=600)
plot(region.pe,region.iugr, main="IUGR vs. PE", col=rgb(0,100,0,10,maxColorValue=255), pch=8)
dev.off()

#Averaged Plots 

png("~/Desktop/1new.PEvs.Control.png", width=600, height=600)
plot(avg.pe,avg.control, main="PE vs. Control", col=rgb(0,60,100,10,maxColorValue=255), pch=8)
dev.off()

png("~/Desktop/new.IUGRvs.Control.png", width=600, height=600)
plot(avg.iugr,avg.control, main="IUGR vs. Control", col=rgb(100,60,0,10,maxColorValue=255), pch=8)
dev.off()

png("~/Desktop/new.PEvs.IUGR.png", width=600, height=600)
plot(avg.pe,avg.iugr, main="IUGR vs. PE", col=rgb(0,100,0,10,maxColorValue=255), pch=8)
dev.off()

