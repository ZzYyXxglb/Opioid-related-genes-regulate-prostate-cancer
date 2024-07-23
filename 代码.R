exp<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header = T,row.names = 1,check.names = F)
sampleID<-t(exp[1,])
sampleID<-as.data.frame(sampleID[-1,])
normal <- 0 
tumour <- 0
for (i in 1:551) {
  num <- as.numeric(substring(sampleID[i,1],14,15))
  if (num %in% seq(1,9)) {
    sampleID[i,2] <- "T"
    tumour <- tumour+1}
  if (num %in% seq(10,29)) {
    sampleID[i,2] <- "N"
    normal <- normal+1}
}
normal
tumour
sampletype<-sampleID
colnames(sampletype)<-c('sampleID','type')
write.table(sampletype,file = '03TCGA_PRAD_sampletype.txt',sep = '\t',col.names = T)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("DESeq2")
require(limma)
require(DESeq2)
gene_exprSet<-read.table(file = '12gene_exprSet_PRAD.txt',sep = '\t',header=T,row.names=1,check.names=F)
dimnames<-list(rownames(gene_exprSet),colnames(gene_exprSet))
gene_exprSet<-matrix(as.numeric(as.matrix(gene_exprSet)),nrow=nrow(gene_exprSet),dimnames=dimnames)
gene_exprSet<-gene_exprSet[rowMeans(gene_exprSet)>2,]
gene_exprSet<-round(gene_exprSet,0) 
group<-sapply(strsplit(colnames(gene_exprSet),"\\-"),"[",4)
group<-sapply(strsplit(group,""),"[",1)
design<-factor(group)
group<-as.data.frame(group)
group<-t(group)
colnames(group)<-colnames(gene_exprSet)
group<-t(group)
dds<-DESeqDataSetFromMatrix(countData = gene_exprSet,colData = group,design = ~group)
dds2<-DESeq(dds)
resultsNames(dds2)
sizefactor<-sizeFactors(dds2)
sizefactor<-as.data.frame(sizefactor)
res <-results(dds2, contrast=c('group',"0","1"))
head(res)
resorder<-res[order(res$padj),]
head(resorder)
resorder<-as.data.frame(resorder)
View(resorder)
write.table(resorder,file = '13DESeq_PRAD_NT.txt',sep = '\t',quote = F,row.names = T)
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
require(dplyr)
DESeq_diff<-DESeq%>%
  dplyr::filter(DESeq$padj<0.05)
DESeq_up<-DESeq_diff%>%
  filter(DESeq_diff$log2FoldChange>0.5)
DESeq_down<-DESeq_diff%>%
  filter(DESeq_diff$log2FoldChange<(-0.5))
DESeq_filt<-rbind(DESeq_up,DESeq_down)
write.table(DESeq_filt,file = '30DESeq_NT_filt.txt',sep = '\t',col.names = T)
gene_exprSet<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header = T,row.names=1,check.names = F)
sample<-read.table(file = '03NTpairTCGA_03.txt',sep = '\t',header = T)
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
require(dplyr)
DESeq_diff<-DESeq%>%
  dplyr::filter(DESeq$padj<0.05)
DESeq_diff<-DESeq_diff[order(DESeq_diff$log2FoldChange,decreasing = T),]
DESeq_up<-DESeq_diff[1:200,]
DESeq_diff<-DESeq_diff[order(DESeq_diff$log2FoldChange,decreasing = F),]
DESeq_down<-DESeq_diff[1:200,]
DESeq_filt<-rbind(DESeq_up,DESeq_down)
require(dplyr)
rt<-gene_exprSet[rownames(DESeq_filt),]
min(rt)
max(rt)
rt<-log2(rt+1-min(rt))
min(rt)
max(rt)
heatmap(as.matrix(rt),Colv= NA, 
        col=colorRampPalette(c('blue','white','red'))(100)) 
heatmap(as.matrix(rt), 
        col=colorRampPalette(c('blue','white','red'))(100)) 

require(dplyr)
diff<-DESeq
allDiff=diff[is.na(diff$padj)==FALSE,]
xMax=max(-log10(allDiff$padj))+1
yMax=10
plot(-log10(allDiff$padj), allDiff$log2FoldChange, xlab="-log10(padj)",ylab="log2FoldChange",
     xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.8)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange>foldChange,]#
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col='orangered3',cex=0.8)
diffSub=allDiff[allDiff$padj<padj & allDiff$log2FoldChange<(-foldChange),]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="steelblue2",cex=0.8)
abline(h=0,lty=2,lwd=1)
opi_DESeq$change<-ifelse(opi_DESeq$log2FoldChange>=1.5&opi_DESeq$padj<0.001,'Up',
                         ifelse(opi_DESeq$log2FoldChange<=-1.5&opi_DESeq$padj<0.001,'Down','NotSignificant'))
table(opi_DESeq$change)
require(ggplot2)
ggplot(opi_DESeq,aes(x=opi_DESeq$log2FoldChange,y=-log10(opi_DESeq$padj)))+
  geom_point()
ggplot(opi_DESeq,aes(x=log2FoldChange,y=-log10(padj)))+
  geom_point(aes(color=change))+
  scale_color_manual(values = c('dodgerblue','gray','firebrick'))
require(clusterProfiler)
library(org.Hs.eg.db)    
library(ggplot2)    
require(dplyr)
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
DESeq_diff<-DESeq%>%
  dplyr::filter(DESeq$padj<0.05)
write.table(DESeq_diff,file = '11DESeq_03_607.txt',sep = '\t',col.names = T,row.names = T)
geneNames<-rownames(DESeq_diff)
gene <-  mapIds(org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL')    
BP.params <- enrichGO(   gene   = gene,    
                         OrgDb  = org.Hs.eg.db,    
                         ont   = "BP"  ,    
                         pAdjustMethod = "BH",    
                         pvalueCutoff  = 0.01,    
                         qvalueCutoff  = 0.01)    
MF.params <- enrichGO(   gene   = gene,    
                         OrgDb  = org.Hs.eg.db,    
                         ont   = "MF"  ,    
                         pAdjustMethod = "BH",    
                         pvalueCutoff  = 0.01,    
                         qvalueCutoff  = 0.01)   
CC.params <- enrichGO(   gene   = gene,    
                         OrgDb  = org.Hs.eg.db,    
                         ont   = "CC"  ,    
                         pAdjustMethod = "BH",    
                         pvalueCutoff  = 0.01,    
                         qvalueCutoff  = 0.01)    
BP.list <- setReadable(BP.params, org.Hs.eg.db, keyType = "ENTREZID") 
MF.list <- setReadable(MF.params, org.Hs.eg.db, keyType = "ENTREZID")     
CC.list <- setReadable(CC.params, org.Hs.eg.db, keyType = "ENTREZID")     
BPresult<-as.data.frame(BP.list@result)
MFresult<-as.data.frame(MF.list@result)
CCresult<-as.data.frame(CC.list@result)
rt<-rbind(MFresult[1:10,],BPresult[1:10,])
BPresult<-filter(BPresult,p.adjust<0.01&qvalue<0.01)
MFresult<-filter(MFresult,p.adjust<0.01&qvalue<0.01)
CCresult<-filter(CCresult,p.adjust<0.01&qvalue<0.01)
a<-colnames(rt)
rt<-rt[,-1]
colnames(rt)<-a[1:10]
ggplot(data=rt)+
  geom_bar(aes(x=Description,y=GeneRatio, fill=p.adjust), stat='identity') + 
  coord_flip() + scale_x_discrete(limits=rt$Description) 
GOresult<-read.table(file = '13GO_04DAVID607chart.txt',sep = '\t',header = T,check.names = F)
require(dplyr)
require(clusterProfiler)
library(org.Hs.eg.db)    
library(ggplot2)    
require(dplyr)
BPresult<-filter(GOresult,Category=='GOTERM_BP_DIRECT'&PValue<0.05)
MFresult<-filter(GOresult,Category=='GOTERM_MF_DIRECT'&PValue<0.05)
CCresult<-filter(GOresult,Category=='GOTERM_CC_DIRECT'&PValue<0.05)
rt<-BPresult[c(1:5,7,15,24,25,62),]
rt<-MFresult[1:10,]
rt<-rbind(MFresult[1:10,],BPresult[c(1:5,7,15,24,25,62),])
ggplot(data=rt)+
  geom_bar(aes(x=Term,y=Count, fill=PValue), stat='identity') + 
  coord_flip() + scale_x_discrete(limits=rt$Term) +ylim(0,50)
require(dplyr)
setwd('D:/Paper/opioid/NT_Deseq')
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
DESeq_diff<-DESeq%>%
  dplyr::filter(DESeq$padj<0.05)
DESeq_diff<-DESeq_diff%>%
  filter(abs(log2FoldChange)>0.5)
geneshot<-read.table(file = '05geneshot_opioid.tsv',sep = '\t',header = T,check.names = F)
opioid<-geneshot[1:200,2]
opioid<-DESeq_diff[opioid,]
opioid<-na.omit(opioid)
opioid<-opioid[order(opioid$log2FoldChange,decreasing = T),]
opioid<-rownames(opioid)[1:7]
gene_exprSet<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header=T,row.names=1,check.names=F)
gene_exprSet<-as.data.frame(gene_exprSet)
rt<-gene_exprSet[opioid,]
rt<-t(rt)
sampletype<-read.table(file = '03TCGA_PRAD_sampletype.txt',sep = '\t',header = T, check.names = F)
rt<-cbind(rownames(rt),rt)
colnames(rt)[1]<-'sampleID'
rt<-as.data.frame(rt)
require(dplyr)
rt<-sampletype%>%
  inner_join(rt,by='sampleID')
write.table(rt,file = '14DESeq_05opinon.txt',sep = '\t',col.names = T)
require(WGCNA)
require(impute)
require(reshape2)
require(stringr)
type<-'signed'
corType<-'pearson' #
corFnc<-ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
gene_exprSet<-read.table(file = '12gene_exprSet_PRAD.txt',sep = '\t',header=T,row.names=1,check.names=F)
dimnames<-list(rownames(gene_exprSet),colnames(gene_exprSet))
gene_exprSet<-matrix(as.numeric(as.matrix(gene_exprSet)),nrow=nrow(gene_exprSet),dimnames=dimnames)
gene_exprSet<-gene_exprSet[rowMeans(gene_exprSet)>2,]
gene_exprSet<-log2(gene_exprSet+1) 
group<-read.table(file = '03TCGA_PRAD_sampletype.txt', sep = '\t', header = T, check.names = F)
require(dplyr)
group<-group%>%
  dplyr::filter(type=='T')
group<-group[,1]
gene_exprSet<-gene_exprSet[,group]
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt', sep = '\t', header = T, check.names = F)
require(dplyr)
DESeq_diff<-DESeq%>%
  dplyr::filter(DESeq$padj<0.05)
DESeq<-rbind(DESeq_up,DESeq_down)
gene_exprSet<-gene_exprSet[row.names(DESeq),]
gene_exprSet<-gene_exprSet[order(rowSums(gene_exprSet),decreasing = T),]
gene_exprSet<-t(gene_exprSet)
nGenes<-ncol(gene_exprSet)
nSamples<-nrow(gene_exprSet)
#ɾ??ȱʧֵ
gsg<-goodSamplesGenes(gene_exprSet, verbose = 3)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
A=adjacency(t(gene_exprSet),type="signed")
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)
thresholdZ.k=-2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
gene_exprSet=gene_exprSet[!remove.samples,]
datTraits=datTraits[!remove.samples,]
write.table(gene_exprSet,file = '23gene_exp_WGCNA.txt', sep = '\t',col.names = T,row.names = T)
sampleTree<-hclust(dist(gene_exprSet),method = 'average')
plot(sampleTree, 
     main = 'Sample clustering to detect outliers', 
     sub = '', xlab = '',cex.lab=0.01)
powers<-c(c(1:10),seq(from=12,to=30,by=2))
sft<-pickSoftThreshold(gene_exprSet,powerVector=powers,verbose=5)
par(mfrow=c(1,2))
cex1 <-0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power <- sft$powerEstimate
power
net = blockwiseModules(gene_exprSet, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0('AS-green-FPKM-TOM'),
                       verbose = 3)
save(net,gene_exprSet,file = '17WGCNA.Rdata')
load(file = '17WGCNA.Rdata')
table(net$colors)
moduleLabels<-net$colors
moduleColors<-labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]],
                    'Module colors',
                    dendroLabels = F,hang=0.03,
                    addGuide = T, guideHang = 0.05)
ADJ<-adjacency(gene_exprSet,power = power)
vis<-exportNetworkToCytoscape(ADJ,edgeFile = 'edg1.txt',nodeFile = 'node.txt',threshold = 0.4)
MEs<-net$MEs
MEs_col<-MEs
colnames(MEs_col)<-paste0('ME',labels2colors(
  as.numeric(str_replace_all(colnames(MEs),'ME',''))
))
MEs_col<-orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col,'Eigengene adjacency heatmap',
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2),plotDendrograms = T,
                      xLabelsAngle = 90)
load(file = '17WGCNA.Rdata')
power<-12
table(net$colors)
moduleLabels<-net$colors
moduleColors<-labels2colors(moduleLabels)
load(net$TOMFiles[1],verbose = T)
TOM<-as.matrix(TOM)
plotTOM<-1-(dissTOM^power)
diag(plotTOM)<-NA
TOMplot(plotTOM,net$dendrograms,moduleColors,
        main='Network heatmap plot, all genes',
        #Colors = gplots::colorpanel(250,'red','orange','lemonchiffon')
)
select = sample(nGenes, size = nSelect); 
selectTOM = dissTOM[select, select]; 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = moduleColors[select]; 
sizeGrWindow(9,9) 
plotDiss = 1-(selectTOM^power); 
diag(plotDiss) = NA; 
TOMplot(plotDiss, 
        selectTree, 
        selectColors, 
        main = "Network heatmap plot, selected genes",
        #Colors = gplots::colorpanel(400,'red','orange','lemonchiffon')
) 

ADJ1<-abs(cor(gene_exprSet,use = 'p'))^6
moduleLabels<-net$colors
ALLdegrees1<-intramodularConnectivity(ADJ1,moduleLabels)
write.table(ALLdegrees1,file = '32intramodular_connectivit.txt',sep = '\t',row.names = T,col.names = T)
genecolor<-read.table(file = '27module_gene.txt',sep = '\t',header = T,check.names = F)
require(dplyr)
ALLdegrees1<-cbind(rownames(ALLdegrees1),ALLdegrees1)
colnames(ALLdegrees1)[1]<-'datSummary'
merge<-inner_join(ALLdegrees1,genecolor,by='datSummary')
colnames(merge)[1]<-'gene'
write.table(merge,file = '32intramodular_connectivit.txt',sep = '\t',row.names = F,col.names = T)

load(file = '17WGCNA.Rdata')
load(net$TOMFiles[1],verbose = T)
TOM<-as.matrix(TOM)
moduleLabels<-net$colors
moduleColors<-labels2colors(moduleLabels)
probes<-colnames(gene_exprSet)
dimnames(TOM)<-list(probes,probes)
cyt<-exportNetworkToCytoscape(TOM,threshold = 0,
                              nodeNames = probes,nodeAttr = moduleColors)
cyt$edgeData
cyt$nodeData
save(cyt,file = '18cytoscape.Rdata')
edgedata<-cyt$edgeData
nodedata<-cyt$nodeData
require(dplyr)
rt<-filter(edgedata,edgedata$weight>0.2)
write.table(rt,file = '06WGCNA_05edgedata_filt.txt',sep = '\t',col.names = T)
write.table(nodedata,file = '06WGCNA_05nodedata.txt',sep = '\t',col.names = T)
=genecolor<-read.table(file = '27module_gene.txt',sep = '\t',header = T,check.names = F)
require(dplyr)
module<-'black'
genelist<-filter(genecolor,datColor==module)
genelist<-genelist[,1]
i<-1
for (i in 1:nrow(rt)) {
  a<-(rt[i,1]%in%genelist&rt[i,2]%in%genelist)
  rt[i,6]<-a
}
require(dplyr)
result<-filter(rt,toAltName=='TRUE')
result[,6]<-'NA'
filename<-paste('06WGCN_06',module,'edgedata.csv',sep = '')
write.csv(result,file = filename)
=require(WGCNA)
require(impute)
require(reshape2)
require(stringr)
type<-'signed'
corType<-'pearson' #
corFnc<-ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
clinical<-read.table(file = '19TCGA_PRAD_phenotype.txt',sep = '\t',header = T,check.names = F)
survival<-read.table(file = '22survival_chen.txt',sep = '\t',header = T,check.names = F)
clinical<-clinical[!duplicated(clinical[,3]),]
gene_exprSet<-read.table(file = '23gene_exp_WGCNA.txt',sep = '\t',header = T,check.names = F)
k<-row.names(gene_exprSet)
k<-substr(k,1,12)
row<-match(k,clinical[,3])
dcli<-clinical[row,]
dcli<-dcli[!duplicated(dcli[,3]),]
row.names(dcli)<-dcli$sampleID
dcli<-dcli[,-1]
dclic<-subset(dcli,select = c('age_at_initial_pathologic_diagnosis','biochemical_recurrence',
                              'clinical_M','clinical_T','gleason_score','pathologic_N','pathologic_T',
                              'psa_value'))
dclic<-cbind(rownames(dclic),dclic)
colnames(dclic)[1]<-'sample'
require(dplyr)
dclic<-dplyr::inner_join(dclic,survival,by='sample')
colnames(dclic)[2]<-'age'
rownames(dclic)<-dclic[,1]
dclic<-dclic[,-1]
write.table(dclic,file = '24clinical_WCGNA.txt',sep = '\t',col.names = T,row.names = T)
for (i in 1:489) {
  if (dclic[i,"biochemical_recurrence"]=='YES') {dclic[i,"biochemical_recurrence"]=1}
  if (dclic[i,"biochemical_recurrence"]=='NO') {dclic[i,"biochemical_recurrence"]=0}
  if (dclic[i,"clinical_M"]=='M0') {dclic[i,"clinical_M"]=0}
  if (dclic[i,"clinical_M"]=='M1') {dclic[i,"clinical_M"]=1}
  if (dclic[i,"clinical_M"]=='M1c') {dclic[i,"clinical_M"]=1}
  if (dclic[i,"clinical_M"]=='M1a') {dclic[i,"clinical_M"]=1}
  if (dclic[i,"clinical_T"]=='T1a') {dclic[i,"clinical_T"]=1}
  if (dclic[i,"clinical_T"]=='T1b') {dclic[i,"clinical_T"]='1'}
  if (dclic[i,"clinical_T"]=='T1c') {dclic[i,"clinical_T"]='1'}
  if (dclic[i,"clinical_T"]=='T2a') {dclic[i,"clinical_T"]='2'}
  if (dclic[i,"clinical_T"]=='T2b') {dclic[i,"clinical_T"]='2'}
  if (dclic[i,"clinical_T"]=='T2c') {dclic[i,"clinical_T"]='2'}
  if (dclic[i,"clinical_T"]=='T2') {dclic[i,"clinical_T"]='2'}
  if (dclic[i,"clinical_T"]=='T3a') {dclic[i,"clinical_T"]='3'}
  if (dclic[i,"clinical_T"]=='T3b') {dclic[i,"clinical_T"]='3'}
  if (dclic[i,"clinical_T"]=='T4') {dclic[i,"clinical_T"]='4'}
  if (dclic[i,"pathologic_N"]=='N0') {dclic[i,"pathologic_N"]='0'}
  if (dclic[i,"pathologic_N"]=='N1') {dclic[i,"pathologic_N"]='1'}
  if (dclic[i,"pathologic_T"]=='T2a') {dclic[i,"pathologic_T"]='2'}
  if (dclic[i,"pathologic_T"]=='T2b') {dclic[i,"pathologic_T"]='2'}
  if (dclic[i,"pathologic_T"]=='T2c') {dclic[i,"pathologic_T"]='2'}
  if (dclic[i,"pathologic_T"]=='T3a') {dclic[i,"pathologic_T"]='3'}
  if (dclic[i,"pathologic_T"]=='T3b') {dclic[i,"pathologic_T"]='3'}
  if (dclic[i,"pathologic_T"]=='T4') {dclic[i,"pathologic_T"]='4'}
  if (dclic[i,"pathologic_T"]=='[Discrepancy]') {dclic[i,"pathologic_T"]=''}
}
write.table(dclic,file = '24clinical_WCGNA.txt',sep = '\t',col.names = T,row.names = T)
dclic<-read.table(file = '24clinical_WCGNA.txt',sep = '\t',header = T,check.names = F)
require(flashClust)
gene_exprSet<-cbind(rownames(gene_exprSet),gene_exprSet)
gene_exprSet[,1]<-substr(gene_exprSet[,1],1,12)
gene_exprSet<-gene_exprSet[!duplicated(gene_exprSet$`rownames(gene_exprSet)`),]
rownames(gene_exprSet)<-gene_exprSet[,1]
gene_exprSet<-gene_exprSet[,-1]
A<-adjacency(t(gene_exprSet),type = 'distance')
k<-as.numeric(apply(A,2,sum))-1

dclic_filt<-dclic[,c(1,2,3,5,6,7,8)]
colnames(dclic_filt)<-c('Age','Biochemincal Recurrence','M stage','Gleason Score','N Stage','T Stage','PSA')
dclic_filt<-dclic_filt[,c(1,6,5,3,2,4,7)]
#write.table(dclic_filt,file = '26clinicfordendrotrait.txt',sep = '\t',col.names = T,row.names = T)
traitColors<-data.frame(numbers2colors(dclic_filt,signed = F,
                                       naColor = 'grey'))
dimnames(traitColors)[[2]]<-names(dclic_filt)
outlierColor<-ifelse(Z.k<thresholdZ.k,'red','blue')
datColor<-data.frame(outlier=outlierColor,traitColors)
plotDendroAndColors(sampleTree,
                    groupLabels = names(datColor),=
                    colors = datColor,=
                    main='Sample dendrogram and trait heatmap',
                    cex.colorLabels = 1,cex.rowText = 0.4,cex.dendroLabels = 0.1
)
dataExp<-gene_exprSet
load(file = '17WGCNA.Rdata')
moduleLabelsAutomatic<-net$colors
moduleLabelsAutomatic<-labels2colors(moduleLabelsAutomatic)
moduleColorsFemale<-moduleLabelsAutomatic
t_rt<-dataExp
nGenes<-ncol(t_rt)
nSample<-nrow(t_rt)
MEs0<-moduleEigengenes(t_rt,moduleColorsFemale)$eigengenes
MEsFemales<-orderMEs(MEs0)
modTraitCor<-cor(MEsFemales,dclic,use = 'p')
modTraitP<-corPvalueStudent(modTraitCor,nSample)
textMatrix<-paste(signif(modTraitCor,2),'\n(',signif(modTraitP,1),')',sep = '')
dim(textMatrix)<-dim(modTraitCor)
par(mar=c(4,8,1,4))
labeledHeatmap(Matrix = modTraitCor,xLabels = colnames(dclic),
               yLabels = names(MEsFemales),colorLabels = F,colors = blueWhiteRed(50),
               ySymbols = names(MEsFemales),textMatrix = textMatrix,setStdMargins = F,
               cex.text = 0.6,zlim=c(-1,1),main=paste('Module-trait relationships'),font.lab.x = 0.5,font.lab.y = 0.5)
blocknumber<-1
datColor<-data.frame(moduleLabelsAutomatic)[net$blockGenes[[blocknumber]],]
datKME<-signedKME(t_rt,MEsFemales)
datSummary<-colnames(dataExp)
datout<-data.frame(datSummary,datColor,datKME)
write.table(datout,file = '27module_gene.txt',sep = '\t',col.names = T,row.names = T)
load(file = '17WGCNA.Rdata')
MEs<-net$MEs
MEs_col<-MEs
colnames(MEs_col)<-paste0('ME',labels2colors(
  as.numeric(str_replace_all(colnames(MEs),'ME',''))
))
MEs_col<-orderMEs(MEs_col)
nSample<-nrow(gene_exprSet)
if (corType=='pearson'){
  geneModuleMembership<-as.data.frame(cor(gene_exprSet,MEs_col,use='p'))
  MMPvalue<-as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples = nrow(gene_exprSet)))
}else{
  geneModuleMembershipA<-bicorAndPvalue(gene_exprSet,MEs_col,robustY=robustY)
  geneModuleMembership<-geneModuleMembershipA$bicor
  MMPValue<-geneModuleMembershipA$p
}
if(corType=='pearson'){
  geneTraitCor=as.data.frame(cor(dataExp,dclic,use = 'p'))
  genetraitP=as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor),nSamples = nrow(dataExp)
  ))
}else{
  geneTraitCorA=bicorAndPvalue(gene_exprSet,dclic,robustY=robustY)
  geneTraitCor=as.data.frame(geneTraitCorA$bicor)
  geneTraitP=as.data.frame(geneTraitCorA$p)
}
table(net$colors)
module<-'red'
pheno<-'RFS_score'
modNames<-substring(colnames(MEs_col),3)

moduleLabels<-net$colors
moduleColors<-labels2colors(moduleLabels)
module_column<-match(module,modNames)
pheno_column<-match(pheno,colnames(dclic))
moduleGenes<-moduleColors==module
sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,module_column]),
                   abs(geneTraitCor[moduleGenes,pheno_column]),
                   xlab = paste('Module Membership in',module,'module'),
                   ylab = paste('Gene significance for',pheno),
                   main = paste('Module Membership vs. Gene significant\n'),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module
)


load(file = '17WGCNA.Rdata')
require(WGCNA)
moduleLabelsAutomatic<-net$colors
moduleLabelsAutomatic<-labels2colors(moduleLabelsAutomatic)
eigengene<-moduleEigengenes(gene_exprSet,moduleLabelsAutomatic)$eigengenes
write.table(eigengene,file = '28WGCNA_eigengene.txt',sep = '\t',col.names = T,row.names = T)
survival<-read.table(file = '22survival_chen.txt',sep = '\t',header = T,check.names = F)
k<-substr(rownames(eigengene),1,15)
eigengene<-cbind(k,eigengene)
eigengene<-eigengene[!duplicated(eigengene$k),]
colnames(eigengene)[1]<-'sample'
require(dplyr)
eigengene_survival<-dplyr::inner_join(survival,eigengene,by='sample')
write.table(eigengene_survival,file = '29eigengene_survival.txt',sep = '\t',col.names = T)
require(survival)
rt<-eigengene_survival
geneset<-colnames(rt)[6:16]
pvalueset<-as.data.frame(colnames(rt)[6:16])
for (i in 6:16) {
  gene<-colnames(rt)[i]
  time<-'os'
  status<-'os_score'
  a<-rt[,gene]<median(rt[,gene])
  time<-as.numeric(rt[,time])
  status<-as.numeric(rt[,status])
  diff<-survdiff(Surv(time, status) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(time, status) ~ a, data = rt)
  summary(fit)    
  plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
       main=paste(gene,'os','(p=', pValue ,")",sep=" "))
  pvalueset[(i-5),2]<-pValue
}

geneset<-colnames(rt)[6:16]
for (i in 6:16) {
  gene<-colnames(rt)[i]
  time<-'RFS'

  status<-'RFS_score'
  a<-rt[,gene]<median(rt[,gene])
  time<-as.numeric(rt[,time])
  status<-as.numeric(rt[,status])
  diff<-survdiff(Surv(time, status) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(time, status) ~ a, data = rt)
  summary(fit)    
  plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
       main=paste(gene,'RFS','(p=', pValue ,")",sep=" "))
  pvalueset[(i-5),3]<-pValue
}
colnames(pvalueset)<-c('gene','OS_Pvalue','RFS_Pvalue')
pvalueset
write.table(pvalueset,file = '31survival_eigengene.txt',sep = '\t',col.names = T,row.names = F)


require(WGCNA)
require(stringr)
require(dplyr)
load(file = '17WGCNA.Rdata')
moduleColors <- labels2colors(net$colors)
ADJ1<-abs(cor(gene_exprSet,use="p"))^6
Alldegrees1<-intramodularConnectivity(ADJ1, moduleColors)
MEs_col<-net$MEs
colnames(MEs_col)<-paste0('ME',labels2colors(
  as.numeric(str_replace_all(colnames(net$MEs),'ME',''))
))
MEs_col<-orderMEs(MEs_col)
datKME<-signedKME(gene_exprSet, MEs_col, outputColumnName="MM.")
write.table(datKME,file = "09hubgene_01MM.txt",sep = '\t',col.names = T,row.names = T)
datKME<-read.table(file = "09hubgene_01MM.txt",sep = '\t',header = T,check.names = F)
dclic<-read.table(file = '24clinical_WCGNA.txt',sep = '\t',header = T,check.names = F)
clinicaltrait<-as.data.frame(dclic[,11]) 
GS_value<-0.25
MM_value<-0.8
rownames(gene_exprSet)<-substr(rownames(gene_exprSet),1,15) 
GS1<-as.numeric(cor(clinicaltrait,gene_exprSet, use="p"))
color<-colnames(datKME)
hub<-matrix(1:2,nrow=1,ncol=2,byrow=T)
colnames(hub)<-c('hubgene','module')
for (i in 1:11) {
  FilterGenes<-GS1>GS_value &abs(datKME[,color[i]]) > MM_value
  x<-datKME[FilterGenes,]
  if(nrow(x)!=0)
  {rt<-as.data.frame(rownames(x))
  rt[,2]<-color[i]
  colnames(rt)<-c('hubgene','module')
  hub<-rbind(hub,rt)}
}
hub<-hub[-1,]
hub[,3]<-'up'
colnames(hub)[3]<-'regulated'
for (i in 1:11) {
  FilterGenes<-((GS1<(-GS_value)) &abs(datKME[,color[i]]) > MM_value)
  x<-datKME[FilterGenes,]
  if(nrow(x)!=0){rt<-as.data.frame(rownames(x))
  rt[,2]<-color[i]
  rt[,3]<-'down'
  colnames(rt)<-c('hubgene','module','regulated')
  hub<-rbind(hub,rt)}
}
write.table(hub,file = '09hubgene_04cMap.txt',sep = '\t',col.names = T,row.names = F)
=require(WGCNA)
require(impute)
require(reshape2)
require(stringr)
load('17WGCNA.Rdata')
table(net$colors)
moduleLabels<-net$colors
moduleColors<-labels2colors(moduleLabels)
which.module<-"black"
rt<-t(scale(gene_exprSet[,moduleColors==which.module]))
heatmap(as.matrix(rt),Rowv= NA,
        col=colorRampPalette(c('green','black','red'))(100)) 

sizeGrWindow(8,7)
which.module<-"black"
ME<-MEs_col[,paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
par(mar=c(5, 4.2, 0, 0.7))
rt<-t(scale(gene_exprSet[,moduleColors==which.module]))
heatmap(as.matrix(rt),Rowv= NA,
        col=colorRampPalette(c('green','black','red'))(100)) 
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="samples")

setwd('D:/Paper/opioid/NT_Deseq')
rt<-read.table(file = '10cMap_02red.txt',sep = '\t',header = T,check.names = F)
colnames(rt)<-c('Compounds','Cell_line','Connective_Score','FDR','Norm_CS','Module')
require(ggplot2)
library(ggthemes)
require(dplyr)
a<-paste(rt[,1],rt[,2],sep='_')
rt[,1]<-a
a<-rev(a)
rt$Compounds<-factor(rt$Compounds,levels = a)#????˳?򣬷?????????״ͼ?ᰴ????ĸ˳??????
ggplot(data = rt)+
  geom_bar(aes(x= Compounds,y=Connective_Score,fill=FDR),stat = 'identity')+
  coord_flip()+
  scale_fill_gradient(low = 'blue',high = 'red')+
  xlab('')+
  ylab("") + 
  theme(axis.text.x=element_text(color="black", size=12), 
        axis.text.y=element_text(color="black", size=12)) + 
  scale_y_continuous(expand=c(0, 0)) + 
  scale_x_discrete(expand=c(0,0))
rt %>%
  mutate(Legend = factor(Compounds, levels = a)) %>%
  ggplot(aes(fill = FDR, y = Connective_Score, x = Compounds)) +
  geom_bar(position = "dodge", stat = "identity")
library(vioplot)   
rt <-read.table("10cMap_02red.txt",sep="\t",header=T,check.names=F)   #??ȡ?????ļ?
colnames(rt)<-c('Compounds','Cell_line','Connective_Score','FDR','Norm_CS','Module')
require(ggplot2)
library(ggthemes)
require(dplyr)
a<-paste(rt[,1],rt[,2],sep='_')
rt[,1]<-a
a<-rev(a)
rt$Compounds<-factor(rt$Compounds,levels = a)
ggplot(rt,aes(x=Compounds,y=Norm_CS,fill=Cell_line))+
  geom_bar(position = 'dodge',stat = 'identity')+
  xlab('Compounds')+
  ylab('Norm_CS')+
  labs(fill='Cell line')+
  coord_flip()+
  theme(axis.text.x=element_text(color="black", size=12), 
        axis.text.y=element_text(color="black", size=12)) + 
  scale_y_continuous(expand=c(0, 0)) + 
  scale_x_discrete(expand=c(0,0))
require(ReactomePA)
require(reactome.db)
ls("package:reactome.db")
keytypes(reactome.db)
columns(reactome.db)
keys(reactome.db, keys ="PATHNAME")
AnnotationDbi::select(reactome.db, keys = c("6794"), columns = c("PATHID","PATHNAME"), keytypes="ENTREZID") 
a<- as.list(reactomePATHID2EXTID)$ "R-HSA-111885"
library( "clusterProfiler" )
library( "org.Hs.eg.db" )
df<-bitr(reagenes, fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
write.csv(df,file = '05Reactome_opioid.csv')
df<-read.csv(file = '05Reactome_opioid.csv')
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
DESeq_diff<-DESeq%>%
  dplyr::filter(DESeq$padj<0.05)
deseq<-DESeq_diff[df[,2],]
deseq<-na.omit(deseq)
fisher<-matrix(1:4,nrow = 2,ncol = 2)
fisher<-as.data.frame(fisher)
colnames(fisher)<-c('up','down')
rownames(fisher)<-c('Opioid-related gene','Non-ORG')
ORGnum<-nrow(deseq)
NORGnum<-6353-ORGnum
up<-filter(deseq,log2FoldChange>0)
write.csv(deseq,file = '05Reactome_03Deseq.csv')
opioid<-df[,2]
deseqgene<-rownames(DESeq_diff)
de<-intersect(opioid,deseqgene)
df<-bitr(de, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
rt<-read.table(file = '29eigengene_survival.txt',sep = '\t', header = T,check.names = F)
require(survival)
univariatecox<-matrix(60,nrow = 10,ncol = 6)
colnames(univariatecox)<-c('module','coef','exp(coef)','se(coef)','z','Pr(>|z|)')
for (i in 6:16) {
  cox<-coxph(Surv(rt$os,rt$os_score)~rt[,i],data=rt)
  summary(cox)
}
colnames(rt)[i]
cox<-coxph(Surv(rt$os,rt$os_score)~rt[,i],data=rt)
summary(cox)
colnames(rt)[i]
cox<-coxph(Surv(rt$RFS,rt$RFS_score)~rt[,i],data=rt)
summary(cox)
i<-i+1
rt<-read.table(file = '29eigengene_survival.txt',sep = '\t', header = T,check.names = F)
clinical<-read.table(file = '25num_clinical_WCGNA.txt',sep = '\t',header = T,check.names = F)
colnames(clinical)[1]<-'sample'
require(dplyr)
merge<-inner_join(rt,clinical[,1:9],by='sample')
merge<-merge[,-c(18,19,20)]
rt<-merge
rt<-na.omit(rt)
rt<-rt[,-c(7:13,15,16)]
colnames(rt)[i]
cox<-coxph(Surv(rt$RFS,rt$RFS_score)~rt[,i],data=rt)
summary(cox)
i<-i+1
os.cox<-coxph(Surv(rt$os,rt$os_score)~rt$MEblack+rt$MEblue+rt$MEbrown+rt$MEgreen+rt$MEmagenta+rt$MEpink+rt$MEpurple+rt$MEred+rt$MEturquoise+rt$MEyellow,data = rt)
rfs.cox<-coxph(Surv(rt$RFS,rt$RFS_score)~rt$MEblack+rt$MEblue+rt$MEbrown+rt$MEgreen+rt$MEmagenta+rt$MEpink+rt$MEpurple+rt$MEred+rt$MEturquoise+rt$MEyellow,data = rt)
summary(os.cox)
summary(rfs.cox)
rt<-read.table(file = '29eigengene_survival.txt',sep = '\t', header = T,check.names = F)
clinical<-read.table(file = '25num_clinical_WCGNA.txt',sep = '\t',header = T,check.names = F)
colnames(clinical)[1]<-'sample'
require(dplyr)
merge<-inner_join(rt,clinical[,1:9],by='sample')
merge<-merge[,-c(18,19,20)]
rt<-merge
rt<-na.omit(rt)
os.cox<-coxph(Surv(rt$os,rt$os_score)~rt$MEblack+
                rt$age+
                rt$gleason_score+
                rt$pathologic_T+
                rt$pathologic_N+
                rt$psa_value,data = rt)
summary(os.cox)
rfs.cox<-coxph(Surv(rt$RFS,rt$RFS_score)~rt$MEblack+
                 rt$MEred+
                 rt$age+
                 #rt$gleason_score+
                 rt$pathologic_T+
                 rt$pathologic_N+
                 rt$psa_value,data = rt)
summary(rfs.cox)
cox<-rfs.cox
i<-i+1
install.packages('survminer')
library(survival)
library(survminer)
ggforest(cox, data = rt)
ggforest(cox,  
         data = rt,  
         main = 'Hazard ratio of PRAD',  
         cpositions = c(0.01, 0.15, 0.35), 
         fontsize = 1, 
         refLabel = 'reference',
         noDigits = 3 
)
require(forestplot)
data<-read.csv(file = '33cox_forest.csv')
colnames(data)[1]<-'name'
forestplot(as.matrix(data[,1:2]),
           mean=data$or2,
           lower=data$lower,
           upper=data$upper,
           title='Hazard Ratio',
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           zero=1,
           clip=c(0.4,15),
           graph.pos=2,
           xticks=c(0.1,1,5,1e+4,1e+9),
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           lwd.ci=0.2, boxsize=0.2,
           ci.vertices=TRUE, ci.vertices.height = 0.15,
           #graph.pos=3, #ΪPvalue????ͼ???ڵ?λ??
           xlog=T)
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
gene_exprSet<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header = T,row.names=1,check.names = F)
sampletype<-read.table(file = '03TCGA_PRAD_sampletype.txt',sep = '\t',header = T,check.names = F)
require(dplyr)
gene_exprSet<-gene_exprSet[rownames(DESeq),]
gene_exprSet<-gene_exprSet[,sampletype[,1]]
write.table(gene_exprSet,file = '14DESeq_03GSEA.txt',sep = '\t',col.names = T,row.names = T)
group<-colnames(gene_exprSet)
group<-as.data.frame(group)
for (i in 1:551) {
  ifelse(substr(group[i,1],14,14)=='0',group[i,2]<-'T',group[i,2]<-'N')
}
write.table(group,file = '14DESeq_04GSEAgroup.txt',sep = '\t',col.names = F,row.names = F)
setwd('D:/Paper/opioid/NT_Deseq')
DESeq<-read.table(file = '13DESeq_PRAD_NT.txt',sep = '\t', header = T, check.names = F)
require(dplyr)
DESeq<-filter(DESeq,padj<0.05&abs(log2FoldChange)>log2FC)
gene<-cbind(rownames(DESeq),DESeq[,2])
colnames(gene)<-c('symbol','log2FC')
gene<-as.data.frame(gene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
geneList<-as.numeric(gene[,2])
names(geneList)<-gene[,1]
geneList<-sort(geneList,decreasing = T)
kegmt<-read.gmt("c2.cp.kegg.v7.4.symbols.gmt") #canonical pathway_KEGG
kegmt<-read.gmt("c2.cp.reactome.v7.4.symbols.gmt") #canonical pathway_Reactome
kegmt<-read.gmt("c5.go.v7.4.symbols.gmt") #GO all
kegmt<-read.gmt("21REACTOME_OPIOID_SIGNALLING.gmt") #Reactome opioid signalling
KEGG<-GSEA(geneList,
           pvalueCutoff = 1,
           seed = F,
           TERM2GENE = kegmt) 
GSEAresult<-KEGG@result
ReactomeGSEAresult<-GSEAresult
write.csv(ReactomeGSEAresult,file = '15GSEA_01opioid.csv')
write.csv(GSEAresult,file = '15GSEA_02KEGG.csv')
library(ggplot2)
dotplot(KEGG)
dotplot(KEGG,color="pvalue")  
dotplot(KEGG,split=".sign")+facet_grid(~.sign) 
dotplot(KEGG,split=".sign")+facet_wrap(~.sign,scales = "free") 
library(enrichplot)
gseaplot2(KEGG,1,color="slateblue2",pvalue_table = T) 
gseaplot2(KEGG,1:10,color="red") 
reactome<-read.csv(file = '05Reactome_opioid.csv')
clinical<-read.table(file = '25num_clinical_WCGNA.txt',sep = '\t',header = T,check.names = F)
gene_exprSet<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header = T,row.names = 1,check.names = F)
reactome<-reactome[,1]
rt<-gene_exprSet[reactome,]
rt<-rbind(substr(colnames(gene_exprSet),1,15),rt)
rt<-t(rt)
colnames(rt)[1]<-'sample'
colnames(clinical)[1]<-'sample'
rt<-as.data.frame(rt)
require(dplyr)
merge<-inner_join(clinical,rt,by='sample')
merge<-merge[!duplicated(merge[,1]),]
rownames(merge)<-merge[,1]
merge<-merge[,-c(1:5,7:13)]
merge<-na.omit(merge)
dimnames<-list(rownames(merge),colnames(merge))
merge<-matrix(as.numeric(as.matrix(merge)),nrow=nrow(merge),dimnames=dimnames)
for (i in 1:nrow(merge)) {
  ifelse(merge[i,1]<=7,merge[i,1]<-1,merge[i,1]<-2)
}
merge<-as.data.frame(merge)
require(dplyr)
group1<-filter(merge,gleason_score==1)
group2<-filter(merge,gleason_score==2)
ttest<-matrix(ncol = 2,nrow = ncol(merge)-1)
for (i in 2:ncol(merge)) {
  a<-group1[,i]
  b<-group2[,i]
  c<-t.test(a,b,paired=F)
  ttest[i-1,1]<-colnames(merge)[i]
  ttest[i-1,2]<-c$p.value
}
ttest<-as.data.frame(ttest)
colnames(ttest)<-c('gene','Pvalue')
write.csv(ttest,file = '05Reactome_08gleason_t.csv')
write.csv(merge,file = '05Reactome_08gleason_exp91.csv')
ttest<-ttest[order(ttest[,2],decreasing = F),]
gene<-ttest[1:12,1]
gene[13]<-'gleason_score'
merge<-merge[,gene]
write.csv(merge,file = '05Reactome_08gleason_exp91.csv')

#ANOVA
anova<-matrix(nrow = ncol(merge)-1,ncol = 2)
colnames(anova)<-c('gene','Pvalue')
for (i in 2:ncol(merge)) {
  b<-oneway.test(merge[,i]~factor(merge[,1]))
  anova[i-1,1]<-colnames(merge)[i]
  anova[i-1,2]<-b$p.value
}
anova[,2]<-as.numeric(anova[,2])
anova<-as.data.frame(anova)
write.csv(anova,file = '05Reactome_04Tstage_anova.csv')
write.csv(merge,file = '05Reactome_05Tstage_exp91.csv')
rt<-filter(anova,Pvalue<0.01)
Tstage_filt<-as.data.frame(cbind(merge[,1],merge[,rt[,1]]))
colnames(Tstage_filt)[1]<-'pT stage'
write.csv(Tstage_filt,file = '05Reactome_06T_exp9.csv')
rt<-read.csv(file = '05Reactome_06T_exp9.csv')
require(dplyr)
table(rt[,2])
rt<-rt[order(rt[,2],decreasing = T),]
colnames(rt)[2]<-'pT'
rt2<-filter(rt,pT=='2')
rt3<-filter(rt,pT=='3')
rt4<-filter(rt,pT=='4')
write.csv(rt2,file = '05Reactome_07T2_exp.csv')
write.csv(rt3,file = '05Reactome_07T3_exp.csv')
write.csv(rt4,file = '05Reactome_07T4_exp.csv')
require(WGCNA)
options(stringsAsFactors = F)
dat1<-read.table(file = '23gene_exp_WGCNA.txt', sep = '\t',header = T,row.names = 1,check.names = F)
dim(dat1)
names(dat1)
datExpr<-dat1
indextrain=c(245:489)
indextest=c(1:244)
nSets = 2
multiExpr = list()
multiExpr[[1]] = list(data = datExpr[indextrain, ])
multiExpr[[2]] = list(data = datExpr[indextest, ])
setLabels = c("Human", "Test")
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)
x = load("AS-green-FPKM-TOM-block.1.RData") 
x = load('17WGCNA.Rdata')
colorTrain<-net$colors
colorTest<-net$colors
colorList = list(colorTrain, colorTest)
names(colorList) = setLabels
system.time( {
  mp = modulePreservation(multiExpr, colorList,
                          referenceNetworks = c(1:2),
                          loadPermutedStatistics = FALSE,
                          nPermutations = 150,
                          verbose = 3)
} )
save(mp, file = "17WGCNA_02modulePreservation.RData")
ref = 1 
test = 2 
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2))
a<-signif(statsZ[, "Zsummary.pres", drop = FALSE],2)
a<-a[-c(1,2),]
a<-as.data.frame(a)
rownames(a)<-c('turquoise','purple','blue','brown','yellow','green','red','black','pink','magenta')
genecolor<-read.table(file = '27module_gene.txt',sep = '\t',header = T,check.names = F)
table(genecolor[,2])
a[,2]<-c(1576,31,446,425,298,283,197,160,157,46)
colnames(a)<-c('Z statistic','Module Size')
require(ggplot2)
plot(log2(a$`Module Size`), a$`Z statistic`, xlab="log2(Module Size)",ylab="Z Statistic",
     xlim=c(0,12),ylim=c(-20,60),yaxs="i",pch=20, cex=2)
for (i in 1:10) {
  b<-as.data.frame(a[i,])
  points(log2(b[,2]),b[,1],pch=20,cex=2,col=rownames(a)[i])
}
abline(h=0,lwd=1)
abline(h=2,lty=2,lwd=1,col='orangered')
abline(h=10,lty=2,lwd=1,col='darkblue')
write.csv(a,file = '17WGCNA_Zstatistics.csv')

#17
setwd('D:/Paper/opioid/NT_Deseq')
MSKCC<-read.table(file = '01MSKCC_PRAD_norm.txt',sep = '\t',header = T,row.names = 1,check.names = F)
TCGA<-read.table(file = '12gene_exprSet_PRAD.txt',sep = '\t',header=T,row.names=1,check.names=F)
MSKCC<-na.omit(MSKCC)
rownames(MSKCC)<-MSKCC[,1]
MSKCC<-MSKCC[,-1]
gene_exprSet<-MSKCC[rownames(TCGA),]
require(WGCNA)
require(impute)
require(reshape2)
require(stringr)
type<-'signed'
corType<-'pearson' #
corFnc<-ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
dimnames<-list(rownames(gene_exprSet),colnames(gene_exprSet))
gene_exprSet<-matrix(as.numeric(as.matrix(gene_exprSet)),nrow=nrow(gene_exprSet),dimnames=dimnames)
gene_exprSet<-gene_exprSet[rowMeans(gene_exprSet)>1,]
gene_exprSet<-na.omit(gene_exprSet)
require(WGCNA)
require(stringr)
require(dplyr)
load(file = '17WGCNA_01MSKCC.Rdata')
moduleColors <- labels2colors(net$colors)
ADJ1<-abs(cor(gene_exprSet,use="p"))^6
Alldegrees1<-intramodularConnectivity(ADJ1, moduleColors)
MEs_col<-net$MEs
colnames(MEs_col)<-paste0('ME',labels2colors(
  as.numeric(str_replace_all(colnames(net$MEs),'ME',''))
))
MEs_col<-orderMEs(MEs_col)
datKME<-signedKME(gene_exprSet, MEs_col, outputColumnName="MM.")
write.table(datKME,file = "09hubgene_05MSKCCMM.txt",sep = '\t',col.names = T,row.names = T)

load(file = '17WGCNA_01MSKCC.Rdata')
moduleLabelsAutomatic<-net$colors
moduleLabelsAutomatic<-labels2colors(moduleLabelsAutomatic)
moduleColorsFemale<-moduleLabelsAutomatic
dataExp<-gene_exprSet
t_rt<-dataExp
nGenes<-ncol(t_rt)
nSample<-nrow(t_rt)
blocknumber<-1
MEs0<-moduleEigengenes(t_rt,moduleColorsFemale)$eigengenes
MEsFemales<-orderMEs(MEs0)
datColor<-data.frame(moduleLabelsAutomatic)[net$blockGenes[[blocknumber]],]
datKME<-signedKME(t_rt,MEsFemales)
datSummary<-colnames(dataExp)
datout<-data.frame(datSummary,datColor,datKME)
write.csv(datout,file = '09hubgene_06MSKCCmodulegene.csv')

MMMSKCC<-read.table(file = '09hubgene_05MSKCCMM.txt',sep = '\t',header = T,check.names = F)
moduleMSKCC<-read.csv(file = '09hubgene_06MSKCCmodulegene.csv')
geneMSKCC<-moduleMSKCC[,c(1,3)]
rownames(geneMSKCC)<-geneMSKCC[,1]
geneMSKCC<-geneMSKCC[rownames(MMMSKCC),]
require(dplyr)
rtmskcc<-cbind(geneMSKCC,MMMSKCC)
rtmskcc<-rtmskcc[,-1]
MMTCGA<-read.table(file = '09hubgene_01MM.txt',sep = '\t',header = T,check.names = F)
moduleTCGA<-read.table(file = '27module_gene.txt',sep = '\t',header = T,check.names = F)
geneTCGA<-moduleTCGA[,c(1,2)]
rownames(geneTCGA)<-geneTCGA[,1]
geneTCGA<-geneTCGA[rownames(MMTCGA),]
require(dplyr)
rtTCGA<-cbind(geneTCGA,MMTCGA)
rtTCGA<-rtTCGA[,-1]


require(dplyr)
a<-c('turquoise','purple','blue','brown','yellow','green','red','black','pink','magenta')
b<-c('black','blue','yellow','green','red','turquoise','brown')
mergeratio<-matrix(ncol = 7,nrow = 10)
for (i in 2:11) {
  colorTCGA<-a[i]
  rt<-filter(rtTCGA,datColor==colorTCGA)
  gene1<-rownames(rt)
  MMcolorTCGA<-paste('MM.',colorTCGA,sep = '')
  genelist1<-rt[,MMcolorTCGA]
  for (x in 2:8) {
    genelist2<-rtmskcc[gene1,x]
    verboseScatterplot(x=abs(genelist1),
                       y=abs(genelist2),
                       ylab = paste('Module Membership in',colnames(rtmskcc)[x],'MSKCC'),
                       xlab = paste('Module Membership in',colorTCGA),
                       main = paste('TCGA module',colorTCGA),
                       cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = colorTCGA
    )
  }}



for (i in 2:11) {
  colorTCGA<-a[i]
  gene1<-filter(rtTCGA,datColor==colorTCGA)
  B<-paste('MM.',colorTCGA,sep = '')
  CTCGA<-as.data.frame(gene1[,B])
}


table(net$colors)
module<-'red'
pheno<-'RFS_score'
modNames<-substring(colnames(MEs_col),3)
moduleLabels<-net$colors
moduleColors<-labels2colors(moduleLabels)
module_column<-match(module,modNames)
pheno_column<-match(pheno,colnames(dclic))
moduleGenes<-moduleColors==module
sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,module_column]),
                   abs(geneTraitCor[moduleGenes,pheno_column]),
                   xlab = paste('Module Membership in',module,'module'),
                   ylab = paste('Gene significance for',pheno),
                   main = paste('Module Membership vs. Gene significant\n'),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module
)
load(file = '17WGCNA.Rdata')
gene<-c('TROAP','EME1','CCNL2','CHTF18')
gene_exprSet<-gene_exprSet[,gene]
eigengene<-gene_exprSet
survival<-read.table(file = '22survival_chen.txt',sep = '\t',header = T,check.names = F)
k<-substr(rownames(eigengene),1,15)
eigengene<-cbind(k,eigengene)
eigengene<-as.data.frame(eigengene)
eigengene<-eigengene[!duplicated(eigengene$k),]
colnames(eigengene)[1]<-'sample'
require(dplyr)
eigengene_survival<-dplyr::inner_join(survival,eigengene,by='sample')
require(survival)
rt<-eigengene_survival
geneset<-colnames(rt)[6:9]
pvalueset<-as.data.frame(colnames(rt)[6:9])
for (i in 6:9) {
  gene<-colnames(rt)[i]
  time<-'os'
  status<-'os_score'
  a<-rt[,gene]<median(rt[,gene])
  time<-as.numeric(rt[,time])
  status<-as.numeric(rt[,status])
  diff<-survdiff(Surv(time, status) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(time, status) ~ a, data = rt)
  summary(fit)    #?鿴??????????
  plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
       main=paste(gene,'os','(p=', pValue ,")",sep=" "))
  pvalueset[(i-5),2]<-pValue
}
geneset<-colnames(rt)[6:9]
for (i in 6:9) {
  gene<-colnames(rt)[i]
  time<-'RFS'
  status<-'RFS_score'
  a<-rt[,gene]<median(rt[,gene])
  time<-as.numeric(rt[,time])
  status<-as.numeric(rt[,status])
  diff<-survdiff(Surv(time, status) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=round(pValue,5)
  fit <- survfit(Surv(time, status) ~ a, data = rt)
  summary(fit)    
  plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
       main=paste(gene,'RFS','(p=', pValue ,")",sep=" "))
  pvalueset[(i-5),3]<-pValue
}
colnames(pvalueset)<-c('gene','OS_Pvalue','RFS_Pvalue')
pvalueset
write.table(pvalueset,file = '31survival_02fourgeneinR_Bmodule.txt',sep = '\t',col.names = T,row.names = F)
write.csv(eigengene_survival,file = '31survival_021fourgeneinR_Bmodule.csv')


#24
setwd('D:/Paper/opioid/NT_Deseq')
genecolor<-read.table(file = '27module_gene.txt',sep = '\t',header = T,check.names = F)
require(dplyr)
red<-filter(genecolor,datColor=='red')
black<-filter(genecolor,datColor=='black')
driver<-read.table(file = '36IntOGen-DriverGenes_PRAD.tsv',sep = '\t',header = T,check.names = F)
colnames(driver)[1]<-'datSummary'
blackdriver<-inner_join(black,driver,by='datSummary')
reddriver<-inner_join(red,driver,by='datSummary')





#28
setwd('D:/Paper/opioid/NT_Deseq')
module<-'black'
module<-paste('ME',module,sep = '')
filename<-paste('38cibersort_group',module,'.csv',sep = '')
sample<-read.csv(file = filename)
sample<-sample[order(sample[,3],decreasing = T),]
rt<-sample[,2]
gene_exprSet<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header = T,row.names = 1,check.names = F)
gene_exprSet<-gene_exprSet[,rt]
write.table(gene_exprSet,file = '38matrix_cibersort.txt',sep = '\t',row.names = T,col.names = T)
source("CIBERSORT.R")
results=CIBERSORT("38cibersort_ref.txt", "38matrix_cibersort.txt", perm=200, QN=TRUE)  #ģ??200?Σ???Pֵ
Results <- read.table("CIBERSORT-Results.txt",sep="\t",header=T,check.names=F)
Results <- read.table("38CIBERSORT.filterGDC.txt",sep="\t",header=T,check.names=F)

#Results <- as.data.frame(Results)
require(dplyr)
require(tidyr)
Results <- Results %>% dplyr::filter(P.value <= 0.05) 
colnames(Results)[1]<-'sample'
Results<-inner_join(Results,sample,by='sample')
Resultsorder<-Results[order(Results[,28],decreasing = FALSE),]
Cibersort.filter <- Resultsorder[,-(24:28)]
save(Cibersort.filter,file = '38CIBERSORT.filter.Rda')
write.table(Cibersort.filter, file ="38CIBERSORT.filter.txt", sep ="\t", row.names =TRUE, col.names =TRUE, quote =TRUE)

#barplot
load(file = '38CIBERSORT.filter.Rda')
data <- Cibersort.filter
data <- t(data)
colnames(data)<-data[1,]
data <- data[-1,]
outpdf<-"barplot.pdf"
col=rainbow(nrow(data),s=0.7,v=0.7)
pdf(outpdf,height=10,width=15)
par(las=1,mar=c(10,4,4,15))
a1 <- barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2 <- axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3) 
dev.off() 

#vioplot
#install.packages("vioplot")
library(vioplot)    
high<-82                
low<-85                                                  
rt <-read.table("38CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F) 
rt<-rt[,-1]
pdf("vioplot.pdf",height=8,width=15)      
par(las=1,mar=c(10,6,3,3))
x<-c(1:ncol(rt))
y<-c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(rt)){
  normalData=rt[1:high,i]
  tumorData=rt[(high+1):(high+low),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'orangered')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'seagreen')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
dev.off()


eigengene<-read.table(file = '28WGCNA_eigengene.txt',sep = '\t',header = T,row.names = 1,check.names = F)
module<-'red'
module<-paste('ME',module,sep = '')
group<-matrix(nrow = 489,ncol = 2)
group[,1]<-rownames(eigengene)
group[,2]<-'ns'
for (i in 1:489) {
  if(eigengene[i,module]>=quantile(eigengene[,module],0.75)){group[i,2]<-'high'}
  if(eigengene[i,module]<=quantile(eigengene[,module],0.25)){group[i,2]<-'low'}
}
colnames(group)<-c('sample',module)
filename<-paste('38cibersort_14group',module,'.csv',sep = '')
write.csv(group,file = filename)

module<-'red'
module<-paste('ME',module,sep = '')
filename<-paste('38cibersort_14group',module,'.csv',sep = '')
sample<-read.csv(file = filename)
sample<-sample[order(sample[,3],decreasing = T),]
rt<-sample[,2]
Results <- read.table("38CIBERSORT.filterGDC.txt",sep="\t",header=T,check.names=F)

require(dplyr)
require(tidyr)
Results <- Results %>% dplyr::filter(P.value <= 0.05) 
colnames(Results)[1]<-'sample'
Results<-inner_join(Results,sample,by='sample')
Resultsorder<-Results[order(Results[,28],decreasing = FALSE),]
Resultsorder<-filter(Resultsorder,Resultsorder[,28]!='ns')
Cibersort.filter <- Resultsorder[,-(24:28)]
save(Cibersort.filter,file = '38CIBERSORT.filter.Rda')
write.table(Cibersort.filter, file ="38CIBERSORT.filter.txt", sep ="\t", row.names =TRUE, col.names =TRUE, quote =TRUE)

#barplot
load(file = '38CIBERSORT.filter.Rda')
data <- Cibersort.filter
data <- t(data)
colnames(data)<-data[1,]
data <- data[-1,]
outpdf<-"barplot.pdf"
col=rainbow(nrow(data),s=0.7,v=0.7)
pdf(outpdf,height=10,width=15)
par(las=1,mar=c(10,4,4,15))
a1 <- barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2 <- axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(data[,ncol(data)]) 
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3) 
dev.off() 

#vioplot
#install.packages("vioplot")
library(vioplot)                    
table(Resultsorder[,28])
high<-44       
low<-36  
rt <-read.table("38CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)
rt<-rt[,-1]
pdf("vioplot.pdf",height=8,width=15) 
par(las=1,mar=c(10,6,3,3))
x<-c(1:ncol(rt))
y<-c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(rt)){
  normalData=rt[1:high,i]
  tumorData=rt[(high+1):(high+low),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'orangered')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'seagreen')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
dev.off()



#31
setwd('D:/Paper/opioid/NT_Deseq')
eigengene<-read.table(file = '28WGCNA_eigengene.txt',sep = '\t',header = T,row.names = 1,check.names = F)
Results <- read.table("38CIBERSORT.filterGDC.txt",sep="\t",header=T,check.names=F)
eigengene<-cbind(rownames(eigengene),eigengene)
colnames(eigengene)[1]<-'SampleID'
require(dplyr)
merge<-inner_join(eigengene[,c(1,2,10)],Results,by="SampleID")
library(ggplot2)
require(ggpubr)
merge<-merge[,-(26:28)]
for (i in 4:ncol(merge)) {
  ggplot(merge,aes_string(x='MEred',y=paste(colnames(merge)[i])))+
    geom_point(color='red')+
    stat_smooth(method = 'lm',se=F,color='red')+
    stat_cor(data = merge,method = 'pearson')
  i<-i+1
}
write.csv(merge,file = '38CIBERSORT_eigengene.csv')



genecolor<-read.table(file = '27module_gene.txt',sep = '\t',header = T,check.names = F)
exp<-read.table(file = '01TCGA_PRAD_norm.txt',sep = '\t',header = T,check.names = F)
require(dplyr)
black<-filter(genecolor,datColor=='black')
rownames(exp)<-exp[,1]
exp<-exp[,-1]
exp<-exp[black[,1],]
survival<-read.table(file = '22survival_chen.txt',sep = '\t',header = T,check.names = F)
rt<-t(exp)
rt<-cbind(rownames(rt),rt)
rt[,1]<-substr(rt[,1],1,15)
rt<-rt[!duplicated(rt[,1]),]
rownames(rt)<-rt[,1]
rt<-rt[,-1]
a<-intersect(rownames(rt),survival[,1])
rt<-cbind(rownames(rt),rt)
colnames(rt)[1]<-'sample'
rt<-as.data.frame(rt)
merge<-inner_join(survival,rt,by="sample")
write.csv(merge,file = '39lasso_surv_black.csv')
require(survival)
pFilter<-0.01
merge<-read.csv(file = '39lasso_surv_black.csv')
rownames(merge)<-merge[,2]
rt<-merge[,-c(1,2,3,4)]#
rt<-na.omit(rt)
sigGenes<-c('os_score','os')#
sigGenes<-c('RFS_score','RFS')#
outTab<-data.frame()
for (i in 3:ncol(rt)) {
  cox<-coxph(Surv(os,os_score)~rt[,i],data = rt)
  coxSummary<-summary(cox)
  coxP<-coxSummary$coefficients[,'Pr(>|z|)']
  outTab<-rbind(outTab,
                cbind(id=colnames(rt)[i],
                      HR=coxSummary$conf.int[,'exp(coef)'],
                      HR.95L=coxSummary$conf.int[,'lower .95'],
                      HR.95H=coxSummary$conf.int[,'upper .95'],
                      pvalue=coxSummary$coefficients[,'Pr(>|z|)']))
  if(coxP<pFilter){sigGenes<-c(sigGenes,colnames(rt)[i])}
}
write.csv(outTab,file = '39lasso_cox_os_black.csv')
uniExp<-rt[,sigGenes]
write.csv(uniExp,file = '39lasso_cox_os_black_02exp.csv')
require(glmnet)
require(survival)
rt<-read.csv(file = '39lasso_cox_os_black_02exp.csv',row.names = 1)
min(rt[,2])
x<-as.matrix(rt[,c(3:ncol(rt))])
y<-as.matrix(Surv(rt$os,rt$os_score))
fit<-glmnet(x,y,family = 'cox',maxit = 1500)
plot(fit,xvar = 'lambda',label = T)
cvfit<-cv.glmnet(x,y,family='cox',maxit=1500)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty='dashed')
coef<-coef(fit,s=cvfit$lambda.min)
index<-which(coef!=0)
actCoef<-coef[index]
lassoGene<-row.names(coef)[index]
lassoSigExp<-rt[,lassoGene]
lassoSigExp<-cbind(id=rownames(lassoSigExp),lassoSigExp)
write.csv(lassoSigExp,file = '39lasso_black_geneexp.csv')
genecolor[lassoGene,]

require(survival)
pFilter<-0.01
merge<-read.csv(file = '39lasso_surv_black.csv')
rownames(merge)<-merge[,2]
rt<-merge[,-c(1,2,3,4)]#
rt<-na.omit(rt)
sigGenes<-c('RFS_score','RFS')#
outTab<-data.frame()
for (i in 3:ncol(rt)) {
  cox<-coxph(Surv(RFS,RFS_score)~rt[,i],data = rt)
  coxSummary<-summary(cox)
  coxP<-coxSummary$coefficients[,'Pr(>|z|)']
  outTab<-rbind(outTab,
                cbind(id=colnames(rt)[i],
                      HR=coxSummary$conf.int[,'exp(coef)'],
                      HR.95L=coxSummary$conf.int[,'lower .95'],
                      HR.95H=coxSummary$conf.int[,'upper .95'],
                      pvalue=coxSummary$coefficients[,'Pr(>|z|)']))
  if(coxP<pFilter){sigGenes<-c(sigGenes,colnames(rt)[i])}
}
write.csv(outTab,file = '39lasso_cox_RFS_black.csv')
uniExp<-rt[,sigGenes]
write.csv(uniExp,file = '39lasso_cox_RFS_black_02exp.csv')

geneexp_MSKCC<-read.csv(file = '40MSKCC_exp.csv',row.names = 1)
clinical_MSKCC<-read.csv(file = '40MSKCC_clinical.csv',row.names = 1)
geneexp_MSKCC<-t(geneexp_MSKCC)
require(dplyr)
clinical_MSKCC<-filter(clinical_MSKCC,SAMPLE_CLASS=='Tumor')
geneexp_MSKCC<-cbind(rownames(geneexp_MSKCC),geneexp_MSKCC)
colnames(geneexp_MSKCC)[1]<-'ID'
clinical_MSKCC<-clinical_MSKCC[,c(1,7,8)]
geneexp_MSKCC<-as.data.frame(geneexp_MSKCC)
merge<-inner_join(clinical_MSKCC,geneexp_MSKCC,by='ID')
write.csv(merge,file = '40MSKCC_merge.csv',row.names = F)
#coxģ?͹???
require(survival)
require(survminer)
rt<-read.csv(file = '39lasso_black_RFS_geneexp.csv',row.names = 1)
multicox<-coxph(Surv(RFS,RFS_score)~.,data = rt)#.???????л???????��??Ϊx
multicox<-step(multicox,direction = 'both')
multicoxsum<-summary(multicox)
outTab<-data.frame()
outTab<-cbind(coef=multicoxsum$coefficients[,'coef'],
              HR=multicoxsum$conf.int[,'exp(coef)'],
              HR.95L=multicoxsum$conf.int[,'lower .95'],
              HR.95H=multicoxsum$conf.int[,'upper .95'],
              pvalue=multicoxsum$coefficients[,'Pr(>|z|)']
)
outTab<-cbind(id=rownames(outTab),outTab)
write.csv(outTab,file = '39lasso_multiCox_blackRFS.csv')
ggforest(multicox,main = 'Harzard Ratio',cpositions = c(0.02,0.22,0.4),fontsize = 0.7,refLabel = 'reference',noDigits = 2)
riskscore<-predict(multicox,type = 'risk',newdata = rt)
coxGene<-rownames(multicoxsum$coefficients)
coxGene<-gsub('`','',coxGene)
outCol<-c("RFS","RFS_score",coxGene)
medianTrainRisk=median(riskscore)
risk=as.vector(ifelse(riskscore>medianTrainRisk,"high","low"))

                
colnames(rtTest)[1:2]<-c('RFS','RFS_score')


rtTest<-na.omit(rtTest)
rtTest<-cbind(rtTest[,colnames(rt)])

riskScoreTest=predict(multicox,type="risk",newdata=rtTest)      #????train?õ?ģ??Ԥ??test??Ʒ????
medianTestRisk<-median(riskScoreTest)
riskTest=as.vector(ifelse(riskScoreTest>medianTestRisk,"high","low"))
riskTest<-cbind(id=rownames(cbind(rtTest[,outCol],riskScoreTest,rtTest[,outCol],riskScore=riskScoreTest,risk=riskTest)))
write.csv(riskTest,file = '39lasso_risk_test.csv',row.names = F)
library(survival)
rt<-read.csv("39lasso_risk_train.csv",row.names = 1)
diff<-survdiff(Surv(RFS, RFS_score) ~risk,data = rt)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(RFS, RFS_score) ~ risk, data = rt)
summary(fit)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
rt=read.csv("39lasso_risk_test.csv",row.names = 1)
diff=survdiff(Surv(RFS, RFS_score) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(RFS, RFS_score) ~ risk, data = rt)
summary(fit)    
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
require(pROC)
rt<-read.csv("39lasso_risk_train.csv",row.names = 1)
roc(rt$RFS_score,rt$riskscore,plot=T)
rt=read.csv("39lasso_risk_test.csv",row.names = 1)
roc(rt$RFS_score,rt$riskScore,plot=T)
require(pheatmap)
rt<-read.csv("39lasso_risk_train.csv",row.names = 1)
rt=rt[order(rt$riskscore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
rt1=log2(rt1+0.01)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         #fontsize_row=11,
         #fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
rt<-read.csv("39lasso_risk_test.csv",row.names = 1)
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
rt1=log2(rt1+0.01)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         #fontsize_row=11,
         #fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
rt<-read.csv("39lasso_risk_train.csv",row.names = 1)  
rt=rt[order(rt$riskscore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskscore"]
line[line>10]=10
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
rt<-read.csv("39lasso_risk_test.csv",row.names = 1)   
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=trainMedianScore,v=lowLength,lty=2)
