#############################
library(stringr)
library(cowplot)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(data.table)
#######################################
args<-commandArgs(T)
if(length(args)!=5){
  print ("Usage:")
  print ("cd work_dir ;Rcript anno.R count_matrix sample max_nFeature_RNA doublet_rate")
  print ("Note:")
  print ("Aim to run single-sample umap")
  q(save = "no", status = 0, runLast = TRUE)
}
setwd(args[1])
mtx <- Read10X(args[2])
object.data = mtx
sample<-args[3]
colnames(object.data)<-paste(colnames(object.data),sample,sep="_")
object_name <- CreateSeuratObject(counts = object.data, project =sample, min.cells = 3, min.features = 200)
a=length(colnames(object_name))
print(sample)
print(a)
#####################################
object_name[["percent.mt"]] <- PercentageFeatureSet(object = object_name, pattern = "^MT-|^mt-")
if (sum(object_name@meta.data$percent.mt)== 0){
  object_name@meta.data$percent.mt[1]<-0.000001
}
object_name <- subset(x = object_name, subset = nFeature_RNA > 800 & nFeature_RNA < as.numeric(args[4]))
a=length(colnames(object_name))
print(a)
object_name <- subset(x = object_name, subset = percent.mt < 10)
a=length(colnames(object_name))
print(a)
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes
object_name=NormalizeData(object_name)
object_name <- FindVariableFeatures(object_name,  selection.method = "vst",nfeatures = 2000,  verbose = FALSE)
object_name <- ScaleData(object_name)
object_name=CellCycleScoring(object_name, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
object_name <- RunPCA(object_name)
pdf("cellcycle.pdf")
DimPlot(object_name,reduction = "pca",group.by= "Phase",split.by = "Phase")   
dev.off()
object_name<- SCTransform(object_name, verbose = FALSE,vars.to.regress=c("percent.mt","S.Score", "G2M.Score"))
object_name <- RunPCA(object_name, features = VariableFeatures(object = object_name))
object_name <- FindNeighbors(object_name, dims = 1:30)
object_name <- FindClusters(object_name, resolution = 0.5)
object_name <- RunUMAP(object_name, dims = 1:30)
sweep.res.list_SM <- paramSweep_v3(object_name, PCs = 1:30,sct=TRUE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- object_name@meta.data$RNA_snn_res.1
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(as.numeric(args[5])*length(object_name@meta.data$orig.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE,sct=TRUE)
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value,sct=TRUE)
object_name@meta.data[,"DF_hi.lo"] <- object_name@meta.data[,10]
object_name@meta.data$DF_hi.lo[which(object_name@meta.data$DF_hi.lo == "Doublet" & object_name@meta.data[,11] == "Singlet")] <- "Doublet_lo"
object_name@meta.data$DF_hi.lo[which(object_name@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
object_name@meta.data$Doublet <-eval(parse(text = paste0("object_name@meta.data$DF.classifications_",pN_value,"_",pK_value,'_',nExp_poi))) 
saveRDS(object_name,"object_name.doublet.rds")
object_Singlet <- SubsetData(object_name,subset.name='Doublet',accept.value='Singlet')
a=length(colnames(object_name))
print(a)
combined_analysis<-a
###################注意去重后，需要重新PCA找高可变基因###
P1=readRDS("/dellfsqd2/P1/combined_analysis.rds")
P1$stim="P1"
P3=readRDS("/dellfsqd2/P2/combined_analysis.rds")
P3$stim="P2"
P3=readRDS("/dellfsqd2/P3/combined_analysis.rds")
P3$stim="P3"
P4=readRDS("/dellfsqd2/P4/combined_analysis.rds")
P4$stim="P4"
P5=readRDS("/dellfsqd2/P5/combined_analysis.rds")
P5$stim="P5"
P6=readRDS("/dellfsqd2/P6/combined_analysis.rds")
P6$stim="P6"
con1=readRDS("/dellfsqd2/con1/combined_analysis.rds")
con1$stim="con1"
con2=readRDS("/dellfsqd2/con2/combined_analysis.rds")
con2$stim="con2"
con3=readRDS("/dellfsqd2/con3/combined_analysis.rds")
con3$stim="con3"
con4=readRDS("/dellfsqd2/con4/combined_analysis.rds")
con4$stim="con4"
con5=readRDS("/dellfsqd2/con5/combined_analysis.rds")
con5$stim="con5"
con6=readRDS("/dellfsqd2/con6/combined_analysis.rds")
con6$stim="con6"
con7=readRDS("/dellfsqd2/con7/combined_analysis.rds")
con7$stim="con7"
con8=readRDS("/dellfsqd2/con8/combined_analysis.rds")
con8$stim="con8"
con9=readRDS("/dellfsqd2/con9/combined_analysis.rds")
con9$stim="con9"
con10=readRDS("/dellfsqd2/con10/combined_analysis.rds")
con10$stim="con10"
con11=readRDS("/dellfsqd2/con11/combined_analysis.rds")
con11$stim="con11"
con12=readRDS("/dellfsqd2/con12/combined_analysis.rds")
con12$stim="con12"
#################################3
bm40k.list=list(con1,con2,con3,con4,con5,con6,con7,con8,con9,con10,con11,con12,P6,P5,P4,P3,P2,P1)
features <- SelectIntegrationFeatures(object.list = bm40k.list)
bm40k.list <- lapply(X = bm40k.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE,vars.to.regress=c("percent.mt","S.Score", "G2M.Score"))
  x <- RunPCA(x, features = features, verbose = FALSE)
})
print("FindIntegrationAnchors")
anchors <- FindIntegrationAnchors(object.list = bm40k.list, reduction = "rpca", dims = 1:20,reference = c(11,12))
print("IntegrateData")
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
print ("aa")
DefaultAssay(object = combined) <- "integrated"
saveRDS(combined, file = "combined.rds")
combined <- ScaleData(object = combined, verbose = FALSE,vars.to.regress=c("percent.mt","S.Score", "G2M.Score"))
print("pca")
pdf("Vlnplot.pdf",width=15,height=10)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, group.by = "stim")
dev.off()
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
pdf("DimHeatmap_30pc.pdf",width=15,height=10)
DimHeatmap(object = combined, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.6)
pdf("align_Clusters.pdf",width=15,height=10)
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
pdf("dimplot_sample.pdf",width=15,height=10)
DimPlot(object = combined, reduction = "umap", split.by = "stim")
dev.off()
pdf("cluster.pdf",width=13,height=12)
DimPlot(object = combined, reduction = "umap",label = TRUE)
dev.off()
saveRDS(combined, file = "combined_analysis.rds")
combined.markers <- FindAllMarkers(object = combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use ="wilcox")
topn <- combined.markers %>% group_by(cluster)
topn <- select(topn, gene, everything())
topn$gene=gsub("-","_",topn$gene)
write.table(topn[topn$p_val_adj <= 0.05,],file="GeneDiffExpFilter.xls",sep="\t",col.names = TRUE,row.names = F,quote=F)
d<-combined
##########################################################Figure2###
feature1<-c("CD3D","CD3G","CD3E","CD8A","CD8B","CD4",
            "SPON2","GNLY","GZMB",
            "CD79A","MS4A1","CD79B",
            "LYZ","CD14","S100A8","S100A9",
            "PPBP",
            "JCHAIN","IRF8",
            "GATA2",
            "SPINK2",
            "NEAT1")
d@meta.data$celltype2<-factor(d@meta.data$celltype2,levels = c("CD4_T","CD8_T","gd_T",
                                                               "NK_T", "B_cell","monocyte",
                                                               "Platete","plasma_cell",
                                                               "HSC","NEAT1_cell","SPINK2_cell"))
pdf("heatmap_celltype2.pdf",width = 30,height = 20)
DoHeatmap(combined, feature =feature1, group.by = "celltype2")+scale_fill_gradientn(colors = c("navy","white","firebrick3"))
dev.off()
col5<-c("#5050FFFF",
        "#BA6338FF",
        "#D60047FF",
        "#802268FF",
        "#749B58FF",
        "#FFC20AFF",
        "#7A65A5FF","#FF1463FF","#5A655EFF",
        "#A9A9A9FF","#837B8DFF")
pdf("umap_celltype1_col_major.pdf",width =6,height = 4)
DimPlot(object =d, pt.size=0.0001,label = F, cols=col5)
dev.off()
###############################Figure 3#############
CD4T.diff <- FindMarkers(d, ident.1 = "pws_CD4_T", ident.2 = "control_CD4_T", verbose = FALSE,logfc.threshold = 0.1)
write.csv(file="Diff_exp_genes_CD4T.csv",CD4T.diff)
CD8T.diff <- FindMarkers(d, ident.1 = "pws_CD8_T", ident.2 = "control_CD8_T", verbose = FALSE,logfc.threshold = 0.1)
write.csv(file="Diff_exp_genes_CD8T.csv",CD8T.diff)
mono.diff <- FindMarkers(d, ident.1 = "pws_monocyte", ident.2 = "control_monocyte", verbose = FALSE,logfc.threshold = 0.1)
write.csv(file="Diff_exp_genes_mono.csv",mono.diff)
NK.diff <- FindMarkers(d, ident.1 = "pws_NK_T", ident.2 = "control_NK_T", verbose = FALSE,logfc.threshold = 0.1)
write.csv(file="Diff_exp_genes_NK.csv",NK.diff)
B.diff <- FindMarkers(d, ident.1 = "pws_B_cell", ident.2 = "control_B_cell", verbose = FALSE,logfc.threshold = 0.1)
write.csv(file="Diff_exp_genes_B.csv",B.diff)
gama.diff <- FindMarkers(d, ident.1 = "pws_gama_T", ident.2 = "control_gama_T", verbose = FALSE,logfc.threshold = 0.1)
write.csv(file="Diff_exp_genes_gama.csv",gama.diff)
mono<-mono.diff[which(mono.diff$avg_log2FC>0.25),]
CD8<-CD8T.diff[which(CD8T.diff$avg_log2FC>0.25),]
B<-B.diff[which(B.diff$avg_log2FC>0.25),]
NK<-NK.diff[which(NK.diff$avg_log2FC>0.25),]
CD4<-CD4T.diff[which(CD4T.diff$avg_log2FC>0.25),]
gama<-gama.diff[which(gama.diff$avg_log2FC>0.25),]
library(org.Hs.eg.db)
library(clusterProfiler)
m<-rownames(mono)
erich.go.ALL = enrichGO(gene =m,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)

result=erich.go.ALL@result
write.csv(result,"mono.csv")
m<-rownames(CD8)
erich.go.ALL = enrichGO(gene =m,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)

result=erich.go.ALL@result
write.csv(result,"CD8.csv")
m<-rownames(CD4)
erich.go.ALL = enrichGO(gene =m,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)

result=erich.go.ALL@result
write.csv(result,"CD4.csv")
m<-rownames(B)
erich.go.ALL = enrichGO(gene =m,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)

result=erich.go.ALL@result
write.csv(result,"B.csv")
m<-rownames(NK)
erich.go.ALL = enrichGO(gene =m,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)
result=erich.go.ALL@result
write.csv(result,"NK.csv")
m<-rownames(gama)
erich.go.ALL = enrichGO(gene =m,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)

result=erich.go.ALL@result
write.csv(result,"gd.csv")
#######################################
w<-c("response to interleukin-1","response to interferon-gamma","interleukin-2 production",
     "I-kappaB kinase/NF-kappaB signaling","positive regulation of inflammatory response",
     "positive regulation of cell activation","stress-activated protein kinase signaling cascade",
     "stress-activated MAPK cascade")
ma2<-ma1[which(ma1$Description %in% w),]
ma2$celltype<-factor(ma2$celltype,levels = c("mono","NK","gama","CD4","CD8","B"))
ma2$Description<-factor(ma2$Description,levels = rev(c("response to interleukin-1","response to interferon-gamma","interleukin-2 production",
                                                       "I-kappaB kinase/NF-kappaB signaling","positive regulation of inflammatory response",
                                                       "positive regulation of cell activation","stress-activated protein kinase signaling cascade",
                                                       "stress-activated MAPK cascade")))

p = ggplot(ma2,aes(celltype,Description))
p1= p + geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_color_gradient(low = 'white',high = "#CD5B45",limits=c(-10,15))+
  theme_bw()
pdf("barplot.pdf",width =7.5,height =5)
p1+scale_size_continuous(range=c(2,8))
dev.off()
###########################################3
monocyte<-rownames(mono)
CD8_T <-rownames(CD8)
CD4_T <-rownames(CD4)
NK_cell <-rownames(NK)
B_cell <-rownames(B)
gama_cell<-rownames(gama)
library(UpSetR)
class(monocyte)
listInput <- list(
  monocyte = c(monocyte), 
  CD8_T = c(CD8_T), 
  CD4_T = c(CD4_T),
  NK_cell=c(NK_cell),
  B_cell=c(B_cell),
  gama_cell=c(gama_cell))
getwd()
saveRDS(listInput,"listInput.rds")
pdf("upset.pdf",width =10,height = 6)
upset(fromList(listInput),nsets = 6,order.by = "freq",main.bar.color = "black",matrix.color="red",set_size.numbers_size=8,set_size.show=10)
dev.off()
############################################
gene2<-c("OSM","IL1B","TNF","VEGFA","IFNG","HES1","IFITM3","EGR1","AIF1","PTGS2",
         "CCL3L1","CXCR4","CCL3","CCL4L2","CCL4","JUN","FOS","DUSP1","NFKBIA","CD69",
         "JUNB","JUND","IER2","IER5")
gene<-unique(gene)
av=AverageExpression(d)
dim(av$RNA)
dd1<-t(av$RNA[gene2,1:22])
class(dd1)
dd1<-as.data.frame(dd1)
dd1$celltype<-sapply(row.names(dd1), function(x) unlist(strsplit(x, "_"))[1]) 
dd1$ID<-sapply(row.names(dd1), function(x) unlist(strsplit(x, "_"))[2]) 
dd2<-t(dd1)
dim(dd2)
dd2<-dd2[-c(341,342),]
unique(colnames(dd2))
dd3<-dd2[,c("control_monocyte","pws_monocyte","control_NK_T","pws_NK_T","control_gama_T","pws_gama_T",
            "control_B_cell","pws_B_cell","control_CD4_T","pws_CD4_T","control_CD8_T","pws_CD8_T")]

data_scale <- as.data.frame(t(apply(dd3,1,scale)))
colnames(data_scale) <- colnames(dd3)

ann_colors = list( sample = c(control="#92C5DE",PWS="#D6604D"), 
                   celltype = c(monocyte ="#E4AF69FF", NK ="#990080FF",
                                gama ="#D60047FF", Bcell="#749B58FF",
                                CD4_T="#5DB1DDFF",CD8_T="#BA6338FF"))
annotation_col <- data.frame(
  celltype = rep(c("monocyte","NK","gama","Bcell","CD4_T","CD8_T"), each = 2),
  sample = c("control","PWS","control","PWS","control","PWS","control","PWS","control","PWS","control","PWS")
)

pdf("gene.pdf",width = 6,height = 6)
pheatmap(data_scale,scale = "row",
         col = colorRampPalette(c('#3B3B3B','white','red4'))(100),
         cluster_row = FALSE,cluster_cols = FALSE,
         border=FALSE,
         gaps_row = c(5,10,15),
         gaps_col = c(2,4,6,8,10),
         annotation_col = annotation_col,
         cellwidth = 15, cellheight = 12,
         fontsize = 8,
         annotation_colors = ann_colors)
dev.off()
##################################################W6上传##
library(corrplot)
cor_data <- cor(w6,method="pearson")
bk <- c(seq(-1.5,-0.1,by=0.01),seq(0,1,by=0.01))
ann_colors = list( group = c("CHF"="#92C5DE",pws="#D6604D"), 
                   celltype = c(monocyte ="#E4AF69FF", NK ="#990080FF",
                                gama ="#D60047FF", Bcell="#749B58FF",
                                CD4_T="#5DB1DDFF",CD8_T="#BA6338FF"))
pheatmap(cor_data,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks=bk,
         annotation_col = m,annotation_colors = ann_colors,
         border=FALSE,treeheight_row=20,treeheight_col=10)
getwd()
write.csv(w6,"w6.csv")
#######################################WGCNA#########
pbmc<-remono_sub
# 加载rds文件
# 细胞数较多，且单细胞表达矩阵稀疏，先讲细胞筛选一下
# 采用pseudocell概念
datadf <- as.matrix(pbmc@assays$RNA@data)
idd1 <- pbmc@meta.data
unique(pbmc@meta.data$stim)
Inter.id1<-cbind(rownames(idd1),idd1$stim)
rownames(Inter.id1)<-rownames(idd1)
colnames(Inter.id1)<-c("CellID","stim")
Inter.id1<-as.data.frame(Inter.id1)
head(Inter.id1)
head(datadf)[1:4,1:4]
Inter1<-datadf[,Inter.id1$CellID]
Inter2<-as.matrix(Inter1)
Inter2[1:4,1:4]
Inter.id1$stim<-factor(Inter.id1$stim)
table(pbmc@meta.data$group3)
##########################################
pseudocell.size = 20
new_ids_list1 = list()
length(levels(Inter.id1$stim))
for (i in 1:length(levels(Inter.id1$stim))) {
  cluster_id = levels(Inter.id1$stim)[i]
  cluster_cells <- rownames(Inter.id1[Inter.id1$stim == cluster_id,])
  cluster_size <- length(cluster_cells)     
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)    
  new_ids_list1[[i]] <- pseudo_ids      
}

new_ids <- unlist(new_ids_list1)
new_ids <- as.data.frame(new_ids)
head(new_ids)
new_ids_length <- table(new_ids)
new_ids_length
new_colnames <- rownames(new_ids)
#######################################################3
gc()
colnames(datadf)  
all.data<-datadf[,as.character(new_colnames)] ###add
all.data <- t(all.data)###add
new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                    list(name=new_ids[,1]),FUN=mean)
head(new.data)[1:4,1:4]
rownames(new.data)<-new.data$name
new.data<-new.data[,-1]
new_ids_length<-as.matrix(new_ids_length)##
short<-which(new_ids_length<10)##
new_good_ids<-as.matrix(new_ids_length[-short,])##
result<-t(new.data)[,rownames(new_good_ids)]
dim(result)
######################################33
type = "signed"  #
corType = "pearson" # 
corFnc = ifelse(corType=="pearson", cor, bicor)
corFnc
maxPOutliers = ifelse(corType=="pearson",1,0.05) 
robustY = ifelse(corType=="pearson",T,F)
dataExpr  <- as.matrix(Cluster1)
dataExpr=t(dataExpr)
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
gsg$goodSamples
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
softPower  = 0.5
cor <- WGCNA::cor
############################################
net = blockwiseModules(dataExpr, power =6, maxBlockSize = nGenes,#nGenes
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("dataExpr", ".tom"),
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)
pdf("modules.pdf",12,6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
mypbmc<- CreateSeuratObject(Cluster1)
design=model.matrix(~0+ mypbmc@meta.data$group2)
colnames(design)=levels(mypbmc@meta.data$group2)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("cor_PWS.pdf",6,6)
par(mar=c(5,10,4,5), oma=c(0, 0, 0, 0))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
############################################################
moduleTraitCor_noFP <- cor(MEs[,-12], mypbmc@meta.data[,7:11], use = "p")
moduleTraitPvalue_noFP = corPvalueStudent(moduleTraitCor_noFP, nSamples);
textMatrix_noFP <- paste(signif(moduleTraitCor_noFP, 2), "\n(", signif(moduleTraitPvalue_noFP, 1), ")", sep = "")
METree = hclust(dist(t(MEs[,1:11])), method = "average")
traitColors=numbers2colors(moduleTraitCor_noFP, signed = TRUE,centered=TRUE,lim=c(-1,1))
pdf("module and trait heatmap.pdf",h=4,w=5)
plotDendroAndColors(METree, traitColors,
                    groupLabels = colnames(mypbmc@meta.data[,7:11]),
                    main = "mudule and trait heatmap")
dev.off()
pdf("all_biaoxing.pdf",width = 6,height = 7)
par(mar = c(5,10, 3, 3)); 
labeledHeatmap(Matrix = moduleTraitCor_noFP, 
               xLabels = names(mypbmc@meta.data[,7:11]), 
               yLabels = names(MEs[,-12]), 
               ySymbols = names(MEs[,-12]), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix_noFP,
               setStdMargins = FALSE, 
               cex.text = 0.65, 
               zlim = c(-1,1), 
               main = paste("Module-trait relationships")) 
dev.off()
##############################Figure 4##################
DefaultAssay(mono_sub) <- "RNA"
Subset_Cells.list <- SplitObject(mono_sub, split.by = "stim")
for (i in 1:length(Subset_Cells.list)) {
  Subset_Cells.list[[i]] <- NormalizeData(Subset_Cells.list[[i]], verbose = FALSE)
  Subset_Cells.list[[i]] <- FindVariableFeatures(Subset_Cells.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}
Integrated <- FindIntegrationAnchors(object.list = Subset_Cells.list, dims = 1:15)
Integrated <- IntegrateData(anchorset = Integrated, dims = 1:15)
DefaultAssay(object = mono_sub1) <- "integrated"
mono_sub=ScaleData(Integrated,verbose = FALSE)
mono_sub <- RunPCA(object = mono_sub, npcs = 30, verbose = FALSE)
mono_sub <- RunUMAP(object = mono_sub, reduction = "pca", dims = 1:15)
mono_sub <- FindNeighbors(object = mono_sub, reduction = "pca", dims = 1:15)
mono_sub <- FindClusters(mono_sub, resolution = 0.4)

mono.markers <- FindAllMarkers(mono_sub, only.pos = TRUE, 
                               min.pct = 0.3, logfc.threshold = 0.25)
write.table(mono.markers,file="mono.markers-0.3-10.txt",sep="\t",quote = ,row.names = F)

mono_sub@meta.data$celltype2<-as.character(mono_sub@meta.data$seurat_clusters)
c=which(mono_sub@meta.data$celltype2=="0")
mono_sub@meta.data[c,"celltype2"]<-"classical monocyte"
c=which(mono_sub@meta.data$celltype2=="1")
mono_sub@meta.data[c,"celltype2"]<-"S100 monocyte"
c=which(mono_sub@meta.data$celltype2=="2")
mono_sub@meta.data[c,"celltype2"]<-"S100 monocyte"
c=which(mono_sub@meta.data$celltype2=="3")
mono_sub@meta.data[c,"celltype2"]<-"IFN monocyte"
c=which(mono_sub@meta.data$celltype2=="4")
mono_sub@meta.data[c,"celltype2"]<-"CD16 monocyte"
c=which(mono_sub@meta.data$celltype2=="5")
mono_sub@meta.data[c,"celltype2"]<-"HLA_1 monocyte"
c=which(mono_sub@meta.data$celltype2=="6")
mono_sub@meta.data[c,"celltype2"]<-"HLA_2 monocyte"
c=which(mono_sub@meta.data$celltype2=="7")
mono_sub@meta.data[c,"celltype2"]<-"T like monocyte"
c=which(mono_sub@meta.data$celltype2=="8")
mono_sub@meta.data[c,"celltype2"]<-"NK like monocyte"
c=which(mono_sub@meta.data$celltype2=="9")
mono_sub@meta.data[c,"celltype2"]<-"macrophages"
c=which(mono_sub@meta.data$celltype2=="10")
mono_sub@meta.data[c,"celltype2"]<-"B like monocyte"
saveRDS(mono_sub,"remono_sub.rds")

color1<-c("#4DBBD5FF","#00A087FF", "#3C5488FF","#E64B35FF","#F39B7FFF","#BB002199",
          "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")

Idents(mono_sub)<-"celltype2"
pdf(paste0("./",result.name,"/momo_celltype2",max(dim.use),"PC.pdf"),width = 10,height = 4)
mono_sub@meta.data$celltype2<-factor(mono_sub@meta.data$celltype2, 
                                     levels=c("classical monocyte","S100 monocyte",
                                              "IFN monocyte","CD16 monocyte",
                                              "HLA_1 monocyte","HLA_2 monocyte",
                                              "T like monocyte","NK like monocyte",
                                              "macrophages","B like monocyte"))

DimPlot(object =mono_sub, pt.size=0.01,label = F, cols=color1)
dev.off()
#######################################
library(remotes)
remotes::install_github("lyc-1995/MySeuratWrappers")
library(MySeuratWrappers)
feature2<-c("LGALS2",
            "S100A12",
            "ISG15",
            "FCGR3A",
            "HLA-DRA",
            "FCER1A",
            "TCF7",
            "CCL5",
            "PF4",
            "CD79A")
pdf("Volin_celltype2.pdf",width = 8,height = 6)
VlnPlot(mono_sub,group.by = "celltype2",
        features = feature2,stacked=T,pt.size=0,cols =color1,direction = "horizontal",x.lab = '', y.lab = '')+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) 
dev.off()
###################################################
######################################################
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
#devtools::install_github("JinmiaoChenLab/cytofkit")
#https://github.com/XinleiChen001/cytofexplorer
devtools::install_github("XinleiChen001/cytofexplorer")
rm(list = ls())#清除内存变量
library(flowCore)
library(Rcpp)
library(cytofkit)
library(igraph)
library(ggplot2)
library(ggthemes)
library(Rtsne)
library(dplyr)
library(cytofexplorer)
projectdir<-"E:/SCT/cyto/CYTOF_paper2"
mergeMethod <- "ceil" 
fixedNum <- 2000    
wdir <-"21_PhenoGraph_tSNE"   #
raw_fcs_dir="01_rawfcs"
wdir<-paste0(projectdir,"/",wdir)
raw_fcs_dir<-paste0(projectdir,"/",raw_fcs_dir)
metadata_dir<-paste0(projectdir,"/02_metadata")
setwd("/SCT/cyto/CYTOF_paper2/21_PhenoGraph_tSNE/")
file_name <- list.files(raw_fcs_dir,pattern='.fcs$', full=TRUE)
combined_data_raw <- cytof_exprsMerge(fcsFiles = file_name,
                                      transformMethod = "none",
                                      mergeMethod =mergeMethod,
                                      fixedNum = fixedNum)
paraname<-colnames(combined_data_raw)
paraname<-sub(".*<","",paraname)
paraname<-sub(">.*","",paraname)
paraname<-sub("-","_",paraname)
colnames(combined_data_raw)<-paraname
File_ID<-sub("_[0-9]*$","",row.names(combined_data_raw))
combined_data_raw<-data.frame(combined_data_raw,File_ID)
head(combined_data_raw)

all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
transform_id=which(all_markers$transform==1)
combined_data_transformed<-combined_data_raw

cytofAsinh <- function(value, cofactor = 5) {
  value <- value-1
  loID <- which(value < 0)
  if(length(loID) > 0)
    value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  value <- value / cofactor
  value <- asinh(value) 
  return(value)
}
combined_data_transformed[, transform_id] <- apply(combined_data_transformed[, transform_id,drop = FALSE],
                                                   2,cytofAsinh)
#checkpoint
head(combined_data_transformed)

###################################33
k=25 
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
PhenoGraph_id=which(all_markers$PhenoGraph==1)
PhenoGraph_input_data=combined_data_transformed[,PhenoGraph_id]
PhenoGraph_result <-as.numeric(membership(Rphenograph(data = PhenoGraph_input_data,k=k)))
#Checkpoint
pdf("PhenoGraph_hist.pdf")
hist(PhenoGraph_result,unique(PhenoGraph_result))
dev.off()
max_iter=1500   
perplexity=30  
theta=0.5      
dims = 2       

if (exists('seed')) set.seed(seed)
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
tSNE_para_id=which(all_markers$tSNE==1)
tSNE_input_data=combined_data_transformed[,tSNE_para_id]
tsne_result <- Rtsne(tSNE_input_data,
                     initial_dims = ncol(tSNE_input_data),
                     pca = FALSE,
                     dims = dims,
                     check_duplicates = FALSE,
                     perplexity=perplexity,
                     max_iter=max_iter,
                     theta=theta)$Y
row.names(tsne_result)<-row.names(combined_data_raw)
colnames(tsne_result)<-c("tsne_1","tsne_2")
pdf("tsne_result.pdf")
plot(tsne_result)
dev.off()

write.csv(unique(File_ID),paste0(metadata_dir,"/all_samples.csv"),row.names = FALSE)
groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors=FALSE)
groups$File_ID<-unique(File_ID)
####groups$Group <- paste(groups$Diet, groups$Tissue_Type, sep="_")
#checkpoint
groups
combined_data_analysed <- cbind(combined_data_transformed,
                                tsne_result,
                                PhenoGraph = PhenoGraph_result)
combined_data_plot <- combined_data_analysed
#checkpoint
head(combined_data_plot)
draw_tsne_figs(combined_data_plot=combined_data_plot,
               groups=groups,
               cluster_color=dif_seq_rainbow,
               cluster_name="PhenoGraph",
               major_cond="Tissue_Type",
               reduction_dm1="tsne_1",
               reduction_dm2="tsne_2")
unique(colnames(groups))
aa<-combined_data_plot
#########################
aa$celltype<-as.character(aa$PhenoGraph)
aa$celltype[which(aa$celltype=="6")]<-"CD4 T cell"
aa$celltype[which(aa$celltype=="8")]<-"CD4 T cell"
aa$celltype[which(aa$celltype=="16")]<-"CD4 T cell"
aa$celltype[which(aa$celltype=="17")]<-"CD4 T cell"
aa$celltype[which(aa$celltype=="18")]<-"CD4 T cell"
aa$celltype[which(aa$celltype=="12")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="13")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="14")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="3")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="23")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="24")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="4")]<-"CD8 T cell"
aa$celltype[which(aa$celltype=="1")]<-"NK cell"
aa$celltype[which(aa$celltype=="5")]<-"NK cell"
aa$celltype[which(aa$celltype=="15")]<-"NK cell"
aa$celltype[which(aa$celltype=="7")]<-"NK cell"
aa$celltype[which(aa$celltype=="2")]<-"B cell"
aa$celltype[which(aa$celltype=="22")]<-"B cell"
aa$celltype[which(aa$celltype=="9")]<-"B cell"
aa$celltype[which(aa$celltype=="11")]<-"CD4-CD8-  T cell"
aa$celltype[which(aa$celltype=="20")]<-"CD4-CD8-  T cell"
aa$celltype[which(aa$celltype=="10")]<-"CD14 monocyte"
aa$celltype[which(aa$celltype=="19")]<-"CD16 monocyte"
aa$celltype[which(aa$celltype=="21")]<-"IM"
unique(aa$celltype)
##################
unique(aa$celltype)
aa$celltype<-factor(aa$celltype,levels = c("NK cell",#"#802268FF"
                                           "B cell",#"#749B58FF"
                                           "CD8 T cell",#"#CE3D32FF"
                                           "CD4 T cell",#"#5050FFFF",
                                           "CD4-CD8-  T cell",#"#D595A7FF
                                           "CD14 monocyte",#"#FFC20AFF"
                                           "CD16 monocyte",#"#CC9900FF",
                                           "IM"))

########################################
col3<-c("#802268FF","#749B58FF","#CE3D32FF","#5050FFFF", "#D595A7FF","#FFC20AFF","#CC9900FF","#CD853F")
pdf("celltype2.pdf", width=12,height=10)
ggplot(aa, aes(x=tsne_1,y=tsne_2,color=celltype))+
  geom_point()+theme_bw()+
  scale_color_manual(values=col3)
dev.off()
##############################################
setwd("E:/SCT/cyto/CYTOF_paper2/")
aa$PhenoGraph <- factor(aa$PhenoGraph)

pdf("tsne_plot2.pdf", width=12,height=10)
ggplot(aa, aes(x=tsne_1,y=tsne_2,color=PhenoGraph))+
  geom_point()+theme_bw()+scale_color_manual(values=c(mycolor(10)))
dev.off()

########################################################
aa<-combined_data_plot
aa$X158Gd_CD16=as.numeric(aa$X158Gd_CD16)

pdf("CD16.pdf",width=8,height=8)
ggplot(aa,mapping = aes(x=tsne_1,
                        y=tsne_2,
                        col=X158Gd_CD16))+geom_point()+
  scale_color_gradient(low="#E3E3E3", high="#EE3B3B")+ 
  theme_bw()


dev.off()

###############################################
pdf("CD3d.pdf",width=8,height=8)
aa$X175Lu_CD3_EQ6=as.numeric(aa$X175Lu_CD3_EQ6)
ggplot(aa,mapping = aes(x=tsne_1,
                        y=tsne_2,
                        col=X175Lu_CD3_EQ6))+geom_point()+
  scale_color_gradient(low="#E3E3E3", high="#EE3B3B")+ 
  theme_bw()
dev.off()
#########################################################33
pdf("CD4.pdf",width=8,height=8)
aa$X155Gd_CD4=as.numeric(aa$X155Gd_CD4)
ggplot(aa,mapping = aes(x=tsne_1,
                        y=tsne_2,
                        col=X155Gd_CD4))+geom_point()+
  scale_color_gradient(low="#E3E3E3", high="#EE3B3B")+ 
  theme_bw()
dev.off()
########################################
pdf("CD8.pdf",width=8,height=8)
aa$X168Er_CD8=as.numeric(aa$X168Er_CD8)
ggplot(aa,mapping = aes(x=tsne_1,
                        y=tsne_2,
                        col=X168Er_CD8))+geom_point()+
  scale_color_gradient(low="#E3E3E3", high="#EE3B3B")+ 
  theme_bw()
dev.off()

###############################################33

pdf("CD56_1.pdf",width=8,height=8)
aa$X176Yb_CD56_EQ7=as.numeric(aa$X176Yb_CD56_EQ7)
ggplot(aa,mapping = aes(x=tsne_1,
                        y=tsne_2,
                        col=X176Yb_CD56_EQ7))+geom_point()+
  scale_color_gradient(low="#E3E3E3", high="#EE3B3B")+ 
  theme_bw()
dev.off()

###############################
aa$X142Nd_CD14_EQ2=as.numeric(aa$X142Nd_CD14_EQ2)
pdf("CD14.pdf",width=8,height=8)
ggplot(aa,mapping = aes(x=tsne_1,
                        y=tsne_2,
                        col=X142Nd_CD14_EQ2))+geom_point()+
  scale_color_gradient(low="#E3E3E3", high="#EE3B3B")+ 
  theme_bw()

dev.off()
#######################################Figure 5#########
library(CellChat)
library(patchwork)
cellchat.control <- createCellChat(object =control@assays$RNA@data, meta =control@meta.data,  group.by ="celltype4")
cellchat.PWS <- createCellChat(object =PWS@assays$RNA@data, meta =PWS@meta.data,  group.by ="celltype4")
#############################
cellchat=cellchat.control 
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.control = cellchat
#################################
cellchat=cellchat.PWS
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.PWS = cellchat
##############################################
cc.list=list(con=cc.control,pws=cc.PWS)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "count")
netVisual_diffInteraction(cellchat,weight.scale = T)
pdf(file = "differ-circle.pdf",width =6,height = 6)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
dev.off()
netVisual_heatmap(cellchat)
pdf(file ="diffeheatmap1.pdf",width =6,height = 6)
netVisual_heatmap(cellchat,measure = "weight")
dev.off()
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "control" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "pws" )

pdf(file ="juleiduizhao.pdf",width =4,height = 4)
netAnalysis_signalingRole_scatter(cc.list[[1]], title = names(cc.list)[1], weight.MinMax = weight.MinMax)+
  scale_y_continuous(limits = c(0,0.09))+
  scale_x_continuous(limits = c(0,0.1))
dev.off()

############################################
num.link <- sapply(cc.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cc.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cc.list[[i]], title = names(cc.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

object.list<-cc.list
###############################################3
dd<-object.list$pws@netP$pathways
unique(dd)
pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("VISFATIN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pdf(file ="VISFATIN.pdf",width =4,height = 4)
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[1]))
dev.off()

pdf(file ="VISFATIN2.pdf",width =4,height = 4)
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[2]))
dev.off()
pathways.show <- c("TNF") 
pdf(file ="TNF1.pdf",width =4,height = 4)
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[1]))
dev.off()
pdf(file ="TNF2.pdf",width =4,height = 4)
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[2]))
dev.off()
##################MR#################
library(mediation)
library(data.table)
Meddata<-fread("./F6MR2.txt")
fitM <- lm(M ~ X,     data=Meddata) 
fitY <- lm(Y ~ X + M, data=Meddata) 
fitMed <- mediation::mediate(fitM, fitY, treat="X", mediator="M",sims=2300)
summary(fitMed)
plot(fitMed)
########################












