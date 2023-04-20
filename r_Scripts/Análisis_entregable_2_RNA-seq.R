#Análisis entregable 2 RNA-seq##################################################
#Librerias
library("AnnotationDbi")
library("org.Mm.eg.db")
library("clusterProfiler")
library("DESeq2")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library(biomaRt)
library("genefilter")
library("ggrepel")
library(pheatmap)
library(ashr)
library(ggVennDiagram)
library("dplyr")
library("ComplexHeatmap")
library("readr")
library("tibble")
library("apeglm")
library("openxlsx")
library(plotly)
library(GeneOverlap)
#Cargamos los datos###############################
##Astrocitos###################################### 
###Matriz de cuentas##############################
cts_astrocytes <- read_table("AnalisisR/OmicosEntregable2/astrocitos/FeatureCounts2/featurecountstotal_STAR.mat")
cts_astrocytes <- cts_astrocytes %>% column_to_rownames(var = "Geneid")
cts_astrocytes<-as.matrix(cts_astrocytes)
###Metadatos######################################
coldata_astrocytes <- read.csv("AnalisisR/OmicosEntregable2/astrocitos/coldata_astrocytes.csv",sep="," ,row.names=1)
coldata_astrocytes <- coldata_astrocytes[,c("Nucleo","Region","Celula")]
coldata_astrocytes$Region <- factor(coldata_astrocytes$Region)
coldata_astrocytes$Celula <- factor(coldata_astrocytes$Celula)
coldata_astrocytes$Nucleo <- factor(coldata_astrocytes$Nucleo)
##Neuronas########################################
cts_Tneurons <- read_table("AnalisisR/OmicosEntregable2/neurons/FeatureCounts_thalamus/featurecountstotal_STAR.mat")
cts_Tneurons <- cts_Tneurons %>% column_to_rownames(var = "Geneid")
cts_Tneurons<-as.matrix(cts_Tneurons)

#Analisis Astrocitos Region_Talamo_vs_Corteza#############################################

dds <- DESeqDataSetFromMatrix(countData = cts_astrocytes,
                              colData = coldata_astrocytes,
                              design = ~ Region)
##PRE-FILTERING################################################################
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rld <- rlog(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
##PCA###########################################################################
data<-plotPCA(rld, intgroup = "Region",returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Region, data=data, label=colnames(rld)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/PCA_CXvsTH.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rld)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Region),                 
                                                                                   alpha = 0.2, 
                                                                                   
                                                                                   show.legend = FALSE, 
                                                                                   
                                                                                   level = 0.75
) 
dev.off()
##Resultado#####################################################################
dds <- DESeq(dds)
res<- results(dds)
#res.ape <- lfcShrink(dds=dds
#                     ,res=res
#                     , coef="Region_Talamo_vs_Corteza"
#                     , type="apeglm")


#Tabla resultado
dim(res)
head(res)
summary(res)
sigs <- na.omit(res)
sigs<- sigs[sigs$padj < 0.1,]
df <- as.data.frame(sigs)
dim(df)

##Anotacion####################################################################
listEnsemblArchives()
ensembl_genes <- useMart('ENSEMBL_MART_ENSEMBL',
                         host =  'https://jan2019.archive.ensembl.org')

mouse <- useDataset("mmusculus_gene_ensembl", ensembl_genes)

# ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = "104")
rownames_ordenados <- sort(rownames(df))
datos_ordenados <- df[rownames_ordenados, ]
values <- rownames(datos_ordenados)
# lista_atributos<-listAttributes(ensembl)
tabla_anotacion<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "gene_biotype","description") ,
                       filters = "ensembl_gene_id",
                       values = rownames(df),
                       mart = mouse)

dim(tabla_anotacion)
dim(df)
# # Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combined <- merge(df, tabla_anotacion %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_final<- unique(df_combined)
rownames(df_final)<-df_final$Row.names
df_final$Row.names<-NULL
head(df_final)

##Heatmap##########################################################
#library("genefilter")
#filtramos los genes para quedarnos con los significativos
df.top <- df_final[(abs(df_final$log2FoldChange) > 0),]
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]

df.top$diffexpressed <- "NO"
df.top$diffexpressed[df.top$log2FoldChange > 0] <- "UP"
df.top$diffexpressed[df.top$log2FoldChange < 0 ] <- "DOWN"
table(df.top$diffexpressed)
#creamos la matriz
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), rownames(coldata_astrocytes)] #sig genes x samples
rownames(mat)<- df.top$mgi_symbol
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- paste0(rld$Region,"-",rld$Nucleo)

#heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Region ]
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/Heatmap_CXvsTH.png",
    width=1000, height = 1000)
heatmap.2(mat.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

##Volcano###################################################################
#seleccionamos los genes expresados diferencialmente
df_final$diffexpressed <- "NO"
df_final$diffexpressed[df_final$log2FoldChange > 0 & df_final$padj < 0.1] <- "UP"
df_final$diffexpressed[df_final$log2FoldChange < 0 & df_final$padj < 0.1] <- "DOWN"
#etiquetamos los genes diferenciados
df_final$delabel <- NA
df_final$delabel[df_final$diffexpressed != "NO"] <- df_final$mgi_symbol[df_final$diffexpressed != "NO"]
table(df_final$diffexpressed)
#representamos con ggplot
#library(ggrepel)
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/Volcano_CXvsTH.png",
    width=1000, height = 1000)
ggplot(data=df_final, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("grey", "blue")) +
  geom_vline(xintercept=c(0, 0), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")
dev.off()

##Enriquecimiento###############################################################
###Talamo#######################################################################
genesTalamo <- df.top[(df.top$log2FoldChange) > 0.322, ]

GO_results <- enrichGO(gene = rownames(genesTalamo),
                       OrgDb = "org.Mm.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP", 
                       pAdjustMethod = "BH",
                       readable = TRUE)
#write.csv(GO_resultsgenes_CautBP@result, "AnalisisR/OmicosEntregable2/astrocitos/GO_resultsgenes_CautBP.csv")
ids_go <- c("GO:0007389","GO:0003002","GO:0021953","GO:0045165","GO:0030902","GO:0009952","GO:0035282","GO:0022037","GO:0021536","GO:0030901")
filtered_result <- GO_results[GO_results@result$ID %in% ids_go, ]
GO_results2 <- GO_results
GO_results2@result <- GO_results@result[ids_go,]

####Barplot############
barplotTh<-ggplot(filtered_result, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradient(low = "blue", high = "lightblue") +
  labs(title = "AS-Th",
       subtitle = "GO BP enrichment analysis",
       x = "Counts")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,hjust = 1, vjust = 0.5,size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/barplot_TH.png",
    width=1000, height = 1000)
barplotTh
dev.off()
####Cnetplot################################
genesTalamoLFC<-genesTalamo$log2FoldChange
names(genesTalamoLFC) <- genesTalamo$mgi_symbol
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/cnetplot_TH.png",
    width=1000, height = 1000)
cnetplot(GO_results2,color.params = list(foldChange=genesTalamoLFC),showCategory = 5)+scale_color_gradient2(name='log2FC', low='darkgreen', high='firebrick')
dev.off()
###Corteza######################################################################
genesCorteza<- df.top[(df.top$log2FoldChange) < -0.322, ]
GO_resultsCort <- enrichGO(gene = rownames(genesCorteza),
                       OrgDb = "org.Mm.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       readable = TRUE)
#write.csv(GO_resultsgenes_CautBP@result, "AnalisisR/OmicosEntregable2/astrocitos/GO_resultsgenes_CautBP.csv")
ids_goCX <- c("GO:0045165","GO:0021537","GO:0021543","GO:0061351","GO:0021872","GO:0048663","GO:0021987","GO:0030900","GO:0021895")
filtered_resultCX <- GO_resultsCort[GO_resultsCort@result$ID %in% ids_goCX, ]
GO_resultsCort2 <- GO_resultsCort
GO_resultsCort2@result <- GO_resultsCort@result[ids_goCX,]

####Barplot############
barplotCX<-ggplot(filtered_resultCX, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradient(low = "black", high = "lightgrey") +
  labs(title = "AS-CX",
       subtitle = "GO BP enrichment analysis",
       x = "Counts")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90,hjust = 1, vjust = 0.5,size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/barplot_CX.png",
    width=1000, height = 1000)
barplotCX
dev.off()
####Cnetplot################################
genesCX_LFC<-genesCorteza$log2FoldChange
names(genesCX_LFC) <- genesCorteza$mgi_symbol
png(filename = "AnalisisR/OmicosEntregable2/astrocitos/cnetplot_CX.png",
    width=1000, height = 1000)
cnetplot(GO_resultsCort2,color.params = list(foldChange=genesCX_LFC),showCategory = 5)+scale_color_gradient2(name='log2FC', low='darkgreen', high='firebrick')
dev.off()