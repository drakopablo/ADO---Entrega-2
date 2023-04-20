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
cts_astrocytes <- read_table("OmicosEntregable2/astrocitos/FeatureCounts2/featurecountstotal_STAR.mat")
cts_astrocytes <- cts_astrocytes %>% column_to_rownames(var = "Geneid")
cts_astrocytes<-as.matrix(cts_astrocytes)
cts_As_CX<-cts_astrocytes[, 1:4]
###Metadatos######################################
coldata_astrocytes <- read.csv("OmicosEntregable2/astrocitos/coldata_astrocytes.csv",sep="," ,row.names=1)
coldata_astrocytes <- coldata_astrocytes[,c("Nucleo","Region","Celula")]
coldata_astrocytes$Region <- factor(coldata_astrocytes$Region)
coldata_astrocytes$Celula <- factor(coldata_astrocytes$Celula)
coldata_astrocytes$Nucleo <- factor(coldata_astrocytes$Nucleo)

##Neuronas########################################
###Matriz de cuentas##############################
cts_CXneurons <- read_table("OmicosEntregable2/neurons/FeatureCountsCortex/featurecountstotal_STAR.mat")
cts_CXneurons <- cts_CXneurons %>% column_to_rownames(var = "Geneid")
cts_THneurons <- read_table("OmicosEntregable2/neurons/FeatureCounts_thalamus/featurecountstotal_STAR_th_ns_multi.mat")
cts_THneurons <- cts_THneurons %>% column_to_rownames(var = "Geneid")

cts_neurons<-cbind(cts_CXneurons,cts_THneurons)

###Metadatos######################################
coldata_neurons <- read.csv("OmicosEntregable2/neurons/coldata_neurons.csv",sep="," ,row.names=1)
coldata_neurons <- coldata_neurons[,c("Nucleo","Region","Celula")]
coldata_neurons$Region <- factor(coldata_neurons$Region)
coldata_neurons$Celula <- factor(coldata_neurons$Celula)
coldata_neurons$Nucleo <- factor(coldata_neurons$Nucleo)
#Analisis Astrocitos Region_Talamo_vs_Corteza#############################################
##Creamos el objeto deseq#################################################################
dds <- DESeqDataSetFromMatrix(countData = cts_astrocytes,
                              colData = coldata_astrocytes,
                              design = ~ Region)
##PRE-FILTERING################################################################
dds
keep <- rowSums(counts(dds)) >= 10 #consiideramos como validos solo los genes con más de 10 counts
dds <- dds[keep,]
rld <- rlog(dds)#Normalización de los datos
##PCA###########################################################################
data<-plotPCA(rld, intgroup = "Region",returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Region, data=data, label=colnames(rld)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "OmicosEntregable2/astrocitos/PCA_CXvsTH.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rld)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Region),                 
                                                                                   alpha = 0.2, 
                                                                                   
                                                                                   show.legend = FALSE, 
                                                                                   
                                                                                   level = 0.75
) 
dev.off()
##Expresión diferencial Resultado#####################################################################
dds <- DESeq(dds)
res<- results(dds)
summary(res)#1823 genes upregulados(Talamo) y 1364 genes downregulados (Cortex)

#Filtramos por padj
sigs <- na.omit(res)
sigs<- sigs[sigs$padj < 0.1,]
df <- as.data.frame(sigs)
dim(df)

##Anotacion con ensembl####################################################################
listEnsemblArchives()
ensembl_genes <- useMart('ENSEMBL_MART_ENSEMBL',
                         host =  'https://jan2019.archive.ensembl.org')

mouse <- useDataset("mmusculus_gene_ensembl", ensembl_genes)

# lista_atributos<-listAttributes(ensembl)
tabla_anotacion<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "gene_biotype","description") ,
                       filters = "ensembl_gene_id",
                       values = rownames(df),
                       mart = mouse)

dim(tabla_anotacion)

# # Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combined <- merge(df, tabla_anotacion %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_final<- unique(df_combined)
rownames(df_final)<-df_final$Row.names
df_final$Row.names<-NULL
head(df_final)
##Visualizacion de los datos###################################################

###Heatmap######################################################################
#filtramos los genes para quedarnos con los significativos
df.top <- df_final[(abs(df_final$log2FoldChange) > 0),]
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]

df.top$diffexpressed <- "NO"
df.top$diffexpressed[df.top$log2FoldChange > 0] <- "UP"
df.top$diffexpressed[df.top$log2FoldChange < 0 ] <- "DOWN"

#En resultado del analisis,el numerador es el Talamo,y por tanto, los genes upregulados son genes específicos de este. Mientras que los genes downregulados son específicos de la corteza
table(df.top$diffexpressed)#Numero de genes expresados diferencialmente

#creamos la matriz
rlog_out <- rlog(dds, blind=FALSE) #normalización de los datos
mat<-assay(rlog_out)[rownames(df.top), rownames(coldata_astrocytes)] #sig genes x samples
rownames(mat)<- df.top$mgi_symbol
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- paste0(rld$Region,"-",rld$Nucleo)

#heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Region ]
png(filename = "OmicosEntregable2/astrocitos/Heatmap_CXvsTH.png",
    width=1000, height = 1000)
heatmap.2(mat.scaled, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

###MAplot###################################################################
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
png(filename = "OmicosEntregable2/astrocitos/MaPlot_CXvsTH.png",
    width=1000, height = 1000)
ggplot(data=df_final, aes(x=-log10(baseMean), y=log2FoldChange, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("grey", "blue")) +
  geom_hline(yintercept=c(0, 0), col="red") 
  #geom_hline(xintercept=-log10(0.1), col="red")
dev.off()

###Volcano######################################################################
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
#En resultado del analisis,el numerador es el Talamo,y por tanto, los genes upregulados son genes específicos de este. Mientras que los genes downregulados son específicos de la corteza
###Talamo#######################################################################
genesTalamo <- df.top[(df.top$log2FoldChange) > 0.322, ] #seleccionamos los genes especificos del talamo

GO_results <- enrichGO(gene = rownames(genesTalamo),
                       OrgDb = "org.Mm.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP", 
                       pAdjustMethod = "BH",
                       readable = TRUE)
#Seleccionamos los procesos biologicos relacionados con el talamo
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
png(filename = "OmicosEntregable2/astrocitos/barplot_TH.png",
    width=1000, height = 1000)
barplotTh
dev.off()
####Cnetplot################################
genesTalamoLFC<-genesTalamo$log2FoldChange
names(genesTalamoLFC) <- genesTalamo$mgi_symbol
png(filename = "OmicosEntregable2/astrocitos/cnetplot_TH.png",
    width=1000, height = 1000)
cnetplot(GO_results2,color.params = list(foldChange=genesTalamoLFC),showCategory = 5)+scale_color_gradient2(name='log2FC', low='darkgreen', high='firebrick')
dev.off()
###Corteza######################################################################
genesCorteza<- df.top[(df.top$log2FoldChange) < -0.322, ]#Seleccionamos los genes específicos de la corteza
GO_resultsCort <- enrichGO(gene = rownames(genesCorteza),
                       OrgDb = "org.Mm.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       readable = TRUE)

#Seleccionamos los procesos biologicos relacionados con la corteza
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
png(filename = "OmicosEntregable2/astrocitos/barplot_CX.png",
    width=1000, height = 1000)
barplotCX
dev.off()
####Cnetplot################################
genesCX_LFC<-genesCorteza$log2FoldChange
names(genesCX_LFC) <- genesCorteza$mgi_symbol
png(filename = "OmicosEntregable2/astrocitos/cnetplot_CX.png",
    width=1000, height = 1000)
cnetplot(GO_resultsCort2,color.params = list(foldChange=genesCX_LFC),showCategory = 5)+scale_color_gradient2(name='log2FC', low='darkgreen', high='firebrick')
dev.off()

#Analisis Neuronas Region_Talamo_vs_Corteza#####################################
##Creamos el objeto deseq#################################################################

ddsN <- DESeqDataSetFromMatrix(countData = cts_neurons,
                              colData = coldata_neurons,
                              design = ~ Region)
##PRE-FILTERING################################################################
ddsN
keepN <- rowSums(counts(ddsN)) >= 10
ddsN <- ddsN[keepN,]
rldN <- rlog(ddsN)#Datos normalizados
##PCA###########################################################################
dataN<-plotPCA(rldN, intgroup = "Region",returnData=TRUE)
percentVar <- round(100 * attr(dataN, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Region, data=dataN, label=colnames(rldN)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "OmicosEntregable2/neurons/PCA_CXvsTH_neurons.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rldN)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Region),                 
                                                                                   alpha = 0.2, 
                                                                                   
                                                                                   show.legend = FALSE, 
                                                                                   
                                                                                   level = 0.75
) 
dev.off()

##Resultado#####################################################################
ddsN <- DESeq(ddsN)
resN<- results(ddsN)
resN.ape <- lfcShrink(dds=ddsN,
                    res=resN,
                    coef="Region_Talamo_vs_Corteza",
                    type="apeglm")


#Tabla resultado
dim(resN)
head(resN)
summary(resN)
sigsN <- na.omit(resN)
sigsN<- sigsN[sigs$padj < 0.1,]
dfN <- as.data.frame(sigsN)
dim(dfN)

##Anotación#####################################################################

# lista_atributos<-listAttributes(ensembl)
tabla_anotacionN<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "gene_biotype","description") ,
                       filters = "ensembl_gene_id",
                       values = rownames(dfN),
                       mart = mouse)

dim(tabla_anotacionN)
dim(dfN)
# # Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combinedN <- merge(dfN, tabla_anotacionN %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_finalN<- unique(df_combinedN)
rownames(df_finalN)<-df_finalN$Row.names
df_finalN$Row.names<-NULL
head(df_finalN)

##Seleccion genes Talamo y Corteza###########################################
#Seleccionamos los 400 genes más expresados diferencialmente de Corteza y Talamo
df.topN <- df_finalN[(abs(df_finalN$log2FoldChange) > 0),]
df.topN <- df.topN[order(df.topN$log2FoldChange, decreasing = TRUE),]
#Al ser el numerador el Talamo, los genes upregulados son genes específicos de este. Mientras que los genes downregulados son específicos de la corteza

df_N_TH <- df.topN[1:400, ]
df_N_CX <- df.topN[(nrow(df.topN)-399):nrow(df.topN), ]

##Diagrama de VENN#############################################################

ven_Ns_vs_AS <-list(As_Th=rownames(genesTalamo), 
            Ns_Th=rownames(df_N_TH), 
            As_Ctx = rownames(genesCorteza),
            Ns_Ctx=rownames(df_N_CX))
#library(ggplot2)

diag_venn_Ns_vs_AS=ggVennDiagram(ven_Ns_vs_AS) + scale_fill_gradient(low="blue",high = "red")
png(filename = "DiagVenn_Ns_vs_AS_MULTIMAPPING.png",
    width=1000, height = 1000)
diag_venn_Ns_vs_AS
dev.off()

library(gplots)
#Identificamos los genes que hacen overlap
ven_intersect <- venn(ven_Ns_vs_AS)
genes_overlap<-attr(ven_intersect,"intersections")
genes_overlap_Th<-genes_overlap$`As_Th:Ns_Th`
genes_overlap_Ctx<-genes_overlap$`As_Ctx:Ns_Ctx`
genes_overlap_Th_Ctx<-c(genes_overlap_Th,genes_overlap_Ctx)

#Heatmap de astrocitos, señalando los genes compartidos entre neuronas y astrocitos
##Data Complex Heatmap#########################################################
tabla_complex<-getBM(attributes =c("ensembl_gene_id","mgi_symbol", "gene_biotype","description") ,
                       filters = "ensembl_gene_id",
                       values = rownames(df.top),
                       mart = mouse)

# # Combinar los dos dataframes en función de los nombres de fila coincidentes en el primer dataframe y los valores de la columna "ensembl_gene_id" del segundo dataframe
df_combined_complex <- merge(df.top, tabla_complex %>% select(mgi_symbol,ensembl_gene_id), by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE, all.y=FALSE)

#eliminamos duplicados
df_final_complex<- unique(df_combined_complex)
rownames(df_final_complex)<-df_final_complex$Row.names
df_final_complex$Row.names<-NULL
head(df_final_complex)
##Complex heatmap AS-TH_VS_AS-CTX###############################################
rlog_out2 <- rlog(dds, blind=FALSE) #get normalized count data from dds object
mat2<-assay(rlog_out)[rownames(df_final_complex), rownames(coldata_astrocytes)] #sig genes x samples
base_mean2 <- rowMeans(mat2)
mat.scaled2 <- t(apply(mat2, 1, scale)) #center and scale each column (Z-score) then transpose
genes_anotacion_m<-match(rownames(mat.scaled2), genes_overlap_Th_Ctx)
rownames(mat.scaled2)<- tabla_complex$mgi_symbol 
colnames(mat.scaled2) <- paste0(rld$Region,"-",rld$Nucleo)

anotaciones_h<- rowAnnotation(foo=anno_mark(at = genes_anotacion_m, 
                                            labels = rownames(mat.scaled2)),
                              width = unit(1, "cm") +
                                max_text_width(rownames(mat.scaled2)))


png(filename = "Heatmap_AS-TH_vs_AS-CTX.png",
    width=600, height = 1000)
Heatmap(mat.scaled2, 
        column_title = "AS-TH_vs_AS-CTX",
        column_title_gp = gpar(fill = "white", col = "black", border = "white"),
        right_annotation = anotaciones_h,
        row_names_gp = gpar(fontsize = 0))
dev.off()

#Genes específicos de los astrocitos de diferentes nucleos del talamo#############################################################

##Datos########################################################################
###Astrocitos
cts_As_TH<-cts_astrocytes[,5:17]
col_As_TH<-coldata_astrocytes[5:17,]

###Neuronas
cts_Ns_TH<-cts_neurons[,5:17]
col_Ns_TH<-coldata_neurons[5:17,]

###Total
cts_TH<-cbind(cts_As_TH,cts_Ns_TH)
col_TH<-rbind(col_As_TH,col_Ns_TH)

##Astrocytes############################################################
##Normalizacion#################################################################
ddsTH <- DESeqDataSetFromMatrix(countData = cts_As_TH,
                              colData = col_As_TH,
                              design = ~ Nucleo)
###PREFILTERING AS##################################################################
ddsTH
keepTH <- rowSums(counts(ddsTH)) >= 10
ddsTH <- ddsTH[keepTH,]
###PCA AS###########################################################################
rld_as<-rlog(ddsTH)
data_as<-plotPCA(rld_as, intgroup = "Nucleo",returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#library("ggplot2")
p_as<-qplot(PC1, PC2, color=Nucleo, data=data_as, label=colnames(rld_as)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "OmicosEntregable2/astrocitos/PCA_NucleosTalamo.png",
    width=1000, height = 1000)
p_as + geom_text(aes(label=colnames(rld_as)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Nucleo),                 
                                                                                   alpha = 0.2, 
                                                                                   
                                                                                   show.legend = FALSE, 
                                                                                   
                                                                                   level = 0.75
) 
dev.off()
####dLGN/MGV Resultado#####################################################################
ddsTH$Nucleo <- relevel(ddsTH$Nucleo, "Medial")
ddsTH <- DESeq(ddsTH)
resTH<- results(ddsTH)
#Tabla resultado
#Nucleo Dorsal vs Medial, genes upregulados son especificos del dLGN y los downregulados son específicos del MGV.
summary(resTH)
sigsTH <- na.omit(resTH)
sigsTH<- sigsTH[sigsTH$padj < 0.1,]
dfTH <- as.data.frame(sigsTH)

####spcdLGN/spcMGV#######################################################################
df.topTH <- dfTH[(abs(dfTH$log2FoldChange) > 0.221),]
df.topTH <- df.topTH[order(df.topTH$log2FoldChange, decreasing = TRUE),]

df.topTH$diffexpressed <- "NO"
df.topTH$diffexpressed[df.topTH$log2FoldChange > 0.221] <- "UP"
df.topTH$diffexpressed[df.topTH$log2FoldChange < 0.221] <- "DOWN"
table(df.topTH$diffexpressed)

spcdLGN<-df.topTH[df.topTH$log2FoldChange > 0.221,]
spcMGV<-df.topTH[df.topTH$log2FoldChange < -0.221,]

#Realizamos el mismo proceso para el VPM
###VPM Resultado#####################################################################
ddsVPM <- DESeqDataSetFromMatrix(countData = cts_As_TH,
                                colData = col_As_TH,
                                design = ~ Nucleo)
####PREFILTERING##################################################################
ddsVPM
keepVPM <- rowSums(counts(ddsVPM)) >= 10
ddsVPM <- ddsVPM[keepVPM,]
####VPMResultado##################################################################
ddsVPM$Nucleo <- relevel(ddsVPM$Nucleo, "Medial")
ddsVPM <- DESeq(ddsVPM)
resVPM<- results(ddsVPM)
#Tabla resultado
summary(resVPM)
sigsVPM <- na.omit(resVPM)
sigsVPM<- sigsVPM[sigsVPM$padj < 0.1,]
dfVPM <- as.data.frame(sigsVPM)
dim(dfVPM)

df.topVPM <- dfVPM[(abs(dfVPM$log2FoldChange) > 0.221),]
df.topVPM<- df.topVPM[order(df.topVPM$log2FoldChange, decreasing = TRUE),]
df.topVPM$diffexpressed <- "NO"
df.topVPM$diffexpressed[df.topVPM$log2FoldChange > 0.221] <- "UP"
df.topVPM$diffexpressed[df.topVPM$log2FoldChange < 0.221] <- "DOWN"
table(df.topVPM$diffexpressed)
spcVPM<-df.topVPM[df.topVPM$log2FoldChange > 0.221,]
##Neurons#######################################################################
cts_Ns_TH<-cts_neurons[,7:17]
col_Ns_TH<-coldata_neurons[7:17,]
##Creacion objeto#################################################################
ddsNs_TH <- DESeqDataSetFromMatrix(countData = cts_Ns_TH,
                               colData = col_Ns_TH,
                               design = ~ Nucleo)
##PRE-FILTERING################################################################
ddsNs_TH
keepNs_TH <- rowSums(counts(ddsNs_TH)) >= 10
ddsNs_TH <- ddsNs_TH[keepNs_TH,]
rldNs_Th <- rlog(ddsNs_TH)
##PCA###########################################################################
dataNs_Th<-plotPCA(rldNs_Th, intgroup = "Nucleo",returnData=TRUE)
percentVar <- round(100 * attr(dataNs_Th, "percentVar"))
#library("ggplot2")
p<-qplot(PC1, PC2, color=Nucleo, data=dataNs_Th, label=colnames(rldNs_Th)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
png(filename = "OmicosEntregable2/neurons/PCA_NucleoTh_neurons.png",
    width=1000, height = 1000)
p + geom_text(aes(label=colnames(rldNs_Th)),size=3, vjust=1.5, hjust=0.5)+ stat_ellipse(geom="polygon", aes(fill = Nucleo),                 
                                                                                    alpha = 0.2, 
                                                                                    
                                                                                    show.legend = FALSE, 
                                                                                    
                                                                                    level = 0.75
) 
dev.off()