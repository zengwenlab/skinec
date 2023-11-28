
##function

DoubletFinder_multisamples <- function(seu_list){
  library(tidyverse)
  library(Seurat)
  library(DoubletFinder)
  library(conflicted)
  conflict_prefer("select", "dplyr")
  for (idx in 1:length(seu_list)) {
    data <- NormalizeData(seu_list[[idx]])
    data <- ScaleData(data, verbose = FALSE)
    data <- FindVariableFeatures(data, verbose = FALSE)
    data <- RunPCA(data, npcs = 40, verbose = FALSE)
    data <- RunUMAP(data, reduction = "pca", dims = 1:30)
    data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
    data <- FindClusters(data, resolution = 0.5)
    sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    ## Homotypic Doublet Proportion Estimate 
    annotations <- data@meta.data$ClusteringResults
    homotypic.prop <- modelHomotypic(annotations)   
    #The sequencing data of the 10X platform is calculated based on an increase of 0.8/100 in the doublet ratio for every 1000 cells added
    nExp_poi <- round((ncol(data)*8*1e-6)*nrow(data@meta.data)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies 
    data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = 0.09, 
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    ## save results
    seu_list[[idx]]$doubFind_res = data@meta.data %>% select(contains('DF.classifications'))
    seu_list[[idx]]$doubFind_score = data@meta.data %>% select(contains('pANN'))
  }
  seu_list
}

########

library(Seurat)


##############1, Direder2022#########################################################

counts_skin1 <- Read10X("~/zhouxin/new/rawdata/Direder2022/skin_7")
counts_keloid1 <- Read10X("~/zhouxin/new/rawdata/Direder2022/keloid1")
counts_keloid2 <- Read10X("~/zhouxin/new/rawdata/Direder2022/keloid2")
counts_keloid3L <- Read10X("~/zhouxin/new/rawdata/Direder2022/keloid3L")
counts_keloid3R <- Read10X("~/zhouxin/new/rawdata/Direder2022/keloid3R")
counts_scar1 <- Read10X("~/zhouxin/new/rawdata/Direder2022/scar1")
counts_scar2 <- Read10X("~/zhouxin/new/rawdata/Direder2022/scar2")
counts_scar3 <- Read10X("~/zhouxin/new/rawdata/Direder2022/scar3")

seu1 <- CreateSeuratObject(counts_skin1, project = "skin1", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_keloid1, project = "keloid1", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_keloid2, project = "keloid2", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_keloid3L, project = "keloid3", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_keloid3R, project = "keloid4", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_scar1, project = "scar1", min.cells = 3, min.features = 200)
seu7 <- CreateSeuratObject(counts_scar2, project = "scar2", min.cells = 3, min.features = 200)
seu8 <- CreateSeuratObject(counts_scar3, project = "scar3", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6, seu7, seu8), add.cell.ids=c("skin1", "keloid1", "keloid2", "keloid3", "keloid4", "scar1", "scar2", "scar3"))

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Direder2022/Direder2022_1.rds")

seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Direder2022/Direder2022_2.rds")



seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin1, c(seu_list$scar3, seu_list$scar2, 
                               seu_list$scar1, seu_list$keloid3, seu_list$keloid4, 
                               seu_list$keloid2, seu_list$keloid1))

saveRDS(seu, file = "~/zhouxin/new/raw_rds/Direder2022/Direder2022_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Direder2022/Direder2022_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Direder2022/Direder2022_5.rds")


###################################2, Tabib2018#############################################

counts_skin2 <- Read10X("~/zhouxin/new/rawdata/Tabib2018/SC1control/filtered_feature_bc_matrix")
counts_skin3 <- Read10X("~/zhouxin/new/rawdata/Tabib2018/SC4control/filtered_feature_bc_matrix")
counts_skin4 <- Read10X("~/zhouxin/new/rawdata/Tabib2018/SC18control/filtered_feature_bc_matrix")
counts_skin5 <- Read10X("~/zhouxin/new/rawdata/Tabib2018/SC32control/filtered_feature_bc_matrix")
counts_skin6 <- Read10X("~/zhouxin/new/rawdata/Tabib2018/SC33control/filtered_feature_bc_matrix")
counts_skin7 <- Read10X("~/zhouxin/new/rawdata/Tabib2018/SC34control/filtered_feature_bc_matrix")

seu1 <- CreateSeuratObject(counts_skin2, project = "skin2", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin3, project = "skin3", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin4, project = "skin4", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin5, project = "skin5", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_skin6, project = "skin6", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_skin7, project = "skin7", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6), add.cell.ids=c("skin2", "skin3", "skin4", "skin5", "skin6", "skin7"))

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Tabib2018/Tabib2018_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)


seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Tabib2018/Tabib2018_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin1, c(seu_list$skin2, seu_list$skin3, 
                               seu_list$skin4, seu_list$skin5, seu_list$skin6, 
                               seu_list$skin7))
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Tabib2018_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/rawdata/raw_rds/Tabib2018_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Tabib2018_5.rds")


###########################3, Vorstandlechner2021##############################################################

counts_skin8 <- Read10X("~/zhouxin/new/rawdata/Vorstandlechner2021/skin_1")
counts_skin9 <- Read10X("~/zhouxin/new/rawdata/Vorstandlechner2021/skin_2")
counts_skin10 <- Read10X("~/zhouxin/new/rawdata/Vorstandlechner2021/skin_3")
counts_scar4 <- Read10X("~/zhouxin/new/rawdata/Vorstandlechner2021/scar_1")
counts_scar5 <- Read10X("~/zhouxin/new/rawdata/Vorstandlechner2021/scar_2")
counts_scar6 <- Read10X("~/zhouxin/new/rawdata/Vorstandlechner2021/scar_3")

seu1 <- CreateSeuratObject(counts_skin8, project = "skin8", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin9, project = "skin9", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin10, project = "skin10", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_scar4, project = "scar4", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_scar5, project = "scar5", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_scar6, project = "scar6", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6), add.cell.ids=c("skin8", "skin9", "skin10", "scar4", "scar5", "scar6"))


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Vorstandlechner2021/Vorstandlechner2021_1.rds")

seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Vorstandlechner2021/Vorstandlechner2021_2.rds")


seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin8, c(seu_list$skin9, 
                               seu_list$skin10, seu_list$scar4, seu_list$scar5, 
                               seu_list$scar6))
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Vorstandlechner2021/Vorstandlechner2021_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Vorstandlechner2021/Vorstandlechner2021_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Vorstandlechner2021/Vorstandlechner2021_5.rds")


#############################4, Shim2022_GSE181297#################################################################

counts_keloid5 <- Read10X("~/zhouxin/new/rawdata/Shim2022/ke01")
counts_keloid6 <- Read10X("~/zhouxin/new/rawdata/Shim2022/ke02")
counts_scar7 <- Read10X("~/zhouxin/new/rawdata/Shim2022/ns02")

seu1 <- CreateSeuratObject(counts_keloid5, project = "keloid5", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_keloid6, project = "keloid6", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_scar7, project = "scar7", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3), add.cell.ids=c("keloid5", "keloid6", "scar7"))

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Shim2022/Shim2022_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Shim2022/Shim2022_2.rds")


seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$keloid5, c(seu_list$keloid6, seu_list$scar7))

saveRDS(seu, file = "~/zhouxin/new/raw_rds/Shim2022/Shim2022_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Shim2022/Shim2022_4.rds")
seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Shim2022/Shim2022_5.rds")

####################################5, SoleBoldo2020_GSE130973###################################################

counts_skin11 <- Read10X("~/zhouxin/new/rawdata/SoleBoldo2020/PRJNA542149/SRR9036396Y1/filtered_feature_bc_matrix")
counts_skin12 <- Read10X("~/zhouxin/new/rawdata/SoleBoldo2020/PRJNA542149/SRR9036397Y2/filtered_feature_bc_matrix")
counts_skin13 <- Read10X("~/zhouxin/new/rawdata/SoleBoldo2020/PRJNA542149/SRR9036398O1/filtered_feature_bc_matrix")
counts_skin14 <- Read10X("~/zhouxin/new/rawdata/SoleBoldo2020/PRJNA542149/SRR9036399O2/filtered_feature_bc_matrix")
counts_skin15 <- Read10X("~/zhouxin/new/rawdata/SoleBoldo2020/PRJNA542149/SRR9036400O3/filtered_feature_bc_matrix")


seu1 <- CreateSeuratObject(counts_skin11, project = "skin11", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin12, project = "skin12", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin13, project = "skin13", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin14, project = "skin14", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_skin15, project = "skin15", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5), add.cell.ids=c("skin11", "skin12", "skin13", "skin14", "skin15"))

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/SoleBoldo2020/SoleBoldo2020_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 10)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/SoleBoldo2020/SoleBoldo2020_2.rds")


seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin11, c(seu_list$skin12, seu_list$skin13, 
                                seu_list$skin14, seu_list$skin15))
saveRDS(seu, file = "~/zhouxin/new/raw_rds/SoleBoldo2020/SoleBoldo2020_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/SoleBoldo2020/SoleBoldo2020_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/SoleBoldo2020/SoleBoldo2020_5.rds")


####################################6, Rustagi2022_GSE182208#####################################################

counts_skin16 <- Read10X("~/zhouxin/new/rawdata/Rustagi2022/SC_NS-1")
counts_skin17 <- Read10X("~/zhouxin/new/rawdata/Rustagi2022/SC_NS-2")
counts_skin18 <- Read10X("~/zhouxin/new/rawdata/Rustagi2022/SC_NS-6")
counts_skin19 <- Read10X("~/zhouxin/new/rawdata/Rustagi2022/Sc_Skin3")
counts_skin20 <- Read10X("~/zhouxin/new/rawdata/Rustagi2022/Sc_Skin4")

seu1 <- CreateSeuratObject(counts_skin16, project = "skin16", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin17, project = "skin17", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin18, project = "skin18", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin19, project = "skin19", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_skin20, project = "skin20", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5), add.cell.ids=c("skin16", "skin17", "skin18", "skin19", "skin20"))


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Rustagi2022/Rustagi2022_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Rustagi2022/Rustagi2022_2.rds")


seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin16, c(seu_list$skin17, seu_list$skin18, 
                                seu_list$skin19, seu_list$skin20))
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Rustagi2022/Rustagi2022_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Rustagi2022/Rustagi2022_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Rustagi2022/Rustagi2022_5.rds")


####################################7, Gao2021##################################################################

library(loomR)
hc.2 <- connect('~/zhouxin/new/rawdata/Gao2021/GSE162183_Raw_gene_counts_matrix_LoomFile.loom', mode = 'r+', skip.validate = TRUE)

gse162183_matrix <- hc.2[["matrix"]][,]
gse162183_matrix <- t(gse162183_matrix)

gene=hc.2$row.attrs$Gene[]
barcode=hc.2$col.attrs$CellID[]
colnames(gse162183_matrix)= barcode
row.names(gse162183_matrix)= gene

seu <- CreateSeuratObject(counts = gse162183_matrix, min.cells = 3, min.features = 200)

table(seu$orig.ident)
seu <- subset(seu, orig.ident %in% c("Ctrl1", "Ctrl2", "Ctrl3"))
seu$orig.ident <- factor(seu$orig.ident, levels = c("Ctrl1", "Ctrl2", "Ctrl3"))

seu <- RenameCells(seu, new.names = gsub(pattern = "Ctrl1", replacement = "skin21", colnames(seu)))
seu <- RenameCells(seu, new.names = gsub(pattern = "Ctrl2", replacement = "skin22", colnames(seu)))
seu <- RenameCells(seu, new.names = gsub(pattern = "Ctrl3", replacement = "skin23", colnames(seu)))

table(seu$orig.ident)

seu$orig.ident <- ifelse(str_detect(colnames(seu),"skin21"), "skin21", ifelse(str_detect(colnames(seu),"skin22"), "skin22", "skin23"))

table(seu$orig.ident)
seu$orig.ident <- factor(seu$orig.ident, levels = c("skin21", "skin22", "skin23"))
Idents(seu) <- "orig.ident"
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Gao2021/Gao2021_1.rds")


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)


seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Gao2021/Gao2021_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin21, c(seu_list$skin22, seu_list$skin23))
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Gao2021/Gao2021_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Gao2021/Gao2021_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Gao2021/Gao2021_5.rds")


###################################8, Wang2021################################################################

count_skin24 <- read.table(gzfile("~/zhouxin/new/rawdata/Wang2021/GSM4815803_counts_scHealthy1.txt.gz"))
count_skin25 <- read.table(gzfile("~/zhouxin/new/rawdata/Wang2021/GSM4815804_counts_scHealthy2.txt.gz"))
count_skin26 <- read.table(gzfile("~/zhouxin/new/rawdata/Wang2021/GSM4815805_counts_scHealthy3.txt.gz"))

seu1 <- CreateSeuratObject(count_skin24, project = "skin24", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(count_skin25, project = "skin25", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(count_skin26, project = "skin26", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3), add.cell.ids=c("skin24", "skin25", "skin26"))
table(Idents(seu))
seu@active.ident
table(seu$orig.ident)

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Wang2021/Wang2021_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Wang2021/Wang2021_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin24, c(seu_list$skin25, seu_list$skin26))
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Wang2021/Wang2021_3.rds")
saveRDS(seu_list, file = file = "~/zhouxin/new/raw_rds/Wang2021/Wang2021_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = file = "~/zhouxin/new/raw_rds/Wang2021/Wang2021_5.rds")

####################################9, Deng2021################################################################

counts_keloid7 <- Read10X("~/zhouxin/new/rawdata/Deng2021/KL1")
counts_keloid8 <- Read10X("~/zhouxin/new/rawdata/Deng2021/KL2")
counts_keloid9 <- Read10X("~/zhouxin/new/rawdata/Deng2021/KL3")
counts_scar8 <- Read10X("~/zhouxin/new/rawdata/Deng2021/NS1")
counts_scar9 <- Read10X("~/zhouxin/new/rawdata/Deng2021/NS2")
counts_scar10 <- Read10X("~/zhouxin/new/rawdata/Deng2021/NS3")

seu1 <- CreateSeuratObject(counts_keloid7, project = "keloid7", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_keloid8, project = "keloid8", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_keloid9, project = "keloid9", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_scar8, project = "scar8", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_scar9, project = "scar9", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_scar10, project = "scar10", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6), add.cell.ids=c("keloid7", "keloid8", "keloid9", "scar8", "scar9", "scar10"))
table(Idents(seu))
seu@active.ident
table(seu$orig.ident)


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Deng2021/Deng2021_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Deng2021/Deng2021_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$keloid7, c(seu_list$keloid8, seu_list$keloid9, seu_list$scar8, seu_list$scar9, seu_list$scar10))
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Deng2021/Deng2021_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Deng2021/Deng2021_4.rds")
seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Deng2021/Deng2021_5.rds")



####################################10, Zou2021################################################################
counts_skin27 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077736Y18/filtered_feature_bc_matrix")
counts_skin28 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077737Y22/filtered_feature_bc_matrix")
counts_skin29 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077738Y23/filtered_feature_bc_matrix")
counts_skin30 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077739M44/filtered_feature_bc_matrix")
counts_skin31 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077740M47/filtered_feature_bc_matrix")
counts_skin32 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077741M48/filtered_feature_bc_matrix")
counts_skin33 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077742O70/filtered_feature_bc_matrix")
counts_skin34 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077743O73/filtered_feature_bc_matrix")
counts_skin35 <- Read10X("~/zhouxin/new/rawdata/Zou2021/HRA000395/HRI077744O76/filtered_feature_bc_matrix")

seu1 <- CreateSeuratObject(counts_skin27, project = "skin27", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin28, project = "skin28", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin29, project = "skin29", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin30, project = "skin30", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_skin31, project = "skin31", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_skin32, project = "skin32", min.cells = 3, min.features = 200)
seu7 <- CreateSeuratObject(counts_skin33, project = "skin33", min.cells = 3, min.features = 200)
seu8 <- CreateSeuratObject(counts_skin34, project = "skin34", min.cells = 3, min.features = 200)
seu9 <- CreateSeuratObject(counts_skin35, project = "skin35", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6, seu7, seu8, seu9), add.cell.ids=c("skin27", "skin28", "skin29", "skin30", 
                                                                                     "skin31", "skin32", "skin33", "skin34",
                                                                                     "skin35"))
table(Idents(seu))
seu@active.ident
table(seu$orig.ident)


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Zou2021/Zou2021_1.rds")
seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Zou2021/Zou2021_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin27, c(seu_list$skin28, seu_list$skin29, 
                                seu_list$skin30, seu_list$skin31, seu_list$skin32, 
                                seu_list$skin33, seu_list$skin34, seu_list$skin35))

saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Zou2021_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/rawdata/raw_rds/Zou2021_4.rds")
seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Zou2021_5.rds")

####################################11, Ahlers2021#############################################################
counts_skin36 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S2/filtered_feature_bc_matrix")
counts_skin37 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S4/filtered_feature_bc_matrix")
counts_skin38 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S5/filtered_feature_bc_matrix")
counts_skin39 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S6/filtered_feature_bc_matrix")
counts_skin40 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S7/filtered_feature_bc_matrix")
counts_skin41 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S9/filtered_feature_bc_matrix")
counts_skin42 <- Read10X("~/zhouxin/new/rawdata/Ahlers2021/S12/filtered_feature_bc_matrix")

seu1 <- CreateSeuratObject(counts_skin36, project = "skin36", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin37, project = "skin37", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin38, project = "skin38", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin39, project = "skin39", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_skin40, project = "skin40", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_skin41, project = "skin41", min.cells = 3, min.features = 200)
seu7 <- CreateSeuratObject(counts_skin42, project = "skin42", min.cells = 3, min.features = 200)


seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6, seu7), add.cell.ids=c("skin36", "skin37", "skin38", "skin39", "skin40", "skin41", "skin42"))
table(Idents(seu))
seu@active.ident
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Ahlers2021_seu.rds")

seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Ahlers2021_1.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)

seu <- merge(seu_list$skin36, c(seu_list$skin37, seu_list$skin38, 
                                seu_list$skin39, seu_list$skin40, seu_list$skin41, 
                                seu_list$skin42))

saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Ahlers2021_2.rds")
saveRDS(seu_list, file = "~/zhouxin/new/rawdata/raw_rds/Ahlers2021_3.rds")
seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Ahlers2021_4.rds")

####################################12, Liu2022#####################################################

counts_keloid10 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112758K007CASE/filtered_feature_bc_matrix")
counts_skin43 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112759K007CTRL/filtered_feature_bc_matrix")
counts_keloid11 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112760K009CASE/filtered_feature_bc_matrix")
counts_skin44 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112761K009CTRL/filtered_feature_bc_matrix")
counts_keloid12 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112762K012CASE/filtered_feature_bc_matrix")
counts_skin45 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112763K012CTRL/filtered_feature_bc_matrix")
counts_keloid13 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112764K013CASE/filtered_feature_bc_matrix")
counts_skin46 <- Read10X("~/zhouxin/new/rawdata/Liu2022/HRA000425/HRS112765K013CTRL/filtered_feature_bc_matrix")


seu1 <- CreateSeuratObject(counts_keloid10, project = "keloid10", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin43, project = "skin43", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_keloid11, project = "keloid11", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin44, project = "skin44", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts_keloid12, project = "keloid12", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts_skin45, project = "skin45", min.cells = 3, min.features = 200)
seu7 <- CreateSeuratObject(counts_keloid13, project = "keloid13", min.cells = 3, min.features = 200)
seu8 <- CreateSeuratObject(counts_skin46, project = "skin46", min.cells = 3, min.features = 200)


seu <- merge(seu1, c(seu2, seu3, seu4, seu5, seu6, seu7, seu8), add.cell.ids=c("keloid10", "skin43", "keloid11", "skin44", 
                                                                               "keloid12", "skin45", "skin13", "skin46"))
table(Idents(seu))
seu@active.ident
table(seu$orig.ident)


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Liu2022/Liu2022_1.rds")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Liu2022/Liu2022_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$keloid10, c(seu_list$skin43, seu_list$keloid11, 
                                  seu_list$skin44, seu_list$keloid12, seu_list$skin45, 
                                  seu_list$keloid13, seu_list$skin46))

saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Liu2022_2.rds")
saveRDS(seu_list, file = "~/zhouxin/new/rawdata/raw_rds/Liu2022_3.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/rawdata/raw_rds/Liu2022_4.rds")


####################################13, Xue2022#####################################################

counts_skin47 <- Read10X_h5(filename = "~/zhouxin/new/rawdata/Xue2022/GSM4115878_SC50raw_feature_bc_matrix.h5")
counts_skin48 <- Read10X_h5(filename = "~/zhouxin/new/rawdata/Xue2022/GSM4115880_SC68raw_feature_bc_matrix.h5")
counts_skin49 <- Read10X_h5(filename = "~/zhouxin/new/rawdata/Xue2022/GSM4115885_SC124raw_feature_bc_matrix.h5")
counts_skin50 <- Read10X_h5(filename = "~/zhouxin/new/rawdata/Xue2022/GSM4115886_SC125raw_feature_bc_matrix.h5")

seu1 <- CreateSeuratObject(counts_skin47, project = "skin47", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts_skin48, project = "skin48", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts_skin49, project = "skin49", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts_skin50, project = "skin50", min.cells = 3, min.features = 200)

seu <- merge(seu1, c(seu2, seu3, seu4), add.cell.ids=c("skin47", "skin48", "skin49", "skin50"))
table(Idents(seu))
seu@active.ident
table(seu$orig.ident)


seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Xue2022/Xue2022_1.rds")
seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 15)
table(seu$orig.ident)
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Xue2022/Xue2022_2.rds")

seu_list <- SplitObject(seu, split.by = "orig.ident")
seu_list <- DoubletFinder_multisamples(seu_list = seu_list)
seu <- merge(seu_list$skin47, c(seu_list$skin48, seu_list$skin49, seu_list$skin50))

saveRDS(seu, file = "~/zhouxin/new/raw_rds/Xue2022/Xue2022_3.rds")
saveRDS(seu_list, file = "~/zhouxin/new/raw_rds/Xue2022/Xue2022_4.rds")

seu <- subset(seu, doubFind_res == "Singlet")
saveRDS(seu, file = "~/zhouxin/new/raw_rds/Xue2022/Xue2022_5.rds")

#############################################

