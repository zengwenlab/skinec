seu <- readRDS("~/Skin_project/seu_skin_all_711.rds")
library(Seurat)
library(CellChat)

options(stringsAsFactors = FALSE)
options(futrue.globlas.Maxsize=2000*1024**3)

seu$celltype1 <- ifelse((seu$celltype %in% c("Epi1","Epi2")), "Epi", 
                        ifelse(seu$celltype %in% c("Fib1","Fib2","Fib3"), "Fib", 
                               ifelse(seu$celltype %in% c("Smc1","Smc2"), "Smc",
                                      ifelse(seu$celltype %in% c("Imm1","Imm2"), "Imm", as.character(seu$celltype)))))
table(seu$celltype1)

seu$celltype1 <- factor(seu$celltype1, levels = c("Epi", "Fib", "Imm", "Smc", "VEc", "Melan", "LEc", "Swan", "Mac", "Gland", "Red", "Mast"))
Idents(seu) <- seu$celltype1
seu <- subset(seu, celltype1!= "Red")
seu$celltype1 <- droplevels(seu$celltype1)
Idents(seu) <- seu$celltype1


seu$seurat_clusters <- NULL
seu$SCT_snn_res.0.45 <- NULL
seu$SCT_snn_res.0.4 <- NULL
seu$SCT_snn_res.0.35 <- NULL
seu$SCT_snn_res.0.3 <- NULL
seu$SCT_snn_res.0.25 <- NULL
seu$SCT_snn_res.0.2 <- NULL
seu$SCT_snn_res.0.15 <- NULL
seu$SCT_snn_res.0.1 <- NULL
seu$SCT_snn_res.0.05 <- NULL
seu$SCT_snn_res.0.05 <- NULL
seu$SCT_snn_res.0 <- NULL
seu$nCount_SCT <- NULL
seu$nFeature_SCT <- NULL
seu$percent.mt <- NULL
seu$G2M.Score <- NULL
seu$S.Score <- NULL
seu$Phase <- NULL
seu@reductions$pca <- NULL
seu@reductions$harmony <- NULL
seu@reductions$umap <- NULL


abdomen <- subset(seu, Location == "abdomen")
armlegs <- subset(seu, Location == "arms&legs")
chest <- subset(seu, Location == "chest")
head <- subset(seu, Location == "head")


###abdomen
data.input.abdomen <- GetAssayData(abdomen, assay = "SCT", slot = "data")
labels <- Idents(abdomen)
identity.abdomen <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.abdomen <- createCellChat(object = data.input.abdomen)
cellchat.seu.abdomen <- addMeta(cellchat.seu.abdomen, meta = identity.abdomen, 
                                meta.name = "labels")
cellchat.seu.abdomen <- setIdent(cellchat.seu.abdomen, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.abdomen@DB <- CellChatDB.use

cellchat.seu.abdomen <- subsetData(cellchat.seu.abdomen)
cellchat.seu.abdomen <- identifyOverExpressedGenes(cellchat.seu.abdomen)
cellchat.seu.abdomen <- identifyOverExpressedInteractions(cellchat.seu.abdomen)
cellchat.seu.abdomen <- projectData(cellchat.seu.abdomen, PPI.human)
cellchat.seu.abdomen <- computeCommunProb(cellchat.seu.abdomen)
cellchat.seu.abdomen <- filterCommunication(cellchat.seu.abdomen, min.cells = 10)
df.net.abdomen <- subsetCommunication(cellchat.seu.abdomen, slot.name = "net")
write.csv(df.net.abdomen, "net.LR.abdomen-sct.csv")


cellchat.seu.abdomen <- computeCommunProbPathway(cellchat.seu.abdomen)
df.netp.abdomen <- subsetCommunication(cellchat.seu.abdomen, slot.name = "netP")
write.csv(df.netp.abdomen, "net.patway.abdomen-sct.csv")

cellchat.seu.abdomen <- aggregateNet(cellchat.seu.abdomen)
cellchat.seu.abdomen <- netAnalysis_computeCentrality(cellchat.seu.abdomen)
saveRDS(cellchat.seu.abdomen, file = "cellchat.abdomen-sct.rds")


####
###armlegs
data.input.armlegs <- GetAssayData(armlegs, assay = "SCT", slot = "data")
labels <- Idents(armlegs)
identity.armlegs <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.armlegs <- createCellChat(object = data.input.armlegs)
cellchat.seu.armlegs <- addMeta(cellchat.seu.armlegs, meta = identity.armlegs, 
                                meta.name = "labels")
cellchat.seu.armlegs <- setIdent(cellchat.seu.armlegs, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.armlegs@DB <- CellChatDB.use

cellchat.seu.armlegs <- subsetData(cellchat.seu.armlegs)
cellchat.seu.armlegs <- identifyOverExpressedGenes(cellchat.seu.armlegs)
cellchat.seu.armlegs <- identifyOverExpressedInteractions(cellchat.seu.armlegs)
cellchat.seu.armlegs <- projectData(cellchat.seu.armlegs, PPI.human)
cellchat.seu.armlegs <- computeCommunProb(cellchat.seu.armlegs)
cellchat.seu.armlegs <- filterCommunication(cellchat.seu.armlegs, min.cells = 10)
df.net.armlegs <- subsetCommunication(cellchat.seu.armlegs, slot.name = "net")
write.csv(df.net.armlegs, "net.LR.armlegs-sct.csv")


cellchat.seu.armlegs <- computeCommunProbPathway(cellchat.seu.armlegs)
df.netp.armlegs <- subsetCommunication(cellchat.seu.armlegs, slot.name = "netP")
write.csv(df.netp.armlegs, "net.patway.armlegs-sct.csv")

cellchat.seu.armlegs <- aggregateNet(cellchat.seu.armlegs)
cellchat.seu.armlegs <- netAnalysis_computeCentrality(cellchat.seu.armlegs)
saveRDS(cellchat.seu.armlegs, file = "cellchat.armlegs-sct.rds")


####
###chest
data.input.chest <- GetAssayData(chest, assay = "SCT", slot = "data")
labels <- Idents(chest)
identity.chest <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.chest <- createCellChat(object = data.input.chest)
cellchat.seu.chest <- addMeta(cellchat.seu.chest, meta = identity.chest, 
                              meta.name = "labels")
cellchat.seu.chest <- setIdent(cellchat.seu.chest, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.chest@DB <- CellChatDB.use

cellchat.seu.chest <- subsetData(cellchat.seu.chest)
cellchat.seu.chest <- identifyOverExpressedGenes(cellchat.seu.chest)
cellchat.seu.chest <- identifyOverExpressedInteractions(cellchat.seu.chest)
cellchat.seu.chest <- projectData(cellchat.seu.chest, PPI.human)
cellchat.seu.chest <- computeCommunProb(cellchat.seu.chest)
cellchat.seu.chest <- filterCommunication(cellchat.seu.chest, min.cells = 10)
df.net.chest <- subsetCommunication(cellchat.seu.chest, slot.name = "net")
write.csv(df.net.chest, "net.LR.chest-sct.csv")


cellchat.seu.chest <- computeCommunProbPathway(cellchat.seu.chest)
df.netp.chest <- subsetCommunication(cellchat.seu.chest, slot.name = "netP")
write.csv(df.netp.chest, "net.patway.chest-sct.csv")

cellchat.seu.chest <- aggregateNet(cellchat.seu.chest)
cellchat.seu.chest <- netAnalysis_computeCentrality(cellchat.seu.chest)
saveRDS(cellchat.seu.chest, file = "cellchat.chest-sct.rds")


####
###head
data.input.head <- GetAssayData(head, assay = "SCT", slot = "data")
labels <- Idents(head)
identity.head <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.head <- createCellChat(object = data.input.head)
cellchat.seu.head <- addMeta(cellchat.seu.head, meta = identity.head, 
                             meta.name = "labels")
cellchat.seu.head <- setIdent(cellchat.seu.head, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.head@DB <- CellChatDB.use

cellchat.seu.head <- subsetData(cellchat.seu.head)
cellchat.seu.head <- identifyOverExpressedGenes(cellchat.seu.head)
cellchat.seu.head <- identifyOverExpressedInteractions(cellchat.seu.head)
cellchat.seu.head <- projectData(cellchat.seu.head, PPI.human)
cellchat.seu.head <- computeCommunProb(cellchat.seu.head)
cellchat.seu.head <- filterCommunication(cellchat.seu.head, min.cells = 10)
df.net.head <- subsetCommunication(cellchat.seu.head, slot.name = "net")
write.csv(df.net.head, "net.LR.head-sct.csv")


cellchat.seu.head <- computeCommunProbPathway(cellchat.seu.head)
df.netp.head <- subsetCommunication(cellchat.seu.head, slot.name = "netP")
write.csv(df.netp.head, "net.patway.head-sct.csv")

cellchat.seu.head <- aggregateNet(cellchat.seu.head)
cellchat.seu.head <- netAnalysis_computeCentrality(cellchat.seu.head)
saveRDS(cellchat.seu.head, file = "cellchat.head-sct.rds")



####
###abdomen
data.input.abdomen <- GetAssayData(abdomen, assay = "RNA", slot = "data")
labels <- Idents(abdomen)
identity.abdomen <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.abdomen <- createCellChat(object = data.input.abdomen)
cellchat.seu.abdomen <- addMeta(cellchat.seu.abdomen, meta = identity.abdomen, 
                                meta.name = "labels")
cellchat.seu.abdomen <- setIdent(cellchat.seu.abdomen, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.abdomen@DB <- CellChatDB.use

cellchat.seu.abdomen <- subsetData(cellchat.seu.abdomen)
cellchat.seu.abdomen <- identifyOverExpressedGenes(cellchat.seu.abdomen)
cellchat.seu.abdomen <- identifyOverExpressedInteractions(cellchat.seu.abdomen)
cellchat.seu.abdomen <- projectData(cellchat.seu.abdomen, PPI.human)
cellchat.seu.abdomen <- computeCommunProb(cellchat.seu.abdomen)
cellchat.seu.abdomen <- filterCommunication(cellchat.seu.abdomen, min.cells = 10)
df.net.abdomen <- subsetCommunication(cellchat.seu.abdomen, slot.name = "net")
write.csv(df.net.abdomen, "net.LR.abdomen-RNA.csv")


cellchat.seu.abdomen <- computeCommunProbPathway(cellchat.seu.abdomen)
df.netp.abdomen <- subsetCommunication(cellchat.seu.abdomen, slot.name = "netP")
write.csv(df.netp.abdomen, "net.patway.abdomen-RNA.csv")

cellchat.seu.abdomen <- aggregateNet(cellchat.seu.abdomen)
cellchat.seu.abdomen <- netAnalysis_computeCentrality(cellchat.seu.abdomen)
saveRDS(cellchat.seu.abdomen, file = "cellchat.abdomen-RNA.rds")


####
###armlegs
data.input.armlegs <- GetAssayData(armlegs, assay = "RNA", slot = "data")
labels <- Idents(armlegs)
identity.armlegs <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.armlegs <- createCellChat(object = data.input.armlegs)
cellchat.seu.armlegs <- addMeta(cellchat.seu.armlegs, meta = identity.armlegs, 
                                meta.name = "labels")
cellchat.seu.armlegs <- setIdent(cellchat.seu.armlegs, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.armlegs@DB <- CellChatDB.use

cellchat.seu.armlegs <- subsetData(cellchat.seu.armlegs)
cellchat.seu.armlegs <- identifyOverExpressedGenes(cellchat.seu.armlegs)
cellchat.seu.armlegs <- identifyOverExpressedInteractions(cellchat.seu.armlegs)
cellchat.seu.armlegs <- projectData(cellchat.seu.armlegs, PPI.human)
cellchat.seu.armlegs <- computeCommunProb(cellchat.seu.armlegs)
cellchat.seu.armlegs <- filterCommunication(cellchat.seu.armlegs, min.cells = 10)
df.net.armlegs <- subsetCommunication(cellchat.seu.armlegs, slot.name = "net")
write.csv(df.net.armlegs, "net.LR.armlegs-RNA.csv")


cellchat.seu.armlegs <- computeCommunProbPathway(cellchat.seu.armlegs)
df.netp.armlegs <- subsetCommunication(cellchat.seu.armlegs, slot.name = "netP")
write.csv(df.netp.armlegs, "net.patway.armlegs-RNA.csv")

cellchat.seu.armlegs <- aggregateNet(cellchat.seu.armlegs)
cellchat.seu.armlegs <- netAnalysis_computeCentrality(cellchat.seu.armlegs)
saveRDS(cellchat.seu.armlegs, file = "cellchat.armlegs-RNA.rds")



####
###chest
data.input.chest <- GetAssayData(chest, assay = "RNA", slot = "data")
labels <- Idents(chest)
identity.chest <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.chest <- createCellChat(object = data.input.chest)
cellchat.seu.chest <- addMeta(cellchat.seu.chest, meta = identity.chest, 
                              meta.name = "labels")
cellchat.seu.chest <- setIdent(cellchat.seu.chest, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.chest@DB <- CellChatDB.use

cellchat.seu.chest <- subsetData(cellchat.seu.chest)
cellchat.seu.chest <- identifyOverExpressedGenes(cellchat.seu.chest)
cellchat.seu.chest <- identifyOverExpressedInteractions(cellchat.seu.chest)
cellchat.seu.chest <- projectData(cellchat.seu.chest, PPI.human)
cellchat.seu.chest <- computeCommunProb(cellchat.seu.chest)
cellchat.seu.chest <- filterCommunication(cellchat.seu.chest, min.cells = 10)
df.net.chest <- subsetCommunication(cellchat.seu.chest, slot.name = "net")
write.csv(df.net.chest, "net.LR.chest-RNA.csv")


cellchat.seu.chest <- computeCommunProbPathway(cellchat.seu.chest)
df.netp.chest <- subsetCommunication(cellchat.seu.chest, slot.name = "netP")
write.csv(df.netp.chest, "net.patway.chest-RNA.csv")

cellchat.seu.chest <- aggregateNet(cellchat.seu.chest)
cellchat.seu.chest <- netAnalysis_computeCentrality(cellchat.seu.chest)
saveRDS(cellchat.seu.chest, file = "cellchat.chest-RNA.rds")


####
###head
data.input.head <- GetAssayData(head, assay = "RNA", slot = "data")
labels <- Idents(head)
identity.head <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.head <- createCellChat(object = data.input.head)
cellchat.seu.head <- addMeta(cellchat.seu.head, meta = identity.head, 
                             meta.name = "labels")
cellchat.seu.head <- setIdent(cellchat.seu.head, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.head@DB <- CellChatDB.use

cellchat.seu.head <- subsetData(cellchat.seu.head)
cellchat.seu.head <- identifyOverExpressedGenes(cellchat.seu.head)
cellchat.seu.head <- identifyOverExpressedInteractions(cellchat.seu.head)
cellchat.seu.head <- projectData(cellchat.seu.head, PPI.human)
cellchat.seu.head <- computeCommunProb(cellchat.seu.head)
cellchat.seu.head <- filterCommunication(cellchat.seu.head, min.cells = 10)
df.net.head <- subsetCommunication(cellchat.seu.head, slot.name = "net")
write.csv(df.net.head, "net.LR.head-RNA.csv")


cellchat.seu.head <- computeCommunProbPathway(cellchat.seu.head)
df.netp.head <- subsetCommunication(cellchat.seu.head, slot.name = "netP")
write.csv(df.netp.head, "net.patway.head-RNA.csv")

cellchat.seu.head <- aggregateNet(cellchat.seu.head)
cellchat.seu.head <- netAnalysis_computeCentrality(cellchat.seu.head)
saveRDS(cellchat.seu.head, file = "cellchat.head-RNA.rds")

cellchat.abdomen <- readRDS("~/Skin_project/cellchat-location/cellchat.abdomen-sct.rds")
cellchat.armlegs <- readRDS("~/Skin_project/cellchat-location/cellchat.armlegs-sct.rds")
cellchat.chest <- readRDS("~/Skin_project/cellchat-location/cellchat.chest-sct.rds")
cellchat.head <- readRDS("~/Skin_project/cellchat-location/cellchat.head-sct.rds")


object.list <- list(abdomen = cellchat.abdomen, 
                    armlegs = cellchat.armlegs,
                    chest = cellchat.chest,
                    head = cellchat.head)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



pdf(file = "abdomen.pdf", height = 10, width = 10)
groupSize_abdomen <- as.numeric(table(cellchat.abdomen@idents))
netVisual_circle(cellchat.abdomen@net$weight, vertex.weight = groupSize_abdomen, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - abdomen", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()


pdf(file = "armlegs.pdf", height = 10, width = 10)
groupSize_armlegs <- as.numeric(table(cellchat.armlegs@idents))
netVisual_circle(cellchat.armlegs@net$weight, vertex.weight = groupSize_armlegs, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - armlegs", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()

pdf(file = "chest.pdf", height = 10, width = 10)
groupSize_chest <- as.numeric(table(cellchat.chest@idents))
netVisual_circle(cellchat.chest@net$weight, vertex.weight = groupSize_chest, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - chest", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()

pdf(file = "head.pdf", height = 10, width = 10)
groupSize_head <- as.numeric(table(cellchat.head@idents))
netVisual_circle(cellchat.head@net$weight, vertex.weight = groupSize_head, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - head", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + 
    colSums(x@net$count)-diag(x@net$count)})

# control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}

pdf(file = "整体互作.pdf", width = 15, height = 10)
patchwork::wrap_plots(plots = gg, nrow = 2)
dev.off()


cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.abdomen@idents))
mat <- cellchat.abdomen@net$count
pdf(file = "abdomen-ec.pdf", height = 15, width = 15)
par(mfrow = c(3,3), xpd=TRUE, mar=c(2,2,2,2))
for (i in 1:nrow(mat)){
  mat2 <- matrix (0, nrow = nrow (mat), ncol = ncol (mat), dimnames = dimnames (mat))
  mat2 [i, ] <- mat[i, ] 
  mat2 <- as.data.frame (mat2) 
  data <- data.frame(mat2 [i, ]) 
  data <- select(data, i, everything())
  datal <- data.frame (mat2[-i, ])
  datal <- select(datal, i, everything ())
  data2 <- rbind(data, datal)
  data2 <- as.matrix(data2)
  netVisual_circle(data2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(data2)[1],
                   color.use = cols_clusters,
                   label.edge = T,
                   edge.label.cex = 1.2,
                   vertex.label.cex = 1.2)
}
dev.off()


cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.chest@idents))
mat <- cellchat.chest@net$count
pdf(file = "chest-ec.pdf", height = 15, width = 15)
par(mfrow = c(3,3), xpd=TRUE, mar=c(2,2,2,2))
for (i in 1:nrow(mat)){
  mat2 <- matrix (0, nrow = nrow (mat), ncol = ncol (mat), dimnames = dimnames (mat))
  mat2 [i, ] <- mat[i, ] 
  mat2 <- as.data.frame (mat2) 
  data <- data.frame(mat2 [i, ]) 
  data <- select(data, i, everything())
  datal <- data.frame (mat2[-i, ])
  datal <- select(datal, i, everything ())
  data2 <- rbind(data, datal)
  data2 <- as.matrix(data2)
  netVisual_circle(data2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(data2)[1],
                   color.use = cols_clusters,
                   label.edge = T,
                   edge.label.cex = 1.2,
                   vertex.label.cex = 1.2)
}
dev.off()


cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.head@idents))
mat <- cellchat.head@net$count
pdf(file = "head-ec.pdf", height = 15, width = 15)
par(mfrow = c(3,3), xpd=TRUE, mar=c(2,2,2,2))
for (i in 1:nrow(mat)){
  mat2 <- matrix (0, nrow = nrow (mat), ncol = ncol (mat), dimnames = dimnames (mat))
  mat2 [i, ] <- mat[i, ] 
  mat2 <- as.data.frame (mat2) 
  data <- data.frame(mat2 [i, ]) 
  data <- select(data, i, everything())
  datal <- data.frame (mat2[-i, ])
  datal <- select(datal, i, everything ())
  data2 <- rbind(data, datal)
  data2 <- as.matrix(data2)
  netVisual_circle(data2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(data2)[1],
                   color.use = cols_clusters,
                   label.edge = T,
                   edge.label.cex = 1.2,
                   vertex.label.cex = 1.2)
}
dev.off()

pdf(file = "内皮调控成其他细胞.pdf", height = 10, width = 13)
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1, 2, 3, 4, 5, 6, 
                                                               7, 8, 9, 10, 11),  
                 comparison = c(1, 2, 3, 4), angle.x = 45)
dev.off()

pathways.show <- "CXCL"
pdf(file = "abdomen-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.abdomen, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "armlegs-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.armlegs, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "chest-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.chest, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "head-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.head, signaling = pathways.show, layout = "chord")
dev.off()


pathways.show <- "VEGF"
pdf(file = "abdomen-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.abdomen, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "armlegs-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.armlegs, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "chest-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.chest, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "head-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.head, signaling = pathways.show, layout = "chord")
dev.off()


library(NMF)
library(ggalluvial)
pdf(file = "abdomen-out.pdf", width = 10, height = 5)
selectK(cellchat.abdomen, pattern = "outgoing")
dev.off()


pdf(file = "abdomen-in.pdf", width = 10, height = 5)
selectK(cellchat.abdomen, pattern = "incoming")
dev.off()


pdf(file = "armlegs-out.pdf", width = 10, height = 5)
selectK(cellchat.armlegs, pattern = "outgoing")
dev.off()


pdf(file = "armlegs-in.pdf", width = 10, height = 5)
selectK(cellchat.armlegs, pattern = "incoming")
dev.off()


pdf(file = "chest-out.pdf", width = 10, height = 5)
selectK(cellchat.chest, pattern = "outgoing")
dev.off()


pdf(file = "chest-in.pdf", width = 10, height = 5)
selectK(cellchat.chest, pattern = "incoming")
dev.off()


pdf(file = "head-out.pdf", width = 10, height = 5)
selectK(cellchat.head, pattern = "outgoing")
dev.off()


pdf(file = "head-in.pdf", width = 10, height = 5)
selectK(cellchat.head, pattern = "incoming")
dev.off()

pdf(file = "abdomen-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.abdomen, pattern = "outgoing", k = 4, 
                              height = 10)
dev.off()


pdf(file = "abdomen-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.abdomen, pattern = "incoming", k = 4, 
                              height = 10)
dev.off()



pdf(file = "armlegs-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.armlegs, pattern = "outgoing", k = 2, 
                              height = 10)
dev.off()

pdf(file = "armlegs-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.armlegs, pattern = "incoming", k = 2, 
                              height = 10)
dev.off()

pdf(file = "chest-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.chest, pattern = "outgoing", k = 2, 
                              height = 12)
dev.off()

pdf(file = "chest-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.chest, pattern = "incoming", k = 2, 
                              height = 12)
dev.off()


pdf(file = "head-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.head, pattern = "outgoing", k = 2, 
                              height = 12)
dev.off()

pdf(file = "head-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.head, pattern = "incoming", k = 2, 
                              height = 12)
dev.off()

# river plot
cellchat.abdomen <- identifyCommunicationPatterns(object = cellchat.abdomen, pattern = "outgoing", k = 4, height = 10)
pdf(file = "abdomen.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.abdomen, pattern = "outgoing")
dev.off()


cellchat.abdomen <- identifyCommunicationPatterns(object = cellchat.abdomen, pattern = "incoming", k = 4, height = 10)
pdf(file = "abdomen.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.abdomen, pattern = "incoming")
dev.off()

cellchat.armlegs <- identifyCommunicationPatterns(object = cellchat.armlegs, pattern = "outgoing", k = 4, height = 10)
pdf(file = "armlegs.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.armlegs, pattern = "outgoing")
dev.off()


cellchat.armlegs <- identifyCommunicationPatterns(object = cellchat.armlegs, pattern = "incoming", k = 4, height = 10)
pdf(file = "armlegs.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.armlegs, pattern = "incoming")
dev.off()


cellchat.chest <- identifyCommunicationPatterns(object = cellchat.chest, pattern = "outgoing", k = 4, height = 10)
pdf(file = "chest.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.chest, pattern = "outgoing")
dev.off()


cellchat.chest <- identifyCommunicationPatterns(object = cellchat.chest, pattern = "incoming", k = 4, height = 10)
pdf(file = "chest.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.chest, pattern = "incoming")
dev.off()


cellchat.head <- identifyCommunicationPatterns(object = cellchat.head, pattern = "outgoing", k = 4, height = 10)
pdf(file = "head.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.head, pattern = "outgoing")
dev.off()


cellchat.head <- identifyCommunicationPatterns(object = cellchat.head, pattern = "incoming", k = 4, height = 10)
pdf(file = "head.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.head, pattern = "incoming")
dev.off()













