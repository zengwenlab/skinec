seu <- readRDS("~/seu_skin_all_711.rds")

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

African <- subset(seu, Race == "African")
Asian <- subset(seu, Race == "Asian")
Caucasian <- subset(seu, Race == "Caucasian")


####African
###African
data.input.African <- GetAssayData(African, assay = "RNA", slot = "data")
labels <- Idents(African)
identity.African <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.African <- createCellChat(object = data.input.African)
cellchat.seu.African <- addMeta(cellchat.seu.African, meta = identity.African, 
                                meta.name = "labels")
cellchat.seu.African <- setIdent(cellchat.seu.African, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.African@DB <- CellChatDB.use

cellchat.seu.African <- subsetData(cellchat.seu.African)
cellchat.seu.African <- identifyOverExpressedGenes(cellchat.seu.African)
cellchat.seu.African <- identifyOverExpressedInteractions(cellchat.seu.African)
cellchat.seu.African <- projectData(cellchat.seu.African, PPI.human)
cellchat.seu.African <- computeCommunProb(cellchat.seu.African)
cellchat.seu.African <- filterCommunication(cellchat.seu.African, min.cells = 10)
df.net.African <- subsetCommunication(cellchat.seu.African)
write.csv(df.net.African, "net.LR.African-RNA.csv")

cellchat.seu.African <- computeCommunProbPathway(cellchat.seu.African)
df.netp.African <- subsetCommunication(cellchat.seu.African, slot.name = "netP")
write.csv(df.netp.African, "net.patway.African-RNA.csv")

cellchat.seu.African <- aggregateNet(cellchat.seu.African)
cellchat.seu.African <- netAnalysis_computeCentrality(cellchat.seu.African)
saveRDS(cellchat.seu.African, file = "cellchat.African-RNA.rds")


####Asian
###Asian
data.input.Asian <- GetAssayData(Asian, assay = "RNA", slot = "data")
labels <- Idents(Asian)
identity.Asian <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.Asian <- createCellChat(object = data.input.Asian)
cellchat.seu.Asian <- addMeta(cellchat.seu.Asian, meta = identity.Asian, 
                              meta.name = "labels")
cellchat.seu.Asian <- setIdent(cellchat.seu.Asian, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.Asian@DB <- CellChatDB.use

cellchat.seu.Asian <- subsetData(cellchat.seu.Asian)
cellchat.seu.Asian <- identifyOverExpressedGenes(cellchat.seu.Asian)
cellchat.seu.Asian <- identifyOverExpressedInteractions(cellchat.seu.Asian)
cellchat.seu.Asian <- projectData(cellchat.seu.Asian, PPI.human)
cellchat.seu.Asian <- computeCommunProb(cellchat.seu.Asian)
cellchat.seu.Asian <- filterCommunication(cellchat.seu.Asian, min.cells = 10)
df.net.Asian <- subsetCommunication(cellchat.seu.Asian)
write.csv(df.net.Asian, "net.LR.Asian-RNA.csv")

cellchat.seu.Asian <- computeCommunProbPathway(cellchat.seu.Asian)
df.netp.Asian <- subsetCommunication(cellchat.seu.Asian, slot.name = "netP")
write.csv(df.netp.Asian, "net.patway.Asian-RNA.csv")

cellchat.seu.Asian <- aggregateNet(cellchat.seu.Asian)
cellchat.seu.Asian <- netAnalysis_computeCentrality(cellchat.seu.Asian)
saveRDS(cellchat.seu.Asian, file = "cellchat.Asian-RNA.rds")



####Caucasian
###Caucasian
data.input.Caucasian <- GetAssayData(Caucasian, assay = "RNA", slot = "data")
labels <- Idents(Caucasian)
identity.Caucasian <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.Caucasian <- createCellChat(object = data.input.Caucasian)
cellchat.seu.Caucasian <- addMeta(cellchat.seu.Caucasian, meta = identity.Caucasian, 
                                  meta.name = "labels")
cellchat.seu.Caucasian <- setIdent(cellchat.seu.Caucasian, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.Caucasian@DB <- CellChatDB.use

cellchat.seu.Caucasian <- subsetData(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- identifyOverExpressedGenes(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- identifyOverExpressedInteractions(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- projectData(cellchat.seu.Caucasian, PPI.human)
cellchat.seu.Caucasian <- computeCommunProb(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- filterCommunication(cellchat.seu.Caucasian, min.cells = 10)
df.net.Caucasian <- subsetCommunication(cellchat.seu.Caucasian)
write.csv(df.net.Caucasian, "net.LR.Caucasian-RNA.csv")

cellchat.seu.Caucasian <- computeCommunProbPathway(cellchat.seu.Caucasian)
df.netp.Caucasian <- subsetCommunication(cellchat.seu.Caucasian, slot.name = "netP")
write.csv(df.netp.Caucasian, "net.patway.Caucasian-RNA.csv")

cellchat.seu.Caucasian <- aggregateNet(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- netAnalysis_computeCentrality(cellchat.seu.Caucasian)
saveRDS(cellchat.seu.Caucasian, file = "cellchat.Caucasian-RNA.rds")


####African
###African
data.input.African <- GetAssayData(African, assay = "SCT", slot = "data")
labels <- Idents(African)
identity.African <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.African <- createCellChat(object = data.input.African)
cellchat.seu.African <- addMeta(cellchat.seu.African, meta = identity.African, 
                                meta.name = "labels")
cellchat.seu.African <- setIdent(cellchat.seu.African, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.African@DB <- CellChatDB.use

cellchat.seu.African <- subsetData(cellchat.seu.African)
cellchat.seu.African <- identifyOverExpressedGenes(cellchat.seu.African)
cellchat.seu.African <- identifyOverExpressedInteractions(cellchat.seu.African)
cellchat.seu.African <- projectData(cellchat.seu.African, PPI.human)
cellchat.seu.African <- computeCommunProb(cellchat.seu.African)
cellchat.seu.African <- filterCommunication(cellchat.seu.African, min.cells = 10)
df.net.African <- subsetCommunication(cellchat.seu.African)
write.csv(df.net.African, "net.LR.African-sct.csv")

cellchat.seu.African <- computeCommunProbPathway(cellchat.seu.African)
df.netp.African <- subsetCommunication(cellchat.seu.African, slot.name = "netP")
write.csv(df.netp.African, "net.patway.African-sct.csv")

cellchat.seu.African <- aggregateNet(cellchat.seu.African)
cellchat.seu.African <- netAnalysis_computeCentrality(cellchat.seu.African)
saveRDS(cellchat.seu.African, file = "cellchat.African-sct.rds")



####Asian
###Asian
data.input.Asian <- GetAssayData(Asian, assay = "SCT", slot = "data")
labels <- Idents(Asian)
identity.Asian <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.Asian <- createCellChat(object = data.input.Asian)
cellchat.seu.Asian <- addMeta(cellchat.seu.Asian, meta = identity.Asian, 
                              meta.name = "labels")
cellchat.seu.Asian <- setIdent(cellchat.seu.Asian, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.Asian@DB <- CellChatDB.use

cellchat.seu.Asian <- subsetData(cellchat.seu.Asian)
cellchat.seu.Asian <- identifyOverExpressedGenes(cellchat.seu.Asian)
cellchat.seu.Asian <- identifyOverExpressedInteractions(cellchat.seu.Asian)
cellchat.seu.Asian <- projectData(cellchat.seu.Asian, PPI.human)
cellchat.seu.Asian <- computeCommunProb(cellchat.seu.Asian)
cellchat.seu.Asian <- filterCommunication(cellchat.seu.Asian, min.cells = 10)
df.net.Asian <- subsetCommunication(cellchat.seu.Asian)
write.csv(df.net.Asian, "net.LR.Asian-sct.csv")

cellchat.seu.Asian <- computeCommunProbPathway(cellchat.seu.Asian)
df.netp.Asian <- subsetCommunication(cellchat.seu.Asian, slot.name = "netP")
write.csv(df.netp.Asian, "net.patway.Asian-sct.csv")

cellchat.seu.Asian <- aggregateNet(cellchat.seu.Asian)
cellchat.seu.Asian <- netAnalysis_computeCentrality(cellchat.seu.Asian)
saveRDS(cellchat.seu.Asian, file = "cellchat.Asian-sct.rds")



####Caucasian
###Caucasian
data.input.Caucasian <- GetAssayData(Caucasian, assay = "SCT", slot = "data")
labels <- Idents(Caucasian)
identity.Caucasian <- data.frame(group = labels, row.names = names(labels)) 
cellchat.seu.Caucasian <- createCellChat(object = data.input.Caucasian)
cellchat.seu.Caucasian <- addMeta(cellchat.seu.Caucasian, meta = identity.Caucasian, 
                                  meta.name = "labels")
cellchat.seu.Caucasian <- setIdent(cellchat.seu.Caucasian, ident.use = "labels") 

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat.seu.Caucasian@DB <- CellChatDB.use

cellchat.seu.Caucasian <- subsetData(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- identifyOverExpressedGenes(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- identifyOverExpressedInteractions(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- projectData(cellchat.seu.Caucasian, PPI.human)
cellchat.seu.Caucasian <- computeCommunProb(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- filterCommunication(cellchat.seu.Caucasian, min.cells = 10)
df.net.Caucasian <- subsetCommunication(cellchat.seu.Caucasian)
write.csv(df.net.Caucasian, "net.LR.Caucasian-sct.csv")

cellchat.seu.Caucasian <- computeCommunProbPathway(cellchat.seu.Caucasian)
df.netp.Caucasian <- subsetCommunication(cellchat.seu.Caucasian, slot.name = "netP")
write.csv(df.netp.Caucasian, "net.patway.Caucasian-sct.csv")

cellchat.seu.Caucasian <- aggregateNet(cellchat.seu.Caucasian)
cellchat.seu.Caucasian <- netAnalysis_computeCentrality(cellchat.seu.Caucasian)
saveRDS(cellchat.seu.Caucasian, file = "cellchat.Caucasian-sct.rds")



cellchat.African <- readRDS("~/cellchat-race/cellchat.African-sct.rds")
cellchat.Asian <- readRDS("~/cellchat-race/cellchat.Asian-sct.rds")
cellchat.Caucasian <- readRDS("~/cellchat-race/cellchat.Caucasian-sct.rds")

object.list <- list(Caucasian = cellchat.Caucasian, 
                    Asian = cellchat.Asian,
                    African = cellchat.African
)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


pdf(file = "Caucasian.pdf", height = 10, width = 10)
groupSize_Caucasian <- as.numeric(table(cellchat.Caucasian@idents))
netVisual_circle(cellchat.Caucasian@net$weight, vertex.weight = groupSize_Caucasian, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - Caucasian", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()



pdf(file = "Asian.pdf", height = 10, width = 10)
groupSize_Asian <- as.numeric(table(cellchat.Asian@idents))
netVisual_circle(cellchat.Asian@net$weight, vertex.weight = groupSize_Asian, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - Asian", 
                 arrow.width = 0.7, arrow.size = 0.7)
dev.off()


pdf(file = "African.pdf", height = 10, width = 10)
groupSize_African <- as.numeric(table(cellchat.African@idents))
netVisual_circle(cellchat.African@net$weight, vertex.weight = groupSize_African, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength - African", 
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

pdf(file = "整体互作.pdf", width = 18, height = 5)
patchwork::wrap_plots(plots = gg, nrow = 1)
dev.off()

cols_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
groupSize <- as.numeric(table(cellchat.Caucasian@idents))
mat <- cellchat.Caucasian@net$count
pdf(file = "Caucasian-ec.pdf", height = 15, width = 15)
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
groupSize <- as.numeric(table(cellchat.African@idents))
mat <- cellchat.African@net$count
pdf(file = "African-ec.pdf", height = 15, width = 15)
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
groupSize <- as.numeric(table(cellchat.Asian@idents))
mat <- cellchat.Asian@net$count
pdf(file = "Asian-ec.pdf", height = 15, width = 15)
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


pdf(file = "内皮调控其他细胞.pdf", height = 10, width = 13)
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1, 2, 3, 4, 5, 6, 
                                                               7, 8, 9, 10, 11),  
                 comparison = c(1, 2, 3), angle.x = 45)
dev.off()



pathways.show <- "CXCL"
pdf(file = "African-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.African, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "Asian-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.Asian, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "Caucasian-CXCL.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.Caucasian, signaling = pathways.show, layout = "chord")
dev.off()



pathways.show <- "VEGF"
pdf(file = "African-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.African, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "Asian-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.Asian, signaling = pathways.show, layout = "chord")
dev.off()

pdf(file = "Caucasian-VEGF.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat.Caucasian, signaling = pathways.show, layout = "chord")
dev.off()


library(NMF)
library(ggalluvial)
pdf(file = "African-out.pdf", width = 10, height = 5)
selectK(cellchat.African, pattern = "outgoing")
dev.off()


pdf(file = "African-in.pdf", width = 10, height = 5)
selectK(cellchat.African, pattern = "incoming")
dev.off()


pdf(file = "Asian-out.pdf", width = 10, height = 5)
selectK(cellchat.Asian, pattern = "outgoing")
dev.off()


pdf(file = "Asian-in.pdf", width = 10, height = 5)
selectK(cellchat.Asian, pattern = "incoming")
dev.off()


pdf(file = "Caucasian-out.pdf", width = 10, height = 5)
selectK(cellchat.Caucasian, pattern = "outgoing")
dev.off()


pdf(file = "Caucasian-in.pdf", width = 10, height = 5)
selectK(cellchat.Caucasian, pattern = "incoming")
dev.off()




pdf(file = "Caucasian-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.Caucasian, pattern = "outgoing", k = 2, 
                              height = 10)
dev.off()


pdf(file = "Caucasian-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.Caucasian, pattern = "incoming", k = 2, 
                              height = 10)
dev.off()

pdf(file = "Asian-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.Asian, pattern = "outgoing", k = 3, 
                              height = 12)
dev.off()


pdf(file = "Asian-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.Asian, pattern = "incoming", k = 2, 
                              height = 12)
dev.off()


pdf(file = "African-pattern-out.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.African, pattern = "outgoing", k = 4, 
                              height = 12)
dev.off()


pdf(file = "African-pattern-in.pdf", width = 10, height = 8)
identifyCommunicationPatterns(object = cellchat.African, pattern = "incoming", k = 3, 
                              height = 12)
dev.off()


# river plot
cellchat.African <- identifyCommunicationPatterns(object = cellchat.African, pattern = "outgoing", k = 4, height = 10)
pdf(file = "African.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.African, pattern = "outgoing")
dev.off()

cellchat.African <- identifyCommunicationPatterns(object = cellchat.African, pattern = "incoming", k = 4, height = 10)
pdf(file = "African.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.African, pattern = "incoming")
dev.off()



cellchat.Asian <- identifyCommunicationPatterns(object = cellchat.Asian, pattern = "outgoing", k = 4, height = 10)
pdf(file = "Asian.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.Asian, pattern = "outgoing")
dev.off()

cellchat.Asian <- identifyCommunicationPatterns(object = cellchat.Asian, pattern = "incoming", k = 4, height = 10)
pdf(file = "Asian.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.Asian, pattern = "incoming")
dev.off()



cellchat.Caucasian <- identifyCommunicationPatterns(object = cellchat.Caucasian, pattern = "outgoing", k = 4, height = 10)
pdf(file = "Caucasian.out.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.Caucasian, pattern = "outgoing")
dev.off()

cellchat.Caucasian <- identifyCommunicationPatterns(object = cellchat.Caucasian, pattern = "incoming", k = 4, height = 10)
pdf(file = "Caucasian.in.river.pdf", width = 10, height = 9)
netAnalysis_river(cellchat.Caucasian, pattern = "incoming")
dev.off()










