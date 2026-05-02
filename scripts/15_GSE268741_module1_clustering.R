#!/usr/bin/env Rscript
# GSE268741 (Gannot 2024 NTS scRNA-seq) 模块 1: 分群 + 子集分选
#
# 单样本 (GSM8298251), 9638 cells x 33696 genes
# 流程: QC -> Norm -> HVG -> Scale -> PCA -> UMAP -> Cluster -> Marker 注释
# 子集: Tac1+(Slc17a6+) 神经元 / Cx3cr1+ 免疫 / 其它神经元
#
# 输出: results/GSE268741/module1_clustering/

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/GSE268741/module1_clustering")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  min.cells = 3, min.features = 200,
  pctmt_max = 15, nFeature_min = 500, nFeature_max = 7500,
  hvg_n = 2000, npc = 30,
  resolution = 0.5,
  marker_thr = 0.5,    # 用于 "expressing" 判断 (log-norm 表达 > 0)
  tac1_min_pct = 0.3   # Tac1+ 神经元簇要求簇内 ≥30% 细胞表达 Tac1
)

# 1. 加载 + 创建 Seurat 对象 ------------------------------------------------
cat(">>> 加载 sparse matrix\n")
M <- readRDS(file.path(ROOT, "processed_data/GSE268741/counts_sparse.rds"))
cat(sprintf("    %d gene x %d cell, nnz=%d\n", nrow(M), ncol(M), length(M@x)))

seu <- CreateSeuratObject(counts = M,
                          project = "GSE268741",
                          min.cells = PARAMS$min.cells,
                          min.features = PARAMS$min.features)
seu[["pct_mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
cat(sprintf(">>> 初始 Seurat 对象: %d cells\n", ncol(seu)))
cat(sprintf("    nFeature mean=%.0f median=%.0f range=%d-%d\n",
            mean(seu$nFeature_RNA), median(seu$nFeature_RNA),
            min(seu$nFeature_RNA),  max(seu$nFeature_RNA)))
cat(sprintf("    pct_mt   mean=%.2f median=%.2f range=%.2f-%.2f\n",
            mean(seu$pct_mt), median(seu$pct_mt),
            min(seu$pct_mt),  max(seu$pct_mt)))

# QC plot (pre-filter)
p_qc <- VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","pct_mt"),
                ncol = 3, pt.size = 0.05) &
  theme(legend.position = "none")
ggsave(file.path(OUT, "QC_pre.pdf"), p_qc, width = 10, height = 4,
       device = cairo_pdf)
ggsave(file.path(OUT, "QC_pre.png"), p_qc, width = 10, height = 4,
       dpi = 200, type = "cairo")

# 2. QC 过滤 ----------------------------------------------------------------
cat(">>> QC 过滤\n")
seu <- subset(seu,
              subset = nFeature_RNA >= PARAMS$nFeature_min &
                       nFeature_RNA <= PARAMS$nFeature_max &
                       pct_mt < PARAMS$pctmt_max)
cat(sprintf("    保留 %d cells (nFeat %d-%d, pct_mt<%g)\n",
            ncol(seu),
            PARAMS$nFeature_min, PARAMS$nFeature_max,
            PARAMS$pctmt_max))

# 3. 标准化 + HVG + Scale + PCA --------------------------------------------
cat(">>> Normalize + HVG + Scale + PCA\n")
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = PARAMS$hvg_n, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = PARAMS$npc, verbose = FALSE)

# 4. UMAP + Cluster ---------------------------------------------------------
cat(">>> UMAP + Cluster\n")
seu <- RunUMAP(seu, dims = 1:PARAMS$npc, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:PARAMS$npc, verbose = FALSE)
seu <- FindClusters(seu, resolution = PARAMS$resolution, verbose = FALSE)
cat(sprintf("    cluster 数: %d\n", length(levels(seu))))

# 保存中间 Seurat 对象 (供后续模块直接读)
saveRDS(seu, file.path(ROOT, "processed_data/GSE268741/seurat_module1.rds"))

# 5. Marker 注释 ------------------------------------------------------------
cat(">>> Marker 注释\n")

# 8 主群预期 marker (依据文献 + 经典经验)
marker_list <- list(
  Neuron      = c("Snap25","Stmn2","Syp","Map2","Syn1","Tubb3","Rbfox3","Eno2"),
  Vglut2_Glu  = c("Slc17a6"),                    # 兴奋性谷氨酸能 (Vglut2)
  Vgat_GABA   = c("Slc32a1","Gad1","Gad2"),     # 抑制性 GABA能
  Tac1        = c("Tac1"),                       # 咳嗽核心标记
  Cx3cr1_Imm  = c("Cx3cr1","Aif1","C1qa","C1qb","C1qc","Csf1r","P2ry12"),
  Astrocyte   = c("Gfap","Slc1a3","Aqp4","Sox9"),
  Oligo       = c("Olig1","Olig2","Mog","Mbp","Plp1"),
  Endothel    = c("Pecam1","Cdh5","Cldn5"),
  Ependyma    = c("Foxj1","Ccdc153"),
  Pericyte    = c("Pdgfrb","Mcam","Rgs5"),
  Schwann     = c("Mpz","Pmp22")
)
all_markers <- unique(unlist(marker_list))
all_markers <- intersect(all_markers, rownames(seu))
cat(sprintf("    Marker 在数据中: %d / %d\n",
            length(all_markers), length(unique(unlist(marker_list)))))

# 簇 x marker 平均表达 + pct
mks_avg <- AverageExpression(seu, features = all_markers,
                             assays = "RNA", layer = "data",
                             group.by = "seurat_clusters")$RNA
expr_pct <- function(seu, gene, cluster_id) {
  cells <- colnames(seu)[seu$seurat_clusters == cluster_id]
  v <- GetAssayData(seu, layer = "data")[gene, cells]
  mean(v > 0)
}
cluster_ids <- levels(seu)
mks_pct <- matrix(NA, nrow = length(all_markers), ncol = length(cluster_ids),
                   dimnames = list(all_markers, cluster_ids))
for (g in all_markers) for (c in cluster_ids) {
  mks_pct[g, c] <- expr_pct(seu, g, c)
}

# 主簇身份打分: 哪类 marker 在该簇内最多基因 mean log-expr 高
cluster_anno <- data.frame(
  cluster = cluster_ids,
  n_cells = as.integer(table(seu$seurat_clusters)[cluster_ids]),
  stringsAsFactors = FALSE
)
for (mn in names(marker_list)) {
  ms <- intersect(marker_list[[mn]], rownames(seu))
  if (length(ms) == 0) {
    cluster_anno[[paste0("avg_", mn)]] <- NA; next
  }
  if (length(ms) == 1) {
    cluster_anno[[paste0("avg_", mn)]] <- as.numeric(mks_avg[ms, , drop = TRUE])
  } else {
    cluster_anno[[paste0("avg_", mn)]] <- as.numeric(colMeans(mks_avg[ms, , drop = FALSE]))
  }
}
# 单 marker 用于子集判定
cluster_anno$pct_Tac1   <- as.numeric(mks_pct["Tac1",   ])
cluster_anno$pct_Cx3cr1 <- as.numeric(mks_pct["Cx3cr1", ])
cluster_anno$pct_Slc17a6<- if ("Slc17a6" %in% rownames(mks_pct)) as.numeric(mks_pct["Slc17a6",]) else NA
cluster_anno$pct_Slc32a1<- if ("Slc32a1" %in% rownames(mks_pct)) as.numeric(mks_pct["Slc32a1",]) else NA

# 自动指派 main_class
class_cols <- grep("^avg_", colnames(cluster_anno), value = TRUE)
class_mat  <- as.matrix(cluster_anno[, class_cols])
cluster_anno$top_class <- sub("^avg_", "", class_cols[apply(class_mat, 1, which.max)])

write.table(cluster_anno, file.path(OUT, "cluster_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n    cluster 注释 (top_class):\n"))
print(cluster_anno[, c("cluster","n_cells","top_class",
                       "pct_Tac1","pct_Cx3cr1",
                       "pct_Slc17a6","pct_Slc32a1")],
      row.names = FALSE)

# 6. 子集定义 ---------------------------------------------------------------
cat("\n>>> 子集定义\n")
expr_data <- GetAssayData(seu, layer = "data")
tac1_pos    <- expr_data["Tac1", ]    > 0
slc17a6_pos <- if ("Slc17a6" %in% rownames(expr_data)) {
  expr_data["Slc17a6", ] > 0
} else {
  rep(FALSE, ncol(seu))
}
cx3cr1_pos  <- expr_data["Cx3cr1", ]  > 0
neuron_score<- colMeans(expr_data[intersect(marker_list$Neuron,
                                            rownames(expr_data)), ] > 0)

# Tac1+ Vglut2+ 神经元: 细胞 Tac1>0 & Slc17a6>0 & 神经元 marker score>=0.5 (≥半数 panel marker 表达)
cell_tac1_glu <- tac1_pos & slc17a6_pos & (neuron_score >= 0.5)
# Cx3cr1+ 免疫: Cx3cr1>0 & 神经元 score 低
cell_cx3cr1   <- cx3cr1_pos & (neuron_score < 0.3)
# 其它神经元: 神经元 score>=0.5 但非 Tac1+
cell_other_neu <- (neuron_score >= 0.5) & !cell_tac1_glu

seu$subset_label <- "Other"
seu$subset_label[cell_tac1_glu]   <- "Tac1+Vglut2+"
seu$subset_label[cell_cx3cr1]     <- "Cx3cr1+_immune"
seu$subset_label[cell_other_neu]  <- "Other_neuron"
cat("    subset_label 计数:\n"); print(table(seu$subset_label))

# 7. 命名 cluster (基于 top_class + 关键 marker) ----------------------------
ident_named <- as.character(seu$seurat_clusters)
for (i in seq_len(nrow(cluster_anno))) {
  cid <- cluster_anno$cluster[i]
  cls <- cluster_anno$top_class[i]
  pct_tac1 <- cluster_anno$pct_Tac1[i]
  pct_cx3 <- cluster_anno$pct_Cx3cr1[i]
  lbl <- cls
  if (cls == "Neuron" && !is.na(pct_tac1) && pct_tac1 >= PARAMS$tac1_min_pct) {
    lbl <- "Neuron_Tac1+"
  } else if (cls == "Cx3cr1_Imm" || (!is.na(pct_cx3) && pct_cx3 > 0.5)) {
    lbl <- "Microglia_Cx3cr1+"
  }
  ident_named[ident_named == cid] <- paste0(lbl, "_c", cid)
}
seu$cluster_label <- ident_named
cat("\n    cluster_label 计数:\n"); print(sort(table(seu$cluster_label), decreasing = TRUE))

saveRDS(seu, file.path(ROOT, "processed_data/GSE268741/seurat_module1.rds"))

# 8. 可视化 -----------------------------------------------------------------
cat("\n>>> 可视化\n")

# 8a. UMAP by cluster
p1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE, label.size = 4) +
  labs(title = sprintf("UMAP by cluster (resolution=%.2f, n_cluster=%d)",
                       PARAMS$resolution, length(levels(seu)))) +
  theme(legend.position = "none")
# 8b. UMAP by top_class
seu$top_class <- cluster_anno$top_class[match(seu$seurat_clusters, cluster_anno$cluster)]
p2 <- DimPlot(seu, group.by = "top_class", label = TRUE, label.size = 3.5,
              repel = TRUE) +
  labs(title = "UMAP by top class")
# 8c. UMAP by subset
p3 <- DimPlot(seu, group.by = "subset_label",
              cols = c("Tac1+Vglut2+" = "#D62728",
                       "Cx3cr1+_immune" = "#1F77B4",
                       "Other_neuron" = "#2CA02C",
                       "Other" = "grey80"),
              order = c("Other","Other_neuron","Cx3cr1+_immune","Tac1+Vglut2+")) +
  labs(title = "UMAP: 子集 (红=Tac1+Vglut2+, 蓝=Cx3cr1+, 绿=其它神经元)")

p_combo <- (p1 | p2) / p3
ggsave(file.path(OUT, "UMAP_overview.pdf"), p_combo,
       width = 14, height = 11, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "UMAP_overview.png"), p_combo,
       width = 14, height = 11, dpi = 200, type = "cairo", limitsize = FALSE)
cat("    UMAP_overview.pdf/png\n")

# 8d. Marker DotPlot
key_markers <- c("Snap25","Map2","Tubb3","Rbfox3",
                 "Slc17a6","Slc32a1","Gad1","Tac1",
                 "Cx3cr1","Aif1","C1qa","Csf1r",
                 "Gfap","Aqp4","Olig1","Mog","Mbp",
                 "Pecam1","Cldn5","Pdgfrb","Mpz")
key_markers <- intersect(key_markers, rownames(seu))
Idents(seu) <- "seurat_clusters"
p_dot <- DotPlot(seu, features = key_markers,
                 cols = c("lightgrey", "#D62728"), dot.scale = 5) +
  RotatedAxis() +
  labs(title = "Marker DotPlot per cluster") +
  theme(axis.text.x = element_text(size = 8))
ggsave(file.path(OUT, "marker_DotPlot.pdf"), p_dot,
       width = max(8, length(key_markers) * 0.35 + 3), height = 5,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "marker_DotPlot.png"), p_dot,
       width = max(8, length(key_markers) * 0.35 + 3), height = 5,
       dpi = 200, type = "cairo", limitsize = FALSE)

# 8e. FeaturePlot for 4 关键 marker
p_fp <- FeaturePlot(seu,
                    features = c("Tac1","Slc17a6","Cx3cr1","Slc32a1"),
                    ncol = 2, order = TRUE,
                    cols = c("lightgrey", "#D62728"))
ggsave(file.path(OUT, "FeaturePlot_keymarkers.pdf"), p_fp,
       width = 10, height = 9, device = cairo_pdf)
ggsave(file.path(OUT, "FeaturePlot_keymarkers.png"), p_fp,
       width = 10, height = 9, dpi = 200, type = "cairo")

# 9. 子集导出 (用于模块 2/4) ------------------------------------------------
seu_tac1   <- subset(seu, cells = colnames(seu)[seu$subset_label == "Tac1+Vglut2+"])
seu_cx3    <- subset(seu, cells = colnames(seu)[seu$subset_label == "Cx3cr1+_immune"])
seu_neu    <- subset(seu, cells = colnames(seu)[seu$subset_label %in%
                                                c("Tac1+Vglut2+","Other_neuron")])
saveRDS(seu_tac1, file.path(ROOT, "processed_data/GSE268741/sub_Tac1_Vglut2.rds"))
saveRDS(seu_cx3,  file.path(ROOT, "processed_data/GSE268741/sub_Cx3cr1.rds"))
saveRDS(seu_neu,  file.path(ROOT, "processed_data/GSE268741/sub_neurons.rds"))
cat(sprintf("\n>>> 子集对象保存:\n"))
cat(sprintf("    Tac1+Vglut2+    : %d cells -> sub_Tac1_Vglut2.rds\n",   ncol(seu_tac1)))
cat(sprintf("    Cx3cr1+_immune  : %d cells -> sub_Cx3cr1.rds\n",        ncol(seu_cx3)))
cat(sprintf("    All neurons     : %d cells -> sub_neurons.rds\n",       ncol(seu_neu)))

# 10. 汇总 -----------------------------------------------------------------
sink(file.path(OUT, "module1_summary.txt"))
cat("=== GSE268741 模块 1: NTS 单细胞分群 + 子集分选 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/15_GSE268741_module1_clustering.R\n")
cat("数据: GSE268741 GSM8298251 (单样本, 9638 cells x 33696 genes)\n\n")

cat("=== 参数 ===\n"); print(PARAMS)

cat(sprintf("\n=== QC ===\n"))
cat(sprintf("初始 cells (CreateSeuratObject): %d\n", length(M@Dim) > 0))
cat(sprintf("QC 后 cells: %d\n", ncol(seu)))

cat(sprintf("\n=== 聚类 ===\n"))
cat(sprintf("PC=%d, resolution=%.2f, n_cluster=%d\n",
            PARAMS$npc, PARAMS$resolution, length(levels(seu))))

cat("\n=== Cluster 注释 ===\n")
print(cluster_anno[, c("cluster","n_cells","top_class",
                       "pct_Tac1","pct_Cx3cr1","pct_Slc17a6","pct_Slc32a1")],
      row.names = FALSE)

cat("\n=== 子集计数 ===\n")
print(table(seu$subset_label))

cat(sprintf("\nTop class 计数:\n"))
print(table(seu$top_class))
sink()
cat("\n", paste(readLines(file.path(OUT, "module1_summary.txt")), collapse="\n"), sep="")

cat(sprintf("\n>>> 输出: %s\n", OUT))
