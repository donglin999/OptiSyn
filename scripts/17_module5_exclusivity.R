#!/usr/bin/env Rscript
# 模块 5: NTS 神经元亚群 4 维度排他性对比
#
# 神经元亚群定义 (基于细胞级别表达):
#   1) Tac1+Vglut2+    : Tac1>0 & Slc17a6>0 & 神经元 score>=0.5
#   2) Vglut2+ only    : Slc17a6>0 & Tac1==0 & 神经元 score>=0.5
#   3) Vgat+_GABA      : Slc32a1>0 & Slc17a6==0 & 神经元 score>=0.5
#   4) Other_neuron    : 神经元 score>=0.5 & 非以上 3 类
#
# 4 维度:
#   D1) 流感信号富集: AddModuleScore(GSE161878 上调 DEG 144)
#   D2) 核心靶点表达: AddModuleScore(58 Mo 严选)
#   D3) 免疫细胞互作强度: 该亚群 vs Cx3cr1+ 免疫的 LR sum_prob
#   D4) 咳嗽相关通路富集: AddModuleScore(感觉 + 痛觉 + 谷氨酸传递 marker)
#
# 排他性结论: 只有 Tac1+ 在 4 维度都显著高于其它亚群
#
# 输出: results/GSE268741/module5_exclusivity/

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
OUT  <- file.path(ROOT, "results/GSE268741/module5_exclusivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. 加载 ------------------------------------------------------------------
cat(">>> 加载\n")
seu <- readRDS(file.path(ROOT, "processed_data/GSE268741/seurat_module1.rds"))
seu_cx3 <- readRDS(file.path(ROOT, "processed_data/GSE268741/sub_Cx3cr1.rds"))
strict <- read.delim(file.path(ROOT,
            "results/Mo_protective/strict/core_targets_strict.tsv"),
            stringsAsFactors = FALSE)
genes58 <- strict$gene_symbol

deg_up_161 <- read.delim(file.path(ROOT, "results/GSE161878/GSE161878_DEG_up.tsv"),
                          stringsAsFactors = FALSE)
flu_up_genes <- unique(deg_up_161$SYMBOL[!is.na(deg_up_161$SYMBOL)])
cat(sprintf("    GSE161878 上调 DEG (主分析剔 D-suffix): %d\n", length(flu_up_genes)))

# 2. 神经元亚群定义 (细胞级别) ---------------------------------------------
cat(">>> 神经元亚群定义\n")
expr_data <- GetAssayData(seu, assay = "RNA", layer = "data")
neuron_panel <- intersect(c("Snap25","Stmn2","Syp","Map2","Syn1","Tubb3",
                            "Rbfox3","Eno2","Uchl1"), rownames(expr_data))
neuron_score <- Matrix::colMeans(expr_data[neuron_panel, ] > 0)

tac1_pos    <- if ("Tac1" %in% rownames(expr_data))
                 (expr_data["Tac1", ] > 0) else rep(FALSE, ncol(seu))
slc17a6_pos <- if ("Slc17a6" %in% rownames(expr_data))
                 (expr_data["Slc17a6", ] > 0) else rep(FALSE, ncol(seu))
slc32a1_pos <- if ("Slc32a1" %in% rownames(expr_data))
                 (expr_data["Slc32a1", ] > 0) else rep(FALSE, ncol(seu))

is_neuron   <- neuron_score >= 0.5
seu$neuron_type <- "Non_neuron"
seu$neuron_type[is_neuron & tac1_pos & slc17a6_pos]                       <- "Tac1+Vglut2+"
seu$neuron_type[is_neuron & slc17a6_pos & !tac1_pos & !slc32a1_pos]       <- "Vglut2+_only"
seu$neuron_type[is_neuron & slc32a1_pos & !slc17a6_pos]                   <- "Vgat+_GABA"
seu$neuron_type[is_neuron & seu$neuron_type == "Non_neuron"]              <- "Other_neuron"

# 4 个亚群的细胞数
neuron_groups <- c("Tac1+Vglut2+","Vglut2+_only","Vgat+_GABA","Other_neuron")
cat("    亚群计数:\n")
print(table(seu$neuron_type))

# 3. AddModuleScore (D1, D2, D4) --------------------------------------------
cat("\n>>> 计算 module score\n")

# D1 流感信号富集 (GSE161878 上调)
flu_in_seu <- intersect(flu_up_genes, rownames(seu))
cat(sprintf("    D1 流感 signature: %d / %d 在 NTS 数据中\n",
            length(flu_in_seu), length(flu_up_genes)))
seu <- AddModuleScore(seu, features = list(flu_in_seu),
                      name = "FluSig_", seed = 42, search = FALSE)
# 列名 FluSig_1
colnames(seu@meta.data)[grepl("^FluSig_", colnames(seu@meta.data))] <- "FluSig"

# D2 核心 58 靶点
g58_in <- intersect(genes58, rownames(seu))
cat(sprintf("    D2 58 靶点: %d / %d 在 NTS 数据中\n",
            length(g58_in), length(genes58)))
seu <- AddModuleScore(seu, features = list(g58_in),
                      name = "Targets58_", seed = 42, search = FALSE)
colnames(seu@meta.data)[grepl("^Targets58_", colnames(seu@meta.data))] <- "Targets58"

# D4 咳嗽 / 感觉神经 / 谷氨酸传递 marker
cough_panel <- c(
  # 痛觉/感觉神经标志
  "Tac1","Calca","Calcb","Trpv1","Trpa1","Scn10a","Scn11a","Piezo2",
  "P2rx3","Mrgprd","Gpr151","Ntrk1",
  # 谷氨酸能 (兴奋性传递, 咳嗽必需)
  "Slc17a6","Slc17a7","Grin1","Grin2a","Grin2b","Gria1","Gria2",
  # 突触传递
  "Syt1","Snap25","Vamp2","Stx1a"
)
cough_in <- intersect(cough_panel, rownames(seu))
cat(sprintf("    D4 咳嗽/感觉神经 panel: %d / %d 在 NTS 数据中\n",
            length(cough_in), length(cough_panel)))
seu <- AddModuleScore(seu, features = list(cough_in),
                      name = "CoughSig_", seed = 42, search = FALSE)
colnames(seu@meta.data)[grepl("^CoughSig_", colnames(seu@meta.data))] <- "CoughSig"

# 4. D3 免疫互作强度: 各神经元亚群 vs Cx3cr1+ 免疫 LR sum_prob ---------------
cat("\n>>> D3 计算各亚群 -> Cx3cr1+ 免疫 LR sum_prob\n")
load(file.path(ROOT, "processed_data/CellChatDB/CellChatDB.mouse.rda"))
ia    <- CellChatDB.mouse$interaction
cmplx <- CellChatDB.mouse$complex
ia_use <- ia[ia$annotation %in% c("Secreted Signaling","Cell-Cell Contact"), ]
expand <- function(x) {
  if (x %in% rownames(cmplx)) {
    s <- as.character(cmplx[x, ]); s[nchar(s) > 0]
  } else x
}
hill <- function(E, K = 0.5, n = 2) E^n / (K^n + E^n)
geo_mean <- function(v) {
  v <- v[!is.na(v) & v > 0]
  if (length(v) == 0) return(NA_real_)
  exp(mean(log(v)))
}

# 受体端 (Cx3cr1+ 免疫) 表达 + rank
cx3_data <- GetAssayData(seu_cx3, assay = "RNA", layer = "data")
cx3_mean <- Matrix::rowMeans(cx3_data)
cx3_pct  <- Matrix::rowMeans(cx3_data > 0)
cx3_rn   <- rank(cx3_mean, ties.method = "average") / length(cx3_mean)
names(cx3_rn) <- names(cx3_mean)

compute_LR_sum <- function(send_cells) {
  s_data <- expr_data[, send_cells, drop = FALSE]
  s_mean <- Matrix::rowMeans(s_data)
  s_pct  <- Matrix::rowMeans(s_data > 0)
  s_rn   <- rank(s_mean, ties.method = "average") / length(s_mean)
  names(s_rn) <- names(s_mean)
  total <- 0; n_LR <- 0; max_p <- 0
  cx3cl1_p <- NA
  for (i in seq_len(nrow(ia_use))) {
    lig <- ia_use$ligand[i]; rec <- ia_use$receptor[i]
    lig_subs <- expand(lig); rec_subs <- expand(rec)
    if (!all(lig_subs %in% names(s_mean))) next
    if (min(s_pct[lig_subs]) < 0.05) next
    if (!all(rec_subs %in% names(cx3_mean))) next
    if (min(cx3_pct[rec_subs]) < 0.05) next
    L_rn <- geo_mean(s_rn[lig_subs])
    R_rn <- geo_mean(cx3_rn[rec_subs])
    P <- hill(L_rn) * hill(R_rn)
    total <- total + P; n_LR <- n_LR + 1
    if (P > max_p) max_p <- P
    if (lig == "Cx3cl1" && rec == "Cx3cr1") cx3cl1_p <- P
  }
  list(sum_prob = total, n_LR = n_LR, max_prob = max_p,
       cx3cl1_cx3cr1_prob = cx3cl1_p)
}

D3_list <- list()
for (gp in neuron_groups) {
  cells <- colnames(seu)[seu$neuron_type == gp]
  if (length(cells) < 5) {
    D3_list[[gp]] <- list(sum_prob = NA, n_LR = 0, max_prob = NA,
                          cx3cl1_cx3cr1_prob = NA)
    cat(sprintf("    %-15s: %d cells (<5), 跳过\n", gp, length(cells)))
    next
  }
  res <- compute_LR_sum(cells)
  D3_list[[gp]] <- res
  cat(sprintf("    %-15s: %d cells, n_LR=%d, sum_prob=%.2f, Cx3cl1->Cx3cr1=%.3f\n",
              gp, length(cells), res$n_LR, res$sum_prob,
              ifelse(is.na(res$cx3cl1_cx3cr1_prob), 0, res$cx3cl1_cx3cr1_prob)))
}

# 5. 汇总表 ----------------------------------------------------------------
m <- seu@meta.data
agg_df <- data.frame()
for (gp in neuron_groups) {
  idx <- which(m$neuron_type == gp)
  if (length(idx) < 5) next
  agg_df <- rbind(agg_df, data.frame(
    neuron_type = gp,
    n_cells     = length(idx),
    FluSig_mean    = mean(m$FluSig[idx]),
    Targets58_mean = mean(m$Targets58[idx]),
    CoughSig_mean  = mean(m$CoughSig[idx]),
    LR_sum_prob       = D3_list[[gp]]$sum_prob,
    LR_n_pairs        = D3_list[[gp]]$n_LR,
    Cx3cl1_Cx3cr1_prob = ifelse(is.na(D3_list[[gp]]$cx3cl1_cx3cr1_prob), 0,
                                D3_list[[gp]]$cx3cl1_cx3cr1_prob),
    stringsAsFactors = FALSE
  ))
}
agg_df <- agg_df[order(-agg_df$Targets58_mean), ]

# 排他性: 4 维度的 rank (1=最高)
rank_desc <- function(x) rank(-x, ties.method = "min")
agg_df$rank_FluSig    <- rank_desc(agg_df$FluSig_mean)
agg_df$rank_Targets58 <- rank_desc(agg_df$Targets58_mean)
agg_df$rank_CoughSig  <- rank_desc(agg_df$CoughSig_mean)
agg_df$rank_LR_sum    <- rank_desc(agg_df$LR_sum_prob)

agg_df$top_in_4D <- (agg_df$rank_FluSig    == 1) +
                    (agg_df$rank_Targets58 == 1) +
                    (agg_df$rank_CoughSig  == 1) +
                    (agg_df$rank_LR_sum    == 1)
write.table(agg_df, file.path(OUT, "exclusivity_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(agg_df, file.path(OUT, "exclusivity_summary.csv"), row.names = FALSE)

cat("\n>>> 4 维度排他性表:\n")
print(agg_df[, c("neuron_type","n_cells","FluSig_mean","Targets58_mean",
                  "CoughSig_mean","LR_sum_prob","Cx3cl1_Cx3cr1_prob",
                  "rank_FluSig","rank_Targets58","rank_CoughSig","rank_LR_sum",
                  "top_in_4D")], row.names = FALSE)

# 6. 统计检验 (Tac1+ vs 其它) -----------------------------------------------
cat("\n>>> 统计检验 (Tac1+Vglut2+ vs 其它每个亚群)\n")
test_res <- data.frame()
ref_idx <- which(m$neuron_type == "Tac1+Vglut2+")
for (gp in setdiff(neuron_groups, "Tac1+Vglut2+")) {
  oth_idx <- which(m$neuron_type == gp)
  if (length(oth_idx) < 5) next
  for (sc in c("FluSig","Targets58","CoughSig")) {
    w <- wilcox.test(m[[sc]][ref_idx], m[[sc]][oth_idx],
                     alternative = "greater")
    test_res <- rbind(test_res, data.frame(
      score = sc,
      vs    = gp,
      Tac1_mean = mean(m[[sc]][ref_idx]),
      other_mean= mean(m[[sc]][oth_idx]),
      p_value = w$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
test_res$p_adj_BH <- p.adjust(test_res$p_value, "BH")
write.table(test_res, file.path(OUT, "wilcox_Tac1_vs_others.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(test_res, row.names = FALSE)

# 7. 可视化 ----------------------------------------------------------------
cat("\n>>> 可视化\n")
sub_meta <- m[m$neuron_type %in% neuron_groups, ]
sub_meta$neuron_type <- factor(sub_meta$neuron_type, levels = neuron_groups)

# 7a. 4 个 Violin (D1/D2/D3 prox/D4)
make_vp <- function(score, lbl, col_pal) {
  ggplot(sub_meta, aes(neuron_type, .data[[score]], fill = neuron_type)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.9) +
    scale_fill_manual(values = col_pal) +
    labs(x = NULL, y = "module score", title = lbl) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}
cols <- c("Tac1+Vglut2+" = "#D62728",
          "Vglut2+_only" = "#FF9933",
          "Vgat+_GABA"   = "#1F77B4",
          "Other_neuron" = "#7F7F7F")
p_d1 <- make_vp("FluSig",    "D1: GSE161878 流感信号富集",         cols)
p_d2 <- make_vp("Targets58", "D2: 58 Mo 严选靶点表达",              cols)
p_d4 <- make_vp("CoughSig",  "D4: 咳嗽/感觉神经/谷氨酸 panel",       cols)

# D3 用条形图 (LR sum_prob)
d3_df <- agg_df[, c("neuron_type","LR_sum_prob","Cx3cl1_Cx3cr1_prob")]
d3_df$neuron_type <- factor(d3_df$neuron_type, levels = neuron_groups)
p_d3 <- ggplot(d3_df, aes(neuron_type, LR_sum_prob, fill = neuron_type)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f", LR_sum_prob)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "总 LR sum_prob",
       title = "D3: 与 Cx3cr1+ 免疫互作强度") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p_combo <- (p_d1 | p_d2) / (p_d3 | p_d4) +
  plot_annotation(title = "模块 5: NTS 神经元亚群 4 维度排他性对比")
ggsave(file.path(OUT, "exclusivity_4D_violin.pdf"), p_combo,
       width = 12, height = 9, device = cairo_pdf)
ggsave(file.path(OUT, "exclusivity_4D_violin.png"), p_combo,
       width = 12, height = 9, dpi = 200, type = "cairo")

# 7b. 雷达 (用 4 D 的 mean 标准化到 [0,1] 直观看 Tac1+ 是否四角都拔尖)
norm01 <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-12)
agg_n <- agg_df
for (sc in c("FluSig_mean","Targets58_mean","CoughSig_mean","LR_sum_prob")) {
  agg_n[[sc]] <- norm01(agg_n[[sc]])
}
# 简化: 用 stacked / faceted bar 替代 radar (ggradar 包未必有)
radar_long <- data.frame(
  neuron_type = rep(agg_n$neuron_type, 4),
  dimension = rep(c("D1_FluSig","D2_Targets58","D3_LR_sum","D4_CoughSig"),
                  each = nrow(agg_n)),
  score_norm = c(agg_n$FluSig_mean, agg_n$Targets58_mean,
                 agg_n$LR_sum_prob, agg_n$CoughSig_mean)
)
radar_long$neuron_type <- factor(radar_long$neuron_type, levels = neuron_groups)
p_heat <- ggplot(radar_long, aes(dimension, neuron_type, fill = score_norm)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", score_norm)), size = 3) +
  scale_fill_gradient(low = "white", high = "#D62728") +
  labs(x = NULL, y = NULL, fill = "norm score",
       title = "4 维度归一化对比 (1 = 最高)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(OUT, "exclusivity_4D_heatmap.pdf"), p_heat,
       width = 7, height = 4, device = cairo_pdf)
ggsave(file.path(OUT, "exclusivity_4D_heatmap.png"), p_heat,
       width = 7, height = 4, dpi = 200, type = "cairo")

# 7c. UMAP 上色 by neuron_type
seu_neurons_only <- subset(seu, cells = colnames(seu)[seu$neuron_type %in% neuron_groups])
p_umap <- DimPlot(seu_neurons_only, group.by = "neuron_type",
                  cols = cols, label = TRUE, label.size = 4) +
  labs(title = "NTS 神经元亚群 (4 类)") +
  theme(legend.position = "right")
ggsave(file.path(OUT, "UMAP_neuron_subtypes.pdf"), p_umap,
       width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(OUT, "UMAP_neuron_subtypes.png"), p_umap,
       width = 8, height = 6, dpi = 200, type = "cairo")

# 7d. 核心靶点 DotPlot (Cx3cr1, Cx3cl1, Gpr151, Gstm1, Klf2, Arrb1) on neurons
core_show <- intersect(c("Cx3cr1","Cx3cl1","Gpr151","Gstm1","Klf2","Arrb1",
                          "Cd209a","Btla","Il16","Tac1","Slc17a6","Slc32a1",
                          "Trpv1","Calca","Scn10a"),
                        rownames(seu_neurons_only))
Idents(seu_neurons_only) <- factor(seu_neurons_only$neuron_type,
                                    levels = neuron_groups)
p_dot <- DotPlot(seu_neurons_only, features = core_show,
                 cols = c("lightgrey","#D62728"), dot.scale = 5) +
  RotatedAxis() +
  labs(title = "核心靶点在神经元亚群的表达") +
  theme(axis.text.x = element_text(size = 8))
ggsave(file.path(OUT, "core_targets_DotPlot.pdf"), p_dot,
       width = max(7, length(core_show) * 0.4 + 2), height = 4,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "core_targets_DotPlot.png"), p_dot,
       width = max(7, length(core_show) * 0.4 + 2), height = 4,
       dpi = 200, type = "cairo", limitsize = FALSE)

# 8. 排他性结论 ------------------------------------------------------------
sink(file.path(OUT, "module5_summary.txt"))
cat("=== 模块 5: NTS 神经元亚群 4 维度排他性验证 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/17_module5_exclusivity.R\n\n")

cat("=== 神经元亚群定义 (细胞级别) ===\n")
print(table(seu$neuron_type))

cat("\n=== 4 维度对比 ===\n")
print(agg_df[, c("neuron_type","n_cells","FluSig_mean","Targets58_mean",
                  "CoughSig_mean","LR_sum_prob","Cx3cl1_Cx3cr1_prob")],
      row.names = FALSE)

cat("\n=== 4 维度排名 (1 = 最高) ===\n")
print(agg_df[, c("neuron_type","rank_FluSig","rank_Targets58",
                  "rank_CoughSig","rank_LR_sum","top_in_4D")],
      row.names = FALSE)

cat("\n=== Wilcoxon Tac1+Vglut2+ vs 其它 (alternative='greater') ===\n")
print(test_res, row.names = FALSE)

cat("\n=== 排他性结论 ===\n")
n_tac1_top <- agg_df$top_in_4D[agg_df$neuron_type == "Tac1+Vglut2+"]
oth_top <- agg_df$top_in_4D[agg_df$neuron_type != "Tac1+Vglut2+"]
cat(sprintf("Tac1+Vglut2+ 在 4 维度中排名第一的次数: %d / 4\n", n_tac1_top))
cat(sprintf("其它亚群最高排名第一次数: %d\n",
            ifelse(length(oth_top) > 0, max(oth_top, na.rm = TRUE), 0)))
if (n_tac1_top == 4) {
  cat(">>> 排他性 [完全成立]: Tac1+Vglut2+ 同时拔尖于 4 维度\n")
} else if (n_tac1_top >= 3) {
  cat(sprintf(">>> 排他性 [3/4 维度成立], 仅 1 维度被其它亚群超过\n"))
} else {
  cat(sprintf(">>> 排他性 [部分成立, %d/4]\n", n_tac1_top))
}
sink()
cat("\n", paste(readLines(file.path(OUT, "module5_summary.txt")), collapse="\n"), sep="")

saveRDS(list(seu_meta = seu@meta.data, agg_df = agg_df,
             test_res = test_res, D3_list = D3_list,
             flu_in_seu = flu_in_seu, g58_in = g58_in,
             cough_in = cough_in),
        file.path(OUT, "module5_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
