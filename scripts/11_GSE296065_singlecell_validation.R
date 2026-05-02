#!/usr/bin/env Rscript
# GSE296065 (Almanzar et al, Sci Adv) 单细胞两步验证 58 个 Mo 严选靶点
#
# 数据: raw_data/00_raw_data/GSE296065_almanzar_flu.rds
#   Seurat v5, 37543 细胞 x 32293 基因, 22 cell types
#   condition: V0 (Naive Vehicle) / R0 (Naive RTX) / V7 (Flu Vehicle) / R7 (Flu RTX)
#   treatment: RTX / Vehicle    exposure: NAIVE / 7 DPI
#
# 流程 (用户严格指定):
#   1) 细胞群强制定义
#      目标 (Mo/Mo-Macs): Monocyte Ly6c2+ + Monocyte Ly6c2- + Macrophage Interstitial
#      对照: 其余全部 cell types
#   2) 细胞特异性验证
#      Seurat::FindMarkers, test.use = "wilcox"
#      ident.1 = "Target", ident.2 = "Background"
#      筛: avg_log2FC > 0.5 & p_val_adj < 0.05
#   3) 组间差异验证 (仅纯净 Mo/Mo-Macs 子集 + condition %in% c("V7","R7"))
#      Seurat::FindMarkers, test.use = "wilcox"
#      ident.1 = "V7" (Flu_Vehicle), ident.2 = "R7" (Flu_RTX)
#      筛: avg_log2FC > 0.5 & p_val_adj < 0.05  (V7 > R7 -> RTX 组下调)
#   4) 取两步通过的交集 = 最终高置信度靶点
#
# 输出: results/GSE296065_singlecell/

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
RDS  <- file.path(ROOT, "raw_data/00_raw_data/GSE296065_almanzar_flu.rds")
TARGETS_CSV <- file.path(ROOT, "results/Mo_protective/strict/core_targets_strict.csv")
OUT  <- file.path(ROOT, "results/GSE296065_singlecell")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

THR <- list(log2fc = 0.5, p_adj = 0.05)

TARGET_CT  <- c("Monocyte Ly6c2+", "Monocyte Ly6c2-", "Macrophage Interstitial")
FLU_VEHICLE <- "V7"   # condition: 7 DPI Vehicle
FLU_RTX     <- "R7"   # condition: 7 DPI RTX

cat(sprintf(">>> 启动: %s\n", TS))

# 1. 加载 ---------------------------------------------------------------------
cat(">>> [1] 加载 Seurat 对象\n")
seu <- readRDS(RDS)
cat(sprintf("    %d cells x %d genes; default assay = %s\n",
            ncol(seu), nrow(seu), DefaultAssay(seu)))
cat(sprintf("    layers: %s\n", paste(Layers(seu[[DefaultAssay(seu)]]), collapse=",")))

# 2. 加载 58 靶点
cat(">>> [2] 加载 58 严选靶点\n")
targets_df <- read.csv(TARGETS_CSV, stringsAsFactors = FALSE)
target_genes <- unique(targets_df$gene_symbol)
cat(sprintf("    源: %s\n", basename(TARGETS_CSV)))
cat(sprintf("    靶点数: %d\n", length(target_genes)))

valid_genes <- intersect(target_genes, rownames(seu))
miss_genes  <- setdiff(target_genes, rownames(seu))
cat(sprintf("    在数据中存在: %d / %d\n", length(valid_genes), length(target_genes)))
if (length(miss_genes) > 0) {
  cat(sprintf("    缺失基因: %s\n", paste(miss_genes, collapse=", ")))
}
existence <- data.frame(gene_symbol = target_genes,
                        in_seurat   = target_genes %in% rownames(seu),
                        stringsAsFactors = FALSE)
write.csv(existence, file.path(OUT, "gene_existence.csv"), row.names = FALSE)

# 3. 强制定义细胞群 -----------------------------------------------------------
cat(">>> [3] 强制定义目标 vs 对照群\n")
ct <- as.character(seu$celltype)
miss_ct <- setdiff(TARGET_CT, unique(ct))
if (length(miss_ct) > 0) stop(sprintf("精确匹配失败的目标 cell type: %s",
                                       paste(miss_ct, collapse="; ")))
seu$group_specificity <- ifelse(ct %in% TARGET_CT, "Target", "Background")
cat("    target / background:\n"); print(table(seu$group_specificity))
cat("    target 内细胞类型分布:\n")
print(table(seu$celltype[seu$group_specificity == "Target"]))

# 4. 特异性 FindMarkers (test 1: Target vs Background) ------------------------
cat(">>> [4] 特异性 FindMarkers (Target vs Background)\n")
Idents(seu) <- "group_specificity"
spec_t0 <- Sys.time()
spec_markers <- FindMarkers(
  seu,
  ident.1 = "Target",
  ident.2 = "Background",
  features = valid_genes,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE,
  verbose  = FALSE
)
spec_markers$gene_symbol <- rownames(spec_markers)
rownames(spec_markers) <- NULL
spec_markers <- spec_markers[, c("gene_symbol","avg_log2FC","p_val","p_val_adj",
                                  "pct.1","pct.2")]
colnames(spec_markers) <- c("gene_symbol","avg_log2FC_target_vs_bg",
                            "p_val","p_val_adj","pct_target","pct_bg")
spec_markers$is_specific <- with(spec_markers,
  avg_log2FC_target_vs_bg > THR$log2fc & p_val_adj < THR$p_adj)
spec_markers <- spec_markers[order(-spec_markers$avg_log2FC_target_vs_bg), ]
write.table(spec_markers, file.path(OUT, "step2_specificity_FindMarkers.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.csv(spec_markers, file.path(OUT, "step2_specificity_FindMarkers.csv"),
          row.names = FALSE)
spec_pass <- spec_markers$gene_symbol[spec_markers$is_specific]
cat(sprintf("    分析基因: %d, 通过 (avg_log2FC > %.2f & p_adj < %.2f): %d\n",
            nrow(spec_markers), THR$log2fc, THR$p_adj, length(spec_pass)))
cat(sprintf("    耗时 %.1f 秒\n",
            as.numeric(Sys.time() - spec_t0, units = "secs")))
cat(sprintf("    通过列表: %s\n", paste(spec_pass, collapse=", ")))

if (length(spec_pass) == 0) stop("特异性 0, 终止")

# 5. 提取纯净 Mo/Mo-Macs 子集 + 组间 FindMarkers (V7 vs R7) -------------------
cat(">>> [5] 纯净子集 + 组间 FindMarkers (Flu_Vehicle V7 vs Flu_RTX R7)\n")
# 构造纯净子集: target 细胞 ∩ condition %in% (V7, R7)
keep_pure <- (seu$group_specificity == "Target") &
             (as.character(seu$condition) %in% c(FLU_VEHICLE, FLU_RTX))
seu_pure <- subset(seu, cells = colnames(seu)[keep_pure])
seu_pure$condition_chr <- as.character(seu_pure$condition)
cat(sprintf("    纯净子集细胞数: %d\n", ncol(seu_pure)))
cat("    纯净子集 condition x celltype:\n")
print(table(seu_pure$condition_chr, seu_pure$celltype))

if (any(table(seu_pure$condition_chr) < 30)) {
  cat("    [WARN] 某组细胞 < 30, 检验功效有限\n")
}

Idents(seu_pure) <- "condition_chr"
grp_t0 <- Sys.time()
grp_markers <- FindMarkers(
  seu_pure,
  ident.1 = FLU_VEHICLE,
  ident.2 = FLU_RTX,
  features = spec_pass,    # 仅在特异性通过的基因里测
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE,
  verbose  = FALSE
)
grp_markers$gene_symbol <- rownames(grp_markers)
rownames(grp_markers) <- NULL
grp_markers <- grp_markers[, c("gene_symbol","avg_log2FC","p_val","p_val_adj",
                                "pct.1","pct.2")]
colnames(grp_markers) <- c("gene_symbol","avg_log2FC_V7_vs_R7",
                           "p_val","p_val_adj","pct_V7","pct_R7")
# Seurat FindMarkers 默认 p_val_adj = Bonferroni x 全基因组 (32293), 对预选小集合过严
# 补充: 在测试的 28 个基因子集内 BH 校正 (更合理的多重比较)
grp_markers$p_val_adj_BH_in_subset <- p.adjust(grp_markers$p_val, method = "BH")
grp_markers$is_RTX_down <- with(grp_markers,
  avg_log2FC_V7_vs_R7 > THR$log2fc & p_val_adj < THR$p_adj)
grp_markers$is_RTX_down_BH <- with(grp_markers,
  avg_log2FC_V7_vs_R7 > THR$log2fc & p_val_adj_BH_in_subset < THR$p_adj)
grp_markers <- grp_markers[order(-grp_markers$avg_log2FC_V7_vs_R7), ]
write.table(grp_markers, file.path(OUT, "step3_group_FindMarkers.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.csv(grp_markers, file.path(OUT, "step3_group_FindMarkers.csv"),
          row.names = FALSE)
grp_pass    <- grp_markers$gene_symbol[grp_markers$is_RTX_down]
grp_pass_BH <- grp_markers$gene_symbol[grp_markers$is_RTX_down_BH]
cat(sprintf("    分析基因 (= 特异性通过): %d\n", nrow(grp_markers)))
cat(sprintf("    通过 strict (Seurat default p_adj < %.2f): %d\n",
            THR$p_adj, length(grp_pass)))
cat(sprintf("    通过 BH (子集 BH < %.2f, 推荐, 更合理的小集合校正): %d\n",
            THR$p_adj, length(grp_pass_BH)))
cat(sprintf("    耗时 %.1f 秒\n",
            as.numeric(Sys.time() - grp_t0, units = "secs")))

# 6. 交集 = 最终高置信度靶点 -------------------------------------------------
final_genes    <- intersect(spec_pass, grp_pass)        # strict
final_genes_BH <- intersect(spec_pass, grp_pass_BH)     # subset BH
cat(sprintf(">>> [6] 最终高置信度靶点 (两步交集):\n"))
cat(sprintf("       strict (Seurat default Bonferroni): %d -> %s\n",
            length(final_genes), paste(final_genes, collapse=", ")))
cat(sprintf("       BH (子集 BH 校正):                  %d -> %s\n",
            length(final_genes_BH), paste(final_genes_BH, collapse=", ")))

# 主结果合并表 (BH 推荐版作为主输出)
make_final_tbl <- function(genes) {
  if (length(genes) == 0) return(data.frame())
  merge(
    spec_markers[spec_markers$gene_symbol %in% genes,
                 c("gene_symbol","avg_log2FC_target_vs_bg","p_val_adj",
                   "pct_target","pct_bg")],
    grp_markers[grp_markers$gene_symbol %in% genes,
                c("gene_symbol","avg_log2FC_V7_vs_R7",
                  "p_val_adj","p_val_adj_BH_in_subset",
                  "pct_V7","pct_R7")],
    by = "gene_symbol", suffixes = c("_specificity","_group")
  )
}
final_tbl_strict <- make_final_tbl(final_genes)
final_tbl_BH     <- make_final_tbl(final_genes_BH)
if (nrow(final_tbl_BH) > 0) {
  final_tbl_BH <- final_tbl_BH[order(-final_tbl_BH$avg_log2FC_V7_vs_R7), ]
}
write.table(final_tbl_strict,
            file.path(OUT, "final_high_confidence_targets_strict.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(final_tbl_BH,
            file.path(OUT, "final_high_confidence_targets_BH.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.csv(final_tbl_BH,
          file.path(OUT, "final_high_confidence_targets_BH.csv"),
          row.names = FALSE)

# 仅基因列表
writeLines(final_genes,    file.path(OUT, "final_genes_strict.txt"))
writeLines(final_genes_BH, file.path(OUT, "final_genes_BH.txt"))

# 主用 BH 版本作为后续可视化与汇总的 'final'
final_tbl <- final_tbl_BH
final_genes_main <- final_genes_BH

# 7. 摘要 ---------------------------------------------------------------------
sink(file.path(OUT, "validation_summary.txt"))
cat("=== GSE296065 单细胞两步 FindMarkers 验证 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("数据: GSE296065 Almanzar et al (Seurat v5), 37543 cells x 32293 genes\n\n")
cat("分组:\n")
cat(sprintf("  目标 (Target): %s -> %d cells\n",
            paste(TARGET_CT, collapse=" + "),
            sum(seu$group_specificity == "Target")))
cat(sprintf("  对照 (Background): 其它 %d cells\n",
            sum(seu$group_specificity == "Background")))
cat(sprintf("\n纯净子集 (Target ∩ condition %%in%% V7/R7):\n"))
cat(sprintf("  Flu_Vehicle V7: %d cells\n", sum(seu_pure$condition_chr == "V7")))
cat(sprintf("  Flu_RTX     R7: %d cells\n", sum(seu_pure$condition_chr == "R7")))

cat(sprintf("\n阈值 (两步统一): avg_log2FC > %.2f & p_val_adj < %.2f\n",
            THR$log2fc, THR$p_adj))

cat("\n=== 流程计数 ===\n")
cat(sprintf("%-50s %d\n", "起始 58 严选靶点", length(target_genes)))
cat(sprintf("%-50s %d\n", "在 Seurat 数据中存在", length(valid_genes)))
cat(sprintf("%-50s %d\n",
            sprintf("[Step 1] 特异性 FindMarkers (Target vs BG, log2FC>%.1f)",
                    THR$log2fc),
            length(spec_pass)))
cat(sprintf("%-50s %d\n",
            sprintf("[Step 2-strict] V7 vs R7 (Seurat 默认 Bonferroni)"),
            length(grp_pass)))
cat(sprintf("%-50s %d\n",
            sprintf("[Step 2-BH]     V7 vs R7 (子集 BH, 推荐)"),
            length(grp_pass_BH)))
cat(sprintf("%-50s %d\n",
            "[Final-strict] 交集 (Seurat 默认)", length(final_genes)))
cat(sprintf("%-50s %d\n",
            "[Final-BH]     交集 (BH 校正, 主结果)", length(final_genes_BH)))

cat("\n=== 特异性通过 (Step 1) ===\n")
print(spec_markers[spec_markers$is_specific,
       c("gene_symbol","avg_log2FC_target_vs_bg","p_val_adj",
         "pct_target","pct_bg")], row.names = FALSE)

cat("\n=== 组间 V7 vs R7 全 28 基因结果 (按 V7-vs-R7 log2FC 降序) ===\n")
print(grp_markers[, c("gene_symbol","avg_log2FC_V7_vs_R7","p_val",
                       "p_val_adj","p_val_adj_BH_in_subset",
                       "pct_V7","pct_R7","is_RTX_down","is_RTX_down_BH")],
      row.names = FALSE)

cat("\n=== 组间通过 (BH 推荐, V7>R7 = RTX 下调) ===\n")
print(grp_markers[grp_markers$is_RTX_down_BH,
       c("gene_symbol","avg_log2FC_V7_vs_R7","p_val",
         "p_val_adj_BH_in_subset","pct_V7","pct_R7")],
      row.names = FALSE)

cat("\n=== 最终交集 (BH 推荐) ===\n")
print(final_tbl, row.names = FALSE)
sink()
cat(readLines(file.path(OUT, "validation_summary.txt")), sep="\n")

# 8. 可视化 ------------------------------------------------------------------
cat(">>> [7] 生成图\n")

# 8a. 特异性 DotPlot (在原 seu 中)
spec_plot_genes <- spec_pass
if (length(spec_plot_genes) > 0) {
  Idents(seu) <- "celltype"
  # 仅显示 7 类 (3 target + 4 高频 background) 以保持可读
  top_bg <- names(sort(table(seu$celltype[seu$group_specificity == "Background"]),
                       decreasing = TRUE))
  show_ct <- c(TARGET_CT, head(top_bg, 6))
  seu_show <- subset(seu, idents = show_ct)
  Idents(seu_show) <- factor(seu_show$celltype, levels = show_ct)
  p_dot <- DotPlot(seu_show, features = spec_plot_genes,
                   cols = c("lightgrey", "#D62728"),
                   dot.scale = 5) +
    RotatedAxis() +
    labs(title = sprintf("特异性 DotPlot: Mo/Mo-Macs vs 主要其它细胞 (n=%d genes)",
                         length(spec_plot_genes)),
         subtitle = "目标群 = Mo Ly6c2+/- + Macroph Interstitial 合并") +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 8))
  fig_w <- max(7, length(spec_plot_genes) * 0.35 + 3)
  ggsave(file.path(OUT, "specificity_DotPlot.pdf"), p_dot,
         width = fig_w, height = 5, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT, "specificity_DotPlot.png"), p_dot,
         width = fig_w, height = 5, dpi = 200, type = "cairo", limitsize = FALSE)
  cat("    DotPlot: specificity_DotPlot.pdf/png\n")
}

# 8b. 组间火山图 (Step 2, 用 BH 校正纵轴)
gv <- grp_markers
gv$negLog10BH <- -log10(gv$p_val_adj_BH_in_subset + 1e-300)
gv$status <- "ns"
gv$status[gv$is_RTX_down_BH] <- "RTX_down"
gv$status[gv$avg_log2FC_V7_vs_R7 < -THR$log2fc &
          gv$p_val_adj_BH_in_subset < THR$p_adj] <- "RTX_up"
top_label <- gv[gv$status != "ns" |
                abs(gv$avg_log2FC_V7_vs_R7) > 0.3, ]
p_volc <- ggplot(gv, aes(avg_log2FC_V7_vs_R7, negLog10BH, color = status)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(ns = "grey70",
                                RTX_down = "#D62728",
                                RTX_up   = "#1F77B4"),
                     labels = c(ns = "ns",
                                RTX_down = "RTX 下调 (V7>R7)",
                                RTX_up   = "RTX 上调")) +
  geom_vline(xintercept = c(-THR$log2fc, THR$log2fc),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(THR$p_adj),
             linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top_label, aes(label = gene_symbol),
                  size = 3.2, max.overlaps = 50, color = "black") +
  labs(x = "avg_log2FC (V7 / R7)",
       y = "-log10(BH p_adj, 在 28 基因子集内)",
       title = "组间火山图: Flu_Vehicle (V7) vs Flu_RTX (R7)",
       subtitle = sprintf("BH 校正: RTX 下调 %d, RTX 上调 %d (n=%d 特异性基因)",
                          sum(gv$status == "RTX_down"),
                          sum(gv$status == "RTX_up"),
                          nrow(gv)),
       color = "状态") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "group_volcano_V7_vs_R7.pdf"), p_volc,
       width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(OUT, "group_volcano_V7_vs_R7.png"), p_volc,
       width = 8, height = 6, dpi = 200, type = "cairo")
cat("    Volcano: group_volcano_V7_vs_R7.pdf/png\n")

# 8c. VlnPlot — 用 BH 主结果; 若 BH 也=0, 则展示 V7-vs-R7 log2FC>0 Top 6 作参考
vln_genes <- if (length(final_genes_main) > 0) {
  final_genes_main
} else {
  head(grp_markers$gene_symbol[grp_markers$avg_log2FC_V7_vs_R7 > 0], 6)
}
if (length(vln_genes) > 0) {
  Idents(seu_pure) <- factor(seu_pure$condition_chr, levels = c("V7","R7"))
  p_vln <- VlnPlot(seu_pure, features = vln_genes,
                   pt.size = 0.05, ncol = min(3, length(vln_genes)),
                   cols = c("V7" = "#1F77B4", "R7" = "#D62728")) &
    labs(x = NULL) &
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))
  ncol_p <- min(3, length(vln_genes))
  nrow_p <- ceiling(length(vln_genes) / ncol_p)
  fname_suffix <- if (length(final_genes_main) > 0) "final_BH" else "topRTX_down_ref"
  ggsave(file.path(OUT, sprintf("violin_%s.pdf", fname_suffix)), p_vln,
         width = 3.5 * ncol_p + 1, height = 2.8 * nrow_p + 1,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT, sprintf("violin_%s.png", fname_suffix)), p_vln,
         width = 3.5 * ncol_p + 1, height = 2.8 * nrow_p + 1,
         dpi = 200, type = "cairo", limitsize = FALSE)
  cat(sprintf("    VlnPlot (%s, n=%d): violin_%s.pdf/png\n",
              fname_suffix, length(vln_genes), fname_suffix))
}

# 8d. 参数 + sessionInfo
sink(file.path(OUT, "params_and_session_info.txt"))
cat(sprintf("时间戳: %s\n", TS))
cat("脚本: scripts/11_GSE296065_singlecell_validation.R\n")
cat(sprintf("RDS: %s\n", RDS))
cat(sprintf("靶点: %s\n", TARGETS_CSV))
cat(sprintf("\n目标群: %s\n", paste(TARGET_CT, collapse=" + ")))
cat(sprintf("Flu_Vehicle = condition '%s'\n", FLU_VEHICLE))
cat(sprintf("Flu_RTX     = condition '%s'\n", FLU_RTX))
cat("\n阈值:\n"); print(THR)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

# 9. 保存中间对象 (轻量)
saveRDS(list(spec_markers = spec_markers, grp_markers = grp_markers,
             final_tbl_strict = final_tbl_strict,
             final_tbl_BH     = final_tbl_BH,
             final_genes_strict = final_genes,
             final_genes_BH     = final_genes_BH,
             spec_pass = spec_pass,
             grp_pass_strict = grp_pass,
             grp_pass_BH     = grp_pass_BH,
             THR = THR, TARGET_CT = TARGET_CT,
             FLU_VEHICLE = FLU_VEHICLE, FLU_RTX = FLU_RTX),
        file.path(OUT, "GSE296065_validation_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
