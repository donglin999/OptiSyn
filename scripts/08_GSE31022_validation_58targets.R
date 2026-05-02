#!/usr/bin/env Rscript
# GSE31022 跨毒株时序验证 58 个 Mo 严选靶点
#
# 输入:
#   processed_data/GSE31022/expr_log2.rds  + sample_meta.rds
#   results/Mo_protective/strict/core_targets_strict.tsv  (58 基因)
#
# 验证逻辑:
#   1) 是否覆盖: 基因在 GSE31022 表达矩阵中存在
#   2) 跨毒株保守 (Cross-strain conserved):
#        任一 Day1-6 vs Ctrl 的 log2FC <= -0.5 且 FDR < 0.10
#        (与 GSE42639 同方向: 下调; 阈值放宽因为跨毒株/不同时间窗)
#   3) 全周期进行性下调 (Progressive down, D1-D6):
#        D1-D6 全周期趋势 slope < 0 且 FDR < 0.10
#        (急性期定义改为 D1-D6 全周期, 涵盖病毒复制顶峰 + 持续损伤期)
#   4) 双标 (湿实验最优先): 同时满足 2 + 3
#
# 输出: results/Mo_protective/strict/GSE31022_validation/

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(RColorBrewer)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/Mo_protective/strict/GSE31022_validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. 数据加载 ----------------------------------------------------------------
expr <- readRDS(file.path(ROOT, "processed_data/GSE31022/expr_log2.rds"))
meta <- readRDS(file.path(ROOT, "processed_data/GSE31022/sample_meta.rds"))
strict <- read.delim(file.path(ROOT,
                     "results/Mo_protective/strict/core_targets_strict.tsv"),
                     stringsAsFactors = FALSE)
genes58 <- strict$gene_symbol
cat(sprintf(">>> Mo 严选 58 靶点 / GSE31022 基因总数 %d\n", nrow(expr)))

# 2. 覆盖度 ------------------------------------------------------------------
covered <- genes58 %in% rownames(expr)
cat(sprintf(">>> 在 GSE31022 中覆盖到的: %d / %d\n", sum(covered), length(genes58)))
miss <- genes58[!covered]
if (length(miss) > 0) {
  cat("    未覆盖: ", paste(miss, collapse=", "), "\n")
}
g_have <- genes58[covered]
E <- expr[g_have, , drop = FALSE]
m <- meta

# 3. limma 模型 1: 各 Day vs Ctrl --------------------------------------------
m$group <- factor(m$group, levels = c("Ctrl","D1","D2","D3","D4","D5","D6"))
designAll <- model.matrix(~0 + m$group)
colnames(designAll) <- levels(m$group)
fitAll <- lmFit(expr, designAll)  # 用全部基因拟合 -> FDR 校正基于全基因
cm <- makeContrasts(
  D1_vs_Ctrl = D1 - Ctrl,
  D2_vs_Ctrl = D2 - Ctrl,
  D3_vs_Ctrl = D3 - Ctrl,
  D4_vs_Ctrl = D4 - Ctrl,
  D5_vs_Ctrl = D5 - Ctrl,
  D6_vs_Ctrl = D6 - Ctrl,
  levels = designAll
)
fit2 <- eBayes(contrasts.fit(fitAll, cm), trend = TRUE)

get_tt <- function(coef) {
  t <- topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  t$gene <- rownames(t); rownames(t) <- NULL
  t[, c("gene","logFC","AveExpr","P.Value","adj.P.Val")]
}
tt <- list(
  D1 = get_tt("D1_vs_Ctrl"),
  D2 = get_tt("D2_vs_Ctrl"),
  D3 = get_tt("D3_vs_Ctrl"),
  D4 = get_tt("D4_vs_Ctrl"),
  D5 = get_tt("D5_vs_Ctrl"),
  D6 = get_tt("D6_vs_Ctrl")
)

# 4. limma 模型 2: 趋势 (D1-D3 急性期, D1-D6 全程) ---------------------------
trend_fit <- function(days_keep) {
  idx <- which(m$day %in% days_keep)
  m_sub <- m[idx, , drop = FALSE]
  E_sub <- expr[, idx, drop = FALSE]
  m_sub$severity <- as.numeric(m_sub$day)
  des <- model.matrix(~ m_sub$severity)
  colnames(des) <- c("Intercept","severity")
  fitS <- eBayes(lmFit(E_sub, des), trend = TRUE)
  tts <- topTable(fitS, coef = "severity", number = Inf, sort.by = "none")
  tts$gene <- rownames(tts); rownames(tts) <- NULL
  colnames(tts)[colnames(tts) == "logFC"] <- "slope"
  tts[, c("gene","slope","P.Value","adj.P.Val")]
}
tt_acute <- trend_fit(days_keep = 1:3)   # D1-D3
tt_full  <- trend_fit(days_keep = 1:6)   # D1-D6

# 5. 组均值 + SEM (用于绘图) -------------------------------------------------
expr_long_for_gene <- function(g) {
  x <- expr[g, ]
  data.frame(gene = g, sample_id = names(x),
             group = m$group[match(names(x), m$sample_id)],
             day   = m$day  [match(names(x), m$sample_id)],
             expr  = as.numeric(x), stringsAsFactors = FALSE)
}

# 6. 标注 --------------------------------------------------------------------
THR <- list(
  conserved_log2fc_max = -0.5,
  conserved_fdr        = 0.10,
  prog_slope_max       = 0,        # D1-D6 全周期 slope < 0
  prog_fdr             = 0.10
)

ann_one <- function(g) {
  if (!(g %in% rownames(expr))) {
    return(data.frame(gene_symbol = g, covered = FALSE,
                      D1_logFC = NA, D1_FDR = NA, D2_logFC = NA, D2_FDR = NA,
                      D3_logFC = NA, D3_FDR = NA, D4_logFC = NA, D4_FDR = NA,
                      D5_logFC = NA, D5_FDR = NA, D6_logFC = NA, D6_FDR = NA,
                      best_logFC_anyDay = NA, best_FDR_anyDay = NA,
                      acute_slope_D1_D3 = NA, acute_FDR_D1_D3 = NA,
                      slope_D1_D6 = NA, FDR_D1_D6 = NA,
                      cross_strain_conserved = FALSE,
                      progressive_down_D1_D6 = FALSE,
                      both = FALSE,
                      stringsAsFactors = FALSE))
  }
  vals <- list()
  for (d in c("D1","D2","D3","D4","D5","D6")) {
    r <- tt[[d]][match(g, tt[[d]]$gene), ]
    vals[[paste0(d,"_logFC")]] <- r$logFC
    vals[[paste0(d,"_FDR")]]   <- r$adj.P.Val
  }
  lfc_v <- c(vals$D1_logFC, vals$D2_logFC, vals$D3_logFC,
             vals$D4_logFC, vals$D5_logFC, vals$D6_logFC)
  fdr_v <- c(vals$D1_FDR, vals$D2_FDR, vals$D3_FDR,
             vals$D4_FDR, vals$D5_FDR, vals$D6_FDR)
  # 任一 Day 满足下调
  conserved <- any(lfc_v <= THR$conserved_log2fc_max &
                   fdr_v < THR$conserved_fdr, na.rm = TRUE)
  best_idx  <- which.min(lfc_v)
  best_lfc  <- lfc_v[best_idx]
  best_fdr  <- fdr_v[best_idx]

  acute_r <- tt_acute[match(g, tt_acute$gene), ]
  full_r  <- tt_full [match(g, tt_full$gene),  ]
  prog_pd <- (!is.na(full_r$slope)) &&
             (full_r$slope < THR$prog_slope_max) &&
             (full_r$adj.P.Val < THR$prog_fdr)

  data.frame(
    gene_symbol = g, covered = TRUE,
    D1_logFC = vals$D1_logFC, D1_FDR = vals$D1_FDR,
    D2_logFC = vals$D2_logFC, D2_FDR = vals$D2_FDR,
    D3_logFC = vals$D3_logFC, D3_FDR = vals$D3_FDR,
    D4_logFC = vals$D4_logFC, D4_FDR = vals$D4_FDR,
    D5_logFC = vals$D5_logFC, D5_FDR = vals$D5_FDR,
    D6_logFC = vals$D6_logFC, D6_FDR = vals$D6_FDR,
    best_logFC_anyDay = best_lfc, best_FDR_anyDay = best_fdr,
    acute_slope_D1_D3 = acute_r$slope, acute_FDR_D1_D3 = acute_r$adj.P.Val,
    slope_D1_D6      = full_r$slope,  FDR_D1_D6        = full_r$adj.P.Val,
    cross_strain_conserved   = conserved,
    progressive_down_D1_D6   = prog_pd,
    both = conserved & prog_pd,
    stringsAsFactors = FALSE
  )
}
ann <- do.call(rbind, lapply(genes58, ann_one))

# 加上 GSE42639 strict 类别字段方便对照
ann$func_class_GSE42639 <- strict$func_class[match(ann$gene_symbol, strict$gene_symbol)]
ann$slope_GSE42639      <- strict$trend_slope[match(ann$gene_symbol, strict$gene_symbol)]

# 排序: 双标 > 仅 conserved > 仅 progressive > 其它; 内排 best_logFC 升 (越负越前)
ann$rank_grp <- with(ann,
  ifelse(both, 1L,
  ifelse(cross_strain_conserved, 2L,
  ifelse(progressive_down_D1_D6, 3L, 4L))))
ann <- ann[order(ann$rank_grp, ann$best_logFC_anyDay, ann$slope_D1_D6), ]
rownames(ann) <- NULL

write.table(ann, file.path(OUT, "GSE31022_validation_58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(ann, file.path(OUT, "GSE31022_validation_58.csv"),
          quote = TRUE, row.names = FALSE)

# 子集表
both_tbl  <- ann[ann$both, ]
cons_tbl  <- ann[ann$cross_strain_conserved, ]
prog_tbl  <- ann[ann$progressive_down_D1_D6, ]
write.table(both_tbl, file.path(OUT, "GSE31022_both_conserved_progressive.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cons_tbl, file.path(OUT, "GSE31022_cross_strain_conserved.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prog_tbl, file.path(OUT, "GSE31022_progressive_down_D1_D6.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n>>> 验证结果分布 (n=%d):\n", nrow(ann)));
cat(sprintf("    覆盖到的:                       %d\n", sum(ann$covered)))
cat(sprintf("    跨毒株保守 (any Day<=-0.5 & FDR<%.2f): %d\n",
            THR$conserved_fdr, sum(ann$cross_strain_conserved)))
cat(sprintf("    全周期 D1-D6 进行性下调 (slope<0 & FDR<%.2f): %d\n",
            THR$prog_fdr, sum(ann$progressive_down_D1_D6)))
cat(sprintf("    双标 (湿实验最优先):            %d\n", sum(ann$both)))

# 7. 流程日志 + 参数 ---------------------------------------------------------
sink(file.path(OUT, "validation_summary.txt"))
cat("=== GSE31022 跨毒株时序验证 (Mo 严选 58 靶点) ===\n")
cat(sprintf("时间: %s\n", TS))
cat("数据集: GSE31022 (C57BL/6 mouse H1N1, lung, D1-D6, 3 Ctrl + 3*6 = 21 样本)\n")
cat("源靶点: results/Mo_protective/strict/core_targets_strict.tsv (n=58)\n\n")

cat("阈值:\n"); print(THR)

cat(sprintf("\n%-50s %d / %d\n", "覆盖到", sum(ann$covered), nrow(ann)))
cat(sprintf("%-50s %d\n", "跨毒株保守 (cross_strain_conserved)",
            sum(ann$cross_strain_conserved)))
cat(sprintf("%-50s %d\n", "全周期 D1-D6 进行性下调 (progressive_down_D1_D6)",
            sum(ann$progressive_down_D1_D6)))
cat(sprintf("%-50s %d\n", "双标 (both, 湿实验最优先)", sum(ann$both)))

cat("\n=== 双标基因 (both) ===\n")
print(both_tbl[, c("gene_symbol","best_logFC_anyDay","best_FDR_anyDay",
                   "slope_D1_D6","FDR_D1_D6",
                   "slope_GSE42639","func_class_GSE42639")],
      row.names = FALSE)

cat("\n=== 仅跨毒株保守 (排除已双标) ===\n")
print(cons_tbl[!cons_tbl$both,
               c("gene_symbol","best_logFC_anyDay","best_FDR_anyDay",
                 "slope_D1_D6","FDR_D1_D6",
                 "func_class_GSE42639")], row.names = FALSE)

cat("\n=== 仅 D1-D6 全周期进行性下调 (排除已双标) ===\n")
print(prog_tbl[!prog_tbl$both,
               c("gene_symbol","slope_D1_D6","FDR_D1_D6",
                 "best_logFC_anyDay","best_FDR_anyDay",
                 "func_class_GSE42639")], row.names = FALSE)
sink()
cat(readLines(file.path(OUT, "validation_summary.txt")), sep="\n")

# 8. 58 panel 折线图 ---------------------------------------------------------
g_plot <- ann$gene_symbol[ann$covered]
plot_df <- do.call(rbind, lapply(g_plot, expr_long_for_gene))

agg <- aggregate(expr ~ gene + group, data = plot_df,
                 FUN = function(v) c(mean = mean(v), sem = sd(v)/sqrt(length(v))))
agg <- data.frame(gene = agg$gene, group = agg$group,
                  mean = agg$expr[, "mean"], sem = agg$expr[, "sem"])
agg$group <- factor(agg$group, levels = c("Ctrl","D1","D2","D3","D4","D5","D6"))
# 顺序: 双标 > conserved > acute > others
agg$gene <- factor(agg$gene, levels = ann$gene_symbol[ann$covered])

# 标签 -> 颜色
gene_label <- with(ann,
  ifelse(both, "Both",
  ifelse(cross_strain_conserved, "Conserved",
  ifelse(progressive_down_D1_D6, "Progressive", "Other"))))
names(gene_label) <- ann$gene_symbol
agg$label <- factor(gene_label[as.character(agg$gene)],
                    levels = c("Both","Conserved","Progressive","Other"))

cols <- c(Both = "#D62728", Conserved = "#1F77B4",
          Progressive = "#2CA02C", Other = "grey50")

p <- ggplot(agg, aes(group, mean, group = gene, color = label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.18, linewidth = 0.4) +
  scale_color_manual(values = cols, name = "validation",
                     drop = FALSE) +
  facet_wrap(~ gene, scales = "free_y", ncol = 8) +
  labs(x = NULL, y = "log2 expr (neqc)",
       title = "GSE31022 (H1N1, mouse lung) D1-D6 时序: 58 靶点验证",
       subtitle = sprintf("Both=同时跨毒株保守+D1-D6全周期进行性下调 (n=%d), Conserved=仅跨毒株(n=%d), Progressive=仅全周期下降(n=%d), Other(n=%d)",
                          sum(ann$both),
                          sum(ann$cross_strain_conserved & !ann$both),
                          sum(ann$progressive_down_D1_D6 & !ann$both),
                          sum(!ann$cross_strain_conserved & !ann$progressive_down_D1_D6))) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"))

n_panel <- length(g_plot)
ncol_p <- 8
nrow_p <- ceiling(n_panel / ncol_p)
fig_h <- max(8, nrow_p * 1.4 + 1.5)
ggsave(file.path(OUT, "GSE31022_58_lineplot.pdf"), p,
       width = 16, height = fig_h, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "GSE31022_58_lineplot.png"), p,
       width = 16, height = fig_h, dpi = 200, type = "cairo", limitsize = FALSE)
cat(">>> 58 panel 折线图: GSE31022_58_lineplot.pdf/png\n")

# 9. 双标基因小图 (湿实验汇报) -----------------------------------------------
if (nrow(both_tbl) > 0) {
  agg_b <- agg[agg$gene %in% both_tbl$gene_symbol, ]
  agg_b$gene <- factor(agg_b$gene, levels = both_tbl$gene_symbol)
  p2 <- ggplot(agg_b, aes(group, mean, group = gene)) +
    geom_line(color = "#D62728", linewidth = 0.8) +
    geom_point(color = "#D62728", size = 1.6) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, linewidth = 0.5, color = "grey25") +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    labs(x = NULL, y = "log2 expr (neqc)",
         title = sprintf("GSE31022 跨毒株保守 + D1-D6 全周期进行性下调 (n=%d)",
                         nrow(both_tbl)),
         subtitle = "组均值 +/- SEM, n=3") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT, "GSE31022_both_lineplot.pdf"), p2,
         width = 10, height = max(4, ceiling(nrow(both_tbl)/4) * 2 + 1),
         device = cairo_pdf)
  ggsave(file.path(OUT, "GSE31022_both_lineplot.png"), p2,
         width = 10, height = max(4, ceiling(nrow(both_tbl)/4) * 2 + 1),
         dpi = 200, type = "cairo")
  cat(sprintf(">>> 双标(%d)小图: GSE31022_both_lineplot.pdf/png\n", nrow(both_tbl)))
}

saveRDS(list(ann = ann, tt = tt, tt_acute = tt_acute, tt_full = tt_full,
             THR = THR, plot_df = plot_df, agg = agg),
        file.path(OUT, "GSE31022_validation.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
