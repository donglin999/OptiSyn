#!/usr/bin/env Rscript
# 分析 1: 5 细胞 × 4 contrasts DEG (limma + BH, FDR<0.05, |log2FC|>=1)
# 分析 2: Mo 剂量梯度 DEG 单调性验证
# 输出: DEG 全表 (long format) + 数量统计表 + 2 张图

suppressPackageStartupMessages({
  library(limma)
  library(illuminaMousev2.db)
  library(ggplot2)
  library(dplyr)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/DEG_5cell_dose")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 加载预处理后的表达矩阵 ----
expr_log2 <- readRDS(file.path(ROOT, "processed_data/GSE42639/expr_log2.rds"))   # gene × sample
meta      <- readRDS(file.path(ROOT, "processed_data/GSE42639/sample_meta.rds"))

stopifnot(all(meta$sample_id %in% colnames(expr_log2)))
expr_log2 <- expr_log2[, meta$sample_id]

cells <- c("Nh","Am","Ne","Mo","Ly")
contrasts_def <- list(
  TX91_vs_Sham    = c("TX91",       "Sham"),
  PR8_0p6_vs_Sham = c("0.6LD50PR8", "Sham"),
  PR8_10_vs_Sham  = c("10LD50PR8",  "Sham"),
  PR8_100_vs_Sham = c("100LD50PR8", "Sham")
)
THRESH_FC  <- 1
THRESH_FDR <- 0.05

# ---- DEA 主函数 ----
limma_dea <- function(cell, treat, control) {
  idx <- which(meta$cell == cell)
  m_sub <- meta[idx, ]; E_sub <- expr_log2[, idx]
  rename <- c(Sham="Sham","TX91"="TX91","0.6LD50PR8"="PR8_0p6",
              "10LD50PR8"="PR8_10","100LD50PR8"="PR8_100")
  m_sub$tname <- factor(rename[m_sub$treat])
  if (!(rename[treat]   %in% levels(m_sub$tname))) return(NULL)
  if (!(rename[control] %in% levels(m_sub$tname))) return(NULL)
  design <- model.matrix(~0 + m_sub$tname)
  colnames(design) <- levels(m_sub$tname)
  cm  <- makeContrasts(contrasts = sprintf("%s - %s",
                          rename[treat], rename[control]), levels = design)
  fit <- eBayes(contrasts.fit(lmFit(E_sub, design), cm), trend = TRUE)
  tt  <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  tt$gene     <- rownames(tt)
  tt$cell     <- cell
  tt$contrast <- sprintf("%s_vs_%s", rename[treat], rename[control])
  tt$sig      <- tt$adj.P.Val < THRESH_FDR & abs(tt$logFC) >= THRESH_FC
  tt$direction <- ifelse(!tt$sig, "ns",
                  ifelse(tt$logFC > 0, "up", "down"))
  tt[, c("gene","cell","contrast","logFC","AveExpr",
         "P.Value","adj.P.Val","sig","direction")]
}

# ---- 全表生成 ----
cat(">>> 5 cells × 4 contrasts limma DEA\n")
all_DEG <- list()
for (c in cells) {
  for (cn in names(contrasts_def)) {
    tc <- contrasts_def[[cn]]
    r  <- limma_dea(c, tc[1], tc[2])
    if (!is.null(r)) all_DEG[[paste(c, cn, sep="__")]] <- r
  }
}
DEG_full <- do.call(rbind, all_DEG)
rownames(DEG_full) <- NULL
write.table(DEG_full, file.path(OUT, "DEG_full_5cell_4contrast.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  全表写入: %d 行 × %d 列\n", nrow(DEG_full), ncol(DEG_full)))

# 显著 DEG 子集
DEG_sig <- DEG_full[DEG_full$sig, ]
write.table(DEG_sig, file.path(OUT, "DEG_significant_5cell_4contrast.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  显著 DEG (FDR<%.2f, |log2FC|>=%d): %d\n",
            THRESH_FDR, THRESH_FC, nrow(DEG_sig)))

# ---- 数量统计表 ----
count_tab <- DEG_full %>%
  group_by(cell, contrast) %>%
  summarise(n_total = sum(sig),
            n_up    = sum(sig & direction == "up"),
            n_down  = sum(sig & direction == "down"),
            .groups = "drop") %>%
  as.data.frame()
count_tab$cell     <- factor(count_tab$cell,     levels = c("Mo","Ne","Nh","Am","Ly"))
count_tab$contrast <- factor(count_tab$contrast,
  levels = c("TX91_vs_Sham","PR8_0p6_vs_Sham","PR8_10_vs_Sham","PR8_100_vs_Sham"))
count_tab <- count_tab[order(count_tab$cell, count_tab$contrast), ]
write.table(count_tab, file.path(OUT, "DEG_count_5cell_4contrast.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n>>> DEG 数量统计表:\n")
print(count_tab)

# ============================================================================
# 分析 1 的图: 5 细胞 × 4 contrasts DEG 数量分组柱状图
# ============================================================================
plot_df <- count_tab %>%
  tidyr::pivot_longer(cols = c(n_up, n_down),
                      names_to = "direction", values_to = "n") %>%
  mutate(direction = factor(direction, levels = c("n_up","n_down"),
                            labels = c("Up","Down")))

p1 <- ggplot(plot_df, aes(cell, n, fill = direction)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(0.9),
            vjust = -0.3, size = 2.8) +
  facet_wrap(~ contrast, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(Up = "#D7263D", Down = "#1B98E0")) +
  labs(x = "细胞类型", y = "显著 DEG 数 (FDR<0.05, |log2FC|>=1)",
       title = "分析 1: 5 细胞 × 4 contrasts 显著 DEG 数量",
       subtitle = "limma + eBayes + BH (vs Sham 基线)",
       fill = "方向") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E8EEF8"))
ggsave(file.path(OUT, "analysis1_5cell_DEG_barplot.pdf"), p1,
       width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "analysis1_5cell_DEG_barplot.png"), p1,
       width = 11, height = 7, dpi = 200, type = "cairo")
cat("\n>>> 分析 1 图: analysis1_5cell_DEG_barplot.{pdf,png}\n")

# ============================================================================
# 分析 2: Mo 单细胞剂量梯度 — 折线图 + 柱状图 + 单调验证
# ============================================================================
mo_tab <- count_tab[count_tab$cell == "Mo", ]
# 把 contrast 转为剂量数值轴
mo_tab$dose_label <- factor(mo_tab$contrast,
                            levels = c("TX91_vs_Sham","PR8_0p6_vs_Sham",
                                       "PR8_10_vs_Sham","PR8_100_vs_Sham"),
                            labels = c("TX91 (轻症)","PR8 0.6LD50",
                                       "PR8 10LD50","PR8 100LD50 (重症)"))
mo_tab$dose_num <- c(0.1, 0.6, 10, 100)[match(as.character(mo_tab$contrast),
                          c("TX91_vs_Sham","PR8_0p6_vs_Sham","PR8_10_vs_Sham","PR8_100_vs_Sham"))]
write.table(mo_tab, file.path(OUT, "Mo_dose_DEG_count.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 单调性验证
mo_long <- mo_tab[order(mo_tab$dose_num), ]
mono_total <- all(diff(mo_long$n_total) >= 0)
mono_up    <- all(diff(mo_long$n_up)    >= 0)
mono_down  <- all(diff(mo_long$n_down)  >= 0)
cat("\n>>> Mo 剂量梯度单调性验证:\n")
cat(sprintf("  总 DEG  数:  %s  (序列: %s)\n", ifelse(mono_total,"✓ 单调递增","✗ 非单调"),
            paste(mo_long$n_total, collapse = " → ")))
cat(sprintf("  Up      数:  %s  (序列: %s)\n", ifelse(mono_up,   "✓ 单调递增","✗ 非单调"),
            paste(mo_long$n_up,    collapse = " → ")))
cat(sprintf("  Down    数:  %s  (序列: %s)\n", ifelse(mono_down, "✓ 单调递增","✗ 非单调"),
            paste(mo_long$n_down,  collapse = " → ")))

# 折线 + 点
mo_plot_df <- mo_tab %>%
  tidyr::pivot_longer(cols = c(n_total, n_up, n_down),
                      names_to = "metric", values_to = "n") %>%
  mutate(metric = factor(metric,
                         levels = c("n_total","n_up","n_down"),
                         labels = c("总 DEG","Up","Down")))

p2a <- ggplot(mo_plot_df, aes(dose_label, n, group = metric, color = metric)) +
  geom_line(size = 1.1) +
  geom_point(size = 3.5) +
  geom_text(aes(label = n), vjust = -1.0, size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c("总 DEG" = "#444444",
                                "Up"     = "#D7263D",
                                "Down"   = "#1B98E0")) +
  labs(x = "感染剂量梯度", y = "显著 DEG 数",
       title = "分析 2: Mo 剂量依赖性 DEG 数量折线",
       subtitle = sprintf("单调性: 总 %s | Up %s | Down %s",
                          ifelse(mono_total,"✓","✗"),
                          ifelse(mono_up,   "✓","✗"),
                          ifelse(mono_down, "✓","✗")),
       color = "类型") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  expand_limits(y = max(mo_plot_df$n) * 1.15)

p2b <- ggplot(mo_plot_df, aes(dose_label, n, fill = metric)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(0.9),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("总 DEG" = "#888888",
                               "Up"     = "#D7263D",
                               "Down"   = "#1B98E0")) +
  labs(x = "感染剂量梯度", y = "显著 DEG 数",
       title = "分析 2: Mo 剂量梯度 DEG 数量柱状",
       fill = "类型") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  expand_limits(y = max(mo_plot_df$n) * 1.10)

# 拼接折线 + 柱状
suppressPackageStartupMessages(library(patchwork))
p2 <- p2a / p2b + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "analysis2_Mo_dose_DEG_lineplot.pdf"), p2,
       width = 9, height = 9, device = cairo_pdf)
ggsave(file.path(OUT, "analysis2_Mo_dose_DEG_lineplot.png"), p2,
       width = 9, height = 9, dpi = 200, type = "cairo")
cat("\n>>> 分析 2 图: analysis2_Mo_dose_DEG_lineplot.{pdf,png}\n")

# ---- Mo 4 组 DEG 全表（仅 Mo, 4 contrasts, 完整 logFC/padj） ----
mo_full <- DEG_full[DEG_full$cell == "Mo", ]
write.table(mo_full, file.path(OUT, "Mo_dose_DEG_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n>>> Mo 4 组 DEG 全表: %d 行\n", nrow(mo_full)))

cat(sprintf("\n>>> 全部输出: %s\n", OUT))
