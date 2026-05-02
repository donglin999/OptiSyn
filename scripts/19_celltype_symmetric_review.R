#!/usr/bin/env Rscript
# 19_celltype_symmetric_review.R
# 在 GSE42639 的 5 个细胞类型 (Nh/Am/Ne/Mo/Ly) 上对称跑同一份"严选"流程,
# 用来 review "为什么选 Mo (而非别的细胞)" 的争议.
#
# 流程 (与 script 05+06 完全一致, 仅改 cell 子集):
#   Step1 趋势 (TX91=0..PR8_100=3):  slope<0 & FDR<0.10
#   Step2 崩塌 (PR8_100 vs TX91):     log2FC<=-0.585 & FDR<0.05
#   Step3 基础表达:                    TX91 mean >= 全基因 TX91 mean 中位
#   Step4 反向假阳性剔除:              剔 PR8_100 vs TX91 log2FC>=1 & FDR<0.05
#   Step5 排已下调:                    剔 TX91 vs Sham FDR<0.05 & log2FC<0
#   StepA 阶梯下降:                    0p6 logFC<=-0.30 & 10 logFC<=-0.585 (FDR<0.05)
#
# 比较点:
#   - 每个 cell 通过 Step1-5 的候选数 (= 564 在 Mo)
#   - 每个 cell 通过 StepA 严选的候选数 (= 65 在 Mo)
#   - 关键基因 (Cx3cr1, Cx3cl1, Klf2, Arrb1, Gstm1, Cd209a) 的 trend_slope 跨 5 cell 对比
#
# 输出: results/celltype_review/

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

set.seed(42)
ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))

OUT  <- file.path(ROOT, "results/celltype_review")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

THR <- list(
  step1_fdr = 0.10, step2_lfc = -0.585, step2_fdr = 0.05,
  step3_quantile = 0.5, step4_lfc_min = 1, step4_fdr = 0.05,
  step5_fdr = 0.05,
  stepA_lfc_0p6 = -0.30, stepA_lfc_10 = -0.585, stepA_fdr = 0.05
)

GROUP_LEVELS <- c("Sham","TX91","PR8_0p6","PR8_10","PR8_100")
SEV_MAP <- c(TX91=0, PR8_0p6=1, PR8_10=2, PR8_100=3)

# 1. 加载
expr <- readRDS(file.path(ROOT, "processed_data/GSE42639/expr_log2.rds"))
meta <- readRDS(file.path(ROOT, "processed_data/GSE42639/sample_meta.rds"))
treat_R <- c(Sham="Sham", TX91="TX91",
             "0.6LD50PR8"="PR8_0p6","10LD50PR8"="PR8_10","100LD50PR8"="PR8_100")
meta$treat <- factor(treat_R[meta$treat], levels = GROUP_LEVELS)

cells <- c("Nh","Am","Ne","Mo","Ly")

# 2. 单细胞跑流程
run_one_cell <- function(cell) {
  idx <- which(meta$cell == cell)
  m <- meta[idx, , drop = FALSE]; E <- expr[, idx, drop = FALSE]
  t_tab <- table(m$treat)

  # 主模型
  designAll <- model.matrix(~0 + m$treat)
  colnames(designAll) <- levels(m$treat)
  fitAll <- lmFit(E, designAll)
  cm <- makeContrasts(
    TX91_vs_Sham    = TX91 - Sham,
    PR8_100_vs_TX91 = PR8_100 - TX91,
    PR8_0p6_vs_TX91 = PR8_0p6 - TX91,
    PR8_10_vs_TX91  = PR8_10  - TX91,
    levels = designAll
  )
  fit2 <- eBayes(contrasts.fit(fitAll, cm), trend = TRUE)
  get_tt <- function(coef) {
    t <- topTable(fit2, coef = coef, number = Inf, sort.by = "none")
    t$gene <- rownames(t); rownames(t) <- NULL
    t[, c("gene","logFC","adj.P.Val")]
  }
  ttSham <- get_tt("TX91_vs_Sham")
  tt100  <- get_tt("PR8_100_vs_TX91")
  tt06   <- get_tt("PR8_0p6_vs_TX91")
  tt10   <- get_tt("PR8_10_vs_TX91")

  # 趋势
  inf_idx <- which(m$treat %in% names(SEV_MAP))
  m_inf <- m[inf_idx,]; E_inf <- E[, inf_idx]
  m_inf$severity <- SEV_MAP[as.character(m_inf$treat)]
  designSev <- model.matrix(~ m_inf$severity)
  colnames(designSev) <- c("Intercept","severity")
  fitSev <- eBayes(lmFit(E_inf, designSev), trend = TRUE)
  ttSev <- topTable(fitSev, coef = "severity", number = Inf, sort.by = "none")
  ttSev$gene <- rownames(ttSev); rownames(ttSev) <- NULL
  colnames(ttSev)[colnames(ttSev) == "logFC"] <- "slope"
  ttSev <- ttSev[, c("gene","slope","adj.P.Val")]

  # TX91 mean + 中位阈
  tx91_idx <- which(m$treat == "TX91")
  tx91_mean <- rowMeans(E[, tx91_idx, drop = FALSE])
  median_thr <- as.numeric(quantile(tx91_mean, THR$step3_quantile))

  # Step 1-5
  S1 <- ttSev$gene[ttSev$slope < 0 & ttSev$adj.P.Val < THR$step1_fdr]
  g_collapse <- tt100$gene[tt100$logFC <= THR$step2_lfc & tt100$adj.P.Val < THR$step2_fdr]
  S2 <- intersect(S1, g_collapse)
  keep3 <- tx91_mean[S2] >= median_thr
  S3 <- S2[keep3]
  g_revUp <- tt100$gene[tt100$logFC >= THR$step4_lfc_min & tt100$adj.P.Val < THR$step4_fdr]
  S4 <- setdiff(S3, g_revUp)
  sham_lfc <- ttSham$logFC[match(S4, ttSham$gene)]
  sham_fdr <- ttSham$adj.P.Val[match(S4, ttSham$gene)]
  already_down <- (sham_fdr < THR$step5_fdr) & (sham_lfc < 0)
  already_down[is.na(already_down)] <- FALSE
  S5 <- S4[!already_down]

  # Step A
  m06 <- match(S5, tt06$gene); m10 <- match(S5, tt10$gene)
  pass_06 <- (tt06$logFC[m06] <= THR$stepA_lfc_0p6) & (tt06$adj.P.Val[m06] < THR$stepA_fdr)
  pass_10 <- (tt10$logFC[m10] <= THR$stepA_lfc_10 ) & (tt10$adj.P.Val[m10] < THR$stepA_fdr)
  pass_06[is.na(pass_06)] <- FALSE; pass_10[is.na(pass_10)] <- FALSE
  S_A <- S5[pass_06 & pass_10]

  list(cell = cell, n_total = nrow(E), n_per_group = t_tab,
       S1 = S1, S2 = S2, S3 = S3, S4 = S4, S5 = S5, S_A = S_A,
       ttSev = ttSev, tt100 = tt100, ttSham = ttSham,
       tx91_mean = tx91_mean)
}

cat(">>> 5 细胞类型对称跑流程\n")
res_list <- list()
for (c in cells) {
  cat(sprintf("    [%s]\n", c))
  res_list[[c]] <- run_one_cell(c)
}

# 3. 流程计数表
flow <- data.frame(
  cell = cells,
  S1_trend  = sapply(cells, function(c) length(res_list[[c]]$S1)),
  S2_collapse = sapply(cells, function(c) length(res_list[[c]]$S2)),
  S3_baseline = sapply(cells, function(c) length(res_list[[c]]$S3)),
  S4_revRm    = sapply(cells, function(c) length(res_list[[c]]$S4)),
  S5_protect  = sapply(cells, function(c) length(res_list[[c]]$S5)),
  StepA_strict = sapply(cells, function(c) length(res_list[[c]]$S_A))
)
flow <- flow[order(-flow$StepA_strict), ]
write.table(flow, file.path(OUT, "celltype_flow_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n>>> 5 细胞类型 5 步 + StepA 计数:\n")
print(flow, row.names = FALSE)

# 4. 关键基因的 trend_slope 跨 5 cell
key_genes <- c("Cx3cr1","Cx3cl1","Klf2","Arrb1","Gstm1","Cd209a","Btla","Il16",
               "Insr","Mbp","Bcl11a","Cd177","Tac1","Slc17a6")
key_tbl <- data.frame(gene = key_genes)
for (c in cells) {
  ts <- res_list[[c]]$ttSev
  key_tbl[[paste0("slope_", c)]] <- ts$slope[match(key_genes, ts$gene)]
  key_tbl[[paste0("FDR_", c)]]   <- ts$adj.P.Val[match(key_genes, ts$gene)]
}
write.table(key_tbl, file.path(OUT, "key_genes_slope_5celltypes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n>>> 关键基因 trend_slope (5 细胞) - 仅显示 slope 列:\n")
print(key_tbl[, c("gene", grep("^slope", colnames(key_tbl), value = TRUE))],
      row.names = FALSE)

# 5. StepA 各细胞的 strict 候选交集
cat("\n>>> StepA 严选候选交集 (5 细胞两两):\n")
SA_list <- lapply(cells, function(c) res_list[[c]]$S_A); names(SA_list) <- cells
for (i in seq_along(cells)) {
  for (j in seq_along(cells)) {
    if (i >= j) next
    a <- cells[i]; b <- cells[j]
    inter <- intersect(SA_list[[a]], SA_list[[b]])
    cat(sprintf("    %s ∩ %s = %d  (%s ...)\n",
                a, b, length(inter),
                paste(head(inter, 5), collapse = ", ")))
  }
}
# 5 cells 共有
common5 <- Reduce(intersect, SA_list)
cat(sprintf("\n  5 细胞共有 StepA 候选: %d\n", length(common5)))
if (length(common5) > 0) cat(sprintf("    %s\n", paste(common5, collapse = ", ")))

# 6. 各细胞 Top 10 严选候选
cat("\n>>> 各细胞 Top 10 strict 候选 (按 trend_slope 升序):\n")
for (c in cells) {
  S <- res_list[[c]]$S_A
  ts <- res_list[[c]]$ttSev
  if (length(S) == 0) {
    cat(sprintf("    [%s] 0 候选\n", c)); next
  }
  s_slope <- ts$slope[match(S, ts$gene)]
  ord <- order(s_slope)
  top10 <- S[ord][1:min(10, length(S))]
  cat(sprintf("    [%s] %d 候选, Top 10: %s\n", c, length(S),
              paste(top10, collapse = ", ")))
}

# 7. 可视化: 5 cell × 14 关键基因 slope 热图
plot_df <- data.frame()
for (c in cells) {
  ts <- res_list[[c]]$ttSev
  for (g in key_genes) {
    v <- ts$slope[match(g, ts$gene)]
    plot_df <- rbind(plot_df, data.frame(cell = c, gene = g, slope = v,
                                          stringsAsFactors = FALSE))
  }
}
plot_df$cell <- factor(plot_df$cell, levels = c("Mo","Ne","Am","Nh","Ly"))
plot_df$gene <- factor(plot_df$gene, levels = key_genes)
cap <- 1.5
plot_df$slope_cap <- pmax(pmin(plot_df$slope, cap), -cap)
p <- ggplot(plot_df, aes(cell, gene, fill = slope_cap)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", slope)),
            size = 3, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-cap, cap),
                       oob = scales::squish, name = "trend_slope") +
  labs(x = NULL, y = NULL,
       title = "关键基因 trend_slope 跨 5 个细胞类型对比",
       subtitle = "slope < 0 = 随感染加重持续下降; 仅 Mo 大量基因显著负 slope") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
ggsave(file.path(OUT, "key_genes_slope_heatmap.pdf"), p,
       width = 6, height = 6, device = cairo_pdf)
ggsave(file.path(OUT, "key_genes_slope_heatmap.png"), p,
       width = 6, height = 6, dpi = 200, type = "cairo")

# 8. 流程计数对比条形图
flow_long <- reshape(flow, varying = c("S1_trend","S2_collapse","S3_baseline",
                                        "S4_revRm","S5_protect","StepA_strict"),
                      v.names = "n", timevar = "step",
                      times = c("S1_trend","S2_collapse","S3_baseline",
                                "S4_revRm","S5_protect","StepA_strict"),
                      direction = "long")
flow_long$step <- factor(flow_long$step,
                          levels = c("S1_trend","S2_collapse","S3_baseline",
                                      "S4_revRm","S5_protect","StepA_strict"))
flow_long$cell <- factor(flow_long$cell, levels = c("Mo","Ne","Nh","Am","Ly"))
p_bar <- ggplot(flow_long, aes(step, n, fill = cell)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(0.9),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Mo"="#D62728", "Ne"="#FF9933", "Nh"="#2CA02C",
                                "Am"="#1F77B4", "Ly"="#9467BD")) +
  labs(x = NULL, y = "通过基因数",
       title = "5 细胞类型在筛选流程各步通过数对比",
       subtitle = "Mo 在所有步骤都显著领先 (除 Step3 基础表达类似)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "celltype_flow_counts_bar.pdf"), p_bar,
       width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(OUT, "celltype_flow_counts_bar.png"), p_bar,
       width = 9, height = 5, dpi = 200, type = "cairo")

# 9. 完整候选名单
for (c in cells) {
  S <- res_list[[c]]$S_A
  if (length(S) > 0) {
    ts <- res_list[[c]]$ttSev
    df <- data.frame(gene = S,
                     trend_slope = ts$slope[match(S, ts$gene)],
                     trend_FDR   = ts$adj.P.Val[match(S, ts$gene)])
    df <- df[order(df$trend_slope), ]
    write.table(df, file.path(OUT, sprintf("strict_candidates_%s.tsv", c)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    writeLines("(no candidates)", file.path(OUT, sprintf("strict_candidates_%s.tsv", c)))
  }
}

# 10. RDS
saveRDS(res_list, file.path(OUT, "celltype_review.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
