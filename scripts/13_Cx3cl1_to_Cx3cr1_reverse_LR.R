#!/usr/bin/env Rscript
# 反向 LR 验证: Vagus (Cx3cl1) -> Mo (Cx3cr1)
#
# 生物学基础:
#   Cx3cl1 (Fractalkine) 是神经元 / 上皮细胞分泌的趋化因子, 经典 vagal sensory
#   neuron 信号分子之一. Cx3cr1 是 Mo / 巨噬 / microglia / NK 高度特异受体.
#   反向 LR (vagus 神经元分泌配体 -> Mo 接收) 比 Mo 释放配体 -> 神经元接收
#   更符合 Cx3cl1-Cx3cr1 通路的真实方向.
#
# 验证内容:
#   [A] Cx3cl1 在 GSE161878 (bulk vagal ganglia) 的表达
#       - CPM 在每个样本里
#       - IAV vs Mock 是否稳定 / 上调 / 下调
#       - 是否在 vagal 表达池 (CPM>1 在 >=20% 样本)
#   [B] Cx3cr1 在 GSE42639 Mo (bulk) + GSE296065 (单细胞) 的表达
#       - Mo TX91 mean (轻症维持)
#       - Mo 严选 trend_slope, FDR (重症崩塌情况)
#       - 单细胞 Mo Ly6c2+ 中的 pct
#   [C] Hill 结合概率
#       L = Cx3cl1 in vagus IAV mean (rank-norm)
#       R = Cx3cr1 in Mo TX91 mean (rank-norm)
#       P = H(L) * H(R)
#   [D] 趋势匹配
#       L_dir: vagus Cx3cl1 IAV vs Mock 方向 (上调 / 下调 / 稳定)
#       R_dir: Mo Cx3cr1 TX91 -> PR8_100 方向 (已知严选: 崩塌)
#       Match: 信号通路是否一致
#   [E] Cx3cr1 下游通路: 趋化 / 抗炎 / 稳态
#       搜 GO term -> Mo 严选 58 命中数 -> 趋势汇总
#
# 输出: results/CellChat/Cx3cl1_to_Cx3cr1_reverse/

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(org.Mm.eg.db)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/CellChat/Cx3cl1_to_Cx3cr1_reverse")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(K = 0.5, n_hill = 2,
               cpm_thresh = 1, sample_frac_min = 0.20)

# ============================================================================
# 1. 加载 GSE161878 (vagal) raw counts -> CPM
# ============================================================================
cnt <- read.delim(file.path(ROOT, "raw_data/00_raw_data/GSE161878_gene_counts.txt"),
                  check.names = FALSE, stringsAsFactors = FALSE)
meta161 <- read.delim(file.path(ROOT,
            "raw_data/00_raw_data/GSE161878_sample_metadata.txt"),
            stringsAsFactors = FALSE)
genes_e <- as.character(cnt$GeneID)
M_raw <- as.matrix(cnt[, -1, drop = FALSE]); rownames(M_raw) <- genes_e
meta161 <- meta161[match(colnames(M_raw), meta161$sample_id), , drop = FALSE]

# 主分析: 剔 D-suffix
keep_main <- !grepl("D$", meta161$sample_id)
M_raw <- M_raw[, keep_main, drop = FALSE]
meta161 <- meta161[keep_main, , drop = FALSE]
lib_size <- colSums(M_raw)
cpm <- t(t(M_raw) / lib_size) * 1e6
log2cpm <- log2(cpm + 1)

ann_g <- read.delim(file.path(ROOT, "processed_data/GSE161878/gene_annotation.tsv"),
                    stringsAsFactors = FALSE)
sym <- setNames(ann_g$SYMBOL, ann_g$ENTREZID)
g_sym <- sym[rownames(log2cpm)]
keep_v <- !is.na(g_sym) & nchar(g_sym) > 0
log2cpm <- log2cpm[keep_v,]; cpm <- cpm[keep_v,]; g_sym <- g_sym[keep_v]
ord <- order(g_sym, -rowMeans(log2cpm))
log2cpm <- log2cpm[ord,]; cpm <- cpm[ord,]; g_sym <- g_sym[ord]
keep_first <- !duplicated(g_sym)
log2cpm <- log2cpm[keep_first,]; cpm <- cpm[keep_first,]
rownames(log2cpm) <- g_sym[keep_first]; rownames(cpm) <- g_sym[keep_first]

iav_idx  <- which(meta161$group == "IAV")
mock_idx <- which(meta161$group == "Mock")
cat(sprintf(">>> GSE161878 主分析 %d 样本 (Mock=%d, IAV=%d)\n",
            ncol(log2cpm), length(mock_idx), length(iav_idx)))

# 表达池
n_samples <- ncol(log2cpm)
n_min <- ceiling(PARAMS$sample_frac_min * n_samples)
detected <- rowSums(cpm > PARAMS$cpm_thresh)
expressed_pool <- rownames(log2cpm)[detected >= n_min]
cat(sprintf(">>> vagal 表达池 (CPM>%g 在 >=%d/%d): %d 基因\n",
            PARAMS$cpm_thresh, n_min, n_samples, length(expressed_pool)))

# ============================================================================
# 2. [A] 验证 Cx3cl1 在迷走神经的表达
# ============================================================================
cat("\n>>> [A] Cx3cl1 在 GSE161878 迷走神经的表达\n")
g_L <- "Cx3cl1"

cx3cl1_in_data <- g_L %in% rownames(cpm)
if (!cx3cl1_in_data) stop(sprintf("%s 不在 GSE161878 矩阵中", g_L))
cx3cl1_cpm  <- cpm[g_L, ]
cx3cl1_log2 <- log2cpm[g_L, ]

A_summary <- data.frame(
  sample_id = colnames(cpm),
  group     = meta161$group,
  CPM       = as.numeric(cx3cl1_cpm),
  log2CPM   = as.numeric(cx3cl1_log2),
  stringsAsFactors = FALSE
)
A_summary <- A_summary[order(A_summary$group, -A_summary$CPM), ]
write.table(A_summary, file.path(OUT, "Cx3cl1_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

L_iav  <- mean(cx3cl1_log2[iav_idx])
L_mock <- mean(cx3cl1_log2[mock_idx])
L_iav_cpm  <- mean(cx3cl1_cpm [iav_idx])
L_mock_cpm <- mean(cx3cl1_cpm [mock_idx])
n_above <- sum(cx3cl1_cpm > PARAMS$cpm_thresh)
in_pool <- g_L %in% expressed_pool

# DEA 信息 (主分析)
dea <- read.delim(file.path(ROOT, "results/GSE161878/GSE161878_DEA_full.tsv"),
                  stringsAsFactors = FALSE)
dea_l <- dea[!is.na(dea$SYMBOL) & dea$SYMBOL == g_L, ]

cat(sprintf("    CPM 范围: %.2f - %.2f (n=%d 样本 CPM>%g)\n",
            min(cx3cl1_cpm), max(cx3cl1_cpm), n_above, PARAMS$cpm_thresh))
cat(sprintf("    Mock 平均 log2CPM: %.3f (CPM=%.2f)\n", L_mock, L_mock_cpm))
cat(sprintf("    IAV  平均 log2CPM: %.3f (CPM=%.2f)\n", L_iav,  L_iav_cpm))
cat(sprintf("    在 vagal 表达池: %s (>=20%% 样本 CPM>%g)\n",
            in_pool, PARAMS$cpm_thresh))
if (nrow(dea_l) > 0) {
  cat(sprintf("    DESeq2 IAV vs Mock: log2FC=%.3f, FDR=%.3g, baseMean=%.1f\n",
              dea_l$log2FC_shrink, dea_l$FDR, dea_l$baseMean))
}

# ============================================================================
# 3. [B] Cx3cr1 在 Mo 细胞的表达
# ============================================================================
cat("\n>>> [B] Cx3cr1 在 Mo (GSE42639 + GSE296065) 的表达\n")
g_R <- "Cx3cr1"

mo_res <- readRDS(file.path(ROOT, "results/Mo_protective/Mo_protective_results.rds"))
mo_E <- mo_res$expr
mo_m <- mo_res$meta
strict <- read.delim(file.path(ROOT,
              "results/Mo_protective/strict/core_targets_strict.tsv"),
              stringsAsFactors = FALSE)

mo_groups <- c("Sham","TX91","PR8_0p6","PR8_10","PR8_100")
B_summary <- data.frame(group = mo_groups,
                        log2_mean = vapply(mo_groups, function(gp) {
                          idx <- which(mo_m$treat == gp)
                          mean(mo_E[g_R, idx])
                        }, numeric(1)),
                        stringsAsFactors = FALSE)
write.table(B_summary, file.path(OUT, "Cx3cr1_Mo_5groups.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(B_summary)

# 严选 58 表里 Cx3cr1 信息
strict_cx3cr1 <- strict[strict$gene_symbol == g_R, ]
cat("\n    Mo 严选 strict 表中的 Cx3cr1:\n")
print(strict_cx3cr1[, c("gene_symbol","trend_slope","trend_FDR",
                        "log2FC_100_vs_TX91","FDR_100_vs_TX91",
                        "func_class")])

# 单细胞 GSE296065 中 Cx3cr1
sc_spec_path <- file.path(ROOT,
       "results/GSE296065_singlecell/step2_specificity_FindMarkers.tsv")
sc_spec <- read.delim(sc_spec_path, stringsAsFactors = FALSE)
sc_cx3cr1 <- sc_spec[sc_spec$gene_symbol == g_R, ]
cat("\n    GSE296065 单细胞 Mo/Mo-Macs Target vs Background:\n")
print(sc_cx3cr1[, c("gene_symbol","avg_log2FC_target_vs_bg","p_val_adj",
                    "pct_target","pct_bg")])

R_TX91 <- B_summary$log2_mean[B_summary$group == "TX91"]
R_PR8_100 <- B_summary$log2_mean[B_summary$group == "PR8_100"]
R_Sham <- B_summary$log2_mean[B_summary$group == "Sham"]

# ============================================================================
# 4. [C] Hill 结合概率 (rank-norm 跨数据集)
# ============================================================================
cat("\n>>> [C] CellChat 风格 Hill 结合概率 (反向: vagus L -> Mo R)\n")
# vagus L 表达分位 (在 vagus 全基因里)
vagus_iav_full <- rowMeans(log2cpm[, iav_idx, drop = FALSE])
vagus_rn  <- rank(vagus_iav_full, ties.method = "average") / length(vagus_iav_full)
names(vagus_rn) <- names(vagus_iav_full)
L_rn <- vagus_rn[g_L]
# Mo R 表达分位
mo_tx91_mean <- rowMeans(mo_E[, mo_m$treat == "TX91"])
mo_rn <- rank(mo_tx91_mean, ties.method = "average") / length(mo_tx91_mean)
names(mo_rn) <- names(mo_tx91_mean)
R_rn <- mo_rn[g_R]

hill <- function(E, K = PARAMS$K, n = PARAMS$n_hill) E^n / (K^n + E^n)
P <- hill(L_rn) * hill(R_rn)
cat(sprintf("    L (vagus Cx3cl1) raw=%.3f log2CPM, rn=%.3f, Hill=%.3f\n",
            L_iav, L_rn, hill(L_rn)))
cat(sprintf("    R (Mo Cx3cr1)    raw=%.3f log2 expr, rn=%.3f, Hill=%.3f\n",
            R_TX91, R_rn, hill(R_rn)))
cat(sprintf("    P_LR = H(L)*H(R) = %.3f\n", P))

# 对比之前 v2 中 Cd209a-Ceacam1 的 prob
cat("    [对比] v2 唯一通讯对 Cd209a-Ceacam1 prob = 0.252\n")

# ============================================================================
# 5. [D] 趋势匹配
# ============================================================================
cat("\n>>> [D] 趋势匹配\n")
L_log2FC_IAV_Mock <- if (nrow(dea_l) > 0) dea_l$log2FC_shrink else NA
L_FDR_IAV_Mock    <- if (nrow(dea_l) > 0) dea_l$FDR else NA
L_dir_str <- if (is.na(L_log2FC_IAV_Mock)) {
  "无 DEA"
} else if (abs(L_log2FC_IAV_Mock) < 0.1) {
  "稳定 (|log2FC|<0.1)"
} else if (L_log2FC_IAV_Mock > 0) {
  "上调"
} else {
  "下调"
}
R_dir_str <- "重症崩塌 (slope<0)"   # 来自 GSE42639 Mo 严选 (slope=-0.796)

# 信号传导一致性: vagus L 在感染中 (稳定 / 上调) + Mo R 在重症崩塌
# = 配体在外周保留, 但接收端关闭 -> "接收端衰退" 模式
mode <- if (is.na(L_log2FC_IAV_Mock)) {
  "L 无方向数据"
} else if (L_log2FC_IAV_Mock >= 0 || abs(L_log2FC_IAV_Mock) < 0.3) {
  "信号源稳定/上调 + 接收端崩塌  -> RECEIVER_LOSS (受体丢失主导)"
} else {
  "信号源下调 + 接收端崩塌        -> SOURCE_RECEIVER_CO_LOSS (双向衰退)"
}
cat(sprintf("    L (Cx3cl1) IAV vs Mock log2FC=%.3f (FDR=%.3g) -> %s\n",
            L_log2FC_IAV_Mock, L_FDR_IAV_Mock, L_dir_str))
cat(sprintf("    R (Cx3cr1) Mo TX91->PR8_100 trend_slope=%.3f (FDR=%.3g) -> %s\n",
            strict_cx3cr1$trend_slope[1], strict_cx3cr1$trend_FDR[1], R_dir_str))
cat(sprintf("    模式: %s\n", mode))

# ============================================================================
# 6. [E] Cx3cr1 下游通路: 趋化 / 抗炎 / 稳态
# ============================================================================
cat("\n>>> [E] Cx3cr1 下游通路检索\n")
go_groups_E <- list(
  Chemotaxis = c(
    "GO:0006935",  # chemotaxis
    "GO:0070098",  # chemokine-mediated signaling pathway
    "GO:0042330",  # taxis
    "GO:0030595",  # leukocyte chemotaxis
    "GO:0050900",  # leukocyte migration
    "GO:0030593"   # neutrophil chemotaxis
  ),
  AntiInflammation = c(
    "GO:0050728",  # negative regulation of inflammatory response
    "GO:0001818",  # negative regulation of cytokine production
    "GO:0050777",  # negative regulation of immune response
    "GO:0032720",  # negative regulation of TNF production
    "GO:0032689",  # negative regulation of IFN-gamma production
    "GO:0002862"   # negative regulation of inflammatory response to antigenic stimulus
  ),
  Homeostasis = c(
    "GO:0042592",  # homeostatic process
    "GO:0001890",  # placental development (有些 -> 排除)
    "GO:0048871",  # multicellular organismal homeostasis
    "GO:0030003",  # cellular cation homeostasis
    "GO:0050801",  # ion homeostasis
    "GO:0048872"   # homeostasis of number of cells
  )
)
go2syms <- function(go_ids) {
  out <- tryCatch({
    suppressMessages(suppressWarnings(
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = go_ids, keytype = "GOALL",
                            columns = c("SYMBOL"))
    ))
  }, error = function(e) {
    suppressMessages(suppressWarnings(
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = go_ids, keytype = "GO",
                            columns = c("SYMBOL"))
    ))
  })
  unique(out$SYMBOL[!is.na(out$SYMBOL)])
}
group_genes <- lapply(go_groups_E, go2syms)
for (nm in names(group_genes)) {
  cat(sprintf("    %-20s GO 覆盖: %d 鼠 SYMBOL\n", nm, length(group_genes[[nm]])))
}

genes58 <- strict$gene_symbol
hits <- lapply(group_genes, function(set) intersect(genes58, set))
cat("\n    58 严选命中:\n")
for (nm in names(hits)) {
  if (length(hits[[nm]]) > 0) {
    cat(sprintf("    [%s] %d -> %s\n", nm, length(hits[[nm]]),
                paste(hits[[nm]], collapse = ", ")))
  } else {
    cat(sprintf("    [%s] 0\n", nm))
  }
}

# 把 hits 整理为表 + 各基因 GSE42639 趋势
hit_tbl <- data.frame()
for (nm in names(hits)) {
  if (length(hits[[nm]]) == 0) next
  for (g in hits[[nm]]) {
    row <- strict[strict$gene_symbol == g, ]
    hit_tbl <- rbind(hit_tbl, data.frame(
      gene_symbol = g,
      pathway = nm,
      trend_slope_GSE42639   = row$trend_slope[1],
      trend_FDR_GSE42639     = row$trend_FDR[1],
      log2FC_100_vs_TX91     = row$log2FC_100_vs_TX91[1],
      FDR_100_vs_TX91        = row$FDR_100_vs_TX91[1],
      func_class             = row$func_class[1],
      stringsAsFactors = FALSE
    ))
  }
}
write.table(hit_tbl, file.path(OUT, "Cx3cr1_downstream_pathways_in_58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n    通路命中明细写入 Cx3cr1_downstream_pathways_in_58.tsv\n"))

# ============================================================================
# 7. 汇总报告
# ============================================================================
sink(file.path(OUT, "summary.txt"))
cat("=== Cx3cl1 -> Cx3cr1 反向 LR 验证 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/13_Cx3cl1_to_Cx3cr1_reverse_LR.R\n\n")

cat("===== [A] Cx3cl1 在 GSE161878 迷走神经表达 =====\n")
cat(sprintf("  CPM 范围 %.2f - %.2f, n=%d/%d 样本 CPM>%g\n",
            min(cx3cl1_cpm), max(cx3cl1_cpm), n_above, n_samples, PARAMS$cpm_thresh))
cat(sprintf("  Mock log2CPM mean: %.3f  ; IAV log2CPM mean: %.3f\n", L_mock, L_iav))
cat(sprintf("  在 vagal 表达池: %s\n", in_pool))
if (nrow(dea_l) > 0) {
  cat(sprintf("  DESeq2 IAV vs Mock: log2FC=%.3f, FDR=%.3g, baseMean=%.1f\n",
              dea_l$log2FC_shrink, dea_l$FDR, dea_l$baseMean))
}

cat("\n===== [B] Cx3cr1 在 Mo 细胞表达 =====\n")
cat("  GSE42639 Mo 五组 log2 mean:\n"); print(B_summary)
cat("\n  Mo 严选 strict 表:\n")
print(strict_cx3cr1[, c("gene_symbol","trend_slope","trend_FDR",
                        "log2FC_100_vs_TX91","FDR_100_vs_TX91","func_class")],
      row.names = FALSE)
cat("\n  GSE296065 单细胞 Target vs Background:\n")
print(sc_cx3cr1[, c("gene_symbol","avg_log2FC_target_vs_bg","p_val_adj",
                    "pct_target","pct_bg")], row.names = FALSE)

cat("\n===== [C] Hill 结合概率 =====\n")
cat(sprintf("  L (vagus Cx3cl1 IAV mean log2CPM=%.3f, rn=%.3f, Hill=%.3f)\n",
            L_iav, L_rn, hill(L_rn)))
cat(sprintf("  R (Mo Cx3cr1 TX91 log2 mean=%.3f, rn=%.3f, Hill=%.3f)\n",
            R_TX91, R_rn, hill(R_rn)))
cat(sprintf("  P_LR = %.3f  (对比 Cd209a-Ceacam1 v2 = 0.252)\n", P))

cat("\n===== [D] 趋势匹配 =====\n")
cat(sprintf("  L (Cx3cl1) vagus IAV vs Mock log2FC=%.3f, FDR=%.3g -> %s\n",
            L_log2FC_IAV_Mock, L_FDR_IAV_Mock, L_dir_str))
cat(sprintf("  R (Cx3cr1) Mo TX91->PR8_100 slope=%.3f, FDR=%.3g -> %s\n",
            strict_cx3cr1$trend_slope[1], strict_cx3cr1$trend_FDR[1], R_dir_str))
cat(sprintf("  模式: %s\n", mode))

cat("\n===== [E] Cx3cr1 下游通路在 Mo 严选 58 中的命中 =====\n")
for (nm in names(hits)) {
  if (length(hits[[nm]]) > 0) {
    cat(sprintf("  %-18s n=%d -> %s\n", nm, length(hits[[nm]]),
                paste(hits[[nm]], collapse = ", ")))
  } else {
    cat(sprintf("  %-18s n=0\n", nm))
  }
}

cat("\n详细 hit 表 (含趋势):\n")
print(hit_tbl, row.names = FALSE)
sink()
cat("\n", paste(readLines(file.path(OUT, "summary.txt")), collapse="\n"), sep="")

# ============================================================================
# 8. 可视化
# ============================================================================
# 8a. Cx3cl1 vagal sample-level boxplot
df_box <- A_summary
df_box$group <- factor(df_box$group, levels = c("Mock","IAV"))
p1 <- ggplot(df_box, aes(group, CPM, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.85) +
  scale_fill_manual(values = c(Mock = "#377EB8", IAV = "#E41A1C")) +
  labs(x = NULL, y = "CPM",
       title = "Cx3cl1 在 GSE161878 迷走神经的表达",
       subtitle = sprintf("Mock mean CPM=%.2f, IAV mean CPM=%.2f, %d/%d 样本 CPM>%g",
                          L_mock_cpm, L_iav_cpm, n_above, n_samples,
                          PARAMS$cpm_thresh)) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "Cx3cl1_vagal_boxplot.pdf"), p1, width = 5, height = 4,
       device = cairo_pdf)
ggsave(file.path(OUT, "Cx3cl1_vagal_boxplot.png"), p1, width = 5, height = 4,
       dpi = 200, type = "cairo")

# 8b. Cx3cr1 Mo 5 组折线
mo_long_data <- list()
for (gp in mo_groups) {
  idx <- which(mo_m$treat == gp)
  for (i in idx) {
    mo_long_data[[length(mo_long_data) + 1L]] <- data.frame(
      group = gp, sample = colnames(mo_E)[i],
      expr = mo_E[g_R, i], stringsAsFactors = FALSE)
  }
}
mo_long <- do.call(rbind, mo_long_data)
mo_long$group <- factor(mo_long$group, levels = mo_groups)
agg_R <- aggregate(expr ~ group, data = mo_long,
                   FUN = function(v) c(mean = mean(v),
                                       sem  = sd(v)/sqrt(length(v))))
agg_R <- data.frame(group = agg_R$group,
                    mean = agg_R$expr[, "mean"],
                    sem  = agg_R$expr[, "sem"])
agg_R$group <- factor(agg_R$group, levels = mo_groups)

p2 <- ggplot(agg_R, aes(group, mean, group = 1)) +
  geom_line(color = "#D62728", linewidth = 0.8) +
  geom_point(size = 2.5, color = "#D62728") +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                width = 0.18, color = "grey25") +
  labs(x = NULL, y = "log2 expr (neqc)",
       title = "Cx3cr1 在 Mo 细胞 5 组剂量响应",
       subtitle = sprintf("trend_slope=%.3f, FDR=%.3g (Mo 严选 strict 主路径)",
                          strict_cx3cr1$trend_slope[1], strict_cx3cr1$trend_FDR[1])) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "Cx3cr1_Mo_5group_line.pdf"), p2, width = 5, height = 4,
       device = cairo_pdf)
ggsave(file.path(OUT, "Cx3cr1_Mo_5group_line.png"), p2, width = 5, height = 4,
       dpi = 200, type = "cairo")

# 8c. 通路命中条形图
if (nrow(hit_tbl) > 0) {
  bar_df <- as.data.frame(table(hit_tbl$pathway))
  colnames(bar_df) <- c("pathway","n_genes")
  bar_df$pathway <- factor(bar_df$pathway,
                            levels = bar_df$pathway[order(-bar_df$n_genes)])
  p3 <- ggplot(bar_df, aes(pathway, n_genes, fill = pathway)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = n_genes), vjust = -0.3) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Mo 严选 58 命中数",
         title = "Cx3cr1 下游通路在 Mo 严选 58 中的命中数") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, "downstream_pathway_hits_bar.pdf"), p3,
         width = 5, height = 4, device = cairo_pdf)
  ggsave(file.path(OUT, "downstream_pathway_hits_bar.png"), p3,
         width = 5, height = 4, dpi = 200, type = "cairo")
}

# 8d. 综合: 反向 LR 简明示意
df_pair <- data.frame(
  side = c("L (vagus)","R (Mo TX91)"),
  gene = c(g_L, g_R),
  log2_expr = c(L_iav, R_TX91),
  rn = c(L_rn, R_rn),
  hill = c(hill(L_rn), hill(R_rn))
)
p4 <- ggplot(df_pair, aes(side, hill, fill = side)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%s\nrn=%.2f\nHill=%.2f",
                                gene, rn, hill)),
            vjust = -0.2, size = 3.5) +
  scale_fill_manual(values = c("L (vagus)" = "#1F77B4",
                                "R (Mo TX91)" = "#D62728")) +
  ylim(0, 1.1) +
  labs(x = NULL, y = "Hill score",
       title = sprintf("反向 LR: vagus %s -> Mo %s, P_LR = %.3f",
                       g_L, g_R, P)) +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(OUT, "reverse_LR_overview.pdf"), p4, width = 6, height = 4,
       device = cairo_pdf)
ggsave(file.path(OUT, "reverse_LR_overview.png"), p4, width = 6, height = 4,
       dpi = 200, type = "cairo")

# 9. 保存
saveRDS(list(
  Cx3cl1_per_sample = A_summary, Cx3cr1_Mo_5groups = B_summary,
  L_iav_log2 = L_iav, L_mock_log2 = L_mock,
  L_log2FC_IAV_Mock = L_log2FC_IAV_Mock, L_FDR_IAV_Mock = L_FDR_IAV_Mock,
  R_TX91 = R_TX91, R_rn = R_rn, L_rn = L_rn, prob = P,
  trend_match_mode = mode,
  hit_tbl = hit_tbl, group_genes = group_genes,
  PARAMS = PARAMS),
  file.path(OUT, "reverse_LR_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
