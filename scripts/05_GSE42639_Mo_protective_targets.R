#!/usr/bin/env Rscript
# GSE42639 Mo 细胞 -> 保护性靶点筛选 (轻症维持 -> 重症崩塌)
#
# 数据: processed_data/GSE42639/expr_log2.rds (limma::neqc, gene x sample)
#       processed_data/GSE42639/sample_meta.rds
# 仅 cell == "Mo", 每组 n=3:
#   Sham / TX91 / PR8_0p6 / PR8_10 / PR8_100
#
# 5 步顺序流程 (不可交换):
#   Step 1 趋势 (TX91=0,0p6=1,10=2,100=3): slope<0 & FDR<0.10
#   Step 2 崩塌 (PR8_100 vs TX91): logFC<=-0.585 & FDR<0.05
#   Step 3 基础表达: TX91 mean >= 全基因 TX91 mean 的中位
#   Step 4 反向假阳性剔除 (PR8_100 vs TX91 logFC>=1 & FDR<0.05)
#   Step 5 排除轻症已下调 (TX91 vs Sham): 留 FDR>=0.05 OR logFC>=0
#
# 输出: results/Mo_protective/

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
OUT  <- file.path(ROOT, "results/Mo_protective")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 阈值集中配置 (写入 params 记录文件)
THR <- list(
  step1_slope_max     = 0,        # slope < 0
  step1_fdr_max       = 0.10,
  step2_log2fc_max    = -0.585,   # FC <= 0.67
  step2_fdr_max       = 0.05,
  step3_quantile      = 0.5,      # >= 全基因 TX91 mean 中位
  step4_log2fc_min    = 1,        # 反向: 剔除 >= 1
  step4_fdr_max       = 0.05,
  step5_fdr_max       = 0.05,     # TX91 vs Sham 显著阈
  step5_log2fc_max    = 0
)

GROUP_LEVELS <- c("Sham","TX91","PR8_0p6","PR8_10","PR8_100")
SEV_MAP <- c(TX91=0, PR8_0p6=1, PR8_10=2, PR8_100=3)

# 1. 数据加载与子集 -----------------------------------------------------------
expr_all <- readRDS(file.path(ROOT, "processed_data/GSE42639/expr_log2.rds"))
meta_all <- readRDS(file.path(ROOT, "processed_data/GSE42639/sample_meta.rds"))

mo_idx <- which(meta_all$cell == "Mo")
m <- meta_all[mo_idx, , drop = FALSE]
E <- expr_all[, mo_idx, drop = FALSE]

treat_R <- c(Sham="Sham", TX91="TX91",
             "0.6LD50PR8"="PR8_0p6",
             "10LD50PR8" ="PR8_10",
             "100LD50PR8"="PR8_100")
m$treat <- factor(treat_R[m$treat], levels = GROUP_LEVELS)

# 强制每组 n=3
tab <- table(m$treat)
cat(">>> Mo 样本分布:\n"); print(tab)
if (!all(tab == 3)) {
  stop(sprintf("分组样本量不满足 n=3: %s",
               paste(names(tab), tab, sep="=", collapse=", ")))
}

cat(sprintf(">>> Mo 表达矩阵: %d 基因 x %d 样本\n", nrow(E), ncol(E)))

# 2. 主线性模型 (~0+treat) ----------------------------------------------------
designAll <- model.matrix(~0 + m$treat)
colnames(designAll) <- levels(m$treat)
fitAll <- lmFit(E, designAll)
cm <- makeContrasts(
  TX91_vs_Sham    = TX91    - Sham,
  PR8_100_vs_TX91 = PR8_100 - TX91,
  PR8_0p6_vs_TX91 = PR8_0p6 - TX91,
  PR8_10_vs_TX91  = PR8_10  - TX91,
  levels = designAll
)
fit2 <- eBayes(contrasts.fit(fitAll, cm), trend = TRUE)

get_tt <- function(coef) {
  t <- topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  t$gene <- rownames(t); rownames(t) <- NULL
  t[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
}
ttSham <- get_tt("TX91_vs_Sham")
tt100  <- get_tt("PR8_100_vs_TX91")
tt06   <- get_tt("PR8_0p6_vs_TX91")
tt10   <- get_tt("PR8_10_vs_TX91")

# 3. 趋势模型 (仅感染组) ------------------------------------------------------
inf_idx <- which(m$treat %in% names(SEV_MAP))
m_inf <- m[inf_idx, , drop = FALSE]
E_inf <- E[, inf_idx, drop = FALSE]
m_inf$severity <- SEV_MAP[as.character(m_inf$treat)]
designSev <- model.matrix(~ m_inf$severity)
colnames(designSev) <- c("Intercept","severity")
fitSev <- eBayes(lmFit(E_inf, designSev), trend = TRUE)
ttSev <- topTable(fitSev, coef = "severity", number = Inf, sort.by = "none")
ttSev$gene <- rownames(ttSev); rownames(ttSev) <- NULL
colnames(ttSev)[colnames(ttSev) == "logFC"] <- "slope"
ttSev <- ttSev[, c("gene","slope","AveExpr","t","P.Value","adj.P.Val","B")]

# 4. TX91 平均表达 ------------------------------------------------------------
tx91_cols <- which(m$treat == "TX91")
tx91_mean <- rowMeans(E[, tx91_cols, drop = FALSE])
names(tx91_mean) <- rownames(E)
median_thr <- as.numeric(quantile(tx91_mean, THR$step3_quantile, na.rm = TRUE))

# 5. 顺序筛选 -----------------------------------------------------------------
n_total <- nrow(E)

# Step 1: 趋势
keep1 <- ttSev$slope < THR$step1_slope_max & ttSev$adj.P.Val < THR$step1_fdr_max
S1 <- ttSev$gene[keep1]
write.table(ttSev[order(ttSev$adj.P.Val), ],
            file.path(OUT, "step1_trend_full.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(ttSev[ttSev$gene %in% S1, ][
              order(ttSev$adj.P.Val[match(S1, ttSev$gene)]), ],
            file.path(OUT, "step1_trend_pass.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf(">>> Step1 趋势 slope<0 & FDR<%.2f: %d / %d\n",
            THR$step1_fdr_max, length(S1), n_total))
if (length(S1) == 0) stop("Step 1 无基因通过 -> 数据或趋势模型异常")

# Step 2: 崩塌
g_collapse <- tt100$gene[tt100$logFC <= THR$step2_log2fc_max &
                         tt100$adj.P.Val < THR$step2_fdr_max]
S2 <- intersect(S1, g_collapse)
write.table(tt100[order(tt100$adj.P.Val), ],
            file.path(OUT, "step2_PR8_100_vs_TX91_full.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(tt100[tt100$gene %in% S2, ],
            file.path(OUT, "step2_pass.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf(">>> Step2 崩塌 logFC<=%.3f & FDR<%.2f (与 Step1 交): %d\n",
            THR$step2_log2fc_max, THR$step2_fdr_max, length(S2)))
if (length(S2) == 0) stop("Step 2 无基因通过 -> 崩塌门槛过严或数据无 100LD50 显著下调")

# Step 3: 基础表达
keep3 <- tx91_mean[S2] >= median_thr
S3 <- S2[keep3]
step3_tbl <- data.frame(
  gene = S2,
  TX91_mean = tx91_mean[S2],
  TX91_mean_threshold = median_thr,
  pass = keep3
)
write.table(step3_tbl[order(-step3_tbl$TX91_mean), ],
            file.path(OUT, "step3_TX91_baseline.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf(">>> Step3 TX91 mean >= %.3f (中位): %d\n",
            median_thr, length(S3)))
if (length(S3) == 0) stop("Step 3 后无基因通过 -> 表达量整体偏低")

# Step 4: 反向假阳性剔除 (PR8_100 vs TX91 显著上调)
g_revUp <- tt100$gene[tt100$logFC >= THR$step4_log2fc_min &
                      tt100$adj.P.Val < THR$step4_fdr_max]
S4 <- setdiff(S3, g_revUp)
n_removed_step4 <- length(S3) - length(S4)
step4_removed <- intersect(S3, g_revUp)
write.table(data.frame(removed_gene = step4_removed,
                       log2FC_100_vs_TX91 = tt100$logFC[match(step4_removed, tt100$gene)],
                       FDR_100_vs_TX91   = tt100$adj.P.Val[match(step4_removed, tt100$gene)]),
            file.path(OUT, "step4_reverse_removed.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf(">>> Step4 反向剔除 logFC>=%.0f & FDR<%.2f: 剔除 %d, 剩 %d\n",
            THR$step4_log2fc_min, THR$step4_fdr_max,
            n_removed_step4, length(S4)))
if (length(S4) == 0) stop("Step 4 后无基因留存")

# Step 5: 排除轻症已下调
sham_lfc <- ttSham$logFC[match(S4, ttSham$gene)]
sham_fdr <- ttSham$adj.P.Val[match(S4, ttSham$gene)]
already_down <- (sham_fdr < THR$step5_fdr_max) & (sham_lfc < THR$step5_log2fc_max)
already_down[is.na(already_down)] <- FALSE
pass5 <- !already_down
S5 <- S4[pass5]
step5_tbl <- data.frame(
  gene = S4,
  log2FC_TX91_vs_Sham = sham_lfc,
  FDR_TX91_vs_Sham   = sham_fdr,
  already_down_in_TX91 = already_down,
  pass_step5 = pass5
)
write.table(step5_tbl[order(step5_tbl$FDR_TX91_vs_Sham), ],
            file.path(OUT, "step5_TX91_vs_Sham_check.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf(">>> Step5 排除轻症已下调: %d (剔 %d)\n",
            length(S5), sum(already_down)))
if (length(S5) == 0) stop("Step 5 后无基因留存 -> 候选可能全部为轻症已下调")

# 6. core_targets ------------------------------------------------------------
g <- S5
core <- data.frame(
  gene_symbol            = g,
  trend_slope            = ttSev$slope[match(g, ttSev$gene)],
  trend_FDR              = ttSev$adj.P.Val[match(g, ttSev$gene)],
  log2FC_100LD50_vs_TX91 = tt100$logFC[match(g, tt100$gene)],
  FDR_100LD50_vs_TX91    = tt100$adj.P.Val[match(g, tt100$gene)],
  log2FC_TX91_vs_Sham    = ttSham$logFC[match(g, ttSham$gene)],
  FDR_TX91_vs_Sham       = ttSham$adj.P.Val[match(g, ttSham$gene)],
  meanExpr_TX91          = as.numeric(tx91_mean[g]),
  pass_step5             = TRUE,
  notes                  = "",
  stringsAsFactors       = FALSE
)
# 排序: slope 升 (越负越前) -> FDR_100 升 -> meanExpr_TX91 降
core <- core[order(core$trend_slope,
                   core$FDR_100LD50_vs_TX91,
                   -core$meanExpr_TX91), ]
rownames(core) <- NULL
write.table(core, file.path(OUT, "core_targets.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.csv(core, file.path(OUT, "core_targets.csv"),
          quote = TRUE, row.names = FALSE)
cat(sprintf(">>> core_targets: %d 基因, 已写出\n", nrow(core)))

# 7. 流程追踪与参数记录 -------------------------------------------------------
sink(file.path(OUT, "flow_summary.txt"))
cat("=== Mo 保护性靶点筛选 (轻症维持 -> 重症崩塌) ===\n")
cat(sprintf("时间: %s\n", TS))
cat(sprintf("数据: GSE42639 Mo 细胞, 5 组 x n=3, neqc 标准化\n\n"))
cat(sprintf("%-60s %s\n", "节点", "基因数"))
cat(sprintf("%-60s %s\n", "----", "------"))
cat(sprintf("%-60s %d\n", "起始 (探针注释合并后)", n_total))
cat(sprintf("%-60s %d\n",
            sprintf("Step1 趋势 slope<0 & FDR<%.2f", THR$step1_fdr_max), length(S1)))
cat(sprintf("%-60s %d\n",
            sprintf("Step2 崩塌 logFC<=%.3f & FDR<%.2f (∩Step1)",
                    THR$step2_log2fc_max, THR$step2_fdr_max), length(S2)))
cat(sprintf("%-60s %d\n",
            sprintf("Step3 TX91 mean >= %.3f (q=%.2f)",
                    median_thr, THR$step3_quantile), length(S3)))
cat(sprintf("%-60s %d\n",
            sprintf("Step4 剔反向上调 logFC>=%.0f & FDR<%.2f",
                    THR$step4_log2fc_min, THR$step4_fdr_max), length(S4)))
cat(sprintf("%-60s %d\n",
            sprintf("Step5 排除轻症已下调 (TX91 vs Sham)"), length(S5)))
cat("\n=== Top 20 (按排序规则) ===\n")
print(head(core, 20), row.names = FALSE)
sink()
cat(readLines(file.path(OUT, "flow_summary.txt")), sep = "\n")

# params + sessionInfo
sink(file.path(OUT, "params_and_session_info.txt"))
cat("=== 参数记录 ===\n")
cat(sprintf("时间戳: %s\n", TS))
cat("脚本: scripts/05_GSE42639_Mo_protective_targets.R\n")
cat(sprintf("随机种子: 42\n"))
cat("\n分组等级 (列序固定):\n")
print(GROUP_LEVELS)
cat("\n严重度编码 (Step1 趋势):\n")
print(SEV_MAP)
cat("\n阈值:\n")
print(THR)
cat(sprintf("\nTX91 mean 阈值 (q=%.2f) = %.6f\n",
            THR$step3_quantile, median_thr))
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

# 8. 中间对象 -----------------------------------------------------------------
saveRDS(list(
  ttSham = ttSham, tt100 = tt100, tt06 = tt06, tt10 = tt10, ttSev = ttSev,
  S1 = S1, S2 = S2, S3 = S3, S4 = S4, S5 = S5,
  core = core, meta = m, expr = E, tx91_mean = tx91_mean,
  THR = THR, GROUP_LEVELS = GROUP_LEVELS, SEV_MAP = SEV_MAP
), file.path(OUT, "Mo_protective_results.rds"))

# 9. 图 1: 剂量响应热图 -------------------------------------------------------
make_heatmap <- function(core, E, m, out_prefix) {
  if (nrow(core) == 0) return(invisible(NULL))
  g <- core$gene_symbol
  M <- E[g, , drop = FALSE]

  # 列顺序: 按分组固定, 组内按 sample_id 排序
  col_ord <- unlist(lapply(GROUP_LEVELS, function(lv) {
    idx <- which(m$treat == lv)
    idx[order(m$sample_id[idx])]
  }))
  M <- M[, col_ord, drop = FALSE]
  m_ord <- m[col_ord, , drop = FALSE]

  # 行 z-score
  Z <- t(scale(t(M)))
  Z[is.na(Z)] <- 0

  # 行聚类 (可选)
  if (nrow(Z) > 2) {
    hc <- hclust(dist(Z), method = "average")
    Z <- Z[hc$order, , drop = FALSE]
  }

  # 长格式
  df <- as.data.frame(Z)
  df$gene <- rownames(Z)
  df_long <- reshape(df, varying = colnames(Z), v.names = "z",
                     timevar = "sample", times = colnames(Z),
                     direction = "long")
  df_long$gene   <- factor(df_long$gene, levels = rev(rownames(Z)))
  df_long$sample <- factor(df_long$sample, levels = colnames(Z))
  df_long$group  <- m_ord$treat[match(df_long$sample, m_ord$sample_id)]

  # 颜色封顶以避免极端值压缩
  cap <- 2.5
  df_long$z_cap <- pmax(pmin(df_long$z, cap), -cap)

  p <- ggplot(df_long, aes(sample, gene, fill = z_cap)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0,
                         limits = c(-cap, cap), oob = scales::squish,
                         name = "row z-score") +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL,
         title = sprintf("Mo 保护性候选 (n=%d) 剂量响应热图", nrow(Z)),
         subtitle = "行标准化 (z-score), 列序固定, 列不聚类") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90", color = NA),
          strip.text = element_text(face = "bold"))

  h_in <- max(3, min(20, nrow(Z) * 0.12 + 2))
  ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = h_in,
         limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 8, height = h_in,
         dpi = 200, limitsize = FALSE, type = "cairo")
  invisible(p)
}
suppressWarnings(suppressMessages({
  if (requireNamespace("scales", quietly = TRUE)) {
    make_heatmap(core, E, m, file.path(OUT, "heatmap_dose_response"))
    cat(">>> 热图: heatmap_dose_response.pdf/png\n")
  } else {
    cat("!! scales 包不可用, 跳过热图\n")
  }
}))

# 10. 图 2: Top20 折线图 (组均值 + SEM, n=3) ----------------------------------
make_lineplot <- function(core, E, m, out_prefix, top_n = 20) {
  if (nrow(core) == 0) return(invisible(NULL))
  g_top <- head(core$gene_symbol, top_n)
  M <- E[g_top, , drop = FALSE]

  long_list <- lapply(g_top, function(gn) {
    x <- M[gn, ]
    data.frame(gene = gn,
               sample_id = names(x),
               group = m$treat[match(names(x), m$sample_id)],
               expr = as.numeric(x),
               stringsAsFactors = FALSE)
  })
  long_df <- do.call(rbind, long_list)

  agg <- aggregate(expr ~ gene + group, data = long_df,
                   FUN = function(v) c(mean = mean(v),
                                       sem  = sd(v)/sqrt(length(v))))
  agg <- data.frame(gene = agg$gene, group = agg$group,
                    mean = agg$expr[, "mean"],
                    sem  = agg$expr[, "sem"])
  agg$group <- factor(agg$group, levels = GROUP_LEVELS)
  agg$gene  <- factor(agg$gene, levels = g_top)

  p <- ggplot(agg, aes(group, mean, group = gene)) +
    geom_line(color = "#377EB8") +
    geom_point(color = "#377EB8", size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, color = "grey30") +
    facet_wrap(~ gene, scales = "free_y", ncol = 5) +
    labs(x = NULL, y = "log2 expr (neqc)",
         title = sprintf("Top %d 保护性候选 - 剂量响应 (组均值 +/- SEM, n=3)",
                         length(g_top)),
         subtitle = "误差线: SEM (n=3 生物学重复)") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(face = "bold"))

  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = 8, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 8, dpi = 200, type = "cairo")
  invisible(p)
}
make_lineplot(core, E, m, file.path(OUT, "top20_lineplot"))
cat(">>> Top20 折线图: top20_lineplot.pdf/png\n")

# 11. 图 3 (附加): 火山图 PR8_100 vs TX91 -------------------------------------
make_volcano <- function(tt100, core, out_prefix,
                         lfc_thr = THR$step2_log2fc_max,
                         fdr_thr = THR$step2_fdr_max) {
  d <- tt100
  d$negLogFDR <- -log10(d$adj.P.Val)
  d$status <- "ns"
  d$status[d$logFC <=  lfc_thr & d$adj.P.Val < fdr_thr] <- "down"
  d$status[d$logFC >= -lfc_thr & d$adj.P.Val < fdr_thr & d$logFC > 0] <- "up"
  d$is_core <- d$gene %in% core$gene_symbol
  topN <- head(core$gene_symbol, 15)

  p <- ggplot(d, aes(logFC, negLogFDR)) +
    geom_point(aes(color = status), size = 0.6, alpha = 0.6) +
    geom_point(data = d[d$is_core, ], color = "#E41A1C", size = 1.5) +
    scale_color_manual(values = c(ns = "grey80", down = "#377EB8", up = "#FB9A99")) +
    geom_vline(xintercept = c(lfc_thr, -lfc_thr), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", color = "grey50") +
    labs(x = "log2FC (PR8_100LD50 vs TX91)", y = "-log10(FDR)",
         title = "PR8_100LD50 vs TX91 火山图",
         subtitle = sprintf("红点: core_targets (n=%d), 标注 Top15", nrow(core))) +
    theme_bw(base_size = 10)

  if (requireNamespace("ggrepel", quietly = TRUE) && length(topN) > 0) {
    library(ggrepel)
    p <- p + geom_text_repel(data = d[d$gene %in% topN, ],
                             aes(label = gene), size = 2.8,
                             max.overlaps = 30, color = "black")
  } else if (length(topN) > 0) {
    p <- p + geom_text(data = d[d$gene %in% topN, ],
                       aes(label = gene), size = 2.8,
                       hjust = -0.1, vjust = -0.2)
  }
  ggsave(paste0(out_prefix, ".pdf"), p, width = 7, height = 6, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 7, height = 6, dpi = 200, type = "cairo")
  invisible(p)
}
make_volcano(tt100, core, file.path(OUT, "volcano_PR8_100_vs_TX91"))
cat(">>> 火山图: volcano_PR8_100_vs_TX91.pdf/png\n")

# 12. 可选 GO 富集 (limma::goana, org.Mm.eg.db) -------------------------------
try_goana <- function(core, universe_genes, out_prefix) {
  if (nrow(core) == 0) return(invisible(NULL))
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    cat("!! org.Mm.eg.db 不可用, 跳过富集\n"); return(invisible(NULL))
  }
  ok <- try({
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    sym2eg <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys = unique(c(core$gene_symbol, universe_genes)),
                                    keytype = "SYMBOL",
                                    columns = c("ENTREZID","SYMBOL"))
    sym2eg <- sym2eg[!is.na(sym2eg$ENTREZID), ]
    de_eg  <- unique(sym2eg$ENTREZID[match(core$gene_symbol, sym2eg$SYMBOL)])
    de_eg  <- de_eg[!is.na(de_eg)]
    uni_eg <- unique(sym2eg$ENTREZID[sym2eg$SYMBOL %in% universe_genes])
    uni_eg <- uni_eg[!is.na(uni_eg)]
    if (length(de_eg) < 5) {
      cat(sprintf("!! 候选基因映射 ENTREZ 后 < 5 (%d), 跳过富集\n", length(de_eg)))
      return(NULL)
    }
    go <- limma::goana(de = de_eg, universe = uni_eg, species = "Mm")
    go <- go[order(go$P.DE), ]
    go$FDR <- p.adjust(go$P.DE, "BH")
    write.table(go, paste0(out_prefix, "_GO_full.tsv"),
                sep="\t", quote=FALSE, row.names=TRUE)
    top_go <- head(go[go$FDR < 0.25, ], 30)
    write.table(top_go, paste0(out_prefix, "_GO_top.tsv"),
                sep="\t", quote=FALSE, row.names=TRUE)
    cat(sprintf(">>> GO 富集完成 (de=%d, FDR<0.25 数=%d)\n",
                length(de_eg), nrow(top_go)))
  }, silent = TRUE)
  if (inherits(ok, "try-error")) cat("!! GO 富集出错: ", conditionMessage(attr(ok,"condition")), "\n")
  invisible(NULL)
}
try_goana(core, universe_genes = rownames(E),
          out_prefix = file.path(OUT, "enrichment"))

cat(sprintf("\n>>> 全部输出位于: %s\n", OUT))
