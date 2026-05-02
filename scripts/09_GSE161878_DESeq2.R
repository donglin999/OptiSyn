#!/usr/bin/env Rscript
# GSE161878 迷走神经感觉神经节 (vagal sensory ganglia) IAV vs Mock 差异分析
#
# 数据:
#   raw_data/00_raw_data/GSE161878_gene_counts.txt    (gene x sample, ENTREZ ID, raw count)
#   raw_data/00_raw_data/GSE161878_sample_metadata.txt
#   24 样本: 12 Mock + 12 IAV (含 AM312/313/314 的 D-suffix variant 各 1)
#
# 方法: DESeq2 (Wald + apeglm LFC shrink)
# 阈值: |log2FC| >= 0.585 (FC >= 1.5) & FDR < 0.05
#
# 备注: AM312/AM312D, AM313/AM313D, AM314/AM314D 来自同一动物 (D = duplicate variant).
#       原 metadata: 11 Mock + 13 IAV (含 3 个 D-suffix 重复).
#       D-suffix = 技术重复 (违反样本独立性), 沿用全部样本会扭曲 DESeq2 方差估计.
#       sanity check 显示: 含 D vs 剔 D 主分析 DEG 仅 32% 一致 -> 主分析切换为剔 D-suffix.
#
#       主分析:   11 Mock + 10 IAV (剔 AM312D/313D/314D)  -> GSE161878_DEG_*.tsv
#       敏感性: 11 Mock + 13 IAV (含 D-suffix)            -> GSE161878_DEG_*_withD.tsv
#
# 输出: results/GSE161878/

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(RColorBrewer)
  library(org.Mm.eg.db)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
OUT  <- file.path(ROOT, "results/GSE161878")
PROC <- file.path(ROOT, "processed_data/GSE161878")
dir.create(OUT,  recursive = TRUE, showWarnings = FALSE)
dir.create(PROC, recursive = TRUE, showWarnings = FALSE)

THR <- list(log2fc_abs = 0.585, fdr = 0.05, low_count_min_per_sample = 10,
            low_count_min_n_samples = 3)

# 1. 读取数据 ----------------------------------------------------------------
cnt_path  <- file.path(ROOT, "raw_data/00_raw_data/GSE161878_gene_counts.txt")
meta_path <- file.path(ROOT, "raw_data/00_raw_data/GSE161878_sample_metadata.txt")

cnt <- read.delim(cnt_path, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(meta_path, stringsAsFactors = FALSE)
stopifnot(colnames(cnt)[1] == "GeneID")
genes <- as.character(cnt$GeneID)
M <- as.matrix(cnt[, -1, drop = FALSE]); rownames(M) <- genes
storage.mode(M) <- "integer"

stopifnot(setequal(colnames(M), meta$sample_id))
meta <- meta[match(colnames(M), meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id

cat(sprintf(">>> raw counts: %d genes x %d samples\n", nrow(M), ncol(M)))
cat(">>> 样本分组分布:\n"); print(table(meta$group))

# 2. ENTREZ -> SYMBOL --------------------------------------------------------
ann <- AnnotationDbi::select(org.Mm.eg.db,
                             keys = genes, keytype = "ENTREZID",
                             columns = c("ENTREZID","SYMBOL","GENENAME"))
ann <- ann[!duplicated(ann$ENTREZID), ]
ann <- ann[match(genes, ann$ENTREZID), ]
write.table(ann, file.path(PROC, "gene_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_no_sym <- sum(is.na(ann$SYMBOL))
cat(sprintf(">>> ENTREZ -> SYMBOL: 缺失 %d / %d (后续保留 ENTREZ 作为标识)\n",
            n_no_sym, length(genes)))

# 3. DESeq2 主分析 (12 vs 12) ------------------------------------------------
run_deseq <- function(M, meta, label) {
  meta$group <- factor(meta$group, levels = c("Mock","IAV"))
  dds <- DESeqDataSetFromMatrix(countData = M, colData = meta, design = ~ group)
  # 低表达过滤
  keep <- rowSums(counts(dds) >= THR$low_count_min_per_sample) >=
          THR$low_count_min_n_samples
  cat(sprintf(">>> [%s] 低表达过滤: 保留 %d / %d (>=%d count 在 >= %d 样本)\n",
              label, sum(keep), nrow(dds),
              THR$low_count_min_per_sample, THR$low_count_min_n_samples))
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)

  # apeglm LFC shrink (推荐, 用于排序与火山图)
  res_shrink <- tryCatch(
    lfcShrink(dds, coef = "group_IAV_vs_Mock", type = "apeglm"),
    error = function(e) {
      cat(sprintf("    [%s] apeglm 失败 (%s), 改用 normal\n", label, conditionMessage(e)))
      lfcShrink(dds, coef = "group_IAV_vs_Mock", type = "normal")
    }
  )
  res_raw <- results(dds, name = "group_IAV_vs_Mock")
  list(dds = dds, res_raw = res_raw, res_shrink = res_shrink)
}

# 主分析: 剔 D-suffix (技术重复, 违反独立性)
keep_main <- !grepl("D$", meta$sample_id)
M_main    <- M[, keep_main, drop = FALSE]
meta_main <- meta[keep_main, , drop = FALSE]
n_mock_main <- sum(meta_main$group == "Mock")
n_iav_main  <- sum(meta_main$group == "IAV")
cat(sprintf("\n>>> ===== 主分析 (%d Mock vs %d IAV, 剔 D-suffix) =====\n",
            n_mock_main, n_iav_main))
main <- run_deseq(M_main, meta_main, "main")

# 4. 输出表 ------------------------------------------------------------------
make_full_tbl <- function(res_shrink, res_raw, ann) {
  d <- as.data.frame(res_shrink)
  d$ENTREZID <- rownames(d)
  d$SYMBOL   <- ann$SYMBOL  [match(d$ENTREZID, ann$ENTREZID)]
  d$GENENAME <- ann$GENENAME[match(d$ENTREZID, ann$ENTREZID)]
  # 也带原始 (未 shrink) logFC 以便对照
  d$log2FC_raw <- res_raw$log2FoldChange[match(d$ENTREZID, rownames(res_raw))]
  d$pvalue_raw <- res_raw$pvalue        [match(d$ENTREZID, rownames(res_raw))]
  d <- d[, c("ENTREZID","SYMBOL","GENENAME","baseMean",
             "log2FoldChange","lfcSE","pvalue","padj",
             "log2FC_raw","pvalue_raw")]
  colnames(d)[colnames(d) == "log2FoldChange"] <- "log2FC_shrink"
  colnames(d)[colnames(d) == "padj"] <- "FDR"
  d <- d[order(d$FDR, -abs(d$log2FC_shrink)), ]
  rownames(d) <- NULL
  d
}
full <- make_full_tbl(main$res_shrink, main$res_raw, ann)
write.table(full, file.path(OUT, "GSE161878_DEA_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(full, file.path(OUT, "GSE161878_DEA_full.csv"),
          quote = TRUE, row.names = FALSE)

# 阈值过滤 (用 shrink LFC)
sig <- !is.na(full$FDR) & full$FDR < THR$fdr & abs(full$log2FC_shrink) >= THR$log2fc_abs
full$sig_status <- "ns"
full$sig_status[sig & full$log2FC_shrink >  0] <- "up"
full$sig_status[sig & full$log2FC_shrink <  0] <- "down"
up_tbl   <- full[full$sig_status == "up",   ]
down_tbl <- full[full$sig_status == "down", ]
write.table(up_tbl,   file.path(OUT, "GSE161878_DEG_up.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(down_tbl, file.path(OUT, "GSE161878_DEG_down.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(up_tbl, down_tbl),
            file.path(OUT, "GSE161878_DEG_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n>>> 显著差异基因 (|log2FC|>=%.3f & FDR<%.2f):\n",
            THR$log2fc_abs, THR$fdr))
cat(sprintf("    上调: %d   下调: %d   合计: %d\n",
            nrow(up_tbl), nrow(down_tbl), nrow(up_tbl) + nrow(down_tbl)))

# 5. 敏感性分析: 含 D-suffix 全部样本 (11 Mock + 13 IAV, 原 metadata) --------
n_mock_s <- sum(meta$group == "Mock")
n_iav_s  <- sum(meta$group == "IAV")
cat(sprintf("\n>>> ===== 敏感性 (%d Mock vs %d IAV, 含 D-suffix 重复) =====\n",
            n_mock_s, n_iav_s))
sens <- run_deseq(M, meta, "withD")
full_s <- make_full_tbl(sens$res_shrink, sens$res_raw, ann)
sig_s <- !is.na(full_s$FDR) & full_s$FDR < THR$fdr &
         abs(full_s$log2FC_shrink) >= THR$log2fc_abs
n_up_s   <- sum(sig_s & full_s$log2FC_shrink >  0, na.rm = TRUE)
n_down_s <- sum(sig_s & full_s$log2FC_shrink <  0, na.rm = TRUE)
# 敏感性产物 (附属)
write.table(full_s, file.path(OUT, "GSE161878_DEA_full_withD.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
up_s   <- full_s[sig_s & full_s$log2FC_shrink >  0, ]; up_s$sig_status <- "up"
down_s <- full_s[sig_s & full_s$log2FC_shrink <  0, ]; down_s$sig_status <- "down"
write.table(up_s,   file.path(OUT, "GSE161878_DEG_up_withD.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(down_s, file.path(OUT, "GSE161878_DEG_down_withD.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("    上调: %d   下调: %d   合计: %d\n",
            n_up_s, n_down_s, n_up_s + n_down_s))
ov <- intersect(c(up_tbl$ENTREZID, down_tbl$ENTREZID),
                full_s$ENTREZID[sig_s])
cat(sprintf("    主 (剔 D) ∩ 含 D 敏感性: %d (主 DEG 总数 %d)\n",
            length(ov), nrow(up_tbl) + nrow(down_tbl)))

# 6. 流程日志 ----------------------------------------------------------------
sink(file.path(OUT, "DEA_summary.txt"))
cat("=== GSE161878 迷走神经感觉神经节 IAV vs Mock 差异分析 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/09_GSE161878_DESeq2.R\n")
cat("方法: DESeq2 (Wald) + apeglm LFC shrinkage\n")
cat(sprintf("阈值: |log2FC| >= %.3f & FDR < %.2f\n", THR$log2fc_abs, THR$fdr))
cat(sprintf("低表达过滤: count >= %d 在 >= %d 样本\n",
            THR$low_count_min_per_sample, THR$low_count_min_n_samples))
cat(sprintf("\n样本: 原 %d (%d Mock + %d IAV); 主分析剔 3 个 D-suffix\n", nrow(meta),
            sum(meta$group == "Mock"), sum(meta$group == "IAV")))

cat(sprintf("\n%-55s %-8s %-8s %s\n", "分析", "上调", "下调", "总计"))
cat(sprintf("%-55s %-8d %-8d %d\n",
            sprintf("主分析 (%d Mock vs %d IAV, 剔 D-suffix)",
                    n_mock_main, n_iav_main),
            nrow(up_tbl), nrow(down_tbl), nrow(up_tbl) + nrow(down_tbl)))
cat(sprintf("%-55s %-8d %-8d %d\n",
            sprintf("敏感性 (%d Mock vs %d IAV, 含 D-suffix)",
                    n_mock_s, n_iav_s),
            n_up_s, n_down_s, n_up_s + n_down_s))
cat(sprintf("\n主分析 ∩ 敏感性分析 DEG 一致性: %d / %d\n",
            length(ov), nrow(up_tbl) + nrow(down_tbl)))

cat("\n=== 上调 Top 30 ===\n")
print(head(up_tbl[, c("SYMBOL","ENTREZID","log2FC_shrink","FDR","baseMean")], 30),
      row.names = FALSE)
cat("\n=== 下调 Top 30 ===\n")
print(head(down_tbl[, c("SYMBOL","ENTREZID","log2FC_shrink","FDR","baseMean")], 30),
      row.names = FALSE)
sink()
cat("\n", paste(readLines(file.path(OUT, "DEA_summary.txt")), collapse="\n"), sep="")

# 7. 可视化 ------------------------------------------------------------------
# 7a. PCA
vsd <- vst(main$dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pca_data$sample_id <- rownames(pca_data)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = sample_id)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, vjust = -0.8, color = "black") +
  scale_color_manual(values = c(Mock = "#377EB8", IAV = "#E41A1C")) +
  labs(x = sprintf("PC1: %d%% variance", percentVar[1]),
       y = sprintf("PC2: %d%% variance", percentVar[2]),
       title = "GSE161878 vagal ganglia PCA (vst)",
       subtitle = sprintf("主分析: %d Mock + %d IAV (剔 3 个 D-suffix 重复)",
                          n_mock_main, n_iav_main)) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "PCA.pdf"), p_pca, width = 7, height = 5,
       device = cairo_pdf)
ggsave(file.path(OUT, "PCA.png"), p_pca, width = 7, height = 5,
       dpi = 200, type = "cairo")
cat(">>> PCA: PCA.pdf/png\n")

# 7b. 火山图
make_volcano <- function(d, out_prefix, title_main) {
  d$negLogFDR <- -log10(d$FDR)
  d$negLogFDR[is.infinite(d$negLogFDR)] <- max(d$negLogFDR[is.finite(d$negLogFDR)], na.rm=TRUE)
  d$status <- "ns"
  d$status[d$sig_status == "up"]   <- "up"
  d$status[d$sig_status == "down"] <- "down"
  top_label <- rbind(head(d[d$sig_status == "up",   ], 12),
                     head(d[d$sig_status == "down", ], 12))
  p <- ggplot(d, aes(log2FC_shrink, negLogFDR)) +
    geom_point(aes(color = status), size = 0.7, alpha = 0.6) +
    scale_color_manual(values = c(ns = "grey80", up = "#E41A1C", down = "#377EB8")) +
    geom_vline(xintercept = c(-THR$log2fc_abs, THR$log2fc_abs),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(THR$fdr),
               linetype = "dashed", color = "grey50") +
    geom_text(data = top_label, aes(label = SYMBOL),
              size = 2.6, hjust = -0.1, vjust = -0.2, color = "black") +
    labs(x = "log2FC (IAV vs Mock, shrink)", y = "-log10(FDR)",
         title = title_main,
         subtitle = sprintf("|log2FC|>=%.3f & FDR<%.2f -> 上 %d / 下 %d",
                            THR$log2fc_abs, THR$fdr,
                            sum(d$sig_status == "up"),
                            sum(d$sig_status == "down"))) +
    theme_bw(base_size = 11)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = 6, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 6,
         dpi = 200, type = "cairo")
}
make_volcano(full, file.path(OUT, "volcano_main"),
             "GSE161878 IAV vs Mock 火山图 (主分析)")
cat(">>> 火山图: volcano_main.pdf/png\n")

# 7c. Top 50 DEG 热图 (vst)
top_n <- min(50, nrow(up_tbl) + nrow(down_tbl))
if (top_n > 1) {
  top_genes <- c(head(up_tbl$ENTREZID, ceiling(top_n/2)),
                 head(down_tbl$ENTREZID, floor(top_n/2)))
  top_genes <- top_genes[top_genes %in% rownames(vsd)]
  M_v <- assay(vsd)[top_genes, , drop = FALSE]
  Z <- t(scale(t(M_v))); Z[is.na(Z)] <- 0
  # 行标签用 SYMBOL
  rl <- ann$SYMBOL[match(rownames(Z), ann$ENTREZID)]
  rl[is.na(rl)] <- rownames(Z)[is.na(rl)]
  rownames(Z) <- make.unique(rl)
  # 列序: Mock 先 -> IAV 后 (用主分析 meta = 剔 D-suffix)
  col_ord <- order(meta_main$group, meta_main$sample_id)
  Z <- Z[, col_ord, drop = FALSE]
  m_ord <- meta_main[col_ord, , drop = FALSE]
  hc <- hclust(dist(Z), method = "average"); Z <- Z[hc$order, ]

  df <- as.data.frame(Z); df$gene <- rownames(Z)
  df_long <- reshape(df, varying = colnames(Z), v.names = "z",
                     timevar = "sample", times = colnames(Z), direction = "long")
  df_long$gene   <- factor(df_long$gene, levels = rev(rownames(Z)))
  df_long$sample <- factor(df_long$sample, levels = colnames(Z))
  df_long$group  <- m_ord$group[match(df_long$sample, m_ord$sample_id)]
  cap <- 2.5; df_long$z_cap <- pmax(pmin(df_long$z, cap), -cap)

  p_hm <- ggplot(df_long, aes(sample, gene, fill = z_cap)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-cap, cap),
                         oob = scales::squish, name = "row z") +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL,
         title = sprintf("GSE161878 Top %d DEG 热图 (vst, 行 z-score)", nrow(Z)),
         subtitle = "上下各取 ~半数, 列按 Mock/IAV 分块") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 6),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90", color = NA),
          strip.text = element_text(face = "bold"))
  h_in <- max(6, nrow(Z) * 0.18 + 2)
  ggsave(file.path(OUT, "heatmap_top50_DEG.pdf"), p_hm,
         width = 9, height = h_in, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT, "heatmap_top50_DEG.png"), p_hm,
         width = 9, height = h_in, dpi = 200, type = "cairo", limitsize = FALSE)
  cat(">>> 热图: heatmap_top50_DEG.pdf/png\n")
}

# 8. 参数 + sessionInfo ------------------------------------------------------
sink(file.path(OUT, "params_and_session_info.txt"))
cat(sprintf("时间戳: %s\n", TS))
cat("脚本: scripts/09_GSE161878_DESeq2.R\n")
cat("数据: GSE161878 vagal sensory ganglia, 24 sample, ENTREZ raw count\n")
cat("方法: DESeq2 + apeglm shrink\n\n")
cat("阈值:\n"); print(THR)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

# 9. 保存中间对象 + vst 矩阵 -------------------------------------------------
saveRDS(list(dds = main$dds, res_raw = main$res_raw,
             res_shrink = main$res_shrink,
             full = full, up_tbl = up_tbl, down_tbl = down_tbl,
             vsd = vsd, meta = meta_main, THR = THR,
             dds_withD = sens$dds, full_withD = full_s),
        file.path(OUT, "GSE161878_DEA.rds"))
# vst + meta 保存为主分析版本 (剔 D-suffix), 下游 (CellChat) 默认读这个
saveRDS(assay(vsd), file.path(PROC, "expr_vst.rds"))
saveRDS(meta_main,   file.path(PROC, "sample_meta.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
