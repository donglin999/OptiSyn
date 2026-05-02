#!/usr/bin/env Rscript
# 用同事方法 (detection.tsv + log2FC + t-test + |log2FC|>=1 + padj<0.05)
# 在 5 个 cell types 上跑全 4 contrasts, 看 "Ly 是不是同事看到的突出细胞"
#
# 同事脚本 main_analysis() 的 4 contrast:
#   TX91_vs_Sham
#   PR8_0.6_vs_Sham
#   PR8_10_vs_Sham
#   PR8_100_vs_Sham

suppressPackageStartupMessages({
  library(ggplot2)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT <- file.path(ROOT, "results/celltype_review/colleague_method_check")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

detection <- read.delim(file.path(ROOT, "processed_data/GSE42639/GSE42639_detection.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
intensity <- read.delim(file.path(ROOT, "processed_data/GSE42639/GSE42639_intensity.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(file.path(ROOT, "processed_data/GSE42639/GSE42639_sample_meta.tsv"),
                   stringsAsFactors = FALSE)
detection$probe_id <- trimws(detection$probe_id)
intensity$probe_id <- trimws(intensity$probe_id)

cells <- c("Nh","Am","Ne","Mo","Ly")
contrasts_def <- list(
  TX91_vs_Sham    = c("TX91",       "Sham"),
  `PR8_0.6_vs_Sham` = c("0.6LD50PR8", "Sham"),
  PR8_10_vs_Sham  = c("10LD50PR8",  "Sham"),
  PR8_100_vs_Sham = c("100LD50PR8", "Sham")
)

# 同事方法: detection.tsv + log2(treat_mean) - log2(control_mean) + Welch t-test + BH
colleague_dea <- function(cell, treat, control, det) {
  ctrl_samp  <- meta$sample_id[meta$cell == cell & meta$treat == control]
  treat_samp <- meta$sample_id[meta$cell == cell & meta$treat == treat]
  if (length(ctrl_samp) < 2 || length(treat_samp) < 2) return(NULL)
  ctrl_M  <- as.matrix(det[, ctrl_samp])
  treat_M <- as.matrix(det[, treat_samp])
  ctrl_mean  <- rowMeans(ctrl_M)
  treat_mean <- rowMeans(treat_M)
  log2fc <- log2(treat_mean + 1e-10) - log2(ctrl_mean + 1e-10)
  pval <- vapply(seq_len(nrow(ctrl_M)), function(i) {
    tryCatch(
      t.test(treat_M[i,], ctrl_M[i,], var.equal = FALSE)$p.value,
      error = function(e) 1.0)
  }, numeric(1))
  padj <- p.adjust(pval, method = "BH")
  sig <- padj < 0.05 & abs(log2fc) >= 1
  list(n_total = sum(sig, na.rm = TRUE),
       n_up   = sum(sig & log2fc > 0,  na.rm = TRUE),
       n_down = sum(sig & log2fc < 0,  na.rm = TRUE),
       log2fc = log2fc, padj = padj)
}

# 我们方法 (limma on neqc-expr): 加载已存好的全细胞 DEA
# (results/celltype_amplitude/celltype_DEA_results.rds 已有 PR8hi_vs_Sham, AnyInf, TX91)
# 但 contrasts 不完全一致, 这里现场跑 limma 一致流程

# 先做一次 neqc 标准化
suppressPackageStartupMessages({
  library(limma)
  library(illuminaMousev2.db)
})
sample_cols <- colnames(intensity)[-1]
int_M <- as.matrix(intensity[, sample_cols]); rownames(int_M) <- detection$probe_id
det_M <- as.matrix(detection[, sample_cols]); rownames(det_M) <- detection$probe_id
storage.mode(int_M) <- "numeric"; storage.mode(det_M) <- "numeric"
elist <- new("EListRaw",
             list(E = int_M, other = list(Detection = det_M),
                  targets = meta,
                  genes = data.frame(ProbeID = detection$probe_id)))
elist_n <- neqc(elist)
expr <- elist_n$E
keep_probe <- rowSums(elist_n$other$Detection < 0.05) >= 2
expr_filt <- expr[keep_probe, ]

ours_dea <- function(cell, treat, control) {
  idx <- which(meta$cell == cell)
  m_sub <- meta[idx, ]; E_sub <- expr_filt[, idx]
  treat_R <- c(Sham="Sham","TX91"="TX91","0.6LD50PR8"="PR8_0p6",
               "10LD50PR8"="PR8_10","100LD50PR8"="PR8_100")
  m_sub$tname <- factor(treat_R[m_sub$treat])
  if (!(treat_R[treat] %in% levels(m_sub$tname))) return(NULL)
  if (!(treat_R[control] %in% levels(m_sub$tname))) return(NULL)
  design <- model.matrix(~0 + m_sub$tname)
  colnames(design) <- levels(m_sub$tname)
  cm <- makeContrasts(contrasts = sprintf("%s - %s",
                       treat_R[treat], treat_R[control]), levels = design)
  fit <- eBayes(contrasts.fit(lmFit(E_sub, design), cm), trend = TRUE)
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  sig <- tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1
  list(n_total = sum(sig, na.rm = TRUE),
       n_up   = sum(sig & tt$logFC > 0,  na.rm = TRUE),
       n_down = sum(sig & tt$logFC < 0,  na.rm = TRUE))
}

# 跑全表
cat(">>> 同事方法 5 cells x 4 contrasts:\n")
coll_grid <- data.frame()
for (c in cells) {
  for (cn in names(contrasts_def)) {
    tc <- contrasts_def[[cn]]
    r <- colleague_dea(c, tc[1], tc[2], detection)
    if (is.null(r)) next
    coll_grid <- rbind(coll_grid, data.frame(
      cell = c, contrast = cn,
      method = "colleague_detection_p",
      n_up = r$n_up, n_down = r$n_down, n_total = r$n_total))
  }
}
print(coll_grid)
write.table(coll_grid, file.path(OUT, "colleague_method_5cells_grid.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n>>> 我们方法 (limma+neqc) 5 cells x 4 contrasts:\n")
ours_grid <- data.frame()
for (c in cells) {
  for (cn in names(contrasts_def)) {
    tc <- contrasts_def[[cn]]
    r <- ours_dea(c, tc[1], tc[2])
    if (is.null(r)) next
    ours_grid <- rbind(ours_grid, data.frame(
      cell = c, contrast = cn,
      method = "ours_neqc_limma",
      n_up = r$n_up, n_down = r$n_down, n_total = r$n_total))
  }
}
print(ours_grid)
write.table(ours_grid, file.path(OUT, "ours_method_5cells_grid.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 合并 + 比较
both <- rbind(coll_grid, ours_grid)
write.table(both, file.path(OUT, "compare_methods_5cells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 排序: 哪个细胞在同事方法下 DEG 总和最多?
cat("\n>>> 同事方法各 cell 总 DEG (4 contrasts 求和):\n")
agg_coll <- aggregate(n_total ~ cell, data = coll_grid, FUN = sum)
agg_coll <- agg_coll[order(-agg_coll$n_total), ]
print(agg_coll)

cat("\n>>> 我们方法各 cell 总 DEG (4 contrasts 求和):\n")
agg_ours <- aggregate(n_total ~ cell, data = ours_grid, FUN = sum)
agg_ours <- agg_ours[order(-agg_ours$n_total), ]
print(agg_ours)

# 可视化: 双方法各 cell DEG 数对比
both$cell <- factor(both$cell, levels = c("Mo","Am","Ne","Nh","Ly"))
both$method <- factor(both$method,
                       levels = c("ours_neqc_limma","colleague_detection_p"),
                       labels = c("我们 (neqc+limma+intensity)",
                                  "同事 (t-test + detection p)"))
p <- ggplot(both, aes(cell, n_total, fill = method)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n_total), position = position_dodge(0.9),
            vjust = -0.3, size = 3) +
  facet_wrap(~ contrast, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("我们 (neqc+limma+intensity)" = "#1F77B4",
                                "同事 (t-test + detection p)" = "#FF7F0E")) +
  labs(x = "Cell type", y = "显著 DEG 数 (|log2FC|>=1, padj<0.05)",
       title = "同事方法 vs 我们方法: 5 cells × 4 contrasts DEG 数对比",
       subtitle = "同事方法 detection-p 数据下哪个细胞最突出?") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "compare_methods_5cells_grid.pdf"), p,
       width = 12, height = 8, device = cairo_pdf)
ggsave(file.path(OUT, "compare_methods_5cells_grid.png"), p,
       width = 12, height = 8, dpi = 200, type = "cairo")
cat(sprintf("\n>>> 输出: %s\n", OUT))
