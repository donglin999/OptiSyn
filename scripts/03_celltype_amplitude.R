#!/usr/bin/env Rscript
# 复核 "Mo 细胞主导免疫应答" 结论：对 5 种细胞分别做 limma 差异分析,
# 比较应答幅度（DEG 数 / 平均 |log2FC| / 方差解释）。
#
# 每细胞类型独立建模:
#   ~0 + treat,  treat ∈ {Sham, TX91, 0.6LD50PR8, 10LD50PR8, 100LD50PR8}
# 三个对比:
#   C1 = TX91   - Sham        (轻症应答)
#   C2 = 100LD50PR8 - Sham    (重症应答, 最高致死剂量)
#   C3 = (TX91+0.6+10+100)/4 - Sham  (总体感染应答)
#
# 输出:
#   results/celltype_amplitude/celltype_DEA_summary.tsv
#   results/celltype_amplitude/{cell}_{contrast}_DEGs_full.tsv
#   results/celltype_amplitude/celltype_DEA_results.rds

suppressPackageStartupMessages({
  library(limma)
})

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
OUT  <- file.path(ROOT, "results/celltype_amplitude")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

expr <- readRDS(file.path(ROOT, "processed_data/GSE42639/expr_log2.rds"))
meta <- readRDS(file.path(ROOT, "processed_data/GSE42639/sample_meta.rds"))
stopifnot(identical(colnames(expr), meta$sample_id))

# 处理因子: 重命名为 R 合法名 (PR8 dose 加前缀以满足 makeContrasts)
treat_levels_raw <- c("Sham", "TX91", "0.6LD50PR8", "10LD50PR8", "100LD50PR8")
treat_levels_R   <- c("Sham", "TX91", "PR8_0p6LD50", "PR8_10LD50", "PR8_100LD50")
names(treat_levels_R) <- treat_levels_raw
meta$treat <- factor(treat_levels_R[meta$treat], levels = treat_levels_R)

cells <- c("Nh", "Am", "Ne", "Mo", "Ly")
FDR_TH <- 0.05
LFC_TH <- 1

run_one_cell <- function(cell) {
  idx <- which(meta$cell == cell)
  m   <- meta[idx, , drop = FALSE]
  E   <- expr[, idx, drop = FALSE]
  treat <- droplevels(m$treat)
  cat(sprintf("\n--- %s : n=%d ; treat counts:\n", cell, length(idx)))
  print(table(treat))

  # 设计矩阵
  design <- model.matrix(~0 + treat)
  colnames(design) <- gsub("^treat", "", colnames(design))

  # 三个对比 (用 makeContrasts；对没有的 level 用 NA 跳过)
  contr_list <- list()
  if (all(c("TX91", "Sham") %in% colnames(design)))
    contr_list$TX91_vs_Sham   <- "TX91 - Sham"
  if (all(c("PR8_100LD50", "Sham") %in% colnames(design)))
    contr_list$PR8hi_vs_Sham  <- "PR8_100LD50 - Sham"
  inf_set <- intersect(c("TX91","PR8_0p6LD50","PR8_10LD50","PR8_100LD50"),
                       colnames(design))
  if (length(inf_set) >= 2 && "Sham" %in% colnames(design)) {
    contr_list$AnyInf_vs_Sham <- paste0("(",
        paste0(inf_set, collapse = " + "),
        ")/", length(inf_set), " - Sham")
  }

  cm <- makeContrasts(contrasts = unlist(contr_list), levels = design)
  colnames(cm) <- names(contr_list)

  fit  <- lmFit(E, design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2, trend = TRUE)

  per_contrast <- list()
  smry <- list()
  for (k in names(contr_list)) {
    tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    tt <- tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
    fn <- file.path(OUT, sprintf("%s_%s_DEGs_full.tsv", cell, k))
    write.table(tt, fn, sep="\t", quote=FALSE, row.names=FALSE)

    sig <- tt[tt$adj.P.Val < FDR_TH & abs(tt$logFC) > LFC_TH, , drop = FALSE]
    smry[[k]] <- data.frame(
      cell      = cell,
      contrast  = k,
      n_total   = nrow(tt),
      n_sig     = nrow(sig),
      n_up      = sum(sig$logFC > 0),
      n_down    = sum(sig$logFC < 0),
      mean_abs_logFC_sig = if (nrow(sig)) mean(abs(sig$logFC)) else NA_real_,
      median_abs_logFC_top500 = median(abs(tt$logFC[order(tt$adj.P.Val)][1:500])),
      stringsAsFactors = FALSE
    )
    per_contrast[[k]] <- tt
  }
  list(summary = do.call(rbind, smry),
       per_contrast = per_contrast,
       fit = fit2)
}

results <- lapply(setNames(cells, cells), run_one_cell)
summary_tbl <- do.call(rbind, lapply(results, `[[`, "summary"))
rownames(summary_tbl) <- NULL
summary_tbl$cell <- factor(summary_tbl$cell, levels = cells)
summary_tbl <- summary_tbl[order(summary_tbl$contrast, summary_tbl$cell), ]

cat("\n========= Summary: 各细胞 × 对比 应答幅度 =========\n")
print(summary_tbl, row.names = FALSE)

write.table(summary_tbl, file.path(OUT, "celltype_DEA_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(results, file.path(OUT, "celltype_DEA_results.rds"))
cat(sprintf("\n>>> 输出: %s\n", OUT))

# 显式判断 Mo 是否主导
cat("\n========= Mo 主导性判断 =========\n")
for (k in unique(summary_tbl$contrast)) {
  sub <- summary_tbl[summary_tbl$contrast == k, ]
  ord <- sub[order(-sub$n_sig), ]
  cat(sprintf("  %s : DEG数排序 = %s ; Mo排第 %d (n_sig=%d)\n",
              k, paste(ord$cell, collapse=" > "),
              which(ord$cell == "Mo"),
              sub$n_sig[sub$cell == "Mo"]))
}
