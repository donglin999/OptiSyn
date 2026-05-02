#!/usr/bin/env Rscript
# 从 results/intersectA/intersectA_results.rds 中恢复 A2 ∩ A3 = 49 基因.
# 仅做查看, 不改 04 脚本.

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
ia   <- readRDS(file.path(ROOT, "results/intersectA/intersectA_results.rds"))

# 旧 RDS 字段: ttS1 ttS2 tt06 tt10 tt100 A1 A2 A3 core core_tbl
A2A3 <- intersect(ia$A2, ia$A3)
stopifnot(length(A2A3) == 49)

# 顺手挂上 TX91 mean (需要 expr / meta)
expr <- readRDS(file.path(ROOT, "processed_data/GSE42639/expr_log2.rds"))
meta <- readRDS(file.path(ROOT, "processed_data/GSE42639/sample_meta.rds"))
mo   <- which(meta$cell == "Mo")
treat_R <- c(Sham="Sham", TX91="TX91",
             "0.6LD50PR8"="PR8_0p6","10LD50PR8"="PR8_10","100LD50PR8"="PR8_100")
m  <- meta[mo, ]; m$treat <- treat_R[m$treat]
E  <- expr[, mo]
tx91_mean <- rowMeans(E[, m$treat == "TX91"])

g <- A2A3
out <- data.frame(
  gene_symbol           = g,
  trend_slope           = ia$ttS2$slope[match(g, ia$ttS2$gene)],
  trend_FDR             = ia$ttS2$adj.P.Val[match(g, ia$ttS2$gene)],
  log2FC_0p6_vs_TX91    = ia$tt06$logFC[match(g, ia$tt06$gene)],
  FDR_0p6_vs_TX91       = ia$tt06$adj.P.Val[match(g, ia$tt06$gene)],
  log2FC_10_vs_TX91     = ia$tt10$logFC[match(g, ia$tt10$gene)],
  FDR_10_vs_TX91        = ia$tt10$adj.P.Val[match(g, ia$tt10$gene)],
  log2FC_100_vs_TX91    = ia$tt100$logFC[match(g, ia$tt100$gene)],
  FDR_100_vs_TX91       = ia$tt100$adj.P.Val[match(g, ia$tt100$gene)],
  log2FC_TX91_vs_Sham   = ia$ttS1$logFC[match(g, ia$ttS1$gene)],
  FDR_TX91_vs_Sham      = ia$ttS1$adj.P.Val[match(g, ia$ttS1$gene)],
  meanExpr_TX91         = as.numeric(tx91_mean[g]),
  stringsAsFactors      = FALSE
)
out <- out[order(out$trend_slope, out$FDR_100_vs_TX91, -out$meanExpr_TX91), ]

OUT <- file.path(ROOT, "results/intersectA/Mo_intersectA_A2_inter_A3_49genes.tsv")
write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf(">>> 49 基因表已保存: %s\n", OUT))
print(out, row.names = FALSE)

# 与新版 core_targets (564) 比对
core_new <- read.delim(file.path(ROOT, "results/Mo_protective/core_targets.tsv"))
shared <- intersect(out$gene_symbol, core_new$gene_symbol)
cat(sprintf("\n与新版 564 候选交集: %d / 49\n", length(shared)))
cat("仅在旧 49 中:\n"); print(setdiff(out$gene_symbol, core_new$gene_symbol))
