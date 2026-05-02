#!/usr/bin/env Rscript
# GSE42639 Mo 细胞 -> 交集 A: 三步严格筛选
#
# 数据: processed_data/GSE42639/expr_log2.rds (neqc 标准化), Mo 细胞 5 组
#       Sham / TX91 / PR8_0.6LD50 / PR8_10LD50 / PR8_100LD50, 每组 3 重复.
#
# Step 1 (轻症特异性上调集 A1):
#   limma TX91 vs Sham,  log2FC ≥ 1 & adj.P.Val < 0.05  -> 取上调.
#
# Step 2 (剂量持续下降集 A2, 核心):
#   仅 4 个感染条件 (Sham 不参与), 严重度连续编码:
#     TX91 = 0,  PR8_0.6LD50 = 1,  PR8_10LD50 = 2,  PR8_100LD50 = 3
#   limma ~ severity (BH校正后) slope<0 & adj.P.Val<0.05.
#
# Step 3 (阶梯式下降验证集 A3):
#   3 组两两差异均通过 (基于 ~0+treat, eBayes 一次拟合):
#     PR8_0.6 - TX91 :  logFC ≤ -0.585 & adj.P.Val < 0.05
#     PR8_10  - TX91 :  logFC ≤ -1     & adj.P.Val < 0.05
#     PR8_100 - TX91 :  logFC ≤ -1     & adj.P.Val < 0.05
#
# Step 4 (核心靶点集): A1 ∩ A2 ∩ A3
#
# 输出:
#   results/intersectA/Mo_step1_TX91_vs_Sham.tsv          完整表
#   results/intersectA/Mo_step2_severity_trend.tsv         完整表
#   results/intersectA/Mo_step3_PR8x_vs_TX91.tsv x3        完整表 (3 个对比合并 long)
#   results/intersectA/intersectA_core_targets.tsv         核心靶点 (含全部所需列)
#   results/intersectA/intersectA_summary.txt              三步流计数 + 文字摘要
#   results/intersectA/intersectA_results.rds              所有 fit/topTable 对象

suppressPackageStartupMessages({ library(limma) })

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
OUT  <- file.path(ROOT, "results/intersectA")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 旧产物清理 (避免历史 step 文件混入)
old <- list.files(OUT, full.names = TRUE)
file.remove(old)

expr <- readRDS(file.path(ROOT, "processed_data/GSE42639/expr_log2.rds"))
meta <- readRDS(file.path(ROOT, "processed_data/GSE42639/sample_meta.rds"))

# 1. 数据准备 -----------------------------------------------------------------
mo_idx <- which(meta$cell == "Mo")
m  <- meta[mo_idx, ]
E  <- expr[, mo_idx]
treat_R <- c(Sham="Sham", TX91="TX91",
             "0.6LD50PR8"="PR8_0p6", "10LD50PR8"="PR8_10",
             "100LD50PR8"="PR8_100")
m$treat <- factor(treat_R[m$treat],
                  levels = c("Sham","TX91","PR8_0p6","PR8_10","PR8_100"))
cat(">>> Mo 样本分布:\n"); print(table(m$treat))

# 2. 主模型 (~0 + treat) ------------------------------------------------------
designAll <- model.matrix(~0 + m$treat)
colnames(designAll) <- levels(m$treat)
fitAll <- lmFit(E, designAll)
cm <- makeContrasts(
  TX91_vs_Sham   = TX91 - Sham,
  PR8_0p6_vs_TX91 = PR8_0p6 - TX91,
  PR8_10_vs_TX91  = PR8_10  - TX91,
  PR8_100_vs_TX91 = PR8_100 - TX91,
  levels = designAll
)
fit2 <- eBayes(contrasts.fit(fitAll, cm), trend = TRUE)

get_tt <- function(coef) {
  t <- topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  t$gene <- rownames(t)
  t[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
}
ttS1 <- get_tt("TX91_vs_Sham")
tt06 <- get_tt("PR8_0p6_vs_TX91")
tt10 <- get_tt("PR8_10_vs_TX91")
tt100 <- get_tt("PR8_100_vs_TX91")

# Step 1: 轻症特异性上调
A1 <- ttS1$gene[ttS1$logFC >= 1 & ttS1$adj.P.Val < 0.05]
write.table(ttS1[order(ttS1$adj.P.Val), ], file.path(OUT, "Mo_step1_TX91_vs_Sham.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# 3. Step 2 趋势检验 ----------------------------------------------------------
sev_map <- c(TX91=0, PR8_0p6=1, PR8_10=2, PR8_100=3)
inf_idx <- which(m$treat %in% names(sev_map))
m_inf   <- m[inf_idx, ]
E_inf   <- E[, inf_idx]
m_inf$severity <- sev_map[as.character(m_inf$treat)]
designSev <- model.matrix(~ m_inf$severity)
colnames(designSev) <- c("Intercept","severity")
fitSev <- eBayes(lmFit(E_inf, designSev), trend = TRUE)
ttS2 <- topTable(fitSev, coef = "severity", number = Inf, sort.by = "none")
ttS2$gene <- rownames(ttS2)
colnames(ttS2)[colnames(ttS2)=="logFC"] <- "slope"
ttS2 <- ttS2[, c("gene","slope","AveExpr","t","P.Value","adj.P.Val","B")]
A2 <- ttS2$gene[ttS2$slope < 0 & ttS2$adj.P.Val < 0.05]
write.table(ttS2[order(ttS2$adj.P.Val), ],
            file.path(OUT, "Mo_step2_severity_trend.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# 4. Step 3 阶梯式两两验证 ----------------------------------------------------
write.table(tt06[order(tt06$adj.P.Val), ],
            file.path(OUT, "Mo_step3a_PR8_0p6_vs_TX91.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(tt10[order(tt10$adj.P.Val), ],
            file.path(OUT, "Mo_step3b_PR8_10_vs_TX91.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(tt100[order(tt100$adj.P.Val), ],
            file.path(OUT, "Mo_step3c_PR8_100_vs_TX91.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

g06  <- tt06 $gene[tt06 $logFC <= -0.585 & tt06 $adj.P.Val < 0.05]
g10  <- tt10 $gene[tt10 $logFC <= -1     & tt10 $adj.P.Val < 0.05]
g100 <- tt100$gene[tt100$logFC <= -1     & tt100$adj.P.Val < 0.05]
A3 <- Reduce(intersect, list(g06, g10, g100))

# 5. Step 4 三重交集 ----------------------------------------------------------
core <- Reduce(intersect, list(A1, A2, A3))

# 核心靶点表 (全字段)
build_row <- function(g) {
  data.frame(
    gene = g,
    A_TX91vsSham_logFC   = ttS1$logFC[match(g, ttS1$gene)],
    A_TX91vsSham_FDR     = ttS1$adj.P.Val[match(g, ttS1$gene)],
    B_severity_slope     = ttS2$slope[match(g, ttS2$gene)],
    B_severity_FDR       = ttS2$adj.P.Val[match(g, ttS2$gene)],
    C_PR8_0p6_vs_TX91_logFC = tt06$logFC[match(g, tt06$gene)],
    C_PR8_0p6_vs_TX91_FDR   = tt06$adj.P.Val[match(g, tt06$gene)],
    C_PR8_10_vs_TX91_logFC  = tt10$logFC[match(g, tt10$gene)],
    C_PR8_10_vs_TX91_FDR    = tt10$adj.P.Val[match(g, tt10$gene)],
    C_PR8_100_vs_TX91_logFC = tt100$logFC[match(g, tt100$gene)],
    C_PR8_100_vs_TX91_FDR   = tt100$adj.P.Val[match(g, tt100$gene)],
    stringsAsFactors = FALSE
  )
}
core_tbl <- if (length(core)) build_row(core) else build_row(character(0))
# 排序: 综合得分 = A_logFC - 平均 PR8 vs TX91 logFC (induced 高 + 跌幅大 → 排前)
if (nrow(core_tbl)) {
  core_tbl$score <- core_tbl$A_TX91vsSham_logFC -
    rowMeans(core_tbl[, c("C_PR8_0p6_vs_TX91_logFC",
                          "C_PR8_10_vs_TX91_logFC",
                          "C_PR8_100_vs_TX91_logFC")])
  core_tbl <- core_tbl[order(-core_tbl$score), ]
}
write.table(core_tbl, file.path(OUT, "intersectA_core_targets.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# 6. 摘要 + RDS ---------------------------------------------------------------
sink(file.path(OUT, "intersectA_summary.txt"))
cat("=== 交集 A 三步严格筛选 (Mo 细胞) ===\n\n")
cat(sprintf("Step 1 轻症特异性上调集 A1 (TX91vsSham logFC>=1 & FDR<0.05): %d\n", length(A1)))
cat(sprintf("Step 2 剂量持续下降集    A2 (slope<0 & FDR<0.05)            : %d\n", length(A2)))
cat(sprintf("Step 3 阶梯式下降验证集  A3 (PR8_0.6 vs TX91 ≤-0.585 & PR8_10/100 vs TX91 ≤-1, 全 FDR<0.05): %d\n", length(A3)))
cat(sprintf("  其中:\n"))
cat(sprintf("    PR8_0.6 vs TX91 满足 ≤-0.585: %d\n", length(g06)))
cat(sprintf("    PR8_10  vs TX91 满足 ≤-1    : %d\n", length(g10)))
cat(sprintf("    PR8_100 vs TX91 满足 ≤-1    : %d\n", length(g100)))
cat(sprintf("\nStep 4 核心靶点集 = A1 ∩ A2 ∩ A3 : %d 基因\n", length(core)))
if (length(core)) {
  cat("\n=== 核心靶点完整表 ===\n")
  print(core_tbl, row.names = FALSE)
} else {
  cat("\n核心靶点为空. 中间集合的两两交集情况:\n")
  cat(sprintf("  A1 ∩ A2 : %d\n", length(intersect(A1,A2))))
  cat(sprintf("  A1 ∩ A3 : %d\n", length(intersect(A1,A3))))
  cat(sprintf("  A2 ∩ A3 : %d\n", length(intersect(A2,A3))))
}
sink()
cat(readLines(file.path(OUT, "intersectA_summary.txt")), sep="\n")

saveRDS(list(
  ttS1 = ttS1, ttS2 = ttS2, tt06 = tt06, tt10 = tt10, tt100 = tt100,
  A1 = A1, A2 = A2, A3 = A3, core = core, core_tbl = core_tbl
), file.path(OUT, "intersectA_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
