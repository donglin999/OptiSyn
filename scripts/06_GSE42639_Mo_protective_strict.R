#!/usr/bin/env Rscript
# 06_GSE42639_Mo_protective_strict.R
#
# 在 564 候选基础上叠加两步严选:
#   Step A 阶梯下降验证:
#     PR8_0p6 vs TX91 log2FC <= -0.3
#     PR8_10  vs TX91 log2FC <= -0.585
#     全部 FDR < 0.05
#   Step B 功能注释白名单 (剔除管家基因 / 代谢酶):
#     仅保留命中以下任一 GO 大类的基因:
#       Secreted          : 分泌型蛋白 / 细胞因子配体
#       MembraneReceptor  : 膜受体 / 跨膜信号受体
#       TF                : 转录因子 / 转录调控子
#       NegRegInflamm     : 炎症 / 免疫反应负向调控因子
#
# 输入: results/Mo_protective/Mo_protective_results.rds
# 输出: results/Mo_protective/strict/

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(RColorBrewer)
  library(org.Mm.eg.db)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
IN_RDS <- file.path(ROOT, "results/Mo_protective/Mo_protective_results.rds")
OUT <- file.path(ROOT, "results/Mo_protective/strict")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(IN_RDS))
res <- readRDS(IN_RDS)
core   <- res$core
tt06   <- res$tt06
tt10   <- res$tt10
tt100  <- res$tt100
ttSham <- res$ttSham
ttSev  <- res$ttSev
m      <- res$meta
E      <- res$expr
GROUP_LEVELS <- res$GROUP_LEVELS

n_start <- nrow(core)
cat(sprintf(">>> 起始候选 (Mo_protective core): %d 基因\n", n_start))

# ============================================================================
# Step A: 阶梯下降
# ============================================================================
A_thr <- list(log2fc_0p6 = -0.3, log2fc_10 = -0.585, fdr = 0.05)

g0 <- core$gene_symbol
m06 <- match(g0, tt06$gene)
m10 <- match(g0, tt10$gene)

a_pass_06 <- (tt06$logFC[m06] <= A_thr$log2fc_0p6) & (tt06$adj.P.Val[m06] < A_thr$fdr)
a_pass_10 <- (tt10$logFC[m10] <= A_thr$log2fc_10 ) & (tt10$adj.P.Val[m10] < A_thr$fdr)
a_pass_06[is.na(a_pass_06)] <- FALSE
a_pass_10[is.na(a_pass_10)] <- FALSE

S_A <- g0[a_pass_06 & a_pass_10]
cat(sprintf(">>> Step A 阶梯下降 (0p6 logFC<=%.2f, 10 logFC<=%.3f, FDR<%.2f): %d -> %d\n",
            A_thr$log2fc_0p6, A_thr$log2fc_10, A_thr$fdr, n_start, length(S_A)))
if (length(S_A) == 0) stop("Step A 后无基因留存")

# ============================================================================
# Step B: 黑名单剔除模式
#   规则:
#     1) 命中白名单任一类   -> 保留 (即使同时命中黑名单, 白优先)
#     2) 仅命中黑名单       -> 剔除
#     3) 都未命中           -> 保留 (默认通过, 不强制功能注释)
# ============================================================================
white_groups <- list(
  ImmuneReg = c(
    "GO:0002376",  # immune system process
    "GO:0050776",  # regulation of immune response
    "GO:0002682",  # regulation of immune system process
    "GO:0045088",  # regulation of innate immune response
    "GO:0002683",  # negative regulation of immune system process
    "GO:0050778"   # positive regulation of immune response
  ),
  InflammReg = c(
    "GO:0006954",  # inflammatory response
    "GO:0050727",  # regulation of inflammatory response
    "GO:0050728",  # negative regulation of inflammatory response
    "GO:0050729",  # positive regulation of inflammatory response
    "GO:0001818",  # negative regulation of cytokine production
    "GO:0001819",  # positive regulation of cytokine production
    "GO:0001817"   # regulation of cytokine production
  ),
  SignalTransduction = c(
    "GO:0007165",  # signal transduction
    "GO:0023052",  # signaling
    "GO:0009966",  # regulation of signal transduction
    "GO:0035556"   # intracellular signal transduction
  ),
  CellCommunication = c(
    "GO:0007154",  # cell communication
    "GO:0010646",  # regulation of cell communication
    "GO:0007267"   # cell-cell signaling
  ),
  Secreted = c(
    "GO:0005576",  # extracellular region
    "GO:0005615",  # extracellular space
    "GO:0005125",  # cytokine activity
    "GO:0048018"   # receptor ligand activity
  ),
  Membrane = c(
    "GO:0005886",  # plasma membrane
    "GO:0009986",  # cell surface
    "GO:0038023",  # signaling receptor activity
    "GO:0004888",  # transmembrane signaling receptor activity
    "GO:0005887"   # integral component of plasma membrane
  ),
  TFReg = c(
    "GO:0003700",  # DNA-binding transcription factor activity
    "GO:0140110",  # transcription regulator activity
    "GO:0006355",  # regulation of DNA-templated transcription
    "GO:0006357"   # regulation of transcription by RNA polymerase II
  )
)

black_groups <- list(
  Ribosome = c(
    "GO:0003735",  # structural constituent of ribosome
    "GO:0005840",  # ribosome
    "GO:0022626",  # cytosolic ribosome
    "GO:0006412"   # translation
  ),
  MitoStructure = c(
    "GO:0006119",  # oxidative phosphorylation
    "GO:0042775",  # mitochondrial ATP synthesis coupled electron transport
    "GO:0008137",  # NADH dehydrogenase (ubiquinone) activity
    "GO:0098798",  # mitochondrial protein-containing complex
    "GO:0005747",  # mitochondrial respiratory chain complex I
    "GO:0005746"   # mitochondrial respirasome
  ),
  CytoskeletalStructure = c(
    "GO:0005200",  # structural constituent of cytoskeleton
    "GO:0005198"   # structural molecule activity
  ),
  CarbMetabolism = c(
    "GO:0006096",  # glycolytic process
    "GO:0006094",  # gluconeogenesis
    "GO:0005975",  # carbohydrate metabolic process
    "GO:0006006"   # glucose metabolic process
  ),
  LipidMetabolism = c(
    "GO:0006631",  # fatty acid metabolic process
    "GO:0008610",  # lipid biosynthetic process
    "GO:0006629",  # lipid metabolic process
    "GO:0019216"   # regulation of lipid metabolic process (调控类含义需要观察, 但仍按用户意图归入"代谢酶")
  ),
  AAMetabolism = c(
    "GO:0006520",  # cellular amino acid metabolic process
    "GO:0009063",  # cellular amino acid catabolic process
    "GO:0008652"   # cellular amino acid biosynthetic process
  )
)

go2syms <- function(go_ids) {
  out <- tryCatch({
    suppressMessages(suppressWarnings(
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = go_ids,
                            keytype = "GOALL",
                            columns = c("SYMBOL"))
    ))
  }, error = function(e) {
    suppressMessages(suppressWarnings(
      AnnotationDbi::select(org.Mm.eg.db,
                            keys = go_ids,
                            keytype = "GO",
                            columns = c("SYMBOL"))
    ))
  })
  unique(out$SYMBOL[!is.na(out$SYMBOL)])
}

cat(">>> Step B 白名单 GO 覆盖小鼠基因数:\n")
white_genes <- lapply(white_groups, go2syms)
for (nm in names(white_genes)) cat(sprintf("    [W] %-22s %d\n", nm, length(white_genes[[nm]])))
cat(">>> Step B 黑名单 GO 覆盖小鼠基因数:\n")
black_genes <- lapply(black_groups, go2syms)
for (nm in names(black_genes)) cat(sprintf("    [B] %-22s %d\n", nm, length(black_genes[[nm]])))

# 贴标签
classify_gene <- function(g) {
  w_hits <- names(white_genes)[vapply(white_genes, function(s) g %in% s, logical(1))]
  b_hits <- names(black_genes)[vapply(black_genes, function(s) g %in% s, logical(1))]
  list(white = w_hits, black = b_hits)
}
labels <- lapply(S_A, classify_gene)
white_str <- vapply(labels, function(x) paste(x$white, collapse=";"), character(1))
black_str <- vapply(labels, function(x) paste(x$black, collapse=";"), character(1))

# 决策: 白优先, 仅黑则剔, 都无则保留
n_white <- nchar(white_str) > 0
n_black <- nchar(black_str) > 0
keep_B <- n_white | (!n_white & !n_black)        # 命中白 或 都未命中 -> 保留
drop_reason_only_black <- (!n_white) & n_black   # 仅命中黑 -> 剔除
S_B <- S_A[keep_B]
S_dropped <- S_A[drop_reason_only_black]
cat(sprintf(">>> Step B 黑名单剔除: %d -> %d  (剔仅命中黑名单的 %d, 默认保留无注释 %d)\n",
            length(S_A), length(S_B), sum(drop_reason_only_black),
            sum(!n_white & !n_black)))
if (length(S_B) == 0) stop("Step B 后无基因留存")

# 写出剔除清单 (备查)
if (length(S_dropped) > 0) {
  drop_tbl <- data.frame(
    gene_symbol = S_dropped,
    black_hits  = black_str[match(S_dropped, S_A)],
    log2FC_100_vs_TX91 = tt100$logFC[match(S_dropped, tt100$gene)],
    FDR_100_vs_TX91   = tt100$adj.P.Val[match(S_dropped, tt100$gene)],
    stringsAsFactors = FALSE
  )
  write.table(drop_tbl, file.path(OUT, "stepB_dropped_blacklist_only.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

class_str <- white_str  # 主表使用白名单类别字符串
black_for_keep <- black_str  # 保留集合中, 同时也命中黑名单的基因 (作参考列, 不剔)

# ============================================================================
# 最终主表
# ============================================================================
get_core <- function(col, gv) core[[col]][match(gv, core$gene_symbol)]
final_tbl <- data.frame(
  gene_symbol           = S_B,
  func_class            = class_str[match(S_B, S_A)],
  also_in_blacklist     = black_for_keep[match(S_B, S_A)],
  trend_slope           = get_core("trend_slope", S_B),
  trend_FDR             = get_core("trend_FDR", S_B),
  log2FC_0p6_vs_TX91    = tt06$logFC[match(S_B, tt06$gene)],
  FDR_0p6_vs_TX91       = tt06$adj.P.Val[match(S_B, tt06$gene)],
  log2FC_10_vs_TX91     = tt10$logFC[match(S_B, tt10$gene)],
  FDR_10_vs_TX91        = tt10$adj.P.Val[match(S_B, tt10$gene)],
  log2FC_100_vs_TX91    = get_core("log2FC_100LD50_vs_TX91", S_B),
  FDR_100_vs_TX91       = get_core("FDR_100LD50_vs_TX91", S_B),
  log2FC_TX91_vs_Sham   = get_core("log2FC_TX91_vs_Sham", S_B),
  FDR_TX91_vs_Sham      = get_core("FDR_TX91_vs_Sham", S_B),
  meanExpr_TX91         = get_core("meanExpr_TX91", S_B),
  stringsAsFactors      = FALSE
)
final_tbl <- final_tbl[order(final_tbl$trend_slope,
                              final_tbl$FDR_100_vs_TX91,
                              -final_tbl$meanExpr_TX91), ]
rownames(final_tbl) <- NULL
write.table(final_tbl, file.path(OUT, "core_targets_strict.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(final_tbl, file.path(OUT, "core_targets_strict.csv"),
          quote = TRUE, row.names = FALSE)

# 各功能类计数 (含重叠)
class_count <- sort(table(unlist(strsplit(final_tbl$func_class, ";"))),
                    decreasing = TRUE)
cat("\n>>> 各功能类基因数 (允许重叠, 空 = 无注释默认保留):\n")
print(class_count)
cat(sprintf("\n无任何白名单注释 (默认保留): %d\n",
            sum(final_tbl$func_class == "")))

# 命中炎症 / 免疫调控的高优先候选
inflam_immune <- final_tbl[grepl("InflammReg|ImmuneReg", final_tbl$func_class), ]
write.table(inflam_immune,
            file.path(OUT, "core_targets_strict_InflammImmune.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf(">>> 命中炎症/免疫调控 (湿实验高优): %d\n", nrow(inflam_immune)))

# ============================================================================
# 流程日志 + 参数
# ============================================================================
sink(file.path(OUT, "flow_summary_strict.txt"))
cat("=== Mo 严选靶点 (564 -> 阶梯下降 -> 黑名单剔除) ===\n")
cat(sprintf("时间: %s\n\n", TS))
cat(sprintf("%-65s %s\n", "节点", "基因数"))
cat(sprintf("%-65s %s\n", "----", "------"))
cat(sprintf("%-65s %d\n", "起始 (Mo_protective core_targets)", n_start))
cat(sprintf("%-65s %d\n",
            sprintf("Step A 阶梯下降 (0.6 logFC<=%.2f, 10 logFC<=%.3f, FDR<0.05)",
                    A_thr$log2fc_0p6, A_thr$log2fc_10), length(S_A)))
cat(sprintf("%-65s %d\n",
            "Step B 黑名单剔除 (核糖体/线粒体结构/胞内结构/糖脂氨基酸代谢)",
            length(S_B)))
cat(sprintf("        其中: 命中白名单 %d, 仅命中黑名单(剔) %d, 都未命中(默认保留) %d\n",
            sum(nchar(class_str[match(S_B, S_A)]) > 0),
            length(S_dropped),
            sum(final_tbl$func_class == "")))
cat("\n白名单各类基因数 (允许重叠):\n")
print(class_count)
cat(sprintf("\n命中炎症/免疫调控: %d (写入 core_targets_strict_InflammImmune.tsv)\n",
            nrow(inflam_immune)))
cat(sprintf("剔除清单 (仅命中黑名单, %d 个): stepB_dropped_blacklist_only.tsv\n",
            length(S_dropped)))
cat("\n=== 严选核心靶点 (前 50) ===\n")
print(head(final_tbl, 50), row.names = FALSE)
sink()
cat("\n", paste(readLines(file.path(OUT, "flow_summary_strict.txt")), collapse="\n"), sep="")

# ============================================================================
# 图: 严选热图 + Top 折线图 + 类别条形
# ============================================================================
make_heatmap <- function(genes, E, m, out_prefix, title_main) {
  if (length(genes) == 0) return(invisible(NULL))
  M <- E[genes, , drop = FALSE]
  col_ord <- unlist(lapply(GROUP_LEVELS, function(lv) {
    idx <- which(m$treat == lv); idx[order(m$sample_id[idx])]
  }))
  M <- M[, col_ord, drop = FALSE]
  m_ord <- m[col_ord, , drop = FALSE]
  Z <- t(scale(t(M))); Z[is.na(Z)] <- 0
  if (nrow(Z) > 2) {
    hc <- hclust(dist(Z), method = "average")
    Z <- Z[hc$order, , drop = FALSE]
  }
  df <- as.data.frame(Z); df$gene <- rownames(Z)
  df_long <- reshape(df, varying = colnames(Z), v.names = "z",
                     timevar = "sample", times = colnames(Z), direction = "long")
  df_long$gene   <- factor(df_long$gene, levels = rev(rownames(Z)))
  df_long$sample <- factor(df_long$sample, levels = colnames(Z))
  df_long$group  <- m_ord$treat[match(df_long$sample, m_ord$sample_id)]
  cap <- 2.5
  df_long$z_cap <- pmax(pmin(df_long$z, cap), -cap)

  p <- ggplot(df_long, aes(sample, gene, fill = z_cap)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-cap, cap),
                         oob = scales::squish, name = "row z-score") +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = title_main,
         subtitle = sprintf("n=%d, 行 z-score, 列序固定, 列不聚类", nrow(Z))) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90", color = NA),
          strip.text = element_text(face = "bold"))
  h_in <- max(3, min(20, nrow(Z) * 0.18 + 2))
  ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = h_in,
         limitsize = FALSE, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 8, height = h_in,
         dpi = 200, limitsize = FALSE, type = "cairo")
}

make_lineplot <- function(genes, E, m, out_prefix, title_main, top_n = 20) {
  g_top <- head(genes, top_n)
  if (length(g_top) == 0) return(invisible(NULL))
  M <- E[g_top, , drop = FALSE]
  long_list <- lapply(g_top, function(gn) {
    x <- M[gn, ]
    data.frame(gene = gn,
               sample_id = names(x),
               group = m$treat[match(names(x), m$sample_id)],
               expr = as.numeric(x), stringsAsFactors = FALSE)
  })
  long_df <- do.call(rbind, long_list)
  agg <- aggregate(expr ~ gene + group, data = long_df,
                   FUN = function(v) c(mean = mean(v), sem = sd(v)/sqrt(length(v))))
  agg <- data.frame(gene = agg$gene, group = agg$group,
                    mean = agg$expr[, "mean"], sem = agg$expr[, "sem"])
  agg$group <- factor(agg$group, levels = GROUP_LEVELS)
  agg$gene  <- factor(agg$gene, levels = g_top)

  p <- ggplot(agg, aes(group, mean, group = gene)) +
    geom_line(color = "#377EB8") +
    geom_point(color = "#377EB8", size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, color = "grey30") +
    facet_wrap(~ gene, scales = "free_y", ncol = 5) +
    labs(x = NULL, y = "log2 expr (neqc)",
         title = title_main,
         subtitle = sprintf("Top %d, 组均值 +/- SEM, n=3", length(g_top))) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(face = "bold"))
  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = 8, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 8, dpi = 200, type = "cairo")
}

make_heatmap(final_tbl$gene_symbol, E, m,
             file.path(OUT, "heatmap_strict"),
             "Mo 严选靶点剂量响应热图")
cat(">>> 热图: heatmap_strict.pdf/png\n")

make_lineplot(final_tbl$gene_symbol, E, m,
              file.path(OUT, "top_lineplot_strict"),
              "Mo 严选靶点 Top 剂量响应折线图",
              top_n = min(20, nrow(final_tbl)))
cat(">>> 折线图: top_lineplot_strict.pdf/png\n")

# 类别条形图 (含 "Unannotated" 列)
class_df <- data.frame(class = names(class_count),
                       n = as.integer(class_count))
n_unann <- sum(final_tbl$func_class == "")
if (n_unann > 0) {
  class_df <- rbind(class_df, data.frame(class = "Unannotated", n = n_unann))
}
class_df <- class_df[order(-class_df$n), ]
class_df$class <- factor(class_df$class, levels = class_df$class)
p_bar <- ggplot(class_df, aes(class, n, fill = class)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = NULL, y = "基因数 (允许重叠)",
       title = "严选靶点功能类别分布 (白名单 + 无注释)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "func_class_bar.pdf"), p_bar, width = 7, height = 4.5,
       device = cairo_pdf)
ggsave(file.path(OUT, "func_class_bar.png"), p_bar, width = 7, height = 4.5,
       dpi = 200, type = "cairo")
cat(">>> 类别分布图: func_class_bar.pdf/png\n")

# 参数 + sessionInfo
sink(file.path(OUT, "params_and_session_info_strict.txt"))
cat("=== 严选参数记录 ===\n")
cat(sprintf("时间戳: %s\n", TS))
cat("脚本: scripts/06_GSE42639_Mo_protective_strict.R\n\n")
cat("Step A 阈值:\n"); print(A_thr)
cat("\nStep B 白名单 GO (保留, 命中即留):\n"); print(white_groups)
cat("\nStep B 黑名单 GO (剔, 仅当未命中白名单):\n"); print(black_groups)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

saveRDS(list(
  S_A = S_A, S_B = S_B, S_dropped = S_dropped,
  final_tbl = final_tbl,
  white_groups = white_groups, white_genes = white_genes,
  black_groups = black_groups, black_genes = black_genes,
  A_thr = A_thr, class_count = class_count
), file.path(OUT, "Mo_protective_strict.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
