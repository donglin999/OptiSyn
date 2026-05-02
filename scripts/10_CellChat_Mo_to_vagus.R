#!/usr/bin/env Rscript
# 跨数据集 CellChat 风格 LR 通讯计算
#
# 配体 L (sender = Mo cell, GSE42639 Mo TX91 baseline):
#   58 严选靶点 ∩ CellChatDB.mouse 'Secreted Signaling' 类的 ligand
#   (Secreted Signaling 已排除 Cell-Cell Contact / ECM-Receptor / 细胞内结构蛋白)
#   也跑一遍 4 双标基因作为更严的 L 列表
#
# 受体 R (receiver = vagal sensory ganglia, GSE161878 IAV upregulated):
#   343 上调 DEG ∩ CellChatDB.mouse 受体基因 (单基因 receptor 或 complex 任一亚基)
#
# 结合概率 (CellChat Hill, 跨数据集 -> 用 rank-normalized 表达):
#   E_L_norm = rank(E_L) / N (Mo TX91 mean, 范围 [0,1])
#   E_R_norm = rank(E_R) / N (vagal IAV mean,    范围 [0,1])
#   complex R: 取所有 subunit 的几何平均 (CellChat 默认)
#   P = E_L^n/(K^n+E_L^n) * E_R^n/(K^n+E_R^n),  K=0.5, n=2
#
# 注意: 这是跨数据集风格的简化版 CellChat (单细胞通讯算法本来设计在单一矩阵内).
# rank-norm 提供了相对可比的基线; P 是相对值, 用于排序/筛选高可信对.
#
# 输出: results/CellChat/

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"
OUT  <- file.path(ROOT, "results/CellChat")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(K = 0.5, n_hill = 2,
               complex_combine = "geometric_mean",
               prob_high_thresh = 0.10)

# 1. 加载数据 ----------------------------------------------------------------
load(file.path(ROOT, "processed_data/CellChatDB/CellChatDB.mouse.rda"))
ia    <- CellChatDB.mouse$interaction
cmplx <- CellChatDB.mouse$complex
cat(sprintf(">>> CellChatDB.mouse: %d interactions\n", nrow(ia)))
cat(">>> annotation 分布:\n"); print(table(ia$annotation))

# Secreted Signaling 子集 (主, 严判)
# 也准备 Secreted + Cell-Cell Contact (放宽: 含膜识别配体如 DC-SIGN/Cd209a)
ia_ss   <- ia[ia$annotation == "Secreted Signaling", , drop = FALSE]
ia_ss_cc <- ia[ia$annotation %in% c("Secreted Signaling","Cell-Cell Contact"),
               , drop = FALSE]
cat(sprintf(">>> Secreted Signaling: %d, +Cell-Cell Contact: %d\n",
            nrow(ia_ss), nrow(ia_ss_cc)))

# Mo 严选 58 + 4 双标
strict <- read.delim(file.path(ROOT,
                       "results/Mo_protective/strict/core_targets_strict.tsv"),
                     stringsAsFactors = FALSE)
genes58 <- strict$gene_symbol
both4   <- read.delim(file.path(ROOT,
                       "results/Mo_protective/strict/GSE31022_validation/GSE31022_both_conserved_progressive.tsv"),
                     stringsAsFactors = FALSE)$gene_symbol

# Mo 表达 (GSE42639 Mo TX91 group mean)
mo_res <- readRDS(file.path(ROOT, "results/Mo_protective/Mo_protective_results.rds"))
mo_E   <- mo_res$expr
mo_m   <- mo_res$meta
tx91_idx <- which(mo_m$treat == "TX91")
mo_tx91_mean <- rowMeans(mo_E[, tx91_idx, drop = FALSE])
names(mo_tx91_mean) <- rownames(mo_E)
cat(sprintf(">>> GSE42639 Mo 矩阵: %d 基因, TX91 mean (n=3) 计算完毕\n", length(mo_tx91_mean)))

# 迷走神经表达 (GSE161878 raw count -> log2(CPM+1) IAV group mean)
# 使用 raw count 而非 DESeq2 过滤后的 vst, 保留低表达受体 (rank-norm 自然给低分)
cnt_path <- file.path(ROOT, "raw_data/00_raw_data/GSE161878_gene_counts.txt")
meta_path <- file.path(ROOT, "raw_data/00_raw_data/GSE161878_sample_metadata.txt")
cnt_raw <- read.delim(cnt_path, check.names = FALSE, stringsAsFactors = FALSE)
vst_m   <- read.delim(meta_path, stringsAsFactors = FALSE)
genes_e <- as.character(cnt_raw$GeneID)
M_raw   <- as.matrix(cnt_raw[, -1, drop = FALSE]); rownames(M_raw) <- genes_e
vst_m   <- vst_m[match(colnames(M_raw), vst_m$sample_id), , drop = FALSE]
iav_idx <- which(vst_m$group == "IAV")

# log2(CPM+1)
lib_size <- colSums(M_raw)
cpm <- t(t(M_raw) / lib_size) * 1e6
log2cpm <- log2(cpm + 1)

ann_g <- read.delim(file.path(ROOT, "processed_data/GSE161878/gene_annotation.tsv"),
                    stringsAsFactors = FALSE)
sym_lookup <- setNames(ann_g$SYMBOL, ann_g$ENTREZID)
g_sym <- sym_lookup[rownames(log2cpm)]
keep_v  <- !is.na(g_sym) & nchar(g_sym) > 0
log2cpm <- log2cpm[keep_v, , drop = FALSE]
g_sym   <- g_sym[keep_v]
# 同基因多 ENTREZ -> 取最大 CPM 行 (代表性表达)
ord <- order(g_sym, -rowMeans(log2cpm))
log2cpm <- log2cpm[ord, ]; g_sym <- g_sym[ord]
keep_first <- !duplicated(g_sym)
log2cpm <- log2cpm[keep_first, ]; rownames(log2cpm) <- g_sym[keep_first]

vagus_iav_mean <- rowMeans(log2cpm[, iav_idx, drop = FALSE])
cat(sprintf(">>> GSE161878 raw->log2CPM: %d SYMBOL, IAV mean (n=%d)\n",
            length(vagus_iav_mean), length(iav_idx)))

# GSE161878 上调 DEG (343)
deg_up <- read.delim(file.path(ROOT, "results/GSE161878/GSE161878_DEG_up.tsv"),
                     stringsAsFactors = FALSE)
up_sym <- unique(deg_up$SYMBOL[!is.na(deg_up$SYMBOL)])
cat(sprintf(">>> GSE161878 上调 DEG SYMBOL: %d\n", length(up_sym)))

# 2. 受体展开 (单基因 + complex 亚基) ----------------------------------------
expand_receptor <- function(rec_name) {
  if (rec_name %in% rownames(cmplx)) {
    subs <- as.character(cmplx[rec_name, ])
    subs <- subs[nchar(subs) > 0]
    return(subs)
  } else {
    return(rec_name)
  }
}

# 3. 候选配体 L (58, 4) --------------------------------------------------------
# 把 ligand 也展开 (少数 ligand 是 complex)
expand_ligand <- function(lig_name) {
  if (lig_name %in% rownames(cmplx)) {
    subs <- as.character(cmplx[lig_name, ])
    subs <- subs[nchar(subs) > 0]
    return(subs)
  } else {
    return(lig_name)
  }
}

build_ia_table <- function(L_pool, ia_pool) {
  rows <- list()
  for (i in seq_len(nrow(ia_pool))) {
    lig  <- ia_pool$ligand[i]
    lig_subs <- expand_ligand(lig)
    if (!any(lig_subs %in% L_pool)) next
    rec  <- ia_pool$receptor[i]
    rec_subs <- expand_receptor(rec)
    rows[[length(rows) + 1L]] <- data.frame(
      ligand_complex = lig,
      ligand_subunits = paste(lig_subs, collapse = ";"),
      receptor_complex = rec,
      receptor_subunits = paste(rec_subs, collapse = ";"),
      n_recep_subunits = length(rec_subs),
      pathway_name = ia_pool$pathway_name[i],
      interaction_name = ia_pool$interaction_name[i],
      annotation = ia_pool$annotation[i],
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# 严判 (Secreted Signaling 唯一)
ia_58_ss <- build_ia_table(genes58, ia_ss)
ia_4_ss  <- build_ia_table(both4,   ia_ss)
# 放宽 (Secreted + Cell-Cell Contact)
ia_58_cc <- build_ia_table(genes58, ia_ss_cc)
ia_4_cc  <- build_ia_table(both4,   ia_ss_cc)

cat(sprintf(">>> 58 ∩ SS ligand: %d LR pairs   ; +CCC: %d LR pairs\n",
            nrow(ia_58_ss), nrow(ia_58_cc)))
cat(sprintf(">>> 4  ∩ SS ligand: %d LR pairs   ; +CCC: %d LR pairs\n",
            nrow(ia_4_ss),  nrow(ia_4_cc)))

# 4. R 上调过滤 + 结合概率 --------------------------------------------------
# rank-normalize 表达
rank_norm <- function(x) {
  r <- rank(x, ties.method = "average")
  r / length(x)
}
mo_rn    <- rank_norm(mo_tx91_mean)
vagus_rn <- rank_norm(vagus_iav_mean)

hill <- function(E, K = PARAMS$K, n = PARAMS$n_hill) {
  E^n / (K^n + E^n)
}

geo_mean <- function(v) {
  v <- v[!is.na(v) & v > 0]
  if (length(v) == 0) return(NA_real_)
  exp(mean(log(v)))
}

# DEG 全表 (R 方向 / 显著性)
deg_full <- read.delim(file.path(ROOT, "results/GSE161878/GSE161878_DEA_full.tsv"),
                       stringsAsFactors = FALSE)
deg_full <- deg_full[!is.na(deg_full$SYMBOL) & nchar(deg_full$SYMBOL) > 0, ]

compute_lr_table <- function(ia_tbl, label) {
  if (nrow(ia_tbl) == 0) return(data.frame())
  out <- list()
  for (i in seq_len(nrow(ia_tbl))) {
    lig_subs <- strsplit(ia_tbl$ligand_subunits[i], ";")[[1]]
    rec_subs <- strsplit(ia_tbl$receptor_subunits[i], ";")[[1]]

    # ligand: Mo TX91 表达 (raw + norm), 多亚基 geo mean
    lig_raw_v <- mo_tx91_mean[lig_subs]
    lig_rn_v  <- mo_rn[lig_subs]
    if (!all(lig_subs %in% names(mo_tx91_mean))) next
    lig_raw <- geo_mean(lig_raw_v)
    lig_rn  <- geo_mean(lig_rn_v)

    # receptor: vagal IAV 表达, 多亚基 geo mean
    rec_raw_v <- vagus_iav_mean[rec_subs]
    rec_rn_v  <- vagus_rn[rec_subs]
    n_in_vagus <- sum(!is.na(rec_raw_v))
    if (n_in_vagus == 0) next
    rec_raw <- geo_mean(rec_raw_v)
    rec_rn  <- geo_mean(rec_rn_v)

    # R 在 GSE161878 中的差异方向/显著性 (取所有 subunit 中 max logFC)
    deg_sub <- deg_full[deg_full$SYMBOL %in% rec_subs, , drop = FALSE]
    if (nrow(deg_sub) > 0) {
      idx_max <- which.max(deg_sub$log2FC_shrink)
      r_max_lfc <- deg_sub$log2FC_shrink[idx_max]
      r_min_fdr <- deg_sub$FDR         [idx_max]
      r_dir_up  <- !is.na(r_max_lfc) && r_max_lfc > 0
    } else {
      r_max_lfc <- NA_real_; r_min_fdr <- NA_real_; r_dir_up <- FALSE
    }

    # R 严判: 全亚基上调
    n_up_sub <- sum(rec_subs %in% up_sym)
    all_up   <- n_up_sub == length(rec_subs)
    any_up   <- n_up_sub >= 1

    # 配体来源信息
    L_slope <- strict$trend_slope[match(lig_subs[1], strict$gene_symbol)]
    L_FDR   <- strict$trend_FDR  [match(lig_subs[1], strict$gene_symbol)]
    L_class <- strict$func_class [match(lig_subs[1], strict$gene_symbol)]

    # CellChat Hill 结合概率
    P <- hill(lig_rn) * hill(rec_rn)

    out[[length(out) + 1L]] <- data.frame(
      L_label = label,
      ligand_complex   = ia_tbl$ligand_complex[i],
      ligand_subunits  = ia_tbl$ligand_subunits[i],
      receptor_complex = ia_tbl$receptor_complex[i],
      receptor_subunits = ia_tbl$receptor_subunits[i],
      n_recep_subunits = ia_tbl$n_recep_subunits[i],
      n_up_subunits    = n_up_sub,
      all_subunits_up  = all_up,
      any_subunit_up   = any_up,
      r_dir_up         = r_dir_up,
      pathway_name = ia_tbl$pathway_name[i],
      annotation   = ia_tbl$annotation[i],
      interaction_name = ia_tbl$interaction_name[i],
      L_expr_Mo_TX91 = lig_raw, L_rn = lig_rn,
      R_expr_vagus_IAV = rec_raw, R_rn = rec_rn,
      L_trend_slope = L_slope, L_trend_FDR = L_FDR, L_func_class = L_class,
      R_max_log2FC = r_max_lfc, R_min_FDR = r_min_fdr,
      prob = P,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  d <- do.call(rbind, out)
  d <- d[order(-d$prob), ]
  rownames(d) <- NULL
  d
}

# 计算: 严判 (Secreted Signaling) + 放宽 (+ Cell-Cell Contact)
cat("\n>>> 计算 LR 表 (R 在迷走神经表达即纳入, 标注上调与显著性)\n")
lr_58_ss <- compute_lr_table(ia_58_ss, "58_SS")
lr_58_cc <- compute_lr_table(ia_58_cc, "58_SS+CCC")
lr_4_ss  <- compute_lr_table(ia_4_ss,  "4_SS")
lr_4_cc  <- compute_lr_table(ia_4_cc,  "4_SS+CCC")
cat(sprintf("    58 SS:        %d  (其中 R 全亚基上调 = %d, R 任一亚基上调 = %d, R 上调方向 = %d)\n",
            nrow(lr_58_ss),
            if (nrow(lr_58_ss)) sum(lr_58_ss$all_subunits_up) else 0,
            if (nrow(lr_58_ss)) sum(lr_58_ss$any_subunit_up)  else 0,
            if (nrow(lr_58_ss)) sum(lr_58_ss$r_dir_up)        else 0))
cat(sprintf("    58 SS+CCC:    %d  (其中 R 全亚基上调 = %d, R 任一亚基上调 = %d, R 上调方向 = %d)\n",
            nrow(lr_58_cc),
            if (nrow(lr_58_cc)) sum(lr_58_cc$all_subunits_up) else 0,
            if (nrow(lr_58_cc)) sum(lr_58_cc$any_subunit_up)  else 0,
            if (nrow(lr_58_cc)) sum(lr_58_cc$r_dir_up)        else 0))
cat(sprintf("    4  SS:        %d\n", nrow(lr_4_ss)))
cat(sprintf("    4  SS+CCC:    %d\n", nrow(lr_4_cc)))

# 全表输出 (主用 SS)
write.table(lr_58_ss, file.path(OUT, "LR_from_58_SecretedSignaling.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr_58_cc, file.path(OUT, "LR_from_58_SS_plus_CellCellContact.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr_4_ss,  file.path(OUT, "LR_from_4both_SecretedSignaling.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr_4_cc,  file.path(OUT, "LR_from_4both_SS_plus_CellCellContact.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(lr_58_ss,   file.path(OUT, "LR_from_58_SecretedSignaling.csv"),
          quote = TRUE, row.names = FALSE)

# 高可信子集 (主选: 严判 SS + R 全亚基上调 FDR<0.05)
filter_high <- function(d) {
  if (nrow(d) == 0) return(d)
  d[d$all_subunits_up & d$prob >= PARAMS$prob_high_thresh, ]
}
filter_mid <- function(d) {
  if (nrow(d) == 0) return(d)
  d[(d$any_subunit_up | d$r_dir_up) & d$prob >= PARAMS$prob_high_thresh, ]
}
high_58 <- filter_high(lr_58_ss)
high_4  <- filter_high(lr_4_ss)
mid_58  <- filter_mid(lr_58_ss)

write.table(high_58, file.path(OUT, "LR_high_confidence_from_58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(high_4,  file.path(OUT, "LR_high_confidence_from_4both.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mid_58, file.path(OUT, "LR_mid_confidence_from_58.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n>>> 高可信 (SS + R 全亚基上调 + prob>=%.2f): from 58 = %d, from 4 = %d\n",
            PARAMS$prob_high_thresh, nrow(high_58), nrow(high_4)))
cat(sprintf(">>> 中可信 (SS + R 任一亚基上调或方向上调 + prob>=%.2f): from 58 = %d\n",
            PARAMS$prob_high_thresh, nrow(mid_58)))

# 5. 流程汇总 ----------------------------------------------------------------
sink(file.path(OUT, "CellChat_summary.txt"))
cat("=== Mo -> vagal sensory ganglia LR 通讯 (CellChat 风格) ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/10_CellChat_Mo_to_vagus.R\n\n")
cat("配体源: GSE42639 Mo 严选 58 (Mo_protective/strict/core_targets_strict.tsv)\n")
cat(sprintf("受体源: GSE161878 主分析 (剔 D-suffix) 上调 DEG %d (|log2FC|>=0.585 & FDR<0.05)\n",
            length(up_sym)))
cat(sprintf("数据库: CellChatDB.mouse  SS=%d, SS+CCC=%d (LR pairs)\n\n",
            nrow(ia_ss), nrow(ia_ss_cc)))
cat("参数:\n"); print(PARAMS)
cat("\n=== 58 候选 (sender), Secreted Signaling 严判 ============\n")
cat(sprintf("LR pair (R 在迷走神经检测到):       %d\n", nrow(lr_58_ss)))
cat(sprintf("  其中 R 全亚基显著上调 (FDR<0.05):  %d\n",
            if (nrow(lr_58_ss)) sum(lr_58_ss$all_subunits_up) else 0))
cat(sprintf("  其中 R 任一亚基显著上调:           %d\n",
            if (nrow(lr_58_ss)) sum(lr_58_ss$any_subunit_up) else 0))
cat(sprintf("  其中 R 表达方向上调 (log2FC>0):    %d\n",
            if (nrow(lr_58_ss)) sum(lr_58_ss$r_dir_up) else 0))
cat(sprintf("高可信 (全亚基上调 + prob>=%.2f):    %d\n",
            PARAMS$prob_high_thresh, nrow(high_58)))
cat(sprintf("中可信 (任一上调或方向上调 + prob):   %d\n", nrow(mid_58)))

cat("\n--- 58 SS 全部 LR (按 prob 降序) ---\n")
if (nrow(lr_58_ss) > 0) {
  print(lr_58_ss[, c("ligand_subunits","receptor_subunits",
                     "pathway_name","prob",
                     "all_subunits_up","any_subunit_up","r_dir_up",
                     "R_max_log2FC","R_min_FDR",
                     "L_trend_slope","L_func_class")],
        row.names = FALSE)
}

cat("\n=== 58 候选 (sender), 放宽: + Cell-Cell Contact ==========\n")
cat(sprintf("LR pair: %d (其中 R 显著上调全亚基: %d)\n",
            nrow(lr_58_cc),
            if (nrow(lr_58_cc)) sum(lr_58_cc$all_subunits_up) else 0))
if (nrow(lr_58_cc) > 0) {
  print(lr_58_cc[, c("annotation","ligand_subunits","receptor_subunits",
                     "pathway_name","prob",
                     "all_subunits_up","any_subunit_up","r_dir_up",
                     "R_max_log2FC","R_min_FDR")],
        row.names = FALSE)
}

cat("\n=== 4 双标 (sender) ======================================\n")
cat(sprintf("SS LR: %d   ;   SS+CCC LR: %d\n", nrow(lr_4_ss), nrow(lr_4_cc)))
if (nrow(lr_4_cc) > 0) {
  print(lr_4_cc[, c("annotation","ligand_subunits","receptor_subunits",
                    "pathway_name","prob","all_subunits_up","any_subunit_up",
                    "r_dir_up","R_max_log2FC","R_min_FDR")],
        row.names = FALSE)
}
sink()
cat(readLines(file.path(OUT, "CellChat_summary.txt")), sep = "\n")

# 6. 可视化 ------------------------------------------------------------------
# 6a. 高可信 LR 网络 (二部图: ligand -> receptor)
draw_network <- function(d, out_prefix, title_main) {
  if (nrow(d) == 0) return(invisible(NULL))
  nodes_l <- unique(d$ligand_complex)
  nodes_r <- unique(d$receptor_complex)
  nl <- length(nodes_l); nr <- length(nodes_r)
  yl <- if (nl > 1) seq(1, 0, length.out = nl) else 0.5
  yr <- if (nr > 1) seq(1, 0, length.out = nr) else 0.5
  pos_l <- data.frame(name = nodes_l, x = 0, y = yl, side = "L (Mo)",
                      stringsAsFactors = FALSE)
  pos_r <- data.frame(name = nodes_r, x = 1, y = yr, side = "R (vagus)",
                      stringsAsFactors = FALSE)
  pos <- rbind(pos_l, pos_r)
  edges <- data.frame(
    x1 = pos_l$x[match(d$ligand_complex, pos_l$name)],
    y1 = pos_l$y[match(d$ligand_complex, pos_l$name)],
    x2 = pos_r$x[match(d$receptor_complex, pos_r$name)],
    y2 = pos_r$y[match(d$receptor_complex, pos_r$name)],
    prob = d$prob,
    pathway = d$pathway_name,
    stringsAsFactors = FALSE
  )
  p <- ggplot() +
    geom_segment(data = edges,
                 aes(x = x1, y = y1, xend = x2, yend = y2,
                     alpha = prob, linewidth = prob),
                 color = "#377EB8") +
    geom_point(data = pos, aes(x, y, color = side), size = 3) +
    geom_text(data = pos_l, aes(x, y, label = name),
              hjust = 1.1, size = 3) +
    geom_text(data = pos_r, aes(x, y, label = name),
              hjust = -0.1, size = 3) +
    scale_color_manual(values = c("L (Mo)" = "#E41A1C",
                                  "R (vagus)" = "#4DAF4A")) +
    scale_linewidth_continuous(range = c(0.3, 1.8)) +
    scale_alpha_continuous(range = c(0.3, 1)) +
    coord_cartesian(xlim = c(-0.4, 1.4), clip = "off") +
    labs(title = title_main,
         subtitle = sprintf("%d LR pairs, 边粗 = prob (Hill K=%g, n=%d)",
                            nrow(d), PARAMS$K, PARAMS$n_hill)) +
    theme_void(base_size = 10) +
    theme(plot.margin = margin(15, 80, 15, 80))
  h <- max(4, max(nl, nr) * 0.25 + 2)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = h,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = h,
         dpi = 200, type = "cairo", limitsize = FALSE)
}

# 主可视: 严判 SS 全部 LR 的二部网络 (按 prob 着色)
draw_network(lr_58_ss, file.path(OUT, "network_58_SS"),
             "Mo 58 SS 配体 -> 迷走神经受体 (按 prob)")
# 放宽 SS+CCC
draw_network(lr_58_cc, file.path(OUT, "network_58_SS_plus_CCC"),
             "Mo 58 SS+CCC 配体 -> 迷走神经受体 (按 prob)")
# 4 双标
draw_network(lr_4_cc,  file.path(OUT, "network_4both_SS_plus_CCC"),
             "Mo 4 双标配体 -> 迷走神经受体 (SS+CCC)")
cat(">>> 网络图: network_*.pdf/png\n")

# 通路条形图
plot_pathway_bar <- function(d, out_prefix, title_main) {
  if (nrow(d) == 0) return(invisible(NULL))
  pw <- aggregate(prob ~ pathway_name, data = d, FUN = max)
  pw <- pw[order(-pw$prob), ]
  top_pw <- head(pw, 20)
  top_pw$pathway_name <- factor(top_pw$pathway_name, levels = rev(top_pw$pathway_name))
  p_pw <- ggplot(top_pw, aes(prob, pathway_name)) +
    geom_col(fill = "#377EB8") +
    geom_text(aes(label = sprintf("%.3f", prob)), hjust = -0.1, size = 3) +
    expand_limits(x = max(top_pw$prob) * 1.15) +
    labs(x = "max prob (within pathway)", y = NULL,
         title = title_main) +
    theme_bw(base_size = 10)
  ggsave(paste0(out_prefix, ".pdf"), p_pw, width = 8, height = 6, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p_pw, width = 8, height = 6,
         dpi = 200, type = "cairo")
}
plot_pathway_bar(lr_58_cc, file.path(OUT, "pathway_topprob_from_58_SSplusCCC"),
                 "Top 20 通路 (Mo 58 SS+CCC)")

saveRDS(list(ia_58_ss = ia_58_ss, ia_58_cc = ia_58_cc,
             ia_4_ss = ia_4_ss,  ia_4_cc = ia_4_cc,
             lr_58_ss = lr_58_ss, lr_58_cc = lr_58_cc,
             lr_4_ss = lr_4_ss, lr_4_cc = lr_4_cc,
             high_58 = high_58, mid_58 = mid_58, high_4 = high_4,
             PARAMS = PARAMS,
             mo_tx91_mean = mo_tx91_mean, mo_rn = mo_rn,
             vagus_iav_mean = vagus_iav_mean, vagus_rn = vagus_rn),
        file.path(OUT, "CellChat_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
