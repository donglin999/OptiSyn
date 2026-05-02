#!/usr/bin/env Rscript
# CellChat 风格 LR 配对 v2: 神经元区室受体池
#
# 配体 L: GSE42639 Mo 严选 58 ∩ CellChatDB SS (+CCC 放宽) ligand
# 受体 R: 重新构建 (新)
#   1. GSE161878 (bulk vagal ganglia) CPM > 1 在 >= 20% 样本 -> "vagal 表达池"
#   2. 排除免疫 marker 黑名单 (经典 + GSE296065 单细胞免疫细胞特异 marker)
#   3. ∩ CellChatDB receptor (单基因 + complex 任一亚基)
#   4. (附) 神经元 marker enrichment: 用 vagal nociceptor / sensory neuron classical marker
#      验证保留下来的基因池里神经元成分 (信息列, 不强制)
# 配对:
#   遍历 SS+CCC LR pairs, L 在 58, R 单基因或 complex 任一亚基在受体池
# 概率:
#   Hill: P = H(E_L) * H(E_R), K=0.5, n=2; 表达用 rank-norm
# 趋势匹配标签:
#   L_dir = -1 (Mo TX91->100 下调, 严选 trend_slope<0)
#   R_dir_IAV_vs_Mock: + (logFC>0) 或 - (logFC<0)
#   trend_match:
#     "co_down"  : L 下调 + R 方向下调 (保护轴共衰退)
#     "mirror_up": L 下调 + R 方向上调 (信号丧失 -> 代偿上调)
#     "no_trend" : R logFC=NA (低表达过滤外, 无方向信息)
#
# 输出: results/CellChat/v2_neuron_compartment/

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(Seurat)
  library(dplyr)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/CellChat/v2_neuron_compartment")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(K = 0.5, n_hill = 2,
               cpm_thresh = 1,
               sample_frac_min = 0.20,   # >= 20% 样本里 CPM>thresh
               prob_high_thresh = 0.10)

# ============================================================================
# 1. 加载 CellChatDB
# ============================================================================
load(file.path(ROOT, "processed_data/CellChatDB/CellChatDB.mouse.rda"))
ia    <- CellChatDB.mouse$interaction
cmplx <- CellChatDB.mouse$complex
ia_ss <- ia[ia$annotation == "Secreted Signaling", , drop = FALSE]
ia_ss_cc <- ia[ia$annotation %in% c("Secreted Signaling","Cell-Cell Contact"), , drop = FALSE]
cat(sprintf(">>> CellChatDB.mouse: SS=%d, SS+CCC=%d\n",
            nrow(ia_ss), nrow(ia_ss_cc)))

# ============================================================================
# 2. 加载 GSE42639 Mo (配体表达)
# ============================================================================
mo_res <- readRDS(file.path(ROOT, "results/Mo_protective/Mo_protective_results.rds"))
mo_E <- mo_res$expr
mo_m <- mo_res$meta
tx91_idx <- which(mo_m$treat == "TX91")
mo_tx91_mean <- rowMeans(mo_E[, tx91_idx, drop = FALSE])
strict <- read.delim(file.path(ROOT,
            "results/Mo_protective/strict/core_targets_strict.tsv"),
            stringsAsFactors = FALSE)
genes58 <- strict$gene_symbol
cat(sprintf(">>> Mo TX91 mean: %d genes\n", length(mo_tx91_mean)))

# ============================================================================
# 3. 加载 GSE161878 raw counts -> CPM, 计算 vagal 表达池
# ============================================================================
cnt_path <- file.path(ROOT, "raw_data/00_raw_data/GSE161878_gene_counts.txt")
meta_path <- file.path(ROOT, "raw_data/00_raw_data/GSE161878_sample_metadata.txt")
cnt_raw <- read.delim(cnt_path, check.names = FALSE, stringsAsFactors = FALSE)
vst_m <- read.delim(meta_path, stringsAsFactors = FALSE)
genes_e <- as.character(cnt_raw$GeneID)
M_raw <- as.matrix(cnt_raw[, -1, drop = FALSE]); rownames(M_raw) <- genes_e
vst_m <- vst_m[match(colnames(M_raw), vst_m$sample_id), , drop = FALSE]

# 主分析: 剔 D-suffix
keep_main <- !grepl("D$", vst_m$sample_id)
M_raw <- M_raw[, keep_main, drop = FALSE]
vst_m <- vst_m[keep_main, , drop = FALSE]

lib_size <- colSums(M_raw)
cpm <- t(t(M_raw) / lib_size) * 1e6
log2cpm <- log2(cpm + 1)

ann_g <- read.delim(file.path(ROOT, "processed_data/GSE161878/gene_annotation.tsv"),
                    stringsAsFactors = FALSE)
sym_lookup <- setNames(ann_g$SYMBOL, ann_g$ENTREZID)
g_sym <- sym_lookup[rownames(log2cpm)]
keep_v <- !is.na(g_sym) & nchar(g_sym) > 0
log2cpm <- log2cpm[keep_v, , drop = FALSE]; cpm <- cpm[keep_v, , drop = FALSE]
g_sym <- g_sym[keep_v]
# 同 SYMBOL 多 ENTREZ: 取 mean log2cpm 最大的代表
ord <- order(g_sym, -rowMeans(log2cpm))
log2cpm <- log2cpm[ord, ]; cpm <- cpm[ord, ]; g_sym <- g_sym[ord]
keep_first <- !duplicated(g_sym)
log2cpm <- log2cpm[keep_first, ]; cpm <- cpm[keep_first, ]
rownames(log2cpm) <- g_sym[keep_first]; rownames(cpm) <- g_sym[keep_first]

# 表达阈值: CPM > 1 (即 log2cpm > log2(2))
n_samples <- ncol(log2cpm)
n_min <- ceiling(PARAMS$sample_frac_min * n_samples)
detected <- rowSums(cpm > PARAMS$cpm_thresh)
expressed_pool <- rownames(log2cpm)[detected >= n_min]
cat(sprintf(">>> GSE161878 主分析 %d 样本 (剔 D-suffix)\n", n_samples))
cat(sprintf(">>> vagal 表达池 (CPM>%g 在 >=%d/%d 样本): %d 基因\n",
            PARAMS$cpm_thresh, n_min, n_samples, length(expressed_pool)))

# 对每个基因取所有样本的 IAV/Mock 平均表达 (用于受体表达分数)
iav_idx <- which(vst_m$group == "IAV")
mock_idx <- which(vst_m$group == "Mock")
vagus_iav_mean  <- rowMeans(log2cpm[, iav_idx,  drop = FALSE])
vagus_mock_mean <- rowMeans(log2cpm[, mock_idx, drop = FALSE])

# ============================================================================
# 4. 免疫细胞 marker 黑名单
# ============================================================================
# 4a. 经典免疫 marker (手动 curate, 鼠 SYMBOL)
immune_classic <- c(
  # Pan-leukocyte
  "Ptprc",
  # Myeloid (Mo / Neu / DC / Macroph)
  "Itgam","Itgax","Itgb2","Itgal","Cd68","Cd14","Csf1r","Adgre1","Mertk",
  "Ly6c1","Ly6c2","Ly6g","S100a8","S100a9","Mpo","Elane","Ngp","Camp",
  "Fcgr1","Fcgr2b","Fcgr3","Fcgr4","Fcer1g","Tyrobp","Fpr1","Fpr2",
  "Ccr2","Ccr5","Cx3cr1","Cd163","Cd86","Cd80","H2-Aa","H2-Ab1","H2-Eb1",
  "Cd74","Ctss","Lyz1","Lyz2","Marco","Ifitm1","Ifitm3",
  # T cells
  "Cd3d","Cd3e","Cd3g","Cd4","Cd8a","Cd8b1","Foxp3","Lat","Lck","Trac","Trbc1","Trbc2",
  # B cells
  "Cd19","Ms4a1","Cd79a","Cd79b","Ighm","Cd22","Pax5",
  # NK / ILC
  "Klrb1c","Ncr1","Klrk1","Gzma","Gzmb","Prf1",
  # Mast / basophil / eosinophil
  "Kit","Fcer1a","Mcpt1","Mcpt8",
  # Complement & immune-restricted
  "C1qa","C1qb","C1qc","C3","C3ar1","C5ar1",
  # Mo/Mφ markers
  "Plac8","Lgals3","Ms4a6c","Ms4a7","Spi1"
)
# 4b. GSE296065 单细胞: 在免疫细胞中高度特异的 marker
# (从 296065 验证脚本结果衍生; 这里轻量做一遍 FindAllMarkers 取 top markers per non-neuronal cluster)
sc_pp_path <- file.path(ROOT, "processed_data/GSE296065_immune_markers.rds")
if (!file.exists(sc_pp_path)) {
  cat(">>> 计算 GSE296065 单细胞免疫 marker (one-shot, 缓存)\n")
  rds_sc <- file.path(ROOT, "raw_data/00_raw_data/GSE296065_almanzar_flu.rds")
  seu <- readRDS(rds_sc)
  Idents(seu) <- "celltype"
  am <- tryCatch({
    FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 1,
                   min.pct = 0.25, return.thresh = 1e-50,
                   verbose = FALSE)
  }, error = function(e) {
    cat("  FindAllMarkers 出错: ", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(am)) saveRDS(am, sc_pp_path) else saveRDS(data.frame(), sc_pp_path)
  rm(seu); gc()
}
sc_immune_full <- readRDS(sc_pp_path)
# 收严: 单细胞 marker 必须 avg_log2FC > 2 & p_val_adj < 1e-100 (高强度特异)
if (nrow(sc_immune_full) > 0) {
  sc_strict <- sc_immune_full[sc_immune_full$avg_log2FC > 2 &
                              sc_immune_full$p_val_adj < 1e-100, ]
  immune_sc <- unique(sc_strict$gene)
  cat(sprintf("    单细胞强免疫 marker (avg_log2FC>2 & p_adj<1e-100): %d\n",
              length(immune_sc)))
  cat(sprintf("    (原 avg_log2FC>1 & p_adj<1e-50 共 %d 个, 阈值过松剔除)\n",
              length(unique(sc_immune_full$gene))))
} else {
  immune_sc <- character(0)
}

immune_blacklist <- unique(c(immune_classic, immune_sc))
cat(sprintf(">>> 免疫黑名单合计: %d 基因 (经典 %d + 单细胞 %d)\n",
            length(immune_blacklist),
            length(immune_classic),
            length(setdiff(immune_sc, immune_classic))))

# ============================================================================
# 5. 神经元正向 marker (信息列, 验证表达池含神经元成分)
# ============================================================================
neuron_markers <- c(
  # Pan-neuronal / structural
  "Snap25","Stmn2","Stmn3","Map2","Tubb3","Nefl","Nefm","Nefh","Gap43",
  "Syn1","Syp","Syt1","Rbfox3","Eno2","Uchl1","Mapt",
  # Vagal / DRG sensory neuron
  "Gpr151","Trpv1","Trpa1","Scn10a","Scn11a","Piezo2","P2rx3","Tac1","Calca",
  "Mrgprd","Trpm8","Ntrk1","Ntrk2","Ntrk3","Ret","Gfra2","Gfra3",
  # Interoceptive / mechano vagal
  "Vipr2","Calcb","Sst","Npy","Cgrp"
)
neuron_in_pool <- intersect(neuron_markers, expressed_pool)
cat(sprintf(">>> 神经元 marker (panel %d) ∩ 表达池: %d -> %s\n",
            length(neuron_markers), length(neuron_in_pool),
            paste(head(neuron_in_pool, 20), collapse=", ")))

# ============================================================================
# 6. 构建受体池 R (神经元相关)
# ============================================================================
# CellChatDB receptor 单基因 + complex 全亚基
expand_x <- function(x) {
  if (x %in% rownames(cmplx)) {
    s <- as.character(cmplx[x, ]); s[nchar(s) > 0]
  } else x
}
all_rec_in_db <- unique(unlist(lapply(unique(ia$receptor), expand_x)))

# 神经元相关受体池 (两版)
R_pool_pre        <- intersect(expressed_pool, all_rec_in_db)
R_pool_neu_strict <- setdiff(R_pool_pre, immune_blacklist)         # 经典 + 单细胞
R_pool_neu_relax  <- setdiff(R_pool_pre, immune_classic)            # 仅经典 (放宽)
R_pool_neu        <- R_pool_neu_strict   # 默认严判
cat(sprintf("\n>>> 受体池构建:\n"))
cat(sprintf("    vagal 表达池                          : %d\n", length(expressed_pool)))
cat(sprintf("    ∩ CellChatDB receptor                 : %d\n", length(R_pool_pre)))
cat(sprintf("    严判 (剔经典+单细胞免疫 marker)        : %d\n", length(R_pool_neu_strict)))
cat(sprintf("    放宽 (仅剔经典免疫 marker)             : %d\n", length(R_pool_neu_relax)))

# === 诊断: 3 个候选受体 (Tnfrsf14/Cd4/Ceacam1) 命运追踪 ===
cat("\n>>> 诊断: 3 候选受体在各阶段命运\n")
diag_recs <- c("Tnfrsf14","Cd4","Ceacam1")
for (g in diag_recs) {
  cpm_v <- if (g %in% rownames(cpm)) cpm[g,] else NA
  log2_v <- if (g %in% rownames(log2cpm)) log2cpm[g,] else NA
  in_pool   <- g %in% expressed_pool
  in_db     <- g %in% all_rec_in_db
  in_classic<- g %in% immune_classic
  in_sc     <- g %in% immune_sc
  in_neuR   <- g %in% R_pool_neu
  cat(sprintf("  [%s] CPM range %s, n_samples_CPM>%g=%d, in_vagal_pool=%s, in_CCDB_rec=%s, in_classic_BL=%s, in_sc_BL=%s, in_neuron_R_pool=%s\n",
              g,
              if (any(!is.na(cpm_v))) sprintf("%.2f-%.2f", min(cpm_v, na.rm=TRUE), max(cpm_v, na.rm=TRUE)) else "NA",
              PARAMS$cpm_thresh,
              if (any(!is.na(cpm_v))) sum(cpm_v > PARAMS$cpm_thresh, na.rm=TRUE) else 0,
              in_pool, in_db, in_classic, in_sc, in_neuR))
}

# 保存被剔除的"免疫受体"列表 (信息)
removed_immune_R <- intersect(R_pool_pre, immune_blacklist)
write.table(data.frame(gene = removed_immune_R,
                       in_classic_blacklist = removed_immune_R %in% immune_classic,
                       in_sc_blacklist      = removed_immune_R %in% immune_sc),
            file.path(OUT, "removed_immune_receptors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# 7. 配体 L
# ============================================================================
L_in_db_ss   <- intersect(genes58, unique(unlist(lapply(unique(ia_ss$ligand), expand_x))))
L_in_db_cc   <- intersect(genes58, unique(unlist(lapply(unique(ia_ss_cc$ligand), expand_x))))
cat(sprintf("\n>>> 配体 (58 ∩ CellChatDB):\n"))
cat(sprintf("    SS:        %d -> %s\n", length(L_in_db_ss),  paste(L_in_db_ss, collapse=", ")))
cat(sprintf("    SS+CCC:    %d -> %s\n", length(L_in_db_cc),  paste(L_in_db_cc, collapse=", ")))

# ============================================================================
# 8. 上调 DEG 信息 (用于 R 方向标签, 不作过滤)
# ============================================================================
deg_full <- read.delim(file.path(ROOT, "results/GSE161878/GSE161878_DEA_full.tsv"),
                       stringsAsFactors = FALSE)
deg_full <- deg_full[!is.na(deg_full$SYMBOL) & nchar(deg_full$SYMBOL) > 0, ]

# ============================================================================
# 9. LR 配对 + Hill prob + 趋势标签
# ============================================================================
mo_rn <- rank(mo_tx91_mean, ties.method = "average") / length(mo_tx91_mean)
names(mo_rn) <- names(mo_tx91_mean)
vagus_rn <- rank(vagus_iav_mean, ties.method = "average") / length(vagus_iav_mean)
names(vagus_rn) <- names(vagus_iav_mean)

hill <- function(E, K = PARAMS$K, n = PARAMS$n_hill) E^n / (K^n + E^n)
geo_mean <- function(v) {
  v <- v[!is.na(v) & v > 0]
  if (length(v) == 0) return(NA_real_)
  exp(mean(log(v)))
}

build_lr <- function(ia_pool, label) {
  rows <- list()
  for (i in seq_len(nrow(ia_pool))) {
    lig <- ia_pool$ligand[i]
    rec <- ia_pool$receptor[i]
    lig_subs <- expand_x(lig)
    rec_subs <- expand_x(rec)

    # L 必须在 58 候选 (任一亚基)
    if (!any(lig_subs %in% genes58)) next
    # L 必须在 Mo 表达
    if (!all(lig_subs %in% names(mo_tx91_mean))) next
    # R 至少一个亚基在神经元受体池
    rec_in_neu <- rec_subs %in% R_pool_neu
    if (!any(rec_in_neu)) next

    # 表达计算
    L_raw <- geo_mean(mo_tx91_mean[lig_subs])
    L_rn  <- geo_mean(mo_rn       [lig_subs])
    rec_subs_have <- rec_subs[rec_subs %in% names(vagus_iav_mean)]
    if (length(rec_subs_have) == 0) next
    R_raw <- geo_mean(vagus_iav_mean [rec_subs_have])
    R_mock_raw <- geo_mean(vagus_mock_mean[rec_subs_have])
    R_rn  <- geo_mean(vagus_rn      [rec_subs_have])

    P <- hill(L_rn) * hill(R_rn)

    # 受体方向 (任一亚基的 max log2FC)
    ds <- deg_full[deg_full$SYMBOL %in% rec_subs, , drop = FALSE]
    if (nrow(ds) > 0) {
      idx_max <- which.max(abs(ds$log2FC_shrink))
      r_lfc <- ds$log2FC_shrink[idx_max]
      r_fdr <- ds$FDR         [idx_max]
    } else {
      r_lfc <- NA_real_; r_fdr <- NA_real_
    }

    # 趋势匹配
    L_slope <- strict$trend_slope[match(lig_subs[1], strict$gene_symbol)]
    L_dir   <- "down"  # 58 严选都是 slope<0
    if (is.na(r_lfc)) {
      trend <- "no_trend"
      R_dir <- NA
    } else if (r_lfc > 0) {
      trend <- "mirror_up"
      R_dir <- "up"
    } else {
      trend <- "co_down"
      R_dir <- "down"
    }

    # 神经元 marker 同表达池存在性
    n_rec_in_neu_pool   <- sum(rec_in_neu)
    n_rec_subs_total    <- length(rec_subs)
    rec_pool_complete   <- n_rec_in_neu_pool == n_rec_subs_total

    rows[[length(rows) + 1L]] <- data.frame(
      L_label = label,
      ligand_complex = lig,
      ligand_subunits = paste(lig_subs, collapse = ";"),
      receptor_complex = rec,
      receptor_subunits = paste(rec_subs, collapse = ";"),
      n_rec_subunits = n_rec_subs_total,
      n_rec_in_neuron_pool = n_rec_in_neu_pool,
      rec_pool_complete = rec_pool_complete,
      pathway_name = ia_pool$pathway_name[i],
      annotation   = ia_pool$annotation[i],
      interaction_name = ia_pool$interaction_name[i],
      L_expr_Mo_TX91   = L_raw, L_rn = L_rn,
      R_expr_vagus_IAV = R_raw, R_expr_vagus_Mock = R_mock_raw, R_rn = R_rn,
      L_trend_slope_GSE42639 = L_slope,
      L_dir = L_dir,
      R_log2FC_IAV_vs_Mock = r_lfc, R_FDR = r_fdr, R_dir = R_dir,
      trend_match = trend,
      prob = P,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  d <- do.call(rbind, rows)
  d <- d[order(-d$prob), ]
  rownames(d) <- NULL
  d
}

cat("\n>>> LR 配对 - 严判 R 池 (经典+单细胞免疫黑名单)\n")
R_pool_neu <- R_pool_neu_strict
lr_ss <- build_lr(ia_ss, "SS_strict")
lr_cc <- build_lr(ia_ss_cc, "SS+CCC_strict")
cat(sprintf("    SS LR pairs:     %d\n", nrow(lr_ss)))
cat(sprintf("    SS+CCC LR pairs: %d\n", nrow(lr_cc)))

cat("\n>>> LR 配对 - 放宽 R 池 (仅经典免疫黑名单, 保留双重身份基因如 Ceacam1)\n")
R_pool_neu <- R_pool_neu_relax
lr_ss_relax <- build_lr(ia_ss, "SS_relax")
lr_cc_relax <- build_lr(ia_ss_cc, "SS+CCC_relax")
cat(sprintf("    SS LR pairs (放宽):     %d\n", nrow(lr_ss_relax)))
cat(sprintf("    SS+CCC LR pairs (放宽): %d\n", nrow(lr_cc_relax)))

# 默认主结果用 strict, 放宽作附属
write.table(lr_ss_relax, file.path(OUT, "LR_v2_SecretedSignaling_relax.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr_cc_relax, file.path(OUT, "LR_v2_SS_plus_CellCellContact_relax.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 写表
write.table(lr_ss, file.path(OUT, "LR_v2_SecretedSignaling.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr_cc, file.path(OUT, "LR_v2_SS_plus_CellCellContact.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(lr_cc, file.path(OUT, "LR_v2_SS_plus_CellCellContact.csv"),
          quote = TRUE, row.names = FALSE)

# 趋势匹配子集
filter_match <- function(d, mode) {
  if (nrow(d) == 0) return(d)
  d[d$trend_match == mode & d$prob >= PARAMS$prob_high_thresh, ]
}
match_codown_ss   <- filter_match(lr_ss, "co_down")
match_mirror_ss   <- filter_match(lr_ss, "mirror_up")
match_codown_cc   <- filter_match(lr_cc, "co_down")
match_mirror_cc   <- filter_match(lr_cc, "mirror_up")
write.table(match_codown_ss, file.path(OUT, "match_co_down_SS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(match_mirror_ss, file.path(OUT, "match_mirror_up_SS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(match_codown_cc, file.path(OUT, "match_co_down_SS_CCC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(match_mirror_cc, file.path(OUT, "match_mirror_up_SS_CCC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n>>> 趋势匹配 (prob>=%.2f):\n", PARAMS$prob_high_thresh))
cat(sprintf("    SS     co_down (L↓ R↓):     %d\n", nrow(match_codown_ss)))
cat(sprintf("    SS     mirror_up (L↓ R↑):   %d\n", nrow(match_mirror_ss)))
cat(sprintf("    SS+CCC co_down:              %d\n", nrow(match_codown_cc)))
cat(sprintf("    SS+CCC mirror_up:            %d\n", nrow(match_mirror_cc)))

# ============================================================================
# 10. 流程汇总
# ============================================================================
sink(file.path(OUT, "summary.txt"))
cat("=== Mo -> vagal 神经元区室 LR 通讯 v2 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("脚本: scripts/12_CellChat_neuron_compartment_LR.R\n\n")

cat("配体源: GSE42639 Mo 严选 58 (results/Mo_protective/strict/core_targets_strict.tsv)\n")
cat("受体源: GSE161878 (bulk vagal ganglia) 神经元区室池 (新)\n")
cat("数据库: CellChatDB.mouse SS / SS+Cell-Cell-Contact\n")

cat(sprintf("\n参数:\n"))
print(PARAMS)

cat("\n=== 受体池构建 ===\n")
cat(sprintf("vagal 表达池 (CPM>%g 在 >=%d/%d 样本): %d\n",
            PARAMS$cpm_thresh, n_min, n_samples, length(expressed_pool)))
cat(sprintf("∩ CellChatDB receptor: %d\n", length(R_pool_pre)))
cat(sprintf("剔除经典+单细胞免疫黑名单 (%d): -%d\n",
            length(immune_blacklist), length(removed_immune_R)))
cat(sprintf("最终神经元相关受体池: %d\n", length(R_pool_neu)))
cat(sprintf("神经元 marker 验证 (panel %d): %d 在表达池\n",
            length(neuron_markers), length(neuron_in_pool)))
cat(sprintf("  含: %s\n", paste(head(neuron_in_pool, 20), collapse=", ")))

cat("\n=== 配体 ===\n")
cat(sprintf("58 ∩ SS ligand:      %d -> %s\n", length(L_in_db_ss),
            paste(L_in_db_ss, collapse=", ")))
cat(sprintf("58 ∩ SS+CCC ligand:  %d -> %s\n", length(L_in_db_cc),
            paste(L_in_db_cc, collapse=", ")))

cat("\n=== LR 配对结果 ===\n")
cat(sprintf("SS LR pairs (R 至少一亚基在神经元池): %d\n", nrow(lr_ss)))
cat(sprintf("SS+CCC LR pairs:                      %d\n", nrow(lr_cc)))

cat("\n=== 高可信通讯对 (prob>=0.10 + 趋势匹配) ===\n")
cat(sprintf("%-30s %s\n", "类别", "数"))
cat(sprintf("%-30s %d\n", "SS co_down (L↓ R↓)",     nrow(match_codown_ss)))
cat(sprintf("%-30s %d\n", "SS mirror_up (L↓ R↑)",   nrow(match_mirror_ss)))
cat(sprintf("%-30s %d\n", "SS+CCC co_down",          nrow(match_codown_cc)))
cat(sprintf("%-30s %d\n", "SS+CCC mirror_up",        nrow(match_mirror_cc)))

cat("\n--- SS 全部 LR (按 prob 降序) ---\n")
if (nrow(lr_ss) > 0) {
  print(lr_ss[, c("ligand_subunits","receptor_subunits","pathway_name",
                  "prob","R_log2FC_IAV_vs_Mock","R_FDR",
                  "trend_match","R_expr_vagus_IAV")],
        row.names = FALSE)
}
cat("\n--- SS+CCC 全部 LR (严判, 按 prob 降序) ---\n")
if (nrow(lr_cc) > 0) {
  print(lr_cc[, c("annotation","ligand_subunits","receptor_subunits",
                  "pathway_name","prob","R_log2FC_IAV_vs_Mock","R_FDR",
                  "trend_match")],
        row.names = FALSE)
}

cat("\n=== 放宽版 R 池 (仅经典免疫 marker 黑名单) ===\n")
cat(sprintf("放宽 SS LR pairs:     %d\n", nrow(lr_ss_relax)))
cat(sprintf("放宽 SS+CCC LR pairs: %d\n", nrow(lr_cc_relax)))
cat("\n--- SS+CCC LR (放宽) ---\n")
if (nrow(lr_cc_relax) > 0) {
  print(lr_cc_relax[, c("annotation","ligand_subunits","receptor_subunits",
                         "pathway_name","prob","R_log2FC_IAV_vs_Mock","R_FDR",
                         "trend_match","R_expr_vagus_IAV","R_expr_vagus_Mock")],
        row.names = FALSE)
}
sink()
cat(readLines(file.path(OUT, "summary.txt")), sep = "\n")

# ============================================================================
# 11. 可视化
# ============================================================================
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
    x1 = pos_l$x[match(d$ligand_complex,  pos_l$name)],
    y1 = pos_l$y[match(d$ligand_complex,  pos_l$name)],
    x2 = pos_r$x[match(d$receptor_complex, pos_r$name)],
    y2 = pos_r$y[match(d$receptor_complex, pos_r$name)],
    prob = d$prob,
    trend = d$trend_match,
    stringsAsFactors = FALSE
  )
  p <- ggplot() +
    geom_segment(data = edges,
                 aes(x = x1, y = y1, xend = x2, yend = y2,
                     color = trend, alpha = prob, linewidth = prob)) +
    geom_point(data = pos, aes(x, y), color = "grey20", size = 3) +
    geom_text(data = pos_l, aes(x, y, label = name),
              hjust = 1.1, size = 3) +
    geom_text(data = pos_r, aes(x, y, label = name),
              hjust = -0.1, size = 3) +
    scale_color_manual(values = c(co_down   = "#1F77B4",
                                   mirror_up = "#D62728",
                                   no_trend  = "grey60")) +
    scale_linewidth_continuous(range = c(0.3, 2)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    coord_cartesian(xlim = c(-0.4, 1.4), clip = "off") +
    labs(title = title_main,
         subtitle = sprintf("%d LR pairs, 蓝=共下调 红=镜像上调 灰=无方向, 边粗=prob",
                            nrow(d))) +
    theme_void(base_size = 10) +
    theme(plot.margin = margin(15, 80, 15, 80))
  h <- max(4, max(nl, nr) * 0.25 + 2)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = h,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = h,
         dpi = 200, type = "cairo", limitsize = FALSE)
}
draw_network(lr_ss, file.path(OUT, "network_v2_SS_strict"),
             "Mo 58 -> 神经元受体 (SS, 严判)")
draw_network(lr_cc, file.path(OUT, "network_v2_SS_CCC_strict"),
             "Mo 58 -> 神经元受体 (SS+CCC, 严判)")
draw_network(lr_ss_relax, file.path(OUT, "network_v2_SS_relax"),
             "Mo 58 -> 神经元受体 (SS, 放宽: 仅经典黑名单)")
draw_network(lr_cc_relax, file.path(OUT, "network_v2_SS_CCC_relax"),
             "Mo 58 -> 神经元受体 (SS+CCC, 放宽)")
cat("\n>>> 网络图: network_v2_*.pdf/png\n")

saveRDS(list(lr_ss = lr_ss, lr_cc = lr_cc,
             lr_ss_relax = lr_ss_relax, lr_cc_relax = lr_cc_relax,
             match_codown_ss = match_codown_ss,
             match_mirror_ss = match_mirror_ss,
             match_codown_cc = match_codown_cc,
             match_mirror_cc = match_mirror_cc,
             expressed_pool = expressed_pool,
             R_pool_neu_strict = R_pool_neu_strict,
             R_pool_neu_relax = R_pool_neu_relax,
             removed_immune_R = removed_immune_R,
             immune_blacklist = immune_blacklist,
             neuron_in_pool = neuron_in_pool,
             PARAMS = PARAMS),
        file.path(OUT, "v2_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
