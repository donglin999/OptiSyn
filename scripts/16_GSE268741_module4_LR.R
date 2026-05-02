#!/usr/bin/env Rscript
# 模块 4: NTS 内 Cx3cr1+ 免疫 ↔ Tac1+Vglut2+ 神经元 LR 互作
#
# 数据: 同一份 GSE268741 NTS scRNA-seq (单数据集, 不存在跨数据集尺度问题)
#   sender / receiver 子集均来自模块 1
#     sub_Tac1_Vglut2.rds  (212 cells, Tac1+ 神经元)
#     sub_Cx3cr1.rds       (215 cells, Cx3cr1+ 免疫)
#
# 算法: CellChat 风格 Hill 结合概率
#   双向计算:
#     [N->I]  Tac1+ 神经元 (sender)  -> Cx3cr1+ 免疫 (receiver)
#     [I->N]  Cx3cr1+ 免疫 (sender) -> Tac1+ 神经元 (receiver)
#   E_L_norm = rank(mean_L_in_sender)   / N_genes (sender 子集内 rank-norm)
#   E_R_norm = rank(mean_R_in_receiver) / N_genes (receiver 子集内 rank-norm)
#   prob = H(E_L) * H(E_R), K=0.5, n=2
#   complex 子单元: 取 geometric mean
#
# 重点关注:
#   1) Cx3cl1 -> Cx3cr1 (Tac1+ 神经元到 Cx3cr1+ 免疫)
#   2) 流感相关炎症介质: IFN-γ (Ifng), TNF, 趋化因子 (Ccl/Cxcl), IL-1β/6, etc
#   3) 神经免疫互作: Tac1 (Substance P) -> Tacr1, Glu (Slc17a6) -> NMDA/AMPA
#
# 输出: results/GSE268741/module4_LR/

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)
TS <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

ROOT <- "/Users/wuxiuxiang/project/dongmei/OptiSyn"

# 启用中文字体 (showtext + macOS STHeiti)
source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
OUT  <- file.path(ROOT, "results/GSE268741/module4_LR")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(K = 0.5, n_hill = 2, prob_high = 0.10,
               min_pct_subset = 0.05)   # 子集内 ≥5% 细胞表达才计入

# 1. 加载子集 ---------------------------------------------------------------
cat(">>> 加载子集\n")
seu_tac <- readRDS(file.path(ROOT, "processed_data/GSE268741/sub_Tac1_Vglut2.rds"))
seu_cx3 <- readRDS(file.path(ROOT, "processed_data/GSE268741/sub_Cx3cr1.rds"))
cat(sprintf("    Tac1+Vglut2+   %d cells x %d genes\n", ncol(seu_tac), nrow(seu_tac)))
cat(sprintf("    Cx3cr1+_immune %d cells x %d genes\n", ncol(seu_cx3), nrow(seu_cx3)))

get_mean_pct <- function(seu) {
  d <- GetAssayData(seu, assay = "RNA", layer = "data")
  list(mean_expr = Matrix::rowMeans(d),
       pct_expr  = Matrix::rowMeans(d > 0))
}
m_tac <- get_mean_pct(seu_tac)
m_cx3 <- get_mean_pct(seu_cx3)

# 2. 加载 CellChatDB --------------------------------------------------------
load(file.path(ROOT, "processed_data/CellChatDB/CellChatDB.mouse.rda"))
ia    <- CellChatDB.mouse$interaction
cmplx <- CellChatDB.mouse$complex
cat(sprintf(">>> CellChatDB.mouse: %d LR pairs (SS=%d, CCC=%d, ECM=%d)\n",
            nrow(ia),
            sum(ia$annotation == "Secreted Signaling"),
            sum(ia$annotation == "Cell-Cell Contact"),
            sum(ia$annotation == "ECM-Receptor")))

# Cell-Cell Contact / Secreted Signaling 都纳入 (神经-免疫互作两类都常见)
ia_use <- ia[ia$annotation %in% c("Secreted Signaling","Cell-Cell Contact"), ]
cat(sprintf(">>> 纳入 SS+CCC LR pairs: %d\n", nrow(ia_use)))

expand <- function(x) {
  if (x %in% rownames(cmplx)) {
    s <- as.character(cmplx[x, ]); s[nchar(s) > 0]
  } else x
}

# 3. 双向 Hill 通讯计算 ----------------------------------------------------
hill <- function(E, K = PARAMS$K, n = PARAMS$n_hill) E^n / (K^n + E^n)
geo_mean <- function(v) {
  v <- v[!is.na(v) & v > 0]
  if (length(v) == 0) return(NA_real_)
  exp(mean(log(v)))
}

build_rank_norm <- function(mean_expr) {
  rn <- rank(mean_expr, ties.method = "average") / length(mean_expr)
  names(rn) <- names(mean_expr)
  rn
}
rn_tac <- build_rank_norm(m_tac$mean_expr)
rn_cx3 <- build_rank_norm(m_cx3$mean_expr)

compute_one_dir <- function(sender_name, receiver_name,
                            mean_S, pct_S, rn_S,
                            mean_R, pct_R, rn_R) {
  out <- list()
  for (i in seq_len(nrow(ia_use))) {
    lig <- ia_use$ligand[i]
    rec <- ia_use$receptor[i]
    lig_subs <- expand(lig)
    rec_subs <- expand(rec)

    # L 在 sender 子集中
    if (!all(lig_subs %in% names(mean_S))) next
    L_pct <- min(pct_S[lig_subs])     # complex: 所有亚基都要 ≥ 阈
    if (L_pct < PARAMS$min_pct_subset) next
    L_mean <- geo_mean(mean_S[lig_subs])
    L_rn   <- geo_mean(rn_S[lig_subs])

    # R 在 receiver 子集中
    if (!all(rec_subs %in% names(mean_R))) next
    R_pct <- min(pct_R[rec_subs])
    if (R_pct < PARAMS$min_pct_subset) next
    R_mean <- geo_mean(mean_R[rec_subs])
    R_rn   <- geo_mean(rn_R[rec_subs])

    P <- hill(L_rn) * hill(R_rn)

    out[[length(out) + 1L]] <- data.frame(
      direction = sprintf("%s -> %s", sender_name, receiver_name),
      ligand_complex = lig,
      ligand_subunits = paste(lig_subs, collapse = ";"),
      receptor_complex = rec,
      receptor_subunits = paste(rec_subs, collapse = ";"),
      annotation = ia_use$annotation[i],
      pathway_name = ia_use$pathway_name[i],
      interaction_name = ia_use$interaction_name[i],
      L_mean_sender  = L_mean,  L_pct_sender  = L_pct,  L_rn = L_rn,
      R_mean_receiver = R_mean, R_pct_receiver = R_pct, R_rn = R_rn,
      prob = P,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  d <- do.call(rbind, out); d[order(-d$prob), ]
}

cat(">>> 计算 [N->I] Tac1+ -> Cx3cr1+\n")
lr_n2i <- compute_one_dir("Tac1+_neuron", "Cx3cr1+_immune",
                          m_tac$mean_expr, m_tac$pct_expr, rn_tac,
                          m_cx3$mean_expr, m_cx3$pct_expr, rn_cx3)
cat(sprintf("    LR pairs: %d\n", nrow(lr_n2i)))

cat(">>> 计算 [I->N] Cx3cr1+ -> Tac1+\n")
lr_i2n <- compute_one_dir("Cx3cr1+_immune", "Tac1+_neuron",
                          m_cx3$mean_expr, m_cx3$pct_expr, rn_cx3,
                          m_tac$mean_expr, m_tac$pct_expr, rn_tac)
cat(sprintf("    LR pairs: %d\n", nrow(lr_i2n)))

# 高可信子集
high_n2i <- lr_n2i[lr_n2i$prob >= PARAMS$prob_high, ]
high_i2n <- lr_i2n[lr_i2n$prob >= PARAMS$prob_high, ]
cat(sprintf(">>> 高可信 (prob>=%.2f): N->I=%d, I->N=%d\n",
            PARAMS$prob_high, nrow(high_n2i), nrow(high_i2n)))

# 写表
write.table(lr_n2i, file.path(OUT, "LR_neuron_to_immune_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lr_i2n, file.path(OUT, "LR_immune_to_neuron_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(high_n2i, file.path(OUT, "LR_neuron_to_immune_high.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(high_i2n, file.path(OUT, "LR_immune_to_neuron_high.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(rbind(high_n2i, high_i2n),
          file.path(OUT, "LR_high_confidence_bidirectional.csv"),
          row.names = FALSE)

# 4. 关键 LR 验证 ---------------------------------------------------------
cat("\n>>> 关键 LR 检验\n")
focus <- function(d, lig, rec, dir_lbl) {
  hit <- d[d$ligand_complex == lig & d$receptor_complex == rec, ]
  if (nrow(hit) > 0) {
    cat(sprintf("    [%s] %s -> %s : prob=%.3f, L_pct=%.2f, R_pct=%.2f, pathway=%s\n",
                dir_lbl, lig, rec, hit$prob[1],
                hit$L_pct_sender[1], hit$R_pct_receiver[1],
                hit$pathway_name[1]))
  } else {
    cat(sprintf("    [%s] %s -> %s : NOT in pool (低 pct 或不在 DB)\n",
                dir_lbl, lig, rec))
  }
}
# Cx3cl1 -> Cx3cr1 (神经元 -> 免疫)
focus(lr_n2i, "Cx3cl1", "Cx3cr1", "N->I")
# Tac1 (Substance P) -> Tacr1
focus(lr_n2i, "Tac1", "Tacr1", "N->I")
# 反向: 免疫释放炎症介质给神经元
for (pair in list(c("Tnf","Tnfrsf1a"), c("Tnf","Tnfrsf1b"),
                  c("Il6","Il6ra"), c("Il1b","Il1r1"),
                  c("Ifng","Ifngr1"), c("Cxcl10","Cxcr3"),
                  c("Ccl2","Ccr2"), c("Apoe","Trem2"),
                  c("Spp1","CD44"))) {
  focus(lr_i2n, pair[1], pair[2], "I->N")
}

# 5. 流感相关炎症介质 LR ---------------------------------------------------
cat("\n>>> 流感相关炎症介质过滤\n")
flu_keywords <- c("CCL","CXCL","IL[0-9]","IFN","TNF","IL1","IL6","TGFB",
                  "TLR","CSF","C1Q","C3","COMPLEMENT","TSLP")
infl_pat <- paste(flu_keywords, collapse = "|")
flu_n2i <- lr_n2i[grepl(infl_pat, lr_n2i$pathway_name, ignore.case = TRUE) |
                  grepl(infl_pat, lr_n2i$ligand_complex, ignore.case = TRUE), ]
flu_i2n <- lr_i2n[grepl(infl_pat, lr_i2n$pathway_name, ignore.case = TRUE) |
                  grepl(infl_pat, lr_i2n$ligand_complex, ignore.case = TRUE), ]
cat(sprintf("    流感炎症介质 N->I: %d, I->N: %d\n", nrow(flu_n2i), nrow(flu_i2n)))
write.table(flu_n2i, file.path(OUT, "LR_flu_inflammation_N_to_I.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(flu_i2n, file.path(OUT, "LR_flu_inflammation_I_to_N.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6. 互作强度量化 (按通路汇总) ---------------------------------------------
pathway_strength <- function(d, lab) {
  if (nrow(d) == 0) return(data.frame())
  agg <- aggregate(prob ~ pathway_name, data = d,
                   FUN = function(v) c(n_LR = length(v),
                                       max_prob = max(v),
                                       sum_prob = sum(v)))
  data.frame(direction = lab, pathway = agg$pathway_name,
             n_LR     = as.integer(agg$prob[, "n_LR"]),
             max_prob = as.numeric(agg$prob[, "max_prob"]),
             sum_prob = as.numeric(agg$prob[, "sum_prob"]),
             stringsAsFactors = FALSE)
}
pw_n2i <- pathway_strength(lr_n2i, "N->I")
pw_i2n <- pathway_strength(lr_i2n, "I->N")
pw_both <- rbind(pw_n2i, pw_i2n)
pw_both <- pw_both[order(-pw_both$sum_prob), ]
write.table(pw_both, file.path(OUT, "pathway_strength.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n>>> Top 20 通路 (按 sum_prob):\n")
print(head(pw_both, 20), row.names = FALSE)

# 7. 摘要日志 -------------------------------------------------------------
sink(file.path(OUT, "module4_summary.txt"))
cat("=== 模块 4: NTS Cx3cr1+ ↔ Tac1+Vglut2+ LR 互作 ===\n")
cat(sprintf("时间: %s\n", TS))
cat("数据: GSE268741 模块 1 子集 (同一数据集, 无跨数据集问题)\n\n")

cat("=== 子集 ===\n")
cat(sprintf("Tac1+Vglut2+ 神经元: %d cells\n", ncol(seu_tac)))
cat(sprintf("Cx3cr1+ 免疫:        %d cells\n", ncol(seu_cx3)))

cat(sprintf("\n=== 参数 ===\n")); print(PARAMS)
cat(sprintf("\nLR 数据库: CellChatDB.mouse SS+CCC, %d pairs\n", nrow(ia_use)))

cat(sprintf("\n=== 双向 LR 通讯计算 ===\n"))
cat(sprintf("[N->I] Tac1+ 神经元 -> Cx3cr1+ 免疫 : %d LR pairs (高可信 %d)\n",
            nrow(lr_n2i), nrow(high_n2i)))
cat(sprintf("[I->N] Cx3cr1+ 免疫 -> Tac1+ 神经元 : %d LR pairs (高可信 %d)\n",
            nrow(lr_i2n), nrow(high_i2n)))

cat("\n=== Top 20 [N->I] (按 prob 降序) ===\n")
if (nrow(lr_n2i) > 0) {
  print(head(lr_n2i[, c("ligand_complex","receptor_complex","pathway_name",
                        "annotation","prob","L_pct_sender","R_pct_receiver")], 20),
        row.names = FALSE)
}

cat("\n=== Top 20 [I->N] (按 prob 降序) ===\n")
if (nrow(lr_i2n) > 0) {
  print(head(lr_i2n[, c("ligand_complex","receptor_complex","pathway_name",
                        "annotation","prob","L_pct_sender","R_pct_receiver")], 20),
        row.names = FALSE)
}

cat("\n=== 关键 LR 检验 ===\n")
key_pairs_n2i <- list(c("Cx3cl1","Cx3cr1"), c("Tac1","Tacr1"))
for (pp in key_pairs_n2i) {
  hit <- lr_n2i[lr_n2i$ligand_complex == pp[1] &
                lr_n2i$receptor_complex == pp[2], ]
  if (nrow(hit) > 0) {
    cat(sprintf("[N->I] %s -> %s: prob=%.3f, L_pct=%.2f, R_pct=%.2f\n",
                pp[1], pp[2], hit$prob[1],
                hit$L_pct_sender[1], hit$R_pct_receiver[1]))
  } else {
    cat(sprintf("[N->I] %s -> %s: NOT in pool\n", pp[1], pp[2]))
  }
}
key_pairs_i2n <- list(c("Tnf","Tnfrsf1a"), c("Tnf","Tnfrsf1b"),
                      c("Il6","Il6ra"), c("Il1b","Il1r1"),
                      c("Ifng","Ifngr1"), c("Cxcl10","Cxcr3"),
                      c("Ccl2","Ccr2"), c("Apoe","Trem2"))
for (pp in key_pairs_i2n) {
  hit <- lr_i2n[lr_i2n$ligand_complex == pp[1] &
                lr_i2n$receptor_complex == pp[2], ]
  if (nrow(hit) > 0) {
    cat(sprintf("[I->N] %s -> %s: prob=%.3f, L_pct=%.2f, R_pct=%.2f\n",
                pp[1], pp[2], hit$prob[1],
                hit$L_pct_sender[1], hit$R_pct_receiver[1]))
  } else {
    cat(sprintf("[I->N] %s -> %s: NOT in pool\n", pp[1], pp[2]))
  }
}

cat("\n=== 流感相关炎症 LR ===\n")
cat(sprintf("N->I (神经元->免疫): %d\n", nrow(flu_n2i)))
if (nrow(flu_n2i) > 0) {
  print(flu_n2i[, c("ligand_complex","receptor_complex","pathway_name",
                    "prob","L_pct_sender","R_pct_receiver")],
        row.names = FALSE)
}
cat(sprintf("\nI->N (免疫->神经元): %d\n", nrow(flu_i2n)))
if (nrow(flu_i2n) > 0) {
  print(flu_i2n[, c("ligand_complex","receptor_complex","pathway_name",
                    "prob","L_pct_sender","R_pct_receiver")],
        row.names = FALSE)
}

cat("\n=== 通路强度量化 (Top 20 sum_prob) ===\n")
print(head(pw_both, 20), row.names = FALSE)
sink()
cat("\n", paste(readLines(file.path(OUT, "module4_summary.txt")), collapse="\n"), sep="")

# 8. 可视化 ---------------------------------------------------------------
cat(">>> 可视化\n")

# 8a. Top 30 LR 双向网络条形图 (按 prob)
top_n2i <- head(lr_n2i, 30)
top_i2n <- head(lr_i2n, 30)
top_n2i$lr <- paste0(top_n2i$ligand_complex, " → ", top_n2i$receptor_complex)
top_i2n$lr <- paste0(top_i2n$ligand_complex, " → ", top_i2n$receptor_complex)
top_n2i$direction <- "N -> I"
top_i2n$direction <- "I -> N"

bar_df <- rbind(top_n2i[, c("lr","prob","direction","pathway_name","annotation")],
                top_i2n[, c("lr","prob","direction","pathway_name","annotation")])
bar_df$lr <- factor(bar_df$lr, levels = unique(bar_df$lr[order(-bar_df$prob)]))
p_bar <- ggplot(bar_df, aes(prob, lr, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("N -> I" = "#D62728", "I -> N" = "#1F77B4")) +
  facet_wrap(~ direction, ncol = 2, scales = "free_y") +
  labs(x = "结合概率 (Hill K=0.5, n=2)", y = NULL,
       title = "Top 30 LR 双向通讯",
       subtitle = "Tac1+Vglut2+ 神经元 ↔ Cx3cr1+ 免疫细胞") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 7))
ggsave(file.path(OUT, "top30_LR_bar.pdf"), p_bar,
       width = 13, height = 9, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "top30_LR_bar.png"), p_bar,
       width = 13, height = 9, dpi = 200, type = "cairo", limitsize = FALSE)

# 8b. Cx3cl1 - Cx3cr1 表达 boxplot
draw_lr_box <- function(seu_S, seu_R, gene, title_main, out_file) {
  d_S <- GetAssayData(seu_S, layer = "data")[gene, ]
  d_R <- GetAssayData(seu_R, layer = "data")[gene, ]
  df <- rbind(
    data.frame(set = "Tac1+_neuron", expr = as.numeric(d_S)),
    data.frame(set = "Cx3cr1+_immune", expr = as.numeric(d_R))
  )
  p <- ggplot(df, aes(set, expr, fill = set)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, outlier.size = 0.5) +
    scale_fill_manual(values = c("Tac1+_neuron" = "#D62728",
                                 "Cx3cr1+_immune" = "#1F77B4")) +
    labs(x = NULL, y = "log-normalized expression",
         title = sprintf("%s 表达 (%s)", gene, title_main)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  ggsave(out_file, p, width = 5, height = 4, device = cairo_pdf)
  ggsave(sub("\\.pdf$",".png", out_file), p, width = 5, height = 4,
         dpi = 200, type = "cairo")
}
if ("Cx3cl1" %in% rownames(seu_tac)) {
  draw_lr_box(seu_tac, seu_cx3, "Cx3cl1", "在两子集中",
              file.path(OUT, "Cx3cl1_violin.pdf"))
}
if ("Cx3cr1" %in% rownames(seu_tac)) {
  draw_lr_box(seu_tac, seu_cx3, "Cx3cr1", "在两子集中",
              file.path(OUT, "Cx3cr1_violin.pdf"))
}

# 8c. 通路 sum_prob 条形图 Top 25 (双向合并)
pw_top <- head(pw_both, 25)
pw_top$pathway <- factor(pw_top$pathway,
                         levels = rev(unique(pw_top$pathway)))
p_pw <- ggplot(pw_top, aes(sum_prob, pathway, fill = direction)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = c("N->I" = "#D62728", "I->N" = "#1F77B4")) +
  labs(x = "通路总互作强度 (sum prob)", y = NULL,
       title = "Top 25 通路互作强度",
       fill = "方向") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(OUT, "pathway_strength_top25.pdf"), p_pw,
       width = 9, height = 8, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "pathway_strength_top25.png"), p_pw,
       width = 9, height = 8, dpi = 200, type = "cairo", limitsize = FALSE)

# 8d. 双向 LR 简化网络图 (Top 15 双向合并)
all_top <- rbind(head(lr_n2i, 15), head(lr_i2n, 15))
all_top$L_node <- paste0(all_top$ligand_complex, "@",
                          ifelse(grepl("^Tac1", all_top$direction), "N", "I"))
all_top$R_node <- paste0(all_top$receptor_complex, "@",
                          ifelse(grepl("immune$", all_top$direction), "I", "N"))
nodes <- unique(c(all_top$L_node, all_top$R_node))
nl <- length(nodes)
node_df <- data.frame(name = nodes,
                      side = sub(".*@","", nodes),
                      stringsAsFactors = FALSE)
node_df$x <- ifelse(node_df$side == "N", 0, 1)
n_per <- table(node_df$side)
node_df$y <- NA
for (sd in names(n_per)) {
  idx <- which(node_df$side == sd)
  ys <- if (length(idx) > 1) seq(1, 0, length.out = length(idx)) else 0.5
  node_df$y[idx] <- ys
}
edge_df <- data.frame(
  x1 = node_df$x[match(all_top$L_node, node_df$name)],
  y1 = node_df$y[match(all_top$L_node, node_df$name)],
  x2 = node_df$x[match(all_top$R_node, node_df$name)],
  y2 = node_df$y[match(all_top$R_node, node_df$name)],
  prob = all_top$prob,
  direction = all_top$direction,
  stringsAsFactors = FALSE
)
p_net <- ggplot() +
  geom_segment(data = edge_df, aes(x = x1, y = y1, xend = x2, yend = y2,
                                    color = direction,
                                    alpha = prob, linewidth = prob)) +
  geom_point(data = node_df, aes(x, y, color = side), size = 4) +
  geom_text(data = node_df, aes(x, y, label = sub("@.*","", name),
                                 hjust = ifelse(x == 0, 1.1, -0.1)),
            size = 3.2, color = "black") +
  scale_color_manual(values = c(
    N = "#D62728", I = "#1F77B4",
    "Tac1+_neuron -> Cx3cr1+_immune" = "#D62728",
    "Cx3cr1+_immune -> Tac1+_neuron" = "#1F77B4")) +
  scale_linewidth_continuous(range = c(0.4, 2)) +
  scale_alpha_continuous(range = c(0.4, 1)) +
  coord_cartesian(xlim = c(-0.4, 1.4), clip = "off") +
  labs(title = "NTS Neuro-Immune LR Network (Top 15 bidirectional)",
       subtitle = "Left=N(Tac1+ neuron), Right=I(Cx3cr1+ immune); red=N->I, blue=I->N") +
  theme_void(base_size = 10) +
  theme(plot.margin = margin(15, 80, 15, 80),
        legend.position = "none",
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
ggsave(file.path(OUT, "LR_network_top15_bidirectional.pdf"), p_net,
       width = 11, height = 9, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "LR_network_top15_bidirectional.png"), p_net,
       width = 11, height = 9, dpi = 200, type = "cairo", limitsize = FALSE)

cat("    可视化完成\n")

# 9. 保存中间对象
saveRDS(list(lr_n2i = lr_n2i, lr_i2n = lr_i2n,
             high_n2i = high_n2i, high_i2n = high_i2n,
             flu_n2i = flu_n2i, flu_i2n = flu_i2n,
             pw_both = pw_both, PARAMS = PARAMS,
             m_tac = m_tac, m_cx3 = m_cx3),
        file.path(OUT, "module4_results.rds"))

cat(sprintf("\n>>> 输出: %s\n", OUT))
