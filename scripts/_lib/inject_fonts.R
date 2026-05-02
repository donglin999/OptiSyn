# 在 11 个画图脚本里注入 source(setup_fonts.R)
files <- c("scripts/05_GSE42639_Mo_protective_targets.R",
           "scripts/06_GSE42639_Mo_protective_strict.R",
           "scripts/08_GSE31022_validation_58targets.R",
           "scripts/09_GSE161878_DESeq2.R",
           "scripts/10_CellChat_Mo_to_vagus.R",
           "scripts/12_CellChat_neuron_compartment_LR.R",
           "scripts/13_Cx3cl1_to_Cx3cr1_reverse_LR.R",
           "scripts/15_GSE268741_module1_clustering.R",
           "scripts/16_GSE268741_module4_LR.R",
           "scripts/17_module5_exclusivity.R",
           "scripts/18_module2_3_enrichment.R")

for (f in files) {
  txt <- readLines(f)
  if (any(grepl("setup_fonts.R", txt, fixed = TRUE))) {
    cat(sprintf("    [skip] %s already has setup_fonts\n", f)); next
  }
  i <- grep("^ROOT[ ]*<-", txt)[1]
  if (is.na(i)) {
    cat(sprintf("    [warn] %s no ROOT line, skip\n", f)); next
  }
  inject <- c(
    "",
    "# 启用中文字体 (showtext + macOS STHeiti)",
    'source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))'
  )
  txt <- append(txt, inject, after = i)
  writeLines(txt, f)
  cat(sprintf("    [ok] %s injected at line %d\n", f, i + 1))
}
