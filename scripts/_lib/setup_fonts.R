# scripts/_lib/setup_fonts.R
# 启用中文字体支持 (macOS 内置 STHeiti Medium / Songti)
# 在画图脚本顶部 source 此文件:
#   source(file.path(ROOT, "scripts/_lib/setup_fonts.R"))
#
# 副作用:
#   - 注册字体 'CJK' 与 'CJKsong'
#   - showtext_auto() 全局开启 (ggsave/grDevices 自动接管)
#   - showtext_opts(dpi = 300) 适配高分辨率 PNG

suppressPackageStartupMessages({
  library(sysfonts)
  library(showtext)
})

# 注册第一个可用 CJK 字体
.font_candidates <- c(
  "/System/Library/Fonts/STHeiti Medium.ttc",
  "/System/Library/Fonts/PingFang.ttc",
  "/Library/Fonts/PingFang.ttc",
  "/System/Library/Fonts/STHeiti Light.ttc"
)
.font_chosen <- NULL
for (.fp in .font_candidates) {
  if (file.exists(.fp)) {
    sysfonts::font_add("CJK", regular = .fp)
    .font_chosen <- .fp
    break
  }
}
if (!is.null(.font_chosen)) {
  cat(sprintf(">>> CJK font registered: %s\n", .font_chosen))
}
.song <- "/System/Library/Fonts/Supplemental/Songti.ttc"
if (file.exists(.song)) sysfonts::font_add("CJKsong", regular = .song)

# showtext 全局接管 (PDF + PNG 渲染时自动用注册的 CJK 字体)
# dpi = 100 匹配 ggsave 默认 dpi (避免字号过大放大);
# ggsave 实际用 dpi=200 时, showtext 会自动按比例适配.
showtext::showtext_auto()
showtext::showtext_opts(dpi = 100)

# 给 ggplot2 默认主题指定 CJK family (避免每个图重复 base_family)
if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(.font_chosen)) {
  ggplot2::theme_set(ggplot2::theme_grey(base_family = "CJK"))
}

# 注意: showtext + cairo_pdf 兼容性差; ggsave 用默认 pdf 设备即可,
#       PNG 用 type="cairo-png" 或默认.
# 旧脚本里 device = cairo_pdf 不需要改, showtext_auto 会拦截 grDevices::pdf,
# 但若仍乱码, 把 device = cairo_pdf 移除 (改默认 pdf), 把 type="cairo" 换成 NULL.
