% OptiSyn 项目开题汇报
% xiuxiang
% 2026-05-03

## 封面

**靶向"外周肺-中枢孤束核"神经免疫互作轴的<br/>流感后咳嗽机制及中药干预深度报告**

项目开题汇报

汇报人: xiuxiang | 日期: 2026-05-03

数据: GSE42639 / GSE31022 / GSE161878 / GSE296065 / GSE268741 + CellChatDB.mouse

## Slide 2 — 临床问题 + 神经免疫学范式转移

**流感后慢性咳嗽 (Post-influenza chronic cough)**: 流感感染急性期过后 >3 周持续咳嗽，临床困境长期未解；单纯抗病毒/镇咳治疗效果有限。

**范式转移** — "气道局部炎症" → "**神经-免疫互作环路**"

| 旧范式 | 新范式 |
|---|---|
| 病毒清除后炎症残留 | TRPV1+ vagal 感觉神经致敏 |
| 局部气道高反应 | NTS Tac1+ 神经元中枢放大 |
| 单一靶点干预 | "外周肺-迷走-脑干"轴 |

**本研究核心**: 系统刻画 Mo (Cx3cr1+) → vagus → NTS Tac1+/Cx3cr1+ 微胶质完整环路，定位中药多靶点干预节点。

## Slide 3 — 核心假设：三段闭环

**总命题**:

```
外周肺 Cx3cr1+ Mo 保护轴失效  →  迷走神经 Cx3cl1 配体投射上传
                              ↓
NTS Tac1+ Vglut2+ 神经元 ←→ NTS Cx3cr1+ 微胶质双向通讯
                              ↓
                       中枢咳嗽信号放大 → 慢性咳嗽病理
```

**信号轴中央**: Cx3cl1-Cx3cr1

- **配体端 (Cx3cl1)**: vagus 节神经元 + NTS Tac1+ 神经元
- **受体端 (Cx3cr1)**: 外周肺 Mo + NTS 微胶质
- **双端崩塌**: Mo 端受体下调 (slope=−0.80) + Am 端配体下调 (slope=−0.60)

## Slide 4 — 三级科学问题

| 假设 | 物理位置 | 待证 | 数据集 |
|---|---|---|---|
| **H1 外周诱发缺口** | 肺髓系 | Cx3cr1+ Mo 数量/功能在感染剂量下崩塌 | GSE42639, GSE31022, GSE296065 |
| **H2 神经传导桥梁** | 迷走神经节 | vagus 端"受体丢失"+ Cx3cl1 配体上传 | GSE161878 + CellChatDB |
| **H3 中枢放大缺口** | 脑干 NTS | Tac1+ 神经元 ↔ Cx3cr1+ 微胶质双向通讯 | GSE268741 |

**三段闭环 → 三个干预节点 → 中药多靶点对应**

## Slide 5 — 总体技术路线

```
[① 外周端]
GSE42639 (5 cells × 5 dose)  →  Mo 严选 564→58→17 高优先
        ↓
GSE31022 (D1-D6 跨毒株时序)  →  4 双标 (Gstm1/Dusp18/Insr/6430548M08Rik)
        ↓
GSE296065 (lung sc + RTX)    →  Cx3cr1 唯一通过 BH 双判

[② 神经传输]
GSE161878 (vagus bulk)       →  144 上 / 10 下 (含 Gpr151 nociceptor)
        +
CellChatDB v3 反向            →  Cx3cl1 → Cx3cr1, prob=0.426

[③ 中枢端]
GSE268741 (NTS sc, 5 模块)   →  Tac1+ ↔ Cx3cr1+ 双向, LR prob=0.602
                                Tac1+ 排他性 4 维度第一
```

## Slide 6 — 数据资源

| 数据集 | 类型 | 用途 | 样本量 | 假设链中的角色 |
|---|---|---|---|---|
| GSE42639 | Illumina BeadChip | 5 细胞 × 5 H1N1 剂量 | 75 (Mo 15) | H1 主筛 |
| GSE31022 | Illumina BeadChip | H1N1 D1-D6 时序 | 21 | H1 跨毒株 |
| GSE161878 | RNA-seq | vagus IAV vs Mock | 21 (剔 D-suffix) | H2 主筛 |
| GSE296065 | scRNA (Seurat) | lung 22 cells + RTX | 37,543 cells | H1 单细胞 + 药理 |
| GSE268741 | scRNA (10x) | NTS 神经-免疫 | 9,638 cells | H3 中枢端 |
| CellChatDB.mouse | 配-受体库 | LR 配对 | 2019 LR pairs | H2 桥接 |

## Slide 7 — 分析 1: 5 细胞 × 4 contrasts DEG (对照实验)

**方法**: limma + eBayes + BH | **阈值**: FDR<0.05 & |log2FC|≥1

| Cell | TX91 | PR8 0.6 | PR8 10 | PR8 100 | **总和** |
|---|---:|---:|---:|---:|---:|
| **Mo** | **373** | **1317** | **1551** | **1788** | **5029** ★ |
| Ne | 139 | 555 | 1099 | 1193 | 2986 |
| Nh | 99 | 520 | 1211 | 979 | 2809 |
| Am | 87 | 422 | 914 | 979 | 2402 |
| Ly | 77 | 408 | 349 | 459 | 1293 |

**Mo 在所有 4 个 contrasts 全部第一**, 领先第 2 名 1.5×–2.7×。

![](../results/DEG_5cell_dose/analysis1_5cell_DEG_barplot.png)

## Slide 8 — 分析 2: Mo 剂量依赖 DEG 单调验证

**目的**: 验证 Mo 应答规模随感染剂量 **严格单调递增** → 锁定"剂量响应感受器"。

| 指标 | TX91 | PR8 0.6 | PR8 10 | PR8 100 | 单调? |
|---|---:|---:|---:|---:|:---:|
| 总 DEG | 373 | 1317 | 1551 | 1788 | ✓ 严格递增 |
| Up | 321 | 783 | 864 | 965 | ✓ 严格递增 |
| Down | 52 | 534 | 687 | 823 | ✓ 严格递增 (16× 跨度) |

**Down 增幅最大** (52→823, 16×) → 与 "重症崩塌" 假设逻辑闭环。

![](../results/DEG_5cell_dose/analysis2_Mo_dose_DEG_lineplot.png)

## Slide 9 — Mo 保护性靶点严选: 5 步 + 严选叠加

**生物学定义**: "TX91 轻症维持 → PR8 剂量上升持续下降 → 100LD50 显著崩塌"

| 步骤 | 模型 / 阈值 | 累计 |
|---|---|---:|
| 起始 | 探针注释 | 16,005 |
| Step 1 趋势 | slope<0 & FDR<0.10 | 844 |
| Step 2 崩塌 | PR8_100 vs TX91 log2FC≤−0.585 | 645 |
| Step 3 基础表达 | TX91 mean ≥ 全基因中位 | 601 |
| Step 4 反向假阳剔除 | — | 601 |
| Step 5 排已下调 | 剔 TX91 vs Sham 已下调 | **564** |
| **Step A 阶梯下降** | PR8_0.6 + PR8_10 严判 | **65** |
| **Step B 黑/白名单** | 剔 7 纯代谢酶 | **58 严选 ★** |

## Slide 10 — 17 高优先候选 + 4 跨毒株双标

**17 高优先候选** (炎症/免疫调控): `Cx3cr1`, `Klf2`, `Arrb1`, `Abcd2`, `Mbp`, `Il16`, `Cd177`, `Rasgrp2`, `Rnf144a`, `Bcl11a`, `Btla`, `Cd209a`, `Rnase2b`, `Ampd3`, `Rasgrp4`, `Dock8`, `Rin3`

**4 双标基因** (GSE31022 D1-D6 跨毒株时序最稳健):

| 基因 | best logFC | slope (D1-D6) | 功能 |
|---|---:|---:|---|
| **Gstm1** | −2.67 | −0.426 | 谷胱甘肽 S-转移酶, 抗氧化 ★ |
| Dusp18 | −1.25 | −0.352 | MAPK 通路负调 |
| Insr | −0.63 | −0.151 | 胰岛素受体 / SignalTransduction |
| 6430548M08Rik | −0.69 | −0.146 | RIKEN cDNA |

**Cx3cr1 / Gstm1 / Dusp18** 三个核心靶点 → 对应莲花清瘟多靶点干预（Slide 18）。

## Slide 11 — 跨系统导火桥梁: vagus DESeq2

**GSE161878 — 迷走神经节 IAV vs Mock**

**关键修复**: AM312D/313D/314D 是同动物技术重复 → 主分析剔除 (11+10) → **144 上 / 10 下**。含 D 版 343/122 仅 32% 一致 → 同动物伪重复污染必须剔除。

**上调主程序**:

- 髓系浸润: `Ly6c2`, `Tyrobp`, `Cybb`, `Mpeg1`, `Ctss`, `Fcgr1/4`
- MHC-II 抗原呈递: `H2-Aa`, `H2-Q7`, `Cd74`
- Type I IFN 抗病毒: `Oasl2`, `Cxcl9`, `Cfb`, `Tlr2`, `Ubd`
- ★ **Gpr151** (vagal nociceptor 标志, FDR=1.7e−4)

**下调 (10 个)**: 神经元结构/突触维持基因 → 神经元稳态损伤。

## Slide 12 — CellChat 通讯 (3 次迭代 → 反向方向修正)

| 版本 | 方向 | 命中 | 状态 |
|---|---|---|---|
| v1 严判 | Mo → Vagus | **0** | 受体 (Tnfrsf14/Cd4) 在 vagus 21/21 样本 CPM<1 |
| v2 放宽 | 神经元区室受体池 | 1 | 弱 |
| **v3 反向 ★** | **Vagal Cx3cl1 → Mo Cx3cr1** | **prob=0.426** | 强命中 |

**关键发现**:

- vagus Cx3cl1 CPM 2.5–8.5, **21/21 样本固有表达** (IAV 不上调, 是固有信号源)
- Mo Cx3cr1 重症崩塌 slope=−0.80 → 受体端故障
- 双数据集独立支持: NTS sc 模块 4 也得到 Cx3cl1→Cx3cr1 prob=0.602

⚠️ **诚实标注**: v1→v3 是 post-hoc refinement, 论文 Methods 必按时间顺序写; 建议补 vagal Cx3cl1 蛋白染色 (IHC) 独立验证。

## Slide 13 — NTS 中枢端: GSE268741 5 模块 (模块 1-3)

**模块 1 — 分群 + 子集分选**

- 输入: 31 xlsx → sparse 33,696 × 9,638
- QC: 6,468 cells 通过, 22 簇
- **关键子集**: Tac1+ Vglut2+ 兴奋性神经元 **212 cells**, Cx3cr1+ 微胶质 **215 cells**, 其它神经元 2,472 cells

**模块 2 — Tac1+ marker 富集**

- FindMarkers Tac1+ vs Others: 506 上 / 310 下
- GO BP: **74 通路 FDR<0.05** (突触兴奋性 / 谷氨酸 / 痛觉传递)

**模块 3 — 跨数据集通路一致性**

- vagal vs Tac1+ marker: 4 共享 secretion 通路 + Lck 共享 signaling 基因

## Slide 14 — NTS 模块 4+5: LR 互作 + Tac1+ 排他性

**模块 4 — 神经-免疫双向 LR ★**

| LR 方向 | 强度 | 含义 |
|---|---:|---|
| **Cx3cl1 → Cx3cr1** (Tac1+ → 微胶质) | **prob=0.602** | 中枢神经向微胶质投射 |
| 神经→免疫 (全) | 227 LR pairs | 双向通讯网络 |
| 免疫→神经 (Tnfrsf1a + Ackr1 反馈) | 135 LR pairs | 微胶质回馈细胞因子 |

**模块 5 — Tac1+ 排他性 4 维度**: LR 互作强度 / 咳嗽通路富集 / Cx3cl1 表达 pct (69.3%) / Wilcox vs 其它 — **全部 Tac1+ 第一**。

⚠️ **单 GSM 限制**: 9,638 cells 来自 1 GSM 样本, 严格意义上属 pseudoreplication; 论文必须明确标注 + **务必寻找第 2 份独立 NTS sc 数据**。

## Slide 15 — 单细胞 RTX 救援: GSE296065

**实验设计漏斗**:

```
58 严选 → 28 通过 specificity (Target vs BG, BH 校正)
         → 1 通过 BH 双判 (V7 vs R7 组间) ★
            └─ Cx3cr1: log2FC=0.86, FDR<0.05
```

**Cx3cr1 单细胞特征**:

| 指标 | Target Mo | Background |
|---|---:|---:|
| pct | **46%** | 11% |
| Specificity log2FC | **2.62** | — |
| RTX 救援 (V7 vs R7) | **+0.86** | — |

**解读**: RTX (TRPV1+ 神经阻断) → 外周 Mo Cx3cr1 表达回升 → 中枢→外周反馈被阻断后保护轴恢复; 28→1 漏斗 → **不能宣称"58 严选都是 RTX 靶点"**, 只 Cx3cr1 是核心枢纽。

## Slide 16 — 完整机制链路 (6 节点闭环)

```
[① 外周肺]
肺 Mo Cx3cr1 受体 (TX91=9.88 → PR8_100=7.50, slope=−0.80)
       ↑                                   ↓ 保护轴失效
       │                                   ↓
[② vagus]                          趋化/抗炎/稳态三轴下游下调
vagus Cx3cl1 (CPM 2.5–8.5,                (58 严选: 7+4+9 命中)
固有表达 21/21)
       │ axonal projection
       ↓
[③ NTS]
NTS Tac1+ Vglut2+ 神经元 (212 cells)
       │ Cx3cl1 释放 (pct 69.3%)
       ↓ prob=0.602
NTS Cx3cr1+ 微胶质 (215 cells)
       │ 反馈 TNF/Cxcl10/Ccl2 (Tnfrsf1a + Ackr1)
       ↓
[④ Tac1+ 突触兴奋性] → 谷氨酸 + 痛觉传递通路 (74 GO BP)
       │ vagal 投射回外周
       ↓
[⑤ 咳嗽病理]    [⑥ 中药节点: 莲花清瘟多靶点干预]
```

## Slide 17 — 候选靶点优先级 (A/B/C 档)

**A 档** (强证据 — 跨多数据集 + 机制清晰):

| 靶点 | 数据集贯穿 | 机制 |
|---|---|---|
| **Cx3cr1** ★★★ | 5 个数据集 | RTX 救援唯一通过, slope=−0.80, 受体端枢纽 |
| **Cx3cl1** ★★★ | vagus + NTS | 21/21 固有表达, NTS Tac1+ pct 69%, LR prob=0.602 |

**B 档** (多数据集证据): Klf2 / Arrb1 / Btla / Il16 / Cd209a — Mo 严选 + GSE31022 同向

**C 档** (探索性): **Gstm1** (跨毒株 4 双标核心, 抗氧化) / Dusp18 / Insr / 6430548M08Rik / **Gpr151** (vagal nociceptor 标志, 与"咳嗽"主线契合)

## Slide 18 — 中药多靶点干预: 莲花清瘟整合调控

**核心病理节点 ↔ 莲花清瘟成分映射**:

| 中药学机制 | 对应肺感染病变节点 | 我们项目命中靶点 |
|---|---|---|
| 抗病毒 + 免疫调控 (黄芩苷/连翘苷) | 病毒清除 + 抗炎调节 | `H2-Aa`, `H2-Q7`, `Oasl2`, `Cxcl9` |
| 抗氧化 (鱼腥草苷/GSH 通路) | 氧化应激缓解 | **Gstm1** ★ (4 双标核心) |
| 抗炎细胞因子 (大黄素/TNF-α 抑制) | TNF/Cxcl 反馈链 | `Tnfrsf1a`, `Ackr1` |
| **神经-胶质细胞通讯调控** | NTS Tac1+ ↔ Cx3cr1+ 微胶质 | **Cx3cr1 / Cx3cl1** ★ |

**核心干预假设**: 莲花清瘟 = 多靶点联合介入 H1+H2+H3 三段闭环; 不是单一抗病毒, 而是同时缓解外周 Mo 保护轴衰退 + vagus 端致敏 + 中枢 NTS 神经-免疫反馈。

## Slide 19 — 后续四大实验模块（湿实验设计）

| 模块 | 核心方法 | 关键评分 / 终点 | 创新性 |
|---|---|---|---|
| **① 体内药效学** | LD10/3 → LD50/3 + WBP 体描记 | 咳嗽次数/强度/吸气幅度 | 标准基线 |
| **② 外周机制验证** | 流式 + scRNA + Western + qPCR | Cx3cr1+ Mo 数量, Cx3cr1/Klf2/Arrb1 表达 | Mo 单细胞验证 |
| **③ 中枢放大验证** | 多色 IHC + Luminex + 磷酸化通路 | NTS Tac1+ ↔ Cx3cr1+ 信号通路, IL-6/TNF | NTS 区域通讯 |
| **④ 化学遗传学体内验证 ★** | **DREADDs + Tac1-Cre + AAV-DIO-hM4Di + i.p. CNO** | 抑制 Tac1+ → 咳嗽是否被阻断 | **本项目最大创新点** |

**模块 ④ — 黄金标准**: AAV1-DIO-hM4D(Gi)-mCherry → Tac1-Cre 小鼠 NTS, i.p. CNO 抑制 Tac1+ 神经元活性 → 直接因果性证明 NTS Tac1+ 是流感咳嗽的中枢驱动节点。

## Slide 20 — 创新性 + 可行性

**创新性 (4 维)**:

1. 跨学科理论框架 — "外周肺-迷走-脑干"轴的完整刻画 (而非局部气道)
2. 数据规模 + 多组学 — 5 个 GEO 数据集 + sc + bulk + LR 互作
3. 化学遗传学因果验证 — DREADDs Tac1+ 抑制是黄金标准, 国内同类研究稀缺
4. 中医-神经免疫学融合 — 莲花清瘟多靶点 vs 现代神经免疫学环路

**可行性 (4 维)**:

1. ✅ 数据已就位 — 5 个 GEO + CellChatDB 全部下载, 流程跑通
2. ✅ 关键预实验已完成 — Mo 严选 58 + Cx3cr1 单细胞 RTX 通过 + NTS LR prob=0.602
3. ✅ 靶点优先级清晰 — A 档 (Cx3cr1/Cx3cl1) 直接进入湿实验验证
4. ⚠️ DREADDs 实验 — 需要 Tac1-Cre 小鼠系 + AAV 病毒制备 + 体描记设备 (需协作单位)

## Slide 21 — 风险与限制 (透明声明)

| # | 风险 | 严重度 | 应对 |
|---|---|:---:|---|
| 1 | **GSE268741 单 GSM** 无生物学重复 | 🔴 | 寻找第 2 份独立 NTS sc 数据 (Allen Brain / 协作组) |
| 2 | CellChat v1→v3 是 post-hoc refinement | 🟠 | 论文如实按时间顺序写 + IHC 蛋白染色独立验证 |
| 3 | RTX 救援 28→1 衰减 | 🟠 | 仅 Cx3cr1 通过, 不宣称"58 全是 RTX 靶点" |
| 4 | 5 cells StepA 候选两两交集=0 | 🟠 | 加做 within-cell rank sensitivity 分析 |
| 5 | GSE161878 D-suffix 是推断 | 🟡 | 联系投稿者确认 + Supplementary 双版本 |
| 6 | TX91 是带毒 baseline | 🟡 | Figure legend 必明确 |
| 7 | NTS Cx3cr1+ pct=100% 是定义同义反复 | 🟡 | "operationally-defined pool" |

**态度**: 主结论稳健 + 异常透明 + 路径清晰 → **A− 总评**

## Slide 22 — 总结 + Q & A

**三句话总结**:

1. **机制链已闭环**: 外周 Mo Cx3cr1 → vagus Cx3cl1 → NTS Tac1+ ↔ Cx3cr1+ 微胶质 → 咳嗽
2. **核心枢纽锁定**: Cx3cl1-Cx3cr1 信号轴 是干预的 single most actionable node
3. **下一步**: DREADDs Tac1+ 抑制 + 莲花清瘟多靶点干预 = 黄金标准 in vivo 验证

**已完成产出**: 6 个数据集预处理 + DEA + 多组学整合; 14 个分析脚本; 36 张图 + 2 份完整报告; GitHub 项目 OptiSyn (commit `4431654`)

**讨论问题**:

1. NTS 单 GSM 限制, 优先用 Allen Brain 还是再做一份独立 sc 实验?
2. DREADDs 实验的 Tac1-Cre 小鼠系是否有现成协作单位?
3. 莲花清瘟"剂量爬升"实验设计 (LD10/3 vs LD50/3) 是否合理?

致谢: 感谢导师指导 + 实验室同仁 + GEO 数据贡献者
