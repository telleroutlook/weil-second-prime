# 研究计划 — weil-second-prime 第二素数窗口

**创建日期：2026-08-08**

> 新程序员接手先读 `HANDOFF.md`（操作交接）。本 PLAN 是战略地图。
> 数学纪律见 `docs/PROOF_CONSTITUTION.md`（从 weil-first 移植）。

---

## 状态快照（一句话真相）

**脚手架阶段。尚无证书。** 窗口 $L \in (\tfrac12\log3,\ \log2) \approx (0.549, 0.693)$，
两素数 $p=2,3$ 均在单跳区，故素数层是真正的**双平移**对象，带交叉耦合
$J_{ij}(\tau_2,\tau_3)$，$\tau_p=\log p/L$。$c_L$ 在窗口右端约 1.82。

---

## 第零编 · 定位（为什么做这个，双目的）

这是从 weil-first 的收尾中得到的**修正战略**（2026-08-08）：

1. **真实数学是宿主。** 第二窗口是新数学（双素数耦合、新代数 $\mathbb{Q}[\tau_2,\tau_3]$、
   新数值陷阱），本身有可发表价值（方法射程扩展）。
2. **proofctl 是随宿主进化的器官，不是独立主线。** 验证工具脱离真实被验证对象就退化成
   "自己测自己"——它能想到的自测恰恰都是它已会防的，发现不了盲区。盲区只有真实研究的
   意外才能暴露。第一窗口的 C10/C11 正是被 weil-first 的真实 bug（copy 证书、S0 漏项）
   逼出来的。第二窗口的新 bug 类（交叉项漏项、不可通约平移的数值陷阱）预期驱动下一轮
   proofctl 条件。

> **诚实边界**：横向扩展方法射程，**不逼近 RH**。墙还在，第二窗口不碰它。

---

## 第一编 · 铁律（方向筛选器，移植自 weil-first）

任何新方向先过这三条（`docs/PROOF_CONSTITUTION.md` PART B/D）：
1. **难度守恒**：任何"让 RH 变简单"的步骤是错觉；合法步骤只无损搬运难度。
2. **禁放缩 (C″)**：RH ⟺ Λ=0，零裕度临界。任何"≥ε 一致下界"论证必错。
3. **叙事抵抗**：支撑精彩叙事的数字要提高验证标准。逐元素比对 artifact 再讲差异。
   区分"过程缺陷"与"结论错误"。

---

## 第二编 · 近期（脚手架 → 第一个探路结论）

| # | 动作 | 状态 |
|---|---|---|
| S1 | 仓库脚手架（README/CLAUDE/PLAN/HANDOFF/骨架目录） | ✅ 完成 |
| S2 | 移植共享机件：archimedean 积分器、interval、ldlt、log_moments、kernel + **单素数极限自检** | ⬜ |
| S3 | 双平移素数层 `legendre_shift_2prime.py`（从 weil-first 原型移入并补全交叉项） | ⬜ |
| S4 | **首要发现动作**：per-sector 素数影响 profile（廉价 mutation 式探针，**certify 级精度**） | ⬜ |
| S5 | 第二窗口 schema + domain + 第一份 pilot 证书骨架 | ⬜ |

> **S4 是最重要的短期动作**（来自 weil-first 实测 steer）：第一窗口 even 扇区素数项
> M2 近乎惰性（归零仅移 pivot 0.003）——真实的"素数 vs 零"张力（RH 关心的）不在 even
> 扇区。故第二窗口**不要对称分配算力**：先 profile 各扇区/各素数/交叉项影响，把硬算力
> 花在真正移动 pivot 处（预期 odd 扇区 / 更大 L / 交叉项 $J(\tau_2,\tau_3)$）。
> 短期测量重定向长期投入。
>
> **精度纪律（来自 weil-first 30× 误差教训）**：S4 的"哪项影响大/小"是算力投向决策，
> 必须 certify 级（Arb interval）确认后才能定论。快扫可做初筛，但任何方向性结论须明确
> 标注精度等级（探路/certify），且有 certify 级数字支撑。探路精度的数字不得下决策
> 判决（PROOF_CONSTITUTION A3）。

> **S2 移植纪律**：
> - 移植来源：`weil-first-prime/src/archimedean/` 的**当前版本**（S_KK-only bug 已于
>   2026-08-08 修复，四项 S0 已验证）。不得携带任何 S_KK-only 残留逻辑。
> - **S2 验收门槛**（不等 S3）：单素数极限自检——weil-second 在 τ₃→0（3 号素数不参与）
>   时，Schur 矩阵各元素必须**逐元素**复现 weil-first 在同 L 点的结果。误差容限：
>   max|C_second − C_first| < 1e-10（float64）。只有通过此检验才标记 S2 为完成。

---

## 第三编 · 中期（有明确路径）

- **数学**：把 Theorem 3 三区间分解验证到 $n=2,3$ 双跳全窗口；给 $\geq3$ 个 $L$ 点
  certify 级严格下界廓线（复用 weil-first 修正后的四项 S0 + 双平移 S2 工具）。
- **proofctl 协同**：交叉项 checker 是否有新的"漏项"盲区？$J(\tau_2,\tau_3)$ 验证需要
  什么新 condition？每暴露一个内核盲区，上游修 proofctl（`~/github/proofctl`）并记入
  方法论论文（weil-first 的 M3《When the Pilot Audits the Tool》的续篇/扩展）。

---

## 第四编 · 硬边界与不可逾越（移植 + 第二窗口特有）

- 结论只能是"$L < \log2$ 第二窗口的有限尺度 Weil 正性"。
- **$L=\log2$ 是硬边界**：$n=2$ 退出单跳区、$n=4$ 进入，Theorem 3 框架需真正扩展。
  **不得外推过 $L=\log2$。**
- **不得**升级为 RH / "接近 RH" / 声称对 RH 有直接推论。
- 任何 certify 断言过 proofctl C01–C11；copy-generator 与漏项 checker 被 C10/C11 封死。

---

## 附 · 代码可信度地图（接手必读，随移植更新）

| 文件 | 状态 | 说明 |
|---|---|---|
| （移植中） | — | S2 完成后填写：archimedean 积分器等来自 weil-first 可信副本 |
| src/prime_layer/legendre_shift_2prime.py | 原型（未认证） | 双平移种子，来自 weil-first 原型；补全交叉项前不得用于证书 |

最危险 bug 模式（继承 + 放大）：漏二阶矩项 **或漏一个素数平移/交叉项** → 残差偏小 →
判据假通过。S0 必须四项；S2 必须含两个素数平移的全部交叉项。

## 附2 · 环境前置

- proofctl 在 `~/github/proofctl`（可修改并发布，需要时先修上游再继续本仓库）；
  `~/bin/proofctl` 为部署副本。github 推送若需代理，用环境变量 `${HTTPS_PROXY:-}`，勿写入文件。
- Python：python-flint(Arb)、numpy；LaTeX 用 tectonic。
- 长任务用 `~/.local/bin/run_and_wait.sh -t <秒> -- <命令>`，前台阻塞，禁裸 `&`。
