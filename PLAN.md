# 研究计划 — weil-second-prime 第二素数窗口

**创建日期：2026-08-08**

> 新程序员接手先读 `HANDOFF.md`（操作交接）。本 PLAN 是战略地图。
> 数学纪律见 `docs/PROOF_CONSTITUTION.md`（从 weil-first 移植）。

---

## 状态快照（一句话真相）

**当前真相（2026-08-18 更新）：k=28 “负 frontier 模”是 raw-GL8 截断误差假信号。**
旧 `submatrix_chain.py` 在 $M_K$ 中使用 `skip_remainder=True`，高阶对角项严重失真：
例如 $M_0(P55,P55)$ raw-GL 中心为 $0.422735$，而 Richardson-Arb focused checkpoint
为 $[0.012416,0.012561]$。对整数见证 $v=3P_{53}+P_{55}$，focused Arb 计算给出
$v^TCv\in[+0.2292033466,+0.3953780542]$，整体为正，因此旧链的
$\lambda_0(k28)=-0.181816$ 不能作为负性或 mode 结构证据。该 focused 结果只否定这个
见证的非正定用法，不证明全矩阵正定。旧 `submatrix_k18..28.json` 全部隔离为 legacy
raw-GL discovery；修正链将写入 `submatrix_rich_k*.json` / `submatrix_rich_row*.npz`。
当前 N=25 Arb 尝试继续从逐 pair checkpoint 恢复，N=27 串行排队。

**方法论发现（C13 候选）：pilot-sign firewall。** raw-center 数值可以用于找方向，但只要
截断误差未进入区间，就不能参与符号叙事或模式命名。哪怕文件位于 `pilots/`，也必须携带
remainder mode 与质量等级；否则一次性能优化会把假符号传播到论文草案。`fit_b0.py` 已
显式 warning 并优先读取 corrected `submatrix_rich_k*` 序列；旧 B0/branch/bootstrap
数值全部撤回，不作为证据。

**历史快照（2026-08-15~16）。** 窗口 $L \in (\tfrac12\log3,\ \log2) \approx (0.549, 0.693)$，
两素数 $p=2,3$ 均在单跳区，故素数层是真正的**双平移**对象，带交叉耦合
$J_{ij}(\tau_2,\tau_3)$，$\tau_p=\log p/L$。

**关键发现（2026-08-15~16）：** L=0.56 odd/N=15 certify 在 η=0.5 下确认非正定（C[0][0]≈-5.4e-3）。
主因：kappa(L=0.56)=2.056（比第一窗口 1.255 大 64%）。**N≈34 估计已撤销**（基于错误
F[0,0]≈6.7e-3，实际 0.119，相差 18×）。新发现：(0,0) 近消去结构（R₀[0,0]≈1.08e-3，
R₂/R₀=6.18），η* 关键区间 (0.887, 6.97)；在 η*=2.49 下 C[0,0](N=15)=+3.10e-3>0。
全矩阵 λ_min(C(η*)) 正定性由 eta_scan_N15.py 确认中（~108min，可恢复）。N=19 certify 进行中。

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
| S2 | 移植共享机件：archimedean 积分器、interval、ldlt、log_moments、kernel + **单素数极限自检** | ✅ 完成 |
| S3 | 双平移素数层 `legendre_shift_2prime.py`（从 weil-first 原型移入并补全交叉项） | ✅ 完成 |
| S4 | **首要发现动作**：per-sector 素数影响 profile（廉价 mutation 式探针，**certify 级精度**） | ✅ 完成 |
| S5 | 第二窗口 schema + domain + 第一份 pilot 证书骨架 | ✅ 完成 |
| S6 | proofctl 全链集成：init + graph/policy + 真 generator(C10) + checker mutation(C11) + release gate | ✅ 完成 |

> **S6 结果（2026-08-08）**：claim `thm-second-cross-structure` 经 proofctl
> replay→verify **ACCEPTED**，`release --dry-run` **PASS**（C01–C11 全绿 + 5 个
> metadata key 全绿，无 blocker）。与第一窗口逼出 C10/C11 内核条件不同，第二窗口未逼出
> 新内核条件，但暴露一个 proofctl **集成鲁棒性 bug**：replay 把 attestation metadata
> 反序列化成 `map[string]string`，任一非字符串值（float 100.0）会使整个 metadata map
> 静默丢弃（连带 obligation_results）→ 释放门以 "meta:X not verified" 阻塞。已**上游修复**
> （proofctl commit 86374a7：改 `map[string]json.RawMessage` + 强制转字符串，真解析错误
> 时告警），重建部署二进制至 `~/bin/proofctl`。记入方法论论文续篇。

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
  **C12 候选（2026-08-16 已提议）：** Bernstein 爆炸属于新盲区类——参数默认值在循环变量
  增大后触发超多项式误差增长，C11 的 mutation catalog（漏项检测）无法捕获。已向 proofctl
  上游提议 T37-13（`proofctl/PLAN.md`）：新增 `blowup_parameter_test` catalog 类型，
  覆盖 `integrate_S_KK`/`integrate_S_VK` 的高-k 路径。

---

## 第三编bis · 战略下调(2026-08-08,两份外包判定被独立裁决证伪后)

> 触发:`docs/OUTSOURCED_VERDICTS_PROOFCTL_ADJUDICATION.md`。dominance 与
> margin-collapse 两份外包报告均被独立复现证伪(dominance 的 Layer-1 特征值符号错误,
> Layer-2 一致下界有反例;margin-collapse 净效应猜想被第4节一票否决)。裁决可信,因其
> **先在 S4 已认证锚点上复现 min-pivot=−0.0197 再报分歧**(宪法 A2)。据此重定战略预期。

**第二窗口真实图景(按 grade 严格标注,不重蹈过度自信):**

| 命题 | grade | 状态 |
|---|---|---|
| 交叉项 J(τ₂,τ₃) 非零(度量A) | **certify(Arb,S4)** | **成立**(不受裁决影响) |
| 交叉项"主导"正定(Δλ≈+0.195) | — | **证伪**(真值 +0.046,N=6 处 −0.015) |
| 交叉项"延缓坍缩"(ΔV>0) | pilot(锚定) | **证伪**(三点 ΔV 强负) |
| L=0.55、L=0.6 非正定 | **certify(min-pivot)** | **成立** |
| even/N=8 采样点 λ_min<0 ∈[−0.282,−0.115] | **pilot(锚定,5 点)** | 成立(**非**"全窗口 certify") |
| 全窗口/odd/更大 N 非正定 | — | **未做**,不得断言 |

**成就上限下调:** 从"交叉项主导第二窗口的新数学"→"交叉项是非零的新结构(度量A),
但在测试的 even/N=8 下**不主导**正定性、**不延缓坍缩**"。这是**更弱但诚实可发表**的结果
(结构存在性 + 负结果),不是失败——方法射程的诚实边界本身有价值。

**论文预期(相应下调):**
- "第二窗口交叉项结构":核心 claim(非零,度量A)仍立,但**删去"主导"措辞**,诚实降级为
  "非零新结构 + even/N=8 下不主导 + 该扇区/尺度非正定的负结果"。
- "第二窗口某点正定(FP-second)":**even/N=8 的 5 个采样点均 λ_min<0,此路基本堵死**;
  仅存的活口是 **odd 扇区或更大 N**(裁决诚实标注未扫)——若要复活此方向,先跑
  `scripts/eig_scan_second_window.py --sector odd --N 7` 与更大 N,拿 certify 级数字,
  不得靠乐观。

**监督者元层面记录(2026-08-08):** dominance 初审只审了"方法严格(Weyl+残差界)"就给高
评价,未独立复现数值,被漂亮外壳骗过;独立裁决 agent 复现后抓到 Layer-1 符号错误。教训:
**方法对 ≠ 代码对 ≠ 数值对(宪法 D 系列);监督者本身也不可免于独立复现的仲裁。**

**H1 pilot 进展（2026-08-15）：** odd 扇区 N-收敛扫描（L=0.6，float）：

| N | d | eig_full | Δλ（交叉项贡献） |
|---|---|---|---|
| 3 | 7 | -0.858 | — |
| 5 | 11 | -0.449 | — |
| 7 | 15 | **-0.247** | **+0.042** |

L=0.6 增量递减（0.409→0.202），几何极限约 −0.045，odd/N=7 非正定。

**左端多点扫描（odd/N=7，float）：**

| L | eig_full | Δλ |
|---|---|---|
| 0.56 | -0.167 | +0.018 |
| 0.575 | -0.217 | +0.024 |
| 0.59 | -0.242 | +0.039 |
| 0.60 | -0.247 | +0.042 |

**L=0.56 N-收敛（float pilot，正确 kappa=2.056，2026-08-16 完整版）：**

| N | d | eig_full | 增量 | 增量比 |
|---|---|---|---|---|
| 7 | 15 | -0.5747 | — | — |
| 9 | 19 | -0.3895 | +0.1852 | — |
| 11 | 23 | -0.2770 | +0.1125 | 0.607 |
| 13 | 27 | -0.1879 | +0.0891 | 0.792 |
| 15 | 31 | -0.1155 | +0.0724 | 0.812 |
| 17 | 35 | -0.0760 | +0.0395 | 0.545 |

**certify 结果（Arb 区间，L=0.56，odd，η=0.5）：**

| N | b_L | certify pivot | 结论 |
|---|---|---|---|
| 15 | +0.136 | (0,0)∈[-5.39e-3,-5.38e-3] | **严格非正定**（pivot 严格负）|
| 17 | +0.256 | (2,2)∈[-9.03e-4,+3.27e-3] | **不定**（pivot 跨零，精度不足）|
| 19 | +0.362 | **进行中**（285/361，~85min剩余） | — |

N=17 pivot(2,2) 中心值约 +1.2e-3（所有 η 测试均为正），但区间半宽 ~2e-3 > 中心值：Arb 级是**精度问题，不是符号问题**。

**N=17 Arb eta 快速扫描（从 checkpoint 运行，每次 6s，η=0.5/1.0/2.49/4.0）：**

| η | pivot(2,2) 区间 | 中心 | 下界 |
|---|---|---|---|
| 0.5 | [-9.03e-4,+3.27e-3] | +1.18e-3 | -9.03e-4 |
| **1.0** | [-7.88e-4,+3.59e-3] | **+1.40e-3** | **-7.88e-4（最优）** |
| 2.49 | [-2.33e-3,+4.90e-3] | +1.28e-3 | -2.33e-3 |
| 4.0 | [-4.24e-3,+6.22e-3] | +0.99e-3 | -4.24e-3 |

**N=17 float λ_min 全矩阵 eta 扫描（从 Arb checkpoint 中心值，4s，2026-08-16）：**
‖R₀‖_F=3.43e-2, ‖R₂‖_F=5.11e-2, k=1.491, Frobenius η*=1.221

| η | λ_min(C) |
|---|---|
| 0.5 | -0.07372 |
| **1.0** | **-0.06278（最优）** |
| 1.221（Frob η*） | -0.06580 |
| 2.49（逐元素 η*） | -0.09861 |
| 4.0 | -0.14591 |

**关键诊断：** N=17 全矩阵 λ_min = -0.063（最优 η=1.0），远非正定。
pivot(2,2) 中心正只反映 3×3 领块行为；剩余 Schur 补携带深度负特征值。
逐元素 η*=2.49 对全矩阵反而更差（λ_min 从 -0.074 恶化到 -0.099）。

**方法边界分析（2026-08-16 更新）：** 完整分析见 `docs/method_boundary.md`。
- **重大更正**：早期 "N≈34 / 24h" 估计基于错误的 F[0,0]≈6.7e-3（实际 0.119，相差18×）
- **正确分析（η=0.5）**：(0,0) pivot 要求 N≥17；其他 pivot 更大 N（certify 进行中）
- **关键发现（2026-08-16）：** R₀[0,0] ≈ 1.08e-3（近消去！），R₂[0,0] = 6.65e-3，R₂/R₀=6.18
- **η 关键区间：** C[0,0] > 0 当且仅当 η ∈ (0.887, 6.97)；η=0.5 在范围外，η≥0.9 即可
- **N=15 eta* 结果：** C[0,0](η*=2.49) = +3.10e-3 > 0！(0,0) pivot 在 N=15 已可为正
- **待确认：** 完整矩阵 λ_min(C(η*)) 是否 > 0（scripts/eta_scan_N15.py，~2h，可恢复）
- 更大 L（如 0.60, 0.62）kappa 更大，方法边界更高，**无改善**（已确认）
- T1 pilot 的"正 λ_∞ 外推"（L=0.62, 0.65）均基于 b_L<0 数据，无效已撤销

注：旧表（错误 kappa）已全部废弃。

**L=0.56 N-收敛（float pilot，使用正确 kappa=2.056，2026-08-15 修正）：**

| N | d | b_L（正确 kappa） | 说明 |
|---|---|---|---|
| 7 | 15 | **-0.573** | b_L<0，Schur 判据不适用 |
| 9 | 19 | **-0.343** | b_L<0，Schur 判据不适用 |
| 11 | 23 | **-0.157** | b_L<0，Schur 判据不适用 |
| 13 | 27 | **+0.000214** | b_L 刚刚为正，矩阵几乎肯定非正定 |
| 15 | 31 | **+0.136** | certify 确认 C[0][0]≈-5.4e-3，**非正定** |
| 17 | 35 | +0.256 | 未测试 |
| 19 | 39 | +0.362 | 未测试 |

**阈值分析：** threshold = c_L(0.56) + kappa(0.56) = 1.835 + 2.056 = 3.891
- H(d) > 3.891 当且仅当 d ≥ 27（N ≥ 13）
- 即使 b_L > 0，C = b_L·F - R_eta 仍需正定才算成立

**~~粗估 C[0][0] 转正所需 N（已废弃，2026-08-16）~~：**
- ~~N=15 时 C[0][0] ≈ -5.4e-3，F[0][0] ≈ 6.7e-3（由 b_L 差值反推）~~
- ~~需额外 Δb_L ≈ 5.4e-3/6.7e-3 ≈ 0.80 使 C[0][0] > 0~~
- ~~目标 b_L ≈ 0.136 + 0.80 = 0.936 → H(d) ≈ 4.83 → d ≈ 70 → **N ≈ 34**（奇）~~
- ~~代价：35×35 矩阵 × ~70s/元 × 1225 元 ≈ **24 小时 certify**，不可行~~

**错误根源：F[0,0]≈6.7e-3 是错的（实际 0.119，差 18×）。正确分析见"方法边界分析"节。**

**~~N=13 零点结果~~（已撤销，2026-08-15）：** ~~第二窗口 L=0.56 odd 扇区在 N=13~~
~~首次出现正值 eig_full=+0.000002~~
该"正值"完全由 eig_scan 使用了错误的第一窗口 kappa=1.255 引起，**无效**。
2. ✅ 建立 2-prime Arb certify 路径（`checker/fp_second/certify_fp_second.py`，2026-08-15）：
   - 全路径：Arb 积分器 → Interval M0/S0 → Fraction 点区间素数层 M2/S2 → Schur 组装 → LDL^T
   - 所有常数（c2, c3, c_L、κ、b_L）均通过 Arb→Fraction-Interval 认证
   - 架构：`build_matrices_iv` + `assemble_schur_iv` + `ldlt_factor`
   - 烟测（N=3 odd）通过：全管道无崩溃，LDL^T 正确拒绝（预期，d=7 时 b_L<0）
   - b_L 区间极紧（宽度 ~6e-64）：d=27/N=13 时 b_lo=+0.000214，d=31/N=15 时 b_lo=+0.136
   - 运行方式：`python3 -m checker.fp_second.certify_fp_second --L 5600 10000 --sector odd --N 15`
3. ✗ certify N=15 运行（5721s）结果：**非正定**（第一 pivot C[0][0]∈[-5.39e-3,-5.38e-3]）
   原因：kappa=2.056（正确）→ b_L=0.136（正确）→ C=0.136·F-R_eta 仍非正定
   eig_scan 的"正值"是假阳性（使用了错误 kappa=1.255）
4. ~~下一步分析：评估 L=0.56 odd 正定化所需最小 N（估计 N≈34，计算不可行）；~~
   ~~考虑路径 B（诚实负结果）或 H3（大-N 变分路径）~~
   **已更新（2026-08-16）**：N≈34 估计基于错误 F[0,0]，已废弃。
   实际方向：η* 优化（C[0,0] 在 N=15 已正，全矩阵扫描进行中）+ N=17/19 certify。

**⚠️ 关键 bug 发现（2026-08-15）：eig_scan 用了错误 kappa**

`scripts/eig_scan_second_window.py` 第 39 行导入了第一窗口常数：
```python
from checker.fp035.recompute_schur import ..., KAPPA_FLOAT, ...
# KAPPA_FLOAT = 1.25528305  ← 第一窗口 L=0.35 的 kappa，对第二窗口完全错误！
```
正确值：`compute_kappa(5600, 10000) = 2.0560`（比 KAPPA_FLOAT 大 0.801）

影响：
- 错误 b_L（N=15 odd）= H(31) - 1.835 - 0 - **1.255** = **0.937**
- 正确 b_L（N=15 odd）= H(31) - 1.835 - 0 - **2.056** = **0.136**
- 差值 0.801 × F[0][0] > 0 → C_wrong[0][0] > 0，但 C_correct[0][0] ≈ **-5.4e-3 < 0**

结论：N=13/N=15 的"首次正值 eig_full=+2e-6"是 **虚假信号**，由错误 kappa 引起。
Certify（使用正确 kappa=2.056）已于 2026-08-15 运行 5721s 并正确拒绝：
```
✗ NOT CERTIFIED: Pivot at (0,0) not strictly positive: [-5.390356e-03, -5.382838e-03]
```

**修复（2026-08-15）**：`eig_scan_second_window.py` 改为对每个 L 调用
`compute_kappa(L_num, L_den, prec=128)`，不再使用第一窗口硬编码常数。
已同时启动修正后的 N=7/9/11 扫描（nohup 后台，pilots/eig_scan_corrected_N*.json）。
结果待出后更新本节。

**certify 耗时估算（N=13/N=15 odd）：**
- n=13(N=13) 或 n=15(N=15) 个基函数 → n² = 169/225 个矩阵元
- 每元 M0+S0 约 70s（depth_2d=4, depth_3d=3）→ N=13 约 3.3h，N=15 约 4.4h
- 注意：kappa bug 已修复，使用 compute_kappa 而非 KAPPA_FLOAT

---



- 结论只能是"$L < \log2$ 第二窗口的有限尺度 Weil 正性"。
- **$L=\log2$ 是硬边界**：$n=2$ 退出单跳区、$n=4$ 进入，Theorem 3 框架需真正扩展。
  **不得外推过 $L=\log2$。**
- **不得**升级为 RH / "接近 RH" / 声称对 RH 有直接推论。
- 任何 certify 断言过 proofctl C01–C11；copy-generator 与漏项 checker 被 C10/C11 封死。

---

## 第四编 · HTF 集成路线（新数学成就，2026-08-15 加入）

> **目标：用 HTF 的字符串图 DSL 和 MPO-DMRG 引擎打开两条当前无法通过纯 python-flint
> 路径到达的数学结论。**  
> 不是重新验证已有结论。每一条 H-task 的成功判据是**从未存在过的 certify 级命题**。

### 数学背景（当前活口）

当前已知（certify 级）：
- even/N=8，五个采样 L 点：$\lambda_\min \in [-0.282, -0.115]$，**全部非正定**。

当前完全未做：
- **odd 扇区任意 N**（PLAN 第三编bis 明确标注"未做，不得断言"）。
- even/odd 扇区 **N≥12** 的任意 L 点（dense Arb 在 N=12 的复杂度为 $O(N^3)$ × Arb 开销，
  单点约 3–5 小时；N=16 不可行）。

两条活口对应两种新数学结论：
1. **Odd 扇区正定性**：若存在 $L$ 使 odd/N≥7 的 min-pivot > 0（certify），则首次得到
   "第二窗口存在正定点"——这是 FP-second 的唯一残存可能。
2. **大 N even 扇区翻转**：若 N=12 even 扇区在某 $L$ 点 min-pivot > 0（certify），则说明
   N=8 的非正定是有限尺度效应而非本质障碍——同样是新的结论。

---

### H1 · Odd 扇区首次 certify 扫描（**已执行，结果为负，2026-08-15**）

**执行结果（N=15 odd，L=0.56，certify 级，5721s）：**
- b_L(d=31) = +0.136（正确 kappa=2.056）> 0 ✓
- LDL^T 第一 pivot C[0][0] ∈ [-5.39e-3, -5.38e-3] < 0 → **非正定，认证失败**
- 根本原因：kappa(L=0.56) = 2.056 远大于第一窗口 1.255，b_L 仍不足以使 C ≻ 0

**kappa 问题的数学含义：**
- 对所有第二窗口 L ∈ (0.549, 0.693)，kappa ≥ 2.017（远大于第一窗口 1.255）
- 阈值 H_d > c_L + kappa ≥ 3.83，要求 N ≥ 13（odd）才能 b_L > 0
- 即使 b_L > 0，C 正定还需 b_L 足够大：粗估 L=0.56 需 N ≈ 34（~24h certify，不可行）

**H1 路径 B 结论（2026-08-15）：**
| 步骤 | 内容 | 结论 |
|---|---|---|
| H1a | float 扫描 odd/N=7 | ✅（b_L < 0 于 N≤11，正值只来自错误 kappa） |
| H1b | certify odd/N=15，L=0.56 | ✅（失败，C[0][0]≈-5.4e-3） |
| H1c | min-pivot > 0 → 发布 claim | ❌ 未达成 |
| H1d | 所有点 min-pivot < 0 → 诚实记录 | **✅ 确认** |

**触发 H3（大-N 路径）或接受路径 B（诚实负结果）。**

---

### H2 · J_ij(τ₂,τ₃) → HTF Box 张量化

**目的：** 将交叉耦合矩阵 $J_{ij}(\tau_2,\tau_3)$ 表达为 HTF 的 `Box` 对象，使其进入
HTF 的类型系统和认证引擎，为 H4 的全管道字符串图铺路。

**数学映射：**
$$J_{ij}(\tau_2,\tau_3) \quad\longrightarrow\quad \texttt{Box}("J\_cross",\; (\texttt{Wire}("i", d_i),),\; (\texttt{Wire}("j", d_j),))$$
其中 $d_i = d_j = $ 当前矩阵维度（由扇区 N 决定），Arb interval 条目由现有
`legendre_shift_2prime.compute_F` 提供。

| 步骤 | 内容 |
|---|---|
| H2a | 在 `src/htf_bridge/` 下新建 `j_cross_box.py`：接受 τ₂, τ₃（Arb），输出 HTF `Box`，内部 tensor 为 Arb interval 矩阵。 |
| H2b | 验收：`contract(J_box, F, mode="certified")` 输出的区间中心与 `compute_F` 结果逐元素一致（max\|center diff\| < 1e-10），区间半径 ≤ python-flint 直接计算的 2 倍。 |
| H2c | 补充 `tests/htf_bridge/test_j_cross_box.py`：单素数极限（τ₃→0）时退化为 `legendre_shift.py` 的 J（与 S2 自检一致的 pattern）。 |

**注意：** H2 不改变任何现有数值路径，纯新增桥接层。若 HTF 桥接出现精度损失，以
python-flint 为 ground truth，不得反向降精度。

---

### H3 · 大-N Schur 矩阵 via HTF MPO-DMRG（打开 N≥12 路径）

**目的：** 当 N=12,16 时 dense Arb 不可行，HTF 的 MPO + 2-site DMRG 提供两侧 variational
界（variational upper bound + Temple lower bound），使大-N 结论成为可能。

**数学结构：**
Schur 矩阵 $C$ 是一个对称正半定算子候选，可以写成 MPO 形式（每个矩阵元是算子积
的迹，适合 1D bond-dimension 截断）。HTF 的 `htf.mpo.MPO` + `htf.mpo.MPODMRGSolver`
对 $C$ 做变分最小化，给出 $\lambda_\min(C)$ 的变分上界；结合 Temple 界（现有
`htf.gap` 机制）得到两侧区间。

| 步骤 | 内容 |
|---|---|
| H3a | **可行性探针（pilot）**：N=8 even 扇区，同时跑 dense Arb 和 HTF DMRG，比较 $\lambda_\min$ 估计。若 DMRG 结果与 Arb 一致（相对误差 < 1%），则 H3 路径可信。 |
| H3b | N=10,12 even 扇区，HTF DMRG 变分估计 $\lambda_\min$（float 精度），确认 N→∞ 趋势（递增还是递减）。 |
| H3c | 若 H3b 显示 N=12 某 L 点 $\lambda_\min > 0$（float），则用 Temple 界（HTF certified 模式）给出两侧认证界，判断是否足够给出 certify 级结论。 |
| H3d | 若 H3c 成功：走 proofctl 链，发布 `thm-second-even-large-N`。 |

**验收门（H3a）：** DMRG 在 N=8 的估计与 Arb 参考值相差 < 1%；bond dimension χ=32
足够（χ 收敛性测试：χ=16 vs χ=32 差异 < 0.1%）。

**H3a pilot 结果（2026-08-16）：** `rayleigh_certificate(C, psi_min)` 在 L=0.56 odd N=5 验证通过：
- HTF Rayleigh midpoint = numpy λ_min（相对误差 = 0，radius=1.11e-16）
- H3a 验收门 **PASS**（见 `pilots/h3_variational_pilot.json`）
- HTF 路径可行；下一步是在 N=8 even 使用 MPO/DMRG 对 HTF dense-Rayleigh 做标定
**诚实边界：** Temple 界在 HTF 中标注 `[heuristic]`（见 HTF CLAUDE.md §4）——
进入 certify 链前必须确认是否升级为 `[engineering]`，不得以 heuristic 界发布 claim。

---

### H4 · 全 Schur 管道字符串图（组合认证层）

**目的：** 将完整计算管道表达为 HTF 字符串图组合，使数值正确性可以由**字符串图合成
定律**推导，而非仅靠独立 checker 抽查。这是方法论贡献，不是重新验证结论。

**管道拓扑：**
```
[Archimedean_A] ─┐
[Archimedean_B] ─┤─→ [Assemble_S0] ─→ [Prime_Layer_J] ─→ [LDLT] ─→ min-pivot
[J_cross_box]  ─┘
```
每个节点是一个 HTF `Box`；`>>` 组合在 composition time 做 wire-dimension 类型检查；
`contract(..., mode="certified")` 给出端到端的 Arb 区间。

| 步骤 | 内容 |
|---|---|
| H4a | 将现有 `src/archimedean/integrator_a.py` 包装为 HTF Box（`archimedean_a_box.py`），接受 L（Arb），输出 $A$ 矩阵（Arb interval tensor）。 |
| H4b | 将 `src/assemble/recompute_schur_2prime.py` 的装配逻辑重表达为字符串图的 `>>` 组合，输出 $C$ 矩阵 Box。 |
| H4c | 将 `src/archimedean/ldlt.py` 包装为 HTF Box（输入对称矩阵，输出 min-pivot Arb ball）。 |
| H4d | 端到端测试：字符串图管道在 L=0.6/N=7/odd 的输出与现有 python 路径逐元素一致（验收同 H2b 精度标准）。 |

**注意：** H4 是增量新增，不删除现有 python 路径。两路并行运行，以现有 python-flint
路径为 ground truth。若存在分歧，以 ground truth 为准并查 HTF 桥接 bug。

---

### H5 · 成就目标与论文映射

| 成就 | 触发条件 | grade |
|---|---|---|
| **FP-second：第二窗口存在 certify 正定点** | H1c 成功 OR H3d 成功 | certify(Arb) |
| **全扇区负结果：even+odd N≤12 非正定** | H1d + H3b 均为负 | certify(Arb) |
| **双素数耦合张量字符串图框架** | H2+H4 完成 | engineering |
| **大-N 变分路径开通** | H3a+H3b 完成 | research（需 H3c 升 certify） |

论文方向（按成就 landing 情况选择）：
- 路径 A（H1c 成功）：*"Second Prime Window Weil Positivity: Certificate-First Proof via Odd-Sector Split-Residual Schur"*
- 路径 B（全负）：*"Finite-Scale Non-Positivity in the Second Prime Window: Honest Boundary of the Split-Residual Method"*（负结果，方法射程诚实边界，仍可发表）
- 路径 C（H3d 成功）：*"Large-N Weil Positivity via Tensor-Network Variational Certification"*（新方法贡献）

> **诚实约束（任何路径）**：结论边界 $L < \log 2$；不声称 RH 或 RH 推论；
> 负结果不等于"第二窗口不存在正定点"（只是"当前方法/尺度未找到"）。

---



| 文件 | 状态 | 说明 |
|---|---|---|
| src/archimedean/{integrator_a,integrator_b,interval,ldlt,log_moments,kernel}.py | ✅ 可信（移植自 weil-first 当前版，四项 S0 已验证） | 逐字移植，tests 通过。 |
| src/archimedean/bernstein.py | ✅ 已修复 P0 bug（2026-08-15）| `bernstein_mk_bound` 旧版 M_f=7/4·a 忽略 Bernstein 椭圆上 P_n(z) 的指数增长；修正后 M_f 包含 (2R)^{n_row} 因子。N=15 负向 certify 有效（bound=8e-8）。Bernstein 可用至 N≤21（n_row≤41，bound=0.31<1）。 |
| src/archimedean/integrator_a.py | ✅ 三项扩展（2026-08-15）| (1) x 方向自适应 subdivision（n_xsub 对称于 n_ysub）；(2) Richardson GL-8/GL-4 mode（use_bernstein=False）：2\*|GL8-GL4|，为 N≥22 提供余量覆盖；(3) skip_remainder=True：float-only path 跳过所有余量计算；(4) GL-8 节点移至 kx 循环外（n_sub 次→1 次）。135 tests 通过。|
| checker/archimedean/{integrate,check_archimedean,replay}.py | ✅ 可信（移植） | 共享 archimedean checker 机件。 |
| checker/fp035/recompute_schur.py | ✅ 可信（移植，作单素数 ground truth）| 正确四项 S0 + min-pivot。skip_remainder=True（float-only path）。 |
| checker/fp_second/certify_fp_second.py | ✅ CLAUDE.md 长任务规范 + --no-bernstein（2026-08-15）| checkpoint/resume/observable/incremental-durable；--no-bernstein 启用 Richardson 余量（N≥22 必需）；use_bernstein 贯穿 build_matrices_iv→integrate_M_K。 |
| src/prime_layer/legendre_shift.py | ✅ 可信（移植） | 单素数 J/E，Fraction 精确算术，25 tests 通过。 |
| src/prime_layer/legendre_shift_2prime.py | ✅ 可信（S3 完成：M2 双平移 + S2 交叉项 F 完整） | `compute_F` 四条独立不变量验证：F(τ,τ)=E、i+j 奇→0、算子交换对称、cross-Cauchy-Schwarz。**无 F=0 静默漏项**（C11）。c3≠0 且无 F_provider 仍 raise。 |
| src/assemble/recompute_schur_2prime.py | 🟡 pilot（float center，非 certify） | 复用 weil-first 四项 S0 装配，仅换素数层。含 S4 profile 开关（include_tau2/3/cross, swap_c2_c3），默认全开。verdict 需 S5 的 Arb certify 路径。 |
| scripts/single_prime_limit_check.py + tests/prime_layer/test_single_prime_limit.py | ✅ 可信 | S2 验收门：c3=0 时两素数层逐元素复现 weil-first，max\|dC\|=0（精确）。 |

> **关键结构事实（2026-08-08，S3 期间发现）**：判据需 b_L = H_d − c_L − L0 − κ > 0，
> 其中 d 是首补自由度 = 2N（偶）/2N+1（奇），**非自由小整数**。第二窗口 c_L≈1.82–2.04
> 比第一窗口（1.365）大，阈值 H_d > c_L+κ 要求 d≥12（L=0.55）到 d≥15（L=0.69），
> 即 **N≥6–8**。这是 handoff "第二窗口积分更重" 警告的具体形式。S4/S5 的 certify 级
> 跑必须 N 足够大使 b_L>0，绝不能用 d=1 快扫下判决（tiny-N 的 pivot 符号是无意义的
> plumbing）。

> **S2 验收结果（2026-08-08）**：单素数极限自检通过。L=0.6 两扇区 N=3，
> 完整装配 Schur 矩阵 C 逐元素 max\|C_second − C_first\| = 0.00e+00（精确复现，
> 因 c3=0 路径重算同一 J/E）。prime-layer 级 + assembled-C 级双重验证。

最危险 bug 模式（继承 + 放大）：漏二阶矩项 **或漏一个素数平移/交叉项** → 残差偏小 →
判据假通过。S0 必须四项；S2 必须含两个素数平移的全部交叉项。

## 附2 · 环境前置

- proofctl 在 `~/github/proofctl`（可修改并发布，需要时先修上游再继续本仓库）；
  `~/bin/proofctl` 为部署副本。github 推送若需代理，用环境变量 `${HTTPS_PROXY:-}`，勿写入文件。
- Python：python-flint(Arb)、numpy；LaTeX 用 tectonic。
- 长任务用 `~/.local/bin/run_and_wait.sh -t <秒> -- <命令>`，前台阻塞，禁裸 `&`。
