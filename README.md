# Structural Readings of Eight Millennium-Type Problems
# 八道千禧年型问题的结构解读

[![verify](https://github.com/zjaizlh8/github_upload/actions/workflows/verify.yml/badge.svg)](https://github.com/zjaizlh8/github_upload/actions/workflows/verify.yml)
![probes](https://img.shields.io/badge/C%20probes-20-blue)
![gates](https://img.shields.io/badge/archive%20gates-4%2F4-brightgreen)
![problems](https://img.shields.io/badge/problems-8-purple)
[![license](https://img.shields.io/badge/license-CC%20BY%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)

```bash
git clone https://github.com/zjaizlh8/github_upload.git && cd github_upload
make run     # 编译 20 支验证器并逐支实跑 / build + run every verifier
make gates   # 4 道落档门禁扫全部 .c / 4 archive gates over all sources
```
Windows 无 `make` 时：`gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o t <any>/code/verify_*.c -lm && ./t`。
CI（GitHub Actions）在每次 push 上跑同一条链：编译 24 个 `.c` → 跑 20 支验证器 → 4 道门禁。

**直达 / direct:** [站点 Pages](https://zjaizlh8.github.io/github_upload/) ·
[[v1.0.2 Windows x64 预编译 24 个 exe](https://github.com/zjaizlh8/github_upload/releases/latest)](https://github.com/zjaizlh8/github_upload/releases/latest) ·
[CI 运行记录](https://github.com/zjaizlh8/github_upload/actions/workflows/verify.yml)

**Author / 作者:** zjaizlh8 (GitHub) · 2026-10-26
**License / 许可:** CC BY 4.0 (text and code / 文本与代码); conjectures themselves are public mathematical statements / 猜想表述本身为公有数学问题。

A set of eight structural readings and machine-verified studies of the seven Clay Millennium
Problems (plus the twin-prime conjecture), built with one shared framework.
用同一套框架对七道 Clay 千禧年问题（外加孪生素数猜想）做成的八份结构解读与机器验证研究。

- **Structural layer / 结构层** (three base objects: the unit `1`, the Cartesian product, the place-value
  power `b^k`; the unique carrier `729 = 3^6`; the four basic operations ≅ Klein four-group)——
  三件基础对象（单位 `1`、笛卡尔积、位值幂 `b^k`）；唯一基底 `729 = 3^6`；四类基本运算 ≅ 克莱因四元群。
- **Spectral method / 谱结构方法** (three readings of one object: intra-layer spectrum
  `X_k = Z_2^k` with multiplicities `C(k,j)`; inter-layer nesting `L_n = X_{n·m}`, `π∘ι = id`;
  closure condition `nm = k`; the generating law `A(N)`, the closure law `Q(N)`, the recurrence
  `x⁶ − 49x⁴ + 504x² − 1296`, radius set `{2,3,6}`, multiplicity law `c_λ = 2m_λ + [λ=β]`)——
  同一对象的三次读数（层内取谱 `X_k = Z_2^k` 重数 `C(k,j)`；层间嵌套 `L_n = X_{n·m}`、`π∘ι = id`；闭合条件
  `nm = k`；生成函数 `A(N)`、闭包律 `Q(N)`、递推 `x⁶ − 49x⁴ + 504x² − 1296`、半径集 `{2,3,6}`、重数律 `c_λ = 2m_λ + [λ=β]`）。
- **Matrix probe / 矩阵探针** (layer/grid/gap census; the onestep domain runs 41 layers / 368 budget / 0 gaps,
  all `D1=0 · S2=1 · CLOSED=1`)——层/格/缺口普查；onestep 域实跑 41 层 / 预算 368 / 缺口 0，全 `D1=0 · S2=1 · CLOSED=1`。

## Directory layout / 目录布局 (one folder per conjecture / 每个猜想一个目录)

```
hardy-littlewood/   Hardy–Littlewood twin-prime asymptotic   HL 孪生渐近      (code + results + paper.md / paper.zh.md)
twin-prime/         Twin primes are infinite                  孪生素数无穷      (code + results + paper.md / paper.zh.md)
riemann-hypothesis/ All non-trivial zeros lie on σ = 1/2      RH 零点在 σ=1/2  (code + results + paper.md / paper.zh.md)
p-vs-np/            P ≠ NP: the permanent has no closed form  P ≠ NP           (code + results + paper.md / paper.zh.md)
hodge-conjecture/   Hodge classes are algebraic               霍奇类代数        (code + results + paper.md / paper.zh.md)
yang-mills/         Yang–Mills existence and mass gap         杨-米尔斯质量隙   (code + results + paper.md / paper.zh.md)
navier-stokes/      Navier–Stokes existence and smoothness    NS 存在性光滑性   (code + results + paper.md / paper.zh.md)
bsd-conjecture/     Birch–Swinnerton-Dyer ord = rank          BSD 阶=秩        (code + results + paper.md / paper.zh.md)
shared/             shared simulator, run logs, matrix CSVs, registry docs / 共用模拟器、运行日志、矩阵 CSV、登记文档
```

Each `paper.md` is the English version of the corresponding Chinese paper (`paper.zh.md`); both
follow the same four-part skeleton: abstract → structural layer → spectral method → matrix probe
→ main theorem (with the load-bearing step written as an explicit **Lemma + honest notes**) →
construction layer (construct / estimate / sufficiency) → numerical tables → reproduction
commands → honest notes.
每篇 `paper.md` 为英文压缩版（按同一节序压缩，篇幅约中文版 1/5；**不等于逐节全译**——中文版另有 §0′ 倒数形式、§0 记号、§5bis 构造层、§5ter/§5quater、§7bis 独立复算、附录 A/B 若干段未进入英文版，差异见 §Known gaps）：摘要 → 结构层 → 谱结构方法 → 矩阵探针
→ 主定理（承重步显式写成 **引理 + 照实**）→ 构造层（构造 / 估计 / 充分性）→ 数值表 → 复现命令 → 照实。

## Status table / 状态表 (honest grading / 诚实分级)

Every main theorem is stated in **conditional form**: *if the load-bearing lemma holds, then the
theorem holds*. The lemma's status (proved in the archive / framework-caliber / open problem) is
what the honest note records — the theorem head and the honest note no longer contradict each other.
每篇主定理一律以**条件定理**形态陈述：*若承重引理成立，则定理成立*。引理的状态（库内已证 / 体系口径 /
公开问题）由照实记录——定理头与照实不再互相矛盾。

| Paper 篇 | Conditional theorem 条件定理 | Load-bearing lemma 承重引理（前提） | Machine witness 机器见证 |
|---|---|---|---|
| HL | if sign-cancellation thm then `π₂(x) ~ 2C₂·Li₂(x)` 若符号相消定理成立则 … | Möbius random series (probe [412], external open step) Möbius 随机级（探针[412]，外部开放步） | `verify_hl.c` (11) / `verify_hl_tech.c` (9), PASS |
| Twin 孪生 | if Lemma 5.1 then `π₂(x) ≥ ⌊log₂x/6⌋ → ∞` | Lemma 5.1 = weak twin-prime form 引理 5.1 = 孪生猜想弱形式 | `verify_twin_count.c` (6) / `verify_twin_tech.c` (10), PASS |
| RH | if Lemma 4.1 then zeros on σ=1/2 | Lemma 4.1 = unique-radius theorem (framework caliber) 唯一半径定理（体系口径） | `verify_rh.c` (12) / `verify_rh_tech.c` (9), PASS |
| P vs NP | if Lemma 4.1 then P ≠ NP | Lemma 4.1 = Valiant VP≠VNP (open; Pólya m=n closed) Valiant VP≠VNP（开放；Pólya m=n 已闭） | `verify_pvsnp.c` / `verify_pvsnp_tech.c` (9), PASS |
| Hodge 霍奇 | if Lemma 4.1 then (p,p) classes algebraic | Lemma 4.1 = AJ torsion, p≥2 (Griffiths 1969 conflict recorded) AJ 挠性 p≥2（Griffiths 1969 冲突已录） | `verify_hodge.c` (10) / `verify_hodge_tech.c` (9), PASS |
| YM | if Lemma 4.1 then mass gap Δ>0 | Lemma 4.1 = λ₁>0 in 4D SU(3) continuum limit (Clay core) 4D SU(3) 连续极限 λ₁>0（Clay 核心） | `verify_ym.c` (9) / `verify_ym_tech.c` (9), PASS |
| NS | if Lemma 6.1 then no blow-up | Lemma 6.1 = reconnection monotonicity (archive 4.65 conjecture) 重联单调（库内 4.65 conjecture） | `verify_ns.c` (10) / `verify_ns_tech.c` (9), PASS |
| BSD | if Lemma 4.1 then ord = rank | Lemma 4.1 = pairing = zero thickness (= BSD itself for general rank) 配对=零点厚度（一般秩即 BSD 本身） | `verify_bsd.c` (13) / `verify_bsd5.c` (10) / `verify_bsd6.c` (10) / `verify_bsd_tech.c` (9), PASS |

**判据数口径 / counting caliber**: 括号里的数字 = 该支**可失败判据的条数**（不是打印了多少行，也不是
"运行没崩"）。凡只打印不判的项都已在本版改成返回 dev 并计入 `RESULT PASS=n FAIL=m`；`EXIT` 由
`FAIL>0` 或 `dev>阈值` 共同决定，因此**任一判据失配都会被 CI 看见**。
Numbers in parentheses = the count of **failable criteria** per witness (not printed lines, not "it ran").
Every print-only item was converted to return a deviation and enter `RESULT PASS=n FAIL=m`; `EXIT`
depends on `FAIL>0 || dev>threshold`, so any mismatch is visible to CI.

**Two caliber notes / 两处口径说明**:
- `verify_rh.c` = **12**（P1–P12；P5/P11/P12 已从"只打印"改为可失败判据），旧稿曾按 9 记。
- `verify_twin_count.c` = **6**（π₂ 五档对公认值 + 结构律"非 6m 中心者恒 1"），论文旧稿未单列此支计数。
- `verify_simlayer.c` = **17**，是全系列**共用**的模拟器见证，不计入任何单篇的主支数。

**What each witness actually verifies / 每支见证实际验的是什么**（防止把「实现精度」读成「定理成立」）：

| Witness 见证 | Actually verified 实际所验 | Not verified 未验 |
|---|---|---|
| `verify_rh.c` | ζ 的复 Euler–Maclaurin 实现精度（N=1000 · K=10 · t≤50）：`max dev = 1.562e−09 < 1e−08` | σ=1/2 上的零点存在性（Lemma 4.1 仍为体系口径） |
| `verify_pvsnp.c` | GF(2) 上 `rank = (n−1)²+1` 对 n=2..8 全中（17 项）、`bad_id=0`、SAT 0/512、UNSAT 0/65536 | `VP ≠ VNP`（Lemma 4.1 仍是公开问题） |
| `verify_twin_count.c` | π₂(x) 五档与公认值逐一相符 + 结构律「非 6m 中心者恒 1（唯 (3,5)）」 | 孪生猜想本身（Lemma 5.1 = 弱形式） |
| 其余各支 | 条件定理的**前提所涉的有限格**（矩阵恒等式、秩表、常数账、复算表） | 前提在一般情形下的成立性 |

All 20 C probes compile with `gcc -Wall -Wextra -O2 -D_GNU_SOURCE` with **zero warnings**,
pass the four project gates, and run with `EXIT=0`. The four gate sources are bundled in
`shared/code/` (`check_onestep.c`, `chk_loops.c`, `chk_res.c`, `chk_scan.c`), so the gate
check is fully reproducible from this repository:
20 支 C 探针全部 `gcc -Wall -Wextra -O2 -D_GNU_SOURCE` 零告警，过四道门禁，实跑 `EXIT=0`。
四道门禁源码随仓提供于 `shared/code/`（`check_onestep.c`、`chk_loops.c`、`chk_res.c`、`chk_scan.c`），
门禁检查可从本仓库完全复现：

```bash
# build the gates once / 编译门禁（一次即可）
for g in check_onestep chk_loops chk_res chk_scan; do
  gcc -Wall -Wextra -O2 -o /tmp/$g shared/code/$g.c
done
# then run all four on every probe; all must exit 0 / 然后对每支探针跑四道，全部 exit 0
for f in */code/verify_*.c shared/code/verify_simlayer.c; do
  for g in check_onestep chk_loops chk_res chk_scan; do /tmp/$g "$f" || echo "GATE FAIL: $f"; done
done
```

Gate semantics: `main` has no top-level loop, every loop is inside a helper function
(`check_onestep`), every loop chain is inside a judge block (`chk_loops`), no hard-coded
`res=` literal in printf (`chk_res`), no tautological scan (`chk_scan`).
门禁含义：`main` 顶层零循环、循环全在辅助函数内（`check_onestep`）；每条循环链都在判项块内（`chk_loops`）；
printf 无写死的 `res=` 字面量（`chk_res`）；无恒真扫描（`chk_scan`）。

Each conjecture directory holds two versions of its paper:
每个猜想目录存放同一论文的两个版本：

| File 文件 | Language 语言 | Role 角色 |
|---|---|---|
| `paper.md` | English 英文 | The abridged English version written for this repository: same section order, ~1/5 length; **not** a full translation — the sections that only exist in the Chinese version are listed in §Known gaps 本仓库英文压缩版：节序一致、篇幅约 1/5；**不是**全译——仅中文版有的段落列在 §Known gaps |
| `paper.zh.md` | Chinese 中文 | The original paper from the `D:\lw\papers\` archive, kept verbatim `D:\lw\papers\` 档案中的中文原件，一字未动 |

The same pairing applies under `shared/docs/`: `registry.md` / `registry.zh.md` and
`reciprocal-ledger.md` / `reciprocal-ledger.zh.md`.
`shared/docs/` 下同式配对：`registry.md` / `registry.zh.md`、`reciprocal-ledger.md` / `reciprocal-ledger.zh.md`。

> **Language note / 语言说明:** the English documents, file names and the matrix CSVs in this
> repository are in English; the Chinese papers are kept as separate `.zh.md` files (never mixed
> inside an English file). The C sources and the `.out.txt` run logs are kept verbatim as machine
> witnesses (their comments and output labels are in Chinese, faithful to the `D:\dpgo` probe
> archive); translating them would break the one-to-one correspondence with the archive files they
> reproduce. The matrix CSVs' layer names are English translations of the archive's registered
> layer names, listed in the same row order.
> 仓库中的英文文档、文件名与矩阵 CSV 为英文；中文论文以独立 `.zh.md` 文件存放（绝不混入英文文件）。
> C 源码与 `.out.txt` 运行日志按机器见证原件保留（注释与输出标签为中文，忠实对应 `D:\dpgo` 探针档案）；
> 翻译它们会破坏与所复现档案文件的一一对应。矩阵 CSV 的层名为档案登记层名的英文翻译，行序不变。

## Known gaps / 已知缺口

本节是 README 正文里 `§Known gaps` 指向的实体清单：**哪些东西只在中文版里、哪些东西本仓没做**。
This is the entity list that `§Known gaps` refers to elsewhere in this README: what exists only in the
Chinese version, and what this repository does not do.

### 1. 英文版不是中文版的全译 / the English `paper.md` is abridged, not a translation

英文版按同一**节序**压缩，篇幅约中文版的 1/5。下列段落**只在 `paper.zh.md` 里**（节数实测）：

| 篇 Paper | `paper.md` 节数 | `paper.zh.md` 节数 | 仅中文版有的节（部分） |
|---|---|---|---|
| BSD | 11 | 49 | §0′ 倒数形式 · §0″ 倒数三项实证 · §5quater（Reg / c_p / #Ш 全账）· 附录 A · 附录 B |
| Hardy–Littlewood | 12 | 43 | §0′ · §5bis（`E_微` 由什么控制）· §7bis（本篇独立复算）· 附录 A · 附录 B |
| Hodge | 12 | 43 | §0′ · 附录 A · 附录 B（与同系列三篇的分工） |
| Navier–Stokes | 10 | 44 | §0′ · §0 分层口径 · §8 数值表 等 |
| P vs NP / RH / Twin / YM | 10–12 | 40+ | §0′ · §5ter / §5quater · §7bis · 附录 A/B |

⟹ 读英文版时**不要**假设它与中文版逐节等价；需要构造层与附录细节时以 `paper.zh.md` 为准。

### 2. 见证的边界 / what the witnesses do **not** verify

见上方「What each witness actually verifies / 每支见证实际验的是什么」表。三处最容易误读的：

- `verify_rh.c`（12 判据）= ζ 的复 Euler–Maclaurin **实现精度** `max dev = 1.562e−09 < 1e−08`；
  **不**验证 σ=1/2 上零点是否存在（引理 4.1 仍为体系口径）。
- `verify_pvsnp.c`（17 判据）= GF(2) 上 `rank = (n−1)²+1`（n=2..8）、`bad_id=0`、SAT `0/512`、`0/65536`；
  **不**验证 `VP ≠ VNP`。
- `verify_twin_count.c`（6 判据）= π₂(x) 五档对公认值 + 结构律「非 6m 中心者恒 1」；**不**验证孪生猜想。

### 3. 本仓没做的事 / not attempted here

1. **英文版全译**：只补了节序与承重句，未逐节翻译。
2. **跨篇统一定理的陈述**：八条路径共享容器 / 对称 / 闭环三件，但**分别陈述**，没有写出一条统一证明
   （接口见 `shared/docs/registry.md` 的「跨篇统一定理：接口留白」）。
3. **预编译件的版本一致性**：`dist/win64/` 由当前源码编出（24 个 exe + 新校验和，就地 `sha256sum -c` 全通过）；
   `dist/win64-v1.0.1/` 是历史发布快照，其 exe 早于本轮的判据改造（例如其中的 `verify_rh.exe`
   只打印 `RESULT max dev = …`，没有 `RESULT PASS=12`）。要对外发布需另开 tag。
4. **规约级门禁的语义审查**：4 道门禁（`check_onestep` / `chk_loops` / `chk_res` / `chk_scan`）只做**结构**检查
   （`main` 直接层无循环、`res` 不是字面常量、判据不空洞等），**不判断**判据选得对不对、量是不是选错。

## How to reproduce / 如何复现

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe hardy-littlewood/code/verify_hl.c
./hl.exe
```

Each `paper.md` §Reproduction lists the exact commands and expected RESULT lines.
每篇 `paper.md` 的 §复现 列出精确命令与期望的 RESULT 行。

## Honest boundary / 照实边界 (per paper, one sentence each / 每篇一句)

- **HL**: pointwise constant and optimal rate of `E(x)` are not derived here (random-series
  theorem for Möbius weights is the external open step). 逐点常数与 `E(x)` 最优速率未给出（Möbius 权随机级定理是外部开放步）。
- **Twin 孪生**: "each layer keeps a genuine pair" is equivalent to a weak form of the twin-prime
  conjecture; Brun's bound is the true-theorem anchor, the antipodal argument is structural.
  「每层保住一对真孪生」等价孪生猜想弱形式；Brun 界为真定理锚，对径论证为结构论证。
- **RH**: the map "zero set ≅ spectral circle" is a structural isomorphism (framework caliber),
  not a standard theorem. 「零点集 ≅ 谱圆」是结构同构映射（体系口径），非标准定理。
- **P vs NP**: the arithmetic-circuit gap is Valiant's VP≠VNP, which remains an open problem;
  the Pólya m=n case is fully closed in the archive (probes 309–311). 算术电路缺口即 Valiant VP≠VNP，
  仍为开放问题；Pólya m=n 情形在库内完全闭合（探针 309–311）。
- **Hodge 霍奇**: torsionness of AJ images for general classes conflicts with Griffiths (1969) and is
  recorded as such; p=1 needs no lemma (Lefschetz 1924); CM abelian varieties are externally proven.
  一般类的 AJ 像挠性与 Griffiths (1969) 冲突，已如实记录；p=1 不需引理（Lefschetz 1924）；CM 阿贝尔簇有外部证明。
- **YM**: λ₁>0 preserved in the 4D SU(3) continuum limit is the Clay core; the winding-lattice
  completeness is a structural claim. 4D SU(3) 连续极限下 λ₁>0 保持即 Clay 核心；绕数格完备性为结构断言。
- **NS**: reconnection monotonicity is the archive's 4.65 conjecture; the paper proves
  "locking ⇒ BKM closure", not the first-principles derivation of locking.
  重联单调性即库内 4.65 conjecture；本文证「锁定 ⟹ BKM 闭合」，非锁定的一线推导。
- **BSD**: the pairing identification *is* BSD itself for general rank; rank 0/1 is proven
  externally (Gross–Zagier 1986 + Kolyvagin 1990); the constant ledger is realized numerically.
  配对识别对一般秩即 BSD 本身；秩 0/1 有外部证明（Gross–Zagier 1986 + Kolyvagin 1990）；常数账已数值落地。

## Acknowledgments / 致谢

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.
zjaizlh8 说结构、指方向；DeepSeek 执行计算；腾讯元宝 hy 模型参与校验（专业版免费使用）——
人机合一，人工智能是人类的朋友。
