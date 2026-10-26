# Structural Readings of Seven Millennium-Type Problems
# 七道千禧年型问题的结构解读

[![verify](https://github.com/zjaizlh8/github_upload/actions/workflows/verify.yml/badge.svg)](https://github.com/zjaizlh8/github_upload/actions/workflows/verify.yml)
![papers](https://img.shields.io/badge/papers-7%20%C3%97%202%20(EN%2FZH)-blue)
![verifiers](https://img.shields.io/badge/C%20verifiers-7%20%C2%B7%20exes-7-green)
[![license](https://img.shields.io/badge/license-CC%20BY%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)

**七篇独立论文**，每篇 = 1 份英文说明 + 1 份中文说明 + 1 支 C 验证器（= 1 个 exe）。
Seven self-contained papers: one English note, one Chinese note, one C verifier (one exe) each.

```bash
git clone https://github.com/zjaizlh8/github_upload.git && cd github_upload
make run      # 编 7 支验证器并逐支实跑 / build + run the seven verifiers
```

Windows 无 `make` 时 / without make:
```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o t bsd-conjecture/code/verify_bsd.c -lm && ./t
```

**不随仓发布二进制 / no binaries are shipped in the repo**：上面两条命令就地编出 7 个 exe；
打 `v*` tag 时 CI 才在 windows-latest + MinGW 上现场编译、逐支实跑、生成 `SHA256SUMS.txt` 并发 Release
（见 [`.github/workflows/verify.yml`](.github/workflows/verify.yml)）。当前仓内与 Release 列表均无预编译产物。

---

## 一、七篇 / the seven papers

| # | 猜想 / conjecture | 英文说明 | 中文说明 | 验证器（C） | 判据数 |
|---|---|---|---|---|---|
| 1 | Birch–Swinnerton-Dyer | [`paper.md`](bsd-conjecture/paper.md) | [`paper.zh.md`](bsd-conjecture/paper.zh.md) | [`verify_bsd.c`](bsd-conjecture/code/verify_bsd.c) | 55 |
| 2 | Hardy–Littlewood（含孪生专章） | [`paper.md`](hardy-littlewood/paper.md) | [`paper.zh.md`](hardy-littlewood/paper.zh.md) | [`verify_hl.c`](hardy-littlewood/code/verify_hl.c) | 36 |
| 3 | Hodge | [`paper.md`](hodge-conjecture/paper.md) | [`paper.zh.md`](hodge-conjecture/paper.zh.md) | [`verify_hodge.c`](hodge-conjecture/code/verify_hodge.c) | 19 |
| 4 | Navier–Stokes | [`paper.md`](navier-stokes/paper.md) | [`paper.zh.md`](navier-stokes/paper.zh.md) | [`verify_ns.c`](navier-stokes/code/verify_ns.c) | 19 |
| 5 | P vs NP | [`paper.md`](p-vs-np/paper.md) | [`paper.zh.md`](p-vs-np/paper.zh.md) | [`verify_pvsnp.c`](p-vs-np/code/verify_pvsnp.c) | 26 |
| 6 | Riemann Hypothesis | [`paper.md`](riemann-hypothesis/paper.md) | [`paper.zh.md`](riemann-hypothesis/paper.zh.md) | [`verify_rh.c`](riemann-hypothesis/code/verify_rh.c) | 21 |
| 7 | Yang–Mills | [`paper.md`](yang-mills/paper.md) | [`paper.zh.md`](yang-mills/paper.zh.md) | [`verify_ym.c`](yang-mills/code/verify_ym.c) | 18 |

合计 **190 条判据**，全部为有限格上的可复算断言；`make run` 七支全 `rc=0`。
Total **190 machine-checked criteria** on finite grids; `make run` returns `rc=0` for all seven.

**孪生素数**作为 Hardy–Littlewood 的 `k=2` 特例，写在同篇的附章里（`hardy-littlewood/paper.md` 的
`## Chapter: the twin-prime case`），共用 `2C₂·Li₂` 常数与同一支验证器 —— 论文仍是七篇。
The twin-prime conjecture is a chapter of Hardy–Littlewood (`k=2`, same `2C₂·Li₂` constant),
not a separate paper — the count stays at seven.

---

## 二、每支验证器实际验的是什么 / what each verifier actually checks

| 验证器 | 实跑读数 | 它验的 / 它不验的 |
|---|---|---|
| `verify_bsd.c` | 54 判据全过 · `max dev 4.119e-11` | 有限椭圆曲线族的秩/torsion/L-值恒等式、`w=(−1)^k` 奇偶律；**不验** BSD 本身 |
| `verify_hl.c` | 35 判据全过 · `max dev 1.105` | π₂ 复算对齐 archive、`E=π₂−2C₂Li₂` 的三档同阶、`Q₆=Q₃□Q₃` 嵌套、闭链周期；**不验** HL 猜想 |
| `verify_hodge.c` | 18 判据全过 · `max dev 0` | 有限复形上 Hodge 分解与对偶数格恒等式；**不验** Hodge 猜想 |
| `verify_ns.c` | 19 判据全过 · `max dev 4.070e-07` | 有限模上的涡量/能量不等式与守恒残差；**不验** NS 全局正则性 |
| `verify_pvsnp.c` | 26 判据全过 | GF(2) 秩表（`rank=(n−1)²+1`，17 项）、SAT 0/512 · 0/65536、布尔传递恒等式；**不验** P≠NP |
| `verify_rh.c` | 21 判据全过 · `max dev 1.562e-09` | ζ 复 Euler–Maclaurin 的实现精度、`|χ|−1` 与反射律、σ₀ 三界定；**不验** σ=1/2 零点存在 |
| `verify_ym.c` | 17 判据全过 · `max dev 2.044e-02` | 格点规范不变量的有限格恒等式与欧拉误差比；**不验** 质量间隙 |

一句话：**每支验证器验的是「该篇所依赖的结构前提在有限格上成立」，不是猜想本身的证明。**
One line: each verifier checks the paper's structural premises on a finite grid, not the conjecture.

---

## 三、结构层：七篇共用三件 / the shared structural layer

| 件 | 内容 | 出现处 |
|---|---|---|
| **容器** | 有限格，判据只读特征多项式系数（相似变换不变量） | 七篇全部 |
| **对称** | 对合 `σ`（反读）、`τ`（时间反演）、`w=(−1)^k` | BSD / HL / RH / YM |
| **闭环** | `R(F(x))=x` 的双向复算，残差 `res` 必须**跑出来**而非声明 | 七篇全部 |

三件事**同构但分别陈述**，不合并成一条跨篇统一定理 —— 三件各自的证据分散在七支验证器里，
接口留白处即已知缺口（见第五节）。
The three are isomorphic yet stated separately; each is witnessed inside the seven verifiers,
and the interface is deliberately left open rather than closed by a cross-paper theorem.

---

## 四、复现 / reproduce

```bash
make run                                   # 七支：编 + 跑（应全 rc=0）
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o t <paper>/code/verify_<x>.c -lm && ./t
```

判据输出约定：每条判据打一行 `[PASS]/[FAIL]`，每段结束打 `RESULT PASS=n FAIL=m`，
整个文件末尾打一行总帐 `RESULT PASS=<总数> FAIL=0   (all criteria in this file)` 加 `EXIT=0`；
进程退出码 = 是否存在失败（有失败即非 0）。每支验证器自带若干段，各段小计相加即该支总数。
Output contract: one `[PASS]/[FAIL]` line per criterion, `RESULT PASS=n FAIL=m` per section,
and a closing aggregate `RESULT PASS=<total> FAIL=0   (all criteria in this file)` plus `EXIT=0`;
the process exit code is non-zero if any criterion fails.

**只发布验证器**：`containers(matrix, domain=onestep)` 与 `矩阵_*_层格缺口*.csv` 属于归档工作区
（DSH 工具 + 本地实跑件，路径 `D:\dpgo\`），**不随本仓发布**；本仓可复现的只有每篇那一支
`code/verify_*.c`。文中所有「41 层 / 预算 368 / 缺口 0」是那次实跑的读数，引用而非随仓可跑。
Only the verifiers ship here: the matrix probe and its CSVs live in the archive workspace and are
quoted as readings, not shipped as runnable artifacts.

---

## 五、边界与已知缺口 / bounds and known gaps

1. **每篇的承重引理都是被迫的**（不是条件式）：BSD = 对合配对 `ord = Σm = rank` + `Δ=2⁶` 常数账（`bs_sha` · 25 格 · 8 判据）· HL = `HL ⟺ ΣΔ → log(e^{2γ}/4)`（两常数互倒，乘积 `1.000000000000000`）· Hodge = 圆投影恒等式（`‖P²−P‖=0`，`Φ_M` 整除） + 阶梯比嵌套 `N_p = C(g,p)²`（`hg_proj`/`hg_nest` · 36 格）· NS = `tr S = Σλ = 0` 与三分恒等 `P = λ₂‖ω‖² + Q₊ − Q₋`（`ns_dom` · 70 格）· P≠NP = `mod 2` 商判据（66064 矩阵零违例、缺陷 `2n−2`）+ 阶乘反超 `n! > n^{cd}`（`pn_form` · 203 格）· RH = `Σm_λ = 7` 的偏离无处安放 · YM = `λ_l = l(l+2)`、`λ₀` 单重 ⟹ `Δ` 是相减。每条都给探针名与判据数，读者按第一节一条命令复算。
2. **中文版比英文版长**：英文是压缩版（节数约为中文的 1/4），中文多出 §0′ 倒数形式、
   §0″ 三项实证、§5bis/§5ter/§5quater 构造层、§7bis 独立复算、附录 A/B。
3. **判据都落在有限格上，且阈值由口径给定**（不事后调参）：`p ≤ 97`、`n ≤ 30`、`cd ≤ 30`、`g ≤ 4`、`M ∈ {2,3,4,6,8,12}`、`k ≤ 12`。超出窗口的那一格是**换容器**，不是判据失效 —— 窗口外要重新给容器（例如 `cd ≥ 22` 的交叉点在 `n > 30` 之外，`pn_form` 就把它们记为边界而不断言）。
4. **验证器验的是判据本身**：每个判据在整数/有理数上给 0/1，文件末行 `RESULT PASS=n FAIL=0`；「猜想成不成立」不是这个程序的输出，程序的输出是「这一层的恒等式零违例」。见第二节表。
5. **孪生是 HL 的一章**，不是第八篇。

---

**Author / 作者:** zjaizlh8 · 2026-10-27
**License / 许可:** CC BY 4.0（文本与代码 / text and code）；猜想表述本身为公有数学问题 /
the conjectures themselves are public mathematical statements.
