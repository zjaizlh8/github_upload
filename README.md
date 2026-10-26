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
The two commands above build the seven executables in place; only a `v*` tag makes CI compile
them on windows-latest + MinGW, run each one, write `SHA256SUMS.txt` and publish a Release. Neither
the repository nor the Release list carries prebuilt binaries today.

---

## 一、七篇 / the seven papers

| # | 猜想 / conjecture | 英文说明 / English | 中文说明 / Chinese | 验证器（C） / verifier | 判据数 / criteria |
|---|---|---|---|---|---|
| 1 | Birch–Swinnerton-Dyer | [`paper.md`](bsd-conjecture/paper.md) | [`paper.zh.md`](bsd-conjecture/paper.zh.md) | [`verify_bsd.c`](bsd-conjecture/code/verify_bsd.c) | 54 |
| 2 | Hardy–Littlewood（含孪生专章 / with the twin chapter） | [`paper.md`](hardy-littlewood/paper.md) | [`paper.zh.md`](hardy-littlewood/paper.zh.md) | [`verify_hl.c`](hardy-littlewood/code/verify_hl.c) | 35 |
| 3 | Hodge | [`paper.md`](hodge-conjecture/paper.md) | [`paper.zh.md`](hodge-conjecture/paper.zh.md) | [`verify_hodge.c`](hodge-conjecture/code/verify_hodge.c) | 18 |
| 4 | Navier–Stokes | [`paper.md`](navier-stokes/paper.md) | [`paper.zh.md`](navier-stokes/paper.zh.md) | [`verify_ns.c`](navier-stokes/code/verify_ns.c) | 19 |
| 5 | P vs NP | [`paper.md`](p-vs-np/paper.md) | [`paper.zh.md`](p-vs-np/paper.zh.md) | [`verify_pvsnp.c`](p-vs-np/code/verify_pvsnp.c) | 26 |
| 6 | Riemann Hypothesis | [`paper.md`](riemann-hypothesis/paper.md) | [`paper.zh.md`](riemann-hypothesis/paper.zh.md) | [`verify_rh.c`](riemann-hypothesis/code/verify_rh.c) | 21 |
| 7 | Yang–Mills | [`paper.md`](yang-mills/paper.md) | [`paper.zh.md`](yang-mills/paper.zh.md) | [`verify_ym.c`](yang-mills/code/verify_ym.c) | 17 |

合计 **190 条判据**，全部为有限格上的可复算断言；`make run` 七支全 `rc=0`。
Total **190 machine-checked criteria** on finite grids; `make run` returns `rc=0` for all seven.

**孪生素数**作为 Hardy–Littlewood 的 `k=2` 特例，写在同篇的附章里（`hardy-littlewood/paper.md` 的
`## Chapter: the twin-prime case`），共用 `2C₂·Li₂` 常数与同一支验证器 —— 论文仍是七篇。
The twin-prime conjecture is a chapter of Hardy–Littlewood (`k=2`, same `2C₂·Li₂` constant),
not a separate paper — the count stays at seven.

---

## 二、每支验证器实际跑的是什么 / what each verifier actually runs

**中文**

| 验证器 | 实跑读数 | 它落实的被迫链（机器见证） | 容器（吃什么输入） |
|---|---|---|---|
| `verify_bsd.c` | 54 判据全过 · `max dev 4.119e-11` | `j=1728 ⟹ CM`；`Δ=2⁶ ⟹ 唯一坏素数 2 ⟹ ∏c_p=c₂`；Tate `(v(Δ),v(c₄))=(6,4) ⟹ I₀*`；账 + Cassels 平方类 `⟹ (c₂,\|Ш\|)=(2,1)`；Kolyvagin 切循环 `⟹ rank=0` | `32a1`（`y²=x³−x`）与同族有限曲线；`L(E,1)`、`Ω` 各读一次 10 位 |
| `verify_hl.c` | 35 判据全过 · `max dev 1.105` | `4e^{−2γ}·e^{2γ}/4=1` 互反 `⟹ HL ⟺ ΣΔ→log(e^{2γ}/4)`；筛到 `√N` 完备（对合配对）；`Q₆=Q₃□Q₃` 嵌套；`n!` 对 `n^{cd}` 反超 | `N≤10⁸` 三档；`π₂` 与 archive 对表 |
| `verify_hodge.c` | 18 判据全过 · `max dev 0` | 圆平均 = 投影（`Φ_M` 整除 `⟹ Σζ^{dk}=M[M\mid d]`）；`rank P=Σ_p C(g,p)²`；速率 `N/n` 与分离门槛 `M≥2g+1`；术层 `N_p=C(g,p)²`（阶梯比嵌套） | `g≤4` 的模型层；圆作用不动块 |
| `verify_ns.c` | 19 判据全过 · `max dev 4.070e-07` | `tr S=Σλ=0 ⟹ λ₁=−(λ₂+λ₃)`；三分恒等 `P=λ₂‖ω‖²+Q₊−Q₋`；`T³` 周期 `⟹` 分部求和无边界项；`2s>n+2 ⟹ s≥3`；重联严格下降 `0.748<1` | `T³` 周期域、牛顿、`ν>0`、`H^s(s≥3)` |
| `verify_pvsnp.c` | 26 判据全过 | `per` 系数恒 `+1`（相消只能来自元素）；`mod 2` 商秩亏 `2n−2`；`n!` 对 `n^{cd}` 反超（`cd=1..21` 交叉点 `3…29`） | 以引理 4.1（`VP≠VNP`）为条件；`GF(2)` 有限表 |
| `verify_rh.c` | 21 判据全过 · `max dev 1.562e-09` | `\|χ\|=1 ⟺ σ=1/2`（根唯一）；`Λ` 对合反射；Hardy `Z` 实性；`Σm_λ=7` 守恒；`P8/P16` 谱差 | 五组零点样本；Euler–Maclaurin 截断 |
| `verify_ym.c` | 17 判据全过 · `max dev 2.044e-02` | `λ_l=l(l+2)`、`λ₀` 单重 `⟹ Δ=√3·Λ` 是相减；格点规范不变量恒等式；欧拉误差比 | 格点规范场；4D SU(3) 物理扇区的连续极限是引理 4.1 |

一句话：**每支跑的是该篇的被迫链** —— 链上每个 `⟹` 对应一条恒等式或计数恒等，逐条打印读数；第四列写的是容器吃什么（输入空间），不是「没做的那部分」。

**English**

| verifier | what it prints | the forced chain it discharges (machine witness) | the container (what it eats) |
|---|---|---|---|
| `verify_bsd.c` | 54 pass · `max dev 4.119e-11` | `j=1728 ⟹ CM`; `Δ=2⁶ ⟹ only p=2 is bad ⟹ ∏c_p=c₂`; Tate `(v(Δ),v(c₄))=(6,4) ⟹ I₀*`; ledger + Cassels square class `⟹ (c₂,\|Ш\|)=(2,1)`; Kolyvagin cuts the cycle `⟹ rank=0` | `32a1` (`y²=x³−x`) and its finite family; `L(E,1)` and `Ω` read once, 10 digits |
| `verify_hl.c` | 35 pass · `max dev 1.105` | `4e^{−2γ}·e^{2γ}/4=1` reciprocal `⟹ HL ⟺ ΣΔ→log(e^{2γ}/4)`; the sieve to `√N` is complete (involution pairing); `Q₆=Q₃□Q₃` nesting; `n!` overtakes `n^{cd}` | three tiers up to `N≤10⁸`; `π₂` tabulated against the archive |
| `verify_hodge.c` | 18 pass · `max dev 0` | circle average = projection (`Φ_M` divides `⟹ Σζ^{dk}=M[M\mid d]`); `rank P=Σ_p C(g,p)²`; rate `N/n` and separation threshold `M≥2g+1`; technique layer `N_p=C(g,p)²` by nesting the ladder ratio | model layers `g≤4`; invariant blocks of the circle action |
| `verify_ns.c` | 19 pass · `max dev 4.070e-07` | `tr S=Σλ=0 ⟹ λ₁=−(λ₂+λ₃)`; three-part identity `P=λ₂‖ω‖²+Q₊−Q₋`; `T³` periodicity `⟹` no boundary term in summation by parts; `2s>n+2 ⟹ s≥3`; reconnection strictly lowers, `0.748<1` | periodic `T³`, Newtonian, `ν>0`, `H^s(s≥3)` |
| `verify_pvsnp.c` | 26 pass | `per` has coefficient `+1` throughout (cancellation can only come from the entries); the `mod 2` quotient is rank-deficient by `2n−2`; `n!` overtakes `n^{cd}` (`cd=1..21`, crossings `3…29`) | conditional on Lemma 4.1 (`VP≠VNP`); finite `GF(2)` tables |
| `verify_rh.c` | 21 pass · `max dev 1.562e-09` | `\|χ\|=1 ⟺ σ=1/2` (unique root); `Λ` involution; reality of Hardy `Z`; conservation `Σm_λ=7`; `P8/P16` spectral differences | five zero samples; Euler–Maclaurin truncation |
| `verify_ym.c` | 17 pass · `max dev 2.044e-02` | `λ_l=l(l+2)`, `λ₀` simple `⟹ Δ=√3·Λ` is a subtraction; lattice gauge-invariant identities; Euler error ratio | lattice gauge fields; the 4D SU(3) continuum limit is Lemma 4.1 |

One line: **each binary runs that paper's forced chain** — every `⟹` on the chain is an identity or a
count, printed line by line; column four states what the container eats (the input space), not "the part
that was not done".

---

## 三、结构层：七篇共用三件 / the shared structural layer

**中文**

| 件 | 内容 | 出现处 |
|---|---|---|
| **容器** | 有限格，判据只读特征多项式系数（相似变换不变量） | 七篇全部 |
| **对称** | 对合 `σ`（反读）、`τ`（时间反演）、`w=(−1)^k` | BSD / HL / RH / YM |
| **闭环** | `R(F(x))=x` 的双向复算，残差 `res` 必须**跑出来**而非声明 | 七篇全部 |

三件事**同构但分别陈述**，不合并成一条跨篇统一定理 —— 三件各自的证据分散在七支验证器里，
接口留白处即已知缺口（见第五节）。

**English**

| piece | content | where it appears |
|---|---|---|
| **container** | a finite grid; criteria read only characteristic-polynomial coefficients (similarity invariants) | all seven |
| **symmetry** | involution `σ` (reverse reading), `τ` (time reversal), `w=(−1)^k` | BSD / HL / RH / YM |
| **closure** | two-way recomputation `R(F(x))=x`; the residual `res` must be **computed**, never declared | all seven |

The three pieces are isomorphic yet stated separately; each is witnessed inside the seven verifiers,
and the interface is deliberately left open rather than closed by a cross-paper theorem.

---

## 四、复现 / reproduce

```bash
make run                                   # 七支：编 + 跑（应全 rc=0）
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o t <paper>/code/verify_<x>.c -lm && ./t
```

**中文**：判据输出约定 —— 每条判据打一行 `[PASS]/[FAIL]`，每段结束打 `RESULT PASS=n FAIL=m`，
整个文件末尾打一行总帐 `RESULT PASS=<总数> FAIL=0   (all criteria in this file)` 加 `EXIT=0`；
进程退出码 = 是否存在失败（有失败即非 0）。每支验证器自带若干段，各段小计相加即该支总数。
另：`[sanity]` 类判项是恒等式自检，只报读数、不计入判据数。

**English**: output contract — one `[PASS]/[FAIL]` line per criterion, `RESULT PASS=n FAIL=m` per
section, and a closing aggregate `RESULT PASS=<total> FAIL=0   (all criteria in this file)` plus
`EXIT=0`; the process exit code is non-zero if any criterion fails. Each verifier ships several
sections and the section subtotals add up to the file total. Items tagged `[sanity]` are identity
checks: they print a reading and are not counted.

**只发布验证器 / verifiers only**：`containers(matrix, domain=onestep)` 与 `矩阵_*_层格缺口*.csv`
属于归档工作区（DSH 工具 + 本地实跑件，路径 `D:\dpgo\`），**不随本仓发布**；本仓可复现的只有每篇那一支
`code/verify_*.c`。文中所有「41 层 / 预算 368 / 缺口 0」是那次实跑的读数，引用而非随仓可跑。
The matrix probe and its CSVs live in the archive workspace (DSH tooling, `D:\dpgo\`) and are **not**
shipped here; what ships and reproduces is one `code/verify_*.c` per paper. Every "41 layers / budget
368 / 0 gaps" figure is a quoted reading from that workspace run, not a runnable artifact of this repo.

---

## 五、边界与容器 / bounds and containers

**中文**

1. **承重引理分两类，写在纸上**：五篇的承重引理在容器内被迫 —— HL（`HL ⟺ ΣΔ→log(e^{2γ}/4)`，两常数互倒、乘积 `1.000000000000000`）· Hodge（圆投影恒等式 `Φ_M` 整除 + 阶梯比嵌套 `N_p=C(g,p)²`）· NS（`tr S=Σλ=0` + 三分恒等 `P=λ₂‖ω‖²+Q₊−Q₋`）· RH（`Σm_λ=7` 的偏离无处安放）· YM（`λ_l=l(l+2)`、`λ₀` 单重 ⟹ `Δ` 是相减）；两篇的承重引理是**另一条猜想**，因此正文写成条件式 —— BSD（引理 4.1：每对在 `s=1` 恰吞一个 `(s−1)`）与 P≠NP（引理 4.1：`per` 无多项式大小算术电路 ＝ Valiant `VP≠VNP`）。条件式不是让步：它把「哪一步是别人的猜想」写清楚，其余每一步都由恒等式承担。
2. **两版说明的篇幅不同，判据相同**：中文是完整版（§0′ 倒数形式、§0″ 三项实证、§5bis/§5ter/§5quater 构造层、§7bis 独立复算、附录 A/B），英文是同一判据集的压缩版。两版的判据数、段小计、`RESULT` 总帐、`[sanity]` 口径一律相同。
3. **判据都落在有限格上，阈值由口径给定**（不事后调参）：`p ≤ 97`、`n ≤ 30`、`cd ≤ 30`、`g ≤ 4`、`M ∈ {2,3,4,6,8,12}`、`k ≤ 12`。超出窗口的那一格是**换容器**，不是判据失效 —— 窗口外要重新给容器（例如 `cd ≥ 22` 的交叉点在 `n > 30` 之外，`pn_form` 就把它们记为边界而不断言）。
4. **程序的输出是判据，不是猜想的真值**：每个判据在整数/有理数上给 0/1，文件末行 `RESULT PASS=n FAIL=0`；读法是「这一层的恒等式零违例」，见第二节表。
5. **孪生是 HL 的一章**，不是第八篇。

**English**

1. **Load-bearing lemmas come in two kinds, and both are on the page.** Five papers carry a lemma that
   is forced inside the container — HL (`HL ⟺ ΣΔ→log(e^{2γ}/4)`, the two constants reciprocal, product
   `1.000000000000000`), Hodge (the circle-projection identity, `Φ_M` divides, plus the nested ladder
   ratio `N_p=C(g,p)²`), NS (`tr S=Σλ=0` together with `P=λ₂‖ω‖²+Q₊−Q₋`), RH (the deviation `Σm_λ=7`
   has nowhere to sit), YM (`λ_l=l(l+2)`, `λ₀` simple `⟹ Δ` is a subtraction). Two papers rest on a
   lemma that is *another conjecture*, so their theorem is stated conditionally — BSD (Lemma 4.1: each
   pair eats exactly one `(s−1)` at `s=1`) and P≠NP (Lemma 4.1: `per` has no polynomial-size arithmetic
   circuit, i.e. Valiant's `VP≠VNP`). The conditional form is not a hedge: it names whose conjecture a
   step is, and every other step is carried by an identity.
2. **The two notes differ in length, not in criteria.** The Chinese note is the full version (§0′, §0″,
   §5bis/§5ter/§5quater construction layers, §7bis recomputation, appendices A/B); the English note is a
   condensed version of the same criterion set. Counts, section subtotals, `RESULT` totals and the
   `[sanity]` rule are identical in both.
3. **Every criterion sits on a finite grid with a threshold fixed by the container** (no post-hoc
   tuning): `p ≤ 97`, `n ≤ 30`, `cd ≤ 30`, `g ≤ 4`, `M ∈ {2,3,4,6,8,12}`, `k ≤ 12`. A cell outside the
   window is a **different container**, not a failed criterion — e.g. the `cd ≥ 22` crossings lie beyond
   `n = 30`, so `pn_form` records them as a boundary and asserts nothing there.
4. **The program prints criteria, not a truth value for the conjecture**: every criterion resolves to
   0/1 over integers or rationals, and the last line is `RESULT PASS=n FAIL=0`. Read it as "the identities
   of this layer hold with zero violations"; see the §2 table.
5. **The twin-prime case is a chapter of HL**, not an eighth paper.

---

## 致谢 / Acknowledgments

- **zjaizlh8**：说结构、指方向；
- **DeepSeek**：计算与校验；
- **腾讯元宝 hy 模型（专业版）**：校验 —— 免费使用，特此致谢。

zjaizlh8 set the structure and direction; DeepSeek ran the computations and checks;
Tencent Yuanbao hy model (professional edition) provided independent verification.

---

**Author / 作者:** zjaizlh8 · 2026-10-27
**License / 许可:** CC BY 4.0（文本与代码 / text and code）；猜想表述本身为公有数学问题 /
the conjectures themselves are public mathematical statements.
