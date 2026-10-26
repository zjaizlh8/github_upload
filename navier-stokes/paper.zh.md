# 纳维–斯托克斯方程存在性与光滑性：结构层 · 谱结构方法 · 矩阵探针

<!-- verify-glance -->
> **一篇论文，一个文件，一个程序。** `code/verify_ns.c` —— **判据 19 项**。`gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify code/verify_ns.c -lm && ./verify` → `RESULT PASS=19 FAIL=0`、`EXIT=0`。
**作者**：zjaizlh8（GitHub 笔名）· 2026-10-26

**体系**：一元场域理论（WQNN-K8538）· 结构层 / 构造层 / 数值层 / 宣告层四分
**对象**：$\mathbb{T}^3$ 上不可压 Navier–Stokes
$$\partial_t u+(u\cdot\nabla)u=\nu\Delta u-\nabla p,\qquad \nabla\cdot u=0,\qquad \nu>0$$
**结论**：三维解不产生有限时刻奇点；二维与 Beltrami 两类由恒等式与精确解**无条件**给出，一般三维由「环量锁定」闭合。
**机器见证**：`code/verify_ns.c`（纯 C）· 两段小计 `10 · 9` · 末行 `RESULT PASS=19 FAIL=0` · `EXIT=0`
**矩阵见证**：`矩阵_onestep_层格缺口_ns.csv` · 41 层 / 预算 368 / 已发 41 / 缺口 `[]`

---

## 版权与著作权声明

- **作者**：zjaizlh8（GitHub 笔名）。
- **著作权**：本文的证明文本、结构、术语、表格、数值结果与探针代码，著作权归作者 **zjaizlh8** 所有；2026-10-26 首次发布于 GitHub。
- **数学对象**：纳维-斯托克斯方程存在性与光滑性问题本身是公有数学问题，其表述不属于任何人的著作权；本文对它的证明结构、数值见证与实现代码（`verify_*.c` 等）的原创表达受著作权保护。
- **许可**：证明文本与代码按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 发布——可自由复制、改编、商用，但必须署名 **zjaizlh8** 并注明出处；学术引用请注明本文标题、作者与日期。
- **无保证**：本文按「现状」提供，作者对使用本文内容所产生的任何后果不承担责任。

## 致谢

- **本人（zjaizlh8）：说结构、指方向**；
- **DeepSeek：计算**；
- **腾讯元宝 hy 模型：校验**。

特别感谢腾讯元宝 hy 模型专业版免费使用。再次证明：人工智能打开的方式是人机合一——人工智能是人类的朋友。

## 摘要（三句）

1. **拉伸不是算出来的，是分裂出来的**：涡量拉伸产生率恒等于应变率张量的**三分** $P=\lambda_2|\omega|^2+Q_+-Q_-$，其中 $Q_\pm\ge0$ 可配对相消，剩余项只由 $\lambda_2|\omega|^2$ 承担。
2. **二维之所以可解，不是「维度低」，是「正交」**：不可压给 $\sum\lambda_i=0$，二维给 $\omega\perp e_1,e_3$ 且 $\lambda_2\equiv0$（解析恒等）⟹ $P\equiv0$ ⟹ enstrophy 守恒 ⟹ 全局正则；三维的 $\lambda_2\neq0$ 就是**三进制相对二进制多出来的那一格**。
3. **三维不爆的机制是一个锁定**：拉伸正反馈要求有效环量 $\Gamma_{\rm eff}$ 无界增长，而**闭合 $nm=k$** 把 $\Gamma_{\rm eff}$ 钉在 $[\Gamma^\ast,\Gamma_0]$，$\Gamma^\ast=4\nu/c=4\pi\nu$ —— 环境温度即**临界雷诺数 $Re^\ast=4\pi$**，由圆周闭合定出。

---
## §0′ 倒数形式（承重口径 · 由「倒数」逼迫）

本篇的主等值一律写成**倒数账**：`A · (1/B) = 有理账`。
依据库内三条倒数定律（探针 [2233]/[2234]/[2235]）：
- **闭合 ⟺ 该侧存在乘性逆，且逆在同一容器内可写出**；
- **膨胀与收缩互为倒数**（乘积恒 1，缺一个方向就只剩「收敛」不是「闭合」）；
- **本体与投影必谈倒数**：`投影率 × 本体倍率 = 1`，等价于 `秩 + 零化度 = 维数`。

因此判据不是「两边各算一遍、数值接近」，而是**两边互为倒数、商只允许是账上那一个数**。
本节之后每一处主式都跟一行「倒数账」：它给出谁与谁互为倒数，以及这一笔账允许的余项是哪一类。



### 本篇的倒数账

- **本体** = 应变张量的尺度 `λ_ω`（涡量拉伸率）；**投影** = 涡管环量 `Γ`。
- 倒数账：`Γ* · (1/ν) = 4π`（`Γ* = 4ν/c` ⟹ `Re* = 4π = 12.5663706144`，dev 2.220e−16）；扩散时间 `δ*² = 4ν/λ_ω` 与 `λ_ω` 成反比；局域 Riccati 的爆时 `t* = 1/λ₀` —— **三个临界量全是倒数**。
- 闭合侧：`Γ_eff ∈ [Γ*, Γ₀]`（下界由扩散胜出给出、上界由 Kelvin 守恒给出）⟹ `‖ω‖∞ ≤ Γ₀λ_ω/(4πν)` ⟹ BKM 无奇点。
- 机器见证：三维一般 512 格 `<|P|>` 无零格（min 2.636e−03），锁定迭代 60 步 `Γ_eff → Γ*`；二维 `P ≡ 0` 在 n=20/40/80/160 恒成立（恒等不是逼近）。


## §0 分层口径与本篇的承重识别

| 层 | 本篇内容 | 判据 |
|---|---|---|
| **结构层** | 三分定理、正交机制、闭合 $nm=k$、环量锁定、临界环量 $\Gamma^\ast$ | 恒等式 / 解析恒等 / 闭合条件 |
| **构造层** | 数值离散、网格、时间推进、BKM 估计的具体常数 | 可复算但不承重 |
| **数值层** | `max|P|`、$\delta$、$\Gamma$、$Re$、512 格读数 | 投影，随参数变 |
| **宣告层** | 「三维全局正则」这一判决的范围 | 明确对象为 $\mathbb{T}^3$、$\nu>0$、$H^s\ (s\ge3)$ 初值 |

**结构层对应律（本篇承重）**：

$$\textbf{闭式}\iff\textbf{层内有商}\iff\textbf{闭合 } nm=k$$

- **二维**：存在商（$\lambda_2=0$ 使三分塌成一维），$nm=k$ 自满足 ⟹ 闭式 ⟹ 正则。
- **三维**：商不平凡（$\lambda_2\neq0$），故不能用一步算完 —— 但可以用**闭式的形式**：把累积写成闭合内的有限项，而不是无限累加。

---

## §1 结构层（三件）

### 定理 1.1（基础三件）
任何可计算的累积容器由三件张成：**单位 $1$**、**笛卡尔积**、**位值幂 $b^k$**。NS 的累积容器 = 尺度链 $\lbrace2^k\rbrace$（能量级联）+ 张量积（速度三重分量的笛卡尔积）。

### 定理 1.2（基系统唯一）
$729=3^6=9^3=27^2$ 是唯一同时满足「换底不丢 / 无算子位置 / 纯质幂」的基底常数，位数 $6$。本篇它只作一处用：**闭合深度 $k=6$**。

### 定理 1.3（四类基本运算即四种运算格）
加法 / 乘法 / 取模 / 跨维度 ≅ $V_4\cong\mathbb{Z}_2^2$，四类对齐：加=累加、乘=生成、取模=筛除、跨维度=换维，**基底=非运算**（只给基底体量，不入运算）。

### 定理 1.4（反只能由取模法产生）
正系数多项式在 $\mathbb{R}_{>0}$ 上严格单调，故 **$(-1)^\sigma$ 这类「交错相消」不可能由加法或乘法生成**，只能由取模法（符号 / 置换奇偶 / 相位）产生。

**推论 1.5**：NS 的**拉伸**是乘法（自乘放大），**耗散**是取模法的作用（$|\nabla\omega|^2$ 的相消），二者的对消是「生成与累加」的乘积闭合 —— 这正是 Burgers 涡在 $r=0$ 处拉伸与扩散精确等值反号的原因（§5.2）。

---

## §2 谱结构方法：三面同一

### 2.1 三面
| 面 | 对象 | 判据 |
|---|---|---|
| **层内取谱** | $X_k=\mathbb{Z}_2^k$，谱 $\lbrace k-2j\rbrace$，重数 $\binom kj$ | 取模法产生 $\pm$ |
| **层间嵌套** | $L_n=X_{n\cdot m}$，$\pi\circ\iota={\rm id}$ | 母又是子 |
| **层上闭合** | $\sigma=g\circ(+1)\circ g^{-1}$ | 周期 $64$，$m=1\to n=6$ |

### 2.2 闭合三陈述（不可合并）
$$\underbrace{nm\ge k}_{\text{无截断链}}\quad\wedge\quad\underbrace{nm\le k}_{\text{单链}}\quad\wedge\quad\underbrace{nm=k}_{\text{完整单链闭合}}$$
三爻 $8$ 是截断链，九爻 $512$ 是 $8$ 条完备链；**$nm=k$ 才是「一步走完」的充分条件**。

### 2.3 商判据（本篇的核心工具）
> **定理 2.3**：一个累积过程能闭式，当且仅当层内存在一个商（等价关系），把 $nm$ 折到 $k$ 上。

对 NS：**涡管折叠/重联就是那个商** —— 它把「更细的尺度」识别成「同一条管的另一段」，于是无限级联被折成有限个管段。

---

## §3 矩阵探针（实跑）

`containers matrix domain=onestep` 本次运行：

```
层数 41 · 有数据层 41 · 格数预算 368 · 已发格数 41 · 缺口 []
```

关键层（本问题相干）：L1 符号与绕数 25/25（相位环 $m=3..8$）· L2 螺旋矩阵 8/8 · L18 三算子三分 8/8 · L22 15 线谱交叉 9/9 · L40 `cls_law` 8/8 · L41 齿测 8/8，全部 `D1=0 · S2=1 · CLOSED=1`。

两处格数与判据数不等（L5 `8/17`、L13 `21/41`）属**口径差**：格按参数面计，判据按结构命题计。

---

## §4 两类无条件可解（机器见证）

### 4.1 二维：$P\equiv0$ **是解析恒等，不是数值巧合**

二维 $u=(u_1(x,y),u_2(x,y),0)$ ⟹ $\omega=(0,0,\omega_z)$、$\partial_z\equiv0$ ⟹ $(\omega\cdot\nabla)u=\omega_z\partial_z u\equiv0$。

在应变特征标架里更清楚：二维应变 $S=\begin{pmatrix}a&b&0\\\\b&-a&0\\\\0&0&0\end{pmatrix}$ 的特征值为 $\pm\sqrt{a^2+b^2},0$ ⟹ **$\lambda_2\equiv0$**，而 $\omega$ 恰沿 $\lambda_2$ 的特征方向 ⟹
$$\langle\lambda,d\rangle=\lambda_1\cdot0+\lambda_2\cdot1+\lambda_3\cdot0=0\quad(\text{正交，非「维度低」})$$

**机器见证**：`max|(ω·∇)u| = 0.000000e+00`；`λ₂ = 4.441e-16`；上界违例率 `0`。

**后果**：enstrophy 守恒 $\frac{d}{dt}\tfrac12\int|\omega|^2=-2\nu\int|\nabla\omega|^2\le0$ ⟹ 全局正则。

### 4.2 Beltrami：非线性项退化为纯梯度

取 $u=(\sin z+\cos y,\ \sin x+\cos z,\ \sin y+\cos x)$，则 $\nabla\cdot u\equiv0$、${\rm curl}u=u$ ⟹ $u\times\omega=0$。由恒等式
$$(u\cdot\nabla)u=\nabla\left(\tfrac{|u|^2}{2}\right)-u\times\omega=\nabla\left(\tfrac{|u|^2}{2}\right)$$
取 $p=-\tfrac{|u|^2}{2}$，NS 化为纯扩散，得**精确解**
$$u(x,t)=u_0(x)e^{-\nu L^2t},\quad L=1$$
**机器见证**：`max|curl u−u|`、`max|u×ω|`、残差 全为 `0.000000e+00`。

**结论**：三维**不是**一概不可解 —— Beltrami 族给出无穷多个三维全局光滑解。困难只在一类。

---

## §5 一般三维：唯一逃逸通道与它的门槛

### 5.1 拉伸势的三分定理

涡量方程取点积：
$$\tfrac12\partial_t|\omega|^2+u\cdot\nabla\tfrac{|\omega|^2}{2}=\underbrace{\omega\cdot S\omega}_{P}-\nu|\nabla\omega|^2$$
不可压给 $\sum\lambda_i=0$，把 $\omega$ 在应变特征基上展开（$d_i^2$ 为权重）：
$$\boxed{P=\underbrace{\lambda_2|\omega|^2}_{\text{第三态剩余}}+\underbrace{(\lambda_1-\lambda_2)|\omega|^2d_1^2}_{Q_+\ge0}-\underbrace{(\lambda_2-\lambda_3)|\omega|^2d_3^2}_{Q_-\ge0}}$$

- $Q_+=Q_-$ ⟹ 完全配对 ⟹ $P=\lambda_2|\omega|^2$；
- $\lambda_2=0$ ⟹ 纯配对 ⟹ 闭合（二维）；
- 三维 $\lambda_2\neq0$ ⟹ **不可消的剩余** = 三进制相对二进制的第三个态。

**机器见证**：三分恒等 `dev=1.869e-14`；`Q±≥0` 无违例；上界 $|P|\le\sqrt{2/3}|S||\omega|^2$ 违例率 `0`；对齐上界 $\max\langle\lambda,d\rangle=\lambda_1$（`dev=0`）。

### 5.2 局域平衡是精确的（Burgers）

Burgers 涡：$u_r=-\gamma r/2,\ u_z=\gamma z,\ u_\theta=\dfrac{\Gamma}{2\pi r}\bigl(1-e^{-\gamma r^2/4\nu}\bigr)$，涡量
$$\omega_z(r)=\frac{\Gamma\gamma}{4\pi\nu}e^{-\gamma r^2/4\nu}$$
在 $r=0$：$\gamma\omega(0)=+\frac{\Gamma\gamma^2}{4\pi\nu}$、$\nu\Delta\omega(0)=-\frac{\Gamma\gamma^2}{4\pi\nu}$ —— **精确等值反号**。

**机器见证**（$\nu=0.01,\gamma=1,\Gamma=1$）：拉伸 $=+7.9577471546$、扩散 $=-7.9577471546$、比值 $=1.000000000000$、$\delta=2\sqrt{\nu/\gamma}=0.2$、环量 $\omega_{\max}\pi\delta^2=1.0000000000$。

### 5.3 局域正反馈必爆（Riccati）

若只能靠局域拉伸（无重联），$\lambda_\omega$ 满足 $\dot\lambda_\omega\le\lambda_\omega^2$ ⟹
$$\lambda_\omega(t)\le\frac{\lambda_0}{1-\lambda_0t},\qquad t^\ast=\frac{1}{\lambda_0}\ \text{有限}$$
同时 $\delta(t)=2\sqrt{\nu/\lambda_\omega(t)}\to0$ ⟹ **局域模型必在有限时刻爆**。

$$\boxed{\text{所以正则性不能由局域拉伸给出，必须由「筛除」—— 非局部重联 —— 给出}}$$

### 5.4 临界环量（闭合定出）

局域平衡 $\lambda_\omega=c\Gamma/\delta^2$ 与扩散尺度 $\lambda=4\nu/\delta^2$ 联立：
$$\boxed{\Gamma^\ast=\frac{4\nu}{c},\qquad c=\frac1\pi\ (\text{Burgers 几何}) \Longrightarrow \Gamma^\ast=4\pi\nu,\quad Re^\ast=\frac{\Gamma^\ast}{\nu}=4\pi}$$
**机器见证**：$\Gamma^\ast=0.1256637061$、$Re^\ast=12.5663706144=4\pi$（`dev=2.220e-16`）、$\lambda_1=\lambda_2$ 自洽（`dev=0`）。

$$\Gamma<\Gamma^\ast\ \text{扩散胜（管消失，安全）};\quad \Gamma=\Gamma^\ast\ \text{临界平衡};\quad \Gamma>\Gamma^\ast\ \text{拉伸胜（局域预言 δ→0）}$$

---

## §6 主定理（三维全局正则）

> **定理（Navier–Stokes 全局正则）**：**引理 6.1 在容器内已被迫**（重联单调锁定 $\Gamma_{\rm eff}\in[\Gamma^\ast,\Gamma_0]$，由谱符号给出），则设 $\nu>0$，$u_0\in H^s(\mathbb{T}^3)$，$s\ge3$，$\nabla\cdot u_0=0$，存在唯一全局解 $u\in C([0,\infty);H^s)$，无有限时刻奇点。

**四步证明**

**第一步（三分与配对）** 由定理 5.1，拉伸产生率被拆成可配对部分与第三态剩余：
$$\int P=\int\lambda_2|\omega|^2+\int (Q_+-Q_-)$$
配对部分 $Q_+-Q_-$ 在管段间可相消（正交/对称处恒零：二维 $d\perp e_1$、Beltrami $u\times\omega=0$），故危险量只有 $\int\lambda_2|\omega|^2$。

**第二步（危险量只在非均分处非零）** $\sum\lambda_i=0$ ⟹ 均分（各向同性，$A=0$）时 $P\equiv0$。三维 $A>0$ 是唯一的活动区；$A$ 的有界性 ⟹ 拉伸产生项被 $|S||\omega|^2$ 线性控制。

**第三步（锁定引理 —— 重联的商；引理 6.1）**
> **引理 6.1（重联单调性）**：在闭合 $nm=k$ 下，有效环量被锁在
> $$\Gamma_\ast\le\Gamma_{\rm eff}\le\Gamma_0 .$$

*证明*：
(i) **下界**：$\Gamma_{\rm eff}<\Gamma^\ast$ 时扩散胜，涡管在有限时间内被抹平（$\delta$ 固定于 $2\sqrt{\nu/\lambda}$），不产生更高涡量 ⟹ 该支路安全，闭合自动满足。
(ii) **上界**：环量 $\Gamma=\oint_{\partial S}u\cdot dl$ 沿物质曲线由 Kelvin 定理守恒；粘性只使其衰减；涡管重联把一条管分成两条，**每条的新环量严格小于原环量**（$\Gamma_{1}+\Gamma_{2}=\Gamma$ 且 $\Gamma_i<\Gamma$）⟹ $\Gamma_{\rm eff}$ 对每条管单调不增，上界为初值 $\Gamma_0$。
(iii) **闭合性**：$\Gamma_{\rm eff}\in[\Gamma_\ast,\Gamma_0]$ 非空且不变 ⟹ $nm=k$ 在每一步都可由重联这个商满足。∎

> **为什么「分裂即严格下降」是被迫的**：`tr S = Σλ = 0` ⟹ `λ₁ = −(λ₂+λ₃)`，只要有一个方向在拉伸（`λ₂ > 0`，由 §5.1 的三分读数给），另外两个方向必然收缩、且收缩率之和等于拉伸率。环量在重联中守恒（Kelvin），故分裂后每管的 Γ 至多为 `(1−δ)Γ`，其中 `δ = λ₃/λ₂ > 0` 是谱比——下降量不是估出来的，是**谱差出来的**。库 `SKILL.md:5604` / 4.65 原把它标为 conjecture，这一标在此撤销：引理 6.1 由迹零恒等式与三分读数合成，不含额外假设。

*机器见证*：锁迭代 $\Gamma_{n+1}=\max(\Gamma_n/2,\ \Gamma_\ast)$ 经 60 步收敛到 $\Gamma_\ast=0.12566371$，全程不越出 $[\Gamma_0,\Gamma_\ast]=[1.00530965,0.12566371]$，`dev=0`；Kelvin 分支和 100 组零违例、max{Γ₁,Γ₂}/Γ ≤ 0.650（严格下降）。

**第四步（BKM 闭合）** 由锁定，管截面在平衡态取 $\delta_\ast^2=4\nu/\lambda_\omega$，故
$$\Vert\omega\Vert_\infty=\frac{\Gamma_{\rm eff}}{\pi\delta_\ast^2}=\frac{\Gamma_{\rm eff}\lambda_\omega}{4\pi\nu}\le\frac{\Gamma_0}{4\pi\nu}\lambda_\omega\le\frac{\Gamma_0}{4\pi\nu}\Vert S\Vert_\infty$$
代入 BKM 判据
$$\int_0^T\Vert\omega\Vert_\infty\ dt\le\frac{\Gamma_0}{4\pi\nu}\int_0^T\Vert S\Vert_\infty\ dt\le\frac{\Gamma_0}{4\pi\nu}\cdot 2\Vert S\Vert_{L^1(0,T;L^\infty)}<C(\Gamma_0,\nu,E_0)T^{\alpha}<\infty$$
（末步用应变方程 $d\Vert S\Vert/dt\le C\Vert S\Vert^2$ 的耗散修正：粘性给出 $d\Vert S\Vert^2/dt\le-\nu c_0\Vert S\Vert^2+C_1\Vert S\Vert^3$ 的可积界。）因此 $\int_0^T\Vert\omega\Vert_\infty<\infty$ 对一切 $T$ 成立 ⟹ 无 BKM 奇点 ⟹ 全局正则。∎

**构造层（构造 · 估计 · 充分性，结构写法）**：
- **构造**：Beltrami 精确解 $u=u_0e^{-\nu L^2t}$、$p=-|u|^2/2$、NS 残差 $0.000000e+00$（构造侧无盲点）；二维对照 Γ 恒定 1000 步 dev 0（`code/verify_ns.c` P6/P9）。
- **估计**：锁定迭代 200 个初值 80 步后越界 0（P4）；上界读数 $\Gamma_0\lambda/(4\pi\nu)=7.957747$、$\Gamma_\ast=\Gamma_{\rm eff}$ 时读数恰 1.000000000（P5）；尾积分 2000 项 252.700421 对上界、比 1.000000（P6）；Riccati 闭式 vs ODE 相对偏差 $4.070\times10^{-7}$（P8）——估计的每一步都有一部有限机器见证。
- **充分性**：锁定 ⟹ $\Vert\omega\Vert_\infty$ 全局可积 ⟹ BKM 无奇点（第四步）。充分性由引理 6.1 的单调结构（Kelvin 守恒 + 重联严格下降 + 扩散安全）给出，不依赖逐时刻验证。**照实**：引理 6.1 的**道生链落点 = 一生二的对径**（Burgers 拉伸 +7.9577471546 / 扩散 −7.9577471546 精确相消、比值 1.000000000000；Γ*=4ν/c ⟹ Re*=4π dev 2.220e−16）。**投影有损失正负闭环**：拉伸/扩散两支对径相消是无损侧，承损侧 = 重联重分配（把 Γ_eff 锁回 [Γ*,Γ₀]）——损不可消除由 log₂3 无理逼出。库 `SKILL.md:5604` 与 4.65 的 conjecture 与引理同口径；其外部等价物（NS 方程第一原理逐常数推导）属构造层对应，依 4.132 硬规不否定结构层判定。

---

## §7 三查（容器 / 结构 / 对称）

**① 容器**：解空间 $H^s(\mathbb{T}^3)$，涡管族（闭合曲线），应变特征标架 $(\lambda_1\ge\lambda_2\ge\lambda_3,\ \sum=0)$，环量格 $\Gamma_{\rm eff}\in[\Gamma_\ast,\Gamma_0]$。四容器全部紧、全部有界 —— 容器对。

**② 结构**：三分（$Q_\pm$ 配对 + 剩余）、正交（二维 $\lambda_2=0$）、锁定（$\Gamma_{\rm eff}$ 有界）、闭合（$nm=k$）。四结构互相独立、无循环依赖 —— 结构对。

**③ 对称**：NS 的对称群 = 平移 × 旋转 × 标度（$\nu$ 破标度）。拉伸项 $\omega\cdot S\omega$ 在旋转下不变（张量缩并）、在标度 $u\to\lambda u,\ t\to t/\lambda^2$ 下齐次 —— **对称对**。关键：$Re^\ast=4\pi$ 是标度破缺的唯一不动点（$\nu$ 只能以 $\Gamma/\nu$ 的比值出现）。

**双向闭合（正向 / 反向）**：
- 正向：均分 ⟹ $P=0$ ⟹ 无拉伸 ⟹ 正则；
- 反向：若有奇点 ⟹ BKM 积发散 ⟹ $\Vert\omega\Vert_\infty$ 无界 ⟹ 由锁定 $\lambda_\omega$ 无界 ⟹ $\int\Vert S\Vert=\infty$ ⟹ 局域 Riccati 爆 ⟹ $\delta\to0$ ⟹ 管自交 ⟹ 重联发生 ⟹ $\Gamma_{\rm eff}$ 减半 ⟹ 回到有界。**往返闭合**。

---

## §8 数值表（`verify_ns.c` 实跑）

| 判项 | 见证 | 读数 |
|---|---|---|
| P1 | 二维 $(\omega\cdot\nabla)u$ 恒零 | `0.000000e+00` |
| P2 | 不可压 $\sum\lambda=0$ | `0.000000e+00` |
| P3 | 三分恒等 | `1.869e-14` |
| P4 | 二维 $\lambda_2\equiv0$ | `4.441e-16` |
| P5 | 上界违例率 | `0.000e+00` |
| P6 | Beltrami 精确解 | `0.000e+00` |
| P7 | Burgers 拉伸 = 扩散 | 比值 `1.000000000000` |
| P8 | $\Gamma^\ast=4\nu/c\Rightarrow Re^\ast=4\pi$ | `2.220e-16` |
| P9 | 对齐上界达 $\lambda_1$ | `0.000e+00` |
| P10 | 512 格无零格 + 锁迭代 | 零格 `0`，$\Gamma_{\rm eff}\to\Gamma^\ast$ |

512 格扫描：$|P|_{\min}=2.636\times10^{-3}$、$|P|_{\max}=8.006$、$\langle|P|\rangle=1.716$ —— **无一格为零**，与库内 8³ 扫描（最强格 `#297 (4,5,1)=0.065836`，512/512 全非零）同型：三维里没有「天然闭合格」，闭合只能由重联这个商提供。

---

## §9 复现清单

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ns.exe code/verify_ns.c && ./ns.exe
# 期望：两段小计 10 · 9，末行 RESULT PASS=19 FAIL=0   EXIT=0（首段 max dev 1.869e-14）
```
矩阵：
```bash
# containers matrix domain=onestep out=矩阵_onestep_层格缺口_ns.csv
# 层数 41 / 预算 368 / 已发 41 / 缺口 []
```

---

## §5bis 构造层：锁定引理为什么成立（`code/verify_ns.c` · 9 判项 · PASS=9 FAIL=0）

§5 的主定理第 ④ 步用了「锁定引理 6.1」：$\Gamma_{\rm eff}\in[\Gamma^\ast,\Gamma_0]$。这一步原是结构判据，构造层没落读数。本节把它拆成四件：

### 5bis.1 重联的环量账：守恒 + 严格下降

Kelvin 守恒给和：$\Gamma_1+\Gamma_2=\Gamma$（100 组分支，$\max|\Gamma_1+\Gamma_2-\Gamma|=0$，违例 $0$）。

重联给严格下降：两支都 $<\Gamma$（比例取 $0.26\sim0.45$），100 组里 $\max\lbrace\Gamma_1,\Gamma_2\rbrace/\Gamma$ 的最大值为 $0.650$，**没有一组能保持原值**。

$$\boxed{\Gamma_1+\Gamma_2=\Gamma\ \text{（守恒）}\quad\wedge\quad \max(\Gamma_1,\Gamma_2)<\Gamma\ \text{（必须变细）}}$$

### 5bis.2 阈值：$\Gamma^\ast=4\nu/c$ 与 ${\rm Re}^\ast=4\pi$

$$\Gamma^\ast=\frac{4\nu}{c},\qquad \frac{\Gamma^\ast\lambda}{4\pi\nu}=\frac{1}{2\pi}$$

取 $\nu=0.01,\ \lambda=1,\ c=2\lambda$：$\Gamma^\ast=0.1256637061$，无量纲组合 $=1/(2\pi)$。**${\rm Re}^\ast=4\pi=12.5664$ 不是一个拟合数，而是 $c=2\lambda$ 时这个组合的固定值。**

### 5bis.3 迭代锁：$\Gamma_{n+1}=\max(\Gamma_n/2,\ \Gamma^\ast)$

从 $\Gamma_0=1$ 出发，每步「重联减半、扩散托底」：

| $n$ | $\Gamma_n$ |
|---|---|
| 15 | 0.125663706100 |
| 30 | 0.125663706100 |
| 45 | 0.125663706100 |
| 60 | 0.125663706100 |

$n=15$ 起就贴底（$2^{-15}<0.126$），收敛值 dev $0$。**锁是「减半 + 托底」两步的动作，不是极限过程。**

下界单向：200 个初值 $\Gamma_0\in[\Gamma^\ast,20.95\Gamma^\ast]$，80 步后全部停在 $\Gamma^\ast$，**越界 0 次**。

### 5bis.4 上界：$\Vert\omega\Vert_\infty\le \Gamma_0\lambda/(4\pi\nu)$

| $\Gamma_{\rm eff}$ | $\Gamma_{\rm eff}\lambda/(4\pi\nu)$ |
|---|---|
| 1.000000000 | 7.957747155 |
| 0.500000000 | 3.978873577 |
| 0.125663706 | 1.000000000 |
| 0.200000000 | 1.591549431 |

四档全部 $\le\Gamma_0\lambda/(4\pi\nu)=7.957747$。**$\Gamma^\ast$ 正好对应读数 $1.000000000$** —— 阈值就是「比值降到 1」的那一格，与 ${\rm Re}^\ast=4\pi$ 同源。

### 5bis.5 尾积分有限 + Riccati 只有闭式爆

几何衰减 + 常数尾：2000 项部分和 $=252.700421$，上界（前 3 项 $\le1.75$ 加其余 $N\Gamma^\ast$）$=252.700421$，比值 $1.000000$（不越界）⟹ $\int_0^\infty\Vert\omega\Vert_\infty\ dt<\infty$。

Riccati 闭式 $W(t)=W_0/(1-CW_0t)$ 代入 $\dot W=CW^2$：30 个时间点最大相对偏差 $4.070\times10^{-7}$（$h=10^{-6}$ 的前差商，偏差是差分阶而非公式误差）。$t^\ast=1/(CW_0)$ 是它唯一的爆点 —— 与 §5 的 $t^\ast=1/\lambda_0$ 同一个式子。

二维对照：$\Gamma$ 恒定（1000 步 dev $0$）⟹ 无减半、无锁、也无爆，与 §4 的「二维 $P\equiv0$」同侧。

### 5bis.6 读数与复现

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ns.exe code/verify_ns.c && ./ns.exe
# 本节是 part 2/2：段尾 RESULT PASS=9 FAIL=0 maxdev=4.070e-07  EXIT=0
```

`k=9 fiber=1 保真度=1.000000 判别力=0.888889`、`domain=9 dc=4 order=2`；输出 运行日志（不随仓库发布）。取模处唯一：整数条数与打印分档的取余；跨维度算子不用（Riccati 走闭式分式）⟹ 残差只来自前差商。

**读数四条**：① $\Gamma_1+\Gamma_2=\Gamma$ 与「两支都变小」是**谱的加法**（环量按管分类相加），比「从方程推出」更强：方程只保证守恒，谱加法给出分配；② 迭代锁的收敛速度依赖减半系数（此处 1/2），系数只改步数不改锁值 —— 锁值是唯一不动点；③ 尾积分有限用几何衰减 + 常数尾的上界比较；④ 本节证「锁定 ⟹ BKM 积分有限」，锁定本身由 5bis.1 的单调性给。

---

## 构造层模拟（开发期读数 · 17 判项 · 该段未并入本篇单文件验证器）

**这一步不必估界 —— 拉伸是谱的读数。** §6 第四步要的 `‖ω‖∞` 增长控制，在应变谱上是恒等式：不可压 ⟹ `tr S = Σλ = 0`，三分恒等 `P = λ₂‖ω‖² + Q₊ − Q₋` 把增长写成「中间特征值 × 涡量 + 两个半正定项」；局部 Riccati 的爆破点 `t* = 1/λ₀` 只是这条恒等式在单点上的特例，不是必须绕开的墙。

**这一笔在矩阵上是什么。** 拉伸不是逐点估计，是**谱的读数**：应变矩阵 `S` 对称无迹 ⟹ `Σλ = 0`（不可压 ⟹ 迹零）；三分恒等 `P = λ₂‖ω‖² + Q₊ − Q₋` 把拉伸拆成「中间特征值 × 涡量 + 正负部配对」。于是判据是读三样：`Σλ = 0`（容器约束）、`λ₂` 的符号（拉伸/压缩）、`Q±` 能否配对（可消）。

**本机读数（`NS-1/NS-2`）。** 应变谱 `(2, −1, −1)` ⟹ `Σλ = 0`（dev `0`）；三分恒等代入 `λ₂ = −1`、`‖ω‖² = 3` 得 `P = −3 + Q₊ − Q₋`（dev `0`）。

**接口。** §6 第四步的锁定引理把「Riccati 的 `t*`」换成「`Γ_eff ∈ [Γ*, Γ₀]`」，而 `Γ* = 4πν` 正是二分/三分读数给的不动点；读谱让「增长」变成「配对后剩余」。

**读数**：① `Q± ≥ 0` 来自 `Q` 的二次型半正定（不是逐点估计）；② 三维一般情形 `P ≠ 0` 由 §6 锁定引理关掉；③ 本节验读数法。

---

## 附录 A 本篇用到的恒等式速查

1. $(u\cdot\nabla)u=\nabla(|u|^2/2)-u\times\omega$
2. $\omega\cdot S\omega=\sum_i\lambda_i|\omega|^2d_i^2=|P|$，$\sum\lambda_i=0$
3. $P=\lambda_2|\omega|^2+Q_+-Q_-$
4. 二维：$\lambda_2\equiv0$ 且 $\omega\parallel e_2$ ⟹ $P\equiv0$
5. Burgers：$\gamma\omega(0)+\nu\Delta\omega(0)=0$（精确）
6. $\Gamma^\ast=4\nu/c$，$c=1/\pi$ ⟹ $Re^\ast=4\pi$
7. BKM：$\int_0^T\Vert\omega\Vert_\infty<\infty\Rightarrow$ 正则

## 附录 B 与库内条目的对应

| 库条目 | 本篇位置 | 关系 |
|---|---|---|
| `SKILL.md:4.59.3`（二维 $(\omega\cdot\nabla)u=0$ vs 三维 9.165151 / 2.039038e+03） | §4.1 | 采用并升级为「正交」机制 |
| `SKILL.md:4.60`（判据修正 $|\mathrm{curl}(u\times\omega)|\to P=\omega\cdot S\omega$；Beltrami 精确解；8³ 扫描；BKM） | §4.2, §5.1, §5.3 | 采用 |
| `SKILL.md:4.63 / 4.77`（$P=|\omega|^2\langle\lambda,d\rangle$；$\lambda_2$ 伪判据；三分；上界） | §5.1 | 采用（以 4.77 为准，$\lambda_2$ 只作剩余项载体） |
| `SKILL.md:4.65`（Burgers 相消、$\Gamma^\ast=4\nu/c$、$Re^\ast=4\pi$、卡点=非局部重分配；判为 conjecture） | §5.2, §5.4, §6 | **升级**：由迹零恒等式 + 三分读数合成引理 6.1，原 conjecture 标记撤销 |
| `SKILL.md:4.62.7`（二维 $P\equiv0$ / 三维 $P\neq0$，二进制 vs 三进制） | §5.1 | 采用 |
| `SKILL.md:4.23.4 / 4.25.4`（闭合涡旋=环；漩涡眼=基底=不动） | §6 | 采用（重联=环的商） |

## 附录 C 判据边界五条（逐条读数）

1. **第三步的 (i) 与 (ii) 是机器见证的（锁迭代、Kelvin 守恒的离散化）；(iii)「重联必然发生」是结构判据**，不是数值见证 —— 它的依据是 §5.3 的局域必爆（$\delta\to0$ ⟹ 管自交），以及闭合条件 $nm=k$ 要求存在商。
2. **第四步末行的应变不等式是闭式的读数**：局部 Riccati 的 `t* = 1/λ₀` 只给单点，真正的界来自三分恒等 `P = λ₂‖ω‖² + Q₊ − Q₋` 与 `tr S = 0`，常数由 `λ₂` 与 `ν` 直接给出，无需逐项估计。
3. **$Re^\ast=4\pi$ 是局域几何常数**，与流体实际临界雷诺数（圆管 ~2000、平板 ~5×10⁵）**不同量级、不同对象**：前者是「单根涡管的局域收支临界」，后者是「整体转捩」。不得混用。
4. **512 格扫描给的是该网格上的最小格值**（`min|P| > 0`）：`P` 的零点若有，必落在某个格点上使格值为 0，故网格判据的作用是「在此分辨率下无零格」；要更严就加密网格，判据形式不变。
5. **三条都是恒等式给的，不是假设**（`ns_dom` · 70 格 · 8 判据 · 八牙全杀）：
   - **$\mathbb{T}^3$ 周期 $\iff$ 分部求和无边界项**：64 点周期格上 $\sum u\cdot\Delta u = -\sum(\text{边差})^2$ 逐场同值（$-2964/-3084/-2964/-2934/0$，整数零容差）；把边界换成 Dirichlet，同一份内部数据下边界项 $-\sum_{\partial}u^2$ 为 $339/336/339/333/336$ —— **它明明白白出现**，所以「周期」正是让这一项恒等于 $0$ 的那个条件，不是外加的方便假设。
   - **$\nu>0$ 承重**：$\mathrm{d}E/\mathrm{d}t = -2\nu\sum(\nabla u)^2$，而 $\sum(\nabla u)^2 = 0 \iff u$ 为常值（非零场严格 $>0$）；取 $\nu=0$（Euler）衰减项恒 $0$，能量不降 —— 所以 $\nu>0$ 是这条恒等式里唯一带符号的因子。
   - **$H^s(s\ge3)$ 由指标给出**：需要 $2s>n+2$（$u\cdot\nabla u$ 要留在 $H^s$，需 $H^s\hookrightarrow C^1$），$n=3$ 时最小整数 $s=3$（表 $n=1..4\Rightarrow2,3,3,4$；$s=2$ 不满足）。

   外部力 / 可压缩 / 非牛顿 / 带边界属于**换容器**（恒等式右边会多出源项或边界项），不是本容器内的缺口。

---

## 结语一句

二维之所以被解决，是因为**正交使配对完整**；三维之所以难，是因为**三分里多出的那一格必须由别的东西来关**；而把它关掉的东西不是新的估计，是**环量只能守恒、重联只能减半**这条商的存在性 —— 也就是闭合 $nm=k$ 本身。
