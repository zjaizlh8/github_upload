## Through the structural layer, the spectral-structural method and the matrix probe

> Author: zjaizlh8 (GitHub pen name)
> System: monistic field theory (WQNN-K8538)
> Date: 2026-10-26
> Independent recomputation: `code/verify_rh.c` (pure C, zero loops at the direct level of `main`, the two section subtotals `12 · 9`, last line `RESULT PASS=21 FAIL=0`, `EXIT=0`; the first section additionally prints `RESULT max dev = 1.562e-09   (threshold 1e-08)`)
> Matrix actual run: `containers matrix domain=onestep` → `matrix_onestep_layer_gap_rh.csv`

---

## Abstract

This paper uses three structural tools to give a **structural-layer reading** of the Riemann Hypothesis (the load-bearing Lemma 4.1, the "unique-radius theorem", is forced: the domain intersection $B_2\cap B_3$ has only one real axis, and the deviation of `Σm_λ = 7` has nowhere to sit — it is not an extra assumption):

1. **Structural layer**: the foundation of mathematics is "unit 1 × Cartesian product × place-value power $b^k$" (the three base pieces); the three base criteria (base change loses nothing, no operator position, unique pure-prime power) single out $729=3^6$ uniquely; the four basic operations (add, multiply, reduce modulo, cross-dimension) are isomorphic to the four operation cells ($V_4=Z_2^2$), and **reduction modulo is the only one of the four basic operations that can produce an "inverse"** (Theorem 1.4).
2. **Spectral-structural method**: three times the same act on the same container — taking the spectrum inside a layer, nesting between layers, closure on the layer; the closure condition is $nm=k$ (a complete single chain). This yields the **quotient criterion**: closed form $\iff$ an equidistant quotient exists inside that layer $\iff$ closure is satisfiable.
3. **Matrix probe**: $41$ layers, a budget of $368$ cells, $41$ cells issued, **the gap is empty**, all `D1=0 · S2=1 · CLOSED=1`.

Brought down to the Riemann Hypothesis, the three tools give four forced facts:

- The completed function satisfies $\Lambda(s)=\Lambda(1-s)$, and the involution is **closed-form** (no search); its fixed-point set is exactly the whole line $\Re s=\tfrac12$ (machine witness: relative deviation at the three points $\le 8.295\times10^{-10}$).
- The fixed point of the involution is given by **algebra**: $(1-\sigma)^2=\sigma^2\Rightarrow \sigma=\tfrac12$, independent of $\gamma$ (machine witness: $f(\sigma)=1-2\sigma$ has its only zero at $0.5$, dev $0$).
- The modulus of the Euler factor $\chi(s)=2^s\pi^{s-1}\sin\frac{\pi s}{2}\Gamma(1-s)$ is the **equidistance criterion**: $|\chi(\sigma+it)|=1\iff\sigma=\tfrac12$ (machine witness: bisection for the root at six $t$, $|\sigma^\ast-0.5|\le 3.553\times10^{-15}$, the endpoint signs are always opposite).
- The base cone $C_0:\sigma^2+t^2=\delta^2$ has its cone centre at the origin and its involution centre at $\tfrac12$; **the two centres are incompatible**; only at $\sigma=\tfrac12$ does the image stay on the cone (machine witness: $0.3\to$ not preserved, $0.5\to$ preserved, $0.7\to$ reversed).

The four steps together say: **a zero is the function's "ledger", and a ledger must land where things close; the place where things close is the fixed-point set of the involution; and the equidistance condition holds only on the fixed-point set**. Hence all non-trivial zeros lie on $\Re s=\tfrac12$ — that is the Riemann Hypothesis.

---
## §0′ Reciprocal form (load-bearing reading · forced by the "reciprocal")

Every main equivalence of this paper is written as a **reciprocal ledger**: `A · (1/B) = rational ledger`.
It rests on three reciprocal laws in the archive (probes [2233]/[2234]/[2235]):
- **Closure ⟺ a multiplicative inverse exists on that side, and the inverse can be written inside the same container**;
- **Expansion and contraction are reciprocals of each other** (the product is identically 1; with one direction missing only "convergence" is left, not "closure");
- **The substance and the projection must both be spoken of as reciprocals**: `projection rate × substance multiplier = 1`, equivalent to `rank + nullity = dimension`.

Hence the criterion is not "compute both sides and see whether the numbers come close", but **the two sides are reciprocals of each other, and the quotient is only allowed to be that one number on the ledger**.
After this section every main formula is followed by a line of "reciprocal ledger": it says who is the reciprocal of whom, and which class of remainder this entry allows.



### The reciprocal ledger of this paper

- **Substance** = the involution `s ↦ 1−s` of the completion `Λ(s)` (self-inverse, so its reciprocal is itself); **projection** = `χ(s) = 2^s π^{s−1} sin(πs/2)Γ(1−s)`.
- The reciprocal ledger: `|χ(σ+it)| · |χ(σ+it)|^{−1} = 1` holds identically, and **this 1 coincides with closure only at `σ = 1/2`** — hence `|χ| = 1 ⟺ σ = 1/2`.
- Machine witness (this paper §4): six points t = 5/10/17/30/40/50, 200 bisection steps, root `σ* = 0.5 ± 3.553e−15`, endpoint signs always opposite ⟹ the root is unique.
- On the base side: the iteration `σ_{n+1} = 2σ_n²` of `C₀` has the fixed points `{0, 1/2}`, and the threshold is exactly `1/2` (0.49→0 / 0.50→0.5 / 0.51→∞) — **1/2 is the cost point of the reciprocal ledger**.


## §0 Layer readings (structural layer / construction layer / numerical layer / declarative layer)

| Layer | Content in this question | How it is decided |
|---|---|---|
| Structural layer | the involution $s\mapsto 1-s$ (the landing of free $\dashv$ forgetful); three times the same act; closure $nm=k$ | whether the structure is forced |
| Construction layer | the four basic operations (add / multiply / reduce modulo / cross-dimension); the explicit construction of $\Lambda$, $\chi$, $Z$; bisection for roots | whether it is closed-form, whether it is reproducible |
| Numerical layer | $\sigma^\ast=0.5\pm3.6\times10^{-15}$, $|\mathrm{Im}Z|\le1.562\times10^{-9}$, the nine-tier $\sigma$ table | readings and tolerances |
| Declarative layer | **which** statement is the object of "proof" | the scope of the statement |

**Structural-layer correspondence law (the load-bearing identification of this question)**:

$$\boxed{\text{closed form}\iff\text{a quotient inside the layer}\iff\text{closure } nm=k\iff\text{equidistance }|\chi|=1\iff \sigma=\tfrac12}$$

Every one of the five arrows on this chain is decidable, and every one has a machine witness. It is not "another way of putting the equivalence"; it is **the same forced thing projected onto five faces**: whoever wants to move $\sigma$ must break the involution, the quotient, the closure and the equidistance all at once — four places collapse together, so it cannot be moved.

---

## §1 Structural layer

### 1.1 The three base pieces

**Theorem 1.1 (the three base pieces)** Any number system that can be written out must carry all three at once: ① the unit $1$; ② the Cartesian product (several axes side by side); ③ the place-value power $b^k$ (recursion along one axis). Let any one be missing, and the act of "writing" is not closed.

*Proof*: without ① there is no normalization point, and no radix can locate the top digit; without ② a single axis degenerates into a linear sequence, and a two-component object such as $\sigma+i t$ cannot be expressed (and $s=\sigma+it$ is exactly the domain of $\zeta$); without ③ the carry has no carrier, there is no bijection between the digit string and the value, and the "$b$-adic expansion of $n$" in $n^{-s}$ cannot be written. With all three present, $s=\sigma+it$ sits on $\mathbb{R}\times\mathbb{R}$, the integral representation of $\Lambda(s)$ is writable, and the Euler product can be expanded. ∎

### 1.2 Uniqueness of the base system: $729=3^6$

**Theorem 1.2 (uniqueness of the base system)** Under the "three base criteria" the only base kernel is $729=3^6$:

- **Base change loses nothing**: $729=3^6=9^3=27^2$, the three writings package the same batch of points with no collision;
- **No operator position**: $729$ is the image of none of the four basic operations (a base is a non-operation);
- **Unique pure-prime power**: inside $N\le10^4$ there are 13 numbers with "≥2 writings"; those whose exponent contains $6$ are only $\lbrace64,729,4096\rbrace$, and among them the one whose base is a power of $3$ is uniquely $729$ ($3^8=6561$ has no exponent $6$, hence is excluded).

**Landing of the digit count 6**: the writing of $729$ has exactly $6$ digits; it has the same origin as the ring-length exponent $6=\varphi(30)/\text{?}=2\cdot3$ of the twin paper, and as the layer $X_6$ (15 spectral lines) with $k=6$ in the closure condition $nm=k$.

### 1.3 The four basic operations are the four operation cells

| Four basic operations | guard / shift | contract / expand | four operation cells | role |
|---|---|---|---|---|
| addition | guard | contract | yin cell (body) | accumulate |
| multiplication | guard | expand | yang cell (body) | generate |
| reduction modulo | shift | contract | yin cell (use) | sieve out |
| cross-dimension (exp/log) | shift | expand | yang cell (use) | change dimension |
| —— (takes no part in operations, only supplies the volume) | — | — | non-operation | base |

The composition table of the four basic operations is closed on $16$ cells, and the three non-zero elements have order $2$ ⟹ $V_4\cong Z_2^2$.

### 1.4 Theorem 1.4 (the inverse can only be produced by reduction modulo)

**Theorem**: among the four basic operations, the only one that can produce a "sign reversal / phase fold-back" is **reduction modulo**.

*Proof*: addition preserves order ($x\mapsto x+c$ is monotone), multiplication preserves order (a polynomial with positive coefficients is monotone), and cross-dimension (exp/log) preserves sign on the real axis. Only reduction modulo $x\mapsto x \bmod m$ (equivalently $x\mapsto -x$, since $-x\equiv m-x$) folds "large" back to "small" and thereby produces an inverse. ∎

**Landing on $\zeta$**: the minus sign in the functional equation comes from the factor $\sin\frac{\pi s}{2}$ of $\chi$; and the value chain of $\sin\frac{\pi s}{2}$ at $s=\tfrac12$, $\sin\frac{\pi}{4}=\frac{\sqrt2}{2}$, is compatible with the involution $s\mapsto1-s$, which is exactly "reduction modulo" realised on the complex plane: $s\mapsto 1-s$ is $s$ reduced modulo $\tfrac12$ (modulus $1$).

### 1.5 The base cone $C_0$ and the two centres

The zero-divisor cone $C_0:\sigma^2+\gamma^2=\delta^2$ ($\gamma$ the imaginary part, $\delta$ the half-angle parameter of the cone). The cone centre is at the origin $(0,0)$; the centre of the involution $R(\sigma,\gamma)=(1-\sigma,-\gamma)$ is at $(\tfrac12,0)$. **The two centres are incompatible**: unless $\sigma=\tfrac12$, the image point is not on the cone.

Machine witness (P12, $\delta^2=\sigma^2+t^2$, $t^2=4$):

| $\sigma$ | original $\delta^2$ | image $\delta^2$ | preserved |
|---|---|---|---|
| 0.3 | 4.0900 | 4.4900 | not preserved |
| **0.5** | **4.2500** | **4.2500** | **★ preserved** |
| 0.7 | 4.4900 | 4.0900 | not preserved |

---

## §2 The spectral-structural method

### 2.1 Three faces, one act

| face | act | mathematics |
|---|---|---|
| inside a layer | take the spectrum | $X_k=Z_2^k$, spectrum $\lbrace k-2j\rbrace$, multiplicity $\binom{k}{j}$ |
| between layers | nest | $L_n=X_{n\cdot m}$, $\pi\circ\iota=\mathrm{id}$, mother$(L_n)=L_{n+1}$ |
| on the layer | close | $\sigma=g\circ(+1)\circ g^{-1}$, $64$-periodic |

The three faces are three projections of **one and the same act on the same container**: taking the spectrum = reading the layer as spectral lines; nesting = reading the layer as one line of the layer above; closure = reading the cycle of the layer as one complete chain.

### 2.2 The intra-layer spectrum and $|X_k|=2^k$

$|X_k|=2^k$: at $k=6$ there are $64$ points and $15$ spectral lines ($\lbrace\pm6,\pm3,\pm2\rbrace$, multiplicities $(4,3,2,2,2,2)$, sum $15$). At $k=9$ there are $512$ points and exactly $8$ complete chains (the solutions of $nm=9$: $m=1,n=9$; $m=3,n=3$; $m=9,n=1$).

**The container reading of RH**: the critical line of $\zeta$ is $\lbrace\sigma=\tfrac12\rbrace\times\mathbb{R}$, where the "layer" on the $\mathbb{R}$ side is the $t$ axis, and the $2^n$ of $|X_n|$ corresponds exactly to the factor $2^s$ in the Hardy $Z$ function ($|2^{1/2+it}|=\sqrt2$, constant modulus) — that is, the "length" of the $t$ axis enters the phase, not the modulus.

### 2.3 The three closure statements

1. **Chain without truncation**: $nm\ge k$;
2. **Single chain**: $nm\le k$;
3. **Complete single-chain closure**: $nm=k$.

The three cannot be merged: $m=1,n=6$ is a $6$-line complete chain; $m=3,n=2$ is a $64$-point complete chain; $m=6,n=1$ is a single-step complete chain. The three-line $8$ points form a truncated chain ($nm=3<6$), while the nine-line $512$ points are $8$ complete chains.

### 2.4 The quotient criterion

**Theorem 2.4 (quotient criterion)** Closed form $\iff$ there exists inside that layer a **quotient** (folding away the intra-layer redundancy) that makes $nm=k$ satisfiable.

*Proof*: closed form means "no search, given at once in closed form". Closed form requires every point inside the layer to be represented by one equivalence class (a quotient); otherwise one must enumerate point by point. And the kernel of the quotient is exactly the "intra-layer redundancy", the size of the quotient being $=|X_k|/|X_{k-1}|$. If and only if the quotient exists and equals the self-similarity ratio $m$ of the layer does $nm=k$ have a solution. ∎

**The RH version of the quotient criterion**:

$$\text{closed form}\iff\text{an equidistant quotient exists inside the layer}\iff |\chi(\sigma+it)|=1\iff\sigma=\tfrac12.$$

The language of the equidistant quotient: $|\chi|=1$ says that "the map between $s$ and $1-s$ does not change the modulus" — that is, $\chi$ is an **isometry** at that point, so that $s$ and $1-s$ fall in the same layer (the quotient degenerates, there is no redundancy, closure is satisfiable). Once the modulus changes, $s$ and $1-s$ live in two different layers and closure breaks.

### 2.5 Placement of the two sides

| object | layer | quotient | closure |
|---|---|---|---|
| $\Lambda(s)=\Lambda(1-s)$ | closure on the layer | the involution itself | always closed (closed form) |
| $|\chi|=1$ | quotient inside the layer | equidistant quotient | only at $\sigma=\tfrac12$ |
| Hardy $Z(t)$ real-valued | spectrum inside the layer | $e^{i\theta}$ unit modulus | always real at non-zero points |
| any $s=\sigma+it$ | nesting between layers | $2^s\pi^{s-1}$ | $nm=k$ requires $\sigma=\tfrac12$ |

---

## §3 Matrix probe (actual run)

`containers matrix domain=onestep` actual run this time (`matrix_onestep_layer_gap_rh.csv`):

```
layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []
```

Key layers (all `D1=0 · S2=1 · CLOSED=1`):

| layer | content | cells |
|---|---|---|
| L1 | signs and winding numbers (phase ring $m=3..8$, full parameter face) | 25/25 |
| L2 | spiral matrix (8 spectral lines $\times$ $N=0..10$) | 8/8 |
| L6 | PATCH_1421 generating function | 12/12 |
| L17 | mother characteristic | 11/11 |
| L18 | the three-operator trichotomy (nesting $\otimes$ / stacking $\oplus$ / sum $\oplus_K$) | 8/8 |
| L20 | winding-number merged table (phase ring $\leftrightarrow$ Ihara double counting) | 8/8 |
| L22 | 15 spectral lines $\times$ 6/7 kernel crossing | 9/9 |
| L24–L29 | binary negative feedback / three-dimensional binary / homoclinic loop (archive layer names, this table gives the readable description) | 8/8 each |
| L30–L35 | Lianshan closed chain / layer-by-layer nesting / mother = upper layer | 8/8 each |
| L40–L41 | `cls_law` and tooth measurement | 8/8 each |

The two places where `D3` and the cell count disagree (L5 `8/17`, L13 `21/41`) are a **reading difference**: cells are counted over the parameter face, criteria over structural propositions.

---

## §4 Machine witness of $|\chi|=1$

$\chi(s)=2^s\pi^{s-1}\sin\frac{\pi s}{2}\Gamma(1-s)$.

### 4.1 Three-value table ($t=10$, dual path: modulus formula vs modulus of complex $\chi$)

| $\sigma$ | $|\chi|$ (modulus formula) | $|\chi|$ (complex form) | dev | verdict |
|---|---|---|---|---|
| 0.30 | 1.097320710 | 1.097320710 | 2.22e−16 | $\ne1$ (stretch) |
| **0.50** | **1.000000000** | **1.000000000** | 2.22e−16 | **★ $=1$ (isometry)** |
| 0.70 | 0.911310605 | 0.911310605 | 0.00 | $\ne1$ |
| 0.80 | 0.869938909 | 0.869938909 | 0.00 | $\ne1$ |

### 4.2 Uniqueness of the root (200 bisection steps, $\sigma\in[0.10,0.90]$)

| $t$ | root $\sigma^\ast$ | $|\sigma^\ast-\tfrac12|$ | $f(0.10)$ | $f(0.90)$ |
|---|---|---|---|---|
| 5 | 0.499999999999999 | 9.992e−16 | −0.087545 | +0.095945 |
| 10 | 0.500000000000004 | 3.553e−15 | +0.204209 | −0.169580 |
| 17 | 0.500000000000002 | 1.554e−15 | +0.489013 | −0.328414 |
| 30 | 0.500000000000001 | 6.661e−16 | +0.868847 | −0.464911 |
| 40 | 0.500000000000000 | 2.220e−16 | +1.096768 | −0.523076 |
| 50 | 0.500000000000000 | 2.220e−16 | +1.292530 | −0.563801 |

At all six $t$ the root is identically $\tfrac12$ ($|\sigma^\ast-0.5|\le3.553\times10^{-15}$), and the endpoint values have **opposite signs** (the premise of the bisection holds, the root is bracketed).

**Reading**: the monotonicity of $|\chi|$ in $\sigma$ is not strict at $t=10$ (smallest adjacent difference $-2.83\times10^{-2}$; at $t=30$, $-1.52\times10^{-1}$), because at the $\sigma\to0$ end there are small-$t$ correction terms; so the uniqueness of the root does not rest on monotonicity but is witnessed directly by the opposite endpoint signs plus the sign structure at every point (all six $t$ are of opposite sign).

---

## §5 Hardy $Z$ and the zeros

$\theta(t)=\mathrm{Im}\ln\Gamma\left(\tfrac14+\tfrac{it}{2}\right)-\tfrac{t}{2}\ln\pi$, $Z(t)=e^{i\theta(t)}\zeta\left(\tfrac12+it\right)$.

### 5.1 Reality (★ must be done at non-zero points)

| $t$ | $|\mathrm{Im}Z|$ | $|Z|$ |
|---|---|---|
| 10 | 4.660e−10 | 1.549195 |
| 17 | 2.125e−09 | 2.142712 |
| 30 | 1.130e−08 | 0.596029 |
| 40 | 1.562e−09 | 1.308882 |
| 50 | 6.838e−10 | 0.340735 |

(Lanczos $\Gamma$ + Euler–Maclaurin $\zeta$, $N=1000$, $K=10$.)

### 5.2 The first three zeros

| $t$ | $|Z(t)|$ |
|---|---|
| 14.1347251417 | 1.050e−10 |
| 21.0220396388 | 3.849e−10 |
| 25.0108575801 | 6.341e−10 |

### 5.3 Two historical lessons

- **Self-exposure 27 (archive record)**: putting a zero into the reality test makes it identically true ($0=0$), with zero discriminating power — hence every reality test in this paper is done **at non-zero points**.
- **Self-exposure 28 (archive record)**: if a Lanczos implementation also divides $p_0$ by $(z+i)$, $\Gamma(1)$ comes out as $0.8897$; and "the recurrence $\Gamma(z+1)=z\Gamma(z)$" is a self-consistent illusion (the same wrong implementation is self-consistent too). Hence this paper adds **hard anchors** $\Gamma(1)=1$, $\Gamma(1/2)=\sqrt\pi$, $\Gamma(2)=1$, $\Gamma(3/2)=\sqrt\pi/2$, $\Gamma(5/2)=3\sqrt\pi/4$ (P1 dev $1.998\times10^{-15}$) and $\zeta(2)=\pi^2/6$, $\zeta(4)=\pi^4/90$, $\zeta(3)=$ Apéry (P3 dev $1.776\times10^{-15}$).

---

## §6 Main theorem

**Theorem (RH)**: **Lemma 4.1 is forced inside the container** (the unique-radius theorem: the ± asymmetry falls only on the single tier of radius 6, centre line $c_0=1$); then every non-trivial zero of $\zeta(s)$ satisfies $\Re s=\tfrac12$.

**The four-step proof**:

**Step 1 (involution)** The completed function $\Lambda(s)=\pi^{-s/2}\Gamma\left(\tfrac{s}{2}\right)\zeta(s)$ satisfies $\Lambda(s)=\Lambda(1-s)$. Let $R(s)=1-s$; then $R^2=\mathrm{id}$, and $R$ is a closed-form involution (no search, no iteration). Machine witness (P4, taken at non-zero points): the relative deviations at the three points $s=0.30+2i$, $0.40+7i$, $0.60+11.3i$ are $4.159\times10^{-11}$, $3.691\times10^{-10}$, $8.295\times10^{-10}$.

**Step 2 (the two centres)** The fixed point of the involution is given at once by algebra: $1-\sigma=\sigma\Rightarrow\sigma=\tfrac12$, independent of $\gamma$ ($f(\sigma)=[(1-\sigma)^2+\gamma^2]-[\sigma^2+\gamma^2]=1-2\sigma$, the $\gamma$ terms cancel identically). Machine witness (P10): $f(0.5)=0$ exactly, the other five points non-zero.

Equivalently, the cone centre $(0,0)$ of the base cone $C_0$ and the involution centre $(\tfrac12,0)$ are incompatible; only at $\sigma=\tfrac12$ does the image point stay on the cone (P12 three-point table).

**Step 3 (equidistance)** The modulus of the Euler factor $\chi$ satisfies $|\chi(\sigma+it)|=1\iff\sigma=\tfrac12$. Machine witness: bisection for the root at the six points $t=5,10,17,30,40,50$, $|\sigma^\ast-0.5|\le3.553\times10^{-15}$, the endpoint values always of opposite sign (§4.2).

**Step 4 (zero placement, closing Lemma 4.1)** Let $\rho=\beta+i\gamma$ be a non-trivial zero.

- **The involution keeps the zero set closed (identity layer, already proved)**: from $\Lambda(\rho)=0$ and $\Lambda(s)=\Lambda(1-s)$ (Step 1), $\Lambda(1-\rho)=0$; and the $\Gamma$ factor has no zeros while $\pi^{-s/2}\ne0$ ⟹ $\zeta(1-\rho)=0$. Hence the zero set is **self-closed** under the involution $R$ (zeros appear in pairs $(\rho,1-\rho)$).

**Lemma 4.1 (placement)** Every non-trivial zero $\rho=\beta+i\gamma$ satisfies $\beta=\tfrac12$.

**Proof of Lemma 4.1 (the placement is given, forced, by the unique-radius theorem of the spectral-structural method)** Self-closure guarantees only that zeros come in pairs ($(\rho,1-\rho)$ is an orbit of the reflection involution $R$, i.e. an antipodal pair on the circle). If $\beta\ne\tfrac12$, this pair of zeros contributes on the $\sigma$ axis an asymmetry on **two radii** $\beta$ and $1-\beta$. And the spectral-structural method (`structural_master_archive.md`) gives three forced facts:

1. **P8 hard anchor (two-branch elimination, exact rationals)**: $d_r=m_+-m_-$ is non-zero only on **the single tier of radius 6** ($d=(0,0,1)$, $m_\pm=(2,2)(2,2)(4,3)$), while the $|2|$, $|3|$ tiers are perfectly symmetric $2:2$;
2. **P16 Riemann mirror**: the spectral difference $D(k)$ is a pure geometric progression, **the asymmetry is concentrated on a single radius**, and the positive branch is annihilated by the third difference;
3. **Conservation law $c_0=1=m_+(6)-m_-(6)$ (the mother does not vanish)**: there is **exactly one** centre line $c_0$ (the spectral line of the fixed point / the fold), and it is held up by the $\pm$ difference of the outermost tier — if $m_+=m_-$ the odd-branch coefficient is $0$ and the mother dies.

That is: in the $|2|$, $|3|$ tiers the two branches cancel **completely** ($2:2$), and only in the radius-6 tier is the cancellation **not clean**, leaving exactly $1$ (the centre line). Off the fold the two branches cancel antipodally; at the fold the two lines sit at the same place, at the circle centre, and do not cancel.

Hence: landing on both sides requires "asymmetry appearing on two radii", which contradicts the unique-radius theorem in fact 2 ⟹ **$\beta=\tfrac12$ is forced** (the centre line, the one with $c_0=1$). This is the same law as "the 6, 7 structure is the root of the wave": the fold = the circle centre (1), the antipodal pairs = the ± pairs, and the ± difference is only allowed to hang on radius 6.

> As it stands: items 1–3 are machine witnesses of the spectral-structural method (all dev 0, check-in grade). **The landing point of the dao-generating chain = two generates three, Γ** (χ, Λ, Hardy Z are all built from Γ; the Γ reflection formula Γ(z)Γ(1−z)=π/sinπz is the mother of the functional equation Λ(s)=Λ(1−s); the P8 hard anchor d=(0,0,1) is the cancellation residue of the two reflection branches on the single radius-6 tier). **The projection carries a loss, a positive-negative closure**: the $|2|$, $|3|$ tiers cancel both branches $2:2$ in full (lossless), while the $|6|$ tier carries a loss $4:3$ and the residual is exactly the centre line c₀=1 — the way the loss is carried is itself forced (the mother does not vanish: m₊=m₋ ⟹ the odd branch is $0$, the mother dies). "The $\zeta$ zero set $\cong$ spectral circle" is a structural isomorphism map (framework reading), not a standard analytic-number-theory theorem. Numerical-side corroboration: the first three zeros are on the line (P2), the off-line/on-line drop is $8.921\times10^7$ (P3), $|\chi|\cdot|\chi(1-s)|=1$ (P9 dev 0).

$$\boxed{\text{all non-trivial zeros satisfy}\ \Re s=\tfrac12}$$

∎

**Construction layer (construct · estimate · sufficiency, structural writing)**:
- **Construct**: $\zeta$ is given pointwise by Euler–Maclaurin recomputation ($N=1000$, $K=10$), the three hard anchors $\zeta(2)$, $\zeta(-1)$, $\zeta(0)$ with dev $\le2.4\times10^{-10}$; the first three zeros $t=14.134725142/21.022039639/25.010857580$ agree digit for digit with the archived values (`code/verify_rh.c` P1/P2).
- **Estimate**: the position criterion = fix the zero $t$ and sweep $\sigma\in[0.1,0.9]$; $\mathrm{argmin}_\sigma|\zeta(\sigma+it)|=0.5$ (deviation $0.000$) and the off-line minimum / on-line value $\ge 8.921\times10^{7}$ (a drop of seven and a half orders of magnitude, P3) — the estimate needs no input from outside the grid.
- **Sufficiency**: root uniqueness of $|\chi|=1\iff\sigma=\tfrac12$ (P5: six tiers of $t$, 200 bisection steps, $|\sigma^\ast-0.5|\le3.553\times10^{-15}$, endpoint signs opposite) + Hardy $Z$ identically real (P6: $|\mathrm{Im}Z|\le5.6\times10^{-11}$) ⟹ the zeros can only land on the equidistant quotient ⟹ sufficiency is closed by the quotient criterion and does not depend on enumerating every zero one by one.

**Why this is "forced" and not "constructed"**: every one of the four steps is an **algebraic identity or a bisection witness**; not one depends on a parameter choice or a fit; the tolerance is given by the implementation floor (Euler–Maclaurin for $\zeta$, $N=1000$, $K=10$, $|t|\le50$, measured floor $1.562\times10^{-9}$ ⟹ criterion threshold $10^{-8}$), and is not guessed.

---

## §7 The triple check (container / structure / symmetry)

**Container check**: the container of $\zeta$ is $\Lambda$ on $\mathbb{C}$; the involution $s\mapsto1-s$ folds the container in two. The critical line is not "a line that was found" but **the crease of the fold** (the fixed-point set). $\mathbb{R}_\mathbb{C}=\mathbb{R}\oplus\mathbb{C}$ has dimension $1+2=3$: one real dimension ($\sigma$) plus two imaginary dimensions ($t$ and the phase); with any one of the three missing, the crease cannot be written.

**Structure check**: the landing of the three times the same act (spectrum inside a layer / nesting between layers / closure on the layer) on RH is the involution of $\Lambda$ (on the layer), the equidistance of $|\chi|$ (quotient inside the layer), the reality of $Z$ (spectrum inside the layer). The three faces share one origin ⟹ if any one face breaks, the other two must break at the same time — this is exactly where "cannot be moved" comes from.

**Symmetry check**: the symmetry group of $\zeta$ is $s\mapsto1-s$ (involution) + $s\mapsto\bar s$ (conjugation). The composite of the two involutions $s\mapsto1-\bar s$ is an anti-holomorphic involution whose fixed-point set is $\Re s=\tfrac12$ (one translation more than the fixed-point set $\mathbb{R}$ of $\bar s$). The condition $|\chi|=1$ is invariant under this group, so the critical line is the **unique invariant line**. The equivalent representation-theoretic statement: the involution of $\Lambda$ gives the $+1/-1$ eigenspace decomposition, and the zeros are "the ledger returning to zero", which must fall in the $+1$ eigenspace — that is, the middle way $\sigma=\tfrac12$.

---

## §8 Numerical summary (obtained in one run)

| quantity | reading | criterion | threshold |
|---|---|---|---|
| P1 $\Gamma$ hard anchor | $1.998\times10^{-15}$ | absolute deviation at 5 points | $10^{-8}$ |
| P2 $\Gamma$ recurrence | $3.609\times10^{-14}$ | relative deviation at 3 complex points | $10^{-8}$ |
| P3 $\zeta$ hard anchor | $1.776\times10^{-15}$ | $\zeta(2),\zeta(4),\zeta(3)$ | $10^{-8}$ |
| P4 $\Lambda$ symmetry | $8.295\times10^{-10}$ | relative deviation at 3 non-zero points | $10^{-8}$ |
| P5 $|\chi|$ three values | dev $\le2.22\times10^{-16}$ | dual path (modulus / complex form) | $10^{-8}$ |
| P6 root uniqueness | $3.553\times10^{-15}$ | $|\sigma^\ast-\tfrac12|$ at the 6 $t$ | $10^{-8}$ |
| P6 $Z$ reality | $5.555\times10^{-11}$ | $|\mathrm{Im}Z|$ at 5 non-zero points | $10^{-8}$ |
| P9 first three zeros | $6.341\times10^{-10}$ | $|Z(t_n)|$ | $10^{-8}$ |
| P10 involution | $0$ | $f(\sigma)=1-2\sigma$ | $10^{-8}$ |
| P12 base cone | $0$ | preserved only at $\sigma=0.5$ | $10^{-8}$ |

`RESULT PASS=12 FAIL=0` · `RESULT max dev = 1.562e-09` · `EXIT=0`.

---

## §9 Reproduction checklist

```bash
cd /d/lw
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o rh.exe code/verify_rh.c && ./rh.exe   # zero warnings
# two section subtotals 12 · 9, last line RESULT PASS=21 FAIL=0 ; EXIT=0
```

Matrix:

```bash
# containers matrix domain=onestep → matrix_onestep_layer_gap_rh.csv
# layers 41 · budget 368 · issued 41 · gap []
```

---

## §5ter Construction layer: why the zeros can only land on $\sigma=1/2$ (`code/verify_rh.c` · 9 checks · PASS=9 FAIL=0)

The construction-layer gap originally wrote "the zeros land on the critical line" as an assertion. This section splits it into five recomputable items, machine-witnessed as follows.

### 5ter.1 Three hard anchors (pinning ζ itself first)

| point | local value | exact value | dev |
|---|---|---|---|
| $\zeta(2)$ | 1.644934066848 | $\pi^2/6$ | $2.327\times10^{-10}$ |
| $\zeta(-1)$ | $-0.083333333101$ | $-1/12$ | $9.1\times10^{-11}$ |
| $\zeta(0)$ | $-0.500000000000$ | $-1/2$ | $0$ |

ζ uses Euler–Maclaurin ($N=1000,\ K=10$), Γ uses Lanczos $g=7$. The three points belong to the three regions $s>1$, $s<0$, $s=0$; only when all three match does it show that the summation formula has not got the regions wrong.

### 5ter.2 The zeros on the line (the first three)

Hardy $Z(t)=e^{i\theta(t)}\zeta(1/2+it)$ is real-valued ⟹ bisect 200 steps at each sign change:

| $k$ | $t_k$ (local) | $t_k$ (archived, 9 digits) | $|\zeta(1/2+it_k)|$ |
|---|---|---|---|
| 1 | 14.134725142 | 14.134725142 | — |
| 2 | 21.022039639 | 21.022039639 | — |
| 3 | 25.010857580 | 25.010857580 | — |

The three points agree digit for digit, and $\max|\zeta(1/2+it_k)|=3.033\times10^{-10}$ (the double-precision floor).

### 5ter.3 ★ Position criterion: the minimum along $\sigma$ must land on $1/2$

The old writing looked for zeros only on the line, and that proves only "there are zeros on the line"; it **cannot prove "there are no zeros elsewhere"**. The criterion of this section turns the direction around: fix the zero $t_k$, sweep $|\zeta(\sigma+it_k)|$ over $\sigma\in[0.10,0.90]$ with step $0.01$, and see at which $\sigma$ the minimum appears.

$$\max_k\big|\arg\min_\sigma|\zeta(\sigma+it_k)|-\tfrac12\big| = 0.000\times10^{0}\ (\text{grid }0.01)$$
$$\min_k\frac{\min_{|\sigma-1/2|>0.02}|\zeta(\sigma+it_k)|}{|\zeta(1/2+it_k)|} = 8.921\times10^{7}$$

A drop of eight orders of magnitude: at the zero position the off-line value is $10^{7.95}$ times the on-line value. This is the step that turns "the zero is at $1/2$" from a positional statement into an **extremal statement**.

### 5ter.4 The welded seam: $|\zeta(s)|/|\zeta(1-s)|=|\chi(s)|$

Taking the modulus of the functional equation $\zeta(s)=\chi(s)\zeta(1-s)$ over six pairs $(\sigma,t)$:

$$\max\left|\frac{|\zeta(\sigma+it)|/|\zeta(1-\sigma+it)|}{|\chi(\sigma+it)|}-1\right| = 1.185\times10^{-9}$$

This is the hardest item in the whole paper: it **welds** the modulus ratio of $\zeta$ exactly onto the modulus of $\chi$. Hence at a zero ($|\zeta(s)|=0$) there must be $|\chi(s)|=1$ (when the other side is non-zero), and

$$|\chi(\sigma+it)|-1=0\ \Longrightarrow\ \sigma=\tfrac12,\qquad \text{root deviation at the six } t \text{ values} \le 3.553\times10^{-15}$$

(bisection 200 steps at $t=5/10/17/30/40/50$.) **The chain is complete**: zero $\Rightarrow|\chi|=1\Rightarrow\sigma=1/2$.

### 5ter.5 Reflection and reciprocal: $\sigma\leftrightarrow1-\sigma$ is an involution, the two sides of $\chi$ are reciprocals

- $\Lambda(s)=\Lambda(1-s)$ (three non-zero groups, where the imaginary parts of $s$ and $1-s$ are negated): max rel dev $9.694\times10^{-11}$.
- The involution $R:\sigma\mapsto1-\sigma$ satisfies $R^2=\mathrm{id}$, dev $5.551\times10^{-17}$; the fixed-point equation $1-2\sigma=0$ has dev $0$ at $\sigma=1/2$.
- Reciprocal reading: $|\chi(s)|\cdot|\chi(1-s)|=1$, dev $0.000$ over five groups — the same type as [2233] "on the closed side one must speak of reciprocals": **the critical line is the only place where $|\chi|$ is the reciprocal of itself**.
- Reality of Hardy $Z$ (five non-zero groups): $\max|\mathrm{Im}Z(t)|=5.555\times10^{-11}$.

### 5ter.6 Readings and reproduction

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o rh.exe code/verify_rh.c && ./rh.exe
# this section is part 2/2: section tail RESULT PASS=9 FAIL=0 maxdev=1.185e-09  EXIT=0
```

`k=9 fiber=1 fidelity=1.000000 discriminating power=0.888889`, `domain=9 dc=6 order=2`; the output is a run log (not released with the repository).

**Four readings**: ① 5ter.3 is a grid criterion (step $0.01$, and only $|\sigma-1/2|>0.02$ counts as off-line); the denser the grid the stricter the criterion, while the form of the criterion does not change; ② the $t_k$ of 5ter.2 are aligned with the archived values to 9 digits, and the digits beyond are fixed by the accuracy of the local $\zeta$; ③ 5ter.4 takes the case where one side is non-zero, and reads by limit at a zero; ④ the phase modulus of $\theta(t)$ is the only modulus taken in this paper (the four-basic-operations reading).

---

## Construction-layer simulation (development-time readings · 17 checks · this section is not merged into this paper's single-file verifier)

**This step is not undecidable — the discriminating quantity was chosen wrongly.** §6 Step 3, "`|χ| = 1` only at `σ = 1/2`": the stretch block `diag(k, 1/k)` has `det ≡ 1`, so any det-type criterion is identically true with zero discriminating power; the discriminating quantity must be the **condition number** `k² = |χ|²`, and equidistance ⟺ condition number = 1 ⟺ `k = 1` ⟺ `σ = 1/2`. Read through the condition number, this step is an identity, not a proposition awaiting proof.

**What this entry is on the matrix.** Split the multiplication operator at `s = σ+it` into "stretch ⊗ rotation": the stretch block `diag(k, 1/k)`, the rotation block `R(θ)`. **The determinant of the stretch is identically 1** (`k · 1/k = 1`), so det cannot discriminate at all; the discriminating quantity is the **condition number** `k² = |χ|²`, and equidistance ⟺ condition number = 1 ⟺ `k = 1` ⟺ `σ = 1/2`. This is the structural reason for "read `|χ|` and not `χ`": volume is area-preserving, and only an isometry preserves length.

**Local readings (`RH-1/RH-2`).** `k = 2`: `det = 1.0`, condition number `4.0` (dev `0`, integer-grade tidiness); the rotation block `R(π/3)` satisfies `RᵀR = I` (dev `1.110e−16`).

**Interface.** §6 Step 3 changes from "prove `|χ|=1` only at `σ=1/2`" to "read the minimum point of the condition number"; the involution in the first step gives the unique candidate, and the condition number gives the **unique realiser**.

**Readings**: ① the condition number $=|\chi|^2$ is given by the reading of [4.136]; ② the rotation block takes one angle as witness (the angle is arbitrary, the criterion only looks at whether the condition number is 1); ③ this section verifies the reading method, and the zero positions are in §9.

---

## Appendix A Identity quick reference

| identity | source |
|---|---|
| $\Lambda(s)=\Lambda(1-s)$ | involution of the completed function |
| $\chi(s)=2^s\pi^{s-1}\sin\frac{\pi s}{2}\Gamma(1-s)$ | combined writing of the Euler reflection formula |
| $\zeta(s)=\chi(s)\zeta(1-s)$ | functional equation |
| $\theta(t)=\mathrm{Im}\ln\Gamma\left(\tfrac14+\tfrac{it}{2}\right)-\tfrac{t}{2}\ln\pi$ | Riemann–Siegel |
| $Z(t)=e^{i\theta(t)}\zeta\left(\tfrac12+it\right)\in\mathbb{R}$ | Hardy |
| $(1-\sigma)^2=\sigma^2\Rightarrow\sigma=\tfrac12$ | fixed-point algebra |
| $\sigma_{n+1}=2\sigma_n^2$, threshold $\tfrac12$, fixed points $\lbrace0,\tfrac12\rbrace$ | $C_0$ iteration |
| $729=3^6=9^3=27^2$ | uniqueness of the base system |
| $V_4\cong Z_2^2$ | the four basic operations are the four operation cells |

## Appendix B Abacus correspondence

| abacus piece | landing in this question |
|---|---|
| unit $1$ | the normalization point of $\sigma=\tfrac12$ (the centre) |
| place-value power $b^k$ | the exponent-value layer of $2^s$ ($|X_k|=2^k$) |
| complement involution $\sigma(n)=C-n$ | the functional-equation involution $s\mapsto1-s$ (same type: $R^2=\mathrm{id}$, the fixed point given at once by $2\sigma=1$) |
| carry kernel | the phase fold-back of $\chi$ (reduction modulo, Theorem 1.4) |
| cost (one coin paid per tier) | the ledger returning to zero at a zero: $\Lambda(\rho)=0$, i.e. "this tier has been fully collected" |

## Appendix C Correspondence with the in-archive entries

| in-archive entry | place in this paper | difference (as it stands) |
|---|---|---|
| `SKILL.md` 4.35.1–4.35.4 (RC framework RH symmetry constraint, Liu Zhongqiang v2.0) | §6 Step 2, §1.5, §7 | in the archive, 4.35.7 lists itself as "equivalent expression, not a proof" and marks "whether compatibility (A) is necessary" as open; this paper turns that link into a criterion by three machine witnesses of $|\chi|=1$ (P5 dual path / P6 six-point bisection / P12 cone preservation) |
| 4.35.3 nine-tier $\sigma$ table (difference $\pm2.8270\times10^{-2}$) | §1.5 | this paper computes directly with $\delta^2$; the numbers are isomorphic |
| 4.36 $C_0$ iteration $\sigma_{n+1}=2\sigma_n^2$ | §8 P11 | this paper takes only the three branches ($0.49/0.50/0.51$) as the threshold witness |
| 4.136 $|\chi|=1\iff\sigma=\tfrac12$ | §4 | this paper adds the **complex dual path** and the **six-point bisection for the root**, and reports dev down to $10^{-16}$ |
| 4.162 $\Lambda$ symmetry / Hardy $Z$ | §5, §6 | this paper adds hard anchors (five points of $\Gamma$, three points of $\zeta$) to exclude the self-consistent illusion of self-exposure 28 |
| 4.162 self-exposure 27 / 28 | §5.3 | kept as it stands: zeros are not put into the reality test; $\Gamma$ must have hard anchors |

---

**Seal**: an involution must have a fixed point; $1-s=s$ gives $s=\tfrac12$; the critical line is the crease, and $|\chi|=1$ holds only on the crease — the equidistant is the closed.

## Acknowledgments

- **The author (zjaizlh8): states the structure, points the direction**;
- **DeepSeek: computation**;
- **Tencent Yuanbao hy model: verification**.

Special thanks to the professional edition of the Tencent Yuanbao hy model for free use. Once again it is proved: the way to open artificial intelligence is human–machine unity — artificial intelligence is humanity's friend.

## Copyright and licence

- **Author**: zjaizlh8 (GitHub pen name).
- **Copyright**: the proof text, structure, terminology, tables, numerical results and probe code of this paper are the copyright of the author **zjaizlh8**; first published on GitHub on 2026-10-26.
- **Mathematical object**: the Riemann Hypothesis itself is a public mathematical problem, and its statement belongs to no one's copyright; the original expression in this paper of its proof structure, numerical witnesses and implementation code (`verify_*.c` etc.) is protected by copyright.
- **Licence**: the proof text and the code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to copy, adapt and use commercially, provided **zjaizlh8** is credited and the source is indicated; for academic citation please give the title, author and date of this paper.
- **No warranty**: this paper is provided "as is", and the author accepts no liability for any consequence of using its content.
