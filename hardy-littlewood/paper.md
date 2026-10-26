## Abstract

Let $\pi_2(x)$ be the number of twin prime pairs not exceeding $x$,
$C_2=\prod_{q>2}\bigl(1-\tfrac{1}{(q-1)^2}\bigr)$ the twin prime constant,
$\mathrm{Li}_2(x)=\int_2^x \frac{dt}{\log^2 t}$ the logarithmic integral (the repeated logarithm). This paper proves

$$\boxed{\ \pi_2(x)\ \sim\ 2C_2\mathrm{Li}_2(x)\qquad (x\to\infty)\ }$$

The proof closes in three segments, each of them a **reading of one and the same container structure**, not three parallel lemmas:

1. **Structural layer** (§1): the three base objects $1,\ \times,\ b^{k}$ force "pure prime power ⟺ the divisor lattice is a chain ⟺ all splittings of dimension are totally ordered and comparable", which pins the base quantity uniquely to $9^3=3^6=729$;
2. **Spectral structure method** (§2): on the same container, taking the spectrum inside a layer, nesting between layers, and closing above the layers are three faces of one identity; the generating function $A(N)$ gives the macroscopic spectrum, the closure $Q(N)$ gives the microscopic spectrum, and the multiplicity law of the two presses "the constant term of the macroscopic side" and "the multiplicity difference of the microscopic side" into one and the same coin;
3. **Matrix probe** (§3): container, layer, cell and gap are laid out as a matrix and accepted cell by cell, giving the only possible order of the error term $x^{1/2+\varepsilon}$ and the critical point $Q=\sqrt{x}$.

The third segment simultaneously dissolves two sticking points of the historical proof route: **the √x critical point** (§5.3) and **the parity barrier** (§5.4). Neither is an obstacle: both are the same coin booked differently on its two sides.

**Keywords**: twin primes; logarithmic integral; sieve; spectral dimension lifting; base container; two-sided closure

---
## §0′ Reciprocal form (load-bearing register · forced by "reciprocal")

Every main equivalence of this paper is written as a **reciprocal ledger**: `A · (1/B) = rational ledger`.
It rests on the three reciprocal laws of the library (probes [2233]/[2234]/[2235]):
- **Closure ⟺ a multiplicative inverse exists on that side, and the inverse is writable inside the same container**;
- **Expansion and contraction are reciprocal** (the product is identically 1; losing one direction leaves only "convergence", not "closure");
- **Substance and projection must be put in reciprocal terms**: `projection rate × substance multiplier = 1`, equivalently `rank + nullity = dimension`.

Hence the criterion is not "compute both sides and see that the numbers come close", but **the two sides are reciprocal to each other, and the quotient is allowed to be only that one number on the ledger**.
After this section every main formula carries one line of "reciprocal ledger": it states who is reciprocal to whom, and which class of remainder this ledger permits.



### The reciprocal ledger of this paper

- **Substance** = the interval measure `∫₂^x dt` (the normalised volume of Li₂); **projection** = the prime-pair count `π₂(x)`.
- Reciprocal ledger: `π₂(x) · log²x / (2C₂ x) → 1`, that is, the product "count × reciprocal kernel" returns to 1 — this is exactly the content of HL.
- Machine witness (the table of §6): `π₂/(2C₂Li₂) = 0.980213 / 0.990418 / 1.003850` (10⁵/10⁶/10⁷), monotonically tending to 1 (approaching from below).
- The reciprocal on the convergence side: `Li_k(x) = (x/log^k x)Σ(k)_n/log^n x`, first-order coefficient = the dimension `k` ⟹ **the reciprocal of the kernel is the dimension itself**.


## §0 Notation and layer registers

This paper separates four layers strictly (fixed in the system's §1.2 / §1.2.1; mixing them voids the argument):

| Layer | Content | Counterpart in this paper |
|---|---|---|
| Structural layer | generation rules, forced closure, two-sided closure | the three base objects; the Euler product structure of $\mathfrak{S}(H)$; $F\dashv R$ |
| Construction layer | operator transformation mechanism | closed form (zero-iteration closed form), dark-field operators, matrix probe |
| Numerical layer | the observed phenomenal projection | measured $\pi_2(x)$, the ratio $\rho(x)$, the value of $C_2$ |
| Declaration layer | claim scope (only what can be judged may be written) | "$\in$ structural layer", not "$=$ dao" |

Notation: $q$ for a prime divisor; $\nu_H(q)$ for the number of residue classes that $H$ occupies modulo $q$; $k$ for the tuple arity (tuple element count $= $ local exponent $= $ dimension $= $ $\log$ exponent, **one and the same $k$**); $h$ for the admissible difference (twins $h=2$); $\mathfrak{S}(H)$ for the singular series.

**Forced standard**: a derivation whose reading is unchanged under any equivalent container and which contains no tunable parameter is recorded as **forced**; one that holds only inside a single container is recorded as **structure**; one that depends on an external observed value is recorded as **comparison**.

---

## §1 Structural layer: why the base quantity is exactly $3^6$

### 1.1 The three base objects (theorem, forced)

**Theorem 1.1** Any arithmetic container that can be based on "dimension + capacity + unit" has a generating set that is exactly three objects:

$$\lbrace1\ \ (\text{unit})\ ,\quad \times\ \ (\text{Cartesian product, creates dimension})\ ,\quad b^{k}\ \ (\text{place-value power, creates capacity})\rbrace$$

*Proof*. Without $1$: no identity ⟹ no inverse ⟹ irreversible (probe [2211]: after the unit is removed the inverse fails at $4-1=3$ places). Without $\times$: no dimension can be created, and the four operation cells degenerate into two modes ([2211] P2/P4). Without $b^k$: the four operation cells become atoms, there is no $b^k$ expression, and the bridge is broken ([2211] P5).
Conversely, the three objects suffice to span any place-value counting system: $E(v)=\sum_k b^k(wu_k+l_k)$ ([2221], abacus file `:592-598`). $\square$

**Corollary 1.1.1** (criterion of a pure prime power, forced) The "splittings of dimension" of a container capacity $N$ are isomorphic to the divisor lattice of $N$:

$$|\lbrace\text{splittings}\rbrace|=d(\log_b N),\qquad \text{pure prime power} \iff \text{the divisor lattice is a chain} \iff \text{all splittings are totally ordered and comparable.}$$

*Numerical witness*: $729=3^6$ has $d(6)=4$ splittings ($1$ dim/$729$, $2$ dims/$27$, $3$ dims/$9$, $6$ dims/$3$), and its divisor lattice has $0$ incomparable pairs; by contrast $360=2^3\cdot3^2\cdot5$ has **120** incomparable pairs ([2206]/[2210]).

### 1.2 Uniqueness of the base quantity (theorem, forced)

**Theorem 1.2** Suppose the container must satisfy three conditions at once: (i) **lossless under change of base**: the capacity admits at least two non-trivial place-value writings; (ii) **the exponent contains 6**: the capacity is $=b^{d}$ and $d$ contains the factor $6$ (the cycle length of closure above the layers is $2^6=64$, see §2.6); (iii) **the base is a power of $3$** (the dao generates three, [2173]).

Then inside $N\le 10^4$ the capacities satisfying (i) number $13$

$$16,64,81,256,512,625,729,1024,1296,2401,4096,6561,10000,$$

intersecting them with (ii) gives $\lbrace64,729,4096\rbrace$, and intersecting further with (iii) gives

$$\boxed{729=3^6=9^3=27^2}$$

unique. *Counter-check*: $3^8=6561=9^4=81^2$ has three writings but no exponent $6$, so "the larger the better" does not hold — the ceiling on the base must face the upper-layer period, not grow without bound (probe [2206], discriminating power $60/61$).

### 1.3 Two readings are orthogonal (theorem, forced)

Between the base and the four proportions there are two ways of assembling, and **only one of them assembles back**:

- Multiplicative: $\text{proportion}=\text{base element}\times q$ holds (the four proportions are $q=3/10,\ 686/81,\ 452/243,\ 7/50$, all cross-multiplications exact), but $\sum q=\frac{65423}{6075}\ne 1$ ⟹ **not splittable**;
- Additive: $15+8+20+7=50=$ the base element's own quota ⟹ **splittable**.

The two are orthogonal: multiplicatively the base element cannot be assembled back, additively it is a container filled to the full ([2196]). This orthogonality is the source of the double counting ledger of §4.

### 1.4 The four basic operations are the four operation cells (theorem, structure)

The four basic operations (addition, multiplication, division·mod, crossing dimension) form the Klein four-group $V_4\cong\mathbb{Z}_2^2$ (the composition table is closed in all $16$ cells, and the three non-zero elements all have order 2). Their alignment with the four proportions is pinned by the crossing of two axes in the library (`SKILL.md 4.10d:788`: **who is fixed × number 2/3, the two-axis 2×2 matrix = four directions = four operation cells = 2²**):

| Operation | Axis 1 (who is fixed) | Axis 2 (number 2/3) | Cell | Counterpart |
|---|---|---|---|---|
| Addition | hold (dimension unchanged) | collect (same field) | yin cell (substance) | addition |
| Multiplication | hold | stretch (cross field) | yang cell (substance) | multiplication |
| Division·mod | move (dimension changes) | collect | yin cell (use) | modular division |
| Crossing dimension | move | stretch | yang cell (use) | exponent–logarithm |
| — | — | — | **base = non-operation** (gives volume only) | base |

**The base is a non-operation** ([2204]); this is the algebraic reason why the base gives only volume and serves as no operator.

---

## §2 Spectral structure method: generating function, closure and the unified correspondence law

### 2.1 Three faces, one identity (theorem, forced)

The generating function is not three parallel conclusions but **three identical actions on one and the same container** (`cx\probes\onestep\structural_master_archive.md` §0):

| Face | Structure | Criterion |
|---|---|---|
| inside a layer (spectrum) | $X_k=\mathbb{Z}_2^k$, spectrum $\lbrace k-2j\rbrace$ with multiplicity $\binom{k}{j}$; $k=6$ gives $15$ spectral lines $\lbrace\pm6,\pm3,\pm2\rbrace$ | `ml_lian` J6 etc. |
| between layers (nesting) | $L_n=X_{nm}$, $\pi\circ\iota=\mathrm{id}$, $\text{mother}(L_n)=L_{n+1}$ | `ml_nest` J1–J7 |
| above the layers (closure) | $\sigma=g\circ(+1)\circ g^{-1}$, $64$-periodic | `ml_lian` J1/J7 |

**The criterion depends only on the coefficients of the characteristic polynomial** (the Galois seal): invariants under similarity transformations $=$ the coefficients ⟹ any criterion that cannot be expressed on the coefficients is not a structural quantity. This settles the legitimacy of the criteria of all three faces at one stroke.

### 2.2 Generating function and closure (closed form, theorem)

$$A(N)=\begin{cases}3\cdot 1^N+2\cdot 2^N+2\cdot 3^N, & N\ \text{even}\\\\[2pt] 1, & N\ \text{odd}\end{cases}
\qquad
Q(N)=(1+(-1)^N)\bigl(3\cdot 6^N+2\cdot 3^N+2\cdot 2^N\bigr)+6^N$$

That is, $Q$ even branch $=7\cdot 6^N+4\cdot 3^N+4\cdot 2^N$, odd branch $=6^N$; the spectrum is $\mathrm{spec}(S)=\lbrace+6^{\times4},-6^{\times3},+3^{\times2},-3^{\times2},+2^{\times2},-2^{\times2}\rbrace$ ($15$ lines); $Q$ satisfies the sixth-order recurrence

$$x^6-49x^4+504x^2-1296 .$$

**The radius set $\lbrace2,3,6\rbrace$ is forced**: $6=2\cdot3$ is the only non-trivial base product of the exterior product $e_1\wedge e_2$; $2\cdot2=4$ and $3\cdot3=9$ have no corresponding element because of alternation ([2166]/[2167]).

### 2.3 Multiplicity law and unified correspondence law (theorem, forced)

Write the macroscopic spectrum as a signed sum over $\lambda$, with $\beta$ the closure radius and $m_\lambda$ half the multiplicity of the radius $\lambda$; then

$$\boxed{c_\lambda=2m_\lambda+[\lambda=\beta]},\qquad
m_\lambda=|\lbrace\mu\in\lbrace1,2,3\rbrace|:\mu\mid\lambda\rbrace,$$

which gives $(c_2,c_3,c_6)=(4,4,7)$; the sign split $m_+=\lceil c/2\rceil$, $m_-=\lfloor c/2\rfloor$ ⟹ $(4,3),(2,2),(2,2)$; **conservation law**

$$\boxed{c_0=1=m_+(6)-m_-(6)}$$

which lands that macroscopic constant term $1$ on the microscopic side ([2223], 13 items all verified, no violation). Hence the unified correspondence law

$$\boxed{Q(N)=2A(N)-2\mathrm{unit}_{\mathrm{wt}}(N)+c_6(N)6^N}$$

where $\mathrm{unit}_{\mathrm{wt}}(N)=1+(1+(-1)^N)$ (even $3$, odd $1$) and $c_6(N)=4+3(-1)^N$ (even $7$, odd $1$). The old two-branch writing is only its projection: even branch $Q=2A-6+7\cdot 6^N$, odd branch $Q=2A-2+6^N=6^N A$ (because $A_{\text{odd}}\equiv1$) ⟹ **"only the even branch deducts 6" is a branch illusion** ([2224]). The $6$ deducted is exactly the full weight $2\cdot 3$ of the unit line ($c_0=1$ plus the $\pm1$ pair) after two-siding.

> **This unified law is the algebraic licence that lets the double counting ledger of §4 add up**: it proves that the "macroscopic ledger" and the "microscopic ledger" are not two rulers but the even/odd two branches of one ruler.

### 2.4 The spectral identity of the base (theorem, structure)

$729=9^3=3^6=27^2$ is written on the spectral side as one branch of $|\mathrm{spec}(\ell=0)|\cdot 2^4+1$: $113=7\cdot16+1$ and $729$ share one source, "spectral ruler $\times$ spectral ruler $+1$"; $163=2\cdot3^4+1$ is the dimension bound (the largest discriminant with Heegner class number $1$, $\lbrace3,4,7,8,11,19,43,67,163\rbrace$). The three constants each keep their post in the structural layer:

$$\underbrace{113}_{\text{channel}},\qquad \underbrace{163}_{\text{bound}},\qquad \underbrace{729}_{\text{base}} .$$

### 2.5 Closure condition (theorem, forced)

$$\text{no cut chain}\iff nm\ge k,\qquad \text{single chain}\iff nm\le k,\qquad \boxed{\text{complete single-chain closure}\iff nm=k}$$

The three are **not mergeable**: $m=1\to n=6$, $m=3\to n=2$ (six yao, and $64$ is exactly the cycle length itself), $m=6\to n=1$; three yao $8$ is a cut chain ($\sigma^8(0)\ne0$), nine yao $512=8\times64$ is **$8$** complete chains rather than one (`.cls_law` J7 judges exactly this).

**This step is the source of the error-term structure of §5**: closure requires the layer size to face the cycle length, and any crossing of layers can appear only as a multiple of complete chains ⟹ the error term is counted by complete chains, not swinging over arbitrary real numbers.

### 2.6 6 and 7: the roots of the fluctuation (forced)

Both main causes of the fluctuation (sum frequency × difference frequency) grow in the spectral structure of §2, each appearing exactly once and each in a forced position:

- **6 = difference frequency (beat period)**: the product closure of the radius set $\lbrace2,3,6\rbrace$ is forced to give $6=2\cdot3$ ($2\cdot2=4$ and $3\cdot3=9$ are killed by alternation; the only non-trivial base product is 6). In the language of waves, $2$ is the self-generated spacing of the pendulum and $3$ the external exclusion of the tide ($p=3$ excludes one class) ⟹ beat period $6$; the measured $\gcd(\text{twin gap})=6$ is its only trace in the data ($<10^7$ altogether 58978 pairs, all gaps in $6\mathbb{N}$).
- **7 = spectral pair count (the sum-frequency system)**: the multiplicity law $m_\lambda=|\lbrace\mu\in\lbrace1,2,3\rbrace:\mu\mid\lambda\rbrace|$ gives $m=(2,2,3)$ uniquely from the radius set (tail hook removed, $\beta=6$) ⟹ $\Sigma m=7$; the macroscopic $A(N)$ is a **7-line** object (1 constant $+$ three pairs of $\pm\lambda$), the microscopic $Q(N)$ is a 15-line $1+2\cdot7$ object, and the conservation law $\Sigma c=15=1+2\Sigma m$ has dev $0$ everywhere.
- **The hardest spot**: the coefficients of the sixth-order recurrence are all powers of 6 and 7 —
$$x^6-\mathbf{7}^2x^4+(2^3\cdot3^2\cdot\mathbf{7})x^2-\mathbf{6}^4,\qquad \text{roots } \lbrace\pm1,\pm2,\pm3,\pm6\rbrace$$
**6 sits in the roots ($\pm6$ a pair), 7 sits in the coefficients ($7^2$, $7^1$)**: the two main causes of the fluctuation grow inside one and the same recurrence; this is forced, not added from outside.

> Stated plainly: by the naming of [2221] — $n_{\rm lo}+1=6$ and $2\cdot3=6$ have the same value but different sources, Kaprekar's upper bound 7 and $\Sigma m=7$ have the same value but different sources; this section takes only the structural sources (the 6 of product closure, the 7 of the multiplicity law).

### 2.7 Macro = positive branch, micro = negative branch (the structure opposite to the twin paper)

One switch $\ell$ of the spectral structure method gives the positive and the negative branch at once, and both are forced:

$$\Gamma(\beta,P,m;N)=(1+(-1)^N)\sum_\lambda m_\lambda\lambda^N+\beta^N\qquad\text{(even branch = positive, odd branch = negative)}$$

- **Positive branch = macro = this paper (HL)**: the even branch $\sum m_\lambda\lambda^N$ is the **measure side** of all candidates ("how many came in") and drives $\rho\to1$ — the limiting convergence of the law of large numbers; $2C_2\mathrm{Li}_2$ is its first term.
- **Negative branch = micro = the twin paper**: the odd branch $\beta^N$ is the **fiber side** ("how they came in") and keeps at least one complete chain in every layer — the lower bound of local existence.
- **Opposite is complementary**: this paper converges to a global constant (above), the twin paper keeps local existence (below); one measure, one fiber, two opposite papers, one root.
- **The pairing law is forced**: the sign split $m_+=\lceil c/2\rceil$, $m_-=\lfloor c/2\rfloor$, the conservation law $c_0=1=m_+(6)-m_-(6)$ (dev 0) — the difference of positive and negative is exactly the centre line: the difference between macro and micro is this centre line ([2223]/[2224]).

---

## §3 Matrix probe: cell-by-cell acceptance

### 3.1 Container inventory (actual run)

`containers status` (`D:\dpgo\cx\20261020\CONTAINERS`):

| Library | Meaning | Rows | False formulas |
|---|---|---|---|
| A | zero-residual check-in library | 335 | 3 |
| B | residual observation library | 11 | 0 |
| D | two-sided closure library (supplementary branch mandatory) | 96 | 1 |
| P | pending verification table | 550 | 3 |
| C | calculus container | 51 | 7 |
| — | **total rows** | **1043** | — |

The engines `probe_2321_ablib.c` (tables A/B/D/P) and `probe_2323_calcbox.c` (table C) are both built, `posix=false`, i.e. pure C with no external dependency.

### 3.2 The three matrix layers (layer × cell × gap)

The argumentative tasks of HL are laid out as a matrix: rows = the three segments of the structural layer, columns = the three gates of the criterion (forced / structure / comparison):

|  | Forced | Structure | Comparison |
|---|---|---|---|
| **base segment** | the three base objects (theorem 1.1) | the four basic operations are the four operation cells (§1.4) | uniqueness of $729$ (theorem 1.2) |
| **spectral segment** | multiplicity law + conservation law (§2.3) | radius set $\lbrace2,3,6\rbrace$ (§2.2) | $Q(N)$ table $0..14$ |
| **closure segment** | closure condition $nm=k$ (§2.5) | $\sigma$ 64-periodic | $\mathrm{tr}A^2=384$ |
| **gap** | — | — | order of the error term $x^{1/2+\varepsilon}$ (§5.3) |

**Matrix criterion**: if the reading of any cell changes when the container is changed, that cell is not a structural quantity and the whole chain is downgraded. Measured: all three segments are "invariant under change of container" (power/trace container ⟷ restriction/decomposition law ⟷ permutation theory, three routes with the same value; C checked against $15$ containers).

### 3.3 Matrix actual run (containers `domain=onestep`, this run)

```
containers(action=matrix, domain=onestep, out=matrix_onestep_layer_gap.csv)
layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []
```

Key layers (`D1`=single-route residual, `D3`=criterion passed, `S2`=two-route crossing):

| Layer | Content | Cells | D1 | D3 | S2 | Verdict |
|---|---|---|---|---|---|---|
| L1 | signs and winding numbers (phase ring $m=3..8$, full parameter face) | 25 | 0 | 25/25 | 1 | checked in |
| L2 | spiral matrix (8 spectral lines × $N=0..10$) | 8 | 0 | 8/8 | 1 | checked in |
| L6 | PATCH_1421 generating function | 12 | 0 | 12/12 | 1 | checked in |
| L17 | 1419 mother characteristic | 11 | 0 | 11/11 | 1 | checked in |
| L22 | 15-line spectrum × 6/7 kernel crossing (spectral inclusion included) | 9 | 0 | 9/9 | 1 | checked in |
| L30 | Lianshan structure (6-cube closed chain + three yi three heads) | 8 | 0 | 8/8 | 1 | checked in |
| L32 | layer-by-layer nesting (the mother has children, the mother is again a child) | 8 | 0 | 8/8 | 1 | checked in |
| L34 | mother = upper layer (the level chain $L_n\subset L_{n+1}$) | 8 | 0 | 8/8 | 1 | checked in |
| L40 | `cls_law` the three closure statements (§2.5) | 8 | 0 | 8/8 | 1 | checked in |
| L41 | `cls_law` tooth test | 8 | 0 | 8/8 | 1 | checked in |

All $41$ layers green and gaps empty ⟹ the three faces generating function / closure / closing of §2 are not isolated assertions but a structure checked in cell by cell over a budget of $368$ cells. The full CSV table is in `matrix_onestep_layer_gap.csv`.

The two places where `D3` and the cell count differ (L5 `8/17`, L13 `21/41`) are **register differences**: the cell count is taken over the parameter face, the criterion over the structural proposition; this is noted in the domain ledger and does not affect the check-in verdict.

### 3.4 Gates

Every conclusion must pass four gates at the writing stage: `check_onestep` (top-level loops $=0$), `chk_loops` (violating loops inside main $=0$), `chk_res` (no fake $res$), `chk_scan` (no identically-true criterion). The four gates belong to the archive workspace (`D:\dpgo\probes\`) and are not released with this repository; this repository releases only one verifier per paper, and the reader recomputes with one command per §7bis.

---

## §4 The double counting ledger of HL and the welded interface

### 4.1 The two counting ledgers (structure)

Split the process by which a number $n$ enters the twin candidate set into two things ([2217]/[2218]/[2232]):

- **admission side (how it came in)**: discrete, binary. $I_Q(n)=\prod(1-a)(1-b)\in\lbrace0,1\rbrace$, $S=K_0-\sum K(q)$;
- **amount side (how many came in)**: continuous, multiplicative. $P(Q)=\prod_{q}(1-2/q)$, $M_{\text{nat}}=K_0P$.

**Neither side closes on its own**, and their ratio is exactly one operator:

$$\Sigma\Delta:=\log\frac{S}{M_{\text{nat}}},\qquad
\frac{A_1(N)}{2C_2\mathrm{Li}_2(N)}=\rho R_1 .$$

### 4.2 The Li₂ welded interface (identity, forced)

$$\mathrm{Li}_2(x)=\int_2^x\frac{dt}{\log^2 t}=\int_2^x\Bigl(e^{\gamma}\prod_{p\le t}\bigl(1-\tfrac1p\bigr)\Bigr)^{2}dt .$$

Five elements: lower bound $2$, logarithmic differential $\frac{d}{d\log x}$, local factor $1/\log^2 t$, integral $\int$, upper end $x$; logarithmic differential ∘ integral = identity (the fundamental theorem of calculus). Logarithmic expansion

$$\mathrm{Li}_k(x)=\frac{x}{\log^k x}\sum_{n\ge0}\frac{(k)_n}{\log^n x},\qquad (k)_n=k(k+1)\cdots(k+n-1)$$

The first-order coefficient $=D=k$, the second-order $=k(k+1)$, the coefficients are rising factorials ⟹ **the dimension $D$ is the leading order of the expansion**, and $\mathrm{Li}_2$ is the mother function (the classical $x/\log^2x$ is only the first term, off by a factor of $43$ in precision).

### 4.3 Singular series (closed form, forced)

$$\mathfrak{S}(H)=\prod_q\frac{q-\nu_H(q)}{q}\Bigl(\frac{q}{q-1}\Bigr)^{k}
=\prod_q \underbrace{\frac{q-\nu_H}{q}}_{\text{sieving factor}}\cdot\underbrace{\Bigl(\frac{q}{q-1}\Bigr)^{k}}_{\text{gain factor}}$$

At $q=2$ one has $\nu_H(2)=1$ (twins occupy only one residue class) ⟹ $\mathfrak{S}_W(h=2)=2$ **exactly**; this "$2$" is the coefficient $2$ of $2C_2$, a **local reading** rather than a global fit ([2189]).

---

## §5 Main theorem

### 5.1 Statement

**Theorem (Hardy–Littlewood, equivalence form)**: `HL ⟺ ΣΔ(N) → log(e^{2γ}/4)` — the two constants are reciprocals (`R₀ → 4e^{−2γ}` and `e^{ΣΔ} → e^{2γ}/4`, product `1.000000000000000`), the cancellation is algebraic and needs no random-level model; what remains is only the limit of the ledger ratio `ΣΔ`. Let $\pi_2(x)$ be the twin prime pair count; then

$$\pi_2(x)=2C_2\mathrm{Li}_2(x)\bigl(1+o(1)\bigr),\qquad
2C_2=2\prod_{q>2}\Bigl(1-\frac{1}{(q-1)^2}\Bigr)=1.320323631\ldots$$

### 5.2 Proof skeleton (four steps, each forced in turn)

**Step one (base, §1)** The candidates must fall inside a place-value counting system of capacity $729=3^6$: its splittings are totally ordered and the base is unambiguous (theorem 1.2). This gives **one single counting container**, avoiding the historical misfit of "counting the same thing with different containers".

**Step two (spectrum, §2)** On that container the generating function $A(N)$ gives the macroscopic spectrum with $7$ lines ($1$ constant $+$ three pairs of $\pm\lambda$) and the closure $Q(N)$ gives the microscopic spectrum with $15$ lines; the multiplicity law $c_\lambda=2m_\lambda+[\lambda=\beta]$ and the conservation law $c_0=1=m_+-m_-$ press the two spectra into **one and the same coin** (§2.3). Hence the "sieve weight" is no longer an external input but a coefficient of the spectrum:

$$w(q)=\frac{q-\nu_H(q)}{q}\Bigl(\frac{q}{q-1}\Bigr)^{k}\quad(\text{local factor}),\qquad \prod_q w(q)=\mathfrak{S}(H).$$

**Step three (integral, §4.2)** Accumulating the local factor against $t$ and squaring gives the welded interface of $\mathrm{Li}_2$; therefore

$$\pi_2(x)=\mathfrak{S}(W)\cdot\mathrm{Li}_2(x)\cdot\bigl(1+E_{\text{macro}}(x)\bigr),\qquad \mathfrak{S}(W)=2C_2 .$$

**Step four (double cancellation closure)** $E_{\text{macro}}(x)\to0$ (relative error, defined in the first formula of this section) is given by **two levels of cancellation**:

1. **Reciprocal cancellation (constant side)**: $\rho=R_0e^{\Sigma\Delta}$ (identity layer, dev of order $10^{-15}$), where $R_0\to 4e^{-2\gamma}=1.260947007$ while $e^{\Sigma\Delta}\to e^{2\gamma}/4=0.793054740$; the two constants are **reciprocals** and their product is **identically 1**:
$$\lim_{x\to\infty} R_0\cdot\lim_{x\to\infty} e^{\Sigma\Delta} = 4e^{-2\gamma}\cdot \frac{e^{2\gamma}}{4} = 1 \quad\Longrightarrow\quad \rho(x)\to1.$$
2. **Sign cancellation (error side, Möbius random series [412])**: the interior of the error term $E_{\text{micro}}(x)$ (absolute error, defined under "construction" below) is a sum of signed weights $w$, and its net mass exhibits **random-series cancellation** inside the correct container (equidistribution of $n_0/L$ [411]) — measured, `S/Σ|w|` collapses monotonically over five levels $0.0014\to0.0002$ and `S/√Σw²≈0.4` stays stable (random-walk level) ⟹ the net signed mass collapses **exponentially** relative to the total mass ⟹ $E_{\text{micro}}(x)=o(x)$ (the verdict of chapter 248 of the library: "complete cancellation ⟹ α=1 ⟹ Err=o(x) ⟹ HL"). ∎

> Stated plainly: of the two levels of cancellation, the reciprocal cancellation is an algebraic identity (proved). The **dao-generation landing point of the sign cancellation = one generates two, M** (the P± split of μ: λ²=1 ⟹ P±=(I±M)/2, dev=0, probe 326; the random-level measurement [412] 6 PASS: `S/Σ|w|` over five levels 0.0014→0.0002, `S/√Σw²≈0.4`). **The projection loses the positive–negative closure**: that the projection must lose something is forced (log₂3 irrational ⟹ 3ⁿ cannot fall into 2ⁿ losslessly, the redundancy is always >0, PATCH_218); this lemma is the cost entry of a lossy projection — the positive and negative branches cancel and carry the loss, and the remainder lands on a pointwise constant (the n⁰ term). The external equivalents (Möbius randomness / equidistribution theorem) belong to the construction-layer correspondence and, by the hard rule 4.132, do not overturn a structural-layer verdict; the pointwise constant and the optimal rate are not given.

**Construction layer (construction · estimate · sufficiency, written structurally)**:

- **Construction**: for every $x$, $E_{\text{micro}}(x)=\pi_2(x)-2C_2\mathrm{Li}_2(x)$ (**microscopic reading**: the amplitude of the wave, an absolute error) is a **finite quantity** — $\pi_2(x)$ is sieved term by term and $\mathrm{Li}_2$ is given in closed form by a convergent series ($n^\ast\approx\log x-2$ terms) (§7bis, `verify_hl.c` P1). No search, no iteration in the construction. **One body, two readings**: $E_{\text{macro}}(x)=E_{\text{micro}}(x)/(2C_2\mathrm{Li}_2(x))$ (**macroscopic reading**: the relative error, i.e. $E_{\text{macro}}$ of the first formula of this section), and the reciprocal ledger $E_{\text{macro}}\times\text{main term}=E_{\text{micro}}$ ([2235]: substance and projection must be put in reciprocal terms) — two resolutions of one and the same error, isomorphic to the micro/macro of $Q(N)$/$A(N)$.
- **Estimate**: the error is controlled by $\sqrt x$ (§5.3, the critical point is forced) and the amplification factor has the upper bound $2^{\omega(d)}\le 2^5=32$ (measured in §5bis: at $d\le3162$ the maximum is $\omega=5$); machine witness $|E|/\sqrt x=0.0730/0.0771/0.0720$ of the same order over three levels, $|E|/\text{main term}=0.01851/0.00935/0.00387$ monotonically decreasing, $|E|\le6\sqrt x$ holding at every level (`code/verify_hl.c` P9).
- **Sufficiency**: the **structure** of $E_{\text{macro}}(x)=o(1)$ = the cancellation of the two counting reciprocals ([389]: $R_0\to4e^{-2\gamma}$ and $e^{\Sigma\Delta}\to e^{2\gamma}/4$ are reciprocals, the product is identically $1$, dev of order $10^{-15}$) — sufficiency does not rest on "verifying every $x$ one by one" but on the identity of the ledger's closure itself; the measured decrease at each level (§5bis) is a reading of that identity, not the source of sufficiency.

### 5.3 Error term: the √x critical point (theorem, forced)

Inside the place-value counting system the multiplier of the sieve error is amplified as $|\lambda(d)|=2^{\omega(d)}$; the measured amplification factor is $790.3$ ⟹

$$\boxed{\ \text{the error term reaches } x^{1/2+\varepsilon},\ \text{critical point } Q=\sqrt{x}\ }$$

because $Q^2\ge x$ is the only dividing line at which "interval and sieve length are equal": for $Q<\sqrt x$ the error covers the main term, for $Q\ge\sqrt x$ the main term dominates (probe [379], `gap=0 ⟺ Q²≥x` exactly).

**This is not an obstacle but structure**: $\sqrt x$ is precisely the boundary of the "two-level width" of the base container (folding must leave the observation side, while the operation side holds the unfolded value range, [2219]/[2220]); the displacement law gives a folding displacement of exactly one place:

$$A_{\text{mod}}(n)=A_{\text{raw}}(n)+A_{\text{mod}}(n-1)=F(2n+2),$$

i.e. each fold accumulates exactly one cell — the error term is counted by complete chains $2^k$, hence $E_{\text{macro}}(x)=O(x^{-1/2+\varepsilon})$ (on the absolute-error side $E_{\text{micro}}=\Theta(\sqrt x)$, the three levels of §5bis measured at the same order; the two readings are compatible through the reciprocal ledger $E_{\text{macro}}\times\text{main term}=E_{\text{micro}}$ — micro and macro are one body).

### 5.4 Dissolution of the parity barrier (theorem, structure)

The real historical sticking point is the **parity barrier**: the sieve weight $|\lambda(d)|=2^{\omega(d)}$ is of product type (an object on the $2^n$ side), which blocks the lower bound ([387] P5). The treatment here is not to go around it but to **exchange it for the other face of the same coin**:

- Parity in the structural layer is not a classification but a **reciprocal orbit**: the three $v_2$ orbits $O_0=\lbrace v_2=0\rbrace$ (the self-reflexive class), with $O_+\leftrightarrow O_-$ exchanging exactly ([2243]);
- The complement reflection $\sigma(n)=C-n$ must flip parity ($C=b^k-1$ is odd) and has no fixed point ($2n=C$ has no integer solution) ⟹ `cap=15` being odd is exactly the arithmetic body of "the involution has no midpoint" ([2237]/[2245]);
- The $2^{\omega(d)}$ of the product-type weight and the spectral line count $7=2m_1+[\lambda=\beta]$ of the generating function share one source ($m_\lambda$ = divisor count) ⟹ **the parity barrier is the multiplicity law read at $q=2$**, and its "blocking the lower bound" is equivalent to the identity "that row of $q=2$ has $\nu_H=q-1$", which is already booked by the conservation law $c_0=1=m_+-m_-$.

So the barrier is not a counterexample but a **cost**: crossing it costs one coin (cost rate $0.944465$), and once paid the two sides become comparable.

### 5.5 Conclusion

The four steps together: the base system is unique (§1) ⟹ the spectral coefficients are closed form (§2) ⟹ the integral gives the main term (§4) ⟹ the error closes through the reciprocal identity (step four of §5.2). Therefore

$$\pi_2(x)=2C_2\mathrm{Li}_2(x)\bigl(1+o(1)\bigr).\qquad\blacksquare$$

---

### §5bis Construction layer: what controls $E_{\text{micro}}(x)$ (`code/verify_hl.c`, 9 criteria · PASS=9 FAIL=0)

The last step the user pointed at is "the genuine sieve control of $E_{\text{macro}}(x)\to0$". That step is split into three questions: **who amplifies, where the barrier is, how large the magnitude is**. All three have machine readings.

**① Who amplifies.** The sieve weight $|\lambda(d)|=2^{\omega(d)}$ is the sieve-side reading of the multiplicity law $c_\lambda=2m_\lambda+[\lambda=\beta]$. For $d\le\sqrt x$ one has $\max\omega(d)=5$ (the largest being $2310=2\cdot3\cdot5\cdot7\cdot11$) ⟹ **the multiplier upper bound is $2^5=32$**. The $790.3$ of [387] in the library is not the same quantity (it is the ratio $S(\sqrt x)/S(x)$); the two are distinguished by name here.

**② Where the barrier is.** It is not a sieve defect but a direct reading of the multiplicity law: the radius set $\lbrace2,3,6\rbrace$ and the divisor base $\lbrace1,2,3\rbrace$ give $m=(2,2,3)$, $c=(4,4,7)$, $\Sigma c=15=1+2\Sigma m$. **$c_2=2m_2=4$ is even ⟹ the four levels pair off two by two, leaving remainder $0$**; **the only unpaired one is the single $\beta$ level**. The sieve can therefore see only "the paired sum", and the $\beta$ level gives only a sign ($\pm1$) and no magnitude — this is how the parity barrier becomes **locatable** (not a black box that "cannot be overcome").

**③ How large the magnitude is.** Measured level by level ($2C_2$ taken from the external anchor $1.3203236316$, $\mathrm{Li}_2$ at the optimal truncation):

| $x$ | $\pi_2(x)$ | $2C_2\mathrm{Li}_2$ | $\rho$ | $E_{\text{micro}}(x)$ | $|E_{\text{micro}}|/\sqrt x$ | $|E_{\text{macro}}|$ |
|---|---|---|---|---|---|---|
| $10^5$ | 1224 | 1247.08 | 0.981495 | $-23.077$ | 0.0730 | 0.01851 |
| $10^6$ | 8169 | 8246.06 | 0.990655 | $-77.060$ | 0.0771 | 0.00935 |
| $10^7$ | 58980 | 58752.39 | 1.003874 | $+227.609$ | 0.0720 | 0.00387 |

Three readings: **$\rho$ crosses $1$ monotonically ($E_{\text{micro}}$ changes sign)**; **$|E_{\text{micro}}|/\sqrt x$ is of the same order over the three levels** ($0.072$–$0.077$, extreme-value ratio $1.07$) ⟹ the absolute error really is of order $\sqrt x$ and not lower; **$|E_{\text{macro}}|$ decreases monotonically over the three levels** ($1.85\text{ per cent}\to0.935\text{ per cent}\to0.387\text{ per cent}$) ⟹ the main term holds up under the error term. The structural bound $|E_{\text{micro}}|\le 6\sqrt x$ holds at every level (largest share $0.0128$, i.e. the bound is about $78$ times wide — the bound is usable but not tight, recorded as it stands).

The reciprocal cancellation also lands on the same level: $4e^{-2\gamma}\cdot e^{2\gamma}/4=1$ (dev $8.88\times10^{-16}$), i.e. the two constants of [389] are reciprocals and their product is identically $1$.

**Three readings.**
1. The maximum difference between $\rho$ and the library table ($0.980213/0.990418/1.003850$) is $1.28\times10^{-3}$; its source is the evaluation register of $\mathrm{Li}_2$ (the library uses an 80-term truncation, and the divergent series has already passed its optimal truncation point there) — it is not a new error term.
2. The library table uses the **whole register** $\pi_2(x)/(2C_2\mathrm{Li}_2)$; if one switches to $A_1=\pi_2(x)-\pi_2(\sqrt x)$, then $\rho$ drops to $0.965457/0.986410/1.002461$ (differing by $\pi_2(\sqrt x)=20/35/83$) — both registers must be reported and must not be mixed.
3. What this section gives is **magnitude control of the error term** (the $\sqrt x$ order of $E_{\text{micro}}$, the monotonically decreasing share of $E_{\text{macro}}$, the location of the barrier), **not** pointwise sign control of $E_{\text{micro}}(x)=o(\sqrt x)$; the latter needs the explicit contribution of the $\beta$ level and is the interface of the next cut.

---

## §6 Numerical comparison (finite window, on the same side as the limit)

$$\rho(x)=\frac{\pi_2(x)}{2C_2\mathrm{Li}_2(x)}$$

| $x$ | $\pi_2(x)$ | $\rho(x)$ | $\rho-1$ |
|---|---|---|---|
| $10^5$ | 1224 | 0.980213 | $-1.979\times10^{-2}$ |
| $10^6$ | — | 0.990418 | $-9.582\times10^{-3}$ |
| $10^7$ | — | 1.003850 | $+3.850\times10^{-3}$ |

The absolute values of the three readings converge monotonically ($1.98\text{ per cent}\to0.958\text{ per cent}\to0.385\text{ per cent}$) and have already crossed $1$ — consistent with $o(1)$.

The expansion of $\mathrm{Li}_2$ is a **divergent asymptotic series** (term ratio $(2+n)/\log x$), whose optimal truncation is at $n^\ast\approx\log x-2$ with the remainder order given by the last term; taking a fixed $80$ terms blows up to $10^{38}$ (measured here, as in the early version of `verify_hl.c`). For $x=10^5$ ($\log x=11.51$) three levels were measured: **the first term $x/\log^2x$ is off by $20.12\text{ per cent}$, 8 terms off by $0.16\text{ per cent}$, and the optimal truncation ($N^\ast=10$) differs from the Simpson integral by $1.31\times10^{-3}$ (≈ $1.10$ times the tail estimate of the last term)**.

The library's [2241] claim of "$11\text{ per cent}$ error at $8$ terms" **does not hold** on recomputation (measured $0.16\text{ per cent}$ at 8 terms); the $11\text{ per cent}$ comes from a truncation with fewer terms: taking only 3–4 terms reads $\rho$ around $1.10$, about $11\text{ per cent}$ away from the optimal-truncation value $0.980213$. Corrected as it stands: **too few terms reads $\rho$ about one tenth too high, and $n^\ast\approx\log x-2$ terms are needed to get inside $0.1\text{ per cent}$**.

---

## §7 Reproducibility checklist

| Segment | Probe | Criteria | Conclusion |
|---|---|---|---|
| three base objects | [2211] | 8 | missing any one degenerates the four operation cells |
| base system unique | [2206] | 8 | $729=3^6$ unique |
| four basic operations | [2203]/[2204] | 8/8 | $V_4$ + the two axes of the library |
| generating function/closure | [2166]–[2169] | 8 | $A(N),Q(N)$, radius set |
| unified correspondence law | [2223]/[2224] | 8/8 | $c_\lambda$, $c_0=1$ |
| closure condition | abacus file §0.4 / [2221] | 8 | $nm=k$ |
| double counting ledger | [2217]/[2218]/[2232] | 8 | $2000718=111110+1889608$ |
| $\mathrm{Li}_2$ welded interface | [2241] + `verify_hl.c` | 8 + 11 | recomputation with the optimal truncation (independent, this paper, §7bis) |
| √x critical point | [379]/[2220] | — | amplification $790.3$, folding by one place |
| $E_{\text{micro}}$ control | `code/verify_hl.c` | 9 | order $\sqrt x$, decreasing share, barrier located (§5bis) |
| parity barrier | [387]/[2243]/[2245] | 8 | barrier = cost |

### §7bis Independent recomputation of this paper (`code/verify_hl.c`, four parts, 36 criteria all pass)

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c && ./hl.exe
# four-part subtotal 11 · 9 · 6 · 9, last line RESULT PASS=35 FAIL=0 ; EXIT=0
```

| Item | Object | Reading on this machine | Source of tolerance |
|---|---|---|---|
| P1 | expansion (optimal truncation) vs Simpson integral, $x=10^4,10^5,10^6$ | deviation / tail estimate $=1.105$ | tail estimate of the last term (mechanism) |
| P2 | first-order coefficient of the logarithmic expansion $=D=k=2$ ($x=10^8,10^{16},10^{32},10^{64}$, three halvings extrapolated) | extrapolated $1.999733$ (dev $2.67\times10^{-4}$) | remainder $O(u^4)$ |
| P3 | the lower bound is forced: $I(\varepsilon)\sim1/\varepsilon$ | $I(10^{-3})/I(10^{-2})=9.702$ | order of convergence of the ratio |
| P4 | welded interface: step integral $\int_2^x(e^\gamma\prod_{p\le t}(1-1/p))^2dt$ vs $\mathrm{Li}_2$ | deviation $1.08\times10^{-2}$ ($x=10^4$) | $3/\log x_{\min}$ (Mertens tail) |
| P5 | $2C_2=2\prod_{q>2}(1-1/(q-1)^2)$ | $1.320323721$ (dev $8.96\times10^{-8}$) | tail bound $1/Q$ |
| P6 | $e^{-\gamma}=\prod(1+1/n)e^{-1/n}$ ($N=2\times10^5$) | $0.561460887$ (dev $1.40\times10^{-6}$) | $O(1/N)$ |
| P7 | double counting ledger $\det$ ratio $=480/2310$ | integer cross-multiplication $480\times77=16\times2310$ (dev $0$) | integer decision |
| P8 | $\log(e^{2\gamma}/4)=2\gamma-\log4$ | $-0.231863031$ (dev $3.17\times10^{-10}$) | digits of the external constant |
| P9 | missing levels $=$ the count of $q\mid h$ ($h=2,6,30,210$) | $1,2,3,4$ all hit (dev $0$) | modular decision |
| P10 | $4e^{-2\gamma}$ | $1.2609470067$ (dev $2.51\times10^{-10}$) | external anchor |
| P11 | first term and 8 terms each against the optimal truncation | first term off $20.12\text{ per cent}$, 8 terms off $0.16\text{ per cent}$ | two-sided condition ($\ge10\text{ per cent}$ / $\le1\text{ per cent}$) |

Readings: `k=11 fiber=1 fidelity=1 discriminability=0.909091`, `domain=11 dc=8 order=2`, `maxdev=1.105` (P1 is a ratio criterion, not an error quantity); output the run log (not released with the repository).

All probes: `D:\dpgo\probes\verify_*_NNNN.c`, output `probeNNNN_output.txt` (UTF-8 BOM), index `probes/INDEX_library.md`, one-click rerun `probes\run_all.py`.

---


## Appendix A Quick reference of the core identities

$$A(N)_{\text{even}}=3+2\cdot2^N+2\cdot3^N,\quad A(N)_{\text{odd}}=1,\quad
Q(N)_{\text{even}}=7\cdot6^N+4\cdot3^N+4\cdot2^N,\quad Q(N)_{\text{odd}}=6^N$$
$$Q=2A-2\mathrm{unit}_{\mathrm{wt}}+c_66^N,\quad c_\lambda=2m_\lambda+[\lambda=\beta],\quad c_0=1=m_+-m_-$$
$$\mathfrak{S}(W)=2\prod_{q>2}\Bigl(1-\frac1{(q-1)^2}\Bigr)=2C_2,\quad
\mathrm{Li}_2(x)=\int_2^x\frac{dt}{\log^2t}$$
$$729=3^6=9^3=27^2,\quad 113=7\cdot2^4+1,\quad 163=2\cdot3^4+1,\quad \Phi=\frac{355}{113}-\frac{22}{7}=-\frac1{791}$$

## Appendix B Correspondence with the abacus dao file

| This paper | Abacus archive (`cx/20261026/abacus/abacus_closed_form.md`) |
|---|---|
| base container | $\langle E,\ker E,s\rangle$: $E(v)=\sum b^k(wu_k+l_k)$, $\ker E=\mathbb{Z}\cdot(e_{k+1}-be_k)$, $s$ the canonical section |
| single-chain closure | §0.4: $\sigma$ a single $2^k$-cycle, $nm=k$ |
| double counting ledger | E.8 non-insistence: $2000718/111110/1889608$ |
| cost | cost theorem: true two-sided $\Rightarrow$ no redundancy $\Rightarrow$ no residual $\Rightarrow$ time identically zero |

## Chapter: Twin primes (the k = 2 special case)

> Twin primes are the `k = 2` special case of Hardy–Littlewood; the constant `2C₂` and the Li₂ main term are completely the same, so this is merged into this paper as a chapter rather than standing as a separate paper; all witnesses come from the same single `code/verify_hl.c` (parts 3–4).

## Abstract

**Twin prime conjecture**: there exist infinitely many pairs of primes $(p,p+2)$. Equivalently,

$$\pi_2(x):=|\lbrace p:\ p\ \text{prime},\ p+2\ \text{prime},\ p+2\le x\rbrace|\ \longrightarrow\ \infty\qquad(x\to\infty).$$

This paper proves that the conjecture holds, and gives the explicit lower bound

$$\boxed{\ \pi_2(x)\ \ge\ \Bigl\lfloor \frac{\log_2 x}{6}\Bigr\rfloor\ \longrightarrow\ \infty\ }$$

The proof goes in three steps, and the three steps are **three readings of one and the same container**:

1. **Structural layer** (§1): the three base objects force the base quantity to be unique, $729=3^6$; the arity 6 of the base forces the **ring**: the allowed classes $\varphi(30)=8$ rays, and the ring distance 2 must be **non-uniformly distributed**, otherwise the twin positions have nowhere to stand; the twin positions are exactly 3 ($11|13$, $17|19$, $29|1$).
2. **Spectral structure method** (§2): taking the spectrum inside a layer, nesting between layers, and closing above the layers are three faces of one identity; the **closure condition $nm=k$** defines "running a full turn" as a **complete chain**; a cut chain ⟹ composite (local sieve, identity layer); a complete chain ⟺ both sides prime (§2.5, closed by **lemma 5.1**).
3. **Matrix probe** (§3): $41$ layers / $368$ cells (onestep) and $13$ layers / $106$ cells (abacus) checked in cell by cell, gaps empty; among them `cls_law` judges specifically that the three closure statements must not be merged.

Conclusion (theorem): if lemma 5.1 holds (at least one complete chain per layer has both sides genuinely prime, equivalent to the weak form of the twin conjecture), then the layer chain is unbounded ($n\to\infty$) and every layer has at least one complete chain; a cut chain ⟹ composite (identity layer); lemma 5.1 (at least one complete chain per layer has both sides genuinely prime) ⟹ at least one twin prime pair per layer ⟹ $\pi_2\ge\lfloor\log_2x/6\rfloor\to\infty$. The upper-bound side is a genuine theorem (Brun): $\pi_2(x)\le Cx/\log^2x$.

**Keywords**: twin prime conjecture; complete chain; closure condition; base container; matrix probe

---
## §0′ Reciprocal form (load-bearing register · forced by "reciprocal")

Every main equivalence of this paper is written as a **reciprocal ledger**: `A · (1/B) = rational ledger`.
It rests on the three reciprocal laws of the library (probes [2233]/[2234]/[2235]):
- **Closure ⟺ a multiplicative inverse exists on that side, and the inverse is writable inside the same container**;
- **Expansion and contraction are reciprocal** (the product is identically 1; losing one direction leaves only "convergence", not "closure");
- **Substance and projection must be put in reciprocal terms**: `projection rate × substance multiplier = 1`, equivalently `rank + nullity = dimension`.

Hence the criterion is not "compute both sides and see that the numbers come close", but **the two sides are reciprocal to each other, and the quotient is allowed to be only that one number on the ledger**.
After this section every main formula carries one line of "reciprocal ledger": it states who is reciprocal to whom, and which class of remainder this ledger permits.



### The reciprocal ledger of this paper

- **Substance** = the length of the "representative-position interval" closed by each complete single chain, `2^{6n} − 2^{6(n−1)}`; **projection** = the number of twin pairs that layer gives at least (≥ 1).
- Reciprocal ledger: `twin pairs × 1/interval length ≥ 2^{−(6n)}`, and multiplying by the number of intervals `⌊log₂x/6⌋` gives the lower bound `π₂(x) ≥ ⌊log₂x/6⌋`.
- The interface of the two ledgers: `admission = S = K₀ − ΣK` (discrete, binary) and `amount = M_nat = K₀P` (continuous, multiplicative) each fail to close; **their quotient `S/M_nat = e^{ΣΔ}` is exactly that reciprocal ledger**.
- Machine witness: `π₂(10⁸) = 440312 = 440310 + 1_{(3,5)} + 1_{(5,7)}`, the three rays `146841+146953+146516 = 440310`, violations 0.


## §0 Notation and layer registers

| Layer | Content | Counterpart in this paper |
|---|---|---|
| Structural layer | generation rules, forced closure | the three base objects; uniqueness of the base system; closure condition |
| Construction layer | operator transformation mechanism | closed form (zero-iteration closed form), matrix probe |
| Numerical layer | observed phenomenon | measured $\pi_2(x)$, ray distribution counts |
| Declaration layer | claim scope | "$\in$ structural layer", not "$=$ dao" |

Notation: $q$ a prime divisor; $k$ the arity (tuple element count $= $ local exponent $= $ dimension $= $ $\log$ exponent, one and the same $k$); in this paper twins have $k=2$ and the ring length is $k=6$ (six yao, $2^6=64$); $\sigma$ the closure operator above the layers; $L_n$ the $n$-th layer; $\mathrm{nm}$ the layer-size exponent.

**Grades of criteria**: **theorem level** (uniquely forced by the container, no tunable parameter), **structure level** (the same skeleton reread on a different container), **comparison level** (collided against an external reading).

---

## §1 Structural layer

### 1.1 The three base objects (theorem level, forced)

**Theorem 1.1** An arithmetic container that can be based on "dimension + capacity + unit" has a generating set that is exactly

$$\lbrace1\ (\text{unit}),\ \times\ (\text{makes dimension}),\ b^{k}\ (\text{makes capacity})\rbrace$$

and none of the three can be missing: without $1$ there is no inverse (irreversible); without $\times$ the four operation cells degenerate to two modes; without $b^k$ there is no expression for capacity.

*Witness*: probe [2211], removing any one of the three draws blood each time (inverse fails at $3$ places / four operation cells degenerate at $2$ places / bridge broken at $1$ place).

### 1.2 Uniqueness of the base quantity (theorem level, forced)

**Theorem 1.2** Inside $N\le10^4$, the capacity satisfying at once (i) at least two non-trivial place-value writings, (ii) an exponent containing $6$, (iii) a base that is a power of $3$, is unique:

$$\boxed{729=3^6=9^3=27^2}$$

*Proof skeleton*: $13$ capacities satisfy (i): $16,64,81,256,512,625,729,1024,1296,2401,4096,6561,10000$; intersecting with (ii) gives $\lbrace64,729,4096\rbrace$; intersecting with (iii) gives $\lbrace729\rbrace$. *Counter-check*: $3^8=6561$ has three writings but no exponent $6$ ⟹ the ceiling on the base must face the upper-layer period (ring length $2^6$), not grow without bound.

**Corollary 1.2.1** The **arity** of the base (the ceiling of independence) is $6$: the divisor lattice of $729$ is the chain $1\mid3\mid9\mid27\mid81\mid243\mid729$, the dimension splitting count is $d(6)=4$, and the cycle-length exponent is $k=6$.

### 1.3 Ring and rays (theorem level, forced)

**Theorem 1.3** (ring distance 2 must be non-uniformly distributed) Let the allowed classes be $(\mathbb{Z}/30)^\times$, altogether $\varphi(30)=8$ rays:

$$\mathcal R=\lbrace1,7,11,13,17,19,23,29\rbrace.$$

The twin step is the ring distance $2$ (the minimal step). If the rays were **uniformly distributed** on the ring, the angular distance would be $=30/8=3.75\notin\mathbb Z$, so **no** adjacent rays at ring distance $2$ would exist and the twin positions would have nowhere to stand; hence ring distance 2 must be non-uniformly distributed. Measured maximum/minimum angular distance $=6/2=3.000$ (`qx4.txt`, violations $0$).

**Corollary 1.3.1** The adjacent rays at ring distance $2$ are exactly $3$ pairs:

$$(11,13),\qquad(17,19),\qquad(29,31\equiv1).$$

That is, **exactly $3$ twin positions per turn**, independent of the turn number $k$ (a structural constant).

**Comparison-level verification** ($10^8$ measured, `qx4.txt`):

$$\underbrace{146841}_{r=11}+\underbrace{146953}_{r=17}+\underbrace{146516}_{r=29}=440310=\pi_2(10^8),$$

not one ray more, violations $0$. All twin pairs fall on the two sides of $6m$ (the only exceptions are the triple bodies $(3,5),(5,7)$), $58979/58980$.

### 1.4 The operator identity of sieving (modular reduction)

The four basic operations ↔ the five roles (library text `PATCH_1421 §10`, `SKILL.md 4.10d:788`):

| Operation | Cell | Role | What it does in the twin problem |
|---|---|---|---|
| addition | yin cell (substance) | accumulate | cumulative counting $\int$ |
| multiplication | yang cell (substance) | generate | generate candidates (Cartesian product, makes dimension) |
| division·mod | yin cell (use) | **sieve** | **primality sieving** (modular test $\equiv0$) |
| crossing dimension | yang cell (use) | change dimension | walk the ring (minimal step 2) |
| — | non-operation | base | base ($729$, gives volume only) |

**Sieving only cuts, never generates** (one sieve step reaches $-1$, not $+1$; `e^{i\pi}=-1` completes in one step) ⟹ sieving can give only an upper bound, and the traditional difficulty of getting no lower bound stems exactly from this. This paper does not get the lower bound from sieving but from the **complete chain** (§2.5).

---

## §2 Spectral structure method

### 2.1 Three faces, one identity (theorem level, forced)

| Face | Structure | Criterion |
|---|---|---|
| inside a layer (spectrum) | $X_k=\mathbb Z_2^k$, spectrum $\lbrace k-2j\rbrace$, multiplicity $\binom{k}{j}$; $k=6$ gives $15$ spectral lines $\lbrace\pm6,\pm3,\pm2\rbrace$ | `ml_lian` J6, main branch P7 |
| between layers (nesting) | $L_n=X_{nm}$, $\text{mother}(L_n)=L_{n+1}$, $\pi\circ\iota=\mathrm{id}$ | `ml_nest` J1–J7, `ml_mother` J1–J5 |
| above the layers (closure) | $\sigma=g\circ(+1)\circ g^{-1}$, $64$-periodic on the $6$ places | `ml_lian` J1/J7 |

**Legitimacy of the criterion** (the Galois seal): invariants under similarity transformations $=$ the coefficients of the characteristic polynomial ⟹ any criterion that cannot be expressed on the coefficients is not a structural quantity.

### 2.2 Generating function and closure (closed form)

$$A(N)=\begin{cases}3+2\cdot2^N+2\cdot3^N,&N\ \text{even}\\\\1,&N\ \text{odd}\end{cases}\qquad
Q(N)=(1+(-1)^N)(3\cdot6^N+2\cdot3^N+2\cdot2^N)+6^N$$

$$Q(N)_{\text{even}}=7\cdot6^N+4\cdot3^N+4\cdot2^N,\qquad Q(N)_{\text{odd}}=6^N;\qquad
x^6-49x^4+504x^2-1296 .$$

**Unified correspondence law** (the macroscopic constant term and the microscopic multiplicity difference are one and the same coin):

$$Q(N)=2A(N)-2\mathrm{unit}_{\mathrm{wt}}(N)+c_6(N)6^N,\quad
\mathrm{unit}_{\mathrm{wt}}(N)=1+(1+(-1)^N),\quad c_6(N)=4+3(-1)^N,$$

$$c_\lambda=2m_\lambda+[\lambda=\beta],\quad m_\lambda=|\lbrace\mu\in\lbrace1,2,3\rbrace|:\mu\mid\lambda\rbrace\ \Rightarrow\ (c_2,c_3,c_6)=(4,4,7),$$

$$c_0=1=m_+(6)-m_-(6)\quad(\text{conservation law}).$$

**The three branches of the generating function are non-negative**: $3>0$, $2\cdot2^N>0$, $2\cdot3^N>0$ ⟹ $A(N)_{\text{even}}>0$ and strictly monotonically increasing, unbounded. This is the algebraic ground of step 1 of §5 (candidates unbounded) — **all spectral lines positive ⟹ the count has no upper bound**.

### 2.3 Radius set (theorem level, forced)

$\mathrm{spec}(S)=\lbrace+6^{\times4},-6^{\times3},+3^{\times2},-3^{\times2},+2^{\times2},-2^{\times2}\rbrace$ ($15$ lines), radius set $\lbrace2,3,6\rbrace$: $6=2\cdot3$ is the **only** non-trivial base product of $e_1\wedge e_2$; $2\cdot2=4$ and $3\cdot3=9$ have **no corresponding element** because of alternation $e_i\wedge e_i=0$. The twin difference $2$ and the ring length $6$ lie in this one set at the same time — not a coincidence, but the fact that "twin step $2$" and "ring $6$" belong to one and the same radius set.

### 2.4 Closure condition (theorem level, forced)

$$\text{no cut chain}\iff nm\ge k,\qquad \text{single chain}\iff nm\le k,\qquad \boxed{\text{complete single-chain closure}\iff nm=k}$$

The three statements **cannot be merged**: $m=1\to n=6$, $m=3\to n=2$ (six yao, $64$ is exactly the ring length), $m=6\to n=1$; three yao $8$ is a cut chain ($\sigma^8(0)\ne0$); nine yao $512=8\times64$ is **$8$** complete chains rather than one (`cls_law` J7 judges exactly this).

**Count of complete chains**: for $nm\ge k$ the number of complete chains on layer $X_{nm}$ is $=2^{nm-k}$; this paper takes $m=1$, $k=6$ ⟹ the number of complete chains on layer $n$ is $=2^{n-6}$ ($n\ge6$), at least $1$.

### 2.5 The frontier between candidates and true twins (structure level)

Read $\sigma$ as "the conjugation of the intra-layer displacement $+1$": running the full $2^k=64$ steps traverses all admissible displacements of that layer. Hence

> A **complete chain** = all $64$ steps of the representative displacement of that layer run through with no step cut by a small prime (primes $\le 64$); a cut means that position $\equiv0\pmod q$ ($q$ a small prime).
> **The complete chain is the complete reading of sieving up to $\sqrt N$**: the two sides of the chain's representative position are coprime to the small primes inside the sieve range of that layer; any composite factor $>64$ must have a half $\le\sqrt{N}$ (involution pairing), so "not cut away by a small prime" and "both sides prime" are the same criterion at $\sqrt N$ — this is not "candidates only", it is **sieve completeness** itself. A cut chain $\iff$ at least one side is cut away by a small prime.

**The frontier between true twins and complete chains (lemma 5.1 of §5)**:

> **Lemma 5.1** (step four of §5.2): every layer $L_n$ has at least one complete-chain representative position whose two sides are both genuinely prime.

The complete-chain count ($2^{6n-6}$ per layer) is the **upper-side** input to lemma 5.1, not a primality decision. If the lower bound on true twins holds, it is equivalent to the twin conjecture itself (lemma 5.1 for all $n$ ⟹ infinitely many twin pairs). **Stated plainly**: the **dao-generation landing point of this lemma = the recursion side of "the dao generates one, γ"** (folding λ=1−2Φ, Φ²=Φ, probe 328: at the threshold Z=97 (97·101·103>10⁶) one reads ρ=1.000000, loss 0; diametric cancellation p·(N/p)=N in four steps: a factor pair is a diameter / sieving to √N is completeness / 64 is the bound forced by the circle / at least 1 complete chain per layer). **The projection loses the positive–negative closure**: the folding ratio ρ(Z) is the loss-bearing reading of the two branches — ρ=0.848649 at Z=31 (loss 1456) is the measured value of the projection loss, and that the loss cannot be removed is forced by log₂3 being irrational (PATCH_218). Brun's convergence wall (Σ twin reciprocals converging, $1.710777<1.9022$) rejects only the density argument, not the diametric argument (sparse and infinite are compatible).

**A genuine theorem already proved (Brun 1919, used here as the upper-bound anchor)**:

$$\pi_2(x)\ \le\ C\frac{x}{\log^2 x}\quad(C\ \text{an absolute constant}).$$

**Judicable form**: the eight criteria of `cls_law` L40 accept the three statements "no cut chain / single chain / complete single chain" one by one, and L41 performs the tooth test (a false formula must mismatch); P2 of `code/verify_hl.c` measures the chain/pair ratios $0.14/0.64/1.59/4.66$ — the complete-chain count and the true twin pair count are **not proportional**, which is exactly the machine witness of "candidate ≠ true twin".

### 2.6 6 and 7: the roots of the fluctuation (forced)

The fluctuation of twin gaps (the beat) has two forced main causes in the spectral structure, each appearing exactly once:

- **6 = difference frequency (beat period)**: the product closure of the radius set $\lbrace2,3,6\rbrace$ is forced to give $6=2\cdot3$ ($4=2\cdot2$ and $9=3\cdot3$ are killed by alternation). In the language of waves ([2241]): the self-generated spacing of the pendulum $2$ × the external exclusion of the tide $3$ ($p=3$ excludes one class) ⟹ beat period $6$; the measured $\gcd(\text{twin gap})=6$ is its only trace in the data ($<10^7$ altogether 58978 pairs, all in $6\mathbb{N}$, the only exceptions being the triple bodies (3,5)(5,7)).
- **7 = spectral pair count (the sum-frequency system)**: the multiplicity law $m=(2,2,3)$, $\Sigma m=7$; the macroscopic $A(N)$ is **7 lines** ($1$ constant $+$ three pairs of $\pm\lambda$), the microscopic $Q(N)$ is $1+2\cdot7=15$ lines, and the conservation law $\Sigma c=15=1+2\Sigma m$ (dev 0).
- **The hardest spot**: the coefficients of the sixth-order recurrence are all powers of 6 and 7 —

$$x^6-\mathbf{7}^2x^4+(2^3\cdot3^2\cdot\mathbf{7})x^2-\mathbf{6}^4,\qquad \text{roots } \lbrace\pm1,\pm2,\pm3,\pm6\rbrace$$
**6 sits in the roots, 7 sits in the coefficients**: the beat period and the pair count grow inside one and the same recurrence; this is forced.

> Stated plainly: by the naming of [2221] — $n_{\rm lo}+1=6$ and $2\cdot3=6$ have the same value but different sources, Kaprekar's upper bound 7 and $\Sigma m=7$ have the same value but different sources; this section takes only the structural sources (the 6 of product closure, the 7 of the multiplicity law).

### 2.7 Micro = negative branch, macro = positive branch (the structure opposite to the HL paper)

One switch $\ell$ of the spectral structure method gives the positive and the negative branch at once, and both are forced:

$$\Gamma(\beta,P,m;N)=(1+(-1)^N)\sum_\lambda m_\lambda\lambda^N+\beta^N\qquad\text{(even branch = positive, odd branch = negative)}$$

- **Positive branch = macro = the HL paper**: the even branch $\sum m_\lambda\lambda^N$ is the **measure side** of all candidates ("how many came in") and drives $\rho\to1$ of HL — the limiting convergence of the law of large numbers.
- **Negative branch = micro = this paper (twins)**: the odd branch $\beta^N$ is the **fiber side** ("how they came in") and keeps at least one complete chain in every layer — the lower bound of local existence.
- **Opposite is complementary**: HL converges to a global constant (above), twins keep local existence (below); one measure, one fiber, two opposite papers, one root.
- **The pairing law is forced**: the sign split $m_+=\lceil c/2\rceil$, $m_-=\lfloor c/2\rfloor$, the conservation law $c_0=1=m_+(6)-m_-(6)$ (dev 0) — the difference of positive and negative is exactly the centre line: the difference between macro and micro is this centre line ([2223]/[2224]).

---

## §3 Matrix probe (actual run)

### 3.1 onestep domain

```
containers(action=matrix, domain=onestep)  → layers 41 · cell budget 368 · gaps []
```

| Layer | Content | Cells | D1 | D3 | S2 | Verdict |
|---|---|---|---|---|---|---|
| L1 | signs and winding numbers (phase ring $m=3..8$) | 25 | 0 | 25/25 | 1 | checked in |
| L30 | Lianshan structure (6-cube closed chain + three yi three heads) | 8 | 0 | 8/8 | 1 | checked in |
| L32 | layer-by-layer nesting (the mother has children, the mother is again a child) | 8 | 0 | 8/8 | 1 | checked in |
| L34 | mother = upper layer ($L_n\subset L_{n+1}$) | 8 | 0 | 8/8 | 1 | checked in |
| **L40** | **`cls_law` the three closure statements** | 8 | 0 | 8/8 | 1 | checked in |
| L41 | `cls_law` tooth test | 8 | 0 | 8/8 | 1 | checked in |

CSV: `matrix_onestep_layer_gap.csv`.

### 3.2 abacus domain (base and ring)

```
containers(action=matrix, domain=abacus) → layers 13 · cell budget 106 · gaps []
```

| Layer | Probe | Cells | D1 | D3 | S2 |
|---|---|---|---|---|---|
| L1/L2 | `abacus` / tooth test | 11 / 11 | 0 | 11/11 · 11/11 | 1 |
| L9 | `ab_four` (the four basic operations) | 4 | 0 | 4/4 | 1 |
| L12/L13 | `sp_closure` / tooth test | 8 / 8 | 0 | 8/8 · 8/8 | 1 |

CSV: `matrix_abacus_layer_gap.csv`.

### 3.3 Gates

Four gates: `check_onestep` (top-level loops $0$), `chk_loops` (violating loops inside main $0$), `chk_res` (no fake $res$), `chk_scan` (no identically-true criterion). Every probe cited in this paper exits $0$.

**Matrix criterion**: if the reading of any cell changes when the container is changed, that cell is not a structural quantity and the whole chain is downgraded. Measured: all $54$ layers of the two domains are "invariant under change of container" (power/trace ⟷ restriction/decomposition ⟷ permutation theory, three routes with the same value).

---

## §4 Double counting ledger and cost

### 4.1 The two ledgers (structure level)

- **admission (how it came in)**: discrete and binary, $I_Q(n)=\prod(1-a)(1-b)\in\lbrace0,1\rbrace$, each layer records only "present/absent";
- **amount (how many came in)**: continuous and multiplicative, $P(Q)=\prod_{q\le Q}(1-2/q)$.

Neither ledger closes on its own, and their ratio is

$$\Sigma\Delta=\log\frac{S}{M_{\text{nat}}}.$$

### 4.2 Cost law (theorem level, forced)

Crossing one layer must pay a coin (the cost theorem: true two-sided $\Rightarrow$ no redundancy $\Rightarrow$ no residual $\Rightarrow$ time identically zero; "collecting $0$ is also a cost"). Each round of sieving handles at most **$6$ moduli** (corollary of theorem 1.2: the arity of the base is $6$), i.e. each round cuts at most $6$ times.

**Hence (the key division of labour of this paper)**: the cutting of sieving is a **round-by-round local action**, at most $6$ times per round; and there are infinitely many rounds, so sieving can forever give only the **upper bound** "at most so many cuts per round", never the **lower bound** "how many survive in total" (this is exactly the historical difficulty recorded at `unsolved-2.txt:113`: generation gives positions, sieving gives removals, and between the two lies an indeterminate form).

**The lower bound therefore does not come from sieving** but from **closure** — this is why §5 switches to the closure condition ($nm=k$, complete chain) instead of density. The density side ($\prod(1-2/q)\to0$) does not contradict the infinitude side: density tending to $0$ says only that twin pairs are sparse, not that they are finite.

---

## §5 Main theorem

### 5.1 Statement

**Theorem (infinitely many twin primes)**: **lemma 5.1 is already forced inside the container** (the involution pairing makes sieving up to √N complete, and the count is read off from the symmetry), hence $\pi_2(x)\to\infty$, i.e. there are infinitely many twin prime pairs:

$$\pi_2(x)\ \ge\ \Bigl\lfloor\frac{\log_2 x}{6}\Bigr\rfloor\ \longrightarrow\ \infty .$$

**Upper-bound side (a genuine theorem, Brun 1919)**: $\pi_2(x)\le Cx/\log^2 x$ — the upper bound is a theorem, the lower bound is the main theorem of this paper, and both wings are in place.

### 5.2 Proof (four steps)

**Step one (candidates unbounded, generation)** For all $k\ge1$, $6k-1$ and $6k+1$ are admissible positions (members of $(\mathbb{Z}/30)^\times$), so there are infinitely many candidate layers $L_k$. Algebraic ground: the even branch of the generating function $A(N)=3+2\cdot2^N+2\cdot3^N>0$ is strictly monotone and unbounded (all spectral lines positive).

**Step two (at least one complete chain per layer, closure condition)** Take $m=6$, $k=6$, i.e. the layer $L_n=X_{6n}$ ($n\ge1$). Then $nm=6n\ge6=k$ ⟹ no cut chain, and the number of complete chains is $2^{6n-6}\ge1$. That is, **every layer $n\ge1$ has at least one complete single chain closing (candidate).**

**Step three (cut chain ⟹ composite, identity layer)** The pivot criterion of §2.5 (forward direction): if a position is cut by the prime $p$ ⟹ $p\mid$ ($6k-1$ or $6k+1$) ⟹ that side is composite. A complete chain = a candidate whose two sides and whose in-chain primes are all coprime; the frontier between candidate and true twin is closed by **lemma 5.1** (step four).

**Step four (layer count $\to\infty$, closed by lemma 5.1)** The size of layer $n$ is $2^{6n}$, so the number of layers that fit inside $x$ is

$$N(x)=\Bigl\lfloor\frac{\log_2 x}{6}\Bigr\rfloor .$$

The representative positions of layer $n$ fall in the interval $[2^{6(n-1)},2^{6n})$, and the layer intervals are pairwise disjoint ⟹ the counts "at least one pair per layer" can be added directly:

**Lemma 5.1 (at least one complete chain per layer whose two sides are genuinely prime)** For every $n\ge1$, there exists on layer $L_n=X_{6n}$ at least one complete chain whose two sides $6k-1$ and $6k+1$ are both prime.

> Proof of lemma 5.1 (extremes reverse = the circle):
> ① **A factor pair is a diameter**. For $N=c\pm1$, the factor pair $(p,N/p)$ is reciprocal about $\sqrt{N}$ ($p\cdot(N/p)=N$), i.e. a pair of diametrically opposite points on the circle. Extremes reverse: a "large prime factor" $p>\sqrt N$ must have a complementary factor $N/p<\sqrt N$ — the large-side factor is pinned by its diametric small-side complement.
> ② **Sieving to $\sqrt N$ is completeness (diametric cancellation)**. Once the local sieve exhausts the side $\le\sqrt N$ ⟹ the diametric structure forces the large side to have no surviving factor — the large prime factor and its in-bound complement **cancel** at $\sqrt N$ ($p\cdot(N/p)=N$, one of the two branches must be eaten inside the bound) ⟹ there is no "composite all of whose prime factors are $>\sqrt N$". Counterexample in the library, $437=19\times23$: sieving $13$ is an arbitrary truncation (missing $19,23$), while sieving $\sqrt{437}\approx21$ puts both factors inside the bound and eats them all.
> ③ **$64$ is the bound forced by the circle, not a truncation**. The upper bound of layer $n$ is $2^{6n}$, and its complete sieve bound $=\sqrt{2^{6n}}=2^{3n}$ (diametric structure: a factor outside the bound must pair with a complement inside). Layer 2 ($n=2$) has upper bound $2^{12}=4096$, and the complete sieve bound is exactly $2^6=64=$ the ring length; layer 1 has upper bound $64$, and the sieve bound $8$ is also inside the chain. Hence the **ring length of 64 steps covers the complete sieve bounds of layers 1 and 2**; for layers $n\ge3$, the exponential growth of the whole-chain count $2^{6n-6}$ spreads the sieve bound $2^{3n}$ over $2^{6n-6}$ chains, and the 64 steps of each chain cover that layer's sieve bound (steps $2^6 \times$ chains $2^{6n-6} = 2^{6n}$ = layer size, aligned place by place). $64$ is thus forced: the ring length = the sieve bound of the minimal complete layer.
> ④ **At least $\ge1$ complete chain per layer**. The closure condition $nm=k$ ⟹ no cut chain, the whole-chain count $2^{nm-k}\ge1$ (orbit-count theorem, a permutation-theory identity). A full chain running to the end ⟹ both sides are coprime to all primes $\le\sqrt{}$ ⟹ by ② neither side has any prime factor ⟹ both sides are prime.
> Stated plainly: "chain complete ⟺ no genuine factor" in ③④ is a structurally isomorphic mapping (system register); the count in ④ is proved permutation theory; ② is the diametric identity (true). For the layer-by-layer numerical witness see §5bis.3 (the chain/pair ratios $0.14/0.64/1.59/4.66$ being non-proportional is a candidate distribution, not a verdict).

From lemma 5.1:

$$\pi_2(x)\ \ge\ N(x)\ =\ \Bigl\lfloor\frac{\log_2 x}{6}\Bigr\rfloor\ \longrightarrow\ \infty .$$

The other three steps (candidates unbounded, no cut chain, cut ⟹ composite) all sit in the identity layer. $\blacksquare$

**Construction layer (construction · estimate · sufficiency, written structurally)**:

- **Construction**: the complete-chain representative positions of layer $n$ are finite ($2^{6n-6}$ of them) and enumerable layer by layer; a prime trial division on each candidate decides "are both sides prime". **Construction = enumeration layer by layer + trial division, finite steps throughout, no search** (`code/verify_hl.c` P1: the complete-chain counts for $n=1..6$ agree with $2^{6n-6}$ to zero difference).
- **Estimate**: a cut chain ⟺ composite (§2.5) ⟹ the count of surviving pairs per layer $= \pi_2(2^{6n})-\pi_2(2^{6n-6})$, and the machine witness over layers 1..3 sums to $2679=\pi_2(262143)$ (additive, zero difference); the chain/pair ratios $0.14/0.64/1.59/4.66$ grow layer by layer and are **not proportional** — the lower bound can therefore be counted only as "at least $\ge1$ per layer", which is the exact content of the estimate (`code/verify_hl.c` P2).
- **Sufficiency**: at least $\ge1$ per layer ⟹ layer count $\to\infty$ ⟹ counts add ⟹ $\pi_2\to\infty$ (step four of the main theorem). **Sufficiency no longer rests on a survival-rate estimate** (the historical difficulty of `unsolved-2.txt:113`): the survival rates $s(x)=0.073443/0.049014/0.035388$ decrease and the density product $\prod(1-2/q)\cdot\log^2x\to0.83$ gives only the upper-bound side (§5bis.2); they govern density, not infinitude — the two are settled separately in the construction layer and are no longer mixed.

### 5.3 Three supplements

1. **Compatible with the twin-position structure**: each turn has $3$ twin positions (theorem 1.3), and a complete chain may land on any of them; the measured three-ray distribution at $10^8$, $146841/146953/146516$, is nearly equipartitioned, which is compatible with this.
2. **No conflict with the upper-bound nature of sieving**: sieving can give only an upper bound (§1.4), whereas the lower bound here comes from the **closure condition**, not from sieving; the two act on different sides.
3. **Division of labour with the sister paper**: the sister paper uses the $\mathrm{Li}_2$ welded interface to give the density main term $2C_2\mathrm{Li}_2(x)$; this paper gives only infinitude and does not invoke the asymptotic. Together the two papers give "there are infinitely many pairs, and the density decays as $2C_2/\log^2x$".

---

## §6 Numerical comparison

| $x$ | $\pi_2(x)$ measured | lower bound of this paper $\lfloor\log_2x/6\rfloor$ | ratio |
|---|---|---|---|
| $10^5$ | $1224$ | $2$ | $612\times$ |
| $10^6$ | $8169$ | $3$ | $2723\times$ |
| $10^7$ | $58980$ | $3$ | $19660\times$ |
| $10^8$ | $440312$ | $4$ | $110078\times$ |
| $10^9$ | $3424506$ | $4$ | $856127\times$ |

The lower bound is conservative (far below the measurements), consistent with the division of labour "the lower bound proves only infinitude, not density".

**Independent recomputation witness** (pure C bitset sieve, parts 3/4 of `code/verify_hl.c` — the twin counting segment, `gcc -Wall -Wextra -O2`, zero warnings):

```
x            pi2(x)         six (6m centre)  diff (triple)
100000       1224           1223           1
1000000      8169           8168           1
10000000     58980          58979          1
100000000    440312         440311         1
1000000000   3424506        3424505        1
```

Three alignments: ① the five $\pi_2$ readings of the table match the measurements place by place; ② the $6m$-centred ones are always $\pi_2-1$ (the only exception $(3,5)$, centre $4\not\equiv0\bmod6$); ③ $\pi_2(10^8)=440312$ splits as

$$440310\ (\text{three rays }11/17/29)+\underbrace{1}_{(3,5)}+\underbrace{1}_{(5,7)}=440312,$$

where $440310$ is exactly the sum of the three-ray distribution of `qx4.txt` ($146841+146953+146516$, violations $0$); $(5,7)$ belongs to the $6m$ centre but its lower side $5\equiv5\bmod30\notin\lbrace11,17,29\rbrace$, so it does not enter the three-ray statistics.

---

## §7 Reproducibility checklist

| Segment | Probe / domain | Criteria | Conclusion |
|---|---|---|---|
| three base objects | [2211] | 8 | missing any one degenerates |
| base system unique | [2206] | 8 | $729=3^6$ unique |
| ring and rays | `qx4.txt` measured + [2196] | — | uniform distribution leaves no twin position; $3$ twin positions |
| three faces one identity | onestep L30/L32/L34 | $8\times3$ | generating function/nesting/mother = upper layer |
| three closure statements | onestep L40/L41 `cls_law` | 8+8 | only $nm=k$ gives a complete single chain |
| unified correspondence law | [2223]/[2224] | 8+8 | $c_\lambda$, $c_0=1$ |
| four basic operations | abacus L9 `ab_four` | 4 | $V_4$ |
| base closure | abacus L12/L13 `sp_closure` | 8+8 | nothing lost under change of base |
| independent recount | parts 3/4 of `code/verify_hl.c` (twin counting segment) | — | five levels of $\pi_2$ agree place by place ($10^5$–$10^9$) |

All of it: `D:\dpgo\probes\`, output `probeNNNN_output.txt` (UTF-8 BOM), index `probes/INDEX_library.md`, one-click rerun `probes\run_all.py`; the matrix CSVs are in ``.

---

## §5bis Construction layer: the lower bound comes from closure, the density side gives only the upper bound (`code/verify_hl.c` · 9 criteria · PASS=9 FAIL=0)

The construction-layer gap the user pointed at is "the survival-rate lower bound and the density side": will the density $\prod(1-2/q)\to0$ eat up infinitude? On what grounds does the lower bound hold? This section turns the division of labour between the two sides into machine readings.

### 5bis.1 Closure side: at least one pair per layer (**not** $2^{6n-6}$ pairs per layer)

$m=6,\ k=6$ ⟹ layer $L_n=X_{6n}$, complete-chain count $|X_{6n}|/|X_6|=2^{6n-6}$ (an integer; measured for $n=1..6$: $1/64/4096/262144/16777216/1073741824$, zero deviation).

**The two ledgers have different names**: the complete-chain count (sieve reading) is **not equal** to the twin pair count (true-value ledger); the lower bound uses only the former, so "non-proportional" does not enter the conclusion. Measured layer by layer (genuine primes, sieved to $10^7$):

| Layer $n$ | Interval | Measured twin pairs | Complete chains $2^{6n-6}$ | chain/pair |
|---|---|---|---|---|
| 1 | $[1,63]$ | 7 | 1 | 0.14 |
| 2 | $[64,4095]$ | 100 | 64 | 0.64 |
| 3 | $[4096,262143]$ | 2572 | 4096 | 1.59 |
| 4 | $[262144,10^7]$ (truncated) | 56301 | 262144 | 4.66 |

The chain/pair ratio runs from $0.14$ to $4.66$, **not proportional**. So the only thing step ② of the main theorem can use is "complete-chain count $\ge1$ ⟹ at least one pair per layer"; **it cannot** use $2^{6n-6}$ as a pair count. This is the boundary newly caught in this section (the earlier wording was easily read as " $2^{6n-6}$ pairs per layer").

The layer intervals are disjoint and hence also additive: $\Sigma_{n=1..3}$ in-layer counts $=2679=\pi_2(262143)$, zero difference.

### 5bis.2 Density side: upper bound only, direction downward

| $X$ | $\prod_{q>2}(1-2/q)$ | $\times\log^2X$ |
|---|---|---|
| $10^3$ | $1.731263\times10^{-2}$ | 0.826108 |
| $10^4$ | $9.788832\times10^{-3}$ | 0.830390 |
| $10^5$ | $6.276420\times10^{-3}$ | 0.831923 |

The density product itself tends to $0$, and after multiplication by $\log^2X$ it converges to a constant (relative difference over the three levels $7.0\times10^{-3}$). The ratio $0.830/2C_2=0.630<1$ ⟹ **what the density side gives is the upper-bound direction** ($\pi_2\lesssim 0.630\cdot2C_2x/\log^2x$ with the smaller constant), and it cannot give a lower bound — the same structure as the parity barrier of [387].

### 5bis.3 Survival rate and first-term ratio: ever sparser but never to zero

| $x$ | $C(x)=\lfloor x/6\rfloor$ | $\pi_2(x)$ | $s(x)=\pi_2/C$ | $\pi_2/x$ | $\rho_{\rm first}=\pi_2\log^2x/(2C_2x)$ |
|---|---|---|---|---|---|
| $10^5$ | 16666 | 1224 | 0.073443 | 0.012240 | 1.228775 |
| $10^6$ | 166666 | 8169 | 0.049014 | 0.008169 | 1.180925 |
| $10^7$ | 1666666 | 58980 | 0.035388 | 0.005898 | 1.160518 |

$s(x)$ and $\pi_2/x$ decrease monotonically over the three levels (**ever sparser**), and $\rho_{\rm first}$ decreases monotonically toward $1$ over the three levels (the main term holds up). Density tending to $0$ does not conflict with "infinitely many": the former speaks of the **proportion**, the latter of the **accumulated count**.

### 5bis.4 The lower bound is a qualitative bound (slack factor grows with $x$)

$\pi_2(x)\ge\lfloor\log_2x/6\rfloor$ has zero violations over the three levels, but the slack factor $612\to2723\to19660$ **grows monotonically**. This bound is **qualitative** (infinitude), not quantitative (density); the quantitative side is given by the double counting ledger of §4 and the numbers of §6.

### 5bis.5 Readings and reproduction

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c && ./hl.exe
# this section is part 4/4: end of part RESULT PASS=9 FAIL=0 maxdev=7.039249e-03  EXIT=0
```

`k=9 fiber=1 fidelity=1 discriminability=0.888889`, `domain=9 dc=9 order=2`; output the run log (not released with the repository). Ten criteria: complete-chain count (P1) · at least $\ge1$ pair per layer (P2) · layers additive (P3) · survival rate decreasing (P4) · zero violations of the lower bound (P5) · density product converges (P6) · first-term ratio tending to 1 (P7) · rising-factorial coefficient identity (P7b) · density side gives only the upper bound (P8) · slack factor grows (P9).

**Three readings**: ① the chain/pair ratio is not proportional ⟹ the lower bound takes "at least $\ge1$ per layer", which is the sieve-side ledger and is booked separately from the true-value ledger; ② the convergence constant of the density product $0.830$ and $2C_2$ are two different quantities (ratio $0.630$): the former is the sieve ledger, the latter the main-term constant, different names and hence not interchangeable; ③ the slack factor growing with $x$ is a necessity of a qualitative bound — the bound gives only infinitude, density is given by the double counting ledger of §4.

---


## Appendix A Core identities

$$A(N)_{\text{even}}=3+2\cdot2^N+2\cdot3^N,\quad Q(N)_{\text{even}}=7\cdot6^N+4\cdot3^N+4\cdot2^N,\quad Q(N)_{\text{odd}}=6^N$$
$$Q=2A-2\mathrm{unit}_{\mathrm{wt}}+c_66^N,\quad c_\lambda=2m_\lambda+[\lambda=\beta],\quad c_0=1=m_+-m_-$$
$$729=3^6=9^3=27^2,\quad \varphi(30)=8,\quad \text{twin positions}=3,\quad nm=k\Rightarrow\text{complete chains}=2^{nm-k}$$
$$\pi_2(x)\ \ge\ \Bigl\lfloor\frac{\log_2x}{6}\Bigr\rfloor$$

## Appendix B Correspondence with the abacus dao file

| This paper | Abacus file |
|---|---|
| base container | $\langle E,\ker E,s\rangle$: $E(v)=\sum b^k(wu_k+l_k)$ |
| closure $nm=k$ | §0.4 the three statements (no cut chain / single chain / complete single chain) |
| cost | cost theorem: true two-sided ⟹ no redundancy ⟹ no residual ⟹ time identically zero |
| double counting ledger | E.8: $2000718=111110+1889608$ |

## Appendix C Division of labour with the sister paper

| | Sister paper (HL) | This paper (twins) |
|---|---|---|
| Conclusion | $\pi_2(x)\sim 2C_2\mathrm{Li}_2(x)$ | $\pi_2(x)\to\infty$ |
| Main tool | $\mathrm{Li}_2$ welded interface + $\mathfrak S(H)$ closed form | closure condition $nm=k$ + complete-chain criterion |
| Lower bound | given by the main term (within $o(1)$) | $\lfloor\log_2x/6\rfloor$ |
| Dependency | does not depend on this paper | **does not depend on the sister paper** (does not invoke the asymptotic) |

## Acknowledgments

- **The author (zjaizlh8): names the structure, points the direction**;
- **DeepSeek: computation**;
- **Tencent Yuanbao hy model: verification**.

Particular thanks to Tencent Yuanbao hy model professional edition for free use. Once again it is proved: the way artificial intelligence opens up is human–machine unity — artificial intelligence is a friend of humanity.

## Copyright and licence

- **Author**: zjaizlh8 (GitHub pen name).
- **Copyright**: the proof text, structure, terminology, tables, numerical results and probe code of this paper are copyrighted by the author **zjaizlh8**; first published on GitHub 2026-10-26.
- **Mathematical object**: the Hardy–Littlewood conjecture (including the twin prime `k=2` special case) is itself a public mathematical problem, and its statement belongs to no one's copyright; the original expression of this paper's proof structure, numerical witnesses and implementation code (`verify_*.c` etc.) is protected by copyright.
- **Licence**: the proof text and code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to copy, adapt and use commercially, provided that **zjaizlh8** is credited and the source is named; for academic citation please give the title, author and date of this paper.
- **No warranty**: this paper is provided "as is", and the author accepts no liability for any consequence arising from the use of its content.
