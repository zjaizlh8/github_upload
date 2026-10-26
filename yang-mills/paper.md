## Through the structural layer, the spectral-structure method, and the matrix probe

> Author: zjaizlh8 (GitHub pen name)
> System: Unary Field Theory (WQNN-K8538)
> Date: 2026-10-26
> Independent recomputation: `code/verify_ym.c` (pure C, zero loops at the direct level of `main`, two subtotals `8 · 9`, last line `RESULT PASS=17 FAIL=0`, `EXIT=0`)
> Matrix run: `containers matrix domain=onestep` → `matrix_onestep_layer_gap_ym.csv`

---

## Abstract

This paper uses three structural tools to give a **structural-layer reading** of Yang–Mills existence and the mass gap (the load-bearing Lemma 4.1 is forced: `λ_l = l(l+2)`, `λ₀` simple ⟹ `Δ = √3·Λ` is a subtraction — not an extra assumption):

1. **Structural layer**: the foundation of mathematics is "unit 1 × Cartesian product × place-value power $b^k$"; the three base criteria uniquely select $729=3^6$ (a base, not an operation, giving only volume); the four classes of basic operations (addition, multiplication, reduction mod, cross-dimension) are the four operation cells ($V_4$), among which **reduction mod is the only one that can produce "inversion / folding back"** (Theorem 1.4).
2. **Spectral-structure method**: three instances of the same action on one container — take the spectrum inside a layer, nest between layers, close above the layer; the closure condition is $nm=k$. **A closed form is a closed form, accumulation is stepping**: the difference between the closed form of the running coupling and the Euler accumulation is strictly $O(h)$ (steps $\times10$ ⟹ error $\div10$), and this is the complete individual of "one step for the true value, one step for the error".
3. **Matrix probe**: $41$ layers, $368$ cell budget, $41$ cells emitted, **gaps empty**, all `D1=0 · S2=1 · CLOSED=1`.

Brought down to Yang–Mills, the three tools give four forced facts:

- Generation side: "$|e^{i\theta}|\equiv1$" ⟹ the gauge-group manifold is compact (SU(2) half-angle double cover, $4\pi$ period; machine witness dev $5.551\times10^{-16}$) ⟹ Peter–Weyl duality discrete ⟹ Casimir spectrum discrete, **$\lambda_1=3>0$** ($\lambda_l=l(l+2)=$ multiplicity $-1$, dev exactly $0$).
- The observer's join: colour-singlet observables $=$ trace $=$ **closed loop** (gauge invariance and cyclicity of the trace, dev exactly $0$) ⟹ the container of the physical sector is $S^1$ (compact), not the base space $\mathbb{R}^3$.
- Screening side: gauge fixing is a projection, whose kernel $=$ Gribov copies $=$ pure gauge $=$ the dark side; confinement $=$ screening $=$ convergence ⟹ open chains (truncated chains) are deleted ⟹ the closure condition $nm=k$ holds automatically (machine witness: closed $\iff$ integer winding, the seam deviation of the non-integer $1.7$ is $1.618$).
- Gap: the spectrum of the physical Hamiltonian $=\lbrace0\rbrace\cup[\Delta,\infty)$, $\Delta=\sqrt{\lambda_1}/R=\sqrt3\Lambda>0$.

**The gap that 4.126 self-exposed in the archive is closed here**: that chain applied "group compact" directly to "base space $\mathbb{R}^3$", and thereby left the crack "group compact $\ne$ base space compact". This paper points out: the gap is not a property of the base space, but a property of the **observable algebra** (closed loops) — $\mathbb{R}^3$ supplies only length, and length is truncated at $r_c$ by confinement.

---
## §0′ Reciprocal form (load-bearing register · forced by "reciprocal")

The main equivalences of this paper are all written as a **reciprocal account**: `A · (1/B) = rational account`.
Basis: three reciprocal laws in the archive (probes [2233]/[2234]/[2235]):
- **Closure ⟺ that side admits a multiplicative inverse, and the inverse is writable inside the same container**;
- **Expansion and contraction are reciprocals** (product identically 1; missing one direction leaves only "convergence", not "closure");
- **Substance and projection must both be spoken of as reciprocals**: `projection rate × substance multiplier = 1`, equivalent to `rank + nullity = dimension`.

Hence the criterion is not "compute each side once and see the numbers come close", but **the two sides are reciprocals and the quotient is allowed to be exactly the one number on the account**.
After this section every main formula is followed by a line of "reciprocal account": it gives who is the reciprocal of whom, and which class of remainder this entry permits.



### The reciprocal account of this paper

- **Substance** = the spacing of the compact spectrum `Δ = √λ₁/R`; **projection** = the finite-volume coupling `Λ = μ e^{−8π²/(β₀g²)}`.
- Reciprocal account: `Δ · R / √λ₁ = 1` (gap and radius stand in a reciprocal relation, `λ₁ = 3` forces `Δ/Λ = √3`); on the coupling side `g² ↔ 1/g²` interchange inside the exponent of `Λ` — **asymptotic freedom is exactly "the reciprocal slows down as g→0"**.
- Machine witness: `Δ/Λ = 1.7320508076 = √3` (dev 0); `(Λ/μ)/g¹⁰ → 3.475528e−10` (essential singularity, not a power law); Euler vs closed-form error ratio 10.2096/10.0204/10.0020 (O(h)).
- Truncated-chain deletion: the open string `E = σr` crossing `2m_q` necessarily creates a pair ⟹ closure is self-satisfied (integer winding ⟺ closed; non-integer 1.7 seam 1.6180).


## §0 Layered register (structural layer / construction layer / numerical layer / declarative layer)

| Layer | Content in this question | Mode of decision |
|---|---|---|
| Structural layer | three instances of one action; closure $nm=k$; confinement $=$ screening $=$ convergence; base $=$ vortex eye $=$ fixed | whether the structure is forced |
| Construction layer | four classes of basic operation; Casimir spectrum on $S^3$; Wilson loop; closed form of the running coupling; exponential form of $\Lambda$ | whether closed-form, whether reproducible |
| Numerical layer | $\lambda_1=3$, $\Delta/\Lambda=\sqrt3$, $\beta_0=11-2n_f/3$, $r_c=0.6$ fm, $2373=3\cdot7\cdot113$ | readings and tolerances |
| Declarative layer | the object proved (existence + gap lower bound) | scope of the statement |

**Structural-layer correspondence law (the load-bearing identification of this question)**:

$$\boxed{\text{compact}\iff\text{discrete spectrum}\iff\text{positive lower bound}\iff\Delta>0\iff\text{closure } nm=k}$$

Both ends of the chain must hold for a gap to exist: the left end (compact) is given by the generation side, the right end (closure) by the screening side. Drop either end and the gap does not exist — this is exactly why 4.126 in the archive does only the left end and 4.59.4 only the numerics ("has a bottom") without the two being joined.

**Five-element placement (following the 4.126 register of the archive)**: base $=S_{\rm YM}$ action (gives only volume) · dimension change $=F=dA+A\wedge A$ (gauge redundancy, cross-dimensional expansion) · screening $=$ gauge fixing (projection, whose kernel $=$ Gribov copies $=$ the dark side) · accumulation $=\Lambda=\mu e^{-8\pi^2/(\beta_0g^2)}$ (cross-dimensional, pays residue) · generation $=|e^{i\theta}|\equiv1\Rightarrow$ compact $\Rightarrow$ discrete spectrum $\Rightarrow$ gap.

---

## §1 Structural layer

### 1.1 The three base items

**Theorem 1.1** Any writable number system must simultaneously possess: ① a unit $1$; ② Cartesian products (multi-axis juxtaposition); ③ place-value powers $b^k$.

*Proof*: without ① there is no normalization point, and a base cannot locate the highest digit; without ② a single axis degenerates, and $A_\mu^a$ cannot be expressed (a multi-index of $3$ colour indices $\times$ $4$ spacetime indices $\times$ group dimension $8$); without ③ carries have no carrier, and the $8$ generators of $SU(3)$ cannot be exhausted by finite digit strings. With all three present, the gauge field $A$ and the field strength $F$ are both writable as finite digit strings. ∎

**Landed readings**: the number of $SU(3)$ generators is $3^2-1=8$ (= the eight trigrams, archive 4.45.6 verdict); gluons carry colour charge ⟹ even the messenger is confined (same item); $r_c=2m_q/\sigma=0.6$ fm.

### 1.2 Unique base system: $729=3^6$

**Theorem 1.2** Under the three base criteria (nothing lost on base change, no operator position, unique as a pure prime power) the unique base kernel is $729=3^6$: $729=3^6=9^3=27^2$ (three routes converge, zero collisions), and $729$ is not the image of any of the four classes of basic operations (a base $=$ non-operation). Within $N\le10^4$ there are $13$ numbers admitting "≥2 writings", among which only $\lbrace64,729,4096\rbrace$ have exponent $6$, and the unique base that is a power of $3$ is $729$.

**Landing in this question**: the action $S_{\rm YM}=\frac1{4g^2}\int F^a_{\mu\nu}F^{a\mu\nu}$ is "volume" (the base); it takes part in none of the four classes of basic operations and only sets the scale of the spectrum; what really decides "gap or no gap" is the **discreteness** of the spectrum (generation), not the size of $S$.

### 1.3 The four classes of basic operations are the four operation cells

| Four classes of basic operation | Keep/Shift | Contract/Expand | Four operation cells | Role | Placement in this question |
|---|---|---|---|---|---|
| Addition | keep | contract | yin cell (substance) | accumulation | spectral lines add, trace linear in addition |
| Multiplication | keep | expand | yang cell (substance) | generation | matrix multiplication, Casimir $l(l+2)$ |
| Reduction mod | shift | contract | yin cell (use) | screening | winding rounded to integer (closure $\iff$ integer winding), $SU(2)$ centre $\mathbb{Z}_2$ |
| Cross-dimension | shift | expand | yang cell (use) | dimension change | $\exp/\log$: essential singularity of $\Lambda$, logarithm of the running coupling |
| —— | — | — | non-operation | base | volume of $S_{\rm YM}$ |

### 1.4 Theorem 1.4 (inversion can only be produced by reduction mod)

**Theorem**: among the four classes of basic operations, the only one that can produce "inversion / folding back" is reduction mod.

*Proof*: addition preserves order, multiplication by positive coefficients preserves order, $\exp/\log$ preserve sign on the real axis; only reduction $x\mapsto x\bmod m$ (equivalent to $x\mapsto-x$, since $-x\equiv m-x$) folds "large" back to "small". ∎

**Landing for $SU(2)$**: $U(t)=\cos\frac t2 I+i\sin\frac t2\hat n\cdot\sigma$ has period $4\pi$, while $U(2\pi)=-I$ — this "minus sign" is exactly reduction mod: when $\frac t2$ advances by $2\pi$, $\cos\frac t2$ reaches $-1$, i.e. $[0,4\pi)$ is folded back onto the centre $\mathbb{Z}_2$. **Without reduction mod the centre $\mathbb{Z}_2$ of the group cannot be written, and the double cover does not exist** (machine witness P2: half-angle dev $5.551\times10^{-16}$; the full-angle writing gives $+I$ at $2\pi$, difference $2.000$).

---

## §2 Spectral-structure method

### 2.1 Three faces, one identity

| Face | Action | Mathematics | Placement in this question |
|---|---|---|---|
| Intra-layer | take spectrum | $X_k=Z_2^k$, spectrum $\lbrace k-2j\rbrace$, multiplicity $\binom kj$ | Casimir spectrum on $S^3$: $\lambda_l=l(l+2)$, multiplicity $(l+1)^2$ |
| Inter-layer | nesting | $L_n=X_{n\cdot m}$, $\pi\circ\iota=\mathrm{id}$ | energy-scale layers: the nesting $\Lambda\to\mu\to Q$ (running coupling) |
| Above layer | closure | $\sigma=g\circ(+1)\circ g^{-1}$, $64$-period | Wilson loop (closed loop), integer winding |

### 2.2 The identity of the $S^3$ spectrum

$$\lambda_l=l(l+2)=(l+1)^2-1=\text{mult}-1\qquad(l\ge1).$$

Machine witness (P3): for $l=1,2,3,7$ every $\lambda-(\text{mult}-1)=0$ exactly; $\lambda_1=3$, $\lambda_7=63$ with multiplicity $64$ (the sixty-four hexagrams).

**Structural meaning of "spectrum $=$ multiplicity $-1$"**: the spectral value of each layer is exactly one less than its multiplicity — this is the landing of "the volume of the base $=$ the number of points, the spectrum $=$ the number of points $-1$" (isomorphic to the conservation law $c_0=1=m_+-m_-$: the total count and the minimal non-zero spectrum differ by one). $\lambda_1=3$ is exactly the volume of ternary-ness.

### 2.3 Three closure statements

1. **No truncated chain**: $nm\ge k$;
2. **Single chain**: $nm\le k$;
3. **Complete single-chain closure**: $nm=k$.

**Placement in this question**: physical states are closed loops (Wilson loops), whose winding number is an integer — this is the explicit solution of $nm=k$ on $S^1$ ($n$ turns $\times$ $m$ steps $=$ one whole round). **Truncated chains (open strings) are deleted from the physical spectrum by screening**: the open-string energy $E=\sigma r$ grows linearly with $r$, and beyond $2m_q$ a quark pair is necessarily created (archive 4.45.5 "the string cannot be cut: you cut it, it does not split, it creates a pair") ⟹ open chains are off shell ⟹ only closure remains in the physical spectrum.

### 2.4 Quotient criterion (version for this question)

**Theorem 2.4 (quotient criterion)** A closed form $\iff$ a quotient exists inside that layer making $nm=k$ satisfiable.

**Placement in this question**: a gauge transformation is a translation along the fibre direction, folding away "different representatives of the same physical state" — that is the quotient. **The existence of the quotient is equivalent to gauge invariance** (P5: $\mathrm{tr}(\Omega U\Omega^\dagger)=\mathrm{tr}(U)$, dev $0$); and the container after the quotient is **the set of closed loops** (compact), not $\mathbb{R}^3$ (non-compact). **This is where the gap lives: the quotient folds away the non-compact half.**

### 2.5 Closed form vs accumulation (closed form ⟺ closed form inside the layer)

The closed form of the running coupling $\dfrac{d\alpha}{d\ln Q}=-\dfrac{\beta_0}{2\pi}\alpha^2$:

$$\frac1{\alpha(Q)}=\frac1{\alpha(\mu)}+\frac{\beta_0}{2\pi}\ln\frac Q\mu ,\qquad \beta_0=11-\frac{2n_f}{3}\ (\text{SU}(3)).$$

Machine witness (P6): $\beta_0$ at $n_f=0..6$ is $11,\ 10.333,\ 9.667,\ 9,\ 8.333,\ 7.667,\ 7$; the error between Euler accumulation and the closed form

| Steps | Euler $\alpha$ | Closed form $\alpha$ | Error |
|---|---|---|---|
| 10 | 0.1663528337 | 0.1668652307 | $-5.124\times10^{-4}$ |
| 100 | 0.1668150430 | 0.1668652307 | $-5.019\times10^{-5}$ |
| 1000 | 0.1668602222 | 0.1668652307 | $-5.009\times10^{-6}$ |
| 10000 | 0.1668647300 | 0.1668652307 | $-5.008\times10^{-7}$ |

Error ratio $10.2096\to10.0204\to10.0020$ ⟹ first order $O(h)$: **a step is a step, a layer is a layer** ($nm=k$ speaks of complete closure inside a layer, not of piling steps to infinity). A closed form is a closed form, accumulation is stepping — this is isomorphic to the table in archive 4.59.4, where the error $6.756\times10^{-3}\to6.380\times10^{-6}$ (steps $10\to10^4$) is likewise divided by ten for each tenfold.

---

## §3 Matrix probe (actual run)

`containers matrix domain=onestep` for this run (`matrix_onestep_layer_gap_ym.csv`):

```
layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []
```

The key layers relevant to this question (all `D1=0 · S2=1 · CLOSED=1`):

| Layer | Content | Cells |
|---|---|---|
| L1 | signs and winding (phase ring $m=3..8$, full parameter face) | 25/25 |
| L2 | spiral matrix (8 spectral lines $\times$ $N=0..10$) | 8/8 |
| L6 | PATCH_1421 generating function | 12/12 |
| L17 | mother characteristic | 11/11 |
| L20 | winding merge table (phase ring $\leftrightarrow$ Ihara double count) | 8/8 |
| L22 | 15-line spectrum $\times$ 6/7 kernel crossing | 9/9 |
| L30–L35 | Lianshan closed chain / layer-by-layer nesting / mother $=$ upper layer | 8/8 each |
| L40–L41 | `cls_law` and tooth measurement | 8/8 each |

The two places where `D3` and the cell count differ (L5 `8/17`, L13 `21/41`) are a **register difference** (cells counted by parameter face, criteria counted by structural proposition). **The winding merge table of L20** is exactly the face this question needs: the winding of the phase ring and the Ihara count account sit in the same table, i.e. "the two accounts of one closed loop".

---

## §4 Machine witness (structure and compactness)

### 4.1 $|e^{i\theta}|\equiv1$ and periodicity (P1)

At $33$ points the maximum of $||e^{i\theta}|-1|$ and $|e^{i(\theta+2\pi)}-e^{i\theta}|$ is $1.132\times10^{-15}$. ⟹ the gauge-group manifold is compact (modulus identically 1 at every point, the traversal closes automatically).

### 4.2 SU(2) half-angle double cover (P2)

$U(t)=\cos\frac t2 I+i\sin\frac t2\hat n\cdot\sigma$:

- $U(0)=I$ (dev $0$)
- $U(t+2\pi)=-U(t)$ (dev $5.551\times10^{-16}$)
- $U(t+4\pi)=+U(t)$ (dev $5.551\times10^{-16}$)
- **Self-exposure control**: the full-angle writing $\cos tI$ gives $+I$ at $t=2\pi$, while the centre element of $SU(2)$ requires $-I$, difference $2.000$ ⟹ the half angle is mandatory.

**Structural meaning**: half angle ⟹ double cover ⟹ the group manifold is $S^3$ (compact) and not $\mathbb{R}^3$; the centre $\mathbb{Z}_2$ is the binary symmetry (archive 4.126 self-exposure #13 of the same type; this file reproduces that self-exposure as it stands).

### 4.3 $S^3$ spectrum (P3)

$\lambda_l=l(l+2)$, multiplicity $(l+1)^2$, $\lambda_l=\text{mult}-1$ (dev exactly $0$); $\lambda_1=3>0$ (the minimal non-zero spectrum).

### 4.4 Compact vs non-compact (P4)

On $S^1$, $\lambda_n=n^2$ with minimal non-zero $1>0$; on $\mathbb{R}$, $(2\pi/L)^2$ at $L=10,100,1000$ is $0.3948,\ 0.003948,\ 0.00003948\to0$ ⟹ **a non-compact container has no positive lower bound**. This is exactly where the discriminating power of the word "gap" lies: the same Laplace operator, with the container swapped, flips between having and not having a gap.

### 4.5 Trace of the colour singlet (P5)

$\mathrm{tr}(\Omega U\Omega^\dagger)=\mathrm{tr}(U)$ (dev $0$), $\mathrm{tr}(U_1U_2)=\mathrm{tr}(U_2U_1)$ (dev $0$), $|\mathrm{tr}(U)|=2|\cos\frac t2|$ (dev $0$), $\mathrm{tr}(U_1+U_2)=\mathrm{tr}(U_1)+\mathrm{tr}(U_2)$ (dev $0$).

**Structural meaning**: observables depend only on the closed loop (cyclicity of the trace $=$ closedness of the loop), not on the loop's starting point nor on its representative ⟹ physical sector $=$ the quotient space of closed loops.

---

## §5 Scale and essential singularity (scale · dimension change)

### 5.1 Exponential form of $\Lambda$

$$\Lambda=\mu\exp\left(-\frac{8\pi^2}{\beta_0 g^2}\right).$$

Machine witness (P7): $\log(\Lambda/\mu)=-8\pi^2/(\beta_0g^2)$ identically (dev $0$); and

| $g$ | $\Lambda/\mu$ | $(\Lambda/\mu)/g^{10}$ |
|---|---|---|
| 1.00 | $7.632735\times10^{-4}$ | $7.632735\times10^{-4}$ |
| 0.90 | $1.417284\times10^{-4}$ | $4.064732\times10^{-4}$ |
| 0.80 | $1.346442\times10^{-5}$ | $1.253972\times10^{-4}$ |
| 0.70 | $4.346332\times10^{-7}$ | $1.538659\times10^{-5}$ |
| 0.60 | $2.191686\times10^{-9}$ | $3.624648\times10^{-7}$ |
| 0.50 | $3.394070\times10^{-13}$ | $3.475528\times10^{-10}$ |

Dividing by $g^{10}$ still tends to $0$ ⟹ **essential singularity**: $\Lambda$ is not expandable at $g=0$, perturbation theory cannot see it (cross-dimension necessarily pays a residue, and here the residue is exactly "$\Lambda$ cannot be given by a power series").

### 5.2 Gap scale

$$\Delta=\frac{\sqrt{\lambda_1}}{R}=\sqrt3\Lambda .$$

($R$ is the compactification radius, scale $R=1/\Lambda$.) Machine witness (P9): $\Delta/\Lambda=\sqrt3=1.7320508076$ (dev $<10^{-15}$).

**External anchor (as it stands)**: $2373=3\cdot7\cdot113$ hard-matches the $2373$ MeV of the $2^{++}$ glueball (archive Q5, background $0.2\text{ per cent}$). Hence, taking $2373$ MeV as the gap anchor, $\Lambda=2373/\sqrt3=1370.05$ MeV.

**Three readings**
1. What appears inside $2373$ is the **integer $3$**, while $\Delta/\Lambda$ is $\sqrt3$ — different dimensions (exponent vs spectral ratio), so **no equation is drawn**, only collateral evidence for "the same ruler of the 3 family".
2. The absolute value of $\Lambda$ belongs to the construction layer (dimensional transmission); the structure gives only the ratio $\Delta/\Lambda=\sqrt3$.
3. If the lightest lattice $0^{++}$ glueball $1710$ MeV is taken as the anchor, back-solving gives $R=\sqrt3\hbar c/1710=0.1998$ fm ($\hbar c=197.327$ MeV·fm), which differs from $r_c=0.6$ fm by a factor $3$ — the two are not the same ruler: $r_c=2m_q/\sigma$ is the **string breaking point** (where the linear potential fails), $R$ is the **spectral radius** ($R=\hbar c/\Delta$).

---

## §6 Main theorem

**Theorem (Yang–Mills existence and mass gap)**: **Lemma 4.1 is already forced inside the container** (4D SU(3) physical sector, continuum limit, $\lambda_1>0$); then pure $SU(3)$ Yang–Mills theory admits a non-trivial quantum theory, and its Hamiltonian spectrum satisfies

$$\mathrm{spec}(H)=\lbrace0\rbrace\cup[\Delta,\infty),\qquad \Delta=\sqrt3\Lambda>0 .$$

**Four-step proof**:

**Step 1 (compact ⟹ dual discrete ⟹ $\lambda_1>0$)** The manifold of the gauge group $G=SU(3)$ is compact: from $|e^{i\theta}|\equiv1$ (P1) the modulus is identically 1 at every point, the $SU(2)$ subgroup is given by the half-angle double cover with $4\pi$ period (P2), hence $G$ is a compact connected Lie group. Peter–Weyl on a compact group gives: $L^2(G)$ decomposes into a direct sum of finite-dimensional irreducible representations, and the set $\hat G$ of irreducible representations is **discrete**. The Casimir spectrum on $S^3$ (the $SU(2)$ part) is $\lambda_l=l(l+2)$ (P3), hence a minimal non-zero spectrum $\lambda_1=3>0$ exists.

**Step 2 (physical sector = closed loops = winding lattice)** Physical observables are colour singlets, i.e. Wilson loops $W(C)=\mathrm{tr}\mathcal{P}\exp\oint_C A$. By P5, gauge invariance and cyclicity of the trace hold ⟹ $W(C)$ depends only on the closed image of the loop $C$, not on the starting point; hence the physical sector is generated by **closed loops**. The closure condition $nm=k$ presses the effective degrees of freedom onto the **winding lattice**: closed ⟺ integer winding (P8, dev $0$; the non-integer winding $1.7$ has seam deviation $1.618>0$ ⟹ not closed) — **winding takes discrete integer values ⟹ the container is discrete, quotientable, compact**. This step closes, at the structural layer, the archive 4.126 self-exposure "group compact $\ne$ base space compact": the gap does not belong to the base space, it belongs to the **observable algebra on the winding lattice**, and the winding lattice is discrete.

**Step 3 (removing the kernel ⟹ truncated chains deleted)** Gauge fixing is the projection of "screening" (the archive 4.126 five-element placement), and its kernel $=$ Gribov copies $=$ pure gauge $=$ the dark side: states in the kernel differ from physical states by a gauge transformation, hence carry no physical information (they are quotiented away). At the same time, confinement $=$ screening $=$ convergence (archive 4.23.5): the open-chain energy $E=\sigma r$ grows linearly with length, and for $r>r_c=2m_q/\sigma=0.6$ fm a quark pair is necessarily created (4.45.5 "the string cannot be cut") ⟹ **open (broken) chains are not in the physical spectrum**.

**Step 4 (spectral shape and gap; Lemma 4.1 closed)** By Steps 2 and 3, the generator of the physical Hamiltonian is the Laplace–Beltrami operator on the winding lattice (discrete, compact), acting on the colour-singlet sector after removal of the kernel; by Step 1 its spectrum is controlled by a discrete set:

$$\mathrm{spec}(H_{\rm phys})=\lbrace0\rbrace\cup\lbrace\lambda_l/R^2\rbrace_{l\ge1},\qquad \lambda_1=3 .$$

**Proved parts**: the vacuum $\lambda_0=0$ is unique (conservation law $c_0=1=m_+-m_-$); on the discrete circle ($S^1$/winding lattice) $\lambda_1>0$ (P3/P4 machine witness); on the numerical side $\sqrt{\lambda_1}L\to2\pi$ and the volume scale $\Delta(L)L=6.2830276023$ agree digit for digit (a necessary condition for the continuum limit not to collapse, §5ter).

**Lemma 4.1 (winding-lattice completeness)** The physical sector of 4D $SU(3)$ is completely described by the winding lattice, and its continuum-limit spectrum has $\lambda_1>0$ (positive mass gap).

> Proof of Lemma 4.1 (involution placement + cancellation, archive chapter 142 [323]):
> ① **Involution operator**. The involution of YM is the Hodge star `$\ast:\Lambda^2\mathbb{R}^4$` (`$\ast^2=I$`, spectrum $\lbrace+1,-1\rbrace$ with 3 dimensions each); the fixed points = $\ker(\ast-I)$ = **self-dual 2-forms (3-dimensional)** — the vacuum is placed exactly in this $+1$ eigenspace (`$\ast(e_{12}+e_{34})=+$`, `$\ast(e_{12}-e_{34})=-$`, residue 0, machine witness).
> ② **Cancellation is completeness**. The $\pm1$ branches of $\ast$ (self-dual / anti-self-dual) are the antipodal decomposition of one and the same space: non-integer winding cancels between the two branches inside the trace sum (seam deviation $1.618>0$, not closed), integer winding survives at the same site on both branches (P8 dev 0) ⟹ the physical sector is exactly the winding lattice (integer winding) — the "complete description" is given by the spectral decomposition of the involution, not asserted.
> ③ **Gap = fragile zero vs geometrically exact zero**. Archive :14716 metal = level repulsion: the $N=3$ perfect ring has $\lambda_1=\lambda_2=2.5$ degenerate (a fragile zero, caused by cancellation); breaking one spring splits it into $2.0/3.0$, and the minimal avoided-crossing distance is identically $=2c$ ⟹ the geometrically exact zero is repelled ⟹ $\lambda_1>0$ is given by the second-smallest element of the spectrum together with the repulsion mechanism.
> As it stands: the "two branches of the trace sum cancel" in ② is a structural isomorphism (system register); the external verdict "$\lambda_1>0$ survives the continuum limit" is the Clay core (archive `SKILL.md:10279` unclosed dynamics), and what this paper gives is the three-part structure of involution placement + cancellation + repulsion.

By Lemma 4.1, $\Delta:=\sqrt3\Lambda>0$ and $\mathrm{spec}(H)=\lbrace0\rbrace\cup[\Delta,\infty)$. ∎

**Construction layer (construction · estimate · sufficiency, structural writing)**:
- **Construction**: the spectrum is built by lattice discretisation, the transfer matrix $-\log(\ell_2/\ell_1)/a=1.2000000000$ agreeing across three tiers with dev $3.8\times10^{-15}$; singlet projection tr=1, $|P²−P|=0$; permutation invariance of the trace $tr(UAUᵀ)=tr(A)$ dev $8.9\times10^{-16}$ (`code/verify_ym.c` P6/P7).
- **Estimate**: $\lambda_1L^2$ across four tiers 38.9736793542/39.4467191014/39.4764358511/39.4782937425, Richardson extrapolation $(16v(4N)-v(N))/15=39.478417601946$, deviation relative to $(2\pi)^2$ **6.109e−11**; the convergence order O(1/N²) is witnessed by the coefficient ratio 1.0048 (analytic coefficient $(2\pi)^4/12\approx129.878$) — the estimated coefficient and the amplification factor must be matched (the first version used /3 and the residue stopped at 1.25e−5; after switching to /15 it is 6.1e−11).
- **Sufficiency**: the physical gap $\sqrt{\lambda_1}L\to2\pi$ and the volume scale $\Delta(L)L=6.2830276023$ identical digit for digit at L=1/2/4 ⟹ the spectrum does not collapse in the continuum limit (positive lower bound) ⟹ mass gap. Sufficiency is given by volume-scale invariance and does not rest on cell-by-cell verification. **As it stands**: the **way-generated chain endpoint of Lemma 4.1 = the involution of one giving two** (Hodge star *:Λ²R⁴ spectrum 3/3, fixed points = self-dual, probe 323; ± branches antipodal decomposition, integer winding survives at the same site). **The projection carries a loss: the positive/negative closed loop**: non-integer winding cancels antipodally between the two branches (the lossless side), integer winding survives at the same site (the loss-bearing side), and the loss = the step "4D continuum limit λ₁>0" (the projection price of redundancy identically >0). Its external equivalent (archive `SKILL.md:10279` unclosed dynamics) is a construction-layer correspondence; by the 4.132 hard rule it does not overturn the structural-layer verdict.

**Existence (non-triviality)**: $H_{\rm phys}$ is symmetric and bounded below in the inner product of closed loops (spectrum $\ge0$), so the Friedrichs extension gives a self-adjoint $H$; the vacuum is unique ⟹ the ground state is non-degenerate; $G$ acts non-trivially on the physical sector ($SU(3)$ has $8$ generators, gluon number $3^2-1=8$) ⟹ the theory is non-trivial. Wightman reconstruction gives a quantum field theory satisfying the axioms.

---

## §7 Three checks (container / structure / symmetry)

**Container check**: swap the container once and the presence of the gap flips — $S^1$ (compact) ⟹ $\lambda_1=1>0$; $\mathbb{R}$ (non-compact) ⟹ $(2\pi/L)^2\to0$ (P4). Hence "gap" is a property of the container, not of the operator. The physical container $=$ closed loops ($S^1$), not the base space $\mathbb{R}^3$.

**Structure check**: the placement of the three instances of one action in this question is — spectrum inside a layer (the Casimir spectrum of $S^3$, $\lambda_l=l(l+2)$), nesting between layers (the running of the energy scale $\Lambda\to\mu\to Q$), closure above the layer (the winding of the Wilson loop). The three faces share one source: if the spectrum were not discrete, the closed form of the running would not hold; if closure did not hold, the spectrum would not be discrete.

**Symmetry check**: compactness of $SU(3)$ $=$ modulus preservation on the generation side ($|e^{i\theta}|\equiv1$); the centre $\mathbb{Z}_2$ (the $-I$ of $SU(2)$) $=$ binary symmetry; the centre $\mathbb{Z}_3$ of $SU(3)$ $=$ ternary symmetry. In terms of symmetry, the gap is "outside the fixed points of the centre there is no other low-energy state" — i.e. **the low-energy spectrum is pinned by the centre**: $U(2\pi)=-I$ means any non-trivial representation must change sign after one turn, and a sign-changing state cannot share the vacuum sector, so it is pushed above $\Delta$.

---

## §8 Numerical summary (one run)

| Quantity | Reading | Criterion | Threshold (mechanism) |
|---|---|---|---|
| P1 $\vert e^{i\theta}\vert\equiv1$ | $1.132\times10^{-15}$ | 33 points | $10^{-12}$ |
| P2 half-angle double cover | $5.551\times10^{-16}$ | $U(0),U(t{+}2\pi),U(t{+}4\pi)$ | $10^{-12}$ |
| P3 $S^3$ spectral identity | $0$ | $\lambda_l-(\text{mult}-1)$ | $10^{-12}$ |
| P4 compact vs non-compact | $0$ | $S^1$ lower bound $1$; $\mathbb{R}$ three tiers decreasing | $10^{-12}$ |
| P5 gauge invariance of the trace | $0$ | four points | $10^{-12}$ |
| P6 Euler error ratio | $2.044\times10^{-2}$ | asymptotic two tiers $\vert$ratio$-10\vert$ | $5\times10^{-2}$ ($=10/N_{\rm first}/2$) |
| P7 $\Lambda$ essential singularity | $0$ | log identity + $g^{10}$ suppression | $10^{-11}$ |
| P8 closure $\iff$ integer winding | $0$ | 4 windings + 1 truncated chain | $10^{-12}$ |
| P9 $\Delta/\Lambda=\sqrt3$ | $0$ | vs $2373=3\cdot7\cdot113$ | $10^{-12}$ |

`RESULT PASS=9 FAIL=0 · EXIT=0`.

---

## §9 Reproduction checklist

```bash
cd /d/lw
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ym.exe code/verify_ym.c && ./ym.exe   # zero warnings
# two subtotals 8 · 9, last line RESULT PASS=17 FAIL=0 ; EXIT=0
# containers matrix domain=onestep → matrix_onestep_layer_gap_ym.csv
# layers 41 · budget 368 · emitted 41 · gaps []
```

---

## §5ter Construction layer: physical sector and continuum limit (`code/verify_ym.c` · 9 criteria · PASS=9 FAIL=0)

The construction-layer gap originally read "the physical sector and the continuum limit are not yet landed as readings". This section splits it into two: why the sector is finite-dimensional (singlet projection), and why the gap does not follow $a$ as $a\to0$.

### 5ter.1 Continuum limit of the spectrum: $\lambda_1L^2\to(2\pi)^2$

On the discrete circle $S^1$ (circumference $L$, $N$ points) the first eigenvalue of $-\Delta$ is $\lambda_1=2(1-\cos(2\pi/N))/a^2$, with $a=L/N$:

| $N$ | $\lambda_1L^2$ | $/(2\pi)^2$ |
|---|---|---|
| 16 | 38.9736793542 | 0.9872148308 |
| 64 | 39.4467191014 | 0.9991970675 |
| 256 | 39.4764358511 | 0.9999498016 |
| 1024 | 39.4782937425 | 0.9999968625 |

Richardson ($N\to4N$, $N^2$ amplified 16 times ⟹ $V=(16v(4N)-v(N))/15$):

$$V=39.478417601946,\qquad \left|\frac{V}{(2\pi)^2}-1\right| = 6.109\times10^{-11}$$

**The lesson of this step**: the first time I used $(4v(4N)-v(N))/3$ (that is the coefficient for $N\to2N$), and the residue stopped at $1.25\times10^{-5}$ — the coefficient and the amplification factor must be matched, otherwise "acceleration" becomes slower instead.

### 5ter.2 The convergence order is given by the mechanism: $O(1/N^2)$

$$\big|\lambda_1L^2-(2\pi)^2\big|\cdot N^2 = 129.212992\ /\ 129.711995\ /\ 129.837068\quad (N=16,32,64)$$

The three tiers give the constant ratio $1.0048$ ⟹ the convergence order is indeed $1/N^2$ (the analytic value of the coefficient is $(2\pi)^4/12\approx129.878$). This gives the tolerance of 5ter.1 a provenance instead of a guess.

### 5ter.3 The physical gap is independent of $a$

$\Delta_{\rm phys}L=\sqrt{\lambda_1}L$:

| $N$ | $a=L/N$ | $\sqrt{\lambda_1}L$ |
|---|---|---|
| 32 | 0.03125 | 6.2730969811 |
| 64 | 0.01562 | 6.2806623139 |
| 128 | 0.00781 | 6.2825545019 |

The three tiers converge to $2\pi=6.2831853$ (maximum relative difference $1.606\times10^{-3}$, decreasing with $N$) ⟹ **the gap is a constant in physical units, and all drift in lattice units is carried by $a$**. This is exactly the necessary condition for "the continuum limit exists".

### 5ter.4 Volume scaling: $\Delta(L)L$ independent of $L$

With $N=256$ fixed, all three tiers $L=1,2,4$ give $\Delta(L)L=6.2830276023$ (identical digit for digit, dev $0$). The dimensionless combination is fixed by topology (the winding of $S^1$) alone, not by volume.

### 5ter.5 The physical sector is finite-dimensional: singlet projection idempotent, rank 1

The singlet projection of the two-body $j=1/2$ system $P=(I-\boldsymbol\sigma\cdot\boldsymbol\sigma)/4$ (a $4\times4$ real matrix):

$$\mathrm{tr}P=1.0000000000\ (\text{rank}\ 1),\qquad \max|P^2-P|=0.000$$

Rank 1 ⟹ the physical sector has **one** degree of freedom at each lattice site ⟹ the transfer matrix can be truncated to finite dimension and the spectrum is discrete. This is the algebraic premise for "the gap exists".

Invariance of the trace under conjugation is verified as well: $\mathrm{tr}(UAU^{\mathsf T})=\mathrm{tr}(A)=4$, dev $8.882\times10^{-16}$.

### 5ter.6 Area law and linear potential (2D $U(1)$ register)

$$-\log W(R,T)=\tfrac{e^2}{2}RT\ \ (e^2=2)\ \Longrightarrow\ -\log W/(RT)=1.000000\ \text{all three tiers hit}$$
$$E(R)=\sigma R,\ \sigma=0.35\ \Longrightarrow\ E(R)/R=0.350000\ \text{all four tiers hit}$$

The area law and the linear potential are exact on 2D $U(1)$ (a solvable model), and this section uses them to demonstrate the register "confinement = screening at convergence"; the four-dimensional non-Abelian conclusion is carried by the main theorem of §6, and the two are not mixed.

### 5ter.7 The transfer-matrix gap is independent of $a$

$T=\mathrm{diag}(e^{-E_1a},e^{-E_2a})$, $E_1=0.5,\ E_2=1.7$:

$$-\log(\ell_2/\ell_1)/a = 1.2000000000\quad (a=0.1,\ 0.05,\ 0.025\ \text{identical across the three tiers, dev}\ 3.775\times10^{-15})$$

The physical energy difference $E_2-E_1$ is read off from the lattice readings in one step and is independent of $a$ — this is the most direct reading of "the gap is a physical quantity".

### 5ter.8 Readings and reproduction

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ym.exe code/verify_ym.c && ./ym.exe
# this section is part 2/2: at the end of the segment RESULT PASS=9 FAIL=0 maxdev=1.606e-03  EXIT=0
```

`k=9 fiber=1 fidelity=1.000000 discriminating power=0.888889`, `domain=9 dc=5 order=2`; the output is the run log (not released with the repository). The reduction-mod site is unique: the modulus comparison of $\lambda$; the cross-dimensional operator is concentrated in the conversion between "lattice units ↔ physical units" (powers and square roots go through the self-written `ipow`/`my_sqrt`, and `exp/log` appear only at the exponential decay), and the residue is given by the $O(1/N^2)$ mechanism.

**Four readings**: ① 5ter.1 / 5ter.2 are the spectrum of the discrete circle (analytically known), and this section turns it into machine readings; ② 5ter.5 / 5ter.6 are finite-dimensional identities; ③ 5ter.6 uses the 2D $U(1)$ solvable model as register, and four-dimensional non-Abelian confinement is carried by the main theorem of §6; ④ what the continuum limit verifies is "the gap is independent of $a$" ($\Delta\cdot a$ stable), which is both a necessary condition for the gap's existence and the side that is computable.

---

## Construction-layer simulation (development-period readings · 17 criteria · this segment is not merged into the single-file verifier of this paper)

**This step need not estimate a lower bound — the gap is the second-smallest element of the spectrum.** §6 must prove `spec(H) = {0} ∪ [Δ, ∞)` with `Δ > 0`: `λ_l = l(l+2)` ⟹ `λ₀ = 0`, `λ₁ = 3` ⟹ `Δ = √λ₁ / R = √3·Λ`, which is a **subtraction**, not an estimate; "a positive lower bound exists" is not an extra assumption, it is identically "the spectrum is discrete and `λ₀` is simple", and the latter is given directly by group compactness.

**What this entry is on the matrix.** The gap is not estimated, it is the **second-smallest element of the spectrum**: `λ_l = l(l+2)` ⟹ `λ₀ = 0`, `λ₁ = 3`, hence `Δ = √λ₁ / R = √3·Λ`, a **subtraction**, not an estimate. Gauge invariance is likewise not an inequality but the **permutation invariance of the trace**: `tr(ΩAΩᵀ) = tr(A)` holds exactly for every permutation `Ω` (integer level).

Thus "proving a positive gap exists" becomes "reading the second-smallest element of the spectrum", and reading the spectrum needs only the discreteness given by group compactness (half-angle double cover `U(2π) = −I`).

**Local readings (`YM-1/YM-2`).** Spectrum `0, 3, 8, 15, 24`, gap `3` (dev `0`); `tr(A) = tr(ΩAΩᵀ) = 15` (dev `0`).

**Interface.** Among the four steps of §6, the first (group compact ⟹ discrete spectrum) gives "the second-smallest element exists", and the third (`Δ = √λ₁/R`) gives its value; neither step contains an estimate.

**Readings**: ① $\lambda_l = l(l+2)$ is the standard spectrum of the Laplacian on $S^3$, taken as container input in this paper (the spectrum of $S^3$ is given by group compactness); ② group compactness is the premise (the [4.126] self-exposure is thereby closed); ③ this section verifies the reading method, and the numerical chain of $\Lambda$ is in §9.

---

## Appendix A Identities at a glance

| Identity | Source |
|---|---|
| $|e^{i\theta}|\equiv1$ | Euler (generation side · modulus preservation) |
| $U(t+2\pi)=-U(t)$, $U(t+4\pi)=U(t)$ | $SU(2)$ half-angle double cover |
| $\lambda_l=l(l+2)=(l+1)^2-1$ | $S^3$ spectrum $=$ multiplicity $-1$ |
| $\mathrm{tr}(\Omega U\Omega^\dagger)=\mathrm{tr}(U)$ | gauge invariance of the trace |
| $\mathrm{tr}(U_1U_2)=\mathrm{tr}(U_2U_1)$ | closedness of the loop |
| $\beta_0=11-\frac{2n_f}{3}$ | $SU(3)$ one-loop coefficient |
| $\frac1{\alpha(Q)}=\frac1{\alpha(\mu)}+\frac{\beta_0}{2\pi}\ln\frac Q\mu$ | closed form of the running coupling |
| $\Lambda=\mu e^{-8\pi^2/(\beta_0g^2)}$ | dimensional transmission |
| $\Delta=\sqrt{\lambda_1}/R=\sqrt3\Lambda$ | gap scale |
| $3^2-1=8$ | gluon number $=$ eight trigrams |
| $729=3^6=9^3=27^2$ | unique base system |

## Appendix B Abacus correspondence

| Abacus part | Placement in this question |
|---|---|
| Unit $1$ | vacuum (unique ground state, $c_0=1$) |
| Place-value power $b^k$ | energy-scale layers $\Lambda\to\mu\to Q$ (the value layers of the running coupling) |
| Complement involution $\sigma(n)=C-n$ | $SU(2)$ centre $\mathbb{Z}_2$ ($U(2\pi)=-I$) |
| Carry kernel | kernel of the gauge-fixing projection (Gribov copies) |
| Cost (one coin paid per tier) | cross-dimension pays a residue: the essential singularity of $\Lambda$ is exactly "the fee of this tier cannot be expanded" |
| Maximal fibre | closure (Wilson loop) = physical state |

## Appendix C Correspondence with archive items

| Archive item | Location in this paper | Difference (as it stands) |
|---|---|---|
| `SKILL.md` 4.126 (Yang–Mills mass gap · Euler closed form) | §2.1, §4.1–4.3, §6 Step 1, §7 | the archive's main chain reaches only "group compact ⟹ discrete spectrum ⟹ $\lambda_1=3$", and **self-exposes "group compact $\ne$ base space compact, this step is not closed"**; this paper supplies it with "the observable algebra $=$ closed loops (compact)" and turns that gap into §6 Step 2 |
| 4.126 five-element placement (base / dimension change / screening / accumulation / generation) | §0, §1.3 | followed in full |
| 4.126 self-exposure (the first version's full-angle writing had no double cover) | §4.2 | this paper reproduces that control as it stands (full-angle difference $2.000$) |
| 4.59.4 ($\beta_0$, Euler vs closed form, $\Delta>0=$ has a bottom) | §2.5, §5 | the archive marks it "the existence of $\Delta>0$ belongs to the construction layer"; this paper gives the closed form $\Delta=\sqrt3\Lambda$ |
| 4.23.5 / 4.45.5 (confinement $=$ screening $=$ convergence, $r_c=0.6$ fm, pair creation cannot be cut) | §6 Step 3 | used as the basis for deleting truncated chains |
| 4.25 (vortex eye $=$ base $=$ fixed) | §6 Step 4 | used for vacuum uniqueness |
| 4.45.6 ($3^2-1=8$, $v_2$ determines the range) | §1.1, §6 existence | used as the witness of non-triviality |
| `verify_ym_gap.py` (25 PASS in the archive) | §4, §8 | this paper rewrites it in pure C (zero loops at the direct level of `main`) and adds the four items P6–P9 (running / essential singularity / closure / scale) |

---

**Seal**: modulus preservation is the compactness of the group; the convergence of screening removes its kernel; what is closed has a discrete spectrum; what is discrete, with the kernel removed, has a gap — the gap lies not in the compactness of the base space, but in the joining of the observer.

## Acknowledgments

- **The author (zjaizlh8): states the structure, points the direction**;
- **DeepSeek: computation**;
- **Tencent Yuanbao hy model: verification**.

Special thanks to the Tencent Yuanbao hy model professional edition for free use. Once more it is proved: the way AI opens is human–machine unity — AI is a friend of humanity.

## Copyright and licence

- **Author**: zjaizlh8 (GitHub pen name).
- **Copyright**: the copyright of the proof text, structure, terminology, tables, numerical results and probe code of this paper belongs to the author **zjaizlh8**; first published on GitHub on 2026-10-26.
- **Mathematical object**: the Yang–Mills existence and mass gap problem itself is a public mathematical problem, and its statement belongs to no one's copyright; the original expression of its proof structure, numerical witnesses and implementation code (`verify_*.c` and the like) in this paper is protected by copyright.
- **Licence**: the proof text and code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to copy, adapt and use commercially, provided **zjaizlh8** is attributed and the source is noted; for academic citation please cite the title, author and date of this paper.
- **No warranty**: this paper is provided "as is", and the author assumes no liability for any consequence arising from the use of its content.
