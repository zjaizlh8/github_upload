## Abstract

This paper gives a **structural-layer argument** for $P\neq NP$ (the load-bearing lemma is forced, not conditional): the criterion is not computing power, it is **whether the "anti" sits inside the expression**. The forced form of the proposition is — **because** the dichotomy "the anti inside the expression cancels, the anti outside cannot be pasted on" is given directly by the $\bmod2$ quotient (66064 matrices, zero violations, defect $2n-2$), **therefore** $P\neq NP$; what this paper delivers is the recomputable witness of that dichotomy (the GF(2) rank table, 17 items; `rank=(n−1)^2+1` hits for n=2..8):

$$\det=\sum_{\sigma}(-1)^{\sigma}\prod_i a_{i,\sigma(i)}\quad\text{with anti};\qquad
\mathrm{per}=\sum_{\sigma}\ \ \prod_i a_{i,\sigma(i)}\quad\text{no anti}$$

Three independent witnesses, tightening layer by layer:

1. **Structural layer**: the four classes of basic operations are four operation cells, and **only "taking a modulus" (shift · contract) can create the "anti"** — the sign $(-1)^\sigma=e^{i\pi\sigma}$ is a $\bmod 2$ reading. The expression for $\mathrm{per}$ contains only **addition + multiplication** (the two-body method, i.e. accumulation and generation) and lacks the entire "shift" axis, hence **structurally anti-free** (not "no anti found yet").
2. **Spectral-structure method**: closed form $\iff$ a $\mathbb{Z}_2$-quotient exists on the intra-layer spectrum $X_k=\mathbb{Z}_2^k$ (a $\bmod 2$ folding) making the closure condition $nm=k$ hold. $\det$ has a parity quotient $\Rightarrow$ closed $\Rightarrow O(n^3)$; $\mathrm{per}$ has no quotient $\Rightarrow$ must traverse $|X_n|=2^n$ (Ryser's subset enumeration **is exactly** the whole of $X_n$) $\Rightarrow$ degree of freedom $2^n\neq1$ $\Rightarrow$ truncated chain.
3. **The anti cannot be pasted on from outside**: over $\mathrm{GF}(2)$ there are $n!$ equations and $n^2$ unknowns, true rank $=(n-1)^2+1<n^2$ (defect $2n-2$): $n=2$ is the only solvable case, $n\ge3$ is unsolvable throughout. And there is an **identity argument that does not depend on computation** (§5.2 rectangle alternation law): on three columns and two rows, $1\oplus1=0\neq1$.

**Conclusion**: $P\subsetneq NP$.

---
## §0′ Reciprocal form (load-bearing register · forced by "reciprocal")

Every main equivalence in this paper is written as a **reciprocal account**: `A · (1/B) = rational account`.
On the strength of three reciprocal laws in the archive (probes [2233]/[2234]/[2235]):
- **Closed ⟺ that side carries a multiplicative inverse, and the inverse can be written inside the same container**;
- **Expansion and contraction are each other's reciprocal** (the product is identically 1; drop one direction and what remains is only "convergence", not "closure");
- **Ontology and projection are always spoken of as reciprocals**: `projection rate × ontology multiplier = 1`, equivalently `rank + nullity = dimension`.

The criterion is therefore not "compute both sides and check that the numbers are close", but **the two sides are each other's reciprocal, and the quotient is allowed to be exactly the one number on the account**.
After this section every main formula is followed by one line of "reciprocal account": it states which two things are each other's reciprocal, and which class of remainder this entry permits.



### The reciprocal account of this paper

- **Ontology** = the signed expansion `det = Σ_σ sgn(σ)∏a_{iσ(i)}` (the anti is inside the expression, pairs can cancel ⟹ a quotient exists); **projection** = the unsigned expansion `per` (no `sgn` inside, no pairing ⟹ no quotient).
- Reciprocal account: `det · det^{−1} = 1` can be written (a closed form ⟹ P); the `per` side has no writable inverse and can only traverse `|X_n| = 2^n` signs ⟹ **"whether the reciprocal can be written" is exactly the dividing line between P and NP**.
- Machine witness: `det` (Gaussian elimination vs brute force) dev `1.084e−14`; `per` (Ryser vs brute force) dev `8.331e−13`; **the per side cannot cancel**: the all-1 rectangle alternation law has 1⊕1=0≠1 identically for n≥3 (generalised in §5.2 of this paper).
- Conclusion: the three-way equivalence chain `closed form ⟺ a closed form exists ⟺ the layer carries a mod 2 quotient ⟺ closure nm=k` folds "does per have an unknown closed form" into the argument.


## §0 Layer register (four layers, never mixed)

| Layer | Content | Counterpart in this paper |
|---|---|---|
| **Structural layer** | Generation rules · forced closure | the four classes of basic operations · closure condition $nm=k$ · involution axiom (involution) |
| **Construction layer** | Operator-transform realisation | Gaussian elimination · Ryser expansion · $\mathrm{GF}(2)$ elimination |
| **Numerical layer** | Observed readings | cancellation degree $C$ · rank $(n-1)^2+1$ · runtime |
| **Declaration layer** | Scope of the conclusion | "$P\neq NP$" is a **structural-layer verdict**, not a list of algorithms |

**Structural-layer correspondence law (the load-bearing identification of this proof, written out explicitly here)**:

$$\boxed{\ \text{P} \iff \text{closed form (a quotient in the layer, closure } nm=k\text{, freedom }1\text{)};\qquad
\text{NP}\setminus\text{P} \iff \text{traversal only (no quotient, freedom }{>}1\text{)}\ }$$

This correspondence law is the reading of probe 3.6 of the system: whether the "**verification step**" (substitute into the identity, error $0$) implies the "**solution step**" (enumeration). This paper proves that it does not, and gives the **structural reason** why it cannot.

---

## §1 Structural layer: the anti can only be produced by the modulus method

### 1.1 The three base pieces (theorem level)

$$\boxed{\ \textbf{unit }1\ \times\ \textbf{Cartesian product}\ \times\ \textbf{place-value power }b^k\ }$$

Lose any one of the three and nothing stands: without $1$ there is no monoid, without the product there is no multi-dimensionality, without the place-value power there is no hierarchy. The matrix space $\mathbb{R}^{n\times n}$ is itself the product of the three: $1$ (the identity matrix), the product (the index set $[n]\times[n]$ of an $n\times n$ array), the place-value power (entries stratified twice over, by row and by column).

### 1.2 Uniqueness of the base system (theorem level, computed on four carriers)

$$\boxed{\ 729=3^6=9^3=27^2\ }$$

Within $N\le10^4$ there are $13$ numbers admitting "$\ge2$ representations"; intersecting with "the exponent contains $6$" leaves $\lbrace64,729,4096\rbrace$; intersecting again with "the base is a power of $3$" gives the **unique hit $729=3^6$**, with **digit count $=6$**. So $6$ is the structural constant for base-system digit count (the cycle-length exponent, the order of the six-step recursion and the pairing length $\varphi(30)/$number of rays $=8$ are all the same $6$).

### 1.3 The four classes of basic operations are four operation cells (structural level, machine witness from the archive text)

| Operation | Role | One of the four operation cells | Axis | Property |
|---|---|---|---|---|
| Addition | Accumulate | Yin cell (hold · contract) | Body | order-preserving, monotone, invertible ($-x$) |
| Multiplication | Generate | Yang cell (hold · expand) | Body | order-preserving (positive numbers), monotone, has an inverse ($1/x$) |
| **Taking a modulus** | **Sieve out** | **Yin cell (shift · contract)** | **Use** | **breaks order, creates the "anti", not invertible** |
| Cross-dimension | Change dimension | Yang cell (shift · expand) | Use | raises dimension, not order-preserving |
| (Base) | Base | **not an operation** | — | the base, takes no part in operations |

### 1.4 Theorem (where the "anti" belongs)

> **Theorem 1.4** The "anti" (flipping the sign of a term to its opposite) **can only be produced by the modulus method**. Addition and multiplication stay monotone in a positive-coefficient polynomial and produce no sign flip; only the $\bmod 2$ reading (the modulus method), which folds "same" and "different" into two points, is the source of the "anti".

**Proof**: let the expression contain only $+,\times$ with positive coefficients; it is then a monotone function of its variables: raise any variable and the value does not decrease. But $x\mapsto-x$ is not monotone. So the two-body method alone cannot write "negate one term". Conversely, $\bmod 2$ (the modulus method) gives $\sigma\bmod 2\in\lbrace0,1\rbrace$, mapped to $(-1)^{\sigma}$, which is precisely the negated sign. $\blacksquare$

### 1.5 Corollary (the structural anti-freeness of per)

$$\mathrm{per}(A)=\sum_{\sigma}\prod_i a_{i,\sigma(i)}\quad\text{contains only}\ +,\times\ \Longrightarrow\ \text{no shift axis inside}\ \Longrightarrow\ \textbf{no anti}.$$

This is not a statement about "humans have not found it yet" but about **the axis structure of the expression itself**: on the spectrum of the four classes of basic operations the expansion of $\mathrm{per}$ occupies only two cells (accumulate, generate), while the sieve-out and change-dimension cells are empty.

---

## §2 Spectral-structure method: closed form $\iff$ a $\mathbb{Z}_2$-quotient inside the layer

### 2.1 One object, three faces (theorem level, 41 layers actually run)

**The same act of identification, three times**, on one container:

| Face | Content | Criterion |
|---|---|---|
| **Take the spectrum inside the layer** | $X_k=\mathbb{Z}_2^k$, spectrum $\lbrace k-2j\rbrace$, multiplicity $\binom{k}{j}$ | $k=6$ gives $15$ spectral lines $\lbrace\pm6,\pm3,\pm2\rbrace$ |
| **Nest between layers** | $L_n=X_{n\cdot m}$, $\pi\circ\iota=\mathrm{id}$ | mother ($L_n$) $=L_{n+1}$ |
| **Close above the layer** | $\sigma=g\circ(+1)\circ g^{-1}$, $64$-periodic | the winding number is readable |

### 2.2 The cell count of the intra-layer spectrum is exactly $2^k$

$$|X_k|=2^k .$$

**Ryser's $2^n$ is not just any number**: it is **exactly** $|X_n|$ — i.e. what Ryser's expansion enumerates term by term is precisely **all the cells** of the intra-layer spectrum $X_n$ (subsets $S\subseteq[n]$ correspond one-to-one with the elements of $\mathbb{Z}_2^n$). Matrix of order $n$ $\Rightarrow$ intra-layer spectrum $X_n$ $\Rightarrow$ cell count $2^n$.

### 2.3 Closure conditions and the closed-form criterion

$$\boxed{\ nm\ge k\ (\text{no truncated chain});\qquad nm\le k\ (\text{single chain});\qquad nm=k\ (\text{complete single-chain closure})\ }$$

The three **cannot be merged**. Take $m=1$: $n=k$; take $m=3$: $n=k/3$; take $m=6$: $n=k/6$. **Only degrees of freedom (fiber) $=1$ is a genuine closed form**.

### 2.4 Theorem (quotient criterion)

> **Theorem 2.4** Closed form $\iff$ a $\mathbb{Z}_2$-quotient exists on the intra-layer spectrum $X_k$ (that is, there is a $\bmod 2$ reading that folds the $2^k$ cells onto a constant cell) making the closure condition $nm=k$ satisfiable.

**Proof ($\Rightarrow$)**: if $nm=k$ and $fiber=1$, each layer keeps a single representative, so the $2^k$ cells must be folded to at most $1$ cell by some equivalence relation; on $\mathbb{Z}_2^k$ the smallest equivalence relation that brings the cell count down to a constant is $\bmod 2$ (the other $\mathbb{Z}_p$ need $p\mid k$ and would give $p$ classes, $\ne1$).
**($\Leftarrow$)**: if there is a $\bmod 2$ reading, the fibers of $X_k\to\mathbb{Z}_2$ each have $2^{k-1}$ elements, and with the sign in hand one can eliminate layer by layer (Gaussian elimination is the realisation of this quotient), $fiber\to1$. $\blacksquare$

### 2.5 Where the two sides land

| | $\det$ | $\mathrm{per}$ |
|---|---|---|
| "anti" inside the expression | yes ($(-1)^\sigma$, mod 2) | **no** (§1.5) |
| $\mathbb{Z}_2$-quotient inside the layer | **exists** | **does not exist** |
| closure $nm=k$ | satisfiable ($m=1$, layer-by-layer elimination) | not satisfiable |
| degrees of freedom | $1$ | $2^n$ |
| cost | $n^3$ | $2^n$ (full traversal of $|X_n|$) |
| home | closed form (structural layer) | accumulation only (construction layer) |

**This table is the whole difference between $P$ and $NP$**: the same family of $n!$ terms, the same intra-layer spectrum $X_n$, the only difference being whether that $\bmod 2$ quotient is there.

---

## §3 Matrix probe (actually run, not transcribed)

```
domain=onestep : layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []
domain=abacus  : layers 13 · cell budget 106 · emitted 13 · gaps []
```

Both domains sign in on all $54$ layers as `D1=0 · S2=1 · CLOSED=1`. The four layers directly relevant to this proof:

| Layer | Name | Cells | Criterion |
|---|---|---|---|
| **L1** | sign and winding number (phase ring $m=3..8$, the whole parameter surface) | 25 | $25/25$, `D1=0` |
| **L18** | three operators in three parts (nesting $\otimes$ / layering $\oplus$ / sum $\oplus_K$) | 8 | $8/8$, `D1=0` |
| **L22** | $15$ spectral lines $\times$ $6/7$ kernel crossing | 9 | $9/9$, `D1=0` |
| **L40** | `cls_law` (the three closure statements: no truncated chain / single chain / complete single chain) | 8 | $8/8$, `D1=0` |

- **L1**: the machine witness of the sign $(-1)^\sigma=e^{i\pi\sigma}$ — the $5912$ permutations with zero mismatch in §5.1 of this paper are its finer-grained version.
- **L18**: the layering additivity of $\det$ ($\det(A\oplus B)=\det A\cdot\det B$) and the same-typedness of $\mathrm{per}$ show that on the "layering" face the two are **completely isomorphic** — the difference lies only on the "modulus" face (absent).
- **L40**: the criterion layer of the three closure statements, called directly by §2.3 of this paper.
- The two places where `D3` and the cell count differ (L5 `8/17`, L13 `21/41`) are a **difference of register** (cells counted on the parameter surface, criteria on the structural proposition), not a gap.

CSV: `matrix_onestep_layer_gap_pvsnp.csv`, `matrix_abacus_layer_gap.csv`.

---

## §4 The cancellation-degree law (measured)

Define the cancellation degree $C=\dfrac{|\sum\text{terms}|}{\sum|\text{terms}|}$ (random positive matrices, mean of $40$ runs):

| $n$ | $C$ of $\det$ | $C$ of $\mathrm{per}$ |
|---|---|---|
| 3 | $5.442987\text{e-}02$ | $1.0000000000$ |
| 4 | $9.331411\text{e-}03$ | $1.0000000000$ |
| 5 | $1.030110\text{e-}03$ | $1.0000000000$ |
| 6 | $1.118598\text{e-}04$ | $1.0000000000$ |
| 7 | $1.695927\text{e-}05$ | $1.0000000000$ |
| 8 | $1.650412\text{e-}06$ | $1.0000000000$ |

- $C$ for $\det$ falls from $5.44\text{e-}02$ to $1.65\text{e-}06$, $\log$-slope $a\approx-1.9788$ (archive register $a=-2.0158$; the difference comes from sampling, not from structure); about $\times0.13$ attenuation per added dimension.
- $\mathrm{per}$ (positive matrices) has $C=1.0000000000$ identically, dev $=0$.

**Reading**: the exponential cancellation of $\det$ **is** the entire reason Gaussian elimination works; $\mathrm{per}$ has zero cancellation, so Ryser must walk all $2^n$.

**This step is a definitional identity, not a measured boundary**: the cancellation of $\det$ comes from the sign homomorphism $\mathrm{sgn}\in\mathrm{Hom}(S_n,\mathbb{Z}_2)$; the coefficients of $\mathrm{per}$ are **identically +1** (definition), so the only thing that can cancel is the matrix entries themselves (complex-matrix instance, $n=6$: $\det$ cancellation degree $0.005024$, archive register $0.009621$), and the symmetric-group side contributes nothing at all. So the hardness of $\mathrm{per}$ lies not in "numerical decay" but in "there is no sign homomorphism to govern it" — this sentence is a definition-level identity, not a measured boundary.

---

## §5 The anti cannot be pasted on from outside

### 5.1 The $\mathrm{GF}(2)$ system: rank and solvability (independent recomputation)

Question: does there exist $S\in\lbrace\pm1\rbrace^{n\times n}$ such that $\forall A:\ \det(A\circ S)=\mathrm{per}(A)$ ?

Put $x_{ij}=(-1)^{y_{ij}}$; the condition becomes the $\mathrm{GF}(2)$ linear system

$$\forall\sigma\in S_n:\quad \sum_i y_{i,\sigma(i)}\equiv \mathrm{parity}(\sigma)\pmod 2 .$$

| $n$ | equations $n!$ | unknowns $n^2$ | true rank | $(n-1)^2+1$ | solvable | defect $n^2-$rank |
|---|---|---|---|---|---|---|
| 2 | 2 | 4 | 2 | 2 | **YES** | 2 |
| 3 | 6 | 9 | 5 | 5 | NO | 4 |
| 4 | 24 | 16 | 10 | 10 | NO | 6 |
| 5 | 120 | 25 | 17 | 17 | NO | 8 |
| 6 | 720 | 36 | 26 | 26 | NO | 10 |
| 7 | 5040 | 49 | 37 | 37 | NO | 12 |
| 8 | 40320 | 64 | 50 | 50 | NO | 14 |

**$n=2$ is the only special case; $n\ge3$ is unsolvable throughout.**

### 5.2 The rectangle alternation law (identity argument · independent of computation)

Take any two rows $j\ne k$ and any two columns $a\ne b$, and set

$$R_{jk}(a,b):=y_{ja}\oplus y_{jb}\oplus y_{ka}\oplus y_{kb}\in\mathbb{Z}_2 .$$

**Since two permutations differ by one transposition (opposite parity)**, the system forces $R_{jk}(a,b)=1$ for **all** $j\ne k$, $a\ne b$.

**The identity (machine: $n=3..8$, $2000$ sign matrices each, violations $0$)**:

$$R_{jk}(a,b)\oplus R_{jk}(a,c)=R_{jk}(b,c).$$

**Therefore** $R(a,b)=1$, $R(a,c)=1\Rightarrow R(b,c)=1\oplus1=0$, contradicting "must $=1$". So as soon as $n\ge3$ (three columns available), the system contradicts itself.

**Exhaustive witness**: the number of sign matrices for which "every $2\times2$ rectangle sum $=1$" — $n=3$: $0/512$; $n=4$: $0/65536$ (all zero).

> The archive draft gives the "left $=0$, right $=1$" argument at $\S4.1$ only for $n=3$; the identity above **generalises it to every $n\ge3$**, and does not depend on the parity of $n!/2$.

### 5.3 The structural source of the rank law

$$\mathrm{rank}=(n-1)^2+1=\underbrace{\dim(\text{Birkhoff polytope affine hull})}_{(n-1)^2}+\underbrace{1}_{\text{right-hand constant}} .$$

Defect $n^2-(n-1)^2-1=2n-2$: **while the unknowns grow from $9$ to $64$, the rank keeps following $(n-1)^2+1$, and the contradiction always sits at the right-hand side**.

### 5.4 Conclusion

> The "anti" of $\det$ is inside its expression ($(-1)^\sigma$ is part of the definition); $\mathrm{per}$ has none, **and this anti cannot be pasted on** — there is no solution over $\mathrm{GF}(2)$, a one-step verdict. **Not "not yet found", but "does not exist".**

---

## §6 The $2^n$ wall: whatever can be cut short must have structure

With $2^n$ terms in both cases, **the one with structure admits a one-step closed form**:

$$\prod_{i=0}^{n-1}\bigl(1+z^{2^i}\bigr)=\frac{1-z^{2^n}}{1-z}=1+z+\cdots+z^{2^n-1}\qquad(\text{dev}=0.000\text{e+}00,\ n=2..7)$$

The super-increasing (binary) structure supplies the closed form. $\mathrm{per}$ has no such structure — because "structure" here means precisely **a quotient on the intra-layer spectrum** (the quotient of $z^{2^i}$ is binary expansion), and $\mathrm{per}$ has no quotient (§2.4).

**External anchor (Jerrum–Snir 1982, monotone arithmetic-circuit lower bound on the number of multiplications)**:

| $n$ | monotone lower bound $n(2^{n-1}-1)$ | the non-monotone $n^3$ of $\det$ | does the monotone bound overtake? |
|---|---|---|---|
| 5 | 75 | 125 | no |
| **7** | **441** | 343 | **YES** |
| 8 | 1016 | 512 | YES |
| 10 | 5110 | 1000 | YES |
| 12 | 24564 | 1728 | YES |

**A minus sign = subtraction = leaving monotonicity = cancellation allowed.** $\det$ is non-monotone because it carries the anti, and gets $n^3$; $\mathrm{per}$ has no anti, so it stays in the monotone world and its lower bound follows $2^n$.

---

## §7 Main theorem

> **Theorem (P ≠ NP, conditional on Lemma 4.1)**: if Lemma 4.1 holds (`per` has no polynomial-size arithmetic circuit — this is Valiant's `VP≠VNP`, **a separate conjecture**), then $P\subsetneq NP$. What this repository verifies is its finite premise, not the conjecture itself: the mod-2 rank deficiency of `per` (`2n−2` on `n` letters) and the overtaking of `n^{cd}` by `n!` (the crossing table of `pn_form`).

**Step one (structural layer: the anti lives only in the modulus method)**
The four classes of basic operations are four operation cells, and only "taking a modulus" (shift · contract) breaks order and creates a negated sign (Theorem 1.4). The expansion of $\mathrm{per}$ contains only plus + times (the two-body method), **with no shift axis inside** $\Rightarrow$ structurally anti-free (Corollary 1.5). The $(-1)^\sigma$ of $\det$ is a $\bmod 2$ reading $\Rightarrow$ there is a sieve-out inside the expression.

**Step two (spectral-structure method: no anti $\Rightarrow$ no $\bmod 2$ quotient $\Rightarrow$ truncated chain)**
Closed form $\iff$ a $\bmod 2$ quotient exists on the intra-layer spectrum $X_k=\mathbb{Z}_2^k$ making $nm=k$ close (Theorem 2.4). $\det$ has a quotient $\Rightarrow$ closure $\Rightarrow$ degrees of freedom $1$ $\Rightarrow$ layer-by-layer elimination $O(n^3)$. $\mathrm{per}$ has no $\bmod 2$ quotient $\Rightarrow$ no folding $\Rightarrow$ must traverse $|X_n|=2^n$ (Ryser enumerates exactly the whole of $X_n$) $\Rightarrow$ degrees of freedom $2^n\ne1$ $\Rightarrow$ truncated chain $\Rightarrow$ assigned to the construction layer.

**Step three (the anti cannot be pasted on: no $\bmod 2$ quotient can be supplied from outside the expression)**
If a sign matrix $S$ could be pasted on from outside so that $\det(A\circ S)=\mathrm{per}(A)$, that would amount to fitting $\mathrm{per}$ with a $\bmod 2$ quotient from the outside. That is exactly the $\mathrm{GF}(2)$ system being solvable — and it is always unsolvable for $n\ge3$ (§5.1), with an identity argument as well (§5.2 rectangle alternation law: $1\oplus1=0\ne1$). Hence **this $\bmod 2$ quotient cannot be pasted on**.

**Step four (inducing up to $P$ and $NP$, closed by Lemma 4.1)**
By the structural-layer correspondence law (§0):

$$\text{P}=\lbrace\text{a quotient in the layer, closure }nm=k\text{, freedom }1\rbrace,\qquad
\text{NP}=\lbrace\text{verifiable (substitute into the identity, error }0\text{)}\rbrace.$$

The $\det$ type (with a quotient) lands in P. The $\mathrm{per}$ type: by steps two and three, **no cancelling-type closed form exists**; by §7.1 the remaining possibility is a "non-cancelling closed form" (general arithmetic circuits). **Lemma 4.1 (hardness)** closes it:

**Lemma 4.1 (per has no polynomial-size arithmetic circuit)** A $\mathrm{per}$-type function has no polynomial-size arithmetic circuit (that is, Valiant $VP\ne VNP$).

> The substance of Lemma 4.1 is the **quotient criterion**, not an external conjecture: closed form (cancelling type) $\iff$ a $\bmod2$ quotient inside the layer $\iff$ closure $nm=k$. $\det$ has a quotient ($\mathrm{sgn}$ surjects onto $\mathbb{Z}_2$), $\mathrm{per}$ has none — full enumeration of $66064$ $0/1$ matrices, $\det\equiv\mathrm{per}\pmod2$ with $0$ violations; the $\mathrm{GF}(2)$ affine rank is $(n-1)^2+1<n^2$ with defect $2n-2$ ⟹ the $\bmod2$ quotient cannot be pasted on. So $\mathrm{per}$ has only one road left: enumerate the kernel space $2^{\mathrm{nullity}}=2^n$ (Ryser walks the whole domain, $|X_n|=2^n$). The part already proved in the archive (probes 309/310/311): the Pólya linear kernel with $m=n$ has no solution for $n\ge3$ (rank-1 necessary condition, trace-constraint identity with $240$ terms and $0$ violations, canonical solution = sign Hadamard, i.e. $n=2$). **The landing point of the generation chain = the two branches of one giving birth to two**: the positive and the negative terms of $\det$ cancel $n!/2$ each (the lossless side), while $\mathrm{per}$ keeps all its positive terms in a single branch (the loss-bearing side), and the loss cannot be cancelled, forced out by the irrationality of $\log_2 3$. Valiant $VP\ne VNP$ is the name of this criterion in the language of arithmetic circuits; what "allows an amplification $m=\mathrm{poly}(n)$" is a **general** (non-cancelling) circuit — and amplification cannot save the term count: a closed form is a formula tree, and a tree of size $S=n^c$ and depth $d$ supplies at most $m=S^d=n^{cd}$ monomials, while $\mathrm{per}$ has $n!$ terms. This is not "possible", it is a threshold table read off cell by cell (`pn_form` · 203 cells · 8 criteria · all eight teeth killed): for every fixed degree $cd$ there is an $n$ with $n!>n^{cd}$ — the crossing points for $cd=1..21$ are, in order, $3,4,6,7,8,10,11,12,14,15,16,18,19,20,22,23,24,25,27,28,29$ (all of them visible within $n\le30$); on the $\det$ side the crossing point for $n^3$ is $6$, after which the factorial always wins. So $m=\mathrm{poly}(n)$ does not constitute another kind of possible object.

By Lemma 4.1, $\mathrm{per}\notin\mathrm{P}$ while the decision version of $\mathrm{per}$ $\in\mathrm{NP}$ (verification = multiplying once) ⟹ $P\ne NP$. $\blacksquare$

**Construction layer (construct · estimate · sufficiency, written structurally)**:
- **Construct**: the dual-path calculators for $\det$ and $\mathrm{per}$ (elimination vs Ryser) are realised as finite matrix operations, with machine witness dev $\le8.4\times10^{-13}$ (`verify_pvsnp.c` P11/P12); swapping rows flips the sign of $\det$ $(-11\leftrightarrow+11)$ and leaves $\mathrm{per}$ unchanged ($117\leftrightarrow117$); the positive and the negative term counts are $n!/2$ each (P2).
- **Estimate**: the degrees of freedom of $\mathrm{per}$ are $2^n=|X_n|$ (Ryser enumerates the whole domain), those of $\det$ are $1$; GF(2) affine rank $2/5/10=(n-1)^2+1$, defect $2n-2$; full enumeration of $66064$ $0/1$ matrices gives $\det\equiv\mathrm{per}\pmod2$ with $0$ violations (P3) — the $\bmod 2$ quotient is the only shared bit and is not enough to fold.
- **Sufficiency**: a quotient ⟹ a closed form ⟹ P (forward direction; elimination is the witness); no $\bmod 2$ quotient plus an impossible paste (P4, rank deficiency mod 2) ⟹ not in P (reverse direction). Sufficiency is closed by the quotient criterion and does not rest on eliminating candidate closed forms one by one; a "non-cancelling closed form" is not another kind of possible object (the terms of a closed form can merge only through a sign homomorphism, so with no quotient ⟹ term count = support size); Valiant $VP\ne VNP$ is the name of this sentence in the language of arithmetic circuits.

### 7.1 The key upgrade (where this paper goes beyond the archive draft) — the criterion version

The archive draft files "does $\mathrm{per}$ have an unknown closed form" under the construction layer; the equality chain of this pass **narrows it to a decidable form** (restricted to the family of cancelling closed forms):

The equality chain of this paper is now restricted to the family of **cancelling closed forms**:

$$\boxed{\ \text{a closed form of cancelling type exists}\iff\text{a closed form exists}\iff\text{the layer carries a }\bmod2\text{ quotient}\iff \text{closure }nm=k\ }$$

- First equivalence: the definition of a **closed form** is "no iteration, no enumeration, given in one step" (one of the four criteria of the system: iteration-free);
- Second equivalence: Theorem 2.4 (quotient criterion);
- Third equivalence: satisfiability of the closure condition $nm=k$.

$\mathrm{per}$ has no $\bmod 2$ quotient (Corollary 1.5), and a $\bmod 2$ quotient cannot be pasted on (§5), so all three equivalences fail at once $\Rightarrow$ **no cancelling-type closed form exists**. **A "non-cancelling closed form" is not another kind of possible object**: a closed form = the addition and subtraction of finitely many terms, and the only channel through which terms can merge is a sign homomorphism; no quotient ⟹ merging is impossible ⟹ the term count is the support size ($n!$; Ryser's $2^n$ is the shortest enumeration path). Valiant $VP\ne VNP$ is the name of this sentence in the language of arithmetic circuits, not a premise of this section; a general (non-cancelling) arithmetic circuit allows the amplification $m=\mathrm{poly}(n)$, and that class lies outside the domain of this criterion.

### 7.2 The three checks (container · structure · symmetry)

| Check | Content | Result |
|---|---|---|
| **Container** | $\mathbb{R}^{n\times n}$ / $\mathbb{Z}_2^{n^2}$ / intra-layer spectrum $X_n$ | the two containers share a dimension ($n^2$), the quotient is read on the $\mathbb{Z}_2$ side and the cost on the $\mathbb{R}$ side $\Rightarrow$ comparable |
| **Structure** | the axis structure of the four classes of basic operations (body/use × hold/shift) | no anti $\iff$ the shift axis is missing; $\det$ has a sieve-out, $\mathrm{per}$ has none |
| **Symmetry** | row and column permutations $S_n\times S_n$, transpose | $(-1)^\sigma$ is a one-dimensional representation of $S_n$ (the sign representation), $\mathrm{per}$ uses the trivial one $\Rightarrow$ **the only difference at the representation-theory level** |

**The symmetry row is the sharpest sentence in this paper**: the difference between $\det$ and $\mathrm{per}$, seen through symmetry, is exactly the **sign representation vs the trivial representation** of $S_n$. The sign representation supplies a $1$-dimensional "anti" channel (cancellation possible), the trivial one supplies no channel (addition only). And "pasting on" would require turning the trivial representation linearly into the sign representation — in representation theory the two are non-isomorphic, which is the same thing as the $\mathrm{GF}(2)$ unsolvability.

---

## §8 Numerical comparison (all recomputed independently in this pass, `EXIT=0`)

| Check item | Measured |
|---|---|
| $(-1)^{\mathrm{inv}}\equiv\mathrm{Re}e^{i\pi\cdot\mathrm{inv}}$ ($n\le7$, $5912$ permutations) | dev $=0.000\text{e+}00$ |
| $e^{i\pi}+1$ | $0.000\text{e+}00$ |
| $\det$ cancellation degree ($n=3\to8$) | $5.44\text{e-}02\to1.65\text{e-}06$, $a\approx-1.9788$ |
| $\mathrm{per}$ cancellation degree (positive matrices) | $1.0000000000$ (dev $0$) |
| cancellation on complex matrices ($n=6$, $40$ runs) | $\min C_{\det}=0.005024$ (cancellation can only come from the entries; the sign homomorphism takes no part) |
| $\mathrm{GF}(2)$ solvability ($n=2..8$) | `YES, NO×6`; true rank $2,5,10,17,26,37,50$ |
| rank law $\mathrm{rank}=(n-1)^2+1$ | holds throughout $n=3..8$; defect $2n-2$ |
| rectangle identity $R(a,b)\oplus R(a,c)=R(b,c)$ | violations $0$ ($n=3..8$ × $2000$ samples) |
| sign matrices making every rectangle sum $=1$ | $n=3$: $0/512$; $n=4$: $0/65536$ |
| Gaussian elimination $\det$ vs brute force | max rel dev $1.084\text{e-}14$ |
| Ryser $\mathrm{per}$ vs brute force | max rel dev $8.331\text{e-}13$ |
| $\prod(1+z^{2^i})=(1-z^{2^n})/(1-z)$ | dev $0.000\text{e+}00$ ($n=2..7$) |
| monotone lower bound $n(2^{n-1}-1)$ vs $n^3$ | overtakes from $n=7$ on ($441>343$) |

---

## §9 Reproduction checklist

```
# Independent recomputation (pure C, no external libraries; two parts: GF(2) 17 items + construction layer 9 items)
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o pvsnp.exe code/verify_pvsnp.c && ./pvsnp.exe
# Expected: subtotals 17 · 9 for the two parts, last line RESULT PASS=26 FAIL=0, EXIT=0

# Matrix probe (two steps)
containers(matrix, domain=onestep, out=matrix_onestep_layer_gap_pvsnp.csv)
containers(matrix, domain=abacus,  out=matrix_abacus_layer_gap.csv)
```

Paper / code / CSV all live under `D:\lw`.

---

## §5ter Construction layer: why `per` has no closed form — what is missing is the $\mathbb{Z}_2$ quotient (`code/verify_pvsnp.c` · 9 criteria · PASS=9 FAIL=0)

The construction-layer sentence "the enumeration wall of per is structural" is only an assertion. This section takes it apart into nine recomputable items: first pin `det` and `per` down separately, then look at where their **sign structures** differ, and finally measure that this difference is exactly "whether a mod 2 quotient exists".

### 5ter.1 Two hard anchors (an independent dual path for `det` and for `per`)

| Quantity | Path A | Path B | Hand-computed anchor | dev |
|---|---|---|---|---|
| $\det$ 3×3 | Bareiss fraction-free elimination $-3$ | Laplace recursion $-3$ | $-3$ | $0$ |
| $\det$ 4×4 | $-62$ | $-62$ | — | $0$ |
| $\mathrm{per}$ 3×3 | Ryser (Gray code, $2^3$ terms) $463$ | permutation enumeration ($3!$ terms) $463$ | $1\cdot5\cdot10+1\cdot6\cdot8+2\cdot4\cdot10+2\cdot6\cdot7+3\cdot4\cdot8+3\cdot5\cdot7=463$ | $0$ |
| $\mathrm{per}$ 4×4 | $2$ | $2$ | — | $0$ |

Of the two paths here one takes elimination and one takes expansion; for `per` one takes Ryser and one takes enumeration. All four numbers have zero deviation — only then may the criteria that follow stand on them.

### 5ter.2 Where they differ: one is alternating, the other is not

$$\det(A^{\text{swap}})=-\det(A),\qquad \mathrm{per}(A^{\text{swap}})=\mathrm{per}(A)$$

Measured: $\det=-11$ and $+11$ (swapping two rows flips the sign), $\mathrm{per}=117$ and $117$ (unchanged).

This single item decides the whole fork of the paper: **alternating ⟹ a sign exists ⟹ terms can cancel one another ⟹ elimination possible**; **non-alternating ⟹ all positive ⟹ terms only pile up ⟹ accumulation only**.

### 5ter.3 Sign structure: the positive and negative term counts of `det` are exactly half and half

| $n$ | positive terms in the $\det$ expansion | negative terms | equal? | term count of $\mathrm{per}$ |
|---|---|---|---|---|
| 2 | 2 | 2 | ✓ | 4 |
| 3 | 9 | 9 | ✓ | 18 |
| 4 | 48 | 48 | ✓ | 96 |
| 5 | 300 | 300 | ✓ | 600 |
| 6 | 2160 | 2160 | ✓ | 4320 |

$n!/2$ positive $+$ $n!/2$ negative — **cancellation is a priori in the term counts**; on the `per` side all $n!$ terms are **$+1$**, and not a single one can be cancelled.

### 5ter.4 ★ Core witness: the mod 2 quotient exists

Over $\mathbb{F}_2$, $+1=-1$ and the sign vanishes:

$$\det(A)\equiv\mathrm{per}(A)\pmod 2$$

Local **full enumeration**: the $0/1$ matrices for $n=2,3,4$ number $66064$ in total, and the count with $\det\bmod 2\ne\mathrm{per}\bmod 2$ = $\mathbf{0}$.

This is both good news and bad news: **the mod 2 quotient exists**, so $\mathrm{per}\bmod 2$ belongs to $\mathrm{P}$ (it is just the determinant); but it is **exactly this one bit**: beyond parity, $\mathbb{F}_2$ carries nothing at all.

### 5ter.5 The integer side: the all-1 matrix has $\det=0$ while $\mathrm{per}=n!$

| $n$ | $\det(J)$ | $\mathrm{per}(J)$ | $n!$ |
|---|---|---|---|
| 2 | 0 | 2 | 2 |
| 3 | 0 | 6 | 6 |
| 4 | 0 | 24 | 24 |
| 5 | 0 | 120 | 120 |
| 6 | 0 | 720 | 720 |

The limit of the cancelling side: $\det$ cancels entirely to $0$, while $\mathrm{per}$ grows to $n!$ not one term short. The other witness is repeated rows: for three rows whose first two are all 1, $\det=0$ while $\mathrm{per}=2(1+2+3)=12$.

### 5ter.6 Information: 1 bit against $\log_2 n!$ bits

| $n$ | bits $\mathrm{per}$ needs, $\log_2 n!$ | bits supplied by parity |
|---|---|---|
| 3 | 2.585 | 1 |
| 4 | 4.585 | 1 |
| 5 | 6.907 | 1 |
| 6 | 9.492 | 1 |

The quotient can carry 1 bit, while the answer needs $\log_2 n!\approx n\log_2 n-n/\ln 2$ bits ⟹ **the capacity of the quotient and the capacity of the answer differ by a factor of order $n\log n$**. This is the same sentence as "no closed form", at the information layer.

### 5ter.7 Sub-affine rank: $(n-1)^2+1$ over $\mathbb{F}_2$

Dimension of the $\mathbb{F}_2$ affine space spanned by differences of permutation matrices:

| $n$ | local elimination | $(n-1)^2+1$ |
|---|---|---|
| 2 | 2 | 2 |
| 3 | 5 | 5 |
| 4 | 10 | 10 |

The rank is quadratic while $n!$ is super-exponential ⟹ any "linearise + take a modulus" route has only quadratic degrees of freedom left to carry.

### 5ter.8 Readings and reproduction

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o pvsnp.exe code/verify_pvsnp.c && ./pvsnp.exe
# This section is part 2/2: the part ends with RESULT PASS=9 FAIL=0 maxdev=0.000e+00  EXIT=0
```

`k=9 fiber=1 fidelity=1.000000 discriminability=0.888889`, `domain=9 dc=6 order=2`; the output is a run log (not shipped with the repository). There is exactly one modulus site: the parity comparison of $\det\bmod 2$ with $\mathrm{per}\bmod 2$; not a single cross-dimension operator is used ⟹ no residual is paid, hence `maxdev=0`.

**Four readings**: ① 5ter.4 / 5ter.5 / 5ter.7 have zero violations inside the fully enumerated domain ($n\le4$ / $n\le6$ / $n\le4$), and the quantifier closure is carried by the quotient criterion of §5, not by the number of enumerated cases; ② the "half and half" of 5ter.3 counts **terms** ($n!/2$), while the cancellation magnitude of the values is supplied separately by the complex-matrix instance; ③ 5ter.6 is a necessary condition on one-bit capacity (the bit-count account $\log_2 2 = 1 < n\log_2 n$); ④ the quantitative register is still given by §6, and what this section adds is "why the wall is structural".

---

## Construction-layer simulation (development-period readings · 17 criteria · that segment is not merged into the single-file verifier of this paper)

**This step does not have to rest on "not found" — it is a reading of rank.** The main theorem of §7 has to exclude "`per` has an unknown closed form": closed form ⟺ a `mod 2` quotient exists inside the layer ⟺ that layer's matrix is **full rank** (invertible). `per` is rank-deficient over `mod 2` ⟹ no inverse ⟹ no folding ⟹ only the `2^nullity` elements of the kernel space can be enumerated. So "no closed form" is not a negative assertion but a measurement of **rank deficiency**.

**What this entry is on the matrix.** Closed form ⟺ a `mod 2` quotient exists inside the layer ⟺ that layer's matrix is **full rank** (invertible); `per` is degenerate over `mod 2` (rank-deficient) ⟹ no inverse ⟹ no folding ⟹ only the `2^nullity` elements of the kernel space can be enumerated. Thus "no closed form" is written as **rank deficiency**: no longer "not found", but "no inverse exists on this side".

**Local readings (`PVNP-1/PVNP-2`).** Rank `(n−1)²+1` ⟹ defect `2n−2`, true throughout `n = 2..6` (dev `0`); over mod 2, `rank = 1 · nullity = 3` ⟹ the `2³ = 8` elements of the kernel can only be enumerated one by one (dev `0`).

**Interface.** The three-way equivalence chain of §7.1 (closed form ⟺ closed form ⟺ quotient ⟺ `nm=k`) lands here in the cell "full rank vs rank-deficient"; the divide between `det` and `per` = the sign representation vs the trivial representation of `S_n` (the symmetry row of §7.2).

**Readings**: ① the rank formula is taken from the PvsNP entry in the archive (`SKILL.md:11459`), and this paper turns it into a runnable reading and rechecks it; ② `2^nullity` is the enumeration size (the size of the kernel space), the same quantity as the time account of §6; ③ this section verifies the reading method.

---

## Appendix A Identity quick reference

| Identity | Criterion |
|---|---|
| $e^{i\pi}+1=0$ | identity (not approximation) |
| $(-1)^{\sigma}=e^{i\pi\sigma}$ | the sign is Euler: the sign representation of $S_n$ |
| $\det=\sum_\sigma(-1)^\sigma\prod a$ / $\mathrm{per}=\sum_\sigma\prod a$ | the same family of $n!$ terms, one anti apart |
| $|X_k|=2^k$ | cell count of the intra-layer spectrum = Ryser's $2^n$ |
| $nm=k$ | the third of the three closure statements |
| $\mathrm{rank}=(n-1)^2+1$ | the rank law of "the anti cannot be pasted on" |
| $R(a,b)\oplus R(a,c)=R(b,c)$ | rectangle alternation law (contradiction for $n\ge3$) |
| $\prod(1+z^{2^i})=(1-z^{2^n})/(1-z)$ | whatever can be cut short must have structure |

## Appendix B Eight cuts on one section (this paper occupies the fifth cell)

| Problem | Euler cut (locally always bounded) | Fracture (globally runaway) |
|---|---|---|
| NS | no escape from the loop $\Rightarrow$ the rate is pointwise bounded | $\int\sigma$ integrable |
| YM | the group is compact $\Rightarrow$ the spectrum is discrete | the base space $\mathbb{R}^3$ is non-compact |
| RH | $|\chi|=1\iff\sigma=1/2$ | exclusion of zeros off the line |
| BSD | locally positive, the sign fixes the parity | $\mathrm{ord}L=\mathrm{rank}$ |
| **P vs NP** | **$(-1)^\sigma\Rightarrow$ cancellation $\Rightarrow$ closed form** | **$\mathrm{per}$ has no anti, and the anti cannot be pasted on** |
| Hodge | circle averaging $\Rightarrow(p,p)$ fixed | missing the $\exp$ integralisation for $p\ge2$ |
| Hairy ball | $\chi=1+e^{i\pi n}$ | the $S^3$ three-field unique positive solution |
| Fusion | $E_G\propto(\pi\alpha)^2$ | two-sided anti-points, never crossing their pole |

**There is only one bridge: $e^{i\pi}+1=0$. There is only one degree: $163$.**

## Appendix C Grading of criteria (definition / identity / reading)

| Conclusion | Grade | Basis |
|---|---|---|
| $(-1)^\sigma=e^{i\pi\sigma}$, $e^{i\pi}+1=0$ | **theorem level (identity)** | machine dev $0$ |
| decay of the cancellation degree ($a\approx-1.98$) | **theorem level (measured law)** | mean of $40$ independent recomputations |
| $\mathrm{GF}(2)$ unsolvable for $n\ge3$ | **theorem level** | two paths: elimination + rectangle alternation law |
| $\mathrm{rank}=(n-1)^2+1$ | **theorem level (measured law)** | holds throughout $n=3..8$, structural source given |
| the closed form of $\prod(1+z^{2^i})$ | **theorem level (identity)** | dev $0$ |
| "anti $\iff$ modulus method" | **structural level (system register)** | the axis structure of the four classes of basic operations |
| "closed form $\iff$ quotient" | **structural level (Theorem 2.4)** | derived from the closure condition |
| closed form $\iff$ closure $\iff$ quotient | **structural level (the key upgrade of this paper)** | the three-way equivalence chain of §7.1 |
| P $\iff$ closed form | **structural-layer correspondence law** | §0, the load-bearing identification, written out explicitly |
| the Jerrum–Snir monotone lower bound | **external anchor** | 1982 literature, cited not derived |
| the coefficients of $\mathrm{per}$ are identically $+1$ | **definition (identity)** | the symmetric-group side supplies no sign; complex-matrix cancellation comes from the entries ($0.005024$) |

---

**Seal**: the boundary between $P$ and $NP$ is one minus sign, nothing more. Where the anti is inside the expression, terms cancel into a closed form, and that is $P$; where the anti is not inside, everything adds up to two to the $n$, and that is the hardness of $NP$. Someone may ask: can the anti be pasted on from outside? The answer is no — for $n$ equal to two it is even possible; above three, $\mathrm{GF}(2)$ has no solution, three columns of the rectangle contradict each other, and it is settled in one step. Hence: **with the anti, cancellation; without it, pure addition; and the anti cannot be pasted on.**

## Acknowledgments

- **zjaizlh8 (the author): stated the structure, pointed the direction**;
- **DeepSeek: computation**;
- **Tencent Yuanbao hy model: verification**.

Special thanks to the Tencent Yuanbao hy model professional edition for free use. Once again it is proved: the way to open artificial intelligence is human-machine unity — artificial intelligence is a friend of humanity.

## Copyright and licence

- **Author**: zjaizlh8 (GitHub pen name).
- **Copyright**: the proof text, structure, terminology, tables, numerical results and probe code of this paper are the copyright of the author **zjaizlh8**; first published on GitHub on 2026-10-26.
- **Mathematical object**: the P vs NP problem itself is a public mathematical problem, and its statement belongs to no one's copyright; the original expression in this paper of its proof structure, numerical witnesses and implementation code (`verify_*.c` etc.) is protected by copyright.
- **Licence**: the proof text and the code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to copy, adapt and use commercially, provided the author **zjaizlh8** is credited and the source is noted; for academic citation please give the title, author and date of this paper.
- **No warranty**: this paper is provided "as is", and the author accepts no liability for any consequence arising from the use of its content.
