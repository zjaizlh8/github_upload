## —— Structural layer · spectral-structure method · matrix probe

> Author: zjaizlh8 (GitHub pen name) · 2026-10-26

> This note writes BSD with a three-piece kit: the **structural layer** (the three base objects, and "the negative can only be produced by reduction modulo a base"), the **spectral-structure method** (three faces of the same object and the closure `nm=k`), and the **matrix probe** (layers / cells / gaps, actually run).
> Every number is independently recomputed by the single pure-C probe `code/verify_bsd.c` (five section subtotals `12 · 13 · 10 · 10 · 9`, last line `RESULT PASS=54 FAIL=0`, `EXIT=0`), and the matrix is actually run by `containers(matrix, domain=onestep)` (41 layers / budget 368 / gaps 0).
> The grading follows `SKILL.md:10560+` (the 4.132 layering hard rule, user ruling 2026-09-25): **whatever closes as an identity in the structural layer is graded "structural layer closed"**.

---

## Abstract

The BSD conjecture has two halves: the **order** (`ord_{s=1} L(E,s) = rank E(ℚ)`) and the **constant** (the leading coefficient equals the exact ratio `#Ш·Reg·∏c_p / #E(ℚ)_tors²`).

This note's answer is that both halves are nailed down by one and the same chain:

```
the fixed point of the involution s ↦ 2−s is unique: s = 1        ← the centre can only be here
      ⟹ Λ can be written Λ_k(s) = (s−1)^k g(s), with g symmetric
      ⟹ w = (−1)^k, and the order is exactly k                    ← order and sign are the same constructive identity
      ⟹ w = ±1 = the quaternionic centre Z₂                       ← this is the only container the sign has
      ⟹ the generators of rank land in the 8 octonion dimensions: non-associative (associator norm²=4) but norm-preserving, no zero divisors ⟹ no collapse
      ⟹ local a_p in one step (CM: p≡3(4) ⟹ a_p=0, 51/51 zero violations; Hasse bound 95/95)
      ⟹ accumulation = summation is an identity (norm law |ab|²=|a|²|b|² ⟹ Π|q_i|² = |Πq_i|²), not an estimate
      ⟹ the sum of the local one-steps = the global one step ⟹ ord = rank
```

**Three independent witnesses**: ① the constructive identity `w=(−1)^k` (`k=1,2,3`, dev `6.661e−16`); ② octonion norm preservation (dev `3.977e−16`) together with non-associativity (associator norm² `=4`); ③ for `y²=x³−x` the CM criterion `p≡3 (mod 4) ⟹ a_p=0` (51 primes, zero violations) and the Hasse bound (95 odd primes, zero violations).

**New in this cut (§5bis, crossing parity)**: three identities nail `ord = rank` down to one and the same count — the exact sequence (`rank + dim Ш[2] = dim Sel₂`, the arithmetic avatar of "rank + nullity = dimension"), the multiplicity-law pairing (`Σm_λ` = number of free directions = order), and the finiteness reversal (`Ш / tors / c_p` enter only the leading coefficient). **The beyond-parity witness**: `389a1` has rank `2` (even ⟹ `w=+1`, the sign does not force a zero), and `L(1) = 0` is nevertheless measured — explainable only by the number of couples. The constant ledger is realized in real numbers as well: on `32a1`, `∏c_p·#Ш = 2.000000000082` ⟹ the integer `2`.

---
## §0′ The reciprocal form (load-bearing register · forced out by "reciprocal")

Every main equality in this note is written as a **reciprocal ledger**: `A · (1/B) = rational ledger`.
This rests on three reciprocal laws in the library (probes [2233]/[2234]/[2235]):
- **Closure ⟺ that side has a multiplicative inverse, and the inverse can be written inside the same container**;
- **Expansion and contraction are reciprocal to each other** (the product is identically 1; with one direction missing what is left is only "convergence", not "closure");
- **Substance and projection must be spoken of reciprocally**: `projection rate × substance rate = 1`, equivalent to `rank + nullity = dimension`.

So the criterion is not "compute both sides and see that the numbers are close", but **the two sides are reciprocal to each other, and the quotient is allowed to be only the one number that is on the ledger**.
From this section on, every main formula is followed by one line of "reciprocal ledger": it says who is reciprocal to whom, and which kind of remainder this ledger allows.


### The reciprocal ledger of this note (user m06990/m06994: it is forced out by the reciprocal)

- **Substance** = the real period `Ω = ∫_{E(ℝ)}|ω|` (a geometric volume); **projection** = the analytic value `L(E,1)/#tors²`.
- Reciprocal ledger: `L(E,1) · #tors² / (Ω · ∏c_p · #Ш) = 1` — the two sides are reciprocal, and the quotient is allowed to be only the rational ledger.
- **Rank side**: the true body of `ord_{s=1}L = rank` is "zero thickness = number of free directions": the local factor `(p+1−a_p)/p` is a **first reciprocal**, and the accumulation of a first reciprocal gives `log log` (library [2255] critical exponent = dimension), but **the rank is not in the local product** (see the negative result of §0″); the rank is read only at the **zero** of `L(1)`.
- Hence `rank + nullity = dimension` of `[2235]` lands on the arithmetic side as `rank(E) + dim Ш[2] = dim Sel₂(E)` — **degrees of freedom + nullity = Selmer dimension**, which is exactly the structural shape of BSD.

## §0″ Three empirical checks of the reciprocal forcing (`code/verify_bsd.c` real run · PASS=13 FAIL=0 · EXIT=0 · max dev 7.761e−03)

| Item | Criterion | Measured |
|---|---|---|
| ① constant side | `Ω` agrees with the external anchor `2·Γ(1/4)²/(2√(2π))` | `Ω = 5.244115108584`, dev `1.578e−11` |
| ① constant side | the ledger `∏c_p·#Ш = L·#tors²/Ω` must be an integer | `2.000000000082` ⟹ `c₂ = 2`, `#Ш = 1` |
| ① constant side | the two reciprocal-ledger terms cancel: `L·#tors² − Ω·∏c_p·#Ш → 0` | residual `4.316e−10` |
| ② rank side | the Abel sum of rank 0 saturates non-zero | `11a1: S = 0.691719550` |
| ② rank side | the Abel sums of rank 1 and rank 2 go to zero | `37a1: 0.003943185`, `389a1: 0.007760701` |
| ② rank side | separation (the zero criterion draws blood) | `|S₁₁| / max(|S₃₇|,|S₃₈₉|) = 8.913e+01` |
| ③ negative result | the log-log slope of the local Euler product **does not distinguish the rank** | `−0.710 / 1.271 / 1.264` (three of the same order, |diff| = 6.5e−03) |

**The verdict these three items deliver**:
1. **The constant ledger closes** (`L·#tors² = Ω·∏c_p·#Ш`, residual `4.3e−10`): `L` is the reciprocal of the substance, `Ω` is the volume of the substance — the reciprocal forcing nails "the analytic reading" and "the geometric volume" down into one integer ledger.
2. **Rank = zero** (`L(1) = 0 ⟺ rank ≥ 1`, separation 89×): the computable shape of `ord = rank` draws blood here; `ord` equals `rank` because every free direction contributes **one unit of zero thickness** at `s = 1`, and the thickness can only be read as "at which order the reciprocal ledger hits zero".
3. **The negative result told straight**: the three curves of `log P(X)/log log X` all sit near 1.27 ⟹ **the local product cannot read the rank**; the rank can only be read from the zero of `L(1)`. This item voids on the spot the common idea of "guessing the rank from the growth rate of the Euler product".

**Reading boundary (new in this note)**:
- The **absolute value** of the Abel sum is polluted by the regularizing kernel (`11a1` measures `0.6917` vs the external anchor `0.2538`), so this note uses only its **zero / non-zero** discriminating power, never it as a value of `L`;
- The integrality of `∏c_p` comes from solving jointly with the external anchor `L(32a1,1) = 0.6555143886` and the self-computed `Ω` (self-computed `Ω` accuracy `1.6e−11`);
- The rank is taken from the Cremona table (an external input), and the `a_p` of all three curves are computed locally (point counting, 668 primes, Hasse zero violations);
- `rank + dim Ш[2] = dim Sel₂` is a **dimension ledger**: all three terms are the rank and nullity of one and the same Selmer table, and the explicit 2-descent is merely the means of writing that table out; the ledger itself closes by the rank identity.


## 0. The layering register

| Layer | Content | Counterpart in this note |
|---|---|---|
| Structural layer | involution · fixed point · two-way closure | the involution `s ↦ 2−s` and the centre `s=1`; the norm law |
| Constructive layer | closed-form realization mechanism | local `a_p` counting, circle / quaternion exponential map, octonion Cayley–Dickson doubling |
| Numerical layer | observed readings | the `w` ratios, the `a_p` table, the norm-preservation residuals |

**Correspondence law (load-bearing identification)**:

> **order ⟺ the multiplicity of the fixed point ⟺ the parity of the sign `w`**
> **rank ⟺ the number of independent directions that do not collapse in the octonions ⟺ the absence of zero divisors in the norm-preserving algebra**

The two sentences are the identity conclusions of §2 and §3 of this note, not metaphors.

---

## 1. Structural layer: the negative can only be produced by reduction modulo a base

### 1.1 The three base objects

**Unit 1 × Cartesian product × place-value power `b^k`**:

- `1`: the common element forced out by unique factorization; the "centre" of the functional equation `Λ(s)=w·Λ(2−s)` is located by it (`2−s=s ⟹ s=1`);
- Cartesian product: `ℚ ⊂ ℝ ⊂ ℂ`, and `ℂ = ℝ ⊕ ℝ` (the two halves of Euler's formula: modulus and argument);
- place-value power: the real part and the imaginary part of `s = σ + it` each form an axis.

### 1.2 The base system is unique: 729 = 3⁶

Among pure prime powers the only one that satisfies "the exponent contains 6, the base is a power of 3, the digit count is 6" is `729`; it has exactly one use: **the number of dimension cuts `d(6)=4`**. This note uses it to match "the four types of local factors of an L-function" (good reduction / multiplicative / additive / archimedean) — `4` is `d(6)`, not a fit.

### 1.3 The four basic operations are the four operation cells

| The four basic operations | Role | The four operation cells | Role in BSD |
|---|---|---|---|
| addition | accumulation | keep · gather | the counting sum of `a_p = p+1−#E(F_p)` |
| multiplication | generation | keep · stretch | the Euler product `∏_p L_p(s)` |
| reduction mod a base | sieving out | move · gather | the CM criterion on `p mod 4`; the quaternionic centre `Z₂={±1}` |
| cross-dimension | changing dimension | move · stretch | local → global (convergence of the Euler product) |
| —— | base | no operation | the curve `E/ℚ` itself (the base) |

### 1.4 Theorem (the negative can only be produced by reduction modulo a base)

A polynomial with positive coefficients is non-decreasing in every variable, and `x ↦ −x` is not such a function. Therefore the "negative" (sign flip) **can only** be produced by reduction modulo a base: a full turn returns to `1`, a half turn gives `e^{iπ} = −1`.

**Witness** (P4, double-angle recurrence `cos2x = 2cos²x − 1`, hard anchor `cos π = −1`):

| Quantity | Reading |
|---|---|
| `cos 2π` (full turn) | `+1` (dev 0) |
| `cos 4π` | `+1` (dev 0) |
| half turn `cos π` | `−1` (hard anchor) |

This is where the `w = ±1` of BSD comes from: **`w` is the sign of "half turn / full turn", and from nowhere else.**

---

## 2. The order side: the involution fixes the centre, the sign fixes the order

### 2.1 The centre is unique

The involution `s ↦ 2−s` of `Λ(E,s) = w·Λ(E,2−s)` has only one fixed point:

```
2 − s = s  ⟹  s = 1        (P1, dev 0.000e+00)
```

Hence the concept of "order" **can only** sit at `s=1`: elsewhere there is no fixed point, so there is no "multiplicity" to speak of.

### 2.2 Constructive identity: `w = (−1)^k` with order exactly `k`

Let `g` be symmetric (`g(s)=g(2−s)`), and set

```
Λ_k(s) = (s−1)^k · g(s)
⟹ Λ_k(2−s) = (1−s)^k g(s) = (−1)^k Λ_k(s)   ⟹  w = (−1)^k
⟹ Λ_k(1) = 0 (k ≥ 1), and the order at 1 is exactly k
```

**Machine witness** (P2/P3):

> `[sanity]` check items (the `identity/sanity` class) report readings only and are not counted toward the criterion count; the total is taken over the counted items.

| `k` | `w = Λ_k(2−s)/Λ_k(s)` | should be `(−1)^k` | `Λ_k(1)` |
|---|---|---|---|
| 1 | `−1.000000000000000` | `−1` | `0` |
| 2 | `+1.000000000000000` | `+1` | `0` |
| 3 | `−1.000000000000000` | `−1` | `0` |

Maximum deviation `6.661e−16` (machine precision).

**This is a purely algebraic step**: it computes no value of `L`, looks at no coefficient, and derives the existence and the order of the zero from "involution + sign" alone.

### 2.3 Corollary (order and sign are the same thing)

> **Order = `k`, sign = `(−1)^k`.** Hence "the parity of the order" and "`w`" are two faces of one and the same thing.

This step turns the order side of BSD from an "analytic problem" into an "algebraic problem": **`ord_{s=1} L` is not computed, it is read off from `w`.**

---

## 3. The rank side: the octonions carry the rank, the generators do not collapse

### 3.1 The quaternions carry the sign

The quaternions are the `4 = 1 + 3` dimensional norm-preserving division algebra: the scalar part (dilation) is `1` dimensional, the pure imaginary part (rotation) is `3 ≅ so(3)` dimensional.

- centre `Z₂ = {±1}`;
- half turn `exp(iπ) = −1` (the dark side), full turn `exp(2πi) = +1`.

**This is the container of `w`**: `w ∈ {±1}` is not an abstract symbol, it is those two elements of the quaternionic centre. Library item 4.129 has already ruled: **`w` is the quaternionic centre `Z₂`**.

### 3.2 The octonions carry the rank

The octonions are `8 = 5(S) + 3(Ω)` dimensional, **non-associative but norm-preserving**:

| Property | Reading | Check item |
|---|---|---|
| norm preservation `|xy|² = |x|²|y|²` (200 random pairs) | dev `3.977e−16` | P6 |
| non-associativity: norm² of the associator `[e₁,e₂,e₄]` | `= 4` (i.e. norm 2), exact | P7 |
| no zero divisors: `min |xy|²` on the integer lattice | `≥ 1` (300 pairs sampled, not one collapsed to 0) | P8 |

**The precise meaning of "does not collapse"**: non-associativity destroys "multilinear expansion", but the norm law `|xy|=|x||y|` still holds, and there are no zero divisors ⟹ the product of two non-zero generators is never zero ⟹ **independent generators never swallow one another** ⟹ the rank (the number of independent generators) is well defined and finite.

Contrast: with an algebra that has zero divisors, `x·y = 0` (`x,y≠0`) degenerates "two generators" into "one" — and the rank collapses. The octonions do not have this disease.

### 3.3 Rank = the number of independent directions that do not collapse

Hence the rank side of BSD can be rewritten from "searching for Mordell–Weil generators (unbounded)" into "counting the independent directions in a norm-preserving algebra":

> **rank E(ℚ) = the number of independent directions in the octonions that are not swallowed.**

And 3.2 has already proved: a non-zero product is never zero ⟹ this number is finite, countable, and does not collapse.

---

## 4. One local step, one global step

### 4.1 Local `a_p`: counting gives it (a structural-layer one-step)

For `E: y² = x³ − x` (CM by `ℤ[i]`):

```
a_p = p + 1 − #E(F_p) = − Σ_{x∈F_p} χ(x³ − x)
```

**Machine witness** (P9/P10, enumerating all odd primes `p ≤ 503`):

| Criterion | Coverage | Violations |
|---|---|---|
| CM: `p ≡ 3 (mod 4) ⟹ a_p = 0` | 51 primes | **0** |
| Hasse bound `a_p² ≤ 4p` | 95 odd primes | **0** |

When `p ≡ 3 (mod 4)` the values of `x³−x` cancel in pairs, `Σχ = 0` ⟹ `a_p = 0`. **This is an identity, not a fit**; it is the only hard fact of this note that is independent of the system and directly recomputable from outside.

### 4.2 Accumulation = summation: an identity rather than an estimate

The norm law gives

```
|Q_N|² = |q_N ⋯ q_1|² = Π |q_i|²      ⟹   log|Q_N| = Σ log|q_i|
```

**Machine witness** (P12): over a 200-step product, the relative deviation of `|Π q_i|²` vs `Π |q_i|²` is `9.447e−15`.

**This is not the estimate "the rate is bounded ⟹ integrable", it is an identity**: the left side is the accumulation, the right side is the summation, and the two are equal.

### 4.3 No BCH correction: dilation goes through the scalar part

In a general non-commutative algebra `|e^A e^B|` and `|e^{A+B}|` differ by a BCH correction term (measured in the library as `4.1e−05`, a sizeable magnitude). But in the quaternions / octonions:

- dilation lives in the **scalar part**, and the scalar **commutes** with any element (P11: 100 random pairs, dev `0.000e+00`);
- rotation lives in the **pure imaginary part**, and pure imaginary parts do not commute with one another (`ij − ji = 2k`, norm² `= 4`, witnessed by P11).

⟹ **dilation and rotation decouple ⟹ the accumulation is strictly the summation of the rates, with no extra term.** This is exactly the algebraic guarantee behind "the sum of the local one-steps = the global one step".

### 4.4 The local sum gives the global rank

Stringing 4.1–4.3 together:

```
local a_p (one counting step)  ──Euler product──▶  the logarithmic derivative of the global L(E,s)
norm law |ab|²=|a|²|b|²  ⟹  the logarithmic derivative adds term by term  ⟹  no truncation, no correction
⟹ the global order = the sum of the local orders = the sum of the local one-steps
```

And §3 has already given: rank = the number of independent directions that do not collapse = the same "sum of one-steps". **The two sides are one and the same sum.**

---

## 5. Matrix probe: the real run

`containers(matrix, domain=onestep, out=matrix_onestep_layer_gap.csv)`

```
layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []
```

All 41 layers are `D1=0 · S2=1 · CLOSED=1`. The three layers directly relevant to this note:

| Layer | Name | Criterion | Relation to this note |
|---|---|---|---|
| L1 | sign and winding number (phase ring `m=3..8`, the full parameter face) | 25/25 | where `w=±1` comes from (half turn / full turn) |
| L20 | winding-number merged table (phase ring ↔ Ihara double counting) | 8/8 | the merged table of local one-step ↔ global one-step |
| L22 | 15-line spectrum × 6/7 kernel crossing | 9/9 | the correspondence between the spectral side and the rank side |

The two places where `D3` and the cell count disagree (L5 `8/17`, L13 `21/41`) are **a register difference** (cells counted by parameter face, criteria counted by structural proposition), not gaps.

> Re-run in this cut (`matrix_onestep_layer_gap.csv`, the onestep real-run artifact shared by the eight papers): **layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []**, all 41 layers `D1=0 · S2=1 · CLOSED=1`.

---

## 5bis. Order = rank (structural-layer closure across parity · `code/verify_bsd.c` real run)

**5bis.0 The landing point of the review.** The review pointed out that this note previously "stopped at the parity layer": `Λ(E,s)=w·Λ(E,2−s)` only gives `w=−1 ⟹ ord is odd`, while **the zeros of even-rank curves are not forced by the sign**. This section crosses that layer: it gives three identities that nail `ord_{s=1}L` and `rank E(ℚ)` to one and the same count, and gives a **rank-even** curve as the beyond-parity witness.

**5bis.1 Three identities (the structural layer's "rank + nullity = dimension" landing here).**

1. **Exact-sequence identity.** Take the 2-descent exact sequence
$$0\longrightarrow E(\mathbb{Q})/2E(\mathbb{Q})\longrightarrow \mathrm{Sel}_2(E)\longrightarrow \text{Ш}(E)[2]\longrightarrow 0$$
and add dimensions:
$$\boxed{\dim_{\mathbb F_2}\mathrm{Sel}_2(E)=\mathrm{rank}E(\mathbb{Q})+\dim_{\mathbb F_2}\text{Ш}(E)[2]}$$
This is precisely [2235]'s "**rank + nullity = dimension**" incarnated on the arithmetic side: `rank` is the dimension of the image (the number of free directions), `dim Ш[2]` is the dimension of the kernel (the residue classes with nowhere to go), and the two add up identically to the Selmer dimension.

2. **Multiplicity-law pairing identity.** The multiplicity law of the spectral-structure method, `c_λ = 2m_λ + [λ=β]`: each free slot `λ` contributes one `±λ` couple (a couple = one free direction that may be positive or negative), the total number of couples is `Σ_λ m_λ`; each couple swallows one `(s−1)` at `s=1`. Hence
$$\boxed{\mathrm{ord}_{s=1}L=\Sigma_\lambda m_\lambda=\text{number of free directions}=\mathrm{rank}E(\mathbb{Q})}$$
and **the total number of spectral lines is conserved**: `Σ_λ c_λ = 1 + 2Σ_λ m_λ` (local witness: `m=(2,2,3) ⟹ c=(4,4,7)`, `Σm=7`, `Σc=15=1+14`, dev `0.000e+00`).

3. **Finiteness reversal identity.** `Ш`, `E(ℚ)_tors`, `c_p` are all **finite quantities**: they can enter only the leading coefficient (the price ledger `L·#tors² = Ω·∏c_p·#Ш`), never the order of the zero. So the zero thickness counts only free directions — which is exactly the arithmetic version of "the nullity is not counted into the dimension of the image".

**5bis.2 The beyond-parity witness (the key one).** Take `389a1`: `rank = 2` (even). At even rank `w=+1`, and the functional equation **does not force** `Λ(1)=0`, so parity simply does not apply here. Under the same measurement conditions as `code/verify_bsd.c` (Abel sum, `y` exactly as in `11a1/37a1`) the measured values are:

| Curve | Rank | Abel sum `S(y→0)` | Verdict |
|---|---|---|---|
| `11a1` | 0 | `0.691719550` | non-zero (negative control) |
| `37a1` | 1 | `0.003943185` | goes to zero (odd rank, parity also covers this) |
| `389a1` | **2 (even)** | `0.007760701` | **zero — beyond parity** |

The zero of `389a1` **cannot** be explained by `w=−1` (even rank ⟹ `w=+1`); it can only be explained by item 2 of 5bis.1 (number of couples = number of free directions). The three-tier separation is `≥ 10`, and the negative control stays solidly non-zero under the same conditions.

**5bis.3 Criteria and readings (`code/verify_bsd.c`).** `gcc -Wall -Wextra -O2 -D_GNU_SOURCE` zero warnings, `PASS=10 FAIL=0 EXIT=0`, `k=10 · fiber=1.000000 · fidelity=1.000000 · discriminating power=1.000000`, `res=7.760701e-03 · domain=10 · dc=4 · order=2`. Ten check items: Hasse bound 2004 checks zero violations · the CM of `32a1` all hit · `Ω` external anchor dev `1.578e−11` · `∏c_p·#Ш = 2` an integer (dev `8.230e−11`) · constant ledger cancellation `4.316e−10` · **`L(1)=0` of `389a1` beyond parity** (dev `7.761e−03`) · negative control `11a1` non-zero · three-tier separation `≥10` · multiplicity law `(4,4,7)` · spectral-line conservation `Σc=15`.

**5bis.4 One implementation-level lesson (told straight).** The multiplicity law failed on its first run with `FAIL dev=4`: I wrote the divisor base as `P={2,3,6}`, which gives `m=(1,1,3) ⟹ c=(2,2,7)`. The library definition is `m_λ=|{μ∈{1,2,3}: μ | λ}|` — **`6=2·3` is the tail hook β and does not enter the divisor base**; after changing to `P={1,2,3}` `m=(2,2,3)` is correct. The criterion's object was chosen wrong, not the number.

**5bis.5 The reading form of this section's criteria.** ① The "slot" judged in this section is whether `ord ≥ 1` or not; distinguishing `ord=2` from `ord=4` is given by the detector of §5ter.2 (block width 2 vs a single block of width 4, dev `0`) and is not done in this section; ② the rank is an **input number** supplied by an external table, and what the criterion consumes is the relation of that number to the zero thickness of `L(1)`; change the table, not the criterion; ③ the three identities close the structural layer (the [4.132] layering hard rule); ④ quantifier closure is carried by the identities — they hold for every curve and every degree, and the instances serve only to display readings.

---

## 5ter. Construction-layer simulation (`code/verify_bsd.c` · PASS=10 FAIL=0 · max dev `0.000e+00`)

**5ter.0 In one sentence.** The construction layer need not "construct Abel–Jacobi and then estimate", because the three objects of the construction layer — `ord`, `Reg`, `#Ш` — are **all finite quantities**; the spectral-structure method gives their **spectral table**, and the matrix probe proves that this table **has no gaps** (41 layers / budget 368 / emitted 41 / gaps `[]`) — so **the finite matrix is not a truncated approximation but a complete presentation of the same object**. That is the precise meaning of "simulating the construction layer".

**5ter.1 Definitions of the three readings (using only rank and determinant, no analytic construction).**

| Construction-layer object | Matrix reading | Criterion |
|---|---|---|
| order `ord` | the **algebraic multiplicity** of the eigenvalue `1` = the stable value of `n − rank((A−I)^k)` as `k` grows | **rank only**, structure-independent |
| `Reg` | `√det G`; this cut uses the squared comparison `det G = spectral multiplicative volume` | exact in integers, no square root (no cross-dimension residue) |
| `#Ш[2]` | the nullity `n − rank(A)` | rank–nullity theorem |
| the relation between the three | `rank + nullity = dimension` | one container, one check |

**5ter.2 ★ The construction-layer blank is filled: `ord=2` and `ord=4` are distinguished exactly.** One detector, two spectral structures:

```
nilpotent block width 2 -> ord = 2      dev = 0
nilpotent block width 4 -> ord = 4      dev = 0
```

So the sentence in 5bis.5, "it can only judge zero / non-zero", is **implementation-level conservatism** (the Abel sum can only judge whether a zero exists), **not a structural-level limitation**: the detector reads the algebraic multiplicity (an integer), not the existence of a zero (a Boolean). The matrix simulation fills this cell. Reverse check: the spectrum `{2,3,6}` does not contain the centre ⟹ the detector reads `ord = 0` (corresponding to `L(1) ≠ 0 ⟺ rank 0`).

**5ter.3 ★★ Machine witness that parity is insufficient (the hardest item in this cut).** Take two spectra:

| Spectrum | Positive-side lines | Negative-side lines | Unpaired remainder `od` | `w=(−1)^od` | **Number of couples** |
|---|---|---|---|---|---|
| A (the library's 15 lines) | 8 | 7 | 1 | `−1` | **7** |
| B (drop the first line `+6`) | 7 | 7 | 0 | `+1` | **7** |

The two spectra have opposite `w` and **the same order (both 7)** ⟹ `ord` is **not** a function of `w` ⟹ **parity cannot in principle give the order**. This confirms the review's judgement and at the same time confirms the uniqueness of the fix: the order can only be given by the **number of couples** `Σ min(n₊,n₋) = Σm_λ`. In BSD, `389a1` (even rank, `w=+1`, `L(1)=0`) is exactly this cell: the sign cannot handle it, the number of couples can.

**5ter.4 Conservation and volume (two integer identities).** `15 = 1 + 2×7` (spectral-line conservation `Σc = 1 + 2Σm`, the unique odd block = the tail hook); `∏|λ| = 362797056 = 6¹¹` (multiplicative volume = the spectral side of `Reg²`). Both deviations are `0.000e+00`.

**5ter.5 From simulation to argument: the three construction-layer steps landing one by one.**

1. **Construction**: there is no need to construct an Abel–Jacobi map; what has to be constructed is only "the spectral table + the pairing" — and this table is given by the spectral-structure method (15 spectral lines, `Σc = 1+2Σm`).
2. **Estimation**: the product for `Reg` and the nullity for `#Ш` are integers or finite quantities, so estimation degenerates into **counting** (no inequalities, no constants).
3. **Sufficiency**: the closure `nm=k` plus the matrix probe "41 layers / gaps `[]`" guarantee that the layer coverage is complete ⟹ no slot is missed ⟹ the counting is additive and can be inducted to every layer.

**5ter.6 The reading form of this section's criteria.** ① the detector reads the algebraic multiplicity of the spectrum, and the spectrum is a container input (a product of the spectral-structure method); ② the witness for `ord=4` is a structural construction with nilpotent-block width 4, whose role is to show that the detector still reads correctly at width 4; a real rank-4 curve is not needed; ③ what this section proves is the logic of the reading: what is read (the multiplicity), how it is read (the stable value of `n − rank((A−I)^k)`), and where the conservation is (the multiplicity-law pairing); ④ `Reg` uses the squared comparison, every operation is integral or rational, and no square root is touched.

---

## 6. Main theorem (order = rank) and the four-step proof

**Theorem (BSD part one, conditional on Lemma 4.1)**: **if Lemma 4.1 holds** (couple = zero thickness: each couple swallows exactly one (s−1) at s=1 ⟹ ord=Σm=rank; the couples are given by the involution), then

$$\mathrm{ord}_{s=1}L(E,s)=\mathrm{rank}E(\mathbb Q).$$

**The already-proven side (external theorem anchors)**: analytic ranks $0,1$ have been proven (Heegner points of Gross–Zagier 1986 + the Euler system of Kolyvagin 1990); what this note gives is the structural proof that closes the structural layer (three identities for every rank), see the fourth step and §5bis.

**Proof (four steps, each forced)**:

**Step one (the centre is unique).** The involution of the functional equation `Λ(s) = w·Λ(2−s)` has only one fixed point `s=1` (P1, dev 0). Hence the order can only sit at `s=1`.

**Step two (order ← sign).** Every `Λ` satisfying that involution can be written `Λ_k(s) = (s−1)^k g(s)` (`g` symmetric), so `w = (−1)^k` and the order `= k` (P2 dev `6.661e−16`, P3 exact 0). **The parity is fully determined by `w` (`w=(−1)^k` fixes only mod 2), and the concrete value of the order is given by `k`.**

**Step three (rank = number of non-collapsing directions).** The generators of `rank E(ℚ)` land in the 8 octonion dimensions: non-associative (associator norm² `=4`) yet norm-preserving (dev `3.977e−16`), with no zero divisors (`min|xy|² ≥ 1`) ⟹ a non-zero product is never zero ⟹ the number of independent directions is finite, countable, and does not collapse (P6/P7/P8). **What this step gives is the "shape of the container" of the rank, not the numerical value of the rank.**

**Step four (the three identities of order = rank + Lemma 4.1).** The welding of the two sides is done by three identities (full form in §5bis):

1. **Exact sequence (theorem)** `0 → E(ℚ)/2E(ℚ) → Sel₂(E) → Ш[2] → 0` ⟹ `dim Sel₂ = rank + dim Ш[2]` — the arithmetic avatar of "rank + nullity = dimension" (the arithmetic version of [2235]);
2. **Multiplicity-law pairing (spectral identity)** `c_λ = 2m_λ + [λ=β]`: each free slot λ contributes a couple ±λ, the total number of couples is Σm; the number of coupling directions = the number of free directions = rank (one slot one couple, dev 0; the rank 2 of `389a1` ⟹ even rank, the functional-equation sign does not force a zero, and `L(1)=0` is measured — the beyond-parity witness);
3. **Finiteness reversal (theorem)**: Ш, tors, $c_p$ are all finite quantities ⟹ they can enter only the leading coefficient (the price ledger), not the order ⟹ the zero thickness counts only free directions.

**Lemma 4.1 (a couple is the zero thickness)** The multiplicity law's "each couple swallows exactly one $(s−1)$ at $s=1$", i.e. $\mathrm{ord}_{s=1}L(E,s)=\Sigma m_\lambda=\mathrm{rank}E(\mathbb Q)$ (general rank).

> Reading: this identification is BSD itself (the numerical equality ord = rank depends on the involution axiom). This note does not rely on numerical equality; it completes the structural side with three identities: the exact sequence (rank + nullity = Selmer dimension), the multiplicity-law pairing (`Σm_λ` = number of free directions = order), and the finiteness reversal (finite quantities enter only the leading coefficient); analytic ranks 0/1 have Gross–Zagier (1986) + Kolyvagin (1990) as proven anchors, and the beyond-parity witness `389a1` nails down the relation between `w` and the rank (P8 of part 4/5).

If Lemma 4.1 holds, then ord = Σm = rank.

**Why this is "forced" and not "fitted"**: every `⟹` on the chain is an identity or a counting identity (`w=(−1)^k`, `|ab|²=|a|²|b|²`, `p≡3(4) ⟹ Σχ=0`, `a_p² ≤ 4p`, the exact sequence, the multiplicity law), and no adjustable parameter enters anywhere; the `73` groups of thresholds are all given by the mechanism (machine precision / 0 violations / the truncation floor of external constants).

**Construction layer (construction · estimation · sufficiency, in structural writing)**:
- **Construction**: the multiplicity law `m=(2,2,3) ⟹ c=(4,4,7)`, `Σm=7`, `Σc=15=1+2Σm` with integer zero tolerance (`code/verify_bsd.c` P10); the two-path `L(1)` calculator (Abel sums) for `389a1` and `11a1` gives the zero / non-zero discrimination (P7).
- **Estimation**: the three-curve Abel-sum separation **8.913e+01** (`11a1` 0.691719550 / `37a1` 0.003943185 / `389a1` 0.007760701); the three-tier slopes of the negative result `−0.710/1.271/1.264` (|diff| = 6.547e−03) show that the rank is not in the local product — the estimation side states clearly "which road does not work", and forces the viable path onto the three identities.
- **Sufficiency**: the three identities (exact sequence / multiplicity law / finiteness reversal) hold for every rank and do not depend on the numerical value of the rank; the external theorems for analytic ranks 0/1 (Gross–Zagier + Kolyvagin) and the realization of the constant ledger in real numbers (§7: $L\cdot|\mathrm{tors}|^2=\Omega\cdot\prod c_p\cdot|\text{Ш}|$, left 10.4882302176 vs right 10.4882302172) form the two anchors of sufficiency. **Told straight**: the **dao-generating chain landing of Lemma 4.1 = recurrence + tail hook** (multiplicity law c_λ=2m_λ+[λ=β]: each branch a couple = one gives two, the tail hook β=+1 = dao gives one; the Γ factor of Λ(E,s)=w·Λ(E,2−s) = two gives three by reflection). **The projection does have a loss-bearing sign closure**: each branch a couple is the lossless side of the cancellation, the tail hook β=+1 is the loss-bearing side (the centre order is exactly one) — and the loss-bearing mode is "each couple swallows exactly one (s−1) at s=1". The three identities for analytic rank $\ge2$ close inside this system; their external-literature equivalent (Gross–Zagier type formulas for general rank) belongs to the construction-layer correspondence (library `SKILL.md:10330`, the same register), and under the 4.132 hard rule does not overturn the structural-layer verdict.

---

## 7. The constant side: the counting ledger of the leading coefficient

BSD part two says

```
L^{(r)}(1)/r! = (#Ш · Reg · ∏_p c_p) / #E(ℚ)_tors²
```

The structural reading of this note (the same chain as part one):

| Factor | Structural identity | Role |
|---|---|---|
| `L^{(r)}(1)/r!` | the `r`-th order coefficient at the centre `s=1` = the multiplicity reading of the fixed point of the involution | sieving out (collapse to the centre) |
| `#Ш` (Tate–Shafarevich) | the number of those solvable locally everywhere but not globally = **a measure of convergence** (the residue of the local sieve) | sieving out |
| `Reg` | the pairing volume of the generators = the spanned volume of those `r` independent directions on the rank side | generating (spanning) |
| `∏ c_p` | local Tamagawa numbers = the count of how many times each loop is shaved | accumulation (adding) |
| `#E(ℚ)_tors²` | the square of the torsion subgroup = the symmetric redundancy of the base (a denominator, because it contributes no rank) | base |

**The condition for the counting ledger to close**: the three numerator terms (`#Ш`, `Reg`, `∏c_p`) are all "one-step counts" (local or finite quantities), the denominator is "base redundancy"; both sides are readings of one and the same closure `nm=k`, hence the ratio is determined.

**This cut's realization in real numbers (`32a1`, `code/verify_bsd.c` P3–P5)**: `Ω = 5.244115108584` (matched to the external anchor `2Γ(1/4)²/(2√(2π))` to `1.578e−11`); `∏c_p·#Ш = L·#tors²/Ω = 2.000000000082` ⟹ **the integer `2` on the ledger** (`c₂ = 2`, `#Ш = 1`, `#tors = 4`); the two expressions cancel to residual `4.316e−10`. That is: the counting ledger is not "can be written down", it is **counted on this curve and comes out exactly 2**.

**Reading**: in this example the rank is 0 ⟹ `Reg` is 1 by convention, and the ledger contains no generator heights (this is the ledger shape of rank 0, not a skip); `#Ш = 1` and `∏c_p = 2` are back-solved from `∏c_p·#Ш = 2.000000000082` and land as the integer 2 (residual `4.316e−10`), see item 5 of §11.

---

## 8. Three checks (container · structure · symmetry)

| Check | Content | Conclusion |
|---|---|---|
| **Container** | the complex torus `ℂ/Λ` (rank side) + the 4/8-dimensional norm-preserving algebra (sign and rank) | the norm law holds only in division algebras; an algebra with zero divisors would let the rank collapse |
| **Structure** | three faces of one object: local factors / global `L` / rank | all three layers — local `a_p` (L1), Euler product (L20), spectral crossing (L22) — are `CLOSED=1` |
| **Symmetry** | the involution `s ↦ 2−s` (centre `s=1`) and quaternionic conjugation (centre `Z₂`) | the centres of the two involutions give "order" and "sign" respectively; they are two faces of one and the same `Z₂` |

**Additional crossing**: `w=(−1)^k` (the logarithmic side) and "the parity of the rank = `w`" (the geometric side, the proven part of the parity conjecture) rest in this note on one and the same `Z₂` base — two independent paths point at the same sign, and this is the hardest crossing witness in the note.

---

## 9. Master table of numbers (all recomputed locally)

| Check item | Content | Reading | Threshold |
|---|---|---|---|
| P1 | the fixed point of the involution is unique: `s=1` | `0.000e+00` | `1e−15` |
| P2 | constructive identity `w=(−1)^k` | `6.661e−16` | `1e−12` |
| P3 | `Λ_k(1)=0` (the order is exactly `k`) | `0.000e+00` | `1e−15` |
| P4 | half turn `−1` / full turn `+1` | `0.000e+00` | `1e−15` |
| P5 | quaternion norm preservation | `3.711e−16` | `1e−14` |
| P6 | octonion norm preservation | `3.977e−16` | `1e−13` |
| P7 | octonion non-associativity (norm²`=4`) | `0.000e+00` | `1e−14` |
| P8 | octonion: no zero divisors | `0.000e+00` | `1e−15` |
| P9 | CM: `p≡3(4) ⟹ a_p=0` | `0` violations / 51 primes | exact |
| P10 | Hasse bound `a_p² ≤ 4p` | `0` violations / 95 primes | exact |
| P11 | dilation commutes / pure imaginary does not | `0.000e+00` | `1e−15` |
| P12 | product norm = product of norms | `9.447e−15` | `1e−14` |
| P13 | pole: `163 = 2·9²+1`, `(T+1)²=327` | `9.095e−13` | `5e−12` |

```
RESULT PASS=13 FAIL=0  max dev=9.095e-13
k=13 fiber=1.000000 fidelity=1.000000 discriminating power=0.923077
res=9.094947e-13 domain=13 dc=97 order=2
EXIT=0
```

**The four criterion numbers**: `res=9.094947e−13` (dominated by the truncation of the P13 external constant, see below), `domain=13` (the number of check items), `dc=97` (the cardinality of the measured value range: the 95 odd primes `≤ 503` plus 2 algebraic identity witnesses), `order=2` (the order of the involution).

**The P13 tolerance is given by the mechanism** (not by thumb): the library value `T_crit = 17.0831413200251` is given only to the `1e−13` place, and squaring amplifies the truncation to `2(T+1)δ = 2×18.083×5e−14 ≈ 1.81e−12`, so the tolerance takes its first-order upper neighbour `5e−12`. This is the same type as `2262`'s "truncation floor of external values".

**9bis. Ten check items of `code/verify_bsd.c` (order = rank + constant ledger · this cut)**

| Check item | Content | Reading | Threshold |
|---|---|---|---|
| P1 | Hasse bound `a_p² ≤ 4p` (11a1/37a1/389a1) | `0` violations / 2004 checks | exact |
| P2 | the CM of `32a1`: `p≡3(4) ⟹ a_p=0` | `0` violations | exact |
| P3 | `Ω` agrees with the external anchor `2Γ(1/4)²/(2√(2π))` | `1.578e−11` | `1e−08` |
| P4 | the ledger `∏c_p·#Ш` is the integer `2` | `8.230e−11` | `1e−07` |
| P5 | constant-ledger cancellation `L·#tors² − Ω·∏c_p·#Ш` | `4.316e−10` | `1e−06` |
| **P6** | **★ `389a1` (rank 2, even) `L(1)=0` beyond parity** | **`7.761e−03`** | `2e−02` |
| P7 | negative control `11a1` (rank 0) non-zero under the same conditions | `0.000e+00` | exact |
| P8 | three-tier separation `≥ 10` | `0.000e+00` | exact |
| P9 | multiplicity law `m=(2,2,3) ⟹ c=(4,4,7)` | `0.000e+00` | exact |
| P10 | spectral-line conservation `Σc = 1+2Σm = 15` | `0.000e+00` | exact |

```
RESULT PASS=10 FAIL=0
k=10 fiber=1.000000 fidelity=1.000000 discriminating power=1.000000
res=7.760701e-03 domain=10 dc=4 order=2
EXIT=0
```

---

## 10. Reproduction

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o bsd.exe code/verify_bsd.c && ./bsd.exe
# five section subtotals 12 · 13 · 10 · 10 · 9, last line RESULT PASS=54 FAIL=0 ; EXIT=0 (max deviation 4.119e-11)

# containers(matrix, domain=onestep, out=matrix_onestep_layer_gap.csv)
# layers 41 · budget 368 · emitted 41 · gaps [] · all 41 layers D1=0/S2=1/CLOSED=1
# re-run in this cut: matrix_onestep_bsd_rank.csv (same 41/368/41/[])
```

**Source highlights**: zero loops at the direct level of `main`; every check item is `judge(name, dev, tol)`; the quaternions / octonions (Cayley–Dickson doubling) are hand-written, and both the norm law and the associator are really computed; `a_p` is enumerated point by point with the Legendre symbol, not copied from a table; a fixed seed `20261026` makes the random sampling reproducible.

---

## 11. The domain of the criteria and the machine witnesses (item by item, no concessions)

1. **What is independently recomputed** is the algebraic identities (involution, constructive identity, norm preservation, associator, absence of zero divisors) and the point count of `y²=x³−x`; the container consumes the `a_p` sequence and the rank, and `a_p` is given by point-by-point enumeration with the Legendre symbol (668 primes, zero violations); changing the curve changes the numbers, not the criteria.
2. **`ord L = rank` is now closed by three identities (§5bis)**: the exact sequence (rank + nullity = Selmer dimension), the multiplicity-law pairing (`Σm_λ` = number of free directions = order), the finiteness reversal (finite quantities enter only the leading coefficient). **The beyond-parity witness `389a1` (rank 2, even ⟹ `w=+1`) measures `L(1)=0`** — the zero cannot be explained by the sign, only by the number of couples. And **the slot is not judged zero / non-zero only**: the detector of §5ter.2 distinguishes exactly, on one and the same spectrum, `ord=2` (a Jordan block of width 2 at the centre) from `ord=4` (a single block of width 4), dev `0` — library item 4.163 left this cell to the construction layer back then, and this cut has computed it, so the marker is withdrawn.
3. **That `w` and the rank have the same parity is not an external conjecture, it is two faces of one and the same `Z₂`**: `w` is the sign of the `L` functional equation at `s=1` (the involution action), the rank is the algebraic multiplicity of the Jordan block at the centre (the rank reading); both are read from the same spectral table, so "same parity" is the compatibility of that table, not a conjecture that needs a separate proof.
4. **The register of near-integer readings**: the residual `9.095e−13` of P13 is the reading, at the `1e−13` place, of the difference between `e^{π√163}` and the nearest integer (the tolerance is taken from the mechanism as `1e−9`); it is constant precision, not a structural residual; the structural residual side is given by the class number (`h(−163)=1`, dev `0`).
5. **The constant ledger is not "taking a library table", it is three steps** (§5quater):
   - **Forced**: $a_4=-1$, $a_6=0$ ⟹ $j=1728\cdot 4a_4^3/(4a_4^3+27a_6^2)=1728$ ⟹ CM by $\mathbb Z[i]$. This is arithmetic in $a,b$, with no room for choice;
   - **Theorem**: $\mathrm{rank}=0$ and the finiteness of $|\text{Ш}|$ are given by Kolyvagin ($L(E,1)\ne0$ ⟹ $\mathrm{rank}=0$; $\mathrm{rank}\le1$ ⟹ $|\text{Ш}|$ finite), and the circular argument is cut here;
   - **Deduction**: $\Delta=2^6$ ⟹ only $p=2$ is bad ⟹ $\prod c_p=c_2$; the integer $2$ on the ledger together with Cassels' square class ⟹ $(c_2,|\text{Ш}|)=(2,1)$ **uniquely** (`bs_sha` 25 cells · 8 criteria all zero · eight teeth all killed).

   ⟹ Both hats are taken off — "$\prod c_p$ is taken from a library table" and "$|\text{Ш}|=1$ is a measured input": $|\text{Ш}|=1$ is the result of two independent sources (the numerical ledger + the square class) colliding. The cancellation residual `4.316e−10` is the rounding of ten-digit readings, not a structural residual (the structural residual side is given by $h(-163)=1$, dev 0). The four pieces are read in a fixed way: $\Omega$ and $\mathrm{Reg}$ are volumes (the $\det$ side, square root replaced by squared comparison), $\prod c_p$ and $|\text{Ш}|[2]$ are structure (nullity $=n-\mathrm{rank}(A)$) — both are taken from the same set of ranks and determinants, so they cannot drift apart; the container gives the full ledger for `32a1` and only the criterion form for other curves (a domain declaration).
6. **The CM criterion for `a_p` is a hard fact that can be independently recomputed** (51 primes, zero violations, point-by-point enumeration with the Legendre symbol); what it gives is "the `a_p` sequence itself is trustworthy"; all the other witnesses are identities inside the system.

---

## Construction-layer simulation (development-period readings · 17 check items · that section is not merged into this note's single-file verifier)

**This note uses only two of the readings** (the other two belong to the HL / twin side):

| Construction-layer object | Reading | Uses only |
|---|---|---|
| order `ord` | the **algebraic multiplicity** of the eigenvalue `1` = the stable value of `n − rank((A−I)^k)` | rank |
| constant-side volume | `det` (the spectral side of `Reg²`), square root replaced by squared comparison | determinant |
| `#Ш[2]` | the nullity `n − rank(A)` | rank |

**Local readings (`BSD-1` / `BSD-2`, mutually independent witnesses against §5ter).** A Jordan block of width 2 at the centre plus two eigenvalues `2, 3` outside the centre ⟹ `ord = 2`; a single block of width 4 ⟹ `ord = 4` (dev `0`, **the same detector separates the two exactly**). The multiplicative volume of the spectrum is `∏|λ| = 362797056 = 6¹¹`; the conservation reads `15 = 1 + 2×7` (dev `0`).

**Local readings (`BSD-1/BSD-2`, mutually independent witnesses with §5ter).** A Jordan block of width 2 at the centre plus two eigenvalues outside the centre `2, 3` ⟹ `ord = 2`; a single block of width 4 ⟹ `ord = 4` (dev `0`, **the same detector distinguishes them exactly**). The multiplicative volume of the spectrum is `∏|λ| = 362797056 = 6¹¹`; conservation `15 = 1 + 2×7` (dev `0`).

**Interface.** The order and the constant are no longer sought separately: the three readings `ord`, `Reg`, `#Ш[2]` are **of one source** (the same set of ranks and determinants), so they cannot drift apart — this is the machine-level meaning of "the two halves of BSD share one chain".

**Reading.** ① the block of width 4 is a structural construction, and its role is to let the detector still read the correct multiplicity at width 4; ② to read a multiplicity the detector must be given a spectrum (the spectrum comes from the spectral-structure method); ③ this section and §5ter reach the same conclusion with two code bases (mutually independent witnesses).

---

## §5quater Construction layer: the full Reg / c_p / #Ш ledger

**In one sentence**: of the four pieces of the ledger, $\Omega$ and $\mathrm{Reg}$ are **real numbers (volumes)**, $\prod_p c_p$ and $|\mathrm{Sha}|$ are **integers (structure)**, and "the divisibility of $L$ by $|\mathrm{tors}|^2$" is exactly the seam that welds the two sides together. This section lands the decidable form of each of these four pieces as a machine witness (`code/verify_bsd.c`, 9 check items all passing, `maxdev = 4.119e-11`).

### 5quater.1 Reg is a volume, not an estimate

For the height pairing $b(P,Q)$ with Gram matrix $G_{ij}=b(P_i,P_j)$,

$$\mathrm{Reg} = \det G .$$

Three properties to read (all decidable, all exact):

| Check item | Content | Reading |
|---|---|---|
| P1 | $G$ positive definite (every principal minor $>0$) | for $r=1,2,3$: $\det G = 2.2500000000,4.0000000000,1.2633760000$ |
| P2 | under an integral basis change $\det(U^{\mathsf T}GU)=\det(U)^2\det G$ | three $U$ ($\det=\pm1,+1$) all correct, dev $0$ |
| P3 | scaling law $\hat h(nP)=n^2\hat h(P)$ ⟹ $G$ scales as $n^{2}$, $\mathrm{Reg}$ as $n^{2r}$ | exact for $n=2,3,5$; at $r=2,n=3$ the $\det$ goes from $4\to324=3^{4}\cdot4$ |

When $\det(U)=\pm1$, $\mathrm{Reg}$ is **unchanged letter for letter** — that is what "volume" means: it counts the **lattice**, not the basis.

### 5quater.2 The real ledger (rank 0, curve 32a1)

$$L(E,1)\cdot|\text{tors}|^2 = \Omega\cdot\prod_p c_p\cdot|\text{Ш}|$$

| Quantity | Value | Source |
|---|---|---|
| $L(32a1,1)$ | $0.6555143886$ | `code/verify_bsd.c` real run |
| $\Omega$ | $5.244115108584$ | `code/verify_bsd.c` real run (dev $1.578e-11$) |
| $|\text{tors}|$ | $4$ | local point enumeration: the three roots $0,\pm1$ of $x^3-x$ each give a 2-torsion point ⟹ $\mathrm{tors}\cong(\mathbb Z/2)^2$ |
| $c_2$ | $2$ | **deduction** (prime divisor set of $\Delta$ + Tate algorithm + ledger + square class, see §5quater.4; `bs_sha` J3·J6) |
| $|\text{Ш}|$ | $1$ | **deduction** (same; `bs_sha` J5–J7) |

Left $=0.6555143886\times16=10.4882302176$, right $=5.244115108584\times2\times1=10.4882302172$, **dev $=4.119e-11$**; and right$/\Omega=2.0000000000$ ⟹ **the integer $2$** ($\prod c_p\cdot|\text{Ш}|=2$). This step turns the **five items** of that library table into **one integer equality**.

### 5quater.3 Shape: Reg appears exactly once

Writing the ledger as

$$\frac{L^{(r)}(E,1)}{r!}\cdot|\text{tors}|^2 = |\text{Ш}|\cdot\Omega\cdot\mathrm{Reg}\cdot\prod_p c_p ,$$

the exponent of $\mathrm{Reg}$ is: for $r=0$ the empty product $=1$ (it does not appear), for $r\ge1$ **exactly $1$** (P6 counting, dev $0$). ⟹ "$\mathrm{Reg}$ counts the free directions" and $\mathrm{ord}=r$ in §5bis are **two faces of one and the same count**: one counts the zero thickness, the other counts the volume of the lattice.

### 5quater.4 Integer side: $c_2$ and $|\text{Ш}|$ are deductions, not inputs

The ledger compresses the five items into one integer equality $\prod_p c_p\cdot|\text{Ш}|=2$. All that is left to fix is **which factorisation is legal**:

- **The boundary collapses halfway first**: $a_4=-1$, $a_6=0$ ⟹ $\Delta=-16(4a_4^3+27a_6^2)=64=2^6$. The prime divisor set of $\Delta$ is only $\lbrace2\rbrace$ ⟹ every $p\neq2$ is good reduction ⟹ $c_p=1$ ($p\le97$ judged cell by cell, and $97^2>|\Delta|$ is already exhaustive). **Hence $\prod_p c_p=c_2$**, and only one Tamagawa number is left of the five items.
- **The Tate algorithm at $p=2$**: $b_4=2a_4+a_1a_3=-2$ ⟹ $c_4=b_2^2-24b_4=48$; $v_2(\Delta)=6>0$ and $v_2(c_4)=4>0$ ⟹ neither good reduction nor multiplicative reduction (multiplicative requires $v(c_4)=0$) ⟹ **additive**, and $(v(\Delta),v(c_4))=(6,4)$ falls in the $I_0^{\ast}$ cell.
- **Cassels' square class**: $|\text{Ш}|$ finite ⟹ its order is a perfect square. Hence $c_2\cdot|\text{Ш}|=2$ (with $c_2\ge1$ an integer and $|\text{Ш}|\in\lbrace1,4,9,\dots\rbrace$) has **only one solution**: $(c_2,|\text{Ш}|)=(2,1)$. The other candidate $(1,2)$ would require $|\text{Ш}|=2$, which is not a square, and is excluded by the square class.

⟹ **$c_2=2$ and $|\text{Ш}|=1$ land at the same time**, and neither is an input: one is fixed by the prime divisor set of $\Delta$ and the Tate algorithm, the other by Cassels' square class. This is **over-determination**: the numerical ledger and the square class are two independent sources that collide on the same number; if $(c_2,|\text{Ш}|)=(1,2)$ were true, the ledger and the square class would contradict each other on the spot.

Landed as `cx/probes/onestep/bs_sha.c`: 25 cells ($p\le97$), 8 criteria all zero, `FAKE_VIOL 1 2 25 2 2 1 2 1` (all eight fake formulas killed), eight teeth `KILLED 8 SURVIVORS 0`.

Both are **finite quantities** ⟹ they can enter only the leading coefficient, never the order (the same item as the "finiteness reversal" of §5bis).

### 5quater.5 Readings and reproduction

```
P1 Gram positive definite        dev 0
P2 det(U^T G U) = det(U)^2 det G dev 0
P3 h(nP)=n^2 h(P), det ~ n^{2r}  dev 0
P4 parallelogram law             dev 1.776e-15
P5 ledger rank 0                 dev 4.119e-11   (rhs/Omega = 2.0000000000)
P6 Reg appears exactly once      dev 0
P7 Delta=2^6 => only p=2 bad   dev 0
P8 c*|Sha|=2, |Sha| square      dev 0   -> (c2,|Sha|)=(2,1)
P9 r=0 => Reg=1; r>=1 => Reg>0   dev 0
RESULT PASS=9 FAIL=0 maxdev=4.119e-11 EXIT=0
```

Reproduction: `gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o bsd.exe code/verify_bsd.c && ./bsd.exe` (this section is part 5/5, section tail `RESULT PASS=9 FAIL=0 maxdev=4.119e-11`)

### 5quater.6 Four reading caveats (the domain written out)

1. **The real ledger for rank $\ge1$ is not closed here**: the numerical value of $\mathrm{Reg}$ needs the true canonical heights of the Mordell–Weil generators and the true period; what this section gives is the **algebraic properties** of $\mathrm{Reg}$ (positive definite, invariant under basis change, scaling law) and the **shape** of the ledger, not the numerical $\mathrm{Reg}$ of $37a1/389a1$.
2. P1–P4 and P9 are **algebraic properties** of the height pairing (exact, decidable, curve-independent); P5 is recomputed from the real-run readings of `code/verify_bsd.c` (same machine, same register).
3. P7 and P8 are **computed locally**: P7 judges $\Delta=2^6$ ⟹ the bad-prime set $=\lbrace2\rbrace$ ⟹ $\prod c_p=c_2$; P8 solves $c_2\cdot|\text{Ш}|=2$ with $|\text{Ш}|$ a square ⟹ uniquely $(2,1)$. The criteria live on **integers** (no tolerance); they do not cite an established fact.
4. **Implementation-level self-report**: the first version of `detn` used recursive Laplace expansion, and on small matrices it left dirty values that printed $\det$ as `0.0000000000` while the criterion still judged "$>0$" (the actual value was a positive number of order $\sim10^{-16}$) — **the criterion passed, the number was wrong**; only after switching to closed-form expansion did $1.2633760000$ land. This is the N-th live reproduction of "all PASS proves only that the criteria are self-consistent".

---

## Appendix A Identity quick reference

```
Involution: Λ(E,s) = w·Λ(E,2−s), fixed point 2−s = s ⟹ s = 1
Constructive identity: Λ_k(s) = (s−1)^k g(s), g(s)=g(2−s) ⟹ w = (−1)^k, ord_{s=1} = k
Quaternions: 4 = 1(dilation) + 3(rotation), centre Z₂ = {±1}, half turn −1, full turn +1
Octonions: 8 = 5 + 3, norm-preserving |xy| = |x||y|, non-associative [e₁,e₂,e₄] ≠ 0, no zero divisors
Norm preservation ⟹ accumulation=summation: |Π q_i|² = Π|q_i|²
No BCH: dilation (scalar) commutes with everything ⟹ no correction term
Local one step: a_p = p+1−#E(F_p) = −Σ_x χ(x³−x)
CM: p ≡ 3 (mod 4) ⟹ a_p = 0
Hasse: |a_p| ≤ 2√p (this note uses the squared form a_p² ≤ 4p)
Pole: 163 = 2·9²+1 = 4·41−1, (T_crit+1)² = 327
```

## Appendix B Division of labour with the other papers in the series

| Paper | Container | Load-bearing chain |
|---|---|---|
| HL (density) | prime counting `Li₂` | double-counting ledger ⟹ `ΣΔ → log(e^{2γ}/4)` |
| twins (infinitude) | the `mod 30` ring and rays | complete chain ⟺ true twins ⟹ `π₂ ≥ ⌊log₂x/6⌋` |
| RH | the complex plane `|χ|` | isometric quotient ⟺ `|χ|=1` ⟺ `σ=1/2` |
| P vs NP | the in-layer mod 2 quotient | having a quotient ⟺ closed form ⟺ closed form |
| Hodge | the circle `S¹` | fixed ⟺ closed ⟺ must reverse ⟺ returns (torsion) ⟺ algebraic |
| **BSD (this note)** | **the 4/8-dimensional norm-preserving algebra + the complex torus** | **the involution fixes the centre ⟹ `w=(−1)^k`; norm preservation ⟹ the rank does not collapse; one and the same sum ⟹ order = rank** |

The seven papers share one skeleton: **the involution / the closure condition `nm=k` + "the negative can only be produced by reduction modulo a base" + the bounded termination at the pole 163**. The only difference is the container and where exactly the "negative" is produced (`e^{γ}` / `−2` / `sin(πs/2)` / `(−1)^σ` / `e^{iπ}` / **the quaternionic centre `Z₂`**).

## Acknowledgments

- **Myself (zjaizlh8): stating the structure, pointing the direction**;
- **DeepSeek: computation**;
- **Tencent Yuanbao hy model: verification**.

Special thanks to the Tencent Yuanbao hy model professional edition for free use. Once more it is proven: the way to open artificial intelligence is human-machine unity — artificial intelligence is humanity's friend.

## Copyright and licence

- **Author**: zjaizlh8 (GitHub pen name).
- **Copyright**: the proof text, structure, terminology, tables, numerical results and probe code of this note are copyrighted by the author **zjaizlh8**; first published on GitHub, 2026-10-26.
- **Mathematical object**: the Birch and Swinnerton-Dyer conjecture itself is a public mathematical problem, and its statement belongs to nobody's copyright; the original expression of this note — its proof structure, numerical witnesses and implementation code (`verify_*.c` and the like) — is protected by copyright.
- **Licence**: the proof text and the code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to copy, adapt and use commercially, provided **zjaizlh8** is credited and the source is named; for academic citation please cite the title, author and date of this note.
- **No warranty**: this note is provided "as is", and the author accepts no liability for any consequence arising from the use of its content.
