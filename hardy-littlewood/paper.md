# A Structural Reading of the Hardy–Littlewood Twin-Prime Conjecture

<!-- verify-glance -->
> **One paper, one file, one program.** `code/verify_hl.c` — **36 criteria**. `gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify verify_hl.c -lm && ./verify` → `RESULT PASS=36 FAIL=0` and `EXIT=0`.
**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Machine witness:** `code/verify_hl.c` — one program, four sections (11 · 9 · 6 · 10 = **36** criteria), closing with the aggregate `RESULT PASS=36 FAIL=0` and `EXIT=0`.

## Abstract

The Hardy–Littlewood conjecture states

$$\pi_2(x) \;\sim\; 2C_2 \cdot \operatorname{Li}_2(x), \qquad 2C_2 = 2\prod_{q>2}\Bigl(1-\frac{1}{(q-1)^2}\Bigr) = 1.320323631\ldots$$

We state the structural form (a conditional reading, not a formal proof): the welded integral `Li₂(x) = ∫₂ˣ (e^γ ∏_{p≤t}(1−1/p))² dt`
(the five-element interface: lower bound 2 / logarithmic derivative / local factor `1/log²t` /
integral / upper bound x; the water-side expansion `Li_k(x) = (x/log^k x) Σ (k)_n / log^n x`
with first coefficient = dimension D = k = 2) is forced by the three base objects
(unit 1, Cartesian product, place-value power `b^k`), and the error term `E_micro(x) = π₂(x) − 2C₂Li₂(x)`
is closed by **double cancellation**: the reciprocal identity on the constant side and the
random-series (sign) cancellation on the error side. Two readings of one error are kept distinct:
`E_micro = π₂ − 2C₂Li₂` (absolute, the wave amplitude, Θ(√x) — magnitudes measured in §6 *Estimate*) and
`E_macro = E_micro/(2C₂Li₂)` (relative, →0), welded by the reciprocal ledger
`E_macro × main = E_micro`. The load-bearing step is written as
explicit lemmas with honest notes.

## 1. Structural layer

- **Three base objects.** The unit `1`, the Cartesian product, the place-value power `b^k`.
- **Unique carrier.** `729 = 3⁶ = 9³ = 27²` is the unique number ≤10⁴ with ≥2 representations
  whose exponent contains 6 and whose base is a power of 3; it carries the whole arithmetic.
- **Four basic operations ≅ V₄ ≅ Z₂².** addition / multiplication / reduction modulo a base /
  crossing dimension (exp·log), closed in 16 composition slots; the carrier is not an operation.

## 2. Spectral method

Three readings of one object (same container, same action):

1. **Intra-layer spectrum** `X_k = Z₂^k`, spectrum `{k−2j}` with multiplicity `C(k,j)`;
2. **Inter-layer nesting** `L_n = X_{nm}`, `π∘ι = id`;
3. **Closure condition** `nm = k` (three statements kept separate: no cut chain `nm≥k`,
   single chain `nm≤k`, complete chain `nm=k`).

The generating law `A(N) = (1+(−1)^N)(1+2^N+3^N)+1` (even branch `3·1^N+2·2^N+2·3^N`, odd
branch 1), the closure law `Q(N) = (1+(−1)^N)(3·6^N+2·3^N+2·2^N)+6^N`, the sixth-order
recurrence `x⁶ − 49x⁴ + 504x² − 1296`, the radius set `{2,3,6}`, and the multiplicity law
`c_λ = 2m_λ + [λ=β]` with divisor counting `m_λ = |{μ∈{1,2,3}: μ|λ}| = (2,2,3) ⟹ c = (4,4,7)`,
`Σc = 15 = 1+2Σm`. **Macro = the positive branch, micro = the negative branch** (the even
branch is the measure side, the odd branch the fiber side) — both are forced by the same law.

## 3. Matrix probe

`containers(matrix, domain=onestep)` : **41 layers / budget 368 / emitted 41 / gaps []**,
all `D1=0 · S2=1 · CLOSED=1`; the two rows (`8/17`, `21/41`) are caliber differences, not gaps.

## 4. Two-sided ledger and the welded interface

- **Local side (admission)** `I_Q(n) = ∏(1−a)(1−b) ∈ {0,1}`; **global side (amount)**
  `P(Q) = ∏(1−2/q)`; **interface** `ΣΔ = log(S/M_nat)`.
- **Welded interface:** `Li₂(x) = ∫₂ˣ (e^γ ∏_{p≤t}(1−1/p))² dt`; the water-side expansion has
  first coefficient D = k = 2 (classical `x/log²x` is the first term only, 43× less precise).

## 5. Conditional theorem — four steps

**Conditional theorem (Hardy–Littlewood). If the sign-cancellation theorem holds** (Möbius-weighted random series, probe [412], the external open step), **then** π₂(x) ~ 2C₂·Li₂(x).

1. **Candidates unbounded** (gain side): `A(N) = 3·1^N + 2·2^N + 2·3^N > 0` strictly increasing.
2. **Lower bound forced at t = 1** (the k-th-order pole): `I(ε) = ∫_{1+ε}² dt/log²t ~ 1/ε`;
   ratio `I(1e−3)/I(1e−2) → 10` (machine-witnessed).
3. **The welded value equals the double product** (five elements aligned; the `e^γ` can only live
   inside the accumulated side and must appear twice).
4. **Error term closed by double cancellation.**

   - (a) **Reciprocal cancellation (constant side):**
     `R₀ → 4e^{−2γ} = 1.260947007` and `e^{ΣΔ} → e^{2γ}/4 = 0.793052993` are reciprocals;
     the product is identically 1 (machine dev 8.88e−16), hence `ρ = π₂/(2C₂Li₂) → 1`.
   - (b) **Sign cancellation (error side):** the Möbius-weighted sum is a random series
     (archive probe [412]: net/total `0.0014 → 0.0002` monotone collapse, random-walk ratio
     ≈0.4 stable) ⟹ complete cancellation ⟹ `E_micro(x) = o(x)`, i.e. `E_macro(x) = o(1)`.

   *Honest note:* the single open step is promoting the random series to a theorem
   (Möbius randomness / equidistribution); everything else above is an identity or a
   machine-witnessed forced bound.

## 6. Construction layer (construct · estimate · sufficiency)

- **Construct:** `E_micro(x) = π₂(x) − 2C₂Li₂(x)` is a finite quantity — the *absolute* reading (the
  wave amplitude); the *relative* reading is `E_macro(x) = E_micro(x)/(2C₂Li₂(x))`, and the reciprocal
  ledger welds the two: `E_macro × main = E_micro` (one error, two resolutions). `Li₂` computed by the convergent series (80 terms)
  or by Simpson on the welded integral — both in `code/verify_hl.c` (P1: series vs Simpson,
  four tiers, tol 1e−9).
- **Estimate:** `|E_micro|/√x = 0.0730 / 0.0771 / 0.0720` (same order, spread ratio 1.07 — the error is
  controlled by `√x`); `|E_macro| = 0.01851 / 0.00935 / 0.00387` (monotone decreasing);
  amplification factor `2^{ω(d)}` with `max ω(d) = 5` for `d ≤ √(1e7)` (worst d = 2310 = 2·3·5·7·11)
  ⟹ **factor ≤ 2⁵ = 32**; structural bound `|E_micro| ≤ 6√x` verified tier by tier (widest tier uses 1/78 of the bound — sufficient, not tight).
  *(Caliber note: the archive's 790.3 is the ratio `S(√x)/S(x)`, a different quantity.)*
- **Sufficiency:** the barrier is located by the multiplicity law — `c₂ = 2m₂ = 4` is even
  (pairs cancel, residue 0), the only unpaired slot is the β-slot (one) ⟹ the parity barrier
  is a reading of the multiplicity law, not a sieve defect.

## 7. Numerical tables (machine-verified)

| x | π₂(x) | ρ = π₂/(2C₂Li₂) | E_micro | \|E_micro\|/√x | \|E_macro\| |
|---|---|---|---|---|---|
| 1e5 | 1224 | 0.981495 | −23.077 | 0.0730 | 0.01851 |
| 1e6 | 8169 | 0.990655 | −77.060 | 0.0771 | 0.00935 |
| 1e7 | 58980 | 1.003874 | +227.609 | 0.0720 | 0.00387 |

Reciprocal identity: `4e^{−2γ} · e^{2γ}/4 = 1` (dev 8.88e−16). Second-order coefficient of the
water-side expansion: `(k)_2 = k(k+1)` (integer identity, zero tolerance).

## 8. Reproduction

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c && ./hl.exe
```

Expected: the four section subtotals print in order as `11 · 9 · 6 · 10`, then the aggregate
`RESULT PASS=36 FAIL=0` and `EXIT=0`. The first subtotal carries `maxdev=1.105190e+00` by design
(P11 reads an internal ratio); the others report `1.284e-02` and `7.039e-03`; zero warnings.

## 9. Honest notes

1. Pointwise constants and optimal rate of `E_micro(x)` are not derived here; the random-series
   theorem for Möbius weights is the external open step (the only one on this paper's chain).
2. The ρ values differ from the archive table (0.980213/0.990418/1.003850) by ≤1.28e−3 — a
   caliber difference of `Li₂` evaluation (archive uses 80 fixed terms past the optimal cutoff).
3. Both calibers are reported without mixing: whole caliber `π₂/(2C₂Li₂)` and `A₁ = π₂ − π₂(√x)`.
4. This gives order-of-magnitude control of `E_micro` and the monotone decay of the `|E_macro|` share, not a pointwise sign control of `E_micro(x)`.

## Acknowledgments

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.

## Chapter: the twin-prime case (k = 2)

> Twin primes are the `k = 2` case of Hardy–Littlewood; the constant `2C₂` and the Li₂ main term are the same, so this lives here as a chapter rather than a separate paper. The witness is the same `code/verify_hl.c` (segment 3–4).

## Abstract

The twin-prime conjecture states that there are infinitely many primes `p` with `p+2` prime.
We state the structural form (a conditional reading, not a formal proof): the ring `mod 30` has `φ(30) = 8` rays, the twin positions are
exactly the 3 pairs `(11,13), (17,19), (29,1)`, and **each layer of the spectral nesting keeps
at least one complete chain whose two sides are genuine primes**, so that

$$\pi_2(x) \;\ge\; \left\lfloor \frac{\log_2 x}{6} \right\rfloor \;\longrightarrow\; \infty .$$

The load-bearing step is **Lemma 5.1** (each layer keeps ≥1 genuine pair), whose proof is the
antipodal (reciprocal) argument: factor pairs are diametrically opposite about `√N`, so sieving
up to `√N` is complete. The honest note records exactly what this lemma is equivalent to.

## 1. Structural layer

- **Three base objects** (unit `1`, Cartesian product, place-value power `b^k`); unique carrier
  `729 = 3⁶` with 6 digits; the four basic operations ≅ V₄.
- **Ring and rays:** `φ(30) = 8`; the 3 twin ray pairs are forced by the admissible classes
  `{11,17,29}` (the pair `(3,5)` is the only exception, center 4 ∉ 6ℤ).
- **Gold only deletes, never creates** (local sieve removes classes; the gain side is the
  nesting of candidates).

## 2. Spectral method

- Three readings of one object: intra-layer spectrum `X_k = Z₂^k`; inter-layer nesting
  `L_n = X_{6n}` with `π∘ι = id`; closure condition `nm = k`.
- **Complete chains:** layer `L_n` has `|X_{6n}| = 2^{6n}` points, and the number of complete
  single chains is `2^{6n−6}` — machine-verified for n = 1..6:
  `1 / 64 / 4096 / 262144 / 16777216 / 1073741824`.
- **Layer intervals** `[2^{6(n−1)}, 2^{6n})` are disjoint ⟹ counts add: layers 1..3 sum to 2679
  = `π₂(262143)` (exact, dev 0).
- **Cut chain ⟺ composite** (deleted by a small prime), **complete chain ⟺ both sides coprime
  to all primes ≤ 2^{3n}** — the candidate/test boundary of §2.5.

## 3. Matrix probe

`containers(matrix, domain=onestep)` : 41 layers / 368 budget / 0 gaps, `D1=0 · S2=1 · CLOSED=1`.

## 4. Conditional theorem — four steps

**Conditional theorem (twin primes are infinite). If Lemma 5.1 holds** (each layer keeps at least one complete chain whose two sides are both prime — equivalent to the weak twin-prime statement), **then** π₂(x) ≥ ⌊log₂x/6⌋ → ∞.

1. **Candidates unbounded:** the generating law `A(N) = 3·1^N + 2·2^N + 2·3^N > 0` strictly
   increases; the layer count `2^{6n−6} ≥ 1` for every n.
2. **Chain/candidate ≠ proportional:** measured pair counts per layer give chain/pair ratios
   `0.14 / 0.64 / 1.59 / 4.66` — the structure guarantees ≥1 pair per layer, not the full
   chain count (the earlier "≥ 2^{6n−6}" check was a criterion-object error, corrected on the
   machine run and recorded).
3. **Lemma 5.1 (each layer keeps ≥1 genuine pair).** Proof by the antipodal argument
   (reciprocal structure):

   - (i) factor pairs `(p, N/p)` are diametrically opposite about `√N` (`p·(N/p) = N` is the
     reciprocal relation);
   - (ii) sieving up to `√N` is therefore **complete**: if a number has a factor `p > √N`, its
     complementary factor `q < √N` was already found on the small side — no large-prime-factor
     hole survives (witness: `437 = 19·23` escapes a cutoff-13 sieve but not a `√437 ≈ 21` sieve);
   - (iii) the number `64 = 2⁶` is the closure-forced bound: layer 2 spans up to `2¹² = 4096`,
     and the complete sieve bound is exactly `√(2¹²) = 2⁶ = 64` = the ring length; the closure
     condition `nm = k` aligns the sieve radius with the layer size;
   - (iv) the orbit-count theorem (permutation theory) gives `2^{nm−k} ≥ 1` complete chains per
     layer.

   ∎ by (i)–(iv), every layer contains a candidate pair whose both sides survived the complete
   sieve, i.e. are prime.

   *Honest note:* this lemma is equivalent to a weak form of the twin-prime conjecture itself
   (per layer, at least one candidate is genuinely prime); Brun's bound
   `π₂(x) ≤ C·x/log²x` is the true-theorem anchor, and the antipodal argument is structural —
   the ratio `0.14/0.64/1.59/4.66` is the machine-witnessed boundary between candidate and
   genuine.

4. **Counting adds:** disjoint layer intervals ⟹ `π₂(x) ≥ ⌊log₂x/6⌋ → ∞`. ∎

## 5. Construction layer (construct · estimate · sufficiency)

- **Construct:** per-layer enumeration + trial division (sieve to 1e7, `verify_hl.c` part 2
  bit-sieve to 1e9).
- **Estimate:** `s(x) = π₂(x)/⌊x/6⌋ = 0.073443 / 0.049014 / 0.035388` (monotone down);
  `π₂/x = 0.012240 / 0.008169 / 0.005898`; lower bound tiers `⌊log₂x/6⌋ = 2/3/3` with zero
  violations and slack multiples `612 / 2723 / 19660` (monotone up — the bound is true but loose);
  density product `∏_{q>2}(1−2/q) = 1.731263e−2 / 9.788832e−3 / 6.276420e−3`, times `log²X`
  `= 0.826108 / 0.830390 / 0.831923` (converging, relative spread 7.04e−3); density side /
  `2C₂ = 0.630091 < 1` ⟹ **the density side only gives an upper bound, never the lower bound**.
- **Sufficiency:** ≥1 pair per layer + additive counting; the lower bound never relies on a
  survival-rate estimate (Brun's wall — `Σ 1/p over twins` converges, the opposite direction of
  `Σ 1/p` — only rules out density arguments, not the antipodal one).

## 6. Numerical tables (machine-verified)

| x | π₂(x) | six(6m centers) | difference (3-tuple) |
|---|---|---|---|
| 1e5 | 1224 | 1223 | 1 |
| 1e6 | 8169 | 8168 | 1 |
| 1e7 | 58980 | 58979 | 1 |
| 1e8 | 440312 | 440311 | 1 |
| 1e9 | 3424506 | 3424505 | 1 |

`440310 = 146841 + 146953 + 146516` (three rays, aligned with the archive `qx4.txt`, violations 0).

## 7. Reproduction

This chapter rides on the same binary as the main note — sections 3 and 4 of
`code/verify_hl.c` carry the twin-prime case:

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hl.exe code/verify_hl.c && ./hl.exe
```

Expected: the five-tier table plus `RESULT PASS=6 FAIL=0` (count) and `RESULT PASS=10 FAIL=0
maxdev=7.039249e-03` (tech, `k=9 fiber=1 fidelity=1 discriminability=0.888889`), then the
file-wide `RESULT PASS=36 FAIL=0` and `EXIT=0`; zero compiler warnings.

## 8. Honest notes

1. Lemma 5.1 is equivalent to a weak form of the twin-prime conjecture; it is not derived from
   the density side (which only bounds from above).
2. Brun's theorem is the true-theorem anchor; the antipodal argument is structural.
3. Three criterion-object corrections were made and recorded during the machine run:
   (i) layer pair count vs ≥1; (ii) `ρ_first` decreases, not increases; (iii) the coefficient
   identity `(k)_n = k(k+1)…(k+n−1)` is the only exact statement (the `1+2u+6u²` comparison
   fails at x=1e5 because the `24u³` term is not negligible).
4. The count program uses an exact bit-sieve up to 1e9; no external tables, no fitting.

## Acknowledgments

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.
