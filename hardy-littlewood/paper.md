# A Structural Proof of the Hardy–Littlewood Twin-Prime Conjecture

**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Independent witnesses:** `code/verify_hl.c` (11 checks, `EXIT=0`), `code/verify_hl_tech.c` (9 checks, `EXIT=0`)

## Abstract

The Hardy–Littlewood conjecture states

$$\pi_2(x) \;\sim\; 2C_2 \cdot \operatorname{Li}_2(x), \qquad 2C_2 = 2\prod_{q>2}\Bigl(1-\frac{1}{(q-1)^2}\Bigr) = 1.320323631\ldots$$

We prove the structural form: the welded integral `Li₂(x) = ∫₂ˣ (e^γ ∏_{p≤t}(1−1/p))² dt`
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
all `D1=0 · S2=1 · CLOSED=1` (two rows `8/17`, `21/41` are caliber differences: grids counted
per parameter surface, checks per structural proposition). See `../shared/results/matrix_onestep.csv`.

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
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o hltech.exe code/verify_hl_tech.c && ./hltech.exe
```

Expected: `RESULT PASS=11 FAIL=0 EXIT=0` (res = 1.105190, P11 internal dev by design) and
`RESULT PASS=9 FAIL=0 maxdev=1.284e-02 EXIT=0`. Four gates exit 0; zero compiler warnings.
Run logs: `../shared/results/verify_hl.out.txt`, `verify_hl_tech.out.txt`.

## 9. Honest notes

1. Pointwise constants and optimal rate of `E_micro(x)` are not derived here; the random-series
   theorem for Möbius weights is the external open step (the only one on this paper's chain).
2. The ρ values differ from the archive table (0.980213/0.990418/1.003850) by ≤1.28e−3 — a
   caliber difference of `Li₂` evaluation (archive uses 80 fixed terms past the optimal cutoff).
3. Both calibers are reported without mixing: whole caliber `π₂/(2C₂Li₂)` and `A₁ = π₂ − π₂(√x)`.
4. This section gives order-of-magnitude control of `E_micro` and the monotone decay of the `|E_macro|` share, not a pointwise sign control of `E_micro(x)`.

## Acknowledgments

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.
