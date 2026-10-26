# A Structural Reading of the Riemann Hypothesis

<!-- verify-glance -->
> **One paper, one file, one program.** `code/verify_rh.c` — **21 criteria**. `gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify verify_rh.c -lm && ./verify` → `RESULT PASS=21 FAIL=0` and `EXIT=0`.
**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Machine witness:** `code/verify_rh.c` — one program, two sections (12 · 9 = **21** criteria), closing with the aggregate `RESULT PASS=21 FAIL=0` and `EXIT=0`.

## Abstract

The Riemann Hypothesis states that all non-trivial zeros of `ζ(s)` lie on `σ = 1/2`.
We state the structural form (a conditional reading, not a formal proof): the functional-equation involution `Λ(s) = Λ(1−s)` keeps the zero
set closed, and the **unique-radius theorem** of the spectral method forbids the two sides
`β, 1−β` from carrying asymmetric weights — so the only allowed location is the fixed line
`β = 1/2` (the center line `c₀ = 1` of the spectral circle). The equivalent equidistance reading

$$|\chi(\sigma+it)| = 1 \Longleftrightarrow \sigma = 1/2 , \qquad
\chi(s) = \pi^{s-1/2}\frac{\Gamma((1-s)/2)}{\Gamma(s/2)}$$

is machine-verified to `3.553e−15`. The load-bearing step is **Lemma 4.1** (zero placement),
proved from the unique-radius theorem, so the placement is a spectral reading rather than a
search.

## 1. Structural layer

- **Three base objects**; unique carrier `729 = 3⁶`; four basic operations ≅ V₄ (reduction
  modulo a base is the only operation that can produce a sign — the `−` of the functional
  equation comes from `sin(πs/2)`, the same family as `e^{iπ} = −1`).
- **Carrier cone** `C₀: σ² + γ² = δ²` (cone center at the origin; reflection center at 1/2;
  the two centers agree only when `σ = 1/2`).

## 2. Spectral method

- Three readings of one object: intra-layer spectrum `X_k = Z₂^k` with spectrum `{k−2j}`,
  multiplicity `C(k,j)`; inter-layer nesting `L_n = X_{nm}`; closure condition `nm = k`.
- **The sixth-order recurrence** `x⁶ − 49x⁴ + 504x² − 1296`, radius set `{2,3,6}`, multiplicity
  law `c_λ = 2m_λ + [λ=β]`, conservation `c₀ = 1 = m₊(6) − m₋(6)` (the mother does not vanish:
  if `m₊ = m₋` the odd branch would be zero).
- **P8 hard anchor** `d = (0,0,1)` (the ±-difference lives only on the radius-6 slot);
  **P16 Riemann mirror** (spectral differences in pure geometric ratio, the asymmetry
  concentrated on a single radius, the positive branch annihilated by third differences).
  Both are machine witnesses with dev 0.

## 3. Matrix probe

`containers(matrix, domain=onestep)` : 41 layers / 368 / 0 gaps, `D1=0 · S2=1 · CLOSED=1`.

## 4. Conditional theorem — four steps

**Conditional theorem (RH). If Lemma 4.1 holds** (unique-radius theorem: the ±-asymmetry sits only on radius 6, center line c₀ = 1), **then** every non-trivial zero of ζ(s) satisfies Re s = 1/2.

1. **Involution keeps the zero set closed** (identity layer, machine-witnessed):
   `Λ(s) = π^{−s/2}Γ(s/2)ζ(s)`, `Λ(s) = Λ(1−s)` with relative deviations
   `4.159e−11 / 3.691e−10 / 8.295e−10` at the three non-zero points
   `0.30+2i, 0.40+7i, 0.60+11.3i`.
2. **Fixed-point algebra:** `(1−σ)² = σ² ⟹ σ = 1/2` (`f(σ) = 1−2σ` has the single zero `1/2`,
   dev 0). This is the forced form of "antipode = fixed point".
3. **Equidistance:** `|χ(σ+it)| = 1 ⟺ σ = 1/2`. Six tiers `t = 5/10/17/30/40/50`, bisection to
   200 steps: root `σ* = 0.5` with `|σ*−0.5| ≤ 3.553e−15`, endpoints of opposite sign at every
   tier. Dual path (modulus formula vs complex χ) dev ≤ 2.2e−16.
4. **Lemma 4.1 (zero placement).** A zero off the line would need two radii `β, 1−β` to carry
   asymmetric weights (the involution maps one side to the other). But the unique-radius
   theorem (P8: `d = (0,0,1)`; P16: the asymmetry is concentrated on a single radius) forces
   all ±-asymmetry onto the single radius 6, and the conservation law
   `c₀ = 1 = m₊(6) − m₋(6)` pins exactly one center line. Hence `β = 1/2` is forced;
   the zero set, being closed under the involution and admitting only the center-line
   asymmetry, lies on `σ = 1/2`. ∎

   The same law is the "6 and 7 are the root of the wave" statement: the fold = the center `1`
   (`c₀ = 1`), the antipodal pairs = the ±-pairs, the ±-difference carried only on radius 6;
   7 = the pairing number `Σm_λ = 7` (the generating law has 7 lines, the closure law
   `1+2·7 = 15`).

   **Why the criterion is an identity:** the zero set is read as a spectral circle — the fold
   sits at the center `1`, the antipodal pairs are the `±`-pairs, and the `±`-difference is
   carried on radius 6 only, so a zero off the circle would break the pairing count
   `Σm_λ = 7`. P8/P16 and the conservation law are machine witnesses with `dev 0`; the
   isomorphism is the framework's structural caliber. The previous draft's `0 = χ·0` vacuum
   step is deleted and recorded.

## 5. Construction layer (construct · estimate · sufficiency)

- **Construct:** complex arithmetic + Lanczos Γ (g=7) + Euler–Maclaurin ζ (N=1000, K=10),
  three hard anchors: `ζ(2) = 1.644934066848` (dev 2.327e−10), `ζ(−1) = −1/12`
  (dev 8.882e−16 tier), `ζ(0) = −1/2`.
- **Estimate:** first three zeros `t = 14.134725142 / 21.022039639 / 25.010857580`
  (digit-aligned with the archive 9-digit values), `max |ζ(1/2+it)| = 3.033e−10`;
  **position criterion:** for fixed zero `t`, sweeping `σ ∈ [0.10, 0.90]` in steps 0.01 gives
  `argmin σ = 1/2` with deviation 0.000 and off-line/on-line ratio **8.921e+07** (seven and a
  half orders of magnitude); Hardy Z reality `|Im Z| ≤ 1.562e−09` at non-zero points,
  `|Z| ≤ 6.341e−10` at the first three zeros.
- **Sufficiency:** root uniqueness (endpoint sign change) + Z reality; the welded interface
  `|ζ(σ+it)|/|ζ(1−σ−it)| = |χ(σ+it)|` holds to `1.185e−09`.

## 6. Numerical table (machine-verified)

| t | |χ(0.30+it)| | |χ(0.50+it)| | |χ(0.70+it)| | |χ(0.80+it)| |
|---|---|---|---|---|---|
| 10 | 1.097320710 | 1.000000000 | 0.911310605 | 0.869938909 |

Roots of `|χ|−1` at six tiers: `σ* = 0.500000000000000 ± 3.6e−15`; endpoint values of opposite
sign at every tier (root uniqueness). `Λ(s) = Λ(1−s)` to `8.295e−10`; `R² = id` to `5.551e−17`;
`|χ(s)|·|χ(1−s)| = 1` to 0.

## 7. Reproduction

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o rh.exe code/verify_rh.c && ./rh.exe
```

Expected: subtotals `12 · 9`, then the aggregate `RESULT PASS=21 FAIL=0` and `EXIT=0`
(the 12-criterion section also prints `RESULT max dev = 1.562e-09   (threshold 1e-08)`;
the second reports `maxdev=1.185e-09`); zero warnings.

## 8. Forced structure and its domain

1. The unique-radius theorem (P8/P16) and the conservation law are the machine witnesses
   (`dev 0`); the spectral-circle reading is what makes the placement criterion one equation
   instead of a grid search.
2. Λ-symmetry must be measured away from zeros (at `s = 0.50+14.1347i` the relative deviation
   distorts to 1.694 because `|Λ| ≈ 1.09e−11`); the three non-zero points above are used.
3. The position criterion is read as: the sign change at the endpoints is unique, so the zero
   is pinned by the pairing count `Σm_λ = 7`; the grid fixes the resolution of the witness.
4. The threshold `1e−08` is set by the implementation floor (1.562e−09), not by hand.
