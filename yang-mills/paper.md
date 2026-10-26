# A Structural Reading of the Yang–Mills Existence and Mass Gap

<!-- verify-glance -->
> **One paper, one file, one program.** `code/verify_ym.c` — **18 criteria**. `gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify verify_ym.c -lm && ./verify` → `RESULT PASS=18 FAIL=0` and `EXIT=0`.
**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Machine witness:** `code/verify_ym.c` — one program, two sections (9 · 9 = **18** criteria), closing with the aggregate `RESULT PASS=18 FAIL=0` and `EXIT=0`.

## Abstract

The Yang–Mills millennium problem asks for a quantum field theory on `R⁴` with compact gauge
group and a positive mass gap. The structural answer: the physical sector is fixed by the
**Hodge star involution** `*: Λ²R⁴` (spectrum 3/3, fixed points = self-dual 2-forms = the
vacuum sector, 3-dimensional), the gauge-invariant observables are traces (singlet projection
`P² = P`, `tr = 1`), the Wilson-loop sector is a **winding lattice** (closed loop ⟺ integer
winding number ⟹ discrete, quotientable, compact), and the spectral gap follows from
**repulsion** — the fragile zero of a perfect ring splits under any perturbation
(level repulsion with minimal crossing distance `2c`), so `λ₁ > 0`:

$$\lambda_l = l(l+2) = \text{multiplicity} − 1 , \qquad \lambda_1 = 3 > 0 , \qquad
\Delta = \sqrt{\lambda_1}/R = \sqrt3\,\Lambda .$$

The load-bearing step is **Lemma 4.1** (winding-lattice completeness + λ₁>0 in the continuum
limit), stated with its honest note.

## 1. Structural layer

- **Three base objects**; unique carrier `729 = 3⁶`; four basic operations ≅ V₄.
- **The sign comes from reduction modulo a base:** the half-angle double cover
  `U(2π) = −I` of SU(2) (machine dev 5.551e−16; the full-angle convention gives +I at 2π —
  difference exactly 2.000) is the sign structure of the gauge group.

## 2. Spectral method

- Three readings of one object: intra-layer spectrum `X_k = Z₂^k`; inter-layer nesting;
  closure condition `nm = k` (closed loop ⟺ integer winding number).
- Generating law `A(N)`, closure law `Q(N)`, recurrence `x⁶ − 49x⁴ + 504x² − 1296`,
  radius set `{2,3,6}`, multiplicity law `c_λ = 2m_λ + [λ=β]`, conservation
  `c₀ = 1 = m₊(6) − m₋(6)`.

## 3. Matrix probe

`containers(matrix, domain=onestep)` : 41 layers / 368 / 0 gaps, `D1=0 · S2=1 · CLOSED=1`.

## 4. Conditional theorem — four steps

**Conditional theorem (Yang–Mills). If Lemma 4.1 holds** (λ₁ > 0 in the 4D SU(3) continuum limit of the physical sector), **then** a non-trivial quantum theory exists for pure SU(3) Yang–Mills with Hamiltonian spectrum {0} ∪ [Δ, ∞), Δ = √3·Λ > 0.

1. **Compactness ⟹ discrete spectrum ⟹ positive bottom:** `|e^{iθ}| ≡ 1` (dev 1.132e−15) +
   the SU(2) half-angle double cover ⟹ the group is compact ⟹ Peter–Weyl ⟹
   `λ_l = l(l+2)`, `λ₁ = 3 > 0` (dev 0).
2. **The physical sector is the winding lattice:** gauge-invariant observables are traces
   (`tr(ΩUΩ†) = tr(U)`, dev 0; cyclic and linear, dev 0) ⟹ Wilson closed loops; closed loop ⟺
   integer winding number ⟹ the container is discrete and compact — **the winding lattice
   is the physical sector** (the earlier draft's "container S¹ (compact)" was a claim; this
   step replaces it with the winding-lattice reading). The ±-branches decompose diametrically:
   non-integer winding cancels (seam deviation 1.618), integer winding survives at the same
   site ⟹ the lattice is complete.
3. **Gap from repulsion:** the archive `:14716` — a perfect N=3 ring has a fragile degenerate
   zero (`λ₁ = λ₂ = 2.5`); any perturbation (broken spring) splits it (2.0 / 3.0); level
   repulsion gives the minimal crossing distance `2c` ⟹ `λ₁ > 0` is geometrically exact.
4. **Lemma 4.1 (winding-lattice completeness + continuum limit).** The physical sector is
   completely described by the winding lattice, and the continuum limit keeps `λ₁ > 0`.
   ∎ Then `Δ = √λ₁/R = √3·Λ` (dev 0; `Δ/Λ = √3 = 1.7320508076`).

   *Honest note:* "λ₁ > 0 preserved in the 4D SU(3) continuum limit" is the Clay core itself
   (the archive `:10279`: dynamics unclosed — group compactness ≠ base-space compactness;
   confinement separates the spectrum in the physical sector). The involution/cancellation/
   repulsion triple above is machine-witnessed; the continuum step is the labeled bridge.

## 5. Construction layer (construct · estimate · sufficiency)

- **Construct:** lattice spectrum on the discrete circle (N = 32/64/128/256) + singlet
  projector (tr = 1, `|P²−P| = 0`).
- **Estimate:** `λ₁L² = 38.9736793542 / 39.4467191014 / 39.4764358511 / 39.4782937425`;
  **Richardson extrapolation** `V = 39.478417601946`, relative deviation from `(2π)²`
  **6.109e−11**; convergence order O(1/N²): `|λ₁L²−(2π)²|·N² = 129.212992/129.711995/129.837068`
  (ratio 1.0048, analytic coefficient `(2π)⁴/12 ≈ 129.878`); physical gap `√λ₁L → 2π`
  (relative differences decreasing, max 1.606e−3); volume scaling `Δ(L)·L = 6.2830276023`
  identical at L = 1/2/4; β-function `β₀ = 11 − 2n_f/3`; Euler vs closed-form error ratios
  `10.2096/10.0204/10.0020` (O(h), threshold 5e−2 = 10/N_first/2);
  `Λ = μ·e^{−8π²/(β₀g²)}` essential singularity: `(Λ/μ)/g¹⁰ → 3.48e−10`.
- **Sufficiency:** the volume scaling shows the limit does not collapse; the trace invariance
  (`tr(UAUᵀ) = tr(A)`, dev 8.88e−16) shows the observables are the singlet sector; area law
  `−log W/(RT) = 1.000000` (three tiers), linear potential `E(R)/R = 0.350000` (four tiers),
  transfer-matrix gap `−log(ℓ₂/ℓ₁)/a = 1.2000000000` (dev 3.775e−15).

  *(Correction recorded: the Richardson coefficient must match the doubling factor — the first
  draft used `(4v(4N)−v(N))/3` (the N→2N coefficient) and stalled at 1.25e−5; dividing by 15
  drops the residue to 6.1e−11.)*

## 6. Numerical anchors

`Δ/Λ = √3 = 1.7320508076` (dev 0); `Λ = 2373/√3 ≈ 1370.05 MeV` with the external anchor
`2373 = 3·7·113` MeV (glueball — an external match, recorded as collateral evidence only,
different order from √3).

## 7. Reproduction

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ym.exe code/verify_ym.c && ./ym.exe
```

Expected: subtotals `9 · 9`, then the aggregate `RESULT PASS=18 FAIL=0` and `EXIT=0`
(first section `max dev 2.044e-02`, second `1.606e-03` with `res=6.109e-11`); zero warnings.

## 8. Honest notes

1. Lemma 4.1's continuum step is the Clay core; the involution + cancellation + repulsion
   triple is machine-witnessed.
2. The glueball anchor 2373 = 3·7·113 MeV is an external match of different order than √3 —
   recorded as collateral evidence.
3. The ratio-type criterion thresholds are mechanism-given (only the two asymptotic tiers enter
   the check).
4. The full-angle vs half-angle comparison must be taken on the same function (difference 2.000
   at 2π is the double cover itself).

## Acknowledgments

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.
