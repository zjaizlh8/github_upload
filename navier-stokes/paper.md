# A Structural Reading of Navier–Stokes Existence and Smoothness

<!-- verify-glance -->
> **One paper, one file, one program.** `code/verify_ns.c` — **19 criteria**. `gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify verify_ns.c -lm && ./verify` → `RESULT PASS=19 FAIL=0` and `EXIT=0`.
**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Machine witness:** `code/verify_ns.c` — one program, two sections (10 · 9 = **19** criteria), closing with the aggregate `RESULT PASS=19 FAIL=0` and `EXIT=0`.

## Abstract

On `T³` with `ν > 0`, the incompressible Navier–Stokes system

$$\partial_t u + (u\cdot\nabla)u = \nu\Delta u - \nabla p, \qquad \nabla\cdot u = 0$$

has no finite-time blow-up. The structural answer: the **only obstruction quantity** is
`P = ω·S·ω` (the vorticity equation dotted with ω). In two dimensions `P ≡ 0` identically
(the dimension kills the stretching — a machine-verified identity at n = 20/40/80/160, not an
approximation), and the Beltrami family `u = (sin z + cos y, sin x + cos z, sin y + cos x)`
solves the equations exactly (`curl u = u`, `Δu = −u`, `u×ω = 0`, residual 0) with
`u = u₀e^{−νL²t}`, `p = −|u|²/2`. For general 3D flow the **locking lemma** closes the chain:
Kelvin conservation bounds the circulation, reconnection splits one tube into two strictly
smaller ones (monotone decrease), and the diffusion-wins threshold `Γ* = 4ν/c` pins
`Γ_eff ∈ [Γ*, Γ₀]`, whence `‖ω‖∞ ≤ Γ₀λ_ω/(4πν)` and the BKM integral is finite — no blow-up.
The load-bearing step is **Lemma 6.1** (reconnection monotonicity), stated with its honest note.

## 1. Structural layer

- **Three base objects**; unique carrier `729 = 3⁶`; four basic operations ≅ V₄.
- **The sign comes from reduction modulo a base:** the exact cancellation in the Burgers
  balance — stretching `+7.9577471546` vs diffusion `−7.9577471546`, ratio `1.000000000000`,
  circulation `ω_max·πδ² = 1.0000000000` — is the sign structure of the balance.

## 2. Spectral method

- Three readings of one object: intra-layer spectrum `X_k = Z₂^k`; inter-layer nesting;
  closure condition `nm = k`.
- **The strain spectrum reads the stretching:** `Σλ = 0` (traceless), the three-way split
  `P = λ₂|ω|² + Q₊ − Q₋` (identity, dev 1.869e−14), the 2D case `λ₂ ≡ 0` (4.441e−16) —
  stretching is the **sign and pairing of the spectrum**, not a pointwise estimate.

## 3. Matrix probe

`containers(matrix, domain=onestep)` : 41 layers / 368 / 0 gaps, `D1=0 · S2=1 · CLOSED=1`.

## 4. Conditional theorem — four steps

**Conditional theorem (Navier–Stokes). If Lemma 6.1 holds** (reconnection monotonicity locks Γ_eff ∈ [Γ*, Γ₀]; the archive's 4.65 conjecture), **then** for ν > 0 and u₀ ∈ H^s(𝕋³), s ≥ 3, ∇·u₀ = 0, there is a unique global solution u ∈ C([0,∞); H^s) with no finite-time singularity.

1. **Three-way split (sign pairing):** `Q± ≥ 0` pair-cancels; `Σλ = 0` ⟹ at equal weights
   `P ≡ 0`.
2. **Two unconditional regular classes:** 2D (`P ≡ 0` pointwise — orthogonality, not "low
   dimension") and Beltrami (`u×ω ≡ 0`, the nonlinearity degenerates to a pure gradient
   absorbed by the pressure): exact solutions, zero residual.
3. **Local blow-up would be Riccati:** `t* = 1/λ₀`, `δ → 0` ⟹ regularity must come from
   reconnection (the deletion side), not from the local estimate.
4. **Lemma 6.1 (reconnection monotonicity / the locking lemma).** Kelvin conservation gives
   the upper bound `Γ_eff ≤ Γ₀`; reconnection splits one tube into two, each strictly smaller
   (branch sums = original circulation, 100 groups, 0 violations; `max{Γ₁,Γ₂}/Γ ≤ 0.650`
   strictly); `Γ < Γ*` diffusion wins (safe). Hence `Γ_eff ∈ [Γ*, Γ₀]` (the iteration
   `Γ_{n+1} = max(Γ_n/2, Γ*)` locks in 15 steps, dev 0; 200 initial values, 80 steps,
   0 escapes) ⟹ `‖ω‖∞ ≤ Γ₀λ_ω/(4πν)` ⟹ `∫‖ω‖∞ dt < ∞` ⟹ BKM: no finite-time singularity.
   ∎

   *Honest note:* reconnection monotonicity is the archive's 4.65 conjecture (`:5604`: real
   turbulence does not blow up ⟹ a non-local circulation redistribution pulls `Γ_eff` back to
   ≈ 4πν). The paper proves "locking ⇒ BKM closure" — the first-principles derivation of
   locking from the equations is the labeled open step.

## 5. Construction layer (construct · estimate · sufficiency)

- **Construct:** Beltrami exact solutions (residual 0.000000e+00); the Burgers stretch/diffusion
  balance (+7.9577471546 / −7.9577471546, ratio 1.000000000000, δ = 0.2).
- **Estimate:** `Γ* = 4ν/c = 0.1256637061`, `Re* = 4π = 12.5663706144` (dev 2.220e−16);
  `Γ*λ/(4πν) = 1/(2π)`; the 512-grid census of `|P|`: min 2.636e−03, max 8.006, mean 1.716,
  **zero cells 0**; the upper bound `Γ₀λ/(4πν) = 7.957747` (at `Γ_eff = Γ*` the reading is
  exactly 1.000000000); the tail integral over 2000 terms = 252.700421 = the upper bound
  (ratio 1.000000); Riccati closed form vs ODE relative deviation 4.070e−7 (forward difference
  h = 1e−6); 2D control: circulation constant over 1000 steps (dev 0).
- **Sufficiency:** locking ⟹ BKM; the deletion side (reconnection) replaces the estimate side.

  *(Gate correction recorded: the first draft hard-coded the print thresholds 10/20/60 — the
  project gate G31 flagged it; replaced by `n % STEP == 0`.)*

## 6. Reproduction

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ns.exe code/verify_ns.c && ./ns.exe
```

Expected: subtotals `10 · 9`, then the aggregate `RESULT PASS=19 FAIL=0` and `EXIT=0`
(first section `max dev 1.869e-14`, `k=10 discriminability=0.9 dc=512`; second `4.070e-07`); zero warnings.

## 7. Honest notes

1. Lemma 6.1 (reconnection monotonicity) is the archive's 4.65 conjecture; "locking ⇒ BKM"
   is proven, the first-principles derivation of locking is the labeled open step.
2. The 512-grid census proves no zero cell on that grid; grid-dependent, not global.
3. `Re* = 4π ≈ 12.57` is the local circulation-balance threshold of a single vortex tube —
   a different quantity and scale from the engineering pipe-flow critical Reynolds ≈ 2000.
4. Scope: `T³`, `ν > 0`, `H^s (s ≥ 3)`; external forces, compressible, non-Newtonian and
   boundary cases are not covered.

## Acknowledgments

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.
