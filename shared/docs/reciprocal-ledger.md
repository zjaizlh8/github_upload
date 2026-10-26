# Reciprocal Ledger — the Eight Papers' Load-Bearing Lines

Every load-bearing step of the eight papers can be written in one reciprocal (or
involutional) form. This file is the single ledger.

## The eight lines

| # | Paper | Reciprocal form | Machine witness |
|---|---|---|---|
| 1 | Hardy–Littlewood | `R₀ · e^{ΣΔ} = 4e^{−2γ} · e^{2γ}/4 = 1` (constants) and Möbius sign cancellation (error side, random-series behavior) | `verify_hl.c` (11 checks, dev 8.88e−16); `verify_hl_tech.c` (9 checks: \|E\|/√x = 0.0730/0.0771/0.0720, 2^ω ≤ 32) |
| 2 | Twin prime | `p · (N/p) = N` — factor pairs are antipodal about √N; sieving to √N is complete; 64 = 2⁶ = the circle-forced bound (ring length of layer 2) | `verify_twin_count.c` (π₂ = 1224/8169/58980/440312/3424506); `verify_twin_tech.c` (10 checks, chains/pairs = 0.14/0.64/1.59/4.66) |
| 3 | Riemann Hypothesis | `\|χ(s)\| · \|χ(1−s)\| = 1`; only at σ = 1/2 does the involution fix the zero set (unique-radius theorem: P8 d=(0,0,1) + P16 asymmetry concentrated on one radius + conservation c₀ = 1) | `verify_rh.c` (12 checks, dev ≤ 1.562e−09); `verify_rh_tech.c` (9 checks, argmin deviation 0, off/on-line ratio 8.921e+07) |
| 4 | P vs NP | `det ≡ per (mod 2)` for all 0/1 matrices; the mod-2 quotient exists only in this one bit — no closure outside it (rectangle alternation 1⊕1 = 0 ≠ 1) | `verify_pvsnp.c` (11 checks); `verify_pvsnp_tech.c` (9 checks, 66064 matrices, 0 violations; defect 2n−2) |
| 5 | Hodge | Circle averaging = unit-root reciprocal sum `Σ e^{i(p−q)·2πk/M} = M·[p≡q mod M]` (cancellation unless p = q); Abel–Jacobi return = reciprocal integrality `kx ≡ 0 mod 1` | `verify_hodge.c` (10 checks, dev 2.998e−15); `verify_hodge_tech.c` (9 checks: 1/12 returns at k = 12, √2/7 never returns, min distance 1.7745e−2) |
| 6 | Yang–Mills | Hodge star involution `* : Λ²R⁴` (spectrum 3/3); integer winding survives, non-integer cancels (antipodal decomposition); level repulsion λ₁ > 0 | `verify_ym.c` (9 checks); `verify_ym_tech.c` (9 checks: Richardson 6.109e−11, coefficient ratio 1.0048, Δ/Λ = √3 dev 0) |
| 7 | Navier–Stokes | Burgers stretching/diffusion exact cancellation (+7.9577471546 vs −7.9577471546); Γ* = 4ν/c, Re* = 4π; reconnection strictly splits circulation | `verify_ns.c` (10 checks, dev 1.869e−14); `verify_ns_tech.c` (9 checks: 100 groups zero violation, max{Γ₁,Γ₂}/Γ ≤ 0.650, tail integral 252.700421 = upper bound) |
| 8 | BSD | Central-point cancellation `Λ(1) = −Λ(1) ⟹ Λ(1) = 0` when w = −1; pairing = zero thickness (each pair eats exactly one (s−1)); multiplicity law `c_λ = 2m_λ + [λ=β]` ⟹ (4,4,7), Σc = 15 = 1 + 2Σm | `verify_bsd.c` (13 checks); `verify_bsd5.c` (10 checks: 389a1 L(1) = 0 at rank 2, beyond parity); `verify_bsd_tech.c` (9 checks: ∏c_p·#Ш = 2.000000000082) |

## Why "reciprocal" is the common name

- Antipode pairs cancel; the fixed point (center line) survives. `c₀ = 1 = m₊(6) − m₋(6)` is the
  single survivor — this is the same law in all eight papers.
- "6 and 7 are the root of the wave": 6 = the period (2·3 = the closed-product radius; ring
  length 64 = 2⁶), 7 = the pairing count (m = (2,2,3) ⟹ Σm = 7; A(N) has 7 lines,
  Q(N) = 1 + 2·7 = 15 lines). The six-order recurrence x⁶ − 7²x⁴ + (2³·3²·7)x² − 6⁴ has
  6 on the roots and 7 on the coefficients.

## Honest boundaries (one per paper)

1. HL: the unique open step is promoting the Möbius random series to a theorem.
2. Twin: Lemma 5.1 (every layer keeps at least one complete chain with both sides prime) is
   equivalent to the weak twin-prime statement; Brun's upper bound is the true-theorem anchor.
3. RH: "zero set ≅ spectral circle" is a structural isomorphism (system-caliber), not a
   standard theorem; the grid test is evidence, not a global proof.
4. P vs NP: the mod-2 obstruction is forced; per having no poly-size arithmetic circuit
   (Valiant, VP≠VNP) is the transfer lemma, open in full generality (Pólya m=n closed).
5. Hodge: Abel–Jacobi torsion for p ≥ 2 conflicts with Griffiths (1969); p = 1 (Lefschetz
   (1,1)) needs no lemma; CM abelian varieties are the external anchor.
6. YM: the continuum limit λ₁ > 0 in the physical sector is the Clay core (archive 4.126/10279).
7. NS: reconnection monotonicity is the archive 4.65 conjecture, not derived from the equation.
8. BSD: "pairing = zero thickness" is BSD itself; Gross–Zagier + Kolyvagin give rank 0/1.
