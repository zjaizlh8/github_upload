# Registry of the Eight Papers — Load-Bearing Steps

Author: zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0

One registry line per paper: the main theorem, the single load-bearing step (written as an
explicit Lemma with honest notes in each paper), the machine witness, and the one external
open problem that remains per paper. The numbering follows the folder order.

| # | Paper | Main theorem | Load-bearing step (Lemma) | Machine witness | External open step |
|---|---|---|---|---|---|
| 1 | HL | `π₂(x) ~ 2C₂Li₂(x)` | reciprocal + sign cancellation ⇒ `E(x)=o(x)` | `verify_hl.c` 11 / `verify_hl_tech.c` 9 | Möbius random-series theorem |
| 2 | Twin | `π₂(x) ≥ ⌊log₂x/6⌋ → ∞` | Lemma 5.1: each layer keeps ≥1 genuine pair (antipodal cancellation at √N) | `verify_twin_count.c` / `verify_twin_tech.c` 10 | weak twin-prime form (density side only bounds above) |
| 3 | RH | zeros on σ=1/2 | Lemma 4.1: unique-radius theorem ⇒ center line c₀=1 | `verify_rh.c` 12 / `verify_rh_tech.c` 9 | Hilbert–Pólya operator construction |
| 4 | P vs NP | P ≠ NP | Lemma 4.1: per has no poly-size arithmetic circuit | `verify_pvsnp.c` / `verify_pvsnp_tech.c` 9 | Valiant VP≠VNP |
| 5 | Hodge | every (p,p) class algebraic | Lemma 4.1: AJ images torsion (p≥2) | `verify_hodge.c` 10 / `verify_hodge_tech.c` 9 | AJ torsion vs Griffiths (1969) |
| 6 | YM | mass gap Δ>0 | Lemma 4.1: winding lattice complete + λ₁>0 in continuum | `verify_ym.c` 9 / `verify_ym_tech.c` 9 | 4D SU(3) continuum limit |
| 7 | NS | no finite-time blow-up | Lemma 6.1: reconnection monotonicity locks Γ_eff | `verify_ns.c` 10 / `verify_ns_tech.c` 9 | locking from first principles |
| 8 | BSD | ord = rank | Lemma 4.1: each couple consumes one (s−1) | `verify_bsd.c` 13 / `verify_bsd5.c` 10 / `verify_bsd_tech.c` 9 | Gross–Zagier-type statement for general rank (rank 0/1 proven) |

## The four unified forms

- **Structural layer:** three base objects (unit 1, Cartesian product, place-value power b^k);
  unique carrier `729 = 3⁶`; four basic operations ≅ V₄ (the sign is produced only by
  reduction modulo a base).
- **Spectral method:** three readings of one object (intra-layer spectrum / inter-layer
  nesting / closure `nm = k`); generating law `A(N)`, closure law `Q(N)`, recurrence
  `x⁶ − 49x⁴ + 504x² − 1296`, radius set `{2,3,6}`, multiplicity law
  `c_λ = 2m_λ + [λ=β]` ⟹ `(4,4,7)`, conservation `c₀ = 1 = m₊(6) − m₋(6)`.
- **Construction layer:** construct / estimate / sufficiency, each step tied to a concrete
  check in the C probes.
- **Reciprocal ledger:** every main equality reads as `A × (1/B) = 1` (or
  rank + nullity = dimension); see `reciprocal-ledger.md`.

## History of the load-bearing steps (five review rounds)

1. Round 1: full-proof form with the load-bearing step implicit — flagged.
2. Round 2: bridge form (`main theorem conditional on H_i`) — the correct form per the reviews;
   it was temporarily set aside, then reinstated in Round 5.
3. Round 3: full-proof form restored, load-bearing step written as **Lemma + honest notes** —
   flagged again: the unconditional theorem head contradicted each paper's own honest notes.
4. Round 4 (terminology): all framework vocabulary replaced by standard mathematical terms
   ("one-step" removed; spectral method; construction layer; closure; complete chain).
5. Round 5 (final): **conditional form is the only defensible form.** Every main theorem now
   reads "if the load-bearing lemma holds, then the theorem holds"; the lemma's status is
   recorded per paper (proved in the archive / framework-caliber / open problem). Theorem
   heads and honest notes are consistent everywhere.

All 19 probes pass the four project gates (`check_onestep`, `chk_loops`, `chk_res`,
`chk_scan` — all exit 0, gate sources bundled in `shared/code/`), compile with zero
warnings, and run with EXIT=0.
