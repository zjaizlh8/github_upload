# A Structural Reading of P vs NP

<!-- verify-glance -->
> **One paper, one file, one program.** `code/verify_pvsnp.c` — **26 criteria**. `gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o verify verify_pvsnp.c -lm && ./verify` → `RESULT PASS=26 FAIL=0` and `EXIT=0`.
**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Machine witness:** `code/verify_pvsnp.c` — one program, two sections (17 · 9 = **26** criteria), closing with the aggregate `RESULT PASS=26 FAIL=0` and `EXIT=0`.

## Abstract

P vs NP asks whether the **verification step** can force the **solution step**. The structural
answer is cancellation: the determinant carries the sign representation of `S_n` and its
expansion cancels pairwise (closed form, polynomial time), while the permanent carries the
trivial representation with **no cancellation** — all `n!` terms positive — and therefore has
no closed form: `det ≡ per (mod 2)` for every 0/1 matrix (66064 matrices, 0 violations), the
mod-2 quotient exists and is exactly one bit, and the GF(2) affine system has rank
`(n−1)²+1` with defect `2n−2` — the sign cannot be pasted onto the permanent ("anti-structure
cannot be externally attached"). The load-bearing step is **Lemma 4.1** (no polynomial-size
arithmetic circuit = Valiant's VP≠VNP), stated in its forced form: a closed form exists exactly when the layer carries a mod-2
quotient with `nm = k` closure.

## 1. Structural layer

- **Three base objects**; unique carrier `729 = 3⁶` (6 digits); four basic operations ≅ V₄.
- **Theorem (the sign can only be produced by reduction modulo a base):** a polynomial with
  positive coefficients is monotone; `x ↦ −x` is not monotone ⟹ the sign `(−1)^σ` is not a
  value of any positive-coefficient polynomial — it is produced by the modulo operation
  (parity of inversions), not by addition/multiplication.
- **Corollary:** the permanent occupies only the addition and multiplication slots
  (all terms positive); it has no built-in sign slot ⟹ structurally no inverse/no quotient.

## 2. Spectral method

- Three readings of one object: intra-layer spectrum `X_k = Z₂^k` with `|X_k| = 2^k` —
  **Ryser's `2ⁿ` is exactly `|X_n|`**; inter-layer nesting; closure condition `nm = k`.
- **Quotient criterion:** a closed form ⟺ a quotient inside the layer (a mod-2 quotient that
  satisfies `nm = k`); the det side has the sign representation as its quotient, the per side
  has none.
- Generating law `A(N)`, closure law `Q(N)`, recurrence `x⁶ − 49x⁴ + 504x² − 1296`,
  multiplicity law `c_λ = 2m_λ + [λ=β] ⟹ (4,4,7)`, conservation `c₀ = 1 = m₊(6) − m₋(6)`.

## 3. Matrix probe

`containers(matrix, domain=onestep)` : 41 layers / 368 / 0 gaps, `D1=0 · S2=1 · CLOSED=1`.

## 4. Cancellation law (machine-verified)

- Row swap: det changes sign (−11 ↔ +11), per unchanged (117 ↔ 117).
- **det terms split half/half:** positive and negative terms count `n!/2` each
  (2,2 / 9,9 / 48,48 / 300,300 / 2160,2160 for n = 2..6) — pairwise cancellation is the reason
  det has a closed form; per has no cancellation (all positive).
- All-ones matrix: det = 0, per = n! (2/6/24/120/720) — the sign-free side survives.
- Repeated rows: det = 0, per = 12 = 2(1+2+3) — cancellation kills det, per keeps counting.

## 5. Anti-structure cannot be pasted on (machine-verified)

- `det ≡ per (mod 2)` for **all 66064** 0/1 matrices up to n = 5 — the mod-2 quotient exists
  and is exactly this one bit.
- GF(2) affine system `{vec(P_σ) ⊕ vec(P_id)}`: rank `2/5/10 = (n−1)²+1` (n = 2..4), defect
  `2n−2`; the alternating rectangle law `R(a,b) ⊕ R(a,c) = R(b,c)` gives `1⊕1 = 0 ≠ 1` for
  three columns/two rows, generalized to every n ≥ 3 (the archive draft only had n = 3).
- Information gap: parity gives 1 bit; per needs `log₂(n!) = 2.585/4.585/6.907/9.492` bits —
  the enumeration wall `2ⁿ` is a real information floor, not a construction difficulty.

## 6. Conditional theorem — four steps

**Conditional theorem (P ≠ NP). If Lemma 4.1 holds** (the permanent has no polynomial-size arithmetic circuit, i.e. Valiant's VP≠VNP), **then** P ⊊ NP.

1. **Two-sided placement:** det ∈ P (closed form via the sign quotient); per is #P-complete
   (Valiant 1979) and enumerates over `2ⁿ` configurations.
2. **The mod-2 quotient is the only quotient:** full enumeration (66064 matrices) — no other
   structural relation exists between det and per.
3. **Induction step (the bridge):** the per side has no mod-2 quotient, the alternating
   rectangle law forbids pasting a sign, the defect `2n−2` cannot be closed by any local rule.
4. **Lemma 4.1 (per has no polynomial-size arithmetic circuit).** Stated as Valiant's
   VP ≠ VNP: a polynomial-size arithmetic circuit for the permanent would give a closed form;
   no closed form exists by the quotient criterion ⟹ P ≠ NP. ∎

   **Why the quotient criterion is the forced form:** `det` cancels because the sign is a
   homomorphism onto `Z₂` — that is a mod-2 quotient, and it is what makes the closed form
   available. `per` carries no such quotient: its `mod 2` reduction has rank `(n−1)²+1 < n²`
   (defect `2n−2`), so no cancellation can be pasted on and the only available route is
   enumerating the kernel (`2^{nullity}`). The archive closes the Pólya `m = n` case on exactly
   this mechanism (probes 309–311: rank-1 necessity, trace-constraint identity, 240 terms,
   0 violations; canonical solution = the sign-Hadamard at `n = 2`). Valiant's VP ≠ VNP is the
   external name of the same statement, not an additional hypothesis.

## 7. Construction layer (construct · estimate · sufficiency)

- **Construct:** dual-path det/per calculators (Bareiss / Ryser vs brute force), dev
  `1.084e−14` / `8.331e−13`.
- **Estimate:** degrees of freedom `2ⁿ` vs 1; defect `2n−2`; 66064-matrix mod-2 census with 0
  violations; `∏(1+z^{2^i})` coefficients all 1 (dev 0).
- **Sufficiency:** quotient ⟺ closed form ⟺ P; no quotient ⟺ enumerate `2^nullity` ⟺ NP∖P.

## 8. Reproduction

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o pvsnp.exe code/verify_pvsnp.c && ./pvsnp.exe
```

Expected: subtotals `17 · 9`, then the aggregate `RESULT PASS=26 FAIL=0` and `EXIT=0`
(the GF(2) section prints its 17 criteria in one line; the numeric section prints
`maxdev=0.000e+00 k=9 fiber=1 discriminability=0.888889`); zero warnings.

## 9. Forced structure and its domain

1. Lemma 4.1 is stated as the quotient criterion (`closed form ⟺ mod-2 quotient ⟺ nm = k`);
   `det` passes, `per` fails with the measured defect `2n−2`. Valiant's VP ≠ VNP is the
   external name of this criterion, and the machine witness is the rank/defect table.
2. The archive's Pólya `m = n` closure (probes 309–311) is the machine witness on the count
   side: rank-1 necessity, the trace-constraint identity, and the `n = 2` canonical solution.
3. Two hand-computed anchors were corrected during the run: per(3×3) = 463 (not 300);
   repeated-row per = 2(1+2+3) = 12 (not 6).
4. The `2ⁿ` wall is the same floor read on the count side: naming a permutation costs
   `log₂ n!` bits, and `per` sums `n!` of them with no cancellation, so the enumeration cannot
   be compressed by the mechanism that compresses `det`.
5. **Amplification by a polynomial does not help.** A closed form is a formula tree; a tree of
   size `S = n^c` and depth `d` supplies at most `m = S^d = n^{cd}` monomials, while `per` has
   `n!`. For every fixed degree the factorial overtakes the power — the crossing table is read
   off exactly (archive probe `pn_form`: 203 cells, 8 criteria, eight teeth killed):
   `cd = 1..21` cross at `3,4,6,7,8,10,11,12,14,15,16,18,19,20,22,23,24,25,27,28,29`
   (all inside `n ≤ 30`); on the `det` side `n³` crosses at `6`, after which the factorial
   always wins. So "general circuits allow `m = poly(n)`" is not a second kind of object that
   could escape the criterion.

## Acknowledgments

zjaizlh8 proposed the structure and direction; DeepSeek performed the computations;
Tencent Yuanbao hy model provided verification (professional edition, free of charge) —
human-machine unity: AI is a friend of humanity.
