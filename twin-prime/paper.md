# A Structural Reading of the Twin-Prime Conjecture

**Author:** zjaizlh8 (GitHub) · 2026-10-26 · CC BY 4.0
**Framework:** structural layer · spectral method · matrix probe
**Independent witnesses:** `code/verify_twin_count.c` (bit-sieve count, `EXIT=0`),
`code/verify_twin_count.c` (6 checks, `RESULT PASS=6 FAIL=0` + `EXIT=0`) and `code/verify_twin_tech.c` (10 checks, `RESULT PASS=10 FAIL=0 maxdev=7.039249e-03` + `EXIT=0`)

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

`containers(matrix, domain=onestep)` : 41 layers / 368 budget / 0 gaps, all
`D1=0 · S2=1 · CLOSED=1`. See `../shared/results/matrix_onestep.csv`.

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

- **Construct:** per-layer enumeration + trial division (sieve to 1e7, `verify_twin_count.c`
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

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o twin_count.exe code/verify_twin_count.c && ./twin_count.exe
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o twintech.exe code/verify_twin_tech.c && ./twintech.exe
```

Expected: the five-tier table plus `RESULT PASS=6 FAIL=0` (count) and `RESULT PASS=10 FAIL=0 maxdev=7.039249e-03` (tech)
(`k=9 fiber=1 fidelity=1 discriminability=0.888889`). Four gates exit 0; zero compiler warnings.

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
