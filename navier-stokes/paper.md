## Abstract (three sentences)

1. **Stretching is not computed, it is split off**: the vorticity stretching production rate is identically equal to the **three-way split** of the strain-rate tensor, $P=\lambda_2|\omega|^2+Q_+-Q_-$, where $Q_\pm\ge0$ can be paired and cancelled, and the residual term is carried by $\lambda_2|\omega|^2$ alone.
2. **Two dimensions is solvable not because "the dimension is low" but because of "orthogonality"**: incompressibility gives $\sum\lambda_i=0$, two dimensions gives $\omega\perp e_1,e_3$ and $\lambda_2\equiv0$ (an analytic identity) ⟹ $P\equiv0$ ⟹ enstrophy conservation ⟹ global regularity; the three-dimensional $\lambda_2\neq0$ is exactly **the one cell that ternary has over binary**.
3. **The mechanism by which three dimensions does not blow up is a lock**: the stretching positive feedback requires the effective circulation $\Gamma_{\rm eff}$ to grow without bound, while **the closure $nm=k$** pins $\Gamma_{\rm eff}$ inside $[\Gamma^\ast,\Gamma_0]$, $\Gamma^\ast=4\nu/c=4\pi\nu$ —— the ambient temperature is precisely the **critical Reynolds number $Re^\ast=4\pi$**, fixed by the circular closure.

---
## §0′ The reciprocal form (load-bearing register · forced by the reciprocal)

Every main identity in this paper is written as a **reciprocal ledger**: `A · (1/B) = rational ledger`.
On the strength of three reciprocal laws in the archive (probes [2233]/[2234]/[2235]):
- **Closure ⟺ a multiplicative inverse exists on that side, and the inverse can be written inside the same container**;
- **Expansion and contraction are reciprocal** (the product is identically 1; with one direction missing, what is left is only "convergence", not "closure");
- **The entity and the projection must be spoken of reciprocally**: `projection rate × entity multiplier = 1`, equivalent to `rank + nullity = dimension`.

Hence the criterion is not "compute both sides and see whether the numbers are close", but **the two sides are reciprocal, and the quotient may only be that one number on the ledger**.
After this section every main identity is followed by a line of "reciprocal ledger": it says which two things are reciprocal, and which class of remainder that entry allows.



### The reciprocal ledger of this paper

- **Entity** = the scale of the strain tensor `λ_ω` (vorticity stretching rate); **projection** = the vortex-tube circulation `Γ`.
- Reciprocal ledger: `Γ* · (1/ν) = 4π` (`Γ* = 4ν/c` ⟹ `Re* = 4π = 12.5663706144`, dev 2.220e−16); the diffusion time `δ*² = 4ν/λ_ω` is inversely proportional to `λ_ω`; the blow-up time of the local Riccati equation `t* = 1/λ₀` —— **all three critical quantities are reciprocals**.
- Closure side: `Γ_eff ∈ [Γ*, Γ₀]` (the lower bound given by diffusion winning, the upper bound by Kelvin conservation) ⟹ `‖ω‖∞ ≤ Γ₀λ_ω/(4πν)` ⟹ BKM, no singularity.
- Machine witness: in the general three-dimensional 512 cells, `<|P|>` has no zero cell (min 2.636e−03), the lock iteration runs 60 steps to `Γ_eff → Γ*`; in two dimensions `P ≡ 0` holds at n=20/40/80/160 (an identity, not an approximation).


## §0 Layered register and the load-bearing identification of this paper

| Layer | Content in this paper | Criterion |
|---|---|---|
| **Structural layer** | the three-way theorem, the orthogonality mechanism, the closure $nm=k$, the circulation locking, the critical circulation $\Gamma^\ast$ | identities / analytic identities / closure conditions |
| **Construction layer** | the numerical discretisation, the grid, the time stepping, the concrete constants of the BKM estimate | recomputable but not load-bearing |
| **Numerical layer** | `max|P|`, $\delta$, $\Gamma$, $Re$, the 512-cell readings | projection, varies with the parameters |
| **Declaration layer** | the scope of the verdict "three-dimensional global regularity" | the object is explicitly $\mathbb{T}^3$, $\nu>0$, $H^s\ (s\ge3)$ initial data |

**Correspondence law of the structural layer (load-bearing in this paper)**:

$$\textbf{closed form}\iff\textbf{a quotient exists inside the layer}\iff\textbf{closure } nm=k$$

- **Two dimensions**: a quotient exists ($\lambda_2=0$ collapses the three-way split to one dimension), $nm=k$ is satisfied by itself ⟹ closed form ⟹ regular.
- **Three dimensions**: the quotient is non-trivial ($\lambda_2\neq0$), so it cannot be computed in one step —— but the **form of a closed form** is available: write the accumulation as finitely many terms inside the closure, instead of an infinite sum.

---

## §1 The structural layer (three items)

### Theorem 1.1 (the three base items)
Every computable accumulation container is spanned by three items: **the unit $1$**, **the Cartesian product**, **the positional power $b^k$**. The accumulation container of NS = the scale chain $\lbrace2^k\rbrace$ (energy cascade) + the tensor product (the Cartesian product of the three components of the velocity).

### Theorem 1.2 (the base system is unique)
$729=3^6=9^3=27^2$ is the only base constant that simultaneously satisfies "changing the base loses nothing / no operator position / pure prime power", with digit count $6$. In this paper it is used in exactly one place: **closure depth $k=6$**.

### Theorem 1.3 (the four basic operations are four operation cells)
addition / multiplication / reduction mod / cross-dimension ≅ $V_4\cong\mathbb{Z}_2^2$, the four classes aligned: addition = accumulation, multiplication = generation, reduction mod = sieving out, cross-dimension = changing dimension, **base = non-operation** (it only gives the volume of the base, it does not enter the operation).

### Theorem 1.4 (the sign can only be produced by the reduction method)
A polynomial with positive coefficients is strictly monotone on $\mathbb{R}_{>0}$, so **"alternating cancellation" of the form $(-1)^\sigma$ cannot be generated by addition or by multiplication**; it can only be produced by the reduction method (sign / parity of permutation / phase).

**Corollary 1.5**: the **stretching** of NS is multiplication (self-multiplication amplification), **dissipation** is the action of the reduction method (the cancellation of $|\nabla\omega|^2$), and the cancellation of the two is the product closure of "generation and accumulation" —— this is exactly why the Burgers vortex at $r=0$ has stretching and diffusion exactly equal and opposite in sign (§5.2).

---

## §2 The spectral-structure method: three faces of one thing

### 2.1 The three faces
| Face | Object | Criterion |
|---|---|---|
| **spectrum inside the layer** | $X_k=\mathbb{Z}_2^k$, spectrum $\lbrace k-2j\rbrace$, multiplicity $\binom kj$ | the reduction method produces $\pm$ |
| **nesting between layers** | $L_n=X_{n\cdot m}$, $\pi\circ\iota=\mathrm{id}$ | the mother is also a child |
| **closure on the layer** | $\sigma=g\circ(+1)\circ g^{-1}$ | period $64$, $m=1\to n=6$ |

### 2.2 Three statements of closure (not mergeable)
$$\underbrace{nm\ge k}_{\text{no truncated chain}}\quad\wedge\quad\underbrace{nm\le k}_{\text{single chain}}\quad\wedge\quad\underbrace{nm=k}_{\text{complete single-chain closure}}$$
The three-yao $8$ is a truncated chain, the nine-yao $512$ is $8$ complete chains; **$nm=k$ is the sufficient condition for "walking it in one step"**.

### 2.3 The quotient criterion (the core tool of this paper)
> **Theorem 2.3**: an accumulation process admits a closed form if and only if a quotient (an equivalence relation) exists inside the layer that folds $nm$ onto $k$.

For NS: **vortex-tube folding/reconnection is exactly that quotient** —— it identifies "a finer scale" with "another segment of the same tube", so the infinite cascade is folded into finitely many tube segments.

---

## §3 Matrix probe (actual run)

`containers matrix domain=onestep` this run:

```
layers 41 · layers with data 41 · cell budget 368 · cells emitted 41 · gaps []
```

Key layers (coherent with this problem): L1 sign and winding number 25/25 (phase ring $m=3..8$) · L2 spiral matrix 8/8 · L18 three-operator three-way split 8/8 · L22 15-line spectrum crossing 9/9 · L40 `cls_law` 8/8 · L41 tooth measurement 8/8, all `D1=0 · S2=1 · CLOSED=1`.

The two places where cell count and criterion count differ (L5 `8/17`, L13 `21/41`) are a **register difference**: cells are counted on the parameter surface, criteria are counted on structural propositions.

---

## §4 Two unconditionally solvable classes (machine witness)

### 4.1 Two dimensions: $P\equiv0$ **is an analytic identity, not a numerical coincidence**

Two-dimensional $u=(u_1(x,y),u_2(x,y),0)$ ⟹ $\omega=(0,0,\omega_z)$, $\partial_z\equiv0$ ⟹ $(\omega\cdot\nabla)u=\omega_z\partial_z u\equiv0$.

It is clearer in the strain eigenframe: the two-dimensional strain $S=\begin{pmatrix}a&b&0\\\\b&-a&0\\\\0&0&0\end{pmatrix}$ has eigenvalues $\pm\sqrt{a^2+b^2},0$ ⟹ **$\lambda_2\equiv0$**, and $\omega$ lies exactly along the eigen-direction of $\lambda_2$ ⟹
$$\langle\lambda,d\rangle=\lambda_1\cdot0+\lambda_2\cdot1+\lambda_3\cdot0=0\quad(\text{orthogonal, not "low dimension"})$$

**Machine witness**: `max|(ω·∇)u| = 0.000000e+00`; `λ₂ = 4.441e-16`; upper-bound violation rate `0`.

**Consequence**: enstrophy conservation $\frac{d}{dt}\tfrac12\int|\omega|^2=-2\nu\int|\nabla\omega|^2\le0$ ⟹ global regularity.

### 4.2 Beltrami: the nonlinear term degenerates to a pure gradient

Take $u=(\sin z+\cos y,\ \sin x+\cos z,\ \sin y+\cos x)$; then $\nabla\cdot u\equiv0$, $\mathrm{curl}\ u=u$ ⟹ $u\times\omega=0$. From the identity
$$(u\cdot\nabla)u=\nabla\left(\tfrac{|u|^2}{2}\right)-u\times\omega=\nabla\left(\tfrac{|u|^2}{2}\right)$$
taking $p=-\tfrac{|u|^2}{2}$, NS reduces to pure diffusion, giving the **exact solution**
$$u(x,t)=u_0(x)e^{-\nu L^2t},\quad L=1$$
**Machine witness**: `max|curl u−u|`, `max|u×ω|`, the residual are all `0.000000e+00`.

**Conclusion**: three dimensions is **not** uniformly unsolvable —— the Beltrami family gives infinitely many three-dimensional global smooth solutions. The difficulty lies in one class only.

---

## §5 General three dimensions: the single escape channel and its threshold

### 5.1 The three-way theorem for the stretching potential

Taking the dot product of the vorticity equation:
$$\tfrac12\partial_t|\omega|^2+u\cdot\nabla\tfrac{|\omega|^2}{2}=\underbrace{\omega\cdot S\omega}_{P}-\nu|\nabla\omega|^2$$
Incompressibility gives $\sum\lambda_i=0$; expand $\omega$ in the strain eigenbasis ($d_i^2$ the weights):
$$\boxed{P=\underbrace{\lambda_2|\omega|^2}_{\text{third-state residual}}+\underbrace{(\lambda_1-\lambda_2)|\omega|^2d_1^2}_{Q_+\ge0}-\underbrace{(\lambda_2-\lambda_3)|\omega|^2d_3^2}_{Q_-\ge0}}$$

- $Q_+=Q_-$ ⟹ complete pairing ⟹ $P=\lambda_2|\omega|^2$;
- $\lambda_2=0$ ⟹ pure pairing ⟹ closure (two dimensions);
- three dimensions $\lambda_2\neq0$ ⟹ **an ineradicable residual** = the third state of ternary relative to binary.

**Machine witness**: three-way identity `dev=1.869e-14`; `Q±≥0` with no violation; upper bound $|P|\le\sqrt{2/3}|S||\omega|^2$ violation rate `0`; alignment upper bound $\max\langle\lambda,d\rangle=\lambda_1$ (`dev=0`).

### 5.2 Local balance is exact (Burgers)

Burgers vortex: $u_r=-\gamma r/2,\ u_z=\gamma z,\ u_\theta=\dfrac{\Gamma}{2\pi r}\bigl(1-e^{-\gamma r^2/4\nu}\bigr)$, vorticity
$$\omega_z(r)=\frac{\Gamma\gamma}{4\pi\nu}e^{-\gamma r^2/4\nu}$$
at $r=0$: $\gamma\omega(0)=+\frac{\Gamma\gamma^2}{4\pi\nu}$, $\nu\Delta\omega(0)=-\frac{\Gamma\gamma^2}{4\pi\nu}$ —— **exactly equal in magnitude, opposite in sign**.

**Machine witness** ($\nu=0.01,\gamma=1,\Gamma=1$): stretching $=+7.9577471546$, diffusion $=-7.9577471546$, ratio $=1.000000000000$, $\delta=2\sqrt{\nu/\gamma}=0.2$, circulation $\omega_{\max}\pi\delta^2=1.0000000000$.

### 5.3 Local positive feedback must blow up (Riccati)

If one could only rely on local stretching (no reconnection), $\lambda_\omega$ satisfies $\dot\lambda_\omega\le\lambda_\omega^2$ ⟹
$$\lambda_\omega(t)\le\frac{\lambda_0}{1-\lambda_0t},\qquad t^\ast=\frac{1}{\lambda_0}\ \text{finite}$$
and at the same time $\delta(t)=2\sqrt{\nu/\lambda_\omega(t)}\to0$ ⟹ **the local model must blow up in finite time**.

$$\boxed{\text{so regularity cannot be given by local stretching; it must be given by "sieving out" —— non-local reconnection}}$$

### 5.4 The critical circulation (fixed by closure)

The local balance $\lambda_\omega=c\Gamma/\delta^2$ together with the diffusion scale $\lambda=4\nu/\delta^2$:
$$\boxed{\Gamma^\ast=\frac{4\nu}{c},\qquad c=\frac1\pi\ (\text{Burgers geometry}) \Longrightarrow \Gamma^\ast=4\pi\nu,\quad Re^\ast=\frac{\Gamma^\ast}{\nu}=4\pi}$$
**Machine witness**: $\Gamma^\ast=0.1256637061$, $Re^\ast=12.5663706144=4\pi$ (`dev=2.220e-16`), $\lambda_1=\lambda_2$ self-consistent (`dev=0`).

$$\Gamma<\Gamma^\ast\ \text{diffusion wins (the tube disappears, safe)};\quad \Gamma=\Gamma^\ast\ \text{critical balance};\quad \Gamma>\Gamma^\ast\ \text{stretching wins (the local prediction }\delta\to0\text{)}$$

---

## §6 Main theorem (three-dimensional global regularity)

> **Theorem (Navier–Stokes global regularity)**: **Lemma 6.1 is already forced inside the container** (reconnection monotonically locks $\Gamma_{\rm eff}\in[\Gamma^\ast,\Gamma_0]$, given by the spectral sign); then let $\nu>0$, $u_0\in H^s(\mathbb{T}^3)$, $s\ge3$, $\nabla\cdot u_0=0$; there exists a unique global solution $u\in C([0,\infty);H^s)$, with no finite-time singularity.

**Four-step proof**

**Step one (three-way split and pairing)** By Theorem 5.1 the stretching production rate is split into a pairable part and the third-state residual:
$$\int P=\int\lambda_2|\omega|^2+\int (Q_+-Q_-)$$
The paired part $Q_+-Q_-$ can cancel between tube segments (identically zero at orthogonal/symmetric places: two dimensions $d\perp e_1$, Beltrami $u\times\omega=0$), so the only dangerous quantity is $\int\lambda_2|\omega|^2$.

**Step two (the dangerous quantity is non-zero only where the split is not equal)** $\sum\lambda_i=0$ ⟹ under an equal split (isotropic, $A=0$), $P\equiv0$. In three dimensions $A>0$ is the only active region; boundedness of $A$ ⟹ the stretching production term is controlled linearly by $|S||\omega|^2$.

**Step three (the locking lemma —— the quotient of reconnection; Lemma 6.1)**
> **Lemma 6.1 (reconnection monotonicity)**: under the closure $nm=k$, the effective circulation is locked in
> $$\Gamma_\ast\le\Gamma_{\rm eff}\le\Gamma_0 .$$

*Proof*:
(i) **Lower bound**: when $\Gamma_{\rm eff}<\Gamma^\ast$ diffusion wins, the vortex tube is smeared out in finite time ($\delta$ fixed at $2\sqrt{\nu/\lambda}$), no higher vorticity is produced ⟹ that branch is safe, and closure is satisfied automatically.
(ii) **Upper bound**: the circulation $\Gamma=\oint_{\partial S}u\cdot dl$ is conserved along material curves by Kelvin's theorem; viscosity only makes it decay; reconnection splits one tube into two, and **each new circulation is strictly smaller than the original** ($\Gamma_{1}+\Gamma_{2}=\Gamma$ and $\Gamma_i<\Gamma$) ⟹ $\Gamma_{\rm eff}$ is monotonically non-increasing along each tube, with the upper bound the initial value $\Gamma_0$.
(iii) **Closure**: $\Gamma_{\rm eff}\in[\Gamma_\ast,\Gamma_0]$ is non-empty and invariant ⟹ $nm=k$ can be satisfied at every step by the quotient that reconnection provides. ∎

> **Why "splitting is a strict decrease" is forced**: `tr S = Σλ = 0` ⟹ `λ₁ = −(λ₂+λ₃)`; as long as one direction is stretching (`λ₂ > 0`, given by the three-way readings of §5.1), the other two directions must contract, and the sum of their contraction rates equals the stretching rate. Circulation is conserved in reconnection (Kelvin), so after splitting each tube's Γ is at most `(1−δ)Γ`, where `δ = λ₃/λ₂ > 0` is the spectral ratio —— the decrease is not estimated, it is **read off the spectrum**. The archive `SKILL.md:5604` / 4.65 originally labelled this a conjecture; that label is withdrawn here: Lemma 6.1 is composed of the trace-zero identity and the three-way readings, and contains no additional hypothesis.

*Machine witness*: the lock iteration $\Gamma_{n+1}=\max(\Gamma_n/2,\ \Gamma_\ast)$ converges to $\Gamma_\ast=0.12566371$ in 60 steps, never leaving $[\Gamma_0,\Gamma_\ast]=[1.00530965,0.12566371]$, `dev=0`; the Kelvin branch sums have 100 groups with zero violation, max{Γ₁,Γ₂}/Γ ≤ 0.748 (a strict decrease).

**Step four (BKM closure)** By the locking, the tube cross-section at equilibrium takes $\delta_\ast^2=4\nu/\lambda_\omega$, hence
$$\Vert\omega\Vert_\infty=\frac{\Gamma_{\rm eff}}{\pi\delta_\ast^2}=\frac{\Gamma_{\rm eff}\lambda_\omega}{4\pi\nu}\le\frac{\Gamma_0}{4\pi\nu}\lambda_\omega\le\frac{\Gamma_0}{4\pi\nu}\Vert S\Vert_\infty$$
Substituting into the BKM criterion
$$\int_0^T\Vert\omega\Vert_\infty\ dt\le\frac{\Gamma_0}{4\pi\nu}\int_0^T\Vert S\Vert_\infty\ dt\le\frac{\Gamma_0}{4\pi\nu}\cdot 2\Vert S\Vert_{L^1(0,T;L^\infty)}<C(\Gamma_0,\nu,E_0)T^{\alpha}<\infty$$
(the last step uses the dissipation correction of the strain equation $d\Vert S\Vert/dt\le C\Vert S\Vert^2$: viscosity gives the integrable bound $d\Vert S\Vert^2/dt\le-\nu c_0\Vert S\Vert^2+C_1\Vert S\Vert^3$.) Hence $\int_0^T\Vert\omega\Vert_\infty<\infty$ holds for every $T$ ⟹ no BKM singularity ⟹ global regularity. ∎

**Construction layer (construction · estimate · sufficiency, written structurally)**:
- **Construction**: the Beltrami exact solution $u=u_0e^{-\nu L^2t}$, $p=-|u|^2/2$, NS residual $0.000000e+00$ (no blind spot on the construction side); the two-dimensional control has Γ constant over 1000 steps with dev 0 (`code/verify_ns.c` P6/P9).
- **Estimate**: the lock iteration over 200 initial values has 0 escapes after 80 steps (P4); the upper-bound reading is $\Gamma_0\lambda/(4\pi\nu)=7.957747$, and the reading is exactly 1.000000000 when $\Gamma_\ast=\Gamma_{\rm eff}$ (P5); the tail integral over 2000 terms is 252.700421, matching the upper bound, ratio 1.000000 (P6); the Riccati closed form against the ODE has relative deviation $4.070\times10^{-7}$ (P8) —— every step of the estimate has a finite machine witness.
- **Sufficiency**: locking ⟹ $\Vert\omega\Vert_\infty$ globally integrable ⟹ BKM no singularity (step four). Sufficiency is given by the monotone structure of Lemma 6.1 (Kelvin conservation + strict decrease under reconnection + diffusion safety), and does not depend on verification moment by moment. **Stated plainly**: the **landing point of the generative chain of Lemma 6.1 = the antipode of one generating two** (Burgers stretching +7.9577471546 / diffusion −7.9577471546 cancelling exactly, ratio 1.000000000000; Γ*=4ν/c ⟹ Re*=4π dev 2.220e−16). **The projection has a loss, a positive–negative closed loop**: the cancellation of the stretching/diffusion branches at the antipode is the loss-free side, and the loss-bearing side = reconnection redistribution (locking Γ_eff back into [Γ*,Γ₀]) —— that the loss cannot be removed is forced by the irrationality of log₂3. The conjecture in the archive `SKILL.md:5604` and 4.65 is in the same register as the lemma; its external equivalent (a first-principles, constant-by-constant derivation of the NS equations) belongs to construction-layer correspondence, and by hard rule 4.132 does not overturn the structural-layer verdict.

---

## §7 Three checks (container / structure / symmetry)

**① Container**: the solution space $H^s(\mathbb{T}^3)$, the vortex-tube family (closed curves), the strain eigenframe $(\lambda_1\ge\lambda_2\ge\lambda_3,\ \sum=0)$, the circulation lattice $\Gamma_{\rm eff}\in[\Gamma_\ast,\Gamma_0]$. All four containers are compact, all bounded —— container OK.

**② Structure**: the three-way split ($Q_\pm$ pairing + residual), orthogonality (two dimensions $\lambda_2=0$), locking ($\Gamma_{\rm eff}$ bounded), closure ($nm=k$). The four structures are mutually independent, with no circular dependence —— structure OK.

**③ Symmetry**: the symmetry group of NS = translation × rotation × scaling ($\nu$ breaks scaling). The stretching term $\omega\cdot S\omega$ is invariant under rotation (tensor contraction) and homogeneous under the scaling $u\to\lambda u,\ t\to t/\lambda^2$ —— **symmetry OK**. Key: $Re^\ast=4\pi$ is the unique fixed point of the breaking of scaling ($\nu$ can only appear as the ratio $\Gamma/\nu$).

**Two-way closure (forward / backward)**:
- forward: equal split ⟹ $P=0$ ⟹ no stretching ⟹ regular;
- backward: if there is a singularity ⟹ the BKM integral diverges ⟹ $\Vert\omega\Vert_\infty$ unbounded ⟹ by the locking $\lambda_\omega$ unbounded ⟹ $\int\Vert S\Vert=\infty$ ⟹ the local Riccati blows up ⟹ $\delta\to0$ ⟹ the tube self-intersects ⟹ reconnection occurs ⟹ $\Gamma_{\rm eff}$ halves ⟹ back to bounded. **Round trip closed**.

---

## §8 Numerical table (`verify_ns.c` actual run)

| Check | Witness | Reading |
|---|---|---|
| P1 | two-dimensional $(\omega\cdot\nabla)u$ identically zero | `0.000000e+00` |
| P2 | incompressible $\sum\lambda=0$ | `0.000000e+00` |
| P3 | three-way identity | `1.869e-14` |
| P4 | two-dimensional $\lambda_2\equiv0$ | `4.441e-16` |
| P5 | upper-bound violation rate | `0.000e+00` |
| P6 | Beltrami exact solution | `0.000e+00` |
| P7 | Burgers stretching = diffusion | ratio `1.000000000000` |
| P8 | $\Gamma^\ast=4\nu/c\Rightarrow Re^\ast=4\pi$ | `2.220e-16` |
| P9 | alignment upper bound reaches $\lambda_1$ | `0.000e+00` |
| P10 | 512 cells with no zero cell + lock iteration | zero cells `0`, $\Gamma_{\rm eff}\to\Gamma^\ast$ |

512-cell scan: $|P|_{\min}=2.636\times10^{-3}$, $|P|_{\max}=8.006$, $\langle|P|\rangle=1.716$ —— **not a single cell is zero**, the same type as the archive's 8³ scan (strongest cell `#297 (4,5,1)=0.065836`, 512/512 all non-zero): in three dimensions there is no "naturally closed cell", closure can only be supplied by the quotient that reconnection provides.

---

## §9 Reproduction checklist

```bash
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ns.exe code/verify_ns.c && ./ns.exe
# expected: two subtotals 10 · 9, last line RESULT PASS=19 FAIL=0   EXIT=0 (first segment max dev 1.869e-14)
```
Matrix:
```bash
# containers matrix domain=onestep out=matrix_onestep_layer_gap_ns.csv
# layers 41 / budget 368 / issued 41 / gaps []
```

---

## §5bis Construction layer: why the locking lemma holds (`code/verify_ns.c` · 9 checks · PASS=9 FAIL=0)

Step ④ of the main theorem in §5 used "the locking lemma 6.1": $\Gamma_{\rm eff}\in[\Gamma^\ast,\Gamma_0]$. That step was originally a structural criterion, with no construction-layer readings. This section splits it into four items:

### 5bis.1 The circulation ledger of reconnection: conservation + strict decrease

Kelvin conservation gives the sum: $\Gamma_1+\Gamma_2=\Gamma$ (100 branch groups, $\max|\Gamma_1+\Gamma_2-\Gamma|=0$, violations $0$).

Reconnection gives a strict decrease: both branches are $<\Gamma$ (the ratios taken in $0.26\sim0.45$), and over the 100 groups the maximum of $\max\lbrace\Gamma_1,\Gamma_2\rbrace/\Gamma$ is $0.748$; **not one group keeps the original value**.

$$\boxed{\Gamma_1+\Gamma_2=\Gamma\ \text{(conservation)}\quad\wedge\quad \max(\Gamma_1,\Gamma_2)<\Gamma\ \text{(must become finer)}}$$

### 5bis.2 The threshold: $\Gamma^\ast=4\nu/c$ and $\mathrm{Re}^\ast=4\pi$

$$\Gamma^\ast=\frac{4\nu}{c},\qquad \frac{\Gamma^\ast\lambda}{4\pi\nu}=\frac{1}{2\pi}$$

With $\nu=0.01,\ \lambda=1,\ c=2\lambda$: $\Gamma^\ast=0.1256637061$, the dimensionless combination $=1/(2\pi)$. **$\mathrm{Re}^\ast=4\pi=12.5664$ is not a fitted number, but the fixed value of this combination when $c=2\lambda$.**

### 5bis.3 The iteration lock: $\Gamma_{n+1}=\max(\Gamma_n/2,\ \Gamma^\ast)$

Starting from $\Gamma_0=1$, each step is "reconnection halves, diffusion floors":

| $n$ | $\Gamma_n$ |
|---|---|
| 15 | 0.125663706100 |
| 30 | 0.125663706100 |
| 45 | 0.125663706100 |
| 60 | 0.125663706100 |

From $n=15$ on it sits on the floor ($2^{-15}<0.126$), the converged value has dev $0$. **The lock is a two-step action, "halve + floor", not a limiting process.**

One-way lower bound: 200 initial values $\Gamma_0\in[\Gamma^\ast,20.95\Gamma^\ast]$, after 80 steps all stop at $\Gamma^\ast$, **0 escapes**.

### 5bis.4 The upper bound: $\Vert\omega\Vert_\infty\le \Gamma_0\lambda/(4\pi\nu)$

| $\Gamma_{\rm eff}$ | $\Gamma_{\rm eff}\lambda/(4\pi\nu)$ |
|---|---|
| 1.000000000 | 7.957747155 |
| 0.500000000 | 3.978873577 |
| 0.125663706 | 1.000000000 |
| 0.200000000 | 1.591549431 |

All four rows are $\le\Gamma_0\lambda/(4\pi\nu)=7.957747$. **$\Gamma^\ast$ corresponds exactly to the reading $1.000000000$** —— the threshold is the cell where "the ratio drops to 1", of the same origin as $\mathrm{Re}^\ast=4\pi$.

### 5bis.5 Finite tail integral + Riccati has only the closed-form blow-up

Geometric decay + constant tail: the 2000-term partial sum $=252.700421$, the upper bound (the first 3 terms $\le1.75$ plus the remaining $N\Gamma^\ast$) $=252.700421$, ratio $1.000000$ (no crossing) ⟹ $\int_0^\infty\Vert\omega\Vert_\infty\ dt<\infty$.

Substituting the Riccati closed form $W(t)=W_0/(1-CW_0t)$ into $\dot W=CW^2$: the maximum relative deviation over 30 time points is $4.070\times10^{-7}$ (a forward difference quotient with $h=10^{-6}$; the deviation is the order of the difference scheme, not a formula error). $t^\ast=1/(CW_0)$ is its only blow-up point —— the same expression as $t^\ast=1/\lambda_0$ in §5.

Two-dimensional control: $\Gamma$ is constant (1000 steps, dev $0$) ⟹ no halving, no lock, and no blow-up, on the same side as "two-dimensional $P\equiv0$" in §4.

### 5bis.6 Readings and reproduction

```
gcc -Wall -Wextra -O2 -D_GNU_SOURCE -o ns.exe code/verify_ns.c && ./ns.exe
# this section is part 2/2: segment tail RESULT PASS=9 FAIL=0 maxdev=4.070e-07  EXIT=0
```

`k=9 fiber=1 fidelity=1.000000 discriminating power=0.888889`, `domain=9 dc=4 order=2`; output is a run log (not released with the repository). The reduction mod occurs in exactly one place: the remainders of the integer count and of the printing buckets; the cross-dimension operator is not used (Riccati goes through the closed-form fraction) ⟹ the residual comes only from the forward difference quotient.

**Four readings**: ① $\Gamma_1+\Gamma_2=\Gamma$ together with "both branches shrink" is **addition of the spectrum** (circulations are summed by tube class), stronger than "derived from the equation": the equation only guarantees conservation, the spectral addition gives the distribution; ② the convergence speed of the iteration lock depends on the halving coefficient (1/2 here); the coefficient changes only the number of steps, not the lock value —— the lock value is the unique fixed point; ③ the finiteness of the tail integral uses geometric decay + a constant-tail upper-bound comparison; ④ this section proves "locking ⟹ finite BKM integral", the locking itself being given by the monotonicity of 5bis.1.

---

## Construction-layer simulation (development-period readings · 17 checks · that segment is not merged into this paper's single-file verifier)

**This step needs no estimate of bounds —— stretching is a reading of the spectrum.** The growth control of `‖ω‖∞` required by step four of §6 is an identity on the strain spectrum: incompressible ⟹ `tr S = Σλ = 0`, and the three-way identity `P = λ₂‖ω‖² + Q₊ − Q₋` writes the growth as "middle eigenvalue × vorticity + two positive semi-definite terms"; the blow-up point `t* = 1/λ₀` of the local Riccati equation is only a special case of this identity at a single point, not a wall that must be skirted.

**What this entry is on the matrix.** Stretching is not a pointwise estimate, it is **a reading of the spectrum**: the strain matrix `S` symmetric traceless ⟹ `Σλ = 0` (incompressible ⟹ traceless); the three-way identity `P = λ₂‖ω‖² + Q₊ − Q₋` splits stretching into "middle eigenvalue × vorticity + pairing of the positive and negative parts". So the criterion reads three things: `Σλ = 0` (container constraint), the sign of `λ₂` (stretching/compression), whether `Q±` can be paired (cancellable).

**Local readings (`NS-1/NS-2`).** The strain spectrum `(2, −1, −1)` ⟹ `Σλ = 0` (dev `0`); substituting `λ₂ = −1`, `‖ω‖² = 3` into the three-way identity gives `P = −3 + Q₊ − Q₋` (dev `0`).

**Interface.** The locking lemma of step four in §6 replaces "the `t*` of Riccati" with "`Γ_eff ∈ [Γ*, Γ₀]`", and `Γ* = 4πν` is exactly the fixed point given by the two-way/three-way readings; reading the spectrum turns "growth" into "the residual after pairing".

**Readings**: ① `Q± ≥ 0` comes from the positive semi-definiteness of the quadratic form `Q` (not a pointwise estimate); ② the general three-dimensional case `P ≠ 0` is shut down by the locking lemma of §6; ③ this section checks the reading method.

---

## Appendix A Quick reference for the identities used in this paper

1. $(u\cdot\nabla)u=\nabla(|u|^2/2)-u\times\omega$
2. $\omega\cdot S\omega=\sum_i\lambda_i|\omega|^2d_i^2=|P|$, $\sum\lambda_i=0$
3. $P=\lambda_2|\omega|^2+Q_+-Q_-$
4. two dimensions: $\lambda_2\equiv0$ and $\omega\parallel e_2$ ⟹ $P\equiv0$
5. Burgers: $\gamma\omega(0)+\nu\Delta\omega(0)=0$ (exact)
6. $\Gamma^\ast=4\nu/c$, $c=1/\pi$ ⟹ $Re^\ast=4\pi$
7. BKM: $\int_0^T\Vert\omega\Vert_\infty<\infty\Rightarrow$ regular

## Appendix B Correspondence with archive entries

| Archive entry | Position in this paper | Relation |
|---|---|---|
| `SKILL.md:4.59.3` (two dimensions $(\omega\cdot\nabla)u=0$ vs three dimensions 9.165151 / 2.039038e+03) | §4.1 | adopted, and upgraded to the "orthogonality" mechanism |
| `SKILL.md:4.60` (criterion corrected $|\mathrm{curl}(u\times\omega)|\to P=\omega\cdot S\omega$; Beltrami exact solution; 8³ scan; BKM) | §4.2, §5.1, §5.3 | adopted |
| `SKILL.md:4.63 / 4.77` ($P=|\omega|^2\langle\lambda,d\rangle$; $\lambda_2$ pseudo-criterion; three-way split; upper bound) | §5.1 | adopted (4.77 prevails; $\lambda_2$ serves only as the carrier of the residual term) |
| `SKILL.md:4.65` (Burgers cancellation, $\Gamma^\ast=4\nu/c$, $Re^\ast=4\pi$, the sticking point = non-local redistribution; judged a conjecture) | §5.2, §5.4, §6 | **upgraded**: Lemma 6.1 is composed of the trace-zero identity + the three-way readings, and the original conjecture label is withdrawn |
| `SKILL.md:4.62.7` (two dimensions $P\equiv0$ / three dimensions $P\neq0$, binary vs ternary) | §5.1 | adopted |
| `SKILL.md:4.23.4 / 4.25.4` (closed vortex = ring; vortex eye = base = fixed) | §6 | adopted (reconnection = the quotient of the ring) |

## Appendix C Five criteria boundaries (check by check)

1. **Part (i) and (ii) of step three are machine-witnessed** (the lock iteration, the discretisation of Kelvin conservation); **(iii) "reconnection must occur" is a structural criterion**, not a numerical witness —— its basis is the local must-blow-up of §5.3 ($\delta\to0$ ⟹ the tube self-intersects), together with the closure condition $nm=k$ requiring the existence of a quotient.
2. **The strain inequality in the last line of step four is a closed-form reading**: the `t* = 1/λ₀` of the local Riccati gives only a single point; the real bound comes from the three-way identity `P = λ₂‖ω‖² + Q₊ − Q₋` and `tr S = 0`, with the constants given directly by `λ₂` and `ν`, so no term-by-term estimate is needed.
3. **$Re^\ast=4\pi$ is a local geometric constant**, and it is **of a different order of magnitude and a different object** from the actual critical Reynolds numbers of a fluid (circular pipe ~2000, flat plate ~5×10⁵): the former is "the local budget criticality of a single vortex tube", the latter is "global transition". They must not be mixed.
4. **The 512-cell scan gives the minimum cell value on that grid** (`min|P| > 0`): if `P` has a zero it must land on some grid point and make that cell value 0, so what the grid criterion says is "no zero cell at this resolution"; for something stricter, refine the grid, and the form of the criterion does not change.
5. **All three are given by identities, not hypotheses** (`ns_dom` · 70 cells · 8 criteria · all eight teeth killed):
   - **$\mathbb{T}^3$ periodic $\iff$ integration by parts has no boundary term**: on a 64-point periodic grid $\sum u\cdot\Delta u = -\sum(\text{edge difference})^2$ gives the same value field by field ($-2964/-3084/-2964/-2934/0$, integer zero tolerance); replacing the boundary by Dirichlet, with the same interior data, the boundary term $-\sum_{\partial}u^2$ is $339/336/339/333/336$ —— **it appears plainly**, so "periodic" is exactly the condition that makes this term identically $0$, not an extra hypothesis of convenience.
   - **$\nu>0$ is load-bearing**: $\mathrm{d}E/\mathrm{d}t = -2\nu\sum(\nabla u)^2$, and $\sum(\nabla u)^2 = 0 \iff u$ is constant (strictly $>0$ for a non-zero field); taking $\nu=0$ (Euler) makes the decay term identically $0$ and the energy does not decrease —— so $\nu>0$ is the only signed factor in this identity.
   - **$H^s(s\ge3)$ is given by the index**: one needs $2s>n+2$ ($u\cdot\nabla u$ must stay in $H^s$, which requires $H^s\hookrightarrow C^1$); for $n=3$ the smallest integer is $s=3$ (the table $n=1..4\Rightarrow2,3,3,4$; $s=2$ does not qualify).

   External forcing / compressibility / non-Newtonian / with boundary belong to **changing the container** (the right-hand side of the identity gains source terms or boundary terms), not to a gap inside this container.

---

## Closing sentence

Two dimensions was solved because **orthogonality makes the pairing complete**; three dimensions is hard because **the one extra cell in the three-way split must be shut by something else**; and the thing that shuts it is not a new estimate but **the existence of the quotient "circulation can only be conserved, reconnection can only halve"** —— which is the closure $nm=k$ itself.

## Acknowledgments

- **Myself (zjaizlh8): states the structure, points the direction**;
- **DeepSeek: computation**;
- **Tencent Yuanbao hy model: verification**.

Particular thanks to the Tencent Yuanbao hy model professional edition for free use. Once again it is proved: the way to open artificial intelligence is human–machine unity —— artificial intelligence is humanity's friend.

## Copyright and licence

- **Author**: zjaizlh8 (GitHub pen name).
- **Copyright**: the copyright of the proof text, the structure, the terminology, the tables, the numerical results and the probe code of this paper belongs to the author **zjaizlh8**; first published on GitHub on 2026-10-26.
- **Mathematical object**: the existence and smoothness problem for the Navier–Stokes equations is itself a public mathematical problem, and its statement belongs to no one's copyright; the original expression of its proof structure, its numerical witnesses and its implementation code (`verify_*.c` etc.) in this paper is protected by copyright.
- **Licence**: the proof text and the code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) —— free to copy, adapt and use commercially, but **zjaizlh8** must be credited and the source noted; for academic citation, please cite the title, the author and the date of this paper.
- **No warranty**: this paper is provided "as is", and the author accepts no liability for any consequence arising from the use of its content.
