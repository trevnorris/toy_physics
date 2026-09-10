# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v3 (round 3, DOCUMENT branch)

VERDICT: NOT-SOUND (3 must-fix + 4 nits); structural layer confirmed right. Verbatim (1-line trace preamble).

---

I'll form an independent view from the governing sources first, then read the v3 spec against that view. Starting with the decision list, c2 record, house-format precedent, N6 disposition, and the profile-class consult.Next I’ll read the c2 house-format spec, the N6 disposition, and the profile-class consult, still without opening the S11c-d spec.The consult file was the prompt, not the three-way outcome. I’ll pull the consult results and the S11c-a/S11c-b physics that S11c-d actually inherits.I’ll stay on the sources of truth: S11c-a power counting and N6, the consults as context only, and the real c2 export keys.I’ll inspect the real c2 export keys (names only) and the conversion/confinement carry-ins, then freeze my independent picture before opening the spec.Independent picture is in place: localized interface, Born in `η` with live `σ_W`, honest c2 import, two kill channels, and `F'(0)` only for the strong edge. Opening the spec now.I’ll check the order bookkeeping, the `sin²` strong-edge counterexample, and the zero-transfer Fourier identity with a short symbolic calculation.The order check and the coupler counterexample match the spec; the form-factor kinematics do not. Writing the review from that.

# Independent physics review — S11c-d SHARED PHYSICS spec (v3)

Artifact: `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`

I formed the required object from `S11c_decisions.md` (N5/N6/N7/N10–N15), the c2 record/disposition, S11c-a bookkeeping, and the real c2 export keys (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, plus bind-closure; **no** increment). Then I read the spec.

**Independent picture of what S11c-d must be.** Localized **thickness** interface (`Δw₁ ≠ 0`), Born in contrast `η`, with `σ_W` and the wave/edge kinematic (`k L_W`) and the form-factor argument (`Q_n L_W`) kept live and not identified with each other. Consume both closed exports as-is. Emit the linear mixing, two distinct N13 channels as computed objects, and both `O(ε²η²)` absolute and `O(η²)` fractional leakage labels on `ΔA`. Strong edge is a **new** construction; this step only supplies weak coefficients. c2 operand DEBT stays unclosed and material. N6 here is insertion-route + one-sided atom corruption, not `∇W→0`.

---

## Must-fix

### F1. Class membership requires a modulus jump `Δm₁ ≠ 0`. That is not the localized-interface object, and it contradicts the spec’s own vertex.

**Spec (§0.1 and §1c):**

> NAME the profile class = a **localized interface** … required on **both** the thickness `w₁` (`Δw₁ ≠ 0`, so the asymptotes `W₋ ≠ W₊`, §2) **and** the modulus `m₁` (`Δm₁ ≠ 0`, the mixer)

> The **interface** condition … is required on **both**: `Δw₁ ≡ w₁(+∞) − w₁(−∞) ≠ 0` … **and** `Δm₁ ≡ m₁(+∞) − m₁(−∞) ≠ 0` (the modulus interface, the mixer).

**Why it is wrong.** N5 names a profile class for the varying thickness. S11c-a already split the two profiles and forbade tying them:

> `w₁` and `m₁` are dimensionless `O(1)` profile functions; `L_W` is an independent length. … `η` and `σ_W` are varied independently … no engine may replace `σ_W` by `η`  
> (`S11c_a_SHARED_PHYSICS.md` §2a)

The spec itself quotes that independence (`§1c`: “⛔ no engine may tie them, ⛔ no `m₁ = m₁[w₁]`”) and then, in `§2`, forbids identifying the object with the modulus subchannel:

> Constitutive mixing is *driven by* `∇μ_R,bg ≠ 0` … ⛔ but this is **not** the whole support — a `μ_R`-only projection **drops** the tilt and advection channels §1a requires and is an **ablation**, not the object.

A thickness interface with independent `m₁` (constant `μ_R`, or a modulus *bump* `Δm₁ = 0`) is still an interface: `W₋ ≠ W₊`, two-asymptote `L±` still exist, and the imported vertex still mixes through tilt `∇w₁` and N4 advection. Requiring `Δm₁ ≠ 0` as **class membership** drops that case — including the slit-like object “`W` jumps, `μ_R` need not.” It also makes `§5c`’s “bump ⇔ `Δm₁ = 0`” mis-label a thickness interface as a bump.

This is not “wrong on a different input.” It changes which profiles are in-class and which mixing channels the class is allowed to see.

**Minimal fix.** Require the interface condition on **`w₁`** (`Δw₁ ≠ 0` ⇒ `W₋ ≠ W₊`). Keep `m₁` independent. Use `Δm₁` / `m̂₁′(0)` as a **per-profile** discriminant of the modulus subchannel, not as a gate on the class. A bump is `Δw₁ = 0` (and analogously `Δm₁ = 0` for the modulus control), not “either jump vanished.”

---

### F2. `§1d` misidentifies the form-factor argument `Q_n L_W` with WKB/sudden kinematics, and calls `Q_n L_W → 0` “maximal conversion.”

**Spec (§1d):**

> `Q_nL_W → 0` is the **zero-transfer / sudden** limit, where `m̂₁′(Q_nL_W) → Δm₁` (the integrated jump, maximal conversion); `|Q_nL_W| ≫ 1` (reached by `L_W → ∞`, i.e. `σ_W → 0`) is the **WKB / adiabatic** regime, where the form factor is suppressed (Riemann–Lebesgue) and conversion is small — the class N5 says does not match a sharp edge.

Same section correctly forbids an extra `σ_W → 0` / `L_W → ∞` expansion of the first-shape-order operator, and correctly forbids `η → O(1)` inside it.

**Why it is wrong.** Three different things are glued together.

**(i) The identity is right; “maximal conversion” is not.**  
If `m₁′` is integrable,

\[
\widehat{m₁′}(Q_n L_W)=\int m₁′(ξ)\,e^{-i Q_n L_W ξ}\,dξ,\qquad \widehat{m₁′}(0)=Δm₁.
\]

That is zero-transfer, not a conversion maximum. For a one-sign (monotone) `m₁′`, \(|\widehat{m₁′}(q)|\le Δm₁\). The spec’s own class includes overshoots (`§3b`: “monotone steps **and** profiles with localized overshoots/wells”). For those, \(|F|\) need not peak at `Q_n=0`. Algebraic counterexample, still an interface (`Δm₁=1`):

\[
m₁′(ξ)=δ(ξ)+δ(ξ-L)-δ(ξ+L),\qquad F(q)=1-2i\sin(qL),\qquad |F|^2=1+4\sin^2(qL).
\]

\(|F(0)|=1\), \(|F(π/(2L))|=\sqrt5>1\). Conversion also depends on `v_out`, channel density, and the rest of the vertex, not only `|F|`.

**(ii) `|Q_n L_W|≫1` is not WKB.**  
Riemann–Lebesgue is the tail of the form factor in **momentum transfer**. WKB/adiabaticity is whether the **wave** is slow on the edge, i.e. `k L_W ≫ 1` (equivalently a further `σ_W→0` expansion). They part:

- A slow interface still has \(\widehat{m₁′}(0)=Δm₁\). Zero-transfer is **not** suppressed by taking `L_W→∞`.
- `|Q_n L_W|≫1` at **finite** `σ_W` is just large-angle / large-transfer Born on a finite-width edge — in-class, not “the WKB class N5 rejects.”
- `|Q_n L_W|≫1` “reached by `L_W→∞`, i.e. `σ_W→0`” is exactly the extra `σ_W→0` operator limit the same paragraph forbids.

**(iii) “Sudden” `Q_n L_W→0` is not `L_W→0`.**  
On the named homotopy (`λ≡η`, **fixed** `L_W`) , `Q_n L_W→0` is small transfer at fixed sharpness — in-domain. `L_W→0` at fixed `η` is `σ_W=η W̄₀/L_W→∞`, which leaves the **first-`σ_W`** truncation the same way `η→O(1)` leaves the first-`η` truncation. The imported kernel has no `∇²W, ∇³W, …`. That is not a target limit.

S11c-a already split contrast from first-jet (`η` vs `σ_W`; `:189-198`) and treated wave amplitude `ε` as a third bookkeeper. The optical sharpness-on-the-wave is `k L_W`, not `Q_n L_W`. `Q_n L_W` should stay live as the **form-factor argument** (do not Taylor-expand `m̂₁′`). It is not the WKB axis.

An engine that implements `§1d` as written can (a) drop large-`Q_n` conversion as “WKB / out of class,” (b) treat `Q_n=0` as the conversion maximum, or (c) take `L_W→0` or `L_W→∞` as the sudden/WKB evaluations of a first-jet operator.

**Minimal fix.** Keep three independent live quantities: `η`, `σ_W`, and kinematics. Split the last:

- `Q_n L_W` = form-factor argument; `m̂₁′(0)=Δm₁` is zero-transfer, not “sudden” and not “maximal conversion.” `|Q_n L_W|≫1` is the form-factor tail at finite `σ_W`, still in-class.
- `k L_W` (or equivalent) = sharpness on the wave. WKB/adiabatic is `k L_W≫1` / an extra `σ_W→0` expansion — forbidden as an additional expansion of the imported operator, not identified with large `Q_n L_W`.
- Do not take `L_W→0` at fixed `η` (`σ_W→∞`).

---

### F3. The bound-pole solve is exempted from re-expansion, but the “stated error limitation” is never stated. That lets a truncated-model pole be claimed as the N13 channel.

**Spec (§2):**

> ⚠ **Exception:** a **bound-pole** spectral solve (`§3b`) is a **resummation** `G = (G₀⁻¹ − V)⁻¹` that the first-order continuum re-expansion **cannot** create — that spectral solve is **exempt** from the continuum re-expansion and carries its own stated error limitation.

`§3b` then binds this object to N13 (bound capture kills the photon as bulk radiation does) and asks for Jost/Evans zeros of the imported operator, with computed absence permitted and the 1D-well theorem correctly refused.

**Why it is wrong.** The limitation is named and not given. What must be stated is the order conflict:

- Continuum Born: `ΔA=O(ελ)`, omitted operator pieces `O(λ²)` ⇒ relative error `O(λ)`, controlled as `λ→0`.
- A weak bound (when one exists) sits at `E\sim λ²`. The imported operator is only first shape order. Omitted `O(λ²)` pieces compete with the whole binding. Existence and location are **not** controlled parent-theory predictions; they are spectrum of the **truncated** model.

`§3b` already forbids a class-wide existence claim and the 1D-well theorem (right: unequal asymptotes, multi-component, `ω`-dependent, possibly non-Hermitian DtN). That does not replace the truncation caveat. Without it, `S11CD_BOUND_MODE_SPECTRAL_TEST` / `S11CD_CONFINEMENT_CONDITION` may be claimed as N13 physics rather than truncated-model spectral data.

**Minimal fix.** State the limitation next to the exemption: Evans/Jost zeros are spectrum of the first-shape-order imported operator; omitted `O(η²,σ_W²)` terms are the same order as weak binding; emit the zeros (including empty) as truncated-model data, ⛔ not as a controlled parent-theory bound channel. Keep the non-Hermitian / physical-sheet tests.

---

## Items that are right (no finding)

**Order bookkeeping (`§3c`, N12).** Symbolic check (linear, homotopy `λ=η` at fixed `L_W`):

```
A_H = K1*epsilon*lambda
J_in = epsilon**2
J_conv = |K1|**2 * epsilon**2 * lambda**2
C = |K1|**2 * lambda**2
```

Converted amplitude `O(ελ)`, absolute flux `O(ε²λ²)`, incident `O(ε²)`, fraction `O(λ²)=O(η²)` with `ε²` cancelled. Emitting both labels, attaching N12 to `ΔA` not to total `A` unless `A_0` is computed zero, and refusing to read `O(ε²λ²)` as the nonlinear program, are correct. `A_0` computed rather than assumed zero matches withdrawn F.

**Strong-edge bridge (`§3d`, N7).** Coupler check:

```
A_H = -I*epsilon*sin(G*eta)
C = sin(G*eta)**2 = G**2*eta**2 - G**4*eta**4/3 + …
C'(0) = 0
C''(0)/2 = G**2 = lim C/eta**2
∂_η(A_H/ε)|_0 = -I*G
C(π/G) = 0  while  G**2 ≠ 0
```

Splitting amplitude slope from fraction coefficient `½C''(0)` is required (`C'(0)=0`). Nonzero Born coefficient is not a lower bound on strong conversion. Order-unity edge is a **new** construction, not a reduction. `O(1)` / lab number withheld. Conditional on a computed-nonzero weak coefficient. This section is honest.

**c2 import (`§1a–§1b`).** Real write-keys; increment/term-origins/§3d rows correctly absent from the consume-set (70-key delta: two closed operators + bind-closure). Per-engine SOUND values vs cross-engine N6 thread only; matched zeros `(0)−(0)`; `R_N6=18/288` schema-unmatched preserved with `R_cov` no-nonzero; operand DEBT 40/76/18 UNADJUDICATED and material; leftover SHAPE not “just thickness”; F/G withdrawn; sign conventions, six §3d items, three Φ/`V`/omitted-block caveats named. No quiet upgrade.

**N6 (`§5a–§5b`).** Insertion-route + one-sided **profile-atom** corruption (tilt `w₁′`, advection on `RHOBR_CONSTANT`, computed absence / ⛔ `A−A` on `RHO4_CONSTANT`). `∇W_bg→0` / `η→0` rejected as corruption. Anchoring not the test. No `Φ` on amplitudes. Uniform regression improved: jets→0 with live asymptotes, coupling computed not asserted to vanish. Kernel-level Eulerian↔material correctly left as c2’s unclosed debt.

**Two kill channels (`§3b`, N13).** Continuum vs bound are distinct emitted objects; bound capture not folded into `J_conv`; confinement is a computed (possibly empty) object; not a Bloch band; 1D-well theorem not invoked. This is the right N13 structure **once F3’s truncation limitation is actually written**.

**House template / N11 / chain.** Pointer inheritance of the c1/S11b rest-frame domain (`|q v/ω|≪1` + boundary-layer/subsonic; large `k c_s0/|ω|` necessary not sufficient); no convective operator; `v_bulk_normal_0` not `v_0`. Fold topology, digest pin including `ledger_fold.py` + c1 + c2, T7 comparator on the mixing amplitude surfacing the DEBT, blind WL, F9/N14 reservations, supplied-vs-computed split: present, not weakened. FORM computed here is an N2-permitted boundary refinement and is not e’s withheld `O(1)` target.

---

## Nits (do not change the verdict by themselves)

1. **`§5d` `∝ k·a`.** N7 already names that structure, and the same paragraph says the FORM is the **full** vertex. Still an expected-shape phrase an engine can type instead of projecting the imported kernel. Prefer “emit the computed projection; do not substitute a typed `k·a`.”
2. **`§3c` flux bilinear** “S11b quadratic energy / c2 traction–slab pairing.” c2’s pairing row is EMIT-only, not in the consume-set (`§1a`). Define `J` from the S11b energy current / the imported closed operator, not from a dropped c2 tag.
3. **`§5a` route 2** “flattening the **interface faces** to material in-plane coordinates” mixes the in-plane `n̂`-edge with S11c-a’s slab-face `w′` map. The atom corruption is the constructible control; route 2 needs the actual chart named (in-plane `x=X+u` at fixed anchoring, or the `w′` face flattening), not “interface faces.”
4. **`§3d` “form-factor node `Q_n L_W = nπ`”** is a box/`sinc` pattern, not class-generic (`tanh`/sech form factors have no real nodes). Keep “the coefficient can vanish (computed node or `k·a` selection)” without `nπ`.

---

## Verdict

**NOT-SOUND.**

Must-fix: **F1** (`Δm₁ ≠ 0` as class membership), **F2** (`Q_n L_W` identified with WKB/sudden and “maximal conversion”), **F3** (bound-pole truncation limitation referenced, not stated).

Those three change which profiles are built, which kinematic limits are treated as in-class, and what the spec may claim about N13. Order bookkeeping, the strong-edge bridge, c2-import honesty, and the N6 insertion control do not.
