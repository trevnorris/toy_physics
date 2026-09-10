# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v5 (round 5, DOCUMENT branch)

VERDICT: NOT-SOUND (1 must-fix F1 = the v5 Fourier fold; 2 nits). Verbatim (1-line trace preamble).

---

I'll start from the governing sources and form an independent view of S11c-d before opening the spec. Reading the review-legs skill and the required source documents first.I have the governing requirements. Next I’ll pull S11c-a power counting, c2 export objects, and related physics so I can form an independent view before opening the spec.I have the requirements. Next I’ll pull S11c-a/S11c-b physics, c2 export keys, and the consult replies so the spec is judged against those, not against itself.I’ll inspect the real c2 export keys and the consult replies next, then lock in an independent physics view before reading the spec itself.I have an independent picture of the required objects. Next I’ll read the full spec and the real c2 export keys, without treating the v5 fold as already correct.I’ll check the real c2 export key set and the parent Fourier convention so the v5 fold can be judged against the imported objects, not against the spec’s wording.I’ll do a mechanical check for whether the imported kernel actually carries a `DiracDelta` identity, since the v5 fold tells engines to reduce against that piece.I confirmed the consumed c2 kernel has no `DiracDelta`, while that identity lives on the c1 DtN kernel. The rest of the review is against the sources, with a short order-counting derivation.

---

# Independent physics review — S11c-d SHARED PHYSICS (the file still titles itself v4)

**Artifact:** `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`  
**Role:** physics spec both blind engines will read.  
**Sources read first:** `S11c_decisions.md` (N5–N7, N10–N15, S11c-d row), `S11c_c2_self_energy_fold.md`, `S11c_c2_SHARED_PHYSICS.md`, `_measurements/S11c_c2_N6_reconcile_disposition.md`, S11c-a bookkeeping, real c2 export keys, c1 `dtn_kernel` display, c2 `kernel_bridge`/`kernel_apply`. Consult used only as context.

**Verdict: NOT-SOUND.** One must-fix (the v5 Fourier fold). Everything else I checked does not change what an engine computes or what the spec may claim.

---

## Must-fix

### F1 — The 3-D→1-D carrier reduction is not computable from the named import, and the parent Fourier convention is still supplied

**Spec location**

> ⚠⚠ **The imported c2 carrier's normalization is NOT supplied — it must be COMPUTED, ⛔ not asserted.** The c2 engine's transform convention is an **unnormalised forward transform with a normalised inverse** (`c1` uses `DiracDelta(k−k′)` with no `(2π)³` coefficient; the `(2π)⁻³` sits on the application/inverse — the c2 self-energy fold's own stated convention, `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`), so `s11cc2FourierW1ProfileHatTransfer` / `s11cc2FourierW1ProfileJetHat{i}` are ⛔ **not** equal to the normalized `(2π)⁻³∫d³y e^{−iQ·y}(·)`. ⛔ **Do not supply a numeric `(2π)`/`δ²(Q_∥)` map for the imported carrier** … ⇒ **each engine COMPUTES the 3-D→1-D reduction factor** by reducing the one-profile insertion **as it actually appears in `s11cc2ClosedCouplingKernel`** against the imported kernel's own flat identity (its `DiracDelta³(k_out−k_in)` piece) … so the tangential `δ²(Q_∥)`, the `2π`, and the dimensional content **fall out of that computation**
>
> — §1c

> **SUPPLIED …** the §1c localized **thickness** interface … short-range-jet domain, **exact reduced Fourier/3-D carrier convention**, and admissibility
>
> — §8

**Why it is wrong**

The fold’s intent is right: do not type a `[L_W/(2π)]δ²` (or `(2π)⁻³∫`) identity for the imported carrier. The method it names is not available on the objects S11c-d actually consumes, and the surrounding sentences still hand both engines a convention that fixes the factor.

1. **The named identity piece is not in the consumed kernel.** Mechanical existence check: `DiracDelta` does not occur in `scripts/S11c_c2_exports.py` (0 matches). The consume-set objects `s11cc2ClosedCouplingKernel` and `s11cc2ClosedSlabOperator` therefore have no `DiracDelta³(k_out−k_in)` flat identity to reduce against.

   The three-delta identity lives on the **c1** two-momentum kernel, not on c2’s closed coupling kernel. From `S11c_c1_exports.py` `dtn_kernel` display:

   > `FLAT_DIAGONAL, omega*rho_m*DiracDelta(-s11cc1_k_input_1 + s11cc1_k_output_1)*DiracDelta(…_2…)*DiracDelta(…_3…)/s11cc1_q_out_output`

   c2 then **strips** those deltas before folding. `S11c_c2_selfenergy_fold_sympy_audit.py`:

   ```text
   deltas = diagonal.atoms(sp.DiracDelta)
   z0out = diagonal.xreplace({d: sp.S.One for d in deltas})
   ```

   So “the insertion as it actually appears in `s11cc2ClosedCouplingKernel`” is an opaque `Function` of `Q`, and “its `DiracDelta³` piece” is a piece of a **different** import (and not of the c2 kernel at all). An engine that follows §1c literally cannot execute the computation.

2. **Even `dtn_kernel` does not encode the 1-D reduction.** c1’s `FIRST_SHAPE` carries `s11cc1_w1_profile_hat_transfer` as a general 3-D Fourier symbol. The interface restriction `f=f(n̂·y/L_W)` and the factor relating `δ³(Q)` to `δ²(Q_∥)δ(Q_n)` are S11c-d geometry plus a Fourier convention, not a coefficient sitting next to a delta in `s11cc2ClosedCouplingKernel`. Pointing at “the imported kernel’s DiracDelta³ piece” does not determine `L_W`, `δ²(Q_∥)`, or the power of `2π`.

3. **The parent convention is still supplied, from one engine’s script, and it is not shared physics.** Neither `S11c_c1_SHARED_PHYSICS.md` nor `S11c_c2_SHARED_PHYSICS.md` states a `(2π)` / `DiracDelta` convention (no `2π` / `DiracDelta` in those specs). The unnormalised-forward / `(2π)⁻³`-on-inverse sentence is copied from a comment in the **SymPy** audit:

   ```text
   # c1 uses DiracDelta(k-k') without a (2*pi)^3 coefficient; the transform
   # convention has an unnormalised forward transform and normalised inverse.
   p0 = integral(phase0 * diagonal * local_source / (2*sp.pi)**3, *limits0)
   ```

   That is one engine’s application convention, not a sibling-spec object. Putting it in the shared spec (i) leaks a numeric convention the fold claimed to withhold — together with the supplied `f̂_red` and `∂_{y_i}f = n̂_i f′/L_W`, it uniquely determines the 3-D→1-D factor without looking at any kernel — and (ii) points the blind Wolfram engine at the other engine’s construction script. There is **no** WL self-energy engine (`S11c_c2_self_energy_fold.md`: “no WL self-energy engine and no self-energy comparator”); WL must re-derive from specs. A spec-level pointer into `S11c_c2_selfenergy_fold_sympy_audit.py` is the opposite of that control.

4. **§8 contradicts §1c.** §1c says the c2 carrier normalization is **not** supplied and must be computed. §8 lists “exact reduced Fourier/3-D carrier convention” as **SUPPLIED** (unfalsifiable in this build). That is the inventory engines use to decide what they may type.

This changes what engines compute: they cannot follow the named kernel identity; they will either invent incommensurate reductions or both reconstruct the factor from the leftover SymPy convention and share it. A wrong shared `(2π)` power is exactly the defect the fold was trying to prevent (“a wrong one (e.g. off by `(2π)³`) is a defect **both engines would share**”).

**Minimal fix**

- Keep `f̂_red` as S11c-d’s **supplied** 1-D convention (that part of §1c is a named object, not an answer).
- Delete the unnormalised/normalised/`(2π)⁻³` sentence and the pointer to `S11c_c2_selfenergy_fold_sympy_audit.py`.
- Do not tell engines to read a `DiracDelta³` piece out of `s11cc2ClosedCouplingKernel` (it is not there).
- Require each engine to **emit** the 3-D→1-D reduction of its **own** closed-kernel Fourier symbols as a computed object with both operands: (i) the 3-D carrier as it actually appears in that engine’s closed coupling kernel, (ii) the reduced 1-D kernel obtained by applying **that same engine’s** Fourier convention — the convention realized in its construction of the two-momentum identity / profile insertion, not a typed `[L_W/(2π)]` or `(2π)² L_W` — to the §1c geometry `f=f(n̂·y/L_W)`. The comparator joins those reduced kernels.
- If the AGREE’d c1 `dtn_kernel` flat identity is used as the convention witness, name **`dtn_kernel`**, not `s11cc2ClosedCouplingKernel`, and still require the 1-D geometry step to be computed.
- In §8: `f̂_red` is SUPPLIED; the 3-D→1-D factor and reconstruction are COMPUTED.

---

## Items that are not findings (checked, then filtered)

### 1. Profile class / regime — sound

Naming a **localized thickness interface** (`Δw₁ ≠ 0`, hence `W₋ ≠ W₊`; short-range retained jets; unequal-asymptote `f` not treated as ordinary `L¹`) is a correct non-global `N5` object. The bump (`Δw₁ = 0`) is kept as a discriminant, not silently substituted for an edge (`∫ w₁′ = Δw₁` when `w₁′ ∈ L¹`). Independent `m₁` (constant / bump / interface) is right: S11c-a already split `w₁` and `m₁`.

The three axes are independent and live:

```text
η          contrast / Born strength
σ_W        first-jet sharpness
s=Q_n L_W, k_a L_W   kinematics (not grades)
```

Forbidding an extra `σ_W→0` expansion / Taylor-in-`s` / `L_W→∞` is the right “do not kill the edge form factor” control. Forbidding `η → O(1)` **inside the first-shape-order imported operator** is also right: those exports are truncated at first shape order in each of `η` and `σ_W` (S11c-a §2a); using them at `η = O(1)` is an incomplete operator, not a non-perturbative solution. `L_W→0` at fixed `η` (`σ_W→∞`) is correctly excluded as leaving the retained model. Periodic→Bloch and WKB are correctly named as **other** `N5` classes, not this S-matrix evaluated on a periodic `w₁`.

### 2. Order bookkeeping — sound (derivation)

Linear theory, homotopy `λ ≡ η` at fixed `L_W` (so `σ_W = κλ` tracks, grades not identified):

\[
A_H(λ)=ε\bigl(a_0+λ a_1+λ^2 a_2+\cdots\bigr),\qquad
J_{H,\mathrm{out}}=B_H(λ)[A_H,A_H],\qquad
J_{T,\mathrm{in}}=ε^2\bigl(j_0+λ j_1+\cdots\bigr).
\]

Quadratic expansion:

\[
\begin{aligned}
J_H^{(0)}&=ε^2 B_0[a_0,a_0],\\
J_H^{(1)}&=ε^2\{B_0[a_0,a_1]+B_0[a_1,a_0]+B_1[a_0,a_0]\},\\
J_H^{(2)}&=ε^2\{B_0[a_1,a_1]+\cdots\}.
\end{aligned}
\]

- If the computed uniform baseline `a_0=0` (and the `B_1[a_0,a_0]` slot dies with it), converted flux starts at \(O(ε^2λ^2)\) and the fraction \(C=J_H/J_{T,\mathrm{in}}\) is \(O(λ^2)=O(η^2)\): the \(ε^2\) cancels in a **linear** theory.
- If `a_0≠0`, there is an \(O(ε^2)\) baseline and an \(O(ε^2λ)\) interference slot; `N12`’s “leakage \(O(ε^2η^2)\)” is then not the leading physical conversion.

The spec emits both the flux object and the fraction, attaches the `N12` labels **only conditionally** on the computed `K_0`/`A_0` disposition, and does not type `K_0=0` from withdrawn F / S11b. That is the honest reading of `N12` after F was withdrawn (`M3`: S11b decoupling is an oracle, not a premise). The induced-field quadratic \(B[\Delta A,\Delta A]\) is separately named so it is not confused with physical \(J_H[A]-J_H[A_0]\). An \(O(ε^2λ^2)\) quadratic of a **linear** vertex is not the excluded nonlinear-intensity program (`N10`/`N12`). No mis-ordered term is labelled as that program.

### 3. Strong-edge bridge — sound

S11c-d emits weak Taylor data (`∂_λ(ΔA/ε)|_0`, total-fraction coefficients, induced-field `λ²` coefficient) and names `C_{\mathrm{strong}}(1)` as what the lab bounds. The coupler counterexample is stated as **not** a slab model and shows \(C=\sin^2(ηG)=0\) at finite `ηG=nπ` with the same baseline-free `N12` orders — so a nonzero Born coefficient is not a lower bound. The order-unity edge is out of scope and a **new** construction (piecewise-uniform / matched / effective interface with explicit undetermined parameters), not a reduction of the truncated kernel. Factorization \(C_{\mathrm{edge}}=C_{\mathrm{interior}} F(ω,ϑ)\) is correctly not assumed. Nothing here is over-claimed.

### 4. c2 import honesty — sound

Per-engine SOUND is restricted to the closed-operator / closed-kernel **values**. Cross-engine content is the N6 covariance thread only (Reading B); matched zeros are `(0)−(0)`, not operand AGREE. Operand DEBT carrier 40 / source 76 / Φ 18 is UNADJUDICATED and **material**. Both `R_N6 = 18/288` (raw census, schema-unmatched) and `R_{\mathrm{cov}}` no-nonzero are preserved. F/G withdrawn; increment EMIT-only. “Representational-difference-UNADJUDICATED” is not upgraded to “just thickness.” Comparator §7 says a d-projection residual does not close upstream families. Real write-keys match the delta (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, field/coefficient/Fourier carriers; no increment / term-origin / §3d rows).

### 5. N6 control — sound for this consumer

Kernel-level N6/N3/N4 belonged to c2 and is named as still-unclosed debt. S11c-d cannot discharge it by rewriting the already-built kernel (no pre-extraction operands, no term-origin rows). The downstream scattering-coordinate covariance + shape-sensitivity probes are labelled as such, not as kernel independence. `∇W_{\mathrm{bg}}→0` / `η→0` rejected as the vacuous uniform limit. Anchoring corruption excluded. `RHO4_CONSTANT` advection absence is emitted as structural absence, ⛔ no `A−A`. Uniform constructions are regressions (§5b), three separate constant backgrounds, not a gradient-coefficient control.

### 6. Two photon-kill channels — sound

Confinement is the computed transverse survival functional \(P_{T,\mathrm{surv}}=J_{T,\mathrm{out}}/J_{T,\mathrm{in}}\) from both reflected and transmitted T blocks — not an assertion, and not “energy stays in the slab.” Continuum `T→H` current and bound spectral overlap are distinct emitted objects. The bound pole is **not** a Bloch band (nonperiodic profile). Class-wide “a weak 1-D well always binds” is correctly **refused** for this class: an interface has unequal asymptotes, and the thickness pencil is multi-component, `ω`-dependent, nonlocal, and possibly non-Hermitian from the outgoing bulk. Existence is a profile-functional computation (possibly empty), with truncated-operator status called out because omitted \(O(η^2,σ_W^2)\) terms compete with weak binding \(E\simλ^2\). That is more accurate than applying the equal-asymptote 1-D theorem to a jump.

### 7. Answer/recipe discipline — sound aside from F1

The object is the complete two-ended open-channel S-matrix at retained background grade, plus the two kill channels and the flux FORM. No insertion sign, no `∝ k·a`, no `K_0=0`, no falsification number, no `O(1)` grating reductio. `η`, `σ_W`, `Q_n L_W`, `k_a L_W`, both densities, both anchorings, and independent `m₁` stay live. The §2 `GKψ` schematic is truncation/organization, not a target value. F1 is the M2 failure on this spec.

### 8. Completeness / house template — sound aside from F1

`N11a` rest-frame domain is carried (large \(k c_{s0}/|ω|\) necessary, not sufficient; no convective operator; no `v_0` alias). Chain topology, `IMPORT_KEYS` deferred to the directive, T7 comparator, blind WL, F9 reservations, no S11c-d card (N1: one family card after S11c-e) match the c2 house pattern. `N15`: no new constitutive constant. FORM at S11c-d is an N2-permitted boundary refinement and is not identified with S11c-e’s withheld strong-edge target.

---

## Nits (do not change computation or claims if F1 is fixed)

- **Title / authorship still say spec v4, Codex-authored, legs = fresh Claude + Grok.** The file is the orchestrator fold of that v4. That is process/G1 bookkeeping, not physics, but it should be updated so the next pairing is not chosen from a stale header.
- **§7 consume-set** omits `s11cc2MiddleMomentum*` / `s11cc2OutgoingNormalMomentum` that exist in the c2 delta. §7 already sends the exact root set to the build directive; not a physics defect if the closure pulls them in.

---

## Verdict

**NOT-SOUND.**

**Must-fix:** **F1** — the v5 Fourier fold still does not give engines an executable, unsupplied 3-D→1-D reduction of the imported carrier.

**Nits:** stale v4/Codex header; middle-momentum keys left to the directive.

No other outstanding item I derived from the sources changes what is computed or what this spec may claim.
