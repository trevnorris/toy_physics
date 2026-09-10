# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v4 (round 4, Codex-authored, DOCUMENT branch)

VERDICT: NOT-SOUND (1 must-fix F1 + 4 nits). Verbatim (1-line trace preamble).

---

I'll form my own view from the governing sources first, then read the v4 spec against that view. Starting with the review skill and the required source documents.Sources are in. Next I’ll pin the `(ε,η,σ_W)` counting, the c2 objects this step actually consumes, and the confinement/conversion premises before opening the spec.I have the grade definitions from S11c-a. Next I’ll read the power-counting and coupling sections in full, then inspect the real c2 export objects.I have a source-grounded view of the objects and their status. Now I’ll read the v4 spec against that view.I’m checking the inherited Fourier convention and a few load-bearing identities against the real c2/c1 sources.The c2 engine comment on Fourier conventions may contradict the spec. I’ll read that code and the c1 transform definition next.I have enough to judge the eight scrutiny items. Checking one Fourier-convention detail in the c1 kernel identity so the reduced-map claim is sourced, not guessed.# Independent physics review — S11c-d SHARED PHYSICS v4

**Artifact:** `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`  
**Role:** Codex-authored physics spec (rule 15). Two blind engines will treat this as authority.  
**Method:** sources first (decision list, c2 record, c2 spec, N6 disposition, real c2 write-keys), then the spec. No working-tree edits.

---

## Source-grounded view (before the spec)

S11c-d must name a **profile class** and emit the matching object, not a generic `ω(k)` (**N5**). The inherited ansatz is

\[
W_{\mathrm{bg}}= \bar W_0[1+\eta w_1(\xi)],\qquad
\sigma_W=\eta\bar W_0/L_W,
\]

with \(\eta\) and \(\sigma_W\) independent via \(L_W\) (S11c-a §2a). A **localized thickness interface** is the right non-global class for the slit-edge endgame: asymptotically constant, \(\int W_{\mathrm{bg}}'=\eta\bar W_0\Delta w_1\neq 0\). A bump is \(\Delta w_1=0\). That is a class gate, not a spectrum.

The consumed c2 objects are the real write-keys `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` (increment is EMIT-only). They are **per-engine SOUND values**; cross-engine content is **N6 covariance (Reading B) only**. Matched vanishing is \((0)-(0)\), not operand agreement. The operand DEBT (carrier 40 / source 76 / Φ 18) is UNADJUDICATED and **material** to gradient-driven mixing. Preserve both \(R_{N6}=18/288\) (census) and \(R_{\mathrm{cov}}\) no-nonzero. F/G are withdrawn.

Power counting in a **linear** theory, with uniform conversion absent:

| object | order |
|---|---|
| coupling / converted amplitude | \(O(\varepsilon\eta)\) (and \(O(\varepsilon\sigma_W)\), mixed \(O(\varepsilon\eta\sigma_W)\)) |
| incident flux | \(O(\varepsilon^2)\) |
| absolute converted flux | \(O(\varepsilon^2\eta^2)\) |
| fraction \(C=J_{\mathrm{conv}}/J_{\mathrm{in}}\) | \(O(\eta^2)\) (\(\varepsilon^2\) cancels) |

An order-unity edge is **not** a reduction of this Born coefficient (**N7**): \(C=\sin^2(\eta G)\) can vanish at finite \(\eta G\) with \(G\neq 0\).

N6 at this step cannot re-derive kernel-level tilt vs advection (no term-origin rows). The genuine control is two coordinate constructions of the **new scattering object**, plus one-sided mutations; \(\nabla W\to 0\) is not a corruption. Confinement is **transverse survival** (**N13**); continuum conversion and bound capture are distinct; a 1D weak-well theorem does **not** apply to unequal asymptotes.

---

## Must-fix

### F1. The spec asserts a 3-D Fourier **carrier** identity that is not established for the imported c2 hats

**Spec (§1c):**

> For every c2 three-dimensional profile carrier, use its actual `(2π)⁻³` convention and the identities
>
> `(2π)⁻³∫d³y exp(−iQ·y) f(ξ) = [L_W/(2π)] δ²(Q_∥) f̂_red(s)`
>
> Thus `s11cc2FourierW1ProfileHatTransfer` maps by the first identity … The comparator uses the reduced per-unit-edge-area coefficient obtained after factoring exactly `[L_W/(2π)]δ²(Q_∥)`, and each engine also emits the reconstruction of the 3-D carrier.
>
> `A_3D(k_out,k_in) ≡ [L_W/(2π)]δ²(Q_∥) A_edge(...)`

The **integral identity is correct** for that definition of a normalized 3-D transform. With \(\xi=\hat n\cdot y/L_W\), \(s=Q_n L_W\), and \(\hat f_{\mathrm{red}}(s)=\int d\xi\, e^{-is\xi}f(\xi)\),

\[
\int d^3y\, e^{-iQ\cdot y}f(\xi)=(2\pi)^2\delta^2(Q_\parallel)\,L_W\,\hat f_{\mathrm{red}}(s),
\]

so \((2\pi)^{-3}\) times the left-hand side is indeed \([L_W/(2\pi)]\delta^2\hat f_{\mathrm{red}}\). That is not the issue.

The issue is the next sentence: it **identifies the named c2 carrier** `s11cc2FourierW1ProfileHatTransfer` with that normalized integral. That is a claim about an imported object, and the c1/c2 sources do not support it.

**c1 identity kernel** (flat piece) is a bare product of `DiracDelta`s, no \((2\pi)^3\):

```210:210:research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py
delta_k = sp.Mul(*(sp.DiracDelta(a - b) for a, b in zip(k_out, k_in)))
```

**c2 application** documents a *different* split — unnormalized **forward** transform, \((2\pi)^{-3}\) only on the **inverse/application**:

```463:466:research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
    # c1 uses DiracDelta(k-k') without a (2*pi)^3 coefficient; the transform
    # convention has an unnormalised forward transform and normalised inverse.
    p0 = integral(phase0 * diagonal * local_source / (2*sp.pi)**3, *limits0)
    p1 = integral(phase1 * off_diagonal * local_source / (2*sp.pi)**3, *limits1)
```

WL’s named hat is `Inactive[FourierTransform][w1Profile[...], spatial, k−k′]` (c1 `.wl:170–172`), i.e. an unevaluated transform whose \(2\pi\)/sign convention was never reduced to \((2\pi)^{-3}\int\). The c1 reconcile only identifies it as “FT of the thickness profile at the momentum transfer.”

If the hat is the unnormalized forward transform, the same geometry gives

\[
\hat w_{3\mathrm{D}}(Q)=(2\pi)^2\delta^2(Q_\parallel)\,L_W\,\hat f_{\mathrm{red}}(s)
=(2\pi)^3\cdot\frac{L_W}{2\pi}\delta^2\hat f_{\mathrm{red}}.
\]

Factoring the spec’s `[L_W/(2π)]δ²` then leaves an extra \((2\pi)^3\) in `A_edge`. That factor is **not** the R1 interior magnitude (**N7**); it is a kinematic-normalization error in the reduced vertex that both engines will share if they implement the asserted map.

The required “reconstruction of the 3-D carrier” does not catch this: `A_3D` is **defined** to be `[L_W/(2π)]δ² A_edge`, so reconstruction against that definition is \(A-A\).

**Why this changes computation / claims.** Every one-profile insertion in the closed kernel is specialized to the 1-D interface by replacing the opaque hat with \(\alpha\,\delta^2(Q_\parallel)\hat f_{\mathrm{red}}\). The spec freezes \(\alpha=L_W/(2\pi)\) as a supplied fact. Wrong \(\alpha\) rescales the Born vertex and the flux-normalized FORM. An error here makes both engines agree.

**Minimal fix.** Keep the definition of \(\hat f_{\mathrm{red}}\) (that is a supplied convention for S11c-d). **Do not** assert that `s11cc2FourierW1ProfileHatTransfer` equals \((2\pi)^{-3}\int\). Require each engine to **compute** the 3-D\(\to\)1-D coefficient from the imported kernel’s own identity (the flat `DiracDelta³` piece) plus the hat/jet insertion as it actually appears in `s11cc2ClosedCouplingKernel`, emit that factor as an object with both operands, and make reconstruction a round-trip against **that** derived map. The comparator should join the reduced kernels so obtained, not a pre-factored `[L_W/(2π)]` coefficient.

---

## Nits (do not change what is computed if F1 is fixed)

1. **`R_N6 = I_E − I_{M→E} = 18/288` (§1b).** The source census is “18 of 288 columns nonzero,” not the algebraic value of the residual (`S11c_c2_N6_reconcile_disposition.md` §1–2). Write it as a count. Engines do not recompute \(R_{N6}\); this is a claim-hygiene nit.

2. **“‘the sectors decouple at uniform background’ is exactly the withdrawn F” (§1b).** F was the c2 **increment** interpretation (closure-induced coupling; withdrawn instruments). S11b’s uniform decoupling is a different object and still stands as prior art. Operationally, computing \(K_0\)/`A_0` rather than typing \(0\) is the right M2/M3 move for **c2’s** closed kernel. Do not let the F label retract S11b.

3. **§3d menu of S11c-e constructions** (piecewise-uniform / matched interior-exterior / effective interface). Fine as a named handoff; slightly recipe-like for a later step. It does not feed this build.

4. **§5c `Emit (w₁′)̂_red(0)=Δw₁`.** This is a Fourier theorem, not a scattering answer. Emit both operands (form factor at zero, and the jump from limits) rather than asserting the equality as a single payload.

---

## The eight scrutiny items

### 1. Profile class / regime — sound (modulo F1’s map)

Naming a **localized thickness interface** (\(\Delta w_1\neq 0\), \(W_-\neq W_+\), short-range jets) is the correct non-global N5 object: a two-ended S-matrix, not \(\omega(k)\). Interface vs bump is stated correctly: \((\widehat{w_1'})_{\mathrm{red}}(0)=\Delta w_1\), and \(\int W_{\mathrm{bg}}'=\eta\bar W_0\Delta w_1\), so \(\int W'\neq 0\) iff \(\Delta w_1\neq 0\). Independent \(m_1\) (constant / bump / interface) as a subchannel discriminant, not a class gate, matches S11c-a’s independent \(w_1,m_1\).

The three quantities are genuinely independent if counted as \(\{\eta,L_W\text{ (or }\sigma_W\text{)},k\}\). The spec is **more precise** than calling \(kL_W\) a grade: \(\eta,\sigma_W\) are expansion bookkeepers; \(s=Q_n L_W\) and \(k_a L_W\) are kinematic arguments and are kept live. They are not the same axis (\(Q=k_{\mathrm{out}}-k_{\mathrm{in}}\) can be small at large \(kL_W\)).

**Weak-gradient would damage the object.** Extra \(\sigma_W\to 0\) at fixed \(\eta\) is \(L_W\to\infty\): the interface delocalizes and the coupling is the vacuous uniform limit (**N6**). The spec forbids that extra expansion, forbids Taylor-dropping \(s\), and separately forbids \(L_W\to 0\) at fixed \(\eta\) (\(\sigma_W\to\infty\), outside the retained rectangle). **Do not set \(\eta\to O(1)\) in the first-shape-order operators** is correct: those operators retain \(\eta^{\le 1},\sigma_W^{\le 1}\); finite contrast is a new construction (§3d).

A representative \((1+\tanh\xi)/2\) is correctly “instance, never the class.”

### 2. Order bookkeeping — sound; both labels belong; \(\varepsilon^2\) cancels

**Derivation.** Incident field \(u_T=\varepsilon\psi_{\mathrm{in}}\). On the physical path \(\lambda\equiv\eta\) at fixed \(L_W\), \(\sigma_W=(W_0/L_W)\lambda\), so the retained vertex is \(V=\lambda V_1+O(\lambda^2)\) with \(V_1\) containing both the \(\eta\) and \(\sigma_W\) first-jet pieces; the mixed grade \(\eta\sigma_W\) is \(O(\lambda^2)\) on this path and **is** retained (first order in each bookkeeper).

If the computed uniform baseline \(A_0=0\) (not assumed; withdrawn F / M3):

\[
A_H=\varepsilon\lambda\, G_0 V_1\psi_{\mathrm{in}}+O(\varepsilon\lambda^2)=O(\varepsilon\eta),
\]
\[
J_{\mathrm{in}}=\varepsilon^2 j_0+O(\varepsilon^2\lambda)=O(\varepsilon^2),
\]
\[
J_H=B[A_H,A_H]=\varepsilon^2\lambda^2\, B[G_0 V_1\psi,G_0 V_1\psi]+O(\varepsilon^2\lambda^3)=O(\varepsilon^2\eta^2),
\]
\[
C=\frac{J_H}{J_{\mathrm{in}}}=\lambda^2\frac{B[\cdots]}{j_0}+O(\lambda^3)=O(\eta^2).
\]

The \(\varepsilon^2\) in numerator and denominator cancel because both fluxes are quadratic in a linear amplitude. Emitting **both** \(O(\varepsilon^2\eta^2)\) (absolute flux) and \(O(\eta^2)\) (fraction) is required: they are different objects. An \(O(\varepsilon^2\lambda^2)\) quadratic of a **linear** vertex is not the N10 nonlinear program; the spec says so.

If \(A_0\neq 0\), \(J_H\) has \(O(\varepsilon^2)\) and \(O(\varepsilon^2\lambda)\) slots before \(O(\varepsilon^2\lambda^2)\). The spec keeps those slots live and attaches the N12 labels only after the computed baseline/interference disposition. That is the right response to withdrawn F and M3 (do not type \(K_0=0\)). It does not let a mis-ordered term read as nonlinear.

The homotopy \(\lambda\equiv\eta\) at fixed \(L_W\) does **not** freeze \(L_W\) in the operators (M3): \(L_W\) remains a live argument; the Taylor is the slit-edge path “contrast at fixed width.” Formal \((\varepsilon,\eta,\sigma_W)\) grades are still emitted separately.

On that path, omitted parent-theory \(\eta^2\) and \(\sigma_W^2\) are also \(O(\lambda^2)\). The spec’s warning that baseline interference can mix in uncomputed second-order amplitude, while \(B_0[a_1,a_1]\) does not need it, is correct Fermi-golden-rule bookkeeping.

### 3. Strong-edge bridge — sound and honest

S11c-d emits weak Taylor data (\(\partial_\lambda(\Delta A/\varepsilon)|_0\), \(C(0)\), \(C'(0)\), \(\tfrac12 C''(0)\), induced-field quadratic coefficient), not \(C_{\mathrm{strong}}(1)\). Lab bounds the strong **fraction**, not an amplitude evaluated at \(\eta=1\).

**Counterexample (baseline-free, linear in \(\varepsilon\)):** a lossless two-mode coupler with integrated coupling \(\eta G\),

\[
A_H=-i\varepsilon\sin(\eta G),\qquad
C=\sin^2(\eta G)=\eta^2 G^2+O(\eta^4).
\]

At \(\eta G=n\pi\), \(C=0\) with Born coefficient \(G\neq 0\). So a nonzero weak coefficient supplies **no** general lower bound on order-unity conversion. The spec labels this “not a model of the slab.” Correct.

The order-unity edge is out of scope and named as a **new** construction whose weak limit must reproduce S11c-d, not as a reduction of the first-jet kernel. Born can miss reconversion, diagonal reflection, and resonance shifts (frequency/angle, not only magnitude). Not over-claimed for this build.

### 4. Honest c2 supply — sound

§1b marks, in order:

- closed operator/kernel **values** as per-engine SOUND only;
- cross-engine content as the **N6 covariance thread only** (Reading B);
- matched zeros as \((0)-(0)\), not operand agreement;
- operand DEBT 40/76/18 UNADJUDICATED, leftover SHAPE uninspected, **material to this consumer**;
- explicit ban on upgrading that DEBT to “just thickness”;
- both \(R_{N6}\) nonzero raw / schema-unmatched **and** \(R_{\mathrm{cov}}\) no-nonzero;
- F and G withdrawn; uniform amplitude is computed `A_0`;
- S11c-b face-force and #90 signs, six §3d items, three Φ/`V`/leakage caveats, no term-origin/increment rows.

Real write-keys match the export (`s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, field/coefficient/Fourier carriers). Comparator: d-agreement does not close upstream operand families; projection can cancel differences. No quiet upgrade to “closed.”

### 5. N6 control — sound given the consume-set

The control is **not** the vacuous uniform limit. \(\nabla W_{\mathrm{bg}}\to 0\) and \(\eta\to 0\) are rejected as mutations. Uniform left/right/reference are a **separate** §5b smoke test with no prescribed vanishing (aligned with the owed c2 §5e clarification).

Kernel-level Eulerian↔material N3/N4 **cannot** be discharged here (no native pre-extraction operands, no term-origin rows). The spec says so and does not pretend a chart rewrite of the same kernel is an independent shape derivation. What it does instead is the right N6 for **this** object: two constructions of the scattering problem under \(x=X+u\) at fixed anchoring, residual in a common Eulerian channel basis, no \(\Phi\) on amplitudes (avoids the c2 §5c category error).

Mutations: reverse \(w_1'\) (tilt/N3) and omit/reverse \(u\cdot\nabla\rho_4/\rho_4\) (advection/N4), labeled as shape-sensitivity not channel isolation. `RHO4_CONSTANT`: computed structural absence, **no \(A-A\)**. `RHOBR_CONSTANT`: live \(\rho_{4\mathrm{D}}=\rho_{\mathrm{br}}/W_{\mathrm{bg}}\). Corrupting an **anchoring** is excluded. \(\Delta\rho\) does not bridge LAB_HELD\(\leftrightarrow\)MATERIAL_ADVECTED.

### 6. Two photon-kill channels — sound (the 1D-well restriction is a correctness fix)

Confinement is the computed survival functional \(P_{T,\mathrm{surv}}=J_{T,\mathrm{out}}/J_{T,\mathrm{in}}\) from **both** reflected and transmitted transverse blocks. Continuum \(T\to H\) and bound spectral overlap are separate tags. “Energy stays in the slab” is not used.

The bound pole is **not** a Bloch band (profile is nonperiodic). Existence is profile-functional: the class contains monotone steps **and** overshoot/wells. Rejecting a class-wide 1D weak-well theorem is correct: that theorem needs equal asymptotes and an attractive self-adjoint well; this operator is multi-component, \(\omega\)-dependent, nonlocal, and possibly non-Hermitian from outgoing \(Z\). Poles are first truncated-model data; promotion needs a remainder bound because weak binding \(E\sim\lambda^2\) competes with omitted \(O(\eta^2,\sigma_W^2)\). Spectral overlap only (no capture rate without a protocol) is honest; combining with continuum loss without a protocol would be the wrong object.

### 7. Answer/recipe discipline — sound except F1

The object is the complete two-ended S-matrix of the imported closed operator, not a derivation-path question. The GKψ schematic is organizational: it **forbids** identifying \(\langle\psi_0|K_1|\psi_0\rangle\) with the distorted-wave first-order result — necessary because a uniform-mode ME around \(W_0\) cannot treat a non-decaying end-value mismatch as a localized Born insertion. N5 requires the expansion to be fixed; this is that contract, not a typed answer.

No expected sign, no \(\propto k\cdot a\), no class-generic form-factor node, no \(A_0=0\), no unconditional \(C=O(\eta^2)\) as a builder target. Falsification numeric bound / \(O(1)\) grating reductio withheld. Varying quantities kept live (\(\eta,\sigma_W,L_W,s,k_a L_W,m_1\), both densities, both anchorings). Prior art (S11b decoupling, lab \(C(1)\)) is oracle, not premise.

F1 is the exception: a supplied numerical \(2\pi\) map is an expected-value leak about an imported carrier.

### 8. Completeness / house template / N11 — sound

Chain: positional `load_model(base, c1, c2)`, own-rows delta, `IMPORT_KEYS` at the directive, `BUILD_INPUT_DIGESTS` includes `ledger_fold.py`. Comparator: frozen T7, three-valued, join by channel labels, no native boolean. Blind WL: imports nothing. Supplied vs computed is explicit. N11a carried in full (rest-frame, \(|q\cdot v_{\mathrm{bulk\,normal\,0}}/\omega|\ll 1\), boundary-layer/subsonic, large \(k c_{s0}/|\omega|\) necessary not sufficient, never alias \(v_0\)). N14/N15 at the right layer. FORM at S11c-d is an N2-permitted boundary refinement; the weak FORM is required for the N12 order-count that N5 placed here.

Two-asymptote organization is necessary: \(\mathcal{L}_\pm\) contain zero-jet end values; expanding everything around \(W_0\) would make the \(\eta\) correction a **non-localized** perturbation. Channels from full end pencils, with computed \(K_0,K_\pm\) and classifiers, is the right object.

---

## Verdict

**NOT-SOUND.**

**Must-fix:** **F1** — do not supply `(2π)⁻³∫ = [L_W/(2π)]δ² f̂_red` as the identity of `s11cc2FourierW1ProfileHatTransfer` / `JetHat*`. Compute the 3-D\(\to\)1-D factor from the imported kernel.

**Nits:** \(R_{N6}\) census wording; F vs S11b decoupling label; S11c-e construction menu; emit both sides of \((\widehat{w_1'})_{\mathrm{red}}(0)=\Delta w_1\).

Nothing else in the eight items changes what an engine computes or what this spec may claim. The profile class, order bookkeeping (including emitting both \(O(\varepsilon^2\eta^2)\) and \(O(\eta^2)\)), strong-edge honesty, c2 DEBT/F-G/Reading B supply, N6-as-downstream-regression, and two-channel confinement construction are the right physics.
