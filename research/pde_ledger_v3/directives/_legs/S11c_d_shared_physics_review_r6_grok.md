# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v6 (round 6, DOCUMENT branch)

VERDICT: NOT-SOUND (1 must-fix F1 = 3rd-layer Fourier: middle-leg hats + dtn_kernel handed + comparator-join not wired; 3 nits). Verbatim.

---

I'll start from the governing sources and form an independent view of what S11c-d must be, then review the spec against that — not against a prior finding list. First I'll load the review-legs procedure and the source-of-truth documents.Sources are loaded. Next I’ll form the physics view from the c2 record, house-format spec, N6 disposition, and consult, then open the v6 spec.I’ll inspect the c2 house-format spec and the real consume-set (closed kernel, Fourier symbols, what c2 actually strips) before opening S11c-d.I’ll inspect the real c2 kernel exports, Fourier conventions, and house-format comparator sections next — still without opening the S11c-d spec.The grep already exposed a few spec lines. I’ll still ground the Fourier/kernel facts in the real c2 exports and parent specs before judging v6.I have enough of an independent physics picture. Next I’ll mechanically inspect the consume-set, then read the spec itself.I have an independent picture of the consume-set and the parent requirements. The review below is from those sources against spec v6, not from a prior finding list.

---

# Independent physics review — `S11c_d_SHARED_PHYSICS.md` v6

## Independent view (formed from the sources, before the spec)

**Object (N5/N2).** S11c-d must name a profile class and emit the matching profile-conditioned linear mixing, not `ω(k)` for generic `W₀(x)`. For the slit-edge endgame the class is a **localized thickness interface**: smooth, short-range jets, asymptotically constant, with a finite thickness jump. That is equivalent to `∫ n̂·∇W_bg dy_n ≠ 0`, not a bump (`=0`). The matching object is a Born/scattering kernel in contrast `η`, with sharpness `σ_W` and kinematics live — not a gradient expansion, not Bloch, not WKB.

**Regime.** The imported operators are already truncated at first shape order in `η` and `σ_W`. Setting `η→O(1)` inside them is using a truncated operator outside its domain. Expanding in `σ_W` / sending `L_W→∞` kills the live form-factor argument `s=Q_n L_W` and leaves the localized-interface family (M3/N14). `η`, `σ_W=ηW̄₀/L_W`, and `kL_W` are independent (vary `η`, `L_W`, and `k` separately). The order-unity slit is a **new** finite-`ΔW` construction, not a reduction of Born (N7).

**Orders (N12/N10).** Linear response: converted amplitude `O(εη)` if the uniform baseline vanishes; converted flux `O(ε²η²)`; incident flux `O(ε²)`; fraction `C=J_conv/J_in=O(η²)` with `ε²` cancelling because the theory is linear. That `ε²` is `|amplitude|²`, not the excluded nonlinear-intensity program. Because c2’s F (uniform decoupling of the *closed* kernel) is **withdrawn**, `K_0=0` is not a supplied premise; the `O(η²)` label on physical `C` is conditional on the computed baseline.

**c2 consume-set (mechanical).** `s11cc2ClosedCouplingKernel` has **0** `DiracDelta`, **550** `s11cc2FourierW1ProfileHatTransfer` applications of which **48** take **middle-leg** arguments `(k_out−k_mid)` / `(k_mid−k_in)`, **526** `Integral`s containing `π^{-3}` (`=8/(2π)³`) and `L_W`, and thousands of *position-space* `w1_profile` jets. Momentum deltas live on c1 `dtn_kernel` `FLAT_DIAGONAL`; c2 strips them and applies the 3-D inverse in the integrals. Fourier hats have dimension `[L]³`. There is **no** WL self-energy engine: WL re-derives the closed kernel. Cross-engine content is N6 Reading B only (`(0)−(0)`, not operand AGREE); operand DEBT 40/76/18 is unadjudicated and material; `R_N6=18/288` and `R_cov` no-nonzero both stand; F/G withdrawn.

**1-D Fourier (derived).** Supplied 1-D convention `f̂_red(s)=∫dξ e^{-isξ}f(ξ)`, `f=(1/2π)∫ds e^{+isξ}f̂_red` is internally consistent. The 3-D→1-D *factor* is convention-dependent (`(2π)² L_W δ²(Q_∥)` vs `L_W/(2π) δ²(Q_∥)` differ by `(2π)³`) and must be computed from each engine’s realized 3-D convention, not supplied. c2’s own emit path is internally ambiguous (normalized inverse in `kernel_apply` vs `/(2π)³` on the *position* integral in `profile_bindings`), which is exactly why the map must not be typed.

---

## Must-fix findings

### F1 — §1a/§1c Fourier reduction is not executable on the real closed kernel (stripped deltas + middle-leg symbols + `dtn_kernel` as a false identity operand)

**Spec location.**

§1a (the consume-set description):

> its Fourier content is on the **momentum transfer** `Q = k_out − k_in` (`s11cc2FourierW1ProfileHatTransfer(−k_input + k_output, …)`), ⛔ not a single unspecified `k`.

and the reachable-import list in the same paragraph includes `dtn_kernel`.

§1c (the v6 control):

> Each engine must **COMPUTE and EMIT the 3-D→1-D reduction of its own closed-kernel Fourier symbols as an object with BOTH operands**: (i) the 3-D carrier as it actually appears in that engine’s own closed coupling kernel, and (ii) the reduced one-dimensional kernel obtained by applying that same engine’s own Fourier convention — the convention realized in its own construction of the **two-momentum identity / profile insertion**, ⛔ not a typed `[L_W/(2π)]` or `(2π)²L_W` map

§7 load-bearing comparator list names S-matrix / currents / conversion / survival / Riesz, and does **not** name the §1c both-operand reduced-kernel residual.

**Why it is wrong.**

Mechanical fact-lookup on the real export row `s11cc2ClosedCouplingKernel` (string counts on `scripts/S11c_c2_exports.py`, no CAS):

| object in the closed kernel | count |
|---|---|
| `DiracDelta` | **0** |
| `s11cc2FourierW1ProfileHatTransfer(...)` | 550 |
| of which arguments involve `s11cc2MiddleMomentum*` | **48** |
| `Integral(` | 526 |
| `Pow(pi` (`π^{-3}` = `1/(2π)³`) | 532 |

c2’s construction (verbatim) is the reason:

```373:374:research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
    deltas = diagonal.atoms(sp.DiracDelta)
    z0out = diagonal.xreplace({d: sp.S.One for d in deltas})
```

```463:466:research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
    # c1 uses DiracDelta(k-k') without a (2*pi)^3 coefficient; the transform
    # convention has an unnormalised forward transform and normalised inverse.
    p0 = integral(phase0 * diagonal * local_source / (2*sp.pi)**3, *limits0)
```

The deltas therefore live on c1 `dtn_kernel` `FLAT_DIAGONAL` (confirmed on `scripts/S11c_c1_exports.py:86`). The closed kernel’s Fourier symbols also include the retained mixed-grade **three-leg** insertion, which c2 itself calls required:

```398:400:research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
    # Second scattering is required at the mixed retained grade.  A three-leg
    # triangular representation evaluates the ordered operator product, with
    # its middle leg integrated below.  There is no single-momentum division.
```

Two concrete wrong-computation paths on **this** consume-set, not a hypothetical other input:

1. **Identity operand is the wrong object.** An engine that realises “the two-momentum identity” by binding the listed `dtn_kernel` (which still carries `DiracDelta(k−k′)` *without* `(2π)³`) infers a 3-D convention, then applies it *again* to a closed kernel whose Integrals already contain `π^{-3}=1/(2π)³`. Unnormalized 3-D forward vs normalized forward differ by `(2π)³` in `ŵ_3D ↔ δ²(Q_∥) ŵ_red`. That factor lands in the Born amplitude and squares in `C`. §1c’s “as it actually appears” cannot cancel this, because the sentence also *names* the identity construction, and §1a *hands* `dtn_kernel`.

2. **Middle-leg hats are not `Q=k_out−k_in`.** §1a’s exhaustive-sounding “Fourier content is on `Q=k_out−k_in`” is false of the payload: 48 hats are `ŵ(k_out−k_mid)` / `ŵ(k_mid−k_in)`. Dropping them (or leaving them unreduced) drops or mis-reduces the retained `ησ_W` second-scattering grade. S11c-a `:195–198` requires first order in *each* background bookkeeper, including the mixed grade.

Separately, §1c says the comparator joins the reduced kernels; §7’s load-bearing join list does not. That is how a v6 control disappears at comparator authoring.

This is not a supplied `(2π)` leak and not a pointer into a `.py` path. It is a consume-set misdescription that makes the leak-free control non-executable.

**Minimal fix.** In §1a and §1c, state as supplied consume-set facts (not as a typed map):

- the closed kernel **does not carry** 3-D momentum deltas; c2 already applied that identity, and the deltas on `dtn_kernel` are **not** the 3-D→1-D operand;
- Fourier symbols appear as the Function carriers **as they actually appear**, including **middle-leg** transfers of the retained three-leg insertion, not only `Q=k_out−k_in`;
- each engine reduces **every** such symbol plus the already-present 3-D Integrals, emitting both operands; ⛔ do not bind `dtn_kernel` to supply the convention;
- §7’s load-bearing comparator residuals include that both-operand reduced-kernel join.

---

## Scrutiny answers (no additional must-fix)

### 1. Profile class / regime (§1c–§1d, N5/N12/N14) — sound

Naming a **localized thickness interface** (`Δw₁ ≠ 0`, hence `W₋≠W₊`) is the correct non-global N5 object. Equivalence to the interface-vs-bump test:

\[
\int_{-\infty}^{\infty}\hat n\cdot\nabla W_{\mathrm{bg}}\,dy_n
= \sigma_W L_W\Delta w_1
= \eta\bar W_0\Delta w_1
= W_+-W_-.
\]

So `Δw₁≠0` iff `∫ W_0′ ≠ 0` along the edge normal. A bump has `Δw₁=0`. The spec states this via `Δw₁` and emits `(w_1′)̂_red(0)` and `Δw₁` as **two operands** (§5c), not as a typed Fourier theorem. Independent `m₁` (constant / bump / interface) is correct: modulus is not the class gate.

The three bookkeepers are genuinely independent: `η` by contrast, `σ_W=ηW̄₀/L_W` by `L_W` at fixed `η`, `kL_W` (and separately `s=Q_n L_W` vs local `k_a L_W`) by kinematics. The spec’s split of transfer `s` vs adiabatic `k_a L_W` is sharper than a single `kL_W`, and is right (WKB is a different N5 class).

Forbidding extra `σ_W→0` / forbidding `η→O(1)` inside the first-shape-order operators is right. An extra `σ_W` Taylor is a gradient expansion: it kills live `s`-dependence of the form factor and leaves the localized-interface family (M3). Evaluating the truncated vertex at `η∼O(1)` is using a first-order operator as if it were the finite-contrast theory — exactly the N7 category error. `L_W→0` at fixed `η` drives `σ_W→∞` outside the retained rectangle; the spec says so.

### 2. Order bookkeeping (§3c, N12) — sound, and more honest than typing `C=O(η²)`

Linear theory, one converted channel, homotopy `λ≡η` at fixed `L_W` (so `σ_W∝λ`):

\[
A_H(λ)=\varepsilon\bigl(a_0+λ a_1+λ^2 a_2\bigr),\qquad
J_{H,\mathrm{out}}(λ)=B_H(λ)[A_H,A_H],\qquad
J_{T,\mathrm{in}}(λ)=\varepsilon^2\bigl(j_0+λ j_1+\cdots\bigr).
\]

Expanding the current bilinear:

\[
\begin{aligned}
J_H^{(0)}&=\varepsilon^2 B_0[a_0,a_0],\\
J_H^{(1)}&=\varepsilon^2\bigl\{B_0[a_0,a_1]+B_0[a_1,a_0]+B_1[a_0,a_0]\bigr\},\\
J_H^{(2)}&=\varepsilon^2\bigl\{B_0[a_1,a_1]+\cdots\bigr\}.
\end{aligned}
\]

If the computed uniform baseline is `a_0=0` (S11b decoupling surviving the fold), then `J_H^{(0)}=J_H^{(1)}=0` and the leading conversion flux is `O(\varepsilon^2λ^2)=O(\varepsilon^2\eta^2)`. Incident flux is `O(\varepsilon^2)`. The fraction

\[
C=\frac{J_{H,\mathrm{out}}}{J_{T,\mathrm{in}}}=O(η^2)
\]

has `ε²` cancelling because both numerator and denominator are quadratic in a **linear** amplitude: `ε→2ε` multiplies both fluxes by 4 and leaves `C` invariant. That is not N10’s nonlinear-intensity program (which would make `C` itself intensity-dependent).

The spec emits both the `O(ε²η²)` flux slot and the `O(η²)` fraction as **N12 labels conditional on the computed baseline/interference disposition**, and emits the induced-field quadratic `B[ΔA,ΔA]` which carries those orders regardless. It does **not** type `K_0=0` from withdrawn F. That is the right contract: N12’s orders are the linear-mixing bookkeeping, not a licence to assume a withdrawn uniform-decoupling result.

### 3. Strong-edge bridge (§3d, N7) — sound

S11c-d emits weak Taylor coefficients (`∂_λ(ΔA/ε)|_0`, `C(0)`, `C'(0)`, `½C''(0)`, `lim C_{\mathrm{ind}}/λ²`) and names `C_{\mathrm{strong}}(1)` as what the lab bounds, withheld. The counterexample is elementary and sufficient. Lossless two-mode coupler, `A_T=\varepsilon\cos(ηG)`, `A_H=-i\varepsilon\sin(ηG)`:

\[
C=\sin^2(ηG)=η^2 G^2+O(η^4),\qquad C=0\text{ at }ηG=nπ\text{ for }G\neq 0.
\]

A nonzero Born coefficient `G` is compatible with **zero** strong conversion. No general positive lower bound exists; evaluating the Born coefficient at `η=1` is identifying weak data with `C_{\mathrm{strong}}(1)`. The order-unity edge is named as a **new** construction (piecewise-uniform / matched interior-exterior / effective interface with explicit undetermined parameters) whose weak limit must reproduce S11c-d, not a resummation of `F'(0)`. The spec does not claim `C_{\mathrm{edge}}=C_{\mathrm{interior}}·F(ω,ϑ)`. Nothing about the strong edge is over-claimed; the S11c-e handoff sentence is a named obligation, not a solved object.

### 4. Honest c2 import (§1b) — sound

Per-engine SOUND is the closed operator/kernel **VALUES**. Cross-engine is N6 covariance Reading B only; matched zeros are `(0)−(0)`, not operand AGREE. `R_N6=18/288` is preserved as a raw/schema-unmatched census; `R_cov` no-nonzero is preserved. Operand DEBT carrier 40 / source 76 / Φ 18 is UNADJUDICATED, leftover SHAPE uninspected, **material** to this consumer. F/G WITHDRAWN (and F is scoped to the increment, not to S11b’s uniform decoupling-as-oracle). Two S11c-b signs, six §3d, three N6 caveats carried. No silent upgrade to “closed” or “just thickness.” Uniform `A_0`/`K_0` is computed, not typed zero.

### 5. N6 control (§5a) — sound as a *downstream* regression; kernel N6 not over-claimed

The genuine N6 (independent shape/coordinate construction of the kernel + one-sided tilt/advection corruption; `∇W_0→0` forbidden; `RHO4_CONSTANT` advection structurally absent, ⛔ no `A−A`; corrupting an anchoring is not the test) **belonged to c2** and remains the §1b debt. S11c-d does not have term-origin rows or native pre-extraction Eulerian/material operands, so it **cannot** discharge kernel-level N6 by a chart rewrite — and it says so. What it runs is a scattering-coordinate covariance regression plus shape-**sensitivity** probes, explicitly not channel-isolation. `∇W_{\mathrm{bg}}→0` / `η→0` rejected; RHO4 absence emitted as structural absence; `Δρ` does not bridge anchorings. Anchoring corruption excluded. This is honest, not a vacuous uniform limit renamed.

### 6. Two photon-kill channels (§3b, N13) — sound

Confinement = survival of the transverse channel, emitted as `P_{T,\mathrm{surv}}=J_{T,\mathrm{out}}/J_{T,\mathrm{in}}` from both reflected and transmitted T-blocks, ⛔ not asserted. Continuum `T→H` and bound spectral overlap are distinct objects. The bound pole is correctly **not** a Bloch band (profile is nonperiodic). The spec is right to **reject** a class-wide “weak 1D well binds” existence claim: an interface has **unequal** asymptotes (a step, not a well that returns to 0), and the thickness pencil is multi-component, `ω`-dependent, nonlocal, and possibly non-Hermitian from the outgoing DtN. Existence is profile-functional and computed; poles are truncated-operator data until a remainder bound promotes them. No capture probability without a protocol. That is stricter than the naive 1D theorem and is the correct physics for *this* class.

### 7. Answer/recipe discipline (M2/M3) — sound on the object; F1 is the convention-executability hole

The object is named (complete two-ended S-matrix of the full asymptotic pencils, two kill channels, flux-normalized FORM). The `G=G_0+λG_1` display is labelled schematic. No expected value/sign of `K_0`, `C`, or a form-factor node is supplied. The `O(1)`/grating reductio is withheld. `∝k·a` is forbidden as a typed shape. Homotopy `λ≡η` at fixed `L_W` is a path through the two-parameter family, not an identification of grades. Varying `σ_W`, `s`, `k_a L_W` stay live. S11b decoupling is treated as M3 oracle, not a closed-kernel premise.

### 8. Completeness / house format — sound except F1’s comparator-join omission

Chain output (positional `load_model`, own-rows delta, `IMPORT_KEYS` at the directive, `BUILD_INPUT_DIGESTS` pins including `ledger_fold.py`), T7 comparator (join by name, paired residuals, three-valued, no native boolean, PRINTS/decides nothing), blind WL (imports nothing), N11 rest-frame (`|q v_{\mathrm{bulk\_normal\_0}}/ω|≪1`, large `kc_{s0}/|ω|` necessary not sufficient, ⛔ `v_0` alias), N14 fresh names, N15 no new constitutive constant — all present. §8 lists `f̂_red` **SUPPLIED** and the 3-D→1-D reduction factor **COMPUTED** with both operands. FORM-at-d vs FORM-at-e is an N2-permitted boundary refinement and does not import the withheld numeric bound. No N11 carry-in silently lost.

---

## Nits (do not change computation or claims)

- **`ξ` name collision.** S11c-a `:175` has vector `ξ≡y/L_W∈R³`; §1c reuses `ξ≡n̂·y/L_W` as a scalar. The derivative maps are rewritten, so this should not compute wrong, but a distinct `ξ_n` would prevent an engine from substituting the scalar into an inherited vector formula.
- **“Complete S-matrix”** means both ends / all open channels of a **first-grade** response (§2). A reader could hear non-perturbative `S`. One word (“Born-complete” / “first-grade complete”) would remove the ambiguity; the surrounding sentences already constrain it.
- **Authorship line** in the header (Codex v6, legs = fresh Claude + Grok) disagrees with the review packet’s “orchestrator-written” label. Process only.

---

## Verdict: **NOT-SOUND**

**Must-fix:** F1 — the §1c Fourier control is leak-free in intent (`f̂_red` supplied, no typed `(2π)` map, no script path) but **not executable on the real `s11cc2ClosedCouplingKernel`**: that kernel has no momentum deltas (they live on `dtn_kernel`, which §1a nonetheless lists), already contains `π^{-3}` Integrals, and carries middle-leg Fourier symbols that §1a does not name. That can change the reduced kernel’s `(2π)` content and the retained mixed-grade second-scattering piece, i.e. what an engine computes and what the spec may claim as a 1-D Born vertex.

Everything else in the eight scrutiny items is sound, including the interface class, the independent grades, the conditional N12 bookkeeping, the `sin²(ηG)` strong-edge honesty, the c2 debt/F/G/`R_N6`/`R_cov` framing, the downstream-only N6 regression, and the two computed photon-kill channels.
