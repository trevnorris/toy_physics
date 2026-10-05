# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v10 (round 10, DOCUMENT branch)

⭐ VERDICT: **SOUND** — no must-fix. The residual F1 seam is closed (universal reduced-representation rule +
per-section qualifiers leave no downstream construction path from unreduced 3-D content). 3 wording-only nits
(§5a Route E / §0.5 / §6 "consumed verbatim") that "change nothing computed". Verbatim.

---

I'll form my own view from the governing sources first, then read the S11c-d spec and report every way it could be wrong. Loading the review skill and the required source documents now.I've got the decision list. Next I'll read the c2 record, c2 spec, N6 disposition, and profile-class consult before opening the S11c-d spec.The consult prompt is context only. Next I’ll read the c2 SHARED PHYSICS house template and the objects S11c-d actually consumes.I'll keep forming the source-grounded view from c2's remaining house-format sections and S11c-a's power-counting/N6 material, still without opening the S11c-d spec.I still need S11c-a’s power counting and the real c2 export keys before I can judge what S11c-d must consume.I have a source-grounded view of the objects, orders, and controls. Next I’ll read the spec against that view, not against prior finding lists.I'll mechanically check the real c2 rows for Fourier hats, integral measures, and momentum deltas so I can judge whether the spec's reduction claims match the import.## Source-grounded view (formed before the spec)

S11c-d must name a **localized-interface** class and emit the **linear** transverse↔thickness mixing of the **closed** c2 operator/kernel: complete two-ended S-matrix, continuum conversion, a **computed** bound-pole object, a transverse **survival** functional, and `(ε,η,σ_W)` bookkeeping. It consumes the real c2 write-keys `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` (70-key mechanical census: **no** increment, term-origin, or §3d rows). Those two rows are **per-engine SOUND only**. Cross-engine content is the N6 covariance thread (Reading B): preserve **both** `R_N6 = 18/288` (raw, schema-unmatched) **and** `R_cov` no-nonzero; matched vanishing zeros are `(0)−(0)`, not operand agreement. The operand DEBT (carrier 40 / source 76 / Φ 18) is UNADJUDICATED and **material**. F and G are WITHDRAWN. N11a is inert. The order-unity slit edge and the `O(1)` reductio stay out of the builder.

Because the named class is one-dimensional in the edge-normal, `𝓛` and everything that feeds it must be the **§1c-reduced** 1-D operator/kernel. A weak well binds in 1-D only under equal asymptotes + an attractive self-adjoint well; an interface with `Δw₁ ≠ 0` is a **step**, so bound existence cannot be asserted class-wide.

---

## 1. Profile class / regime (§1c–§1d, N5/N12/N14)

**Localized interface is the right non-global object.** N5 forbids `ω(k)` for generic `W₀(x)` and requires a named class. A slit edge is a thickness jump, not a same-asymptote defect. The class gate

```text
Δw₁ ≡ w₁(+∞) − w₁(−∞) ≠ 0  ⇒  W₋ ≠ W₊
```

is the FTC identity `∫ w₁' = Δw₁`, so **interface vs bump** is stated correctly. Independent `m₁` (constant / bump / interface) is also correct: S11c-a already treats `w₁` and `m₁` as independent; `Δm₁` is a modulus-subchannel discriminant, not membership.

**The three (really four) live quantities are independent as bookkeepers.** S11c-a supplies `σ_W ≡ η W̄₀/L_W` with `η` and `L_W` independent, and `ε` independent of both. Incident `k` is a further kinematic axis. The spec keeps

- `η` — zero-jet contrast / Born strength,
- `σ_W` — first-jet sharpness,
- `s = Q_n L_W` — transfer argument of a form factor,
- `k_a L_W` — local wavelength / gap control,

and correctly refuses to infer WKB from large `|s|`. They are not four free grades of the **imported** operator (that operator is already truncated at first `η` and first `σ_W`); they are independent **parameters of the scattering object**. That is the right distinction.

**Forbidding a further `σ_W→0` / `L_W→∞` expansion is right, and would damage the object.** At fixed `η` and `Δw₁ ≠ 0`, `L_W→∞` replaces a finite-width scattering kernel by an infinitely slow two-media mismatch (the WKB class, which N5 assigns a different object). Expanding away `s` dependence replaces `f̂_red(s)` by its zero-transfer moment and kills the profile-conditioned kernel. Forbidding `L_W→0` at fixed `η` is the other wall: it drives `σ_W→∞` outside the first-jet truncation. The lab slit is that limit **and** `η→O(1)`; both are correctly out of scope.

**Do not put `η→O(1)` into first-shape-order operators: correct.** Those operators were built at first shape order. Using them at order-unity contrast is not a “reduction” of this kernel; it is a different construction (§3d).

No finding.

---

## 2. Order bookkeeping (§3c, N12) — derivation

Linear theory, one incident transverse amplitude of order `ε`, first-shape-order vertex of order `(η, σ_W)`:

\[
A_T \sim \varepsilon,\qquad
K = K_0 + \eta K_\eta + \sigma_W K_\sigma + \eta\sigma_W K_{\eta\sigma}+\cdots
\]

Born converted amplitude \(A_H \sim \varepsilon\,K\). Energy current is a **bilinear** on amplitudes, so homogeneous of degree 2:

- \(J_\mathrm{in}\sim\varepsilon^2\),
- if \(K_0=0\) (uniform sectors decouple), \(A_H=O(\varepsilon\eta)\) on the physical path \(\sigma_W\propto\eta\), hence **absolute** converted flux \(J_\mathrm{conv}=O(\varepsilon^2\eta^2)\),
- **fraction** \(C=J_\mathrm{conv}/J_\mathrm{in}=O(\eta^2)\). The \(\varepsilon^2\) **cancels because both fluxes are quadratic in the same linear fields**.

That cancellation is the signature of a **linear** vertex. A nonlinear-intensity program would make \(C=C(\varepsilon)\). Emitting **both** the \(O(\varepsilon^2\eta^2)\) flux label and the \(O(\eta^2)\) fraction label is therefore required, and they are not interchangeable.

If the computed uniform baseline \(K_0\) (or \(a_0\)) does **not** vanish — which this spec must not assume, because c2 F is WITHDRAWN and S11b’s uniform zero is an M3 oracle, not a typed \(K_0=0\) — then

\[
A_H\sim\varepsilon(a_0+\lambda a_1+\cdots),\qquad
J_H\sim\varepsilon^2\big(|a_0|^2 + 2\lambda\,\mathrm{Re}(a_0^*a_1)+\lambda^2|a_1|^2+\cdots\big),
\]

so the physical \(C\) has an \(O(1)\) slot and an \(O(\lambda)\) interference slot **before** \(O(\lambda^2)\). The spec’s split — always emit component grades and the induced-field quadratic \(B[\Delta A,\Delta A]\) with \(O(\varepsilon^2\lambda^2)\) / \(O(\lambda^2)\) labels; attach N12’s “\(O(\varepsilon\eta)\) coupling, \(O(\varepsilon^2\eta^2)\) leakage, \(O(\eta^2)\) fraction” to **physical** \(C_{T\to H}\) only after the computed baseline/interference disposition — is the honest version of N12, not a retreat from it.

It also blocks the N10 misread: “An `O(ε²λ²)` quadratic observable of a linear vertex is not the excluded nonlinear-light program.”

On the named path \(\lambda\equiv\eta\) at fixed \(L_W\), mixed grade \(\eta\sigma_W\) is \(O(\lambda^2)\) **and retained**; pure \(\eta^2\) / double insertion are the same order as omitted operator terms. The spec flags that truncation rather than promoting the displayed \(\lambda^2\) coefficient. Correct.

No finding.

---

## 3. Strong-edge bridge (§3d, N7) — derivation

S11c-d emits only weak Taylor data \(\partial_\lambda(\Delta A_H/\varepsilon)|_0\), \(C(0)\), \(C'(0)\), \(\tfrac12 C''(0)\), and \(\lim C_{H,\mathrm{induced}}/\lambda^2\), and states that the lab bounds \(C_\mathrm{strong}(1)\). The order-unity edge is named as a **new** S11c-e construction whose weak limit must reproduce S11c-d, not as a reduction of this kernel.

Counterexample (lossless two-mode coupler, **not** the slab): \(\mathrm{i} A_T'=\eta G_0 A_H\), \(\mathrm{i} A_H'=\eta G_0 A_T\), \(A_T(0)=\varepsilon\), \(A_H(0)=0\). Then \(A_H=-i\varepsilon\sin(\eta G)\) with \(G=G_0\times\mathrm{length}\), and

\[
C=\sin^2(\eta G)=\eta^2 G^2+O(\eta^4),\qquad C=0\text{ at }\eta G=n\pi\ (n\ge 1).
\]

Born coefficient \(G\neq 0\), yet exact conversion vanishes at finite coupling. So a nonzero weak coefficient supplies **no** general lower bound on strong-edge conversion, and substituting the Born value at \(\eta=1\) is exactly the identification N7 forbids.

Nothing about the strong edge is over-claimed. The “if it cannot be established, a conditional constraint on edge-response parameters, not a shape-independent exclusion” sentence is the right under-claim.

No finding.

---

## 4. Honest c2 import (§1b)

Checked against the c2 record, the N6 disposition, and a mechanical key census of `scripts/S11c_c2_exports.py` (70 keys; `s11cc2ClosedSlabOperator` / `s11cc2ClosedCouplingKernel` present; increment/term-origin/§3d rows absent):

| Required mark | Spec |
|---|---|
| Closed operator/kernel **VALUES** per-engine SOUND only | yes |
| Cross-engine = N6 covariance thread only, Reading B | yes |
| Matched vanishing zeros are `(0)−(0)`, not operand AGREE | yes |
| Both `R_N6 = 18/288` (raw, schema-unmatched) **and** `R_cov` no-nonzero | yes |
| DEBT carrier 40 / source 76 / Φ 18 UNADJUDICATED, leftover SHAPE uninspected, **material** | yes |
| ⛔ not “just thickness” | yes, verbatim |
| F and G WITHDRAWN; increment EMIT-only, not imported | yes |
| Two S11c-b face-force / #90 signs; six §3d items; three N6 premise caveats | yes |
| S11b uniform decoupling stands as M3 oracle; do not type `K_0=0` | yes |
| d projections do not close upstream families | §7, yes |

No silent upgrade. No finding.

---

## 5. N6 control (§5a)

Kernel-level N6 (independent Eulerian vs material **construction** of the off-diagonal kernel, then one-sided tilt/advection corruption) was c2’s control and remains the carried DEBT. S11c-d does not import native pre-extraction operands or term origins, so a chart rewrite of the **same** closed kernel would not be an independent route. The spec says so and does not pretend to discharge it.

What it does instead is a **downstream scattering-coordinate covariance regression** plus shape-**sensitivity** probes (reverse `w₁'` / omit an identifiable `u·∇ρ_4D` factor), explicitly **not** labelled as isolated N3/N4 channels. That is the honest control given the consume-set.

- Vacuous uniform limit is **not** the control; it is a separate §5b regression that “does not validate a gradient coefficient, sign, or parity.”
- `∇W_bg→0` and `η→0` are rejected as mutations.
- Anchoring corruption is excluded; `Δρ` does not bridge `LAB_HELD↔MATERIAL_ADVECTED`.
- `RHO4_CONSTANT`: emit computed structural absence of the advection factor, ⛔ no `A−A`. `RHOBR_CONSTANT` keeps live `∇ρ_4D`.

No finding.

---

## 6. Two photon-kill channels (§3b, N13)

Confinement is **survival of the transverse channel**, emitted as the computed functional

\[
P_{T,\mathrm{surv}}[a_T]=J_{T,\mathrm{out}}/J_{T,\mathrm{in}}
\]

with both reflected and transmitted transverse blocks, ⛔ not asserted. Continuum \(T\to H\) conversion and bound **spectral overlap** are distinct payloads. No “energy stays in the slab.” Not a Bloch band (profile is nonperiodic, stated independently of the pole-set disposition).

**The spec is right to refuse the slogan “a weak 1-D well binds.”** That theorem needs (i) 1-D, (ii) **equal** asymptotes, (iii) an attractive well, (iv) a self-adjoint operator. This class has \(W_-\neq W_+\) (a **step**, not a well); a well appears only for profiles with overshoot; the thickness pencil is multi-component, \(\omega\)-dependent, nonlocal, and possibly non-Hermitian from the outgoing bulk response. Existence is therefore profile-functional and may be empty. Capture **rate** is withheld because no preparation/switching protocol is supplied; spectral overlap is the object that can actually be computed. That is N13, not a dodge.

Poles are first truncated-model spectrum of the first-shape-order operator; promotion needs a remainder bound competing with \(E\sim\lambda^2\). Correct, and it does not leak a pole location.

No finding.

---

## 7. Answer / recipe discipline (M2/M3)

The named objects are the two-asymptote S-matrix, currents, conversion fraction, survival functional, pole set + Riesz data. The homotopy schematic is labelled as logic, not a target. No expected sign/value. Falsification numeric bound / grating reductio withheld. `(1+\tanh\xi)/2` is an instance, never the class. `f̂_red` is a supplied 1-D convention, not a 3-D `(2π)` map. Varying quantities (`η`, `σ_W`, `s`, `k_a L_W`, live end values, live density maps) stay live. The physical path `λ≡η` at fixed `L_W` is a **named homotopy**, not a freeze of the independent grades.

No finding.

---

## 8. Completeness / F1 seam / deferral

House template is intact: honesty section, supplied-vs-computed, bind-closure own-rows delta, T7 comparator, blind Wolfram, N9 denylist cut, N11a rest-frame (`|q\,v_\mathrm{bulk\_normal\_0}/\omega|\ll 1`, large \(kc_{s0}/|\omega|\) necessary not sufficient, never `v_0`), N14 fresh names, N15 moments as scattering data not new constitutive constants.

**F1 seam (the v10 question).** The governing rule is

> Wherever §3, §5, or §6 refers to "the imported closed operator", "the closed operator", or "consumed verbatim" in constructing `𝓛`, its **currents, modes, normalization, resolvent, poles, S-matrix, conversion fraction, survival functional, or the §5a routes**, it denotes the engine's **§1c-reduced** closed operator (and reduced coupling kernel), ⛔ **never** the unreduced 3-D rows.

Local qualifiers match that rule at the load-bearing sites: §0.2 (after §1c reduction of **both** rows), §2 (`𝓛` from reduced rows; unreduced rows are reduction operands **only**), §3a (current and `J` on reduced `𝓛_e`, twice), §3b (pole solve of reduced `𝓛(ω)`), §5a Route M (reduced), §5b (recompute the **reduced** operator), §6 (modes/currents from reduced pencils), §7 (join reduced rows **and** the object built from them, ⛔ “not a side reduction while `𝓛` uses unreduced content”).

Mechanical check of the real rows supports the reduction premise: **no** `DiracDelta`/`KroneckerDelta` in either closed payload (c2 did strip 3-D momentum deltas); both rows **do** carry `Integral`, `Pow(pi, Integer(-3))`, transfer hats, jet hats, `s11cc2MiddleMomentum*`, `s11cc1_k_{in,out}put_*`, and `s11cc2Y{1,2,3}`. Qualitative census in §1a is consistent with the files; the exact per-symbol census remains directive-side.

**Deferral remains legitimate.** The spec names the object (both-operand 3-D→1-D reduction of every convention-bearing element of **both** rows), supplies only `f̂_red`, forbids typing `[L_W/(2π)]δ²(Q_∥)`, forbids an `A−A` reconstruction, and puts the leak-safe control in the comparator join of **reduced** operators/kernels. That is the same pattern as deferring `IMPORT_KEYS` to the build directive. The directive must freeze the **element census and reduction procedure**, not the `(2π)` arithmetic **result** — otherwise both engines would agree on a leaked factor. That is a handoff constraint on the next artifact, not a defect in this spec.

Nothing the v10 fold touched newly breaks the physics: currents, modes, poles, S-matrix, conversion, and survival are all pinned to reduced `𝓛`. The leftover wording nits below do not reopen a construction path.

FORM at d vs e: N2 allows spec-stage boundary refinement; d’s weak-contrast FORM is not e’s withheld order-unity target. N13 confinement **interpretation** remains e; d emits the computed survival object. Acceptable.

---

## Findings

**Must-fix:** none. Nothing outstanding changes what an engine computes or what the spec may claim.

**Nits** (wording only; a careful engine following §1c+§2+the universal rule already does the right thing):

1. **§5a Route E is not locally qualified, unlike Route M.**
   Spec: “Route E constructs the §2–§3 scattering problem directly in Eulerian `x` coordinates.” Route M: “rewrites in `X` the §2–§3 scattering problem built on the **reduced** closed operator.”
   Why it is only a nit: §2 already defines that scattering problem as `𝓛(y_n;ω,k_∥)` after reduction, and the universal rule names “the §5a routes.” An engine that ignored those and built Route E from unreduced 3-D `x` would make the covariance residual incomparable.
   **Minimal fix:** say Route E constructs the **reduced** `𝓛(y_n;ω,k_∥)` in Eulerian \(y_n\).

2. **§0.5 sits outside the universal rule’s §§3–6 scope.**
   Spec: “The scattering-coordinate covariance regression on the **imported closed kernel**.”
   **Minimal fix:** “on the §1c-reduced closed kernel.”

3. **§6 category mix-up.**
   Spec: “the **reduced** closed operator and its **reduced** off-diagonal extraction **are the c2 exports consumed verbatim**.”
   The c2 exports are the unreduced 3-D rows; the reduced rows are **computed** here. Intent is clear from §1c+§2; the sentence can be read as “skip reduction.”
   **Minimal fix:** “the construction operands of `𝓛` are the **computed §1c reductions** of the c2 exports, consumed without further constitutive modification.”

---

## Verdict: **SOUND**

The v10 universal reduced-representation rule, together with the §0.2 / §3a / §3b / §5a / §5b / §6 / §7 qualifiers, leaves **no construction seam** that still feeds `𝓛`, flux normalization, the S-matrix, conversion, survival, or the pole set from unreduced 3-D content. The Fourier deferral is still the right split (object + leak-safe control here; census in the directive). Profile class, independent grades, order bookkeeping (both \(O(\varepsilon^2\eta^2)\) and \(O(\eta^2)\), N10-safe), strong-edge honesty, c2 import status, N6 demotion, and the two photon-kill channels all match the sources. The three nits do not change what is computed or may be claimed.
