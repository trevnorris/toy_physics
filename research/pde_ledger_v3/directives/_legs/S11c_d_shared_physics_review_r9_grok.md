# Grok (grok-4.6 high) review — S11c-d SHARED PHYSICS spec v9 (round 9, DOCUMENT branch)

VERDICT: NOT-SOUND (1 must-fix F1: 𝓛 bound to reduced rows, but J/S/C/survival/pole/§5a still say
"imported closed operator" -> mixed-convention path; fix = one operand rule everywhere). 3 nits. Deferral endorsed.

---

I'll form my own view from the governing sources first, then read the S11c-d spec against that view. Starting with the decisions, c2 record, house template, and N6 disposition.Next I’ll read the c2 record, house template, N6 disposition, and profile-class consult so I have an independent picture before opening the spec.I have the c2 record and N6 disposition. Next I’ll read the c2 house-format spec, then the S11c-d spec.I have the spec. Next I’ll pull the c2 house-format remainder, the S11c-a profile definition, and the actual c2 Fourier carriers so I can check the fold against the real rows.I’ll check S11c-a’s profile definition, the real c2 Fourier export keys, and every place the spec builds `𝓛` or currents, looking for mixed 3-D/1-D conventions.## Independent view (from the sources, before the spec)

S11c-d must name a **localized thickness interface** (`Δw₁ ≠ 0` ⇔ `∫ W_bg′ ≠ 0`), emit the **linear mixing** of the c2 closed operator/kernel as a two-asymptote scattering object, keep `(ε, η, σ_W)` and the kinematic axes live, treat c2 values as **per-engine SOUND** and the cross-engine content as **N6 Reading B plus an unadjudicated operand DEBT**, emit **two distinct photon-kill channels**, and give only weak Taylor data plus a named strong-edge obligation. The v9 fold’s direction is right: `𝓛` has to be built from **reduced** rows, and every convention-bearing 3-D piece of both rows has to be reduced. The census/`(2π)` arithmetic may stay in the build directive.

---

## Must-fix

### F1 — Reduced-operand binding stops at `𝓛`; currents / flux-normalized `S` / covariance can still mix 3-D and 1-D conventions

**Spec locations.**

§2 does bind `𝓛`:

> After each engine performs the §1c 3-D→1-D reduction of **both** imported rows — the closed operator and the closed coupling kernel — the object `𝓛(y_n;ω,k_∥)` is constructed from that engine's **full reduced closed slab operator**, with its **reduced closed coupling kernel** … The reduced rows are the construction operands of `𝓛`, ⛔ not the imported unreduced rows and ⛔ not merely comparator objects.

§7 repeats that for `𝓛` and the S-matrix:

> the `𝓛` and complete flux-normalized S-matrix/mixing amplitudes that each engine constructs from those reduced rows … ⛔ not a side reduction while `𝓛` uses unreduced content.

The same document then builds the **emitted mixing object** from the **unreduced** import.

§3a (twice):

> Derive `J` from the S11b quadratic energy current, evaluated on the **imported closed operator** …
>
> `J_{ab}^{(e)} ≡` the polarized S11b normal energy-current bilinear `𝓙_n[l_a,r_b]` **evaluated with the imported closed operator**.

Modes and the pencil derivative are taken from reduced `𝓛`:

> compute right modes `r_a` and adjoint/left modes `l_a` of `𝓛_e^full(ω,k_n,k_∥)` …
> `N_{ab}^{(e)} ≡ ⟨l_a, (∂_ω𝓛_e^full) r_b⟩` …
> Also emit the identity relating the derived current, `∂_{k_n}𝓛_e^full`, and `∂_ω𝓛_e^full` for the engine's Fourier convention.

§2’s pole exemption:

> The pole solve of §3b uses the resolvent of the retained first-shape-order **imported operator**

while §3b’s pencil is `𝓛(ω)` from §2.

§5a:

> Route E constructs the §2–§3 scattering problem directly in Eulerian `x` coordinates. Route M rewrites the **same imported closed operator** in `X`.

§6:

> The closed operator and its off-diagonal extraction are the c2 exports **consumed verbatim**. Compute `K₀,K₋,K₊`; define modes and currents from the full asymptotic block pencils.

§0 item 2 still sources the two-asymptote basis from the unreduced keys `s11cc2ClosedSlabOperator` / `s11cc2ClosedCouplingKernel`, with no “after reduction.”

§1c already forbids using a flux `δ²` strip as a substitute for reducing the carrier:

> This flux step is distinct from, and ⛔ not a substitute for, the computed carrier reduction.

**Why this is wrong.**

c2’s closed rows are 3-D Fourier objects. The PY fold applies the inverse with `1/(2π)³` and integrates `d³y`, `d³k_out`, `d³k_in`, and `d³k_middle` (`S11c_c2_selfenergy_fold_sympy_audit.py`, `kernel_apply`), after stripping c1’s unnormalized `DiracDelta(k−k′)`. Hats have dimension `[L]³`. The 1-D convention supplied in §1c is a different object:

```text
f̂_red(s) ≡ ∫ dξ e^{−isξ} f(ξ) ,   f(ξ) = (1/2π) ∫ ds e^{+isξ} f̂_red(s) .
```

If modes/`N`/`∂_ω𝓛` live in the reduced 1-D convention and `J` is the S11b current of the unreduced 3-D operator, then:

1. `l_a, r_a` and `𝓙_n[l_a,r_b]` do not act in the same Fourier representation.
2. The required identity among `J`, `∂_{k_n}𝓛`, and `∂_ω𝓛` cannot hold across a 3-D vs 1-D `(2π)`/`δ²(Q_∥)` mismatch.
3. Flux-normalized `S`, `C_{T→H}`, and `P_{T,surv}` all pass through `J`. Those are the N5/N12/N13 objects, not `𝓛` as a block matrix.
4. Evaluating 3-D `J` and then stripping `δ²(Q_∥)` at the flux step is exactly the substitute §1c forbids.
5. §5a covariance of “§2–§3 scattering” on the unreduced operator is a different residual from covariance of reduced `𝓛`.
6. §6 “consumed verbatim” plus §0’s unreduced keys give an engine a licensed path to skip reduction when assembling currents and the S-matrix.

So: **`𝓛` is bound to the reduced rows. The objects S11c-d actually emits are not.** Two engines can agree on mixed-convention fluxes, or disagree because the spec offers both readings. That is a spec error, not a builder choice.

The fold did not newly break the profile class, N6 honesty, or the `f̂_red` deferral. It left the current/S/covariance/method sentences on the pre-fold “imported closed operator” wording.

**Minimal fix.** One operand rule, applied everywhere that constructs mixing data:

- 3-D closed rows = operands of the §1c reduction only.
- `𝓛`, modes, `N`, `J`, flux-normalized `S`, `C`, survival, the pole pencil, and the §5a routes are all built from the **same reduced operator/kernel**.
- Replace “evaluated on / consumed verbatim / imported closed operator” in §3a, §2’s pole sentence, §5a, and §6 with “the reduced closed operator that is the construction operand of `𝓛`.”
- Keep the S11b current *formula* and the closed/nonlocal bulk content; change only the Fourier representative they are evaluated on.
- In §0 item 2, say the two-asymptote basis is obtained **after** the §1c reduction of those two rows.

The §1c element principle on the **rows** can stay. It already covers hats at every argument and every `d³y` / `d³k` measure, with the closure “no convention-bearing element left in a 3-D convention.” That principle is not the hole; the hole is objects derived after the rows.

---

## Nits (do not change the verdict by themselves)

**N1.** The §1c including-list names hats and measures, not the plane-wave phases `e^{ik·y}` or overall `(2π)` weights that sit in those integrands. The closure sentence covers them, and c2 puts `/(2π)³` and the phases inside the integrals, so a faithful reduction of the measure reduces them. Optional: add “including the phases and `(2π)` weights in those integrands.”

**N2.** §1a lists PY Fourier write-keys while deferring the per-row census. The spec already calls this house chain-wiring. It is not a typed `(2π)` map. Leave it unless it starts being read as the census.

**N3.** §5b “recompute the full operator, `K`” does not say “reduced.” A translation-invariant uniform background has no profile hats, but it still has a field Fourier convention. Point it at the same reduced representative as `K₀`.

---

## Scrutiny (no further must-fix)

**1. Profile class / regime.** Naming a localized interface — asymptotically constant, short-range retained jets, `Δw₁ ≠ 0` hence `W₋ ≠ W₊` — is the right non-global N5 object. Interface vs bump is `Δw₁ ≠ 0` vs `= 0`, which is `∫ ∂_{y_n} W_bg = W̄₀ η Δw₁ ≠ 0` vs `= 0`. Correct. `m₁` is independent and not a class gate.

The three axes are independent in the S11c-a sense (`S11c_a_SHARED_PHYSICS.md:189–191`): `η` and `σ_W = η W̄₀/L_W` are independent because `L_W` is an independent length; `s = Q_n L_W` and `k_a L_W` are kinematics, not grades. An extra `σ_W → 0` / `L_W → ∞` expansion would kill the first-jet vertex and/or send `|s| → ∞` at fixed `Q_n` (Riemann–Lebesgue on `(f′)̂_red`, and the WKB class, which N5 forbids). `L_W → 0` at fixed `η` drives `σ_W → ∞`, outside the imported first-shape-order truncation. Forbidding both limits is right. Setting `η → O(1)` inside that truncation is the strong edge, a new construction, not a reduction.

**2. Order bookkeeping.** Linear response, incident `ψ_T ∼ ε`, vertex from the retained operator `K = K_0 + λ K_1 + ⋯` along the physical path `λ ≡ η` at fixed `L_W` (so `σ_W ∝ λ` tracks, formal grades still separate):

\[
A_H = \varepsilon\bigl(a_0 + \lambda a_1 + \lambda^2 a_2 + \cdots\bigr),
\]
\[
J_{H,\mathrm{out}} = B[A_H,A_H]
= \varepsilon^2 B_0[a_0,a_0]
+ \varepsilon^2\lambda\bigl(B_0[a_0,a_1]+B_0[a_1,a_0]+B_1[a_0,a_0]\bigr)
+ \varepsilon^2\lambda^2\bigl(B_0[a_1,a_1]+B_0[a_0,a_2]+\cdots\bigr).
\]

Incident flux `J_{\mathrm{in}} = \varepsilon^2(j_0 + \cdots)`. If the computed baseline `a_0 = 0` (S11b uniform decoupling as M3 oracle, not a typed `K_0 = 0`):

\[
J_{H,\mathrm{out}} = O(\varepsilon^2\lambda^2) = O(\varepsilon^2\eta^2),\qquad
C = J_{H,\mathrm{out}}/J_{\mathrm{in}} = O(\lambda^2) = O(\eta^2).
\]

The `ε²` cancels because both fluxes are quadratic in a **linear** amplitude. That is not the N10 nonlinear program (vertex depending on intensity, `A_H ∼ ε³`, `C` intensity-dependent). Emitting both labels is correct if they sit on absolute flux vs fraction. The spec attaches N12’s “`O(εη)` coupling, `O(ε²η²)` leakage, `O(η²)` fraction” only after the computed baseline/interference disposition, and it emits `A_0`/`K_0` rather than typing `K_0 = 0` from withdrawn F. That is the right refinement of N12 given M3 and withdrawn F.

**3. Strong-edge bridge.** Honest: S11c-d emits weak Taylor data (`∂_λ(ΔA_H/ε)|_0`, `C(0)`, `C'(0)`, `½C''(0)`, induced-field quadratic coefficient) and names `C_{\mathrm{strong}}(1)` as the lab object. Order-unity edge is out of scope and a **new** construction.

Counterexample, lossless two-mode coupler with integrated coupling `ηG`:

\[
\frac{dA_H}{d\zeta} = -i\,(\eta g)\,A_T,\quad
A_H = -i\varepsilon\sin(\eta G),\quad
C = \sin^2(\eta G) = \eta^2 G^2 + O(\eta^4).
\]

Born coefficient `G^2 ≠ 0`, but `C = 0` at finite `ηG = nπ`. A nonzero weak coefficient gives **no** lower bound on strong-edge conversion. The spec does not overclaim a shape-independent exclusion.

**4. Honest c2 import.** §1b matches the c2 record and the Path-B disposition: per-engine SOUND values; cross-engine = N6 covariance thread only; matched covariance zeros are `(0)−(0)`, not operand AGREE; DEBT carrier 40 / source 76 / Φ 18 UNADJUDICATED and material; leftover SHAPE uninspected; ⛔ not “just thickness”; both `R_{N6} = 18/288` (raw census, schema-unmatched) and `R_{\mathrm{cov}}` no-nonzero preserved; F/G withdrawn (c2 increment only; S11b uniform decoupling stands as oracle, `K_0` computed). No quiet upgrade found.

**5. N6 control.** Not the vacuous uniform limit. Genuine kernel-level N6/N3/N4 stay c2’s unclosed debt. §5a is a downstream scattering-coordinate regression and says so. `∇W_{\mathrm{bg}} → 0` / `η → 0` rejected; anchoring corruption excluded; `RHO4_CONSTANT` structural absence with no `A−A`; `RHOBR_CONSTANT` keeps live `∇ρ_{4D}`. Tilt and advection probes are shape-sensitivity, not isolation — correct, because the import has no term-origin rows.

**6. Two photon-kill channels.** Confinement = transverse-channel survival, emitted as `P_{T,\mathrm{surv}}`. Continuum `T→H` and bound spectral overlap are distinct. The spec correctly **does not** invoke the 1-D weak-well theorem: unequal asymptotes, multi-component / nonlocal / `ω`-dependent, possibly non-Hermitian from the outgoing bulk. Existence is computed and may be empty. Not a Bloch band. No capture probability without a protocol. Right.

**7. Answer/recipe discipline.** The object is the complete two-ended S-matrix at first background grade, not a derivation-path question. Power-counting contract is supplied; numeric `O(1)` / grating reductio withheld. `f̂_red` is a supplied 1-D convention, not a typed map for imported hats. No M3 freeze of `η`, `σ_W`, or `s`; the `λ ≡ η` homotopy is a named physical path and does not identify formal grades. Prior art (`K_0 = 0`) is not used as a premise.

**8. Completeness / house template.** N11a rest-frame domain is present (`|q_{\mathrm{out}}·v_{\mathrm{bulk\_normal\_0}}/ω| ≪ 1`, large `k c_{s0}/|ω|` necessary not sufficient, no `v_0` alias). Chain output, T7 comparator, blind WL, supplied-vs-computed, F9/`N14` names, bind-closure `load_model` all present. FORM at d is an N2-permitted boundary refinement; the withheld bound stays orchestrator-side.

**Fold-specific.** The Fourier-reduction approach remains legitimate: name the object, supply `f̂_red`, put a leak-safe CONTROL in the spec, defer the element census and `(2π)` bookkeeping to the build directive (WL will not share PY dummy names). The element principle on the **rows** is complete enough as a spec control. The fold is incomplete in the sense of F1, not illegitimate.

---

## Verdict: **NOT-SOUND**

**Must-fix:** F1 — bind `J`, flux-normalized `S`/`C`/survival, the pole pencil, and §5a covariance to the same reduced rows that are already the construction operands of `𝓛`; kill the leftover “imported closed operator / consumed verbatim” construction path.

**Nits:** N1–N3 above.

Nothing else outstanding changes what an engine computes or what the spec may claim.
