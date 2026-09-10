# Codex-sol (gpt-5.6-sol xhigh) review — S11c-d SHARED PHYSICS spec v3 (round 3, DOCUMENT branch)

VERDICT: NOT-SOUND (10 must-fix); structural layer confirmed sound. Raw transcript (~2.4 MB) trimmed to the review body.

---

The spec is **NOT-SOUND**. The core profile choice and strong-edge reasoning are sound, and the c2 status is reported honestly, but ten issues can change the constructed scattering object, its normalization, or what the result may claim.

## Must-fix findings

### 1. The distorted-wave reduction assumes the withdrawn uniform-decoupling result

**Spec text.** Section 2 says:

> “uniform-mode Born and the two-asymptote distorted-wave construction coincide (`G = G₀ + O(η,σ_W)`, `K = O(first shape)` …)”  
> — [S11c_d_SHARED_PHYSICS.md:234](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:234)

But the same spec correctly says:

> “The uniform-background amplitude is a computed object `A_0`, not an assumed zero.”  
> — [S11c_d_SHARED_PHYSICS.md:220](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:220)

The real closed-kernel export explicitly carries grades

```text
((1,0,0), (1,0,1), (1,1,0), (1,1,1))
```

including an `(ε¹,η⁰,σ_W⁰)` component; see [S11c_c2_exports.py:11](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11). Whether its apparent cancellations establish uniform decoupling is precisely the withdrawn F question:

> “whether ‘the genuine closure-induced coupling decouples’ holds is exactly the F question … do not assert it here.”  
> — [S11c_c2_self_energy_fold.md:166](/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:166)

**Why wrong.** Write

\[
G=G_0+\lambda G_1+\cdots,\quad
K=K_0+\lambda K_1+\cdots,\quad
\psi=\psi_0+\lambda\psi_1+\cdots .
\]

Then

\[
G K\psi
=G_0K_0\psi_0
+\lambda\bigl(G_1K_0\psi_0+G_0K_1\psi_0+G_0K_0\psi_1\bigr)+\cdots .
\]

The uniform-mode matrix element retains only the middle term. Coincidence therefore requires the unearned condition `K₀=0`. Likewise, if the constant asymptotic off-diagonal blocks `K_±` are nonzero, the correct asymptotic channels are eigenchannels of the full block operators, not of the diagonal `L_±` used at [line 226](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:226).

**Minimal fix.** Decompose and compute `K₀`, `K₋`, and `K₊` first. Use the full asymptotic block operators to define channels. Permit the diagonal-sector/uniform-mode simplification only conditionally if the relevant computed off-diagonal baselines vanish.

---

### 2. `ΔA=A−A₀` is not necessarily the gradient amplitude, and its squared norm is not the physical excess flux

**Spec text.**

> “The leakage-relevant object is the induced amplitude `ΔA ≡ A − A₀` (the part sourced by the gradient).”  
> — [S11c_d_SHARED_PHYSICS.md:308](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:308)

It then assigns

> “the induced converted flux `O(ε²λ²)` … `C=J_conv/J_in=O(λ²)`”  
> — [S11c_d_SHARED_PHYSICS.md:311](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:311)

while expressly declining to assume `A₀=0`.

**Why wrong.** The exported multigrade permits

\[
A=\epsilon\bigl(a_{00}+\eta a_{10}+\sigma_W a_{01}
+\eta\sigma_Wa_{11}\bigr).
\]

Thus `A−A₀` includes the zero-jet contrast term `ηa₁₀`; it is not necessarily “the part sourced by the gradient.”

More importantly, for flux pairing \(B\), along \(\sigma_W=c\lambda\),

\[
J[A]=\epsilon^2\!\left[
 B(a_0,a_0)
 +2\lambda\operatorname{Re}B(a_0,a_1)
 +\lambda^2B(a_1,a_1)+\cdots
\right].
\]

Therefore the physical total conversion can contain baseline `O(ε²)` and interference `O(ε²λ)` terms. Although

\[
B(\Delta A,\Delta A)=O(\epsilon^2\lambda^2),
\]

it is neither \(J[A]\) nor generally \(J[A]-J[A_0]\). The `ε²` cancellation correctly gives an `O(λ²)` fraction only when the relevant baseline amplitude vanishes, is orthogonal, or the observable is explicitly defined as the self-flux of a coherently subtracted field.

**Minimal fix.**

- Separate zero-jet (`η`) and first-jet (`σ_W`) amplitude components.
- Emit the total flux expansion, including baseline and interference.
- Rename \(B(\Delta A,\Delta A)/J_{\rm in}\) as an induced-field quadratic form, not the physical total conversion fraction.
- Apply the mandated `O(ε²η²)`/`O(η²)` physical-leakage labels conditionally on the computed baseline/interference disposition.

The basic N12 counting itself is correct for a genuinely `O(εη)` converted field:

\[
A_{\rm conv}=O(\epsilon\eta),\quad
J_{\rm conv}\propto |A_{\rm conv}|^2=O(\epsilon^2\eta^2),\quad
J_{\rm in}=O(\epsilon^2),\quad
C=O(\eta^2).
\]

The defect is identifying that quantity with the observable while `A₀` remains unresolved.

---

### 3. The interface Fourier object lacks the decay, distribution, and normalization premises needed by both engines

**Spec text.**

> “smooth, asymptotically constant, with a finite integrated jump”  
> — [S11c_d_SHARED_PHYSICS.md:38](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:38)

and

> “one reduced 1-D Fourier convention along `n̂`”  
> — [S11c_d_SHARED_PHYSICS.md:176](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:176)

No convention is actually supplied.

**Why wrong.**

1. Asymptotic constancy and a finite conditional jump do not ensure \(w_1',m_1'\in L^1\), which the invoked Riemann–Lebesgue argument requires. The exported kernel also contains higher profile jets, so their necessary decay/domain must be stated.
2. The full profiles \(w_1,m_1\) do not decay for an interface with unequal asymptotes. Their Fourier transforms are distributions, not ordinary localized form factors.
3. The actual c2 carrier uses a three-dimensional `(2π)⁻³` convention. For a profile uniform parallel to the edge,

\[
\frac1{(2\pi)^3}\int d^3y\,L_W\partial_n f\,
 e^{-iQ\cdot y}
=
\frac{L_W}{2\pi}\delta^{(2)}(Q_\parallel)
\int d\xi\,e^{-i(Q_nL_W)\xi}f'(\xi).
\]

The spec does not state whether the tangential delta and `L_W/(2π)` are retained, stripped, or converted to a per-unit-edge-area amplitude. Blind engines can therefore differ by dimensional and \(2\pi\) factors.

**Minimal fix.** Require sufficient short-range decay, at least for every profile derivative consumed by the retained kernel; define the exact reduced transform; specify the mapping from every c2 three-dimensional Fourier carrier; and state the per-unit-area/tangential-delta convention. Treat zero-jet step profiles in coordinate space or after an explicit asymptotic subtraction rather than as ordinary \(L^1\) Fourier functions.

---

### 4. The `Q_nL_W` “sudden/WKB/maximal” claims are too strong, and the stated `nπ` nodes are not generic

**Spec text.**

> “`Q_nL_W→0` is the zero-transfer / sudden limit … integrated jump, maximal conversion; `|Q_nL_W|≫1` … is the WKB / adiabatic regime”  
> — [S11c_d_SHARED_PHYSICS.md:203](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:203)

and

> “a form-factor node `Q_nL_W=nπ` with `Δm₁≠0`”  
> — [S11c_d_SHARED_PHYSICS.md:339](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:339)

**Why wrong.**

- `Q_nL_W→0` can arise from `Q_n→0` at any width; it is not by itself a sudden-edge limit. Taking `L_W→0` at fixed contrast would also drive `σ_W` large, outside the imported first-shape expansion.
- At zero transfer, \(\widehat{m_1'}(0)=\Delta m_1\) under the specified unnormalized transform, but this does not prove maximal total conversion. The full vertex, flux weights, other profile channels, and nonmonotone derivatives can move the maximum elsewhere.
- Large transfer suppresses an \(L^1\) derivative form factor, but “WKB” requires a local wavelength/gap adiabaticity condition; `Q_nL_W` alone is insufficient, especially near a mode crossing.
- `nπ` zeros depend on profile form. For the spec’s representative \(m_1=(1+\tanh\xi)/2\),

\[
\widehat{m_1'}(s)=\frac{\pi s}{2\sinh(\pi s/2)},
\]

which has no nonzero real `nπ` nodes.

**Minimal fix.** Call `Q_nL_W→0` only the zero-transfer limit. State that large `|Q_nL_W|` suppresses this particular derivative-form-factor contribution under the new decay assumptions. Define WKB separately through local modal wave numbers/gaps. Replace `Q_nL_W=nπ` by “a zero of the computed form factor, if one exists.”

---

### 5. The admissibility paragraph freezes the wrong density in one branch

**Spec text.**

> “Name which quantities vary (`W_bg`, `μ_R,bg`, `ρ_br,bg⁰` — both density representatives …)”  
> — [S11c_d_SHARED_PHYSICS.md:179](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:179)

**Source text.**

> `RHO4-CONSTANT`: `ρ_4D,bg⁰` is constant and `ρ_br,bg⁰=ρ_4D,bg⁰W_bg`;  
> `RHOBR-CONSTANT`: `ρ_br,bg⁰` is constant and `ρ_4D,bg⁰=rho_br/W_bg`.  
> — [S11c_a_SHARED_PHYSICS.md:210](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:210)

**Why wrong.** `ρ_br,bg⁰` does not vary in the `RHOBR_CONSTANT` branch; `ρ_4D,bg⁰` does. The present wording omits the latter varying field exactly where the N4 advection probe depends on its gradient.

**Minimal fix.** State the branchwise maps verbatim and carry both density asymptotes/gradients live.

---

### 6. Section 5a does not provide a constructible independent N6 route or separable corruption sites

**Spec text.** It first says:

> “S11c-d does not re-derive the material closed operator … no face-normal carrier factory is imported.”  
> — [S11c_d_SHARED_PHYSICS.md:400](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:400)

but then requires:

> “the SAME insertion built after flattening the interface faces to material in-plane coordinates”  
> — [S11c_d_SHARED_PHYSICS.md:407](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:407)

and separate mutations of the tilt and advective-density “atoms” in the imported kernel.

The governing requirement is:

> derive by direct graph linearization, derive again after flattening into material coordinates, transform exactly, compare, then corrupt one route only.  
> — [S11c_decisions.md:94](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:94)

**Why wrong.** A coordinate rewrite of the same already-constructed imported kernel is not an independent shape/coordinate derivation of its N3/N4 content. Moreover, §1a correctly admits that no c2 term-origin rows are exported. The final kernel contains common `w₁′` factors but no exported provenance that isolates “tilt” from `u·∇ρ₄/ρ₄`; flipping `w₁′` flips multiple channels. The SymPy route therefore has neither the independent material carrier nor a constructible channel-specific mutation site.

The correct details—rejecting `∇W→0`, not corrupting an anchoring, and emitting structural absence rather than `A−A` for `RHO4_CONSTANT`—are present, but they do not make the proposed routes independent.

**Minimal fix.** Either:

- import/reconstruct explicit Eulerian and native-material pre-extraction operands with channel provenance and an exact coordinate/Jacobian map, then mutate those source operands; or
- relabel this as a downstream scattering-coordinate covariance regression, remove the claim that it discharges N6/N3/N4 independence, and carry the c2 N6 debt unchanged.

---

### 7. The jet-zero regression cannot retain unequal asymptotes

**Spec text.**

> “sets the profile jets to zero … while retaining live `η` and the arbitrary constant asymptotes `W_±`, `μ_±`”  
> — [S11c_d_SHARED_PHYSICS.md:444](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:444)

**Derivation.** If \(w_1'(\xi)=0\) everywhere, then \(w_1\) is one constant and

\[
\Delta w_1=\int_{-\infty}^{\infty}w_1'(\xi)d\xi=0.
\]

The same holds for `m₁`. A single smooth profile cannot have zero jets and independently retained unequal left/right asymptotes.

**Minimal fix.** Run two separate uniform regressions: one with the left constant background everywhere and one with the right constant background everywhere. In each run the two ends are equal. Alternatively retain a discontinuous jump, but then its derivative is a delta and the jets are not zero.

---

### 8. The supposedly canonical S-matrix, flux normalization, and confinement object remain noncanonical

**Spec text.**

> “one named incident end (or the complete channel matrix …)”  
> — [S11c_d_SHARED_PHYSICS.md:259](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:259)

and

> “the amplitude carries the explicit `√(v_out/v_in)` factor”  
> — [S11c_d_SHARED_PHYSICS.md:261](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:261)

and

> “`J_conv` = the continuum (thickness / bulk) outgoing transverse-channel flux”  
> — [S11c_d_SHARED_PHYSICS.md:319](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:319)

**Why wrong.**

- “One end or full matrix” permits different emitted objects and labels.
- For a multi-component, frequency-dependent, potentially non-Hermitian operator, field-amplitude flux is not generally just group velocity. It requires the mode’s bilinear current/energy normalization, often involving left/right modes and `∂ωL`.
- A converted thickness/bulk channel is not an “outgoing transverse-channel flux.” The incident flux is transverse; the converted flux uses the converted channel’s own current or the bulk-loss functional.
- N13 requires survival of the transverse channel, but the spec emits only a vague `S11CD_CONFINEMENT_CONDITION`, not the reflected-plus-transmitted transverse survival block/fraction. The source requirement is explicit:  
  > “confinement … means survival of the transverse polarization channel.”  
  > — [S11c_decisions.md:130](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:130)

**Minimal fix.** Require the complete left/right channel matrix; supply or derive an explicit modal flux bilinear for every asymptotic channel; define \(C_{T\to H}=J_{H,\mathrm{out}}/J_{T,\mathrm{in}}\); and emit the transverse survival functional from all outgoing transverse reflection/transmission channels.

---

### 9. The bound-channel output is method-dependent, optional, and overclaims the truncated spectral solve

**Spec text.**

> “emit a profile-functional Jost / Evans determinant and its zeros”  
> — [S11c_d_SHARED_PHYSICS.md:280](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:280)

and

> “`S11CD_CAPTURE_PROTOCOL (or spectral-overlap-only)`”  
> — [S11c_d_SHARED_PHYSICS.md:297](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:297)

Section 2 also promises that the resummed solve “carries its own stated error limitation,” but §3b supplies no such truncation limitation.

**Why wrong.**

- An Evans/Jost determinant is defined only up to multiplication by a nonvanishing analytic factor: \(D_2(\omega)=h(\omega)D_1(\omega)\). Its zeros agree, but its literal value and derivative need not. It is therefore not a canonical comparator payload unless its normalization is fixed.
- The imported operator is nonlocal; an ordinary finite-dimensional Evans construction is not automatically available. A Fredholm determinant would itself need compactness/trace-class premises.
- “Protocol or overlap” allows the two engines to compute different physical objects. A stationary bound state has no outgoing flux, as the spec correctly notes, so a capture probability cannot exist without a fixed preparation/switching/damping protocol.
- Resumming \(L^{(1)}=L_0+\lambda V_1\) does not turn it into an all-orders operator. Poles found near thresholds can move, appear, or disappear under the omitted terms. They are poles of the retained operator unless a separation-from-threshold/error condition is supplied.

**Minimal fix.** Make the canonical object the pole set of the retained resolvent plus normalized Riesz residue/projector, with \(\langle l,\partial_\omega L\,r\rangle=1\) where applicable. Choose spectral overlap only unless a concrete capture protocol is supplied. State the retained-operator truncation error and the threshold domain in which a pole classification may be promoted to a physical claim.

The spec is otherwise right not to invoke the scalar weak-1D-well theorem class-wide, and right that this is not a Bloch band.

---

### 10. The d comparator cannot “surface” the c2 carrier/source/Φ debt merely by comparing the final amplitude

**Spec text.**

> “the comparator SURFACES the DEBT and the §1b representation questions”  
> — [S11c_d_SHARED_PHYSICS.md:533](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:533)

**Source text.**

> “`R_N6` itself, the channels, the reconcile-engine sources, and the guards are ENTIRELY UNMATCHED … a SCHEMA non-join.”  
> — [S11c_c2_N6_reconcile_disposition.md:70](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_disposition.md:70)

**Why wrong.** A residual on a projected final scattering amplitude can reveal a downstream disagreement, but it cannot separately expose the 40 carrier, 76 source, and 18 Φ operand residuals. Projection and channel summation may cancel them. Agreement of the amplitude would close only that particular projected d object, not the upstream operands.

**Minimal fix.** Say the d amplitude remains conditional on and propagates the debt. Claim direct surfacing only if d emits separately defined, common-basis projections of those operand families with a sound join.

## Nits

- The citation for "`w₁` and `m₁` are independent" points to S11c-a `:190`, which discusses independence of `η` and `σ_W`; the relevant profile definition begins at [S11c_a_SHARED_PHYSICS.md:171](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:171).
- “Both residuals” at [S11c_d_SHARED_PHYSICS.md:423](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:423) has no clear antecedent because the naturality residual immediately above is optional.

## Checks that passed

- A localized interface as a derivative-localized, unequal-asymptote object is the correct non-global N5 choice; `∫W_bg′≠0` versus zero for a bump is the correct topological distinction, once the transform convention and decay conditions are fixed.
- Contrast `η`, sharpness grade `σ_W`, and kinematic argument `Q_nL_W` are distinct. The prohibition on taking an additional `σ_W→0` construction limit is appropriate; setting `η=O(1)` inside the truncated c2 operator is invalid.
- The strong-edge section correctly withholds finite-contrast conversion. The counterexample is valid:
  \[
  \sin^2(\eta G)=\eta^2G^2+O(\eta^4),\qquad
  \sin^2(n\pi)=0,
  \]
  so a nonzero Born coefficient supplies no positive strong-edge lower bound.
- The c2 per-engine/cross-engine status, `(0)−(0)` caveat, `18/288` raw residual, `R_cov` no-nonzero result, operand counts, uninspected shape, withdrawn F/G, and materiality of the debt are preserved honestly.
- The `RHO4_CONSTANT` structural-absence rule, rejection of `∇W→0` as N6, exclusion of anchoring corruption, N11 rest-frame condition, chain topology, blind-Wolfram rule, and withholding of the numeric/`O(1)` falsification target are correct.

**Verdict: NOT-SOUND — must-fix findings 1–10 remain.**
