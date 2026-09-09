# Codex-sol (gpt-5.6-sol xhigh) review — S11c-d SHARED PHYSICS spec v1 (DOCUMENT branch)

VERDICT: NOT-SOUND (9 must-fix). Raw transcript (~1.3 MB, echoed source files) trimmed to the review body.

---

NOT-SOUND. Nine must-fix findings survive. No files were modified.

## Must-fix findings

1. **The spec imports c2 objects that do not exist in the real export.**

   Spec: [§1a:78–85](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:78) supplies the closed operator “with `...TERM_ORIGINS` and `...PARITY_BLOCKS`,” the coupling kernel “with `...TERM_ORIGINS`,” and “the self-energy increment … and the six §3d re-adjudication objects.”

   Source: the c2 record says publication deliberately “[dropped] the increment to EMIT-only, [keeping] both closed operators” ([c2 record:59–61](/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:59)). The real delta contains only two high-level derived rows, `'s11cc2ClosedCouplingKernel'` and `'s11cc2ClosedSlabOperator'` ([export:11–12](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11)); there are no c2 term-origin, parity, increment, or six re-adjudication rows.

   This will either fail binding or tempt a builder to reconstruct supposedly supplied provenance.

   **Minimal fix:** supply the exact write keys `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` only. Treat everything else as step-record provenance, derive it anew as an S11c-d output, or first publish it through a separately reviewed upstream delta.

2. **The withdrawn c2 uniform-decoupling interpretation is quietly reinstated.**

   Spec first says correctly that F is withdrawn ([§1b:125–128](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:125)), but later asserts:

   > “With a uniform background the two sectors decouple identically”  
   > “the off-diagonal mixing must vanish”

   ([§2:186–188](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:186), [§5b:363–367](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:363)).

   Source: c2 explicitly says the corresponding “must vanish” wording is imprecise, the replacement is unsettled, and whether the closure-induced coupling decouples is exactly the withdrawn F question ([c2 record:166–171](/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:166)). Consistently, the real closed-kernel payload reports grades including `(ε¹,η⁰,σ⁰)` and `(ε¹,η⁰,σ¹)`, not only an `η¹` block ([export:11](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11)).

   **Minimal fix:** make the uniform specialization a literal computed output with no prescribed zero. Do not assert support only on `∇μ_R≠0` until the exported zero-background components are identified and adjudicated. Any subtraction of a uniform piece must itself be a named computed construction, not an assumed physical zero.

3. **The localized profile is incompletely and inconsistently bound to the inherited profiles.**

   Spec introduces

   > `W₀(x) = W̄₀[1+η f(x/L_W)]`

   and later says the varying quantities are “`W₀(x)`, `μ_R(x)`, `ρ_br⁰(x)`” ([§1c:134–149](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:134)).

   Source: the inherited ansatz instead has two independent profiles,

   > `W_bg = W_0[1+η w₁]`, `μ_R,bg = μ_R[1+η m₁]`

   with `W_0` and `mu_R` remaining constant ledger keys ([S11c-a:171–190](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:171)). The target neither binds `f=w1_profile` nor specifies the localized class or relation for `m1_profile`, even though it says the conversion is driven by `∇μ_R`. Its carrier list actually names constant `mu_R` but omits `mu_R_bg` and `m1_profile` ([§1a:86–88](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:86)).

   Two engines can consequently choose different `μ_R,bg` interfaces, or accidentally turn the constant `W_0`/`mu_R` keys into fields.

   **Minimal fix:** use the inherited fresh names and explicitly state `f≡w1_profile`; separately define `m1_profile`—either a stated function of `f` or an independently localized profile. Spell out both density representatives and their two asymptotic values. Reserve `W_0` and `mu_R` strictly for constants.

4. **The `η`/`σ_W` bookkeeping is incompatible with the advertised orders and truncation.**

   Spec calls `η`, `σ_W`, and `kL_W` “independent live grades,” then projects the response as `[...]_{εη}`, requires the operator/resolvent to be re-expanded to `(η≤1,σ≤1)`, yet emits a quadratic `O(η²)` observable ([§1d:160–180](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:160), [§2:193–209](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:193), [§3c:248–267](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:248)).

   Source: the inherited rule says `η` and `σ_W` must be multigraded independently and “no engine may … assign a common order” ([S11c-a:189–198](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:189)). The real kernel contains the four grades `(1,0,0)`, `(1,0,1)`, `(1,1,0)`, `(1,1,1)` ([export:11](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11)). Thus a literal `εη` projection can discard the `εσ_W` first-jet channel.

   Symbolically, after excluding the unresolved uniform piece,

   \[
   A_H=\epsilon(\eta a+\sigma_W b+\eta\sigma_W c),
   \]

   so

   \[
   J_{\rm conv}=\epsilon^2\!\left(
   \eta^2Q_{\eta\eta}
   +2\eta\sigma_WQ_{\eta\sigma}
   +\sigma_W^2Q_{\sigma\sigma}+\cdots\right),
   \qquad
   C=J_{\rm conv}/J_{\rm in},
   \]

   and the `ε²` cancels, but the background order is not simply `η²` while `η` and `σ_W` are independent. It becomes `O(η²)` only on a specified fixed-shape contrast homotopy where `σ_W∝η`.

   The spec also cannot apply an `η≤1,σ≤1` output truncation to a deliberately quadratic flux.

   **Minimal fix:** distinguish:

   - formal operator multigrades `(ε,η,σ_W)`;
   - the kinematic parameter `kL_W`—not a grade;
   - a one-parameter physical contrast homotopy, e.g. `(η,σ_W)→λ(η,σ_W)` at fixed `L_W`.

   Report amplitude grades explicitly, then allow quadratic observable grades. State `A=O(ελ)` and `C=O(λ²)` on that homotopy.

   The warning against `η=O(1)` is correct. The blanket “weak-gradient is forbidden” language is not: the imported operator is already a first-`σ_W` shape expansion. What must be forbidden is taking an additional `σ_W→0`/WKB limit or expanding away the full `kL_W` form factor.

5. **The interface scattering and flux normalization are not defined sufficiently to compute `C`.**

   Spec requests “outgoing-channel projections” and a “flux-normalized” fraction ([§3a:219–227](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:219), [§3c:250–267](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:250)) but supplies no flux bilinear form, unit-flux channel normalization, incidence side, threshold prescription, or left/right asymptotic modes.

   This matters especially because `W_+\ne W_-`: the two ends have different diagonal operators, wave numbers, and flux velocities. A step-like zero-jet profile is not a short-range perturbation of one global uniform operator. Consequently, the “plain uniform Born versus distorted-wave Born” choice cannot simply be deferred as equivalent ([§2:207–209](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:207)).

   Source: N7 stresses that the dimensionless observable is not fixed until “what couples to what” and its normalization are fixed ([decisions:106–112](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:106)).

   **Minimal fix:** define the two asymptotic operators `L_-` and `L_+`, incoming side(s), open outgoing channels, outgoing/Jost prescription, and the energy or symplectic flux pairing. Either normalize modes to unit flux or include the explicit `J_out/J_in` factor. For a genuine interface, make the two-asymptote distorted basis the construction; uniform Born may remain only a controlled expansion away from thresholds.

6. **The bound-mode claim misapplies the 1D weak-well theorem and leaks an existence result.**

   Spec states:

   > “a weak attractive well binds a mode in 1D, so this is a real photon-kill channel at small `η`”

   ([§3b:237–239](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:237)).

   The theorem applies to a localized attractive scalar self-adjoint potential. For

   \[
   H=-\partial_x^2-\eta V,\qquad V\ge0,\quad V\in L^1,
   \]

   one finds

   \[
   \kappa=\frac{\eta}{2}\int V\,dx+O(\eta^2),\qquad
   E_b=-\kappa^2
      =-\frac{\eta^2}{4}\left(\int V\,dx\right)^2+O(\eta^3).
   \]

   The chosen profile is instead an interface with different limits, not a localized well. A monotone scalar step generically has no state below both continua. Moreover, the actual thickness operator is multi-component, frequency-dependent, nonlocal, and potentially non-Hermitian because of the outgoing bulk response; no attractive-sign or self-adjointness premise is supplied.

   There is a second contradiction: the new bound pole arises from resumming

   \[
   G=(G_0^{-1}-\eta V)^{-1}.
   \]

   Re-expanding `G` only through `O(η)` as required by [§2:201–203](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:201) cannot create that new pole. Plain uniform Born has no bound channel at all.

   Finally, in stationary conservative scattering, a true bound state carries no asymptotic outgoing flux; a “capture rate” requires a wave-packet/switching protocol, damping, or a resonance width. A pole residue and coupling alone are not a conversion probability.

   **Minimal fix:** emit a conditional pole/existence test with no prescribed existence or sign. Distinguish true bound poles from second-sheet resonances. Define the residue with left/right modes and the nonlinear-eigenvalue normalization involving `∂_ωL_H`; define a capture protocol or report only spectral overlap/coupling. If a resummed first-order effective operator is used, explicitly exempt that spectral solve from the continuum Born re-expansion and state its error limitation. “Not a Bloch band” is correct if a pole is actually found.

7. **`F′(0)` is the wrong weak coefficient for a conversion fraction.**

   Spec says S11c-d delivers “`F′(0)` (the `O(η)` slope of the conversion)” while the lab bounds `F(1)` ([§3d:272–285](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:272)).

   Its own order table says the conversion fraction is `O(η²)`. If

   \[
   A_H/\epsilon=\eta a_1+O(\eta^2),
   \qquad
   C(\eta)=\gamma |A_H/\epsilon|^2
           =\gamma|a_1|^2\eta^2+\cdots,
   \]

   then

   \[
   C'(0)=0,\qquad
   \lim_{\eta\to0}\frac{C}{\eta^2}
   =\frac12C''(0)=\gamma|a_1|^2.
   \]

   The spec’s counterexample proves the error directly:

   \[
   C=\sin^2(\eta G),\quad C'(0)=0,\quad \tfrac12C''(0)=G^2.
   \]

   If `F` instead denotes amplitude, `F′(0)` is legitimate, but the lab bounds the flux fraction `|F(1)|²` with flux factors—not `F(1)`.

   **Minimal fix:** emit separately the amplitude coefficient
   `∂η(A_H/ε)|₀` and the fraction coefficient `lim C/η² = C''(0)/2`. Name the strong observable `C_strong(1)`. The counterexample and the “no positive lower bound” conclusion are otherwise correct.

8. **The d-level N6 control omits the covariance comparison that made c2’s Reading B valid.**

   Spec defines only

   > `S11CD_REP_INVARIANCE_RESIDUAL = route1 − route2`

   and says no separate transform is applied ([§5a:329–360](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:329)).

   Source: c2 says the native material operand is **not** the fully `Φ`-transformed image, that the raw `R_N6` is nonzero `18/288`, and that full frame-change faithfulness is checked separately by `R_cov` ([c2 record:145–155](/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:145)). The governing disposition likewise says the matched covariance zeros are `(0)−(0)`, not operand agreement ([disposition:49–59](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_disposition.md:49)).

   The target correctly preserves those facts in §1b, but its own control reverts to the insufficient raw difference and provides no response-level transformation law for the source and asymptotic channel bases.

   **Minimal fix:** emit both a raw component residual analogous to `R_N6` and a separate d-level covariance/naturality residual after transforming the source, incoming/outgoing modes, and flux pairing under `Φ`. Do not treat raw equality as N6 closure. Keep the carrier/source/Φ operand debt explicitly unclosed.

   The tilt/advection one-sided probes, `RHO4_CONSTANT` structural absence with no `A−A`, rejection of `∇W₀→0`, and exclusion of anchoring corruption are otherwise correctly stated.

9. **N15 is applied at the wrong layer and invites new upstream physics.**

   Spec tells d to:

   > “emit any new gradient-of-background invariants as RESULTS (new constants if they appear) — e.g. the edge form factor / integrated-jump moment”

   ([§1c:154–156](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:154), [§4:314](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:314)).

   Source: N15 assigns that constitutive construction to the variable-coefficient operator stage: inherit uniform invariants and emit newly admitted gradient invariants rather than smuggling them into the energy ([decisions:143–148](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:143)). S11c-d is supposed to consume c2’s operator/kernel verbatim; it cannot invent a new constitutive constant without changing upstream physics. A Fourier form factor or integrated jump is profile/scattering data, not a new local constitutive invariant.

   **Minimal fix:** emit profile moments and form factors derived from the imported kernel, with no new constitutive constants. If the scattering calculation exposes a missing local invariant, record it as an upstream N15 debt rather than adding it inside d.

## Nits

- [Line 3](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:3) calls d the “fourth sub-step.” It is the fourth lettered seam, but after the c1/c2 split it is the fifth build unit. Clarify the label.

- In the counterexample, [line 280](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:280) says `C → 0 at finite ηG`; more precisely, `C=0` at finite nonzero `ηG=nπ`.

- [§5c:373–378](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:373) asks for an arbitrary full “edge minus bump” object. That difference is not canonical without choosing two profiles. The invariant discriminant is the zero-transfer limit:
  \(\widehat{W'}(0)=W_+-W_-\), nonzero for an interface and zero for a bump.

## Items that checked out

The interface-versus-bump distinction itself is correct; excluding `η=O(1)` inside the first-shape operator is correct; the `ε²` cancellation in a properly normalized linear conversion fraction is correct; the two-mode counterexample does prove that a nonzero Born coefficient supplies no strong-edge lower bound; and the finite-contrast edge is correctly named as a new downstream construction. The c2 Path-B counts/status, withdrawn F/G status, N11 rest-frame condition, blind-Wolfram rule, chain topology, comparator discipline, and withholding of the numeric/O(1) falsification target are all preserved correctly.

**Verdict: NOT-SOUND — must-fix findings 1–9 remain.**
