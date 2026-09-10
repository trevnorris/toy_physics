# Codex-sol (gpt-5.6-sol xhigh) review — S11c-d SHARED PHYSICS spec v2 (round 2, DOCUMENT branch)

VERDICT: NOT-SOUND (11 must-fix). Raw transcript (~379 KB, echoed source files) trimmed to the review body.
Reviewed text = spec v2 (quoted by line below); v2 was overwritten by v3 in-tree without a separate baseline commit.

---

I found 11 must-fix issues. Several are fold-introduced internal contradictions: the spec preserves c2’s unresolved uniform component, but then silently assumes it vanishes when defining the mixing support and conversion order.

## Must-fix findings

1. The selected class is not actually imposed on the thickness profile

Spec §1c says:

> “the localized-interface condition … is imposed on the mixer `m₁` … and named for `w₁`,” while `w₁` “may itself be an interface or a co-present tilt.”  
> — [S11c_d_SHARED_PHYSICS.md:157](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:157>)

But N5 requires S11c-d to name a profile class for the varying thickness \(W_0(x)\):

> “S11c-d must name a profile class … localized→Born/scattering kernel.”  
> — [S11c_decisions.md:83](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:83>)

For the inherited ansatz,
\[
\int_{-\infty}^{\infty}\partial_yW_{\rm bg}\,dy
=\bar W_0\eta\,[w_1(+\infty)-w_1(-\infty)].
\]
Thus allowing \(w_1\) to be a bump permits zero thickness jump and \(W_-=W_+\), contradicting §2’s supplied assertion:

> “the interface has different asymptotes `W₋ ≠ W₊`.”  
> — [S11c_d_SHARED_PHYSICS.md:202](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:202>)

Minimal fix: require \(\Delta w_1\neq0\) for the thickness interface. If the modulus mixer must also be an interface, separately require \(\Delta m_1\neq0\), while retaining functional independence \(m_1\not\equiv m_1[w_1]\).

2. The interface geometry, Fourier normalization, momentum transfer, and WKB limit are not correctly fixed

The spec inherits vector \(y,\xi\), writes component derivatives, but then uses one-dimensional notions \(m_1(\pm\infty)\), \(\int\partial_\xi d\xi\), and a Jost problem without selecting an interface normal:

> “`ξ ≡ y/L_W` … `∂_{yᵢ}W_bg` … `m₁(−∞) ≠ m₁(+∞)`.”  
> — [S11c_d_SHARED_PHYSICS.md:152](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:152>)

That leaves the two engines free to choose different reductions and Fourier conventions. The real c2 Fourier object uses the full transfer
\[
\mathbf Q=\mathbf k_{\rm out}-\mathbf k_{\rm in},
\]
not a single unspecified \(k\): its exported binding is of the form  
`s11cc2FourierW1ProfileHatTransfer(-k_input + k_output, …)`.  
— [S11c_c2_exports.py:11](</var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11>)

The spec instead repeatedly supplies \(\widehat{m_1'}(kL_W)\), and states:

> “sending it to `q=0` is WKB.”  
> — [S11c_d_SHARED_PHYSICS.md:191](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:191>)

That is backwards. For the representative \(m_1=(1+\tanh\xi)/2\),
\[
\widehat{m_1'}(s)=\frac{\pi s}{2\sinh(\pi s/2)}.
\]
It tends to the integrated jump at \(s=QL_W\to0\); the adiabatic/WKB limit at fixed nonzero \(Q\) is \(L_W\to\infty\), hence \(|QL_W|\to\infty\), where conversion is suppressed.

Minimal fix: define a unit normal \(\hat n\), scalar \(\xi=\hat n\!\cdot y/L_W\), tangential momentum conservation, \(Q_n=k_{{\rm out},n}-k_{{\rm in},n}\), and one reduced one-dimensional Fourier convention. Then write the form factor as \(\widehat{m_1'}(Q_nL_W)\). State that \(Q_nL_W\to0\) is the zero-transfer/sudden limit and \(|Q_nL_W|\gg1\) is the WKB/adiabatic regime.

3. The c2 import manifest still contains a nonexistent key and ambiguously consumes a non-exported object

The spec names:

> “the field carriers `s11cc2FieldTheta`, `s11cc2FieldeW`, …”  
> — [S11c_d_SHARED_PHYSICS.md:92](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:92>)

The real key is lowercase `theta`:

> `'s11cc2Fieldtheta': …`  
> — [S11c_c2_exports.py:19](</var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:19>)

The spec also says the operator, kernel, and “increment VALUES” are followed by:

> “These are the operands S11c-d consumes.”  
> — [S11c_d_SHARED_PHYSICS.md:108](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:108>)

But its own §1a says no self-energy-increment row exists, and c2 records that the increment was dropped from publication/export:

> “drop the increment to EMIT-only, keep both closed operators.”  
> — [S11c_c2_self_energy_fold.md:56](</var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:56>)

Minimal fix: replace `s11cc2FieldTheta` with `s11cc2Fieldtheta`, and state explicitly that d consumes only the two exported operator/kernel rows plus actual carrier rows; the increment value is established audit provenance, not an import operand.

4. The c2 `R_N6` result is placed under a cross-engine heading without its essential “unmatched” qualification

Under:

> “CROSS-ENGINE, dual-engine confirmed,”

the spec preserves:

> “`R_N6 … nonzero (18/288)` AND `R_cov` no-nonzero.”  
> — [S11c_d_SHARED_PHYSICS.md:112](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:112>)

The governing disposition says:

> “`R_N6` itself … [is] ENTIRELY UNMATCHED” and “there is NO direct cross-engine comparison of `R_N6` itself.”  
> — [S11c_c2_N6_reconcile_disposition.md:70](</var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_disposition.md:70>)

Minimal fix: label \(18/288\) explicitly as the per-engine SymPy raw result. Only the Reading-B covariance vanishing statement and specified control/premise subset are dual-engine confirmed; the raw \(R_{N6}\) object was schema-unmatched.

5. The spec reinstates the withdrawn uniform-decoupling claim and drops non-\(\nabla\mu_R\) channels

The supplied framing says:

> “mixing … supported only where `∇μ_R,bg ≠ 0`.”  
> — [S11c_d_SHARED_PHYSICS.md:202](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:202>)

But §1a correctly describes the imported full vertex as containing:

> “tilt `∇w₁`, modulus-gradient `∇m₁`, and N4 advection channels; ⛔ not `∇w₁` alone.”  
> — [S11c_d_SHARED_PHYSICS.md:89](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:89>)

Since \(w_1\) and \(m_1\) are independent, tilt/advection support need not coincide with \(\nabla\mu_R\). More importantly, “supported only” logically supplies zero mixing in the uniform limit, while the governing c2 record says F is withdrawn:

> “F/G conclusions are withdrawn … the increment VALUES stand.”  
> — [S11c_c2_self_energy_fold.md:198](</var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:198>)

Minimal fix: define support from the full imported kernel—tilt, modulus-gradient, and advection—and remove every “only where \(\nabla\mu_R\neq0\)” statement unless explicitly naming the modulus subchannel. Keep the uniform component as computed.

6. The conversion orders are asserted for the total amplitude even though its uniform component is explicitly unresolved; the required \(\eta\)-labels are also omitted

The spec says:

> “the amplitude is `O(ελ)`, the absolute converted flux `O(ε²λ²)` … `C=O(λ²)`,” and “do not report a single `O(εη)`/`O(ε²η²)` collapse.”  
> — [S11c_d_SHARED_PHYSICS.md:290](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:290>)

But it also says the imported kernel has \(\eta\)-zero grades and that the uniform amplitude must be computed. In general,
\[
A(\lambda)=\epsilon(A_0+\lambda A_1+\cdots)
\]
gives
\[
J_{\rm conv}\propto\epsilon^2\!\left(
|A_0|^2+2\lambda\operatorname{Re}A_0^*A_1
+\lambda^2|A_1|^2+\cdots\right).
\]
Thus the quoted orders hold only if \(A_0=0\), exactly the withdrawn/unadjudicated premise.

When \(A_0=0\), the requested physical fixed-\(L_W\) counting is straightforward:
\[
A_{\rm conv}=\epsilon\eta a_1,\quad
J_{\rm conv}=v_{\rm out}\epsilon^2\eta^2|a_1|^2,\quad
J_{\rm in}=v_{\rm in}\epsilon^2,
\]
so
\[
C=\frac{v_{\rm out}}{v_{\rm in}}\,\eta^2|a_1|^2=O(\eta^2).
\]
This is exactly N12:

> “coupling is `O(εη)`; leakage probability/rate `O(ε²η²)`.”  
> — [S11c_decisions.md:123](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:123>)

Minimal fix: emit the full multigrade and \(A_0\). Apply the \(O(\epsilon\eta)\), \(O(\epsilon^2\eta^2)\), and \(O(\eta^2)\) labels either conditionally after computed \(A_0=0\), or to a separately named induced amplitude \(\Delta A=A(\lambda)-A(0)\). Retain \(\lambda\) only as the joint physical homotopy explanation, not as a replacement for N12’s required \(\eta\)-labels.

7. The continuum S-matrix object is not canonical enough for blind comparison

The spec says:

> “name the INCOMING side,”

rather than selecting one or emitting both, and permits amplitudes:

> “normalized to unit flux (or carrying the explicit `J_out/J_in` factor).”  
> — [S11c_d_SHARED_PHYSICS.md:208](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:208>) and [S11c_d_SHARED_PHYSICS.md:240](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:240>)

These alternatives are not the same amplitude object. If \(a_{\rm raw}\) is field-normalized, the flux-normalized amplitude generally contains \(\sqrt{v_{\rm out}/v_{\rm in}}\). Likewise left and right incidence produce different reflection/transmission blocks.

Minimal fix: either emit the complete channel matrix with both incident ends, or fix one end explicitly. Choose one amplitude normalization and branch labels `(incoming end, outgoing end, channel, ω, k_parallel)`; emit raw amplitudes separately if desired.

8. The bound channel is still underdetermined and internally inconsistent

The spec correctly rejects applying the weak-well theorem blindly, but then demands:

> “solve the thickness-diagonal spectral problem … PERMIT COMPUTED ABSENCE.”  
> — [S11c_d_SHARED_PHYSICS.md:258](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:258>)

No concrete profile is supplied—indeed §1c says a specific shape is only a representative. A profile class with fixed asymptotic jump contains both monotone steps and profiles with added localized overshoots/wells; they can have different pole sets. Therefore existence/absence cannot be computed class-wide.

There are three further defects in the same object:

- “below both continua” is insufficient for a true bound state when the closed operator contains a non-Hermitian outgoing-bulk self-energy. A true bound pole must be a physical-sheet, normalizable, zero-width pole with every radiation channel closed.
- The statement  
  > “`Not a Bloch band` holds iff a pole is found”  
  is false. A localized interface is nonperiodic regardless of whether its pole set is empty. If found, the pole is a localized discrete/interface state, not a Bloch band.
- §3b correctly says a pole residue alone is not a capture probability and permits “spectral-overlap-only,” but §3c nevertheless defines  
  > “`J_conv` … plus bound capture.”  
  — [S11c_d_SHARED_PHYSICS.md:265](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:265>) and [S11c_d_SHARED_PHYSICS.md:285](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:285>)

Minimal fix: either supply a concrete profile and capture protocol, or emit a profile-functional Jost/Evans determinant plus left/right residue and spectral coupling, without claiming computed existence. Add a bound contribution to \(J_{\rm conv}\) only when a specified switching/wave-packet/damping protocol produces a dimensionally valid capture rate.

9. §5a does not implement the N6 independent route, and its covariance formula risks double-applying \(\Phi\)

The spec says:

> “Both routes construct … from the imported CLOSED operator/kernel” and explicitly rejects “a fresh face level-set/graph linearization.”  
> — [S11c_d_SHARED_PHYSICS.md:365](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:365>)

N6 requires exactly:

> “derive … by direct level-set/graph linearization, derive it again after flattening faces into material coordinates … then corrupt one route only.”  
> — [S11c_decisions.md:94](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:94>)

The spec also incorrectly identifies the direct graph route with c2’s open `close(extract)` ordering. c2 defines `close(extract)` only as an ablation:

> “`close(extract(SLAB))` … is the §5a ablation, ⛔ not a construction route.”  
> — [S11c_c2_SHARED_PHYSICS.md:150](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:150>)

The real c2 delta exports no N6 route operands; its N6 constructions live in separate diagnostic streams. Consequently the declared d consume-set cannot construct the material route or source-level tilt corruption.

Finally, route 2 is said to be “already in the common Eulerian face basis,” but the next formula uses `route1 − Φ(route2)`. c2 explicitly says no separate final transform is applied to the already-mapped increment:

> “NO separate `T`/final pullback on the increment.”  
> — [S11c_c2_self_energy_fold.md:145](</var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md:145>)

Minimal fix: supply/import the necessary face/source primitives and independently build
\[
\operatorname{extract}(\operatorname{close}SLAB_E),\qquad
\operatorname{extract}(\operatorname{close}SLAB_M)
\]
at fixed anchoring and density. Map the material covector once inside its native builder. Define the separate covariance/naturality residual at the source, mode, and flux-pairing level; do not write a second global \(\Phi\) on an already Eulerian-basis output.

10. The uniform regression is vacuous for the gradient-independent term it claims to test

The spec sets:

> “`W_bg→W̄₀` (`η→0`, `σ_W→0`)”  
> and calls this a smoke test for a forbidden gradient-independent term.  
> — [S11c_d_SHARED_PHYSICS.md:406](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:406>)

N6 says the uniform regression should catch a forbidden gradient-independent term:

> “useful smoke test for a forbidden gradient-independent term.”  
> — [S11c_decisions.md:94](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:94>)

But a bad term \(K_{\rm bad}=c\,\eta\) is gradient-independent and is erased automatically by \(\eta\to0\), whatever \(c\) is.

Minimal fix: set profile jets to zero while retaining arbitrary constant zero-jet values—test both asymptotic uniform backgrounds \(W_\pm,\mu_\pm\), or constant \(w_1,m_1\) with live \(\eta\). Keep the \(\eta=\sigma_W=0\) reference test as an additional regression.

11. Two builder-facing target answers remain leaked

The spec supplies:

> “the leading amplitude coefficient … `(nonzero, the Born vertex content)`.”  
> — [S11c_d_SHARED_PHYSICS.md:303](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:303>)

That is neither guaranteed nor permitted as a target answer. The coefficient can vanish for \(k\!\cdot a=0\), a selection rule, or a form-factor zero. The spec also says:

> “Perturb the FORM … and require the mixing/leakage to move.”  
> — [S11c_d_SHARED_PHYSICS.md:416](</var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:416>)

M2 requires:

> “The spec says what to compute—never what anything equals, is expected, or was measured.”  
> — [CLAUDE.md:72](</var/projects/toy_physics/CLAUDE.md:72>)

Minimal fix: remove “nonzero” and make the strong-edge statement conditional on a nonzero computed coefficient. For the form ablation, emit baseline, altered-form operand, and residual; leave whether it moved to adjudication.

## Checks that passed

- The localized-interface versus bump distinction is physically correct once applied to the actual thickness profile and given an explicit one-dimensional transform: \(\int f'=\Delta f\).
- Keeping \(\eta\), \(\sigma_W\), and the wave-scale variable separately visible is correct. The imported operator is already first order in \(\sigma_W\); an additional WKB limit would discard the finite-width scattering content. Setting \(\eta=O(1)\) inside an operator truncated at first shape order is not controlled.
- Conditional on zero uniform amplitude, the amplitude/absolute-flux/fraction orders are \(O(\epsilon\eta)\), \(O(\epsilon^2\eta^2)\), and \(O(\eta^2)\); the \(\epsilon^2\) cancellation is correct and is not nonlinear-light physics.
- The strong-edge bridge is otherwise honest. For
  \[
  A_H=-i\epsilon\sin(\eta G),\qquad C=\sin^2(\eta G),
  \]
  \(C'(0)=0\), \(\tfrac12C''(0)=G^2\), yet \(C=0\) at \(\eta G=n\pi\). Thus a nonzero Born amplitude coefficient supplies no positive lower bound at finite contrast.
- The c2 carrier/source/\(\Phi\) debt, the `(0)−(0)` limitation, withdrawn F/G status, sign debts, and premise caveats are otherwise carried honestly.
- The tilt and advection corruption choices, fixed-anchoring rule, `RHO4_CONSTANT` structural absence, rejection of `∇W→0` as the N6 corruption, and prohibition on corrupting an anchoring are correct.
- Rejecting the generic 1D weak-well theorem for this multi-component, unequal-asymptote, possibly non-Hermitian interface is correct.
- N11’s rest-frame qualification and the chain/comparator/blind-Wolfram/supplied-versus-computed structure are preserved.
- No numeric lab bound or `O(1)` reductio is leaked.

Nits: none separate from the must-fix items above.

**Verdict: NOT-SOUND.** Outstanding findings: 1–11. I made no working-tree changes.
