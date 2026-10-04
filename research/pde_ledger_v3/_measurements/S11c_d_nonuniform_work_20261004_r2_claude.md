NEEDS REVISION FOR THIS SOURCE-BOUND NONUNIFORM WORK CONSISTENCY METHOD

Most of the method is sound. Two concrete gaps remain that the finite plan does not close, plus a few smaller required edits. I did not execute anything, and I read only part of the packet (coverage below).

## What passes as stated
- **Kinetic input.** `native-kinetic-definition.txt` has T = ε²[ρbr·|u_t|² + μ_W·(∂_e(W_bg·e_local)·e_t)²]/2. `native-local-thickness-definition.txt` gives e_local = W0·eW/W_bg. Together these give a thickness rate of W0·e_t. The plan (§2) states this correctly, with constant μ_W and no W_bg² inertia. This matches S11c_b_SHARED_PHYSICS.md §1c, where T = ½ρ_br⁰|u_t|² + ½μ_W(∂_tδW)².
- **Density.** `native-density-definition.txt` gives RHO4_CONSTANT as rho4 = ρ_br/W0 and Σ0 = ρ_br·W_bg/W0, with derivatives. The plan's operand join matches it.
- **Chemical frame.** `mu_specific_raw` and `affinity_raw` bind μθ/rhobr_bg_exact before the shape projection, and the affinity subtracts pressure/ρ_m. The plan (§2) correctly keeps this distinct from constant ρ_br. It also keeps the projected reciprocal separate from physical inverse density, so epsilon enters once.
- **Mass and constraint records.** `native-virtual-mass-definition.txt` and `native-mass-evolution-definition.txt` are independent definitions. The material route's density-gradient/Jacobian terms and the Eulerian `BACKGROUND_ADVECTION` term are not in the homogeneous constraint. §3 correctly keeps the physical density rate independent. The plan's reading of the chemical/mass-rate work as the difference from the homogeneous virtual variation matches `operator_from_density`, which keeps `ADVECTIVE_MASS_OPERAND` separate and puts the θ row in the mass-evolution slot.
- **Face-memory algebra.** I checked these by hand:
  - χj = j²/Λ + d_t[τj²/(2Λ)] follows from τj_t + j = Λχ.
  - The phasor mean Re(A_mem)|χ|²/2 with A = Λ/(1−iωτ) is consistent.
  - The split (P+Xχ)V + μj = Pv + χj + (Xχ)V follows from μ = χ + P/ρ and v = V + j/ρ.
  - The acoustic constant 1/(4πρ_mω) follows from P = −ρφ_t and a mean of Re(P*v)/2.
- **Separation of objects.** The split between P_rect, formal Q_port and the actual survival asymptotic is explicit. Reflected transverse flux is kept as survival. The threshold remainder and G2/end-current debts are named as deferred.

## Blocker 1: the background support is a free symbol, and the plan has no step for it
- In `S11c_b_brane_operator_sympy_audit.py`, `admissibility_support` (≈l.4081) is built from `f_hold_u_*_0`, `f_hold_theta_0`, `f_hold_e_W_0` and `t_hold_*` (l.339–345). These are `PREMISE` symbols, not a supplied law.
- `S11c_b_SHARED_PHYSICS.md` §2b says the residual (operator operand − support) is "not asserted zero". It also says a nonzero residual is an admissible outcome. The registry carries `f_hold_theta_0` symbolically (7 hits in `S11c_b_exports.py`).
- `S11CB_ADMISSIBILITY_RESIDUAL` appears nowhere in plan.txt or guide.txt. Plan §3 and §4 treat only an "independently supplied support law" as external power.
- No such law exists, so a nonzero background residual is a named, unclassifiable work term. It multiplies the field, so it is linear in the second-order fields and sits in the ε² work balance.
- Required change:
  - Add the ADMISSIBILITY_OPERATOR_OPERAND, SUPPORT_OPERAND and RESIDUAL registry entries for LAB_HELD/RHO4 as an explicit work-classification class.
  - State the retained-grade condition under which the residual is zero. Otherwise carry it as an UNSUPPLIED obstruction, in the same ledger as 20/02, with its coefficient on the second-order fields.
  - Do not take "held" support as zero by name.

## Blocker 2: the face-centre degree of freedom ζ_c is absent from the work ledger
- The source treats ζ_c as an independent face DOF with velocity `zeta_c_t`. Its generalized row, `CENTER_FACE_GENERALIZED_ROW`, comes from `face_generalized_force_rows` (`native-face-work-row-definition.txt`). SHARED_PHYSICS §1a forbids setting ζ_c = 0.
- plan.txt never mentions ζ_c. §3 pairs "five rows" (u×3, θ, e_W). §5 uses V_s = W0·eW_t/2 for LAB_HELD.
- T has no ζ_c inertia, so net face power through ζ_c_t must be balanced by the centre row or appear as a work mismatch. A one-sided or both-face pairing that omits this row can be wrong by that term.
- Required change:
  - Add ζ_c_t × the centre row to the §3 pairing and classification.
  - Say whether ζ_c is slaved or constrained by the saved system in the LAB_HELD/RHO4 specialization. If so, give the source operand that does this.
  - Check that V_s on both faces is consistent with that.

## Smaller required edits (finite checks)
1. **Two orderings.** §4 expands in the defect parameter λ (f0 + λf1 + λ²f2). The work identity also has the amplitude ordering ε². State that time-averaged work uses only the linear-in-ε fields, and that rectified/second-harmonic ε² fields are excluded only because Blocker 1's residual vanishes or is carried. Today this rests on the linearised rows.
2. **Truncation consistency.** Kinetic rows use `rhobr` exactly via `density_pair` (the original is RHO4: ρ_br·W_bg/W0). The stored rows are truncated at `STRONG_ROW_JET_DEPTH`. State that any mismatch between the exact kinetic grade and the truncated stored grade is a named order defect, not absorbed into the residual.
3. **Shape projection.** `mu_specific_raw` multiplies by `source.parameter` and applies `shape(...)`. The plan should name that parameter and projection as part of the normalization join, so controls can move it.
4. **Control applicability.** §7's controls need an applicable nonzero operand each. Add the ζ_c/centre-row control and a support-residual control. An inapplicable control counts as untested, as the plan already says.

## Not blockers (deferred as the guide states)
- The finite-defect survival asymptotic, the threshold remainder and the memory-norm remainder.
- G2/second-order end currents.
- Evaluation of the face integrals.
- The face-map jet-order and L² list in §5. Listing it as a duty is adequate for a method.

## Optional suggestions
- Present the classification as a table keyed by origin: stored, kinetic, mass-evolution, face virtual work, FACE_FLUX, support residual, ζ_c.
- Record the five LAB_HELD/RHO4 term-origin registry keys by name so the worker can cite them.

## Coverage and uncertainty
- **Read fully:** guide.txt, plan.txt, review-prompt.md, and all nine `views/native-*-definition.txt` files.
- **Read in part:**
  - S11c_b_SHARED_PHYSICS.md: §§1c, 2a, 2b, 3a–3d, and the opening/scope lines.
  - S11c_b audit: `admissibility_support`, `task_admissibility`, and the ζ_c bindings by grep.
- **Not read:**
  - The packet index, pin-join views, local-cell tables and native receiving/flux evidence.
  - The S11c-a physics and audit beyond the supplied view definitions.
  - The large exports (grep hits only, no decoding).
- **Uncertainty:** I did not decode the registry value of the admissibility residual. I only saw that the support symbols are free. Blocker 2 rests on source-text reading, and it could be answered if the saved system already slaves ζ_c. Neither point was executed or numerically verified.