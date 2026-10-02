Both verdicts are positive. I found no recipe fault, only the uncomputed objects the method says it will build and test. Everything below is hand source and algebra checking from the staged input only. I executed nothing, and no status flag or zero guard was treated as acceptance.

## Hand checks by item

**1. Ordered transfer and normalization: holds as stated.**
- **Profile hats:** native `shape_source` (`native-c1.py:591-599`) gives h_hat = η·W0·ŵ/2 and tilt_hat = σ·jet_hat/2. I re-derived the transforms from ∫tanh·e^{-itx} = -iπ/sinh(πt/2) and ∫sech²·e^{-iτx} = πτ/sinh(πτ/2). Both give ŵ = δ/2 + L/(4i·sinh(πLt/2)) and jet_hat = L·A, matching `profile-transform-return.json`.
- **Measure:** `profile_bindings` (`native-c2.py:850-861`) uses a (2π)^-3 forward transform. A function of y alone therefore gives exactly δ(Δ_edge)·δ(Δ_edge)·f̂_1D, with no extra 2π.
- **Middle integral:** `kernel_apply` p2 (`:469-470`) uses the same source/(2π)^3 as p1 and a plain d³m. The product of two normalized hats integrated over m is a plain convolution, and the two edge delta sets collapse to one. The Jacobian is +1.
- **Units:** reduced to 1D, the direct slot ∫dt·C·ĥ·ĵ has dimension M L⁻¹T⁻¹. That equals z01 and a·z01·z12·dm, both counted separately.
- **Saved k=0 integrand:** the method's raw integrand at k=0, Q=1/10 equals the saved integrand exactly, including sign and the 25√3/8 factor. I computed this by hand, not by executing anything. I could not find a missing factor or term.

**2. Contact and PV: holds on its stated domain.**
- **Generic coefficient:** I re-solved the saved boundary equations from `boundary-before-guards.json` and reproduced the mixed coefficient and the saved `A11`. The velocity and pressure expressions match the Taylor expansion of n·∇φ at z=h. The saved `mixed` is therefore the full ησ coefficient, not a k=0 selection.
- **Factorization:** C = H·B follows from N = -H·qi² + k·qi·(qh-qi) + k·qh·(qo-qs) and the two squared-depth identities. I get exactly the stated B+, including the signs of the second and third terms.
- **Contact:** C(k;0,Q) = 0 when qi, qo ≠ 0, because qh=qi and qs=qo at H=0. The contact is zero, and t·PV(1/t)=1 gives the stated ordinary integrand.
- **Domain wording:** "regular external depths" should read qi≠0 and qo≠0. Only then is B analytic at t=0 and the sums qh+qi=2qi and qs+qo=2qo nonzero there. On the first-quadrant sheet, qh+qi and qs+qo vanish elsewhere only when both terms vanish, meaning coincident branch points. That is a pole exclusion to list, not a flaw.
- **Not done:** no distribution identity is inferred at qi=0 or qo=0, and the method correctly refuses to.

**3. Sheet and face signs: holds.**
- **Extension:** it matches `native-c2.py:505` and `:582`, with f·q on the same q.
- **Normal:** `native-geometry.py:844-850` has h0=f·W_bg/2, so the lab slope is f·∇W/2. The normal is then (-f·∇h_lab, f)/√, whose tangential part is -∇W/2 for both faces.
- **Lower face:** in outward variables the boundary equations are literally the upper ones, with e^{iqh} appearing for both faces. The derived lower C is therefore the same N. The jet slot carries i·f·qo, and the minus consumer coefficient has the opposite sign (`consumer-unit-joins.json`), so f cancels between them.
- **Native linear check:** native `shape_source` and `dtn_first_kernel` are face-independent in outward variables.
- **Rigid-tilt test:** at first order, the slope coefficient kωρ/q² cancels the H→0 derivative of the height coefficient. That supports the slope sign convention with q unflipped.
- **Mirror:** I did not use reciprocity as a check. The independent h and s slots are not a surface, so it does not apply to the mixed coefficient alone.

**4. Closure and consumers: holds.**
- **Direct term:** with M = I + a·Z3, (M⁻¹Z)_{02} = (z02 − a·z01·z12/M11)/(M00·M22). The direct term therefore carries R(kout)·C·R(kin) once, and the first-shape iteration is unchanged.
- **Zero grade:** z0 = ρω/q matches the saved diagonals √3 and 3/2.
- **Saved referenceMixed:** at the saved point it equals slot·R(kout)·R(kin) with a trace diagonal of 1. `jetMixed` is i·qo times that.
- **Grazing:** explicit 1/q cancellation is necessary but insufficient, as the method says. N, qh, the trace inverse and i·f·qo remain.
- **Consumers:** the zero-grade jet consumers are zero (`zero-grade-jet-consumers.json`), so the new (1,1) jet slot contributes nothing in the rectangle. The pressure slot contributes.
- **trace_three[0,2]=0:** this rests on the saved trace being height-only. The saved evidence is plus-face only, so the minus join is still owed.

**5. Applicability: holds.**
- **Branch locations:** t = -k±Κ and t = k+Q±Κ are correct. qs is a cusp in the numerator. 1/qh is an integrable inverse-square-root singularity.
- **t_b formula:** |t_b| = |qi²/(Κ+k)| is exact. The ≈ qi²/(2|k|) form holds only for fixed nonzero k near |k|.
- **Tails:** the integrand decays like A(t)·ĵ, so tails are fine.
- **Left as acceptance work:** the closed two-sided limit, the integrable bound, and the exchange with the middle integral.

## Non-blocking conditions for implementation

1. **Slot hats and depths.** `fourier_profiles` (`native-c2.py:429-435`) puts one (a−b) argument on every hat. The slot needs ĥ(m−kin) and ĵ(kout−m) with distinct arguments. `outgoing_spectral` (`:439-452`) has no symbol for qs = q(k+kout−m). It also selects its root by `roots[-1]`, so the explicit rule needs an equality check against the native values at the saved points.
2. **Selected-point scale.** At the saved k=0 selection, input grazing puts both branch points at ±Κ = ±qi, which is linear in qi, not qi². Use exact roots for the partition, as the method already says. The output-side qs branch has the analogous scale qo²/(Κ+k+Q).
3. **Closure-input insertion.** The saved `closure-input.json` insertion (a scalar at [0,2]) is a whole-D factor-join control. It cannot be the slot representation, since `kernel_apply` already integrates m.

## Verdicts

**ORDERED CONTACT/SHEET/FACE RULES: SUPPORTED WITH STATED DOMAINS**
Domains: qi≠0 and qo≠0 for the contact and factorization, nongrazing for the raw output, and the positive root on the shared sheet.

**NEXT INCREMENT METHOD: CLEAR FOR BOUNDED IMPLEMENTATION**
No finite matrix changes, and no exact-match claim.

## Next action

Build the symbolic task-local object under the existing guard. Persist these before any check:
- the general-k C and C=H·B residuals;
- the lower-face boundary residuals;
- the times-two, wrong-sheet, wrong-face and slope-omission controls;
- the saved k=0 integrand reproduction and the R-factor reproduction, each with its exact residual.

Then do the closed two-sided limit and local bound. If an exact-match piece lacks a supported integration prescription, record the obstruction and stop.