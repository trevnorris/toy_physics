NEEDS REVISION

**Coverage.** I read `plan.txt`, `guide.txt`, `views/C3.json`, `receiving/outgoing-domain.json`, `receiving/grazing-plus-determinant-and-domains.json`, `views/force-A.json`, `views/force-B.json`, `views/receiving-chart.json` and `views/receiving-dual.json`. I did not open `grazing-minus`, the full five-field operator, the JL/G0/K1/B/DB views, the native files, the affine face files, the pressure assembly, or the step and path identity matrices. I also did not open the other ~1,800 artifacts. All algebra below was done by hand on the supplied C3 and force text. Nothing was executed, and I did not decode any hash or opaque object.

**Sound parts**

- **Evenness.** Every occurrence of `l` in the saved C3 is `l**2`. The reduction to s = p² − q² is therefore a syntactic fact, and the two-ray map is correct. The propagating ray is q = t on [0, p], and the evanescent ray is q = ir with r ≥ 0. Both grazing points are the shared endpoint q = 0.
- **Denominators.** All five distinct denominators are scalar multiples of (3+10i)q + 3i. Its only root is q = −3i/(3+10i) = −(30+9i)/109 = −β. That point has Re q < 0 and Im q < 0, so it lies on neither ray. The plan's "check the actual factors" is satisfied by this one factor, and the denominators can be shown nonzero on both rays without any root search. The chart and dual denominators l² + 1/20 and l²/20 + 1/400 are strictly positive.
- **Determinant step.** The ray determinant, the Re/Im gcd with Bezout, Sturm counts including the algebraic endpoint p, and the comparison with the saved threshold value form a legitimate finite-algebra test. The saved threshold determinant is −1341/10 + 184758i/625, which is nonzero. A root on a ray would be a real obstruction and should stop the route.
- **Growth bound.** Rows 2 and 3 of C3 are rational in q with numerators of degree 3 and denominators of degree 1. Entries therefore grow like r² on the evanescent ray, and adj/det is polynomially bounded if det is not identically zero. A degree drop in det only changes the polynomial exponent.
- **Cross-flux.** The C₋₊ and C₊₋ contractions are the right check. In a non-dissipative homogeneous medium a conserved flux forces the e^{±2ipx} cross term to vanish. With memory dissipation it may not, so a nonzero value is a physical result. It must stop the additive-weights claim, as the plan says.

**Blocking findings: the source-class and decay argument fails on the live operands**

1. **The step is not removed by the Q-identity.** `Q·Ahat = −i·(A')hat` applies to the combination l·F₀ + F₁. That is the `receivingDotDiagnostic`, and the artifact marks it as not a loss projection. The receiving rows are the saved dual applied to the force, and the dual has no (l−p) factor.
   - Column A, etaForce = (iA/5)(1, −p, 0, 0, 0). Dual row 0 gives −2p·(iA/5). Dual rows 1 and 2 are proportional to (1+4lp)/(l²+1/20) and (l−p/5)/(l²+1/20), times iA/5. Row 1 is nonzero at l = p. Row 2 is (4p/5)/(p²+1/20) there, which is also nonzero.
   - `localRight` is (9i/5)(1, −p, 0, 0, 0). So the dual-projected force does not decay at the right end.
   - Consequence: Ahat contributes both δ(l−p) and a principal-value 1/(l−p) term to each receiving row. These sit at l = p, which is q = 0, the grazing point. This is not a rare corner: p² = 119/20 is also the bulk constant in q² = 119/20 − l². The incident momentum and the grazing threshold coincide numerically, so the step singularity and the √ cusp of C⁻¹ occur at the same point.
   - The plan's claims that the weighted Fourier fields are "L1 and L2" and that Riemann–Lebesgue decay holds as x → either end are false for the step part. 1/Q is neither L1 nor L2. The cusp term Q^{-1/2} is L1 but not L2.
   - The correct picture is a non-decaying driven plane wave at the right, about e^{ipx}·C(p,0)⁻¹·(dual·F_right). The slow remainder is a √-cusp tail, roughly x^{-1/2}.
   - The route needs this split made explicitly: asymptotic right-end particular solution plus a remainder driven by A′ and the other localized profiles. Only the remainder can be put in the decay class.

2. **The U_B zero-forcing claim is contradicted by the visible operands.** The B etaForce is (0, iA/10, −iA/5, 0, 0). It is proportional to chart column 0, (0, 1/10, −1/5, 0, 0). Dual row 0 = [0, 2, −4, 0, 0] gives iA, which is nonzero. Dual rows 1 and 2 give zero. `receivingDotDiagnostic` is 0, but that is the k·F contraction, not the dual rows.
   - Unless the force index convention or the pairing of the dual differs from what I read, B drives receiving field 0 with a full step.
   - So "its receiving scalar/longitudinal response is zero" is unsupported until the plan states which field chart column 0 is. The plan calls the three fields longitudinal displacement, theta and e_W, and I could not match column 0 to any of them from the supplied text. The sigma and pressure additions to B also have not been shown to cancel this term.
   - The plan lists this check as to-do. It should be moved ahead of any source-class argument, because it decides whether the B column has a receiving response at all.

3. **The missing physical assumption is that the perturbation is localized.** η switches on in the right half-slab (left end 0, right end 9). The first-order response there is a non-decaying, grazing-coupled field whose radiated flux contains ∫ Re q·|ψ|² dl with |ψ|² ~ 1/Q². That integral diverges at q → 0. Equivalently, a perturbed uniform right half-slab leaks at a constant rate over an infinite length. So an unnormalized "total transverse deficit" from this step is infinite or ill-defined.
   - Either the observable must be a per-length or amplitude-renormalized quantity, or the η step needs a physical cutoff that the plan does not supply. This is a physical choice and should go to the user.
   - Only the sigma and direct-chemical profiles are localized and Schwartz. The η-step force is not.

**Other items**

- Uniqueness "a.e. in weighted L2" is vacuous for the step part, since that part is outside the class. It also does not exclude non-decaying homogeneous solutions at the threshold.
- The plan should state that the K1 = 0 and G0 end-flux results used a matched non-decaying incident-plus-step amplitude, and say how the right-end particular wave in finding 1 enters T1 and the matching.
- The plan should confirm from the original operands that the bulk wavenumber in q² is p. Here ω/cs = √6 ≈ 2.4495, not p ≈ 2.4393. A mismatch in the grazing identification would change everything at q = 0.
- The "direct pressure once" and "affine face" requirements are plausible, but I did not inspect the pressure assembly or the face files. I make no finding on them.

**Smallest sound next step.** Keep the C-block ray-determinant work and the cross-contraction work as written. Both are small and can legitimately stop the route. Before any adjugate or source-class claim, add these two items:
1. A hand-verifiable dual-row projection of both force columns, including the right-end asymptote and the B column, with chart column 0 identified.
2. A step-plus-remainder decomposition with explicit δ and principal-value treatment at l = p, and a decision from the user on how the infinite-length leakage is normalized.

I do not certify a regularity or source-class result, leakage, or a power identity. No power identity is required before the determinant test, because the determinant test is an independent algebraic fact. It is not enough on its own to support the stated decay conclusions.