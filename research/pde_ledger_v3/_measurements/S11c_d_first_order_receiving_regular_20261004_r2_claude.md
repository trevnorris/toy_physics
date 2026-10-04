CLEAR FOR THIS REAL-AXIS RECEIVING AND CROSS-FLUX METHOD

I found no unsound step in the method as written. The clearance is conditional on the requirements and coverage limits below. All checks are hand arithmetic on the supplied views. I ran nothing and decoded nothing opaque.

## What I checked by hand

**Chart and dual (section 0).**
- In `receiving-chart.json` the columns are t, b(l), k(l), e_θ, e_W, as the plan labels them.
- In `receiving-dual.json`, row 2 is (l, 1/5, 1/10)/K2, matching R3 row 0, with K2 = l²+1/20.
- I multiplied out D5·E5 = I for the first three rows against the three chart columns, and it holds.
- `complete-block-matrix.json` shows exact zero off-diagonal blocks between {0,1} and {2,3,4}, so the two-way decoupling is real.
- The transverse symbol is (3/2)(l²−p²), so the transverse poles are only at l=±p.

**Incident columns and the projection table.**
- Factoring i/5 out of U_A gives (1, −5p), as stated.
- I recomputed all five D5 rows on both columns. Every entry of the hand table matches, and rows 3 and 4 vanish.
- U_A = −2ip·t + 4i·b(p), and U_B = i·t. The K2 = 6 reduction at l=p reproduces the 4i coefficient.
- `eta-A` and `eta-B` equal A(x)U with A = 5/2 + 9T/2 + 2T².
- The right-end forces equal 9U, and the left-end forces are 0.
- The identity Q·Ahat = −i(A')hat and A' = (9/2+4T)(1−T²)/L both check, as does ∫A' = 9.

**Section 3a gate.**
- `B.json` gives B_A(l) = −2ip·t + 4i·b(l) and B_B = i·t exactly. Both lie in the transverse chart span for every l.
- So R3(l)B(l) ≡ 0, and the gate d/dl[R3·B] = 0 at l=p is satisfied analytically.
- Numerically, R3(p)B'(p) = −i/30 (using `Bprime.json`) and R3'(p)U_A = +i/30. They cancel.
- The plan is right that a frozen R3(p)B' gives −i/30·δ_p ≠ 0. This is the longitudinal admixture from the rotating polarization, and it is cancelled by the δ' term.
- The sign of the δ' jet checks: the transform of x·e^{ipx} is 2π·i·δ'.
- This gate is an identity of the chart construction, not an independent physical test. The worker should still restore it with full maps. It is cheap, and the dropped-term control is meaningful.

**Source class.**
- U_B: the saved `receivingDotDiagnostic` is identically 0 for eta plus sigma/10. Rows 3 and 4 are zero for sigma and eta. The pressure `source01Envelope` is 0. The three-row forcing is therefore zero by direct reading, not only by a face flag.
- U_A:
  - The sigma rows 3 and 4 vanish at T=±1.
  - The `source01Envelope` polynomial vanishes at T=±1.
  - Row 2 of sigma-A vanishes at T=±1.
  - Together these are consistent with Schwartz localization.
  - I did not verify sigma rows 0 and 1 at the ends.
  - `localRight` and `localLeft` list only the 9U eta step, which is consistent with this.

**C3 block.**
- Every entry depends on l only through l², so evenness is evident.
- Every denominator in `C3.json` is a constant times (q + (30+9i)/109). No other factor appears, and its real part is positive on both rays.
- All coefficients are Gaussian rational, since p² = 119/20. The exact gcd/Sturm plan is therefore feasible, with the interval endpoint p algebraic.
- Because s = p²−q², det C is a function of q alone. The two rays are exactly the real segment q∈[0,p] and the imaginary half-line.

**Cross current.**
- I contracted the upper 3×3 of `JL.json` at kl = −kr. I used B(−p) from the saved B(l) formula, which is transverse to k(−p). I took all four combinations of U_A and U_B. Every entry is 0, with exact cancellation between the −7/50 and −7/100 terms.
- J(p,p) contracted with U reproduces G0[1,1] = 9p/80, so the convention joins.
- This is only a hand preview. The plan still requires the actual contraction, sign-mutation sensitivity, and the phase record. The result is unexecuted evidence, not a certification.
- If the contraction is zero, the end flux has no interference term. That does not by itself make the survival deficit a sum of independent weights.

## Physics points the worker must carry

1. **Degenerate point.** The slab transverse wave speed squared is 3/2, which equals cs². So the incident transverse wavenumber l=p sits exactly at the exterior grazing branch point q=0. The forced source content at Q=0 (value 9/(5K2)) meets the square-root cusp of C3. This is the real reason finite threshold determinants matter. Even if C3⁻¹ is finite there, the response has cusp-type, algebraic x-decay. L¹/L² and Riemann–Lebesgue still hold. Nothing stronger should be claimed.
2. **Uniqueness.** It rests on the outgoing branch selection and the weighted-L² class. It does not extend to distributions supported at l=±p, where C is not smooth. The plan states this correctly.
3. **Large-l behavior.** The entries grow like l², and the l² coefficient matrix may be singular at q→∞. The degree and leading-coefficient records are therefore essential. Polynomial growth follows only because det is rational in r and not identically zero.
4. **Chart factor.** E3(l) carries k(l), which grows like l. The weight in the L¹ claim must include it.

## Why this is the smallest sound next step, and what it will not give

- Regularity of C3 and the U_A/U_B source split are prerequisites for any finite radiated or dissipated power. No power identity needs to precede them, and I see no concrete obstruction to this order.
- A full pass still gives only the first-order outgoing three-field response to U_A. The O(λ²) deficit needs either second-order transverse amplitudes or a power identity that converts the second-order end flux into first-order quantities. K1=0 only kills the O(λ) flux.
- The power identity is therefore the next derivation, not an optional one. It has to include:
  - the exterior acoustic flux,
  - memory dissipation,
  - the mass-rate and chemical work,
  - and the LAB_HELD prescription work.
- A static modulation does no net work only if the modulated operator is energy-conserving apart from the memory part. That is an assumption to be derived, not asserted.
- The nondecaying secular term δ_p·x·U in the transverse sector means first-order perturbation theory is nonuniform in x. The plan keeps it, and any later second-order step must as well.

## Coverage limits

I read these in full:
- `plan.txt`, `guide.txt`, `review-prompt.md`
- the incident-columns, eta, sigma, right-force, chart, dual, block-matrix, C3, B, Bprime, JL, G0 and K1 views
- `force-A.json` and `force-B.json`

I read only the first 945 of 4401 lines of `matching/pressure-transverse-absence.json`. These I did not open:
- `receiving/pressure-receiving-assembly.json` (too large)
- the affine-face-column files
- the outgoing-domain and grazing-determinant records
- `source-contracts/`, `frozen-engine.py`, `original/native/LEFT-original.json`, and the work/energy maps
- the 781 source and 725 receiving artifact families

I therefore have not verified that direct chemical pressure enters the affine faces exactly once. I also have not verified the sign, scaling and outgoing-branch conventions against the originals, or that the force rows are in the same equation basis as the left-multiplication by D5.

The `-p` incident column has no saved first-order forcing. This is not needed for the left-incident survival question, but it blocks any claim about a right-incident or reflected-wave source.

The determinant, adjugate, root and cross-contraction results are all unexecuted. Any root, nonzero residue or nonzero cross term must stop the route.