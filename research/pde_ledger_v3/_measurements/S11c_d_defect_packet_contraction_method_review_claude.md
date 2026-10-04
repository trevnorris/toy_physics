CLEAR FOR THIS BOUNDED FINITE-WINDOW CONTRACTION METHOD

The clearance covers the algebra and domain reordering of J and the three ADDED direct terms only. It is not a numerical, build, author or result clearance. I found no fatal method error. I did not open the saved templates, the 544 addresses, the unit records or the geometry files, so everything that depends on them is listed below as a build obligation.

## Verified by hand against the source

I checked `kernel_components` (`S11c_d_defect_packet_inner_lib.py:53-60`) and `middle`, `q`, `profile` and `exact_cuts` in the same file. I also checked the section 2-3 statements of `input/packet-action-method.md`.

- **Factors:**
  - `profile(t)` is `5t/(2 sinh 5πt)`, which equals A(z) with L=10, and A(0)=1/(2π).
  - `q` is `sqrt(5.95-p²)` on the real branch and `i·sqrt(p²-5.95)` outside. That matches the first-quadrant outgoing branch with κ²=595/100.
  - J's denominator is exactly E(l)E(k+t)E(k)·q(k+t)·(q(k+t)+q(k)).
  - `pref` equals Cd·A1·A2/(E(k)E(l)), with Cd = (-iμ)WL/(4i).
  - The reference lines for J and for Cd/(E·E) hold, so no extra resolvent is inserted.
- **J, height and quadratic (m=k+t):**
  - (k+t)(2k+t) = m(m+k). The method's BJ therefore follows, with m·C0T giving the m·k⁰ part.
  - k(t+2k) = k(m+k), so BD_height = C2T + m·C1T, over q(m).
  - The quadratic term is q(k)²/q(m), which gives X02T.
  - A(t)A(l-k-t) = A(m-k)A(l-m).
  - The window |m-k| ≤ T falls on k only. The l range stays [-K,K], and the profile A(l-m) is unwindowed in the source.
- **Reflected (m=l-t):**
  - 2l-t = l+m.
  - A(t)A(l-k-t) = A(l-m)A(m-k).
  - The k integral is full, so it gives X1(m). The window |l-m| ≤ T falls on l, so it gives I_T(m).
  - The orientation argument is right: dt = -dm, and reversing the limits gives a positive measure. I_T(m) is symmetric in the sense needed, so the two windows act on different variables, as the method says.
- **Wrong-root mutant:**
  - With q(l-t) replaced by q(k+t) and m=k+t, 2l-t = 2l-m+k.
  - The numerator is k(2l+k-m), which gives 2·X1T·Y1C + (X2T - m·X1T)·Y0C. The Y factors are unclipped and the window sits on k.
  - This matches section 5 exactly. It is correctly a different substitution from the correct reflected one.
- **Outer range and wings:**
  - I_T(m) is non-empty exactly for |m| ≤ K+T. It is the full [-K,K] for |m| ≤ T-K (95 for both (27,122) and (29,124)).
  - Clipping occurs on 95 < |m| ≤ K+T (149, or 153 for the enlargement).
  - Both ±(T-K) kinks and the outer endpoints must be cuts. The method lists ±K±T.
- **Bounds:**
  - For q in the closed first quadrant, |q1+q2| ≥ (|q1|+|q2|)/√2. Since Re β>0 and Im β>0, |q+β| ≥ |β|, which is stronger than the stated bound.
  - The ratios q(k)/(q(m)+q(k)) are bounded by √2, so the J and height majorants are C/|q(m)|.
  - The reflected 1/(q(m)+q(l)) is locally integrable even where both depths vanish, since 1/(√|a|+√|b|) is integrable in 2D.
  - The original 3-D integrand has only these integrable singularities, so absolute Fubini on the finite domain is sound for these ordinary terms.
- **Scope statements:** the PV/H, flat, contact and slope exclusions and the "no unsubtracted-PV interchange" disclaimer are consistent with `packet-action-method.md` lines 107-135.

## Build obligations, to be met before or inside the certificate build

1. **Separability join.** The reordering holds only if each address's integrand is exactly Y(l)·N_f(l)·[kernel]·X(k). Any additional (l,k)-dependent coefficient in the saved Jwhole/Dwhole templates would break it. The runtime join must check this per template, with the normal sign at l (+iq, -iq, or 1 for pressure) as stated. I did not read the 20 templates.
2. **Persistence.** Persist the old and new full expressions and their residuals for all four primitives separately, including the polynomial identities above. The domain certificate must record forward and inverse maps, orientation, windows and the transported collision labels.
3. **Collision transport.** The old surfaces k+l=0, k+l=±2κ and l=k come from the two-root and profile coincidences, and they need explicit per-primitive treatment. In the m-form the J, height and quadratic terms have no q(l-t) route, and m=k and m=l are separate profile points. The reflected term has no q(k+t) route. Each surface must be shown absent or retained, not dropped.
4. **Mutant control witnesses.** The wrong-root mutant needs its own original-mutant integrand and domain join, with symbolic witnesses chosen before results. The normal-slot mutation needs a non-flat applicable address.

## Remaining numerical-readiness limits (not approvals)

- **Inner regularity.** The inner functions CrT(m) and YrCT(m) are bounded but not analytic at m=±κ. The ratio of square roots near the double zero may carry non-analytic terms beyond a pure √ behavior. Square-substitution cuts are necessary but may not give the assumed convergence, and this is unestablished. The inner integrals also have a narrow-kernel scale from A, whose poles are at ±0.2i. The A-pole scale should therefore appear in the resolution cuts.
- **Error propagation.** The inner error multiplies by Y0(m)/q(m) and the outer weights, with an inverse-square-root outer singularity. The allocation scheme, request reuse and storage/cost are unapproved, as the method states. I give no runtime forecast.
- **Sharing.** The contractions separate into X-only (C, X1T, X02T) and Y-only (Y0, Y0C, Y1C) pieces. Sharing across addresses is allowed only under the exact-match rule the method states. Preserve all 544 addresses and the zeros regardless.
- **Fubini scope.** The compact-parameter Fubini argument is not uniform in K and T. The tail obligations are unchanged.

## Optional suggestion

State the central window |m| ≤ T-K and the two wing ranges explicitly in the certificate, as a check on the clipped-wing control.

## Coverage

**Read:** `input/method.md`, `input/evidence-guide.md`, `original-source/S11c_d_defect_packet_inner_lib.py` (lines 25-140), and `packet-action-method.md` (lines 60-180). I also saw the file listing, which was truncated at 100 of 1548 entries.

**Not read:**
- The saved templates, adapters, field/unit/geometry JSON and tail files.
- `source-contracts.json`.
- `pressure-readiness-method.md`.
- The other original-source files, including the symbolic `kernel_components` call at `S11c_d_defect_packet_inner.py:139`.

I gave no numerical value. The provenance, unit, 544-address and collision-label claims remain exactly as uncertain as the inherited evidence leaves them.