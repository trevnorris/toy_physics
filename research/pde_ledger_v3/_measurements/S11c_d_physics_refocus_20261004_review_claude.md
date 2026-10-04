I've read the packet and the key native lines. I'm writing up the verdict now, with one anomaly in the LEFT uniform record that I couldn't interpret.

# Verdict: CLEAR FOR THIS PHYSICS REFOCUS AND FIRST-ORDER SOURCE PLAN

This covers only the priority and bounded source-plan preparation. It does not cover a worker, a receiving-field calculation, a grazing expansion, a loss result, calibration or a leakage number. The clearance depends on the conditions below.

## What I read
- **Read fully:** `proposal.txt`, `guide.txt`, `review-prompt.md`, `packet-index.json`, `refocus-inspection.json`, the amendment draft, and the loss assessment.
- **Read in part:**
  - `SHARED_PHYSICS.md` §1c and §3a–3d (lines 183–302 and 415–680).
  - Native source lines 3040–3100 and 3290–3424.
  - Address 8346, input header and return header (grades, row, field).
  - `point-01-LEFT.json`: header, current Gram and eigenvalues, one residual block and the bulk-current block.
  - `defect_packet_action_method.md`, first 60 lines.
  - `near_unity_uniform_plan.md`, by grep only.
- **Not read:**
  - The other 19 address files; I relied on the inspection census for them, which was not independently re-checked.
  - `point-01-RIGHT.json`, `original-numeric-context.json`, the `evidence/physical/` files and the contracted numeric method.
  - Most of the native source (about 545 KB).

## 1. Target
- **L_T is the right target.** It matches `SHARED_PHYSICS.md` §3b (`P_T,surv`) and the amendment §1. It keeps all reflected and transmitted transverse power as surviving. It excludes forward extinction, decay rates and bound capture.
- **Double counting is correctly guarded.** The amendment (§2) says to count non-transverse boundary flux plus dissipation inside one control volume. Survival deficit is not a third mechanism.
- **Weak coefficient vs strong-edge number.** A fixed-input weak coefficient is distinct from the unresolved strong-edge number. §3d states this, and its `sin²(ηG)` counterexample shows a Born coefficient gives no lower bound at η=1. The proposal respects that (lines 16–24, 103–105).

## 2. Order argument
- **The logic is sound and conditional.**
  - **Leading term:** if every loss amplitude has a0=0 and the power form is regular, then P_loss = λ²B0[a1,a1] with a1 = C0Ψ1 + C1Ψ0. That needs L10, L01 and the direct maps, not Ψ2 or L11.
  - **Why grade11 can be parked:** along η=λ, σ_W=λ/10, the grade-11 operator is λ². It enters only Ψ2 and the a0–a2 interference, which vanishes when a0=0. I checked this against §3c's J^(2).
  - **Matches the inspection:** all 20 addresses are grade 11 (I confirmed this directly for 8346), with source and consumer grade 00 and input e_W. The packets are not transverse modes (action method §1).
- **Parking is justified only conditionally.** It is not a global discard. It holds if the conditions below hold, and the proposal says so (lines 44–50).
- **No wave nonlinearity is involved.** The grades are bookkeeping for a linear vertex, consistent with §3c's last sentence.
- **Where the inference fails:**
  1. **a0 ≠ 0:** if any loss channel has a nonzero reference amplitude or drive, B0[a0,a2] enters at λ². Net-zero baseline power is not enough. The amendment §2 states this correctly.
  2. **Vanishing first-order projected source:** L1Ψ0 and C1Ψ0 could be nonzero as entries yet project to zero on the loss sector. That only shows the leading λ² loss is absent. Then L11, L1Ψ1 and the preserved mixed work return as the next order, so grade-11 work is not obsolete.
  3. **Grazing:** I flag this as my own analysis, not something the packet shows.
     - LEFT has depth exactly 0 (`point-01-LEFT.json`, stratum EXACT_GRAZING). The bulk depth root therefore has a branch point exactly at the incident momentum.
     - The step's 1/(Q_n ± i0) pole also sits at Q_n=0, which is the same place.
     - So the loss-projected source must vanish at Q_n=0 fast enough to cancel that pole, against the square-root density of states.
     - If it does not, B0[a1,a1] diverges and the regular expansion fails. Higher orders then cannot be dismissed.
     - The saved finite uniform fields do not bound this.
  4. **Truncation:** the retained operator drops O(η²,σ²) terms. This does not matter at λ² if a0=0, but it does matter for any transverse-sector quantity and for the pole-binding statement in §3b.

## 3. Smallest next ingredient
- **First: the actual projected first-order forcing and face sources.** These are L1Ψ0 and C1Ψ0 on the LEFT incident doublet, including end matching.
- **Why it comes before any inverse.** It decides whether anything is driven, whether a0=0 holds, and whether the Q_n=0 residue vanishes. No inverse is meaningful before that.
- **Unverified item.** `point-01-LEFT.json` around line 1761 has a "reproject" object (route `OMIT_FINITE_NATIVE_COEFFICIENT_THEN_REPROJECT_FULL_FIVE_ROW_SOURCE`). Its residual has a nonzero entry in row index 3, imag 0.2553, with movement 0.0416. I could not tell from the part I read whether this is a face drive or only a reprojection residual. Resolve it before relying on "zero selected face drives" (proposal lines 111–112).
- **Fourier receiving problem.** It plausibly avoids full mixed-kernel action. With the step profile, the source is a single-Q_n function times a pointwise 5×5 system, plus one output-momentum power integral. It is a proposal only. The uniform pencil cannot be inverted naively on the propagating pole or at the grazing branch point. The proposal already says so (lines 74–78).

## 4. Power balance
- **The native forms are insufficient as they stand.**
  - **Closure:** it is solved homogeneously. At `native:3070–3085`, `constant = closure.subs(amplitude,0)` is solved for the exterior amplitude with no external inhomogeneity.
  - **Row-power maps:** `PLUS/MINUS_ROW_POWER_MAP` (3315–3318) pair only the five slab residuals.
  - **Finite balance:** `finite_residual` (3399) contains only that `source_power`. Face kinematic, closure and exterior inhomogeneities are absent.
- **Source-specific extensions needed:**
  1. Work from inhomogeneous face kinematics and closure, and from the exterior drive.
  2. Retained induced work and direct exterior reconstruction if sources are eliminated into the slab rows, without double counting.
  3. Comparison against independently reconstructed face, bulk and lateral power using `C0Ψ1 + C1Ψ0`.
  4. Any external work of the held profile (LAB_HELD). Nonnegativity of the loss is not assumed.

## 5. Priority recommendation
1. Park the grade-11 J/direct recovery as operator-validation work. Do not discard it.
2. Build the projected first-order forcing and face sources on the LEFT doublet, including cross terms and the 10 and 01 grades. Resolve the row-3 residual and test a0=0 in each loss channel.
3. Examine the Q_n=0 residue against the grazing branch point.
4. Only then pose the Fourier receiving problem and add the inhomogeneous boundary work to the power balance.

**Physical blocker if the coefficient cannot be produced:** there is no established regular first-order loss-sector response at exact grazing with a step-end source. It is a missing physical premise, not a computing cost. The calibration and strong-edge number is a separate missing datum.