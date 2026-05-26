---
unit_id: 024
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 024

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:131-183`: the previous direct-rational `zResp`/`nResp` Series block (old lines 121-150) is removed and replaced with a per-pair conservative 2x2 matrix-inverse derivation. New `mPair = {{omegaU^2 - omega^2, -rPair}, {-rPair, omegaW^2 - omega^2}}`, `coupling = {gU, gW}`, `zFromMatrix = First[First[{coupling}.Inverse[mPair].Transpose[{coupling}]]]`, `nFromMatrix = ((Inverse[mPair].coupling)[[2]])^2`. The closed-form references `qRef, hRef, pRef, deltaRef, sRef, zRefRational, nRefRational` are introduced and two new anchor assertions appear: `expectZero["Z_full from matrix inverse matches paper rational", zFromMatrix - zRefRational]` and `expectZero["N_full from matrix inverse matches paper rational", nFromMatrix - nRefRational]`. The downstream `zSeries`/`nSeries` are now derived from `zFromMatrix`/`nFromMatrix` rather than from a pre-chosen rational. Coefficient assertions Z0–Z4, N0–N4, D0–D4 are restructured to reference `qRef/deltaRef`, `pRef^2/deltaRef^2`, etc.
- `mathematica/...wl:252-290` (Section V): `laneRatio[lam_]` is replaced with a `pAFull[lam_]`/`pAAt0[lam_]`/`pAFirst[lam_]` quotient-rule derivation using `D[pAFull[lam], eps] /. eps -> 0` rather than `Series[..., {eps, 0, 1}]`. The lane assertions are reframed as "P20/P21/P22 first-order from quotient rule" against `1*p1Closed`, `(1/2)*p1Closed`, `(-1)*p1Closed`. The downstream defect map (`pbarChk`, `aPChk`, `bPChk`) is rebuilt from `pAAt0 + eps*pAFirst`.

**Assessment:**
F1 is the central sign-resolution remediation for this stage. The directive specified `+rPair` off-diagonals, but Codex blocked F1/F2 in iter1 because the prescribed matrix inverts with the wrong-signed mixed term and the anchor assertions would have failed. After the documented codex remediation pass (session 019e4c8d), F1 was re-applied with `-rPair` off-diagonals, with the deviation block citing the `+R_l A_l W_l` Lagrangian term at `paper/parts/part01_parent_geometry.tex:956` and the matched mixed signs at `paper/stages/stage_024.tex:108,113`. Both anchor assertions now pass (`Z_full from matrix inverse matches paper rational = 0`, `N_full from matrix inverse matches paper rational = 0` on Mathematica log lines 46-49). The N projection uses `(Inverse[mPair].coupling)[[2]]` rather than the directive's literal `coupling.Inverse[mPair].{0,1}`; these are algebraically equivalent (both pick the W-component of the response vector) and the equivalence is borne out by the passing anchor.

The Section V quotient-rule replacement is faithful: it uses `D[..., eps]` as required, derives `pAFirst[lam]` symbolically and assembles the first-order expansion without reusing SymPy's variable names. All eight downstream assertions pass (Mathematica log lines 103-119).

Collateral edits beyond the directive's text, but documented in the Applied block and well-motivated:
1. Memoized canonical sphere moments `i4`/`i6` (Wolfram diff lines 5-22) with index-sort canonicalization to fix the performance hang introduced by the F3/F4-augmented Section IV table — this is a performance fix, not a logic change; the sphere integrals are unchanged values, just cached.
2. Section IV symbol-context reset (`ClearAll[gU, gW, rPair, omegaU, omegaW, mPair, coupling, zFromMatrix, ..., gUiso, gWiso, rrIso, omU, omW]` at .wl:227-230) before the `m20` table — necessary because F3/F4 leaked symbols would otherwise contaminate the spherical Integrate of Y20 against the basis. Section IV's M-M_target check still passes (Mathematica log line 90-91).

Both collateral edits are necessary plumbing for the F1/F4 additions to coexist with Section IV; they do not alter assertion semantics. No regressions.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:261-279`: inserted a new SymPy matrix-inverse anchor block immediately after the existing `P = sp.simplify(OmegaU**2 * GW + Rr * GU)` definition (and before the `Zresp = sp.expand(...)` series extraction). New code defines `M_pair = sp.Matrix([[OmegaU**2 - omega**2, -Rr], [-Rr, OmegaW**2 - omega**2]])`, `coupling = sp.Matrix([GU, GW])`, `eta_response = M_pair.inv() * coupling`, `Z_from_matrix = sp.simplify((coupling.T * eta_response)[0, 0])`, `N_from_matrix = sp.simplify(eta_response[1, 0] ** 2)`, with `Zresp_target = (Q - H*omega**2)/(Delta - S*omega**2 + omega**4)` and `Nresp_target = (P - GW*omega**2)**2/(Delta - S*omega**2 + omega**4)**2`. Two new `expect_zero` calls: `"Z_from_matrix - Zresp_target (physics anchor)"` and `"N_from_matrix - Nresp_target (physics anchor)"`.
- `scripts/...py:300-303`: the existing `subbanner("III.2 — One-pair Maxwell/mixed witness and D(omega)")` is followed by a three-line explanatory `print(...)` block tagging the section as non-tautological given the anchor.

**Assessment:**
Same sign resolution as F1: the directive specified `+Rr` off-diagonals on `M_pair` and Codex iter1 blocked because the assertion would have residual `4*Rr*GU*GW/Delta`. After the remediation pass, F2 was applied with `-Rr` off-diagonals, matching F1. Both anchor assertions pass (SymPy log lines 81-82). The existing `Zresp = sp.expand(...)` / `Nresp = sp.expand(...)` / `Dresp` blocks are unchanged in shape, so the original Z0–Z4, N0–N4, D0–D4 coefficient assertions still print zero (log lines 90-98). The added anchor is non-tautological: it tests that the algebraic combination `Q - H*omega^2` over `Delta - S*omega^2 + omega^4` is the actual matrix-inverse response, so a typo in `Q`, `H`, `P`, `Delta`, or `S` would now produce a non-zero residual.

The explanatory print block is informational only and matches the directive's prescription.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:191`: the `equal_lane = {x20: x0, x21: x0, x22: x0}` definition and the three `expect_zero` calls for `xbar - x0`, `a_x on equal lanes`, `b_x on equal lanes` are removed. The `subbanner("II.1 — Equal-lane data imply exact grouped isotropy")` is retained as a structural marker but the body is empty.
- `scripts/...py:201-213`: new arbitrary-mix reassembly block `mix = {x20: p, x21: q, x22: rr}` with three new `expect_zero` calls (`x20 reassembled - p`, `x21 reassembled - q`, `x22 reassembled - rr`) inserted after the II.2 witness block. All three pass (SymPy log lines 65-67).
- `scripts/...py:233`: the C_alpha self-equality `expect_zero("C_alpha - lambda_B,alpha I_etaalpha", ...)` call is removed; `C1 = sp.simplify(lamB1 * Ieta1)` and `C2 = sp.simplify(lamB2 * Ieta2)` definitions remain (still needed by `Bresp`).
- `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:89-103`: equivalent deletions of `xbar - x0`, `a_x on equal lanes`, `b_x on equal lanes`, and the C_alpha self-equality on line 102. Added mirror reassembly block with `xbarMix, axMix, bxMix, x20Re, x21Re, x22Re` and three new `expectZero` calls; all pass (Mathematica log lines 30-35).

**Assessment:**
Both tautologies are gone from both engines. The replacement reassembly block is non-tautological: the identities `xbar + 4*ax = x20`, `xbar - ax + bx = x21`, `xbar - ax - bx = x22` test the inverse of the `(1,2,2)/5`-weighted decomposition map — they would fail under a typo in any of the linear-combination coefficients. The auditor's algebra (directive lines 273-275) for the reassembly coefficients is correct (verified independently: `xbar + 4ax = (x20+2x21+2x22)/5 + (2x20-x21-x22)*4/10 = (10*x20)/10 = x20`). No collateral edits beyond the directive scope.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:318-385`: new top-level function `lane_collapse_check()` inserted between `isotropic_overlap_moments` and `axisymmetric_splitting`. Defines lane-decorated couplings `GU_{20,21,22}`, `GW_{20,21,22}`, `Rr_{20,21,22}`, `Cc_{20,21,22}` and inner helpers `per_lane_D` and `per_lane_N` matching the directive's prescription verbatim. Builds `D_20, D_21, D_22, N_20, N_21, N_22` and substitutes `iso` to verify three `D` collapse identities, three `N` collapse identities (subsections III.5.1, III.5.2), and a lane-breaking witness (III.5.3) that asserts a nonzero linear-in-delta coefficient of `D_20 - D_21` under `GU_20 -> GU_iso + delta`.
- `scripts/...py:484`: `lane_collapse_check()` added to `__main__` between `isotropic_overlap_moments()` and `axisymmetric_splitting()`.
- `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:185-223`: mirror block inserted before the `SECTION IV — AXISYMMETRIC SPLITTING MATRIX` banner. Uses `perLaneD[gUA_, gWA_, rrA_, cA_]`, `perLaneN[gUA_, gWA_, rrA_]`, with `isoSubs` and `breakSubs` substitution rules. Same six collapse assertions and one lane-breaking witness with `Print["linear coefficient of delta in (D_20 - D_21) = ", fmt[coeffDelta]]` plus an `If[TrueQ[coeffDelta === 0], fail[...], pass[...]]`.

**Assessment:**
All six collapse assertions pass on both engines (SymPy log lines 107-116; Mathematica log lines 72-83). The lane-breaking witness produces a non-trivial coefficient on both sides:
- SymPy (log line 121): `2*(GU_iso*Omega_W**2 - GU_iso*omega**2 + GW_iso*Rr_iso) / (-Omega_U**2*Omega_W**2 + Omega_U**2*omega**2 + Omega_W**2*omega**2 + Rr_iso**2 - omega**4)`
- Mathematica (log line 84): `(2*gUiso*(omega - omW)*(omega + omW) - 2*gWiso*rrIso) / ((omega - omU)*(omega + omU)*(omega - omW)*(omega + omW) - rrIso^2)`

These are algebraically equal after expanding `(omega - omW)(omega + omW) = omega^2 - omW^2`. Both engines confirm the witness is linear-in-delta and nonzero, so the III.5 collapse assertions are demonstrably non-tautological. No collateral edits beyond what the directive prescribed.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 81: `Z_from_matrix - Zresp_target (physics anchor) = 0`
- Line 82: `N_from_matrix - Nresp_target (physics anchor) = 0`
- Line 121: `linear coefficient of delta in (D_20 - D_21) = 2*(GU_iso*Omega_W**2 - GU_iso*omega**2 + GW_iso*Rr_iso)/(...)`
- Line 194: `# exit_code: 0`

Section II.1 banner is printed but body is empty (the equal-lane tautologies were deleted as required). The transcript flows directly to the II.2 witnesses, then the reassembly block, then III.1 (B_n), then the F2 anchor lines, then III.2, then the new III.5 section. All structurally correct.

**Mathematica:** exit=0. Notable lines:
- Line 46-47: `Z_full from matrix inverse matches paper rational = 0` / `PASS: Z_full from matrix inverse matches paper rational`
- Line 48-49: `N_full from matrix inverse matches paper rational = 0` / `PASS: N_full from matrix inverse matches paper rational`
- Line 84-85: `linear coefficient of delta in (D_20 - D_21) = (2*gUiso*(omega - omW)*(omega + omW) - 2*gWiso*rrIso)/(...)` / `PASS: Lane-breaking witness: collapse check is non-tautological`
- Line 127: `# exit_code: 0`

All 38 PASS-labeled assertions pass. The Section IV M-M_target sphere-integral check still produces the 5x5 zero matrix (log line 90), confirming the F1/F4 symbol-context resets did not regress Section IV.

**Output freshness:** The saved canonical outputs (`scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt` mtime 2026-05-21 16:03; `mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt` mtime 2026-05-21 16:53) are STALE relative to the May 26 exec logs and post-fix script mtimes (.py 2026-05-26 00:23; .wl 2026-05-26 00:28). The exec logs at `redteam/exec_logs/` are the post-fix runs and both report exit 0, so the verification itself stands, but the canonical `output/*.txt` mirrors need to be regenerated (orchestrator should refresh via `$RT exec-sympy 024` / `$RT exec-mathematica 024` to update MANIFEST.yaml).

## Material-change assessment

`material_change`: false.

Stage 024's load-bearing closed forms (`B_n`, `Z_n`, `N_n`, `D_n`, `M^(20)`, `P_1` formula, `b=3a`) are unchanged. The matrix-inverse anchor (F1, F2) shows that the existing rational form is the actual response from the conservative 2x2 matrix with `-R` off-diagonals — but the coefficient identities the script asserts (Z0 = Q/Delta, etc.) are exactly the same expressions as before. The lane-decorated `D_{A,n}, N_{A,n}` construction in III.5 (F4) is a new internal check that consumes the same closed forms; it does not alter any output. Section IV's M^(20) matrix is computed by the same triple_overlap path and still produces the diagonal `(sqrt(5)/(7*sqrt(pi))) * (1, 1/2, 1/2, -1, -1)`. Section V's quotient-rule derivation produces the same `P_1 = (N_1*D_0 - N_0*D_1)/D_0^2` formula as before.

Downstream units (025+) consume `B_n, Z_n, N_n, D_n, P_1, kappa_*`; none of these symbolic results changed. The sign convention on the conservative matrix off-diagonal (`-R`) was made explicit in this stage but the paper rational `Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2` (positive cross term) is unchanged — the `-R` matrix inverts to `+2 G_U G_W R / Delta` cross term, which is exactly the paper's `Q/Delta` at omega=0. Therefore downstream stages see the same numerics.

## Side observations (non-blocking)

1. **Output txt freshness gap:** The canonical `output/*.txt` files are 5 days older than the script files. The orchestrator should refresh these (via direct `python3 ... > output/...txt` and `math -script ... > output/...txt`, or via the single-seat `$RT exec-*` calls per the no-parallel rule) before the MANIFEST.yaml is consulted for stage 024. The post-fix run output is captured in `redteam/exec_logs/` but not mirrored to the canonical output path.

2. **Section II.1 empty banner:** After the F3 deletions, the SymPy `subbanner("II.1 — Equal-lane data imply exact grouped isotropy")` prints a banner with no assertions following it. This is cosmetically odd (a section header with no body), but does not affect correctness. A future cleanup could either delete the empty banner or rename II.1 to describe what II.2 + the reassembly block are now testing.

3. **Sign convention deviation (F1/F2):** The deviation from the literal directive prescription (`-Rr` instead of `+Rr` off-diagonal) is well-documented in the Applied blocks and grounded in `paper/parts/part01_parent_geometry.tex:956`. Under a strict scripts-only contract, the verifier cannot independently confirm the paper citation. However, the anchor assertions now pass — meaning the matrix-inverse derivation does land on the paper's stated rational closed form — which is the operationally relevant outcome. Stage 023 had the same `Mblock` sign question and was blocked_unfixable; resolving it for stage 024 via the same matrix structure suggests stage 023 should be re-audited under the now-resolved convention as a follow-up. Flag for orchestrator.

4. **Memoization of i4/i6:** The canonical-sort memoization preserves the integration result (the integrals are symmetric in their indices, so canonical-sort is mathematically safe). This is a performance-only change.

## Verdict justification

All four findings are `resolved`: F1 closed via Section III matrix-inverse derivation and Section V quotient-rule replacement in Mathematica (with documented `-R` sign resolution), F2 closed via the SymPy matrix-inverse anchor block, F3 closed via deletion of tautologies and addition of falsifiable reassembly witnesses, and F4 closed via the new `lane_collapse_check` / III.5 section in both engines. Both engines exit 0, all assertion lines print zero residuals, and the lane-breaking witness produces a documented non-zero linear-in-delta coefficient on both sides. The sign deviation in F1/F2 is documented and operationally validated by the passing anchors. No regressions and no material changes to downstream-consumed symbolic outputs. The canonical `output/*.txt` files need regeneration but the substantive verification stands on the post-fix exec logs.
