---
unit_id: 068
batch: III.3
created_at: 2026-05-26T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-26T12:57:00-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive - unit 068

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead - skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 - tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:43-68`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:28-43`

**Issue:**
Three "derivations" at the top of each script are tautological:

(a) `W_res = C^2 W_wall` is constructed by relabeling `(C*A_in)^2 -> C^2*W_wall` then compared to `C^2*W_wall`. In Mathematica, Solve produces the exact RHS by definition.
(b) `P_res = 1/C_res^2` is constructed by writing `1/Cres**2` then compared to `1/Cres**2`.
(c) The "consistency check" `P_res*C_res^2 - 1 = 0` is literally `(1/Cres**2)*Cres**2 - 1 = 1 - 1 = 0`; the symbol `Pres` is not even in the LHS.

**Required change:**

### SymPy edits (file: `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`)

Replace lines 44-68 (everything between the `# 1. Resonance-corrected wall figure...` header and the `# 2. Exact threshold translation...` header, exclusive of section-2's leading comment) with the following block. This anchors `W_res = C^2 W_wall` in the matched-branch gain decomposition from notes section 1, derives `P_res = 1/C_res^2` from the required-wall-figure amplification, and adds a numeric anchor against the paper card's quoted value.

Before:
```python
# ---------------------------------------------------------------------------
# 1. Resonance-corrected wall figure (derived, not postulated)
# ---------------------------------------------------------------------------
#
# Premise: a wall-incident wave of complex amplitude A_in encounters a
# transmission coefficient C(r) (real, positive at the resonance branch
# we consider).  Power scales as |amplitude|^2, so the transmitted wall
# power figure is |C|^2 * W_wall.  We carry C2 := C^2 as the audit symbol.

A_in, C_sym = sp.symbols("A_in C", positive=True, real=True)
A_trans = C_sym * A_in
# Power in the incident wave ~ |A_in|^2 normalises to W_wall; transmitted power:
W_wall_def = A_in**2
W_res_derived = sp.simplify((A_trans**2).subs(A_in**2, W_wall) * 1)
# Re-express |C|^2 as the audit symbol C2:
W_res_derived = sp.simplify(W_res_derived.subs(C_sym**2, C2))
expect_zero("W_res - C2 * W_wall", W_res_derived - C2 * W_wall)

# At the resonance point r = r_*, C^2 -> C_res^2 by definition of the
# resonance branch; the penalty factor is the *inverse* of the amplification.
Cres = sp.symbols("C_res", positive=True, real=True)
Pres_derived = sp.simplify(1 / Cres**2)
expect_zero("P_res - 1/C_res^2", Pres_derived - 1/Cres**2)
# And consistency with the audit-level symbol P_res introduced as positive:
expect_zero("P_res*C_res^2 - 1", (1/Cres**2)*Cres**2 - 1)
```

After:
```python
# ---------------------------------------------------------------------------
# 1. Resonance-corrected wall figure (derived from matched-branch gain)
# ---------------------------------------------------------------------------
#
# Premise (notes section 1): the matched-branch gain is
#   G_match = rho_star * g_phi^2 * N_phiphi / (m * c_s^2 * K_X),
# and the wall figure of merit is W_wall = kappa * G_match.
# The profile-family gain at coherence C^2 is G_res = C^2 * G_match,
# so W_res = kappa * G_res. We verify that W_res = C^2 * W_wall by
# building each side from these independent component symbols.

rho_star, g_phi, N_phiphi, m_s, c_s, K_X, kappa = sp.symbols(
    "rho_star g_phi N_phiphi m c_s K_X kappa", positive=True, real=True
)
Gmatch_expr = rho_star * g_phi**2 * N_phiphi / (m_s * c_s**2 * K_X)
Wwall_expr = kappa * Gmatch_expr
Gres_expr  = C2 * Gmatch_expr
Wres_expr  = kappa * Gres_expr
expect_zero("W_res - C2 * W_wall (from gain decomposition)",
            sp.simplify(Wres_expr - C2 * Wwall_expr))
# Also confirm the matched limit C2 -> 1 collapses W_res to W_wall:
expect_zero("W_res(C2->1) - W_wall (matched limit)",
            sp.simplify(Wres_expr.subs(C2, 1) - Wwall_expr))

# ---------------------------------------------------------------------------
# 1b. Resonance penalty factor (derived from ratio of required wall figures)
# ---------------------------------------------------------------------------
#
# The required wall figure on the matched branch is W_wall = Pe_req/Delta.
# On the profile family at coherence C^2, it is W_wall = Pe_req/(C^2 Delta).
# The amplification of the *required* wall figure from matched to profile
# is therefore 1/C^2.  At the resonance point C^2 -> C_res^2 (carried in
# from stage 067), this amplification factor is P_res = 1/C_res^2.

Cres, Delta_sym_pre = sp.symbols("C_res Delta", positive=True, real=True)
W_required_match = sp.solve(
    sp.Eq(sp.Symbol("Wm", positive=True) * Delta_sym_pre, Pe_req),
    sp.Symbol("Wm", positive=True),
)[0]
W_required_prof = sp.solve(
    sp.Eq(C2 * sp.Symbol("Wp", positive=True) * Delta_sym_pre, Pe_req),
    sp.Symbol("Wp", positive=True),
)[0]
Pres_from_ratio = sp.simplify((W_required_prof / W_required_match).subs(C2, Cres**2))
expect_zero("P_res - 1/C_res^2 (from required-wall-figure ratio)",
            Pres_from_ratio - 1/Cres**2)

# Numeric anchor: the paper card states P_res = 1.005612487760576 and
# C_res^2 = 0.994418836451529 (carried from stage 067). Verify the link.
Cres_sq_numeric = sp.Float("0.994418836451529", 20)
Pres_numeric    = 1 / Cres_sq_numeric
Pres_paper      = sp.Float("1.005612487760576", 20)
numeric_residual = sp.Abs(Pres_numeric - Pres_paper)
print(f"P_res numeric residual = {numeric_residual}")
if numeric_residual >= sp.Float("1e-12", 20):
    raise AssertionError(
        f"P_res numeric anchor failed: |1/C_res^2 - paper P_res| = {numeric_residual}"
    )
```

Note: the new code introduces local symbols `Wm`, `Wp` inside the two `sp.solve` calls; they are scoped to those lines and do not collide with the existing `W_match_sym`, `W_prof_sym` used in section 2.

### Mathematica edits (file: `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`)

Replace lines 33-43 (the comment block plus the W_res and P_res derivations) with a non-parallel derivation. See F2 for the broader transliteration concern; the below specifically removes the tautologies and uses Mathematica idioms (`FullSimplify`, `Reduce`) rather than mimicking SymPy's `Solve` choreography.

Before:
```mathematica
(* Premise: transmission coefficient C maps wall amplitude Ain to C*Ain.   *)
(* Power = amplitude^2 ; W_wall normalises to Ain^2 ; transmitted power:   *)
(* |C * Ain|^2 = C^2 * W_wall. Use Reduce to extract the C^2 coefficient.  *)
WresRule = First@Solve[Wres == (C2)*Wwall, Wres];
WresDerived = Wres /. WresRule;
expectZero["W_res - C2 * W_wall", WresDerived - C2*Wwall];

(* P_res = 1/C_res^2 derived as the inverse of the amplification *)
PresFromCres = First@Solve[Pres*Cres^2 == 1, Pres];
PresDerived = Pres /. PresFromCres;
expectZero["P_res - 1/C_res^2", PresDerived - 1/Cres^2];
```

After:
```mathematica
(* W_res derived from matched-branch gain decomposition (notes section 1): *)
(*   G_match = rho_star g_phi^2 N_phiphi / (m c_s^2 K_X),                  *)
(*   W_wall  = kappa G_match,    G_res = C^2 G_match,    W_res = kappa G_res *)
Clear[rhoStar, gPhi, Nphi, mPart, cs, KX, kappaW, GmatchExpr, WwallExpr, GresExpr, WresExpr];
$Assumptions = $Assumptions && Element[{rhoStar, gPhi, Nphi, mPart, cs, KX, kappaW}, Reals] &&
  rhoStar > 0 && gPhi > 0 && Nphi > 0 && mPart > 0 && cs > 0 && KX > 0 && kappaW > 0;
GmatchExpr = rhoStar*gPhi^2*Nphi / (mPart*cs^2*KX);
WwallExpr  = kappaW*GmatchExpr;
GresExpr   = C2*GmatchExpr;
WresExpr   = kappaW*GresExpr;
expectZero["W_res - C2 * W_wall (gain decomposition)",
  FullSimplify[WresExpr - C2*WwallExpr, Assumptions -> $Assumptions]];
expectZero["W_res(C2->1) - W_wall (matched limit)",
  FullSimplify[(WresExpr /. C2 -> 1) - WwallExpr, Assumptions -> $Assumptions]];

(* P_res derived from ratio of required wall figures at resonance C2 -> Cres^2. *)
PresFromRatio =
  FullSimplify[
    ((PeReq/(C2*Delta1)) / (PeReq/Delta1)) /. C2 -> Cres^2,
    Assumptions -> $Assumptions && Delta1 > 0
  ];
expectZero["P_res - 1/C_res^2 (required-wall-figure ratio)",
  FullSimplify[PresFromRatio - 1/Cres^2, Assumptions -> $Assumptions]];

(* Numeric anchor: paper card states P_res = 1.005612487760576 and             *)
(* C_res^2 = 0.994418836451529 (carried from stage 067).                       *)
With[{CresSqNum = SetPrecision[0.994418836451529, 20],
      PresPaperNum = SetPrecision[1.005612487760576, 20]},
  PresNumResidual = Abs[1/CresSqNum - PresPaperNum];
  Print["P_res numeric residual = ", fmt[PresNumResidual]];
  If[!TrueQ[PresNumResidual < 10^-12],
    fail["P_res numeric anchor", PresNumResidual]
  ];
];
```

**Verification:**
After Codex applies, the verifier runs `redteam exec-sympy 068` and `redteam exec-mathematica 068` and confirms:
1. Both scripts still exit 0.
2. The new printed lines `W_res - C2 * W_wall (from gain decomposition) = 0`, `W_res(C2->1) - W_wall (matched limit) = 0`, `P_res - 1/C_res^2 (from required-wall-figure ratio) = 0`, and `P_res numeric residual = ...` (a value below 1e-12) appear in both transcripts.
3. The line `expect_zero("P_res*C_res^2 - 1", (1/Cres**2)*Cres**2 - 1)` no longer appears anywhere.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`
- summary: Replaced the tautological W_res and P_res checks with gain-decomposition, required-wall-ratio, and numeric-anchor checks in both scripts.
- deviation: Mathematica check labels were aligned to the directive's verification inventory by adding "from" to two labels.

## F2 - mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:46-84`

**Issue:**
After F1's edits replace the top-of-script W_res / P_res blocks, the remaining sections of the `.wl` (matched-branch threshold via `WmatchSol = First@Solve[Wmatch*DeltaSym == PeReq, Wmatch]`, profile-family threshold via `WprofSol = First@Solve[C2*Wprof*DeltaSym == PeReq, Wprof]`, and "Way B" band width via `gap /. First@Solve[WsuffMatch + gap == Pres*WsuffMatch, gap]`) still mirror the SymPy script's algebra line-by-line. The leftover `gapSym;` on line 74 confirms the port-from-Python provenance. The two engines must derive the result along non-parallel symbolic paths.

**Required change:**

Make the following targeted edits to `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`:

(i) **Remove the stray `gapSym;` line (current line 74).** Delete the entire line.

(ii) **Replace the matched-branch and profile-family threshold derivations (current lines 45-53)** with `Reduce`-based derivations:

Before:
```mathematica
(* Matched-branch threshold from W_match * Delta = Pe_req *)
WmatchSol = First@Solve[Wmatch*DeltaSym == PeReq, Wmatch];
WfailMatch = FullSimplify[(Wmatch /. WmatchSol) /. DeltaSym -> Deltainf];
WsuffMatch = FullSimplify[(Wmatch /. WmatchSol) /. DeltaSym -> Delta0];

(* Profile-family threshold from C2 * W_prof * Delta = Pe_req *)
WprofSol = First@Solve[C2*Wprof*DeltaSym == PeReq, Wprof];
WfailRes = FullSimplify[(Wprof /. WprofSol) /. DeltaSym -> Deltainf];
WsuffRes = FullSimplify[(Wprof /. WprofSol) /. DeltaSym -> Delta0];
```

After:
```mathematica
(* Matched-branch threshold from Reduce[W*Delta == PeReq && W > 0, W].       *)
(* Use Reduce rather than Solve to keep the positivity premise explicit.     *)
WfailMatch = First[Cases[
    Reduce[Wmatch*Deltainf == PeReq && Wmatch > 0 && Deltainf > 0 && PeReq > 0, Wmatch, Reals],
    HoldPattern[Wmatch == rhs_] :> rhs, Infinity
  ]];
WsuffMatch = First[Cases[
    Reduce[Wmatch*Delta0 == PeReq && Wmatch > 0 && Delta0 > 0 && PeReq > 0, Wmatch, Reals],
    HoldPattern[Wmatch == rhs_] :> rhs, Infinity
  ]];
WfailMatch = FullSimplify[WfailMatch, Assumptions -> $Assumptions];
WsuffMatch = FullSimplify[WsuffMatch, Assumptions -> $Assumptions];

(* Profile-family threshold from Reduce[C2*W*Delta == PeReq && W > 0, W].   *)
WfailRes = First[Cases[
    Reduce[C2*Wprof*Deltainf == PeReq && Wprof > 0 && Deltainf > 0 && PeReq > 0 && C2 > 0, Wprof, Reals],
    HoldPattern[Wprof == rhs_] :> rhs, Infinity
  ]];
WsuffRes = First[Cases[
    Reduce[C2*Wprof*Delta0 == PeReq && Wprof > 0 && Delta0 > 0 && PeReq > 0 && C2 > 0, Wprof, Reals],
    HoldPattern[Wprof == rhs_] :> rhs, Infinity
  ]];
WfailRes = FullSimplify[WfailRes, Assumptions -> $Assumptions];
WsuffRes = FullSimplify[WsuffRes, Assumptions -> $Assumptions];
```

(iii) **Replace the "Way B" band-width Solve calls (current lines 73-76)** with a `Series`-expansion derivation:

Before:
```mathematica
(* Way B: Solve WsuffMatch + gap == Pres*WsuffMatch for gap. *)
gapSym;
successBandB = gap /. First@Solve[WsuffMatch + gap == Pres*WsuffMatch, gap];
failureBandB = gap /. First@Solve[WfailMatch + gap == Pres*WfailMatch, gap];
```

After:
```mathematica
(* Way B: first-order Series expansion of Pres*Wmatch around Pres = 1. *)
(* The first-order coefficient in (Pres - 1) is the band width.        *)
successSeries = Series[Pres*WsuffMatch, {Pres, 1, 1}];
failureSeries = Series[Pres*WfailMatch, {Pres, 1, 1}];
successBandB = FullSimplify[Coefficient[Normal[successSeries] - WsuffMatch, Pres - 1, 0]
                            + Coefficient[Normal[successSeries] - WsuffMatch, Pres - 1, 1]*(Pres - 1),
                            Assumptions -> $Assumptions];
failureBandB = FullSimplify[Coefficient[Normal[failureSeries] - WfailMatch, Pres - 1, 0]
                            + Coefficient[Normal[failureSeries] - WfailMatch, Pres - 1, 1]*(Pres - 1),
                            Assumptions -> $Assumptions];
```

This re-derives the band width as a linear-order Taylor coefficient rather than as a Solve for `gap`. The final expressions should still equal `PeReq*(Pres - 1)/Delta0` and `PeReq*(Pres - 1)/Deltainf` respectively. (If the Series-and-Coefficient form returns the same expression directly, the FullSimplify will collapse it.)

**Verification:**
After Codex applies, the verifier confirms:
1. `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl` no longer contains the substring `gapSym`.
2. The script no longer contains `First@Solve[Wmatch*DeltaSym ==` (replaced by Reduce-and-Cases).
3. The script no longer contains `First@Solve[WsuffMatch + gap` (replaced by Series).
4. The script still exits 0 and the final printed expressions for `WfailRes`, `WsuffRes`, `successBandA`, `failureBandA`, `successBandB`, `failureBandB` agree with the SymPy outputs at the algebraic level (the comparison expectZero calls on `WfailRes*C2 - WfailMatch` etc. all still pass).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`
- summary: Replaced the Mathematica threshold derivations with Reduce-and-Cases forms and removed the Python-port gap solve path.
- deviation: F3 supersedes the interim Series-based band-width block with the requested C-form/P-form replacement.

## F3 - insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:113-134`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:68-84`

**Issue:**
The "two ways" band-width check is not actually two independent computations. Way A and Way B both arithmetically collapse to `(Pres - 1) * Wmatch` and do not exercise `Wsuff_res / Wfail_res` independently. The comparison fails to be a cross-check.

**Required change:**

### SymPy edits (file: `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`)

After F1 + F2 are applied, the band-width block (currently SymPy lines 113-134) should be replaced with an alternative that routes through `Cres` instead of just substituting `Pres -> 1/Cres^2`. This makes the cross-check sensitive to the `Pres = 1/Cres^2` relation itself.

Before (current lines 113-134):
```python
# ---------------------------------------------------------------------------
# 3. Profile-sensitive band widths, computed two independent ways
# ---------------------------------------------------------------------------
#
# Way A: as the difference of the F1-derived profile and matched thresholds,
#        evaluated at the resonance point C^2 = 1/P_res.
success_band_widthA = sp.simplify((Wsuff_res - Wsuff_match).subs(C2, 1 / Pres))
failure_band_widthA = sp.simplify((Wfail_res - Wfail_match).subs(C2, 1 / Pres))

# Way B: by solving Wsuff_res = Pres * Wsuff_match for the gap directly.
#        We pose the equation and let Solve extract the difference symbolically.
gap_sym = sp.symbols("gap", real=True)
success_band_widthB = sp.solve(sp.Eq(Wsuff_match + gap_sym, Pres * Wsuff_match), gap_sym)[0]
failure_band_widthB = sp.solve(sp.Eq(Wfail_match + gap_sym, Pres * Wfail_match), gap_sym)[0]

print("Success-side band width (A) =", success_band_widthA)
print("Failure-side band width (A) =", failure_band_widthA)
print("Success-side band width (B) =", success_band_widthB)
print("Failure-side band width (B) =", failure_band_widthB)

# Independent-derivation cross-check: the two ways must agree.
expect_zero("success band A vs B", success_band_widthA - success_band_widthB)
expect_zero("failure band A vs B", failure_band_widthA - failure_band_widthB)
```

After:
```python
# ---------------------------------------------------------------------------
# 3. Profile-sensitive band widths, two routes through different symbols
# ---------------------------------------------------------------------------
#
# Way A (C-form): compute the difference at the resonance point using Cres
#   directly. Wsuff_res evaluated at C^2 = Cres^2 is Pe_req/(Cres^2 Delta_0).
# Way B (P-form): compute the difference using Pres as the audit symbol.
# Both should equal (Pres - 1)*Wsuff_match after the relation Pres = 1/Cres^2.
# Any perturbation of that relation breaks one form but not the other, so the
# cross-check is now sensitive to the Pres = 1/Cres^2 link.

# Cres is defined in section 1b above; reuse it here.
Wsuff_res_C   = Pe_req / (Cres**2 * Delta0)
Wfail_res_C   = Pe_req / (Cres**2 * Deltainf)
success_band_widthA = sp.simplify(Wsuff_res_C - Wsuff_match)
failure_band_widthA = sp.simplify(Wfail_res_C - Wfail_match)

success_band_widthB = sp.simplify((Pres - 1) * Wsuff_match)
failure_band_widthB = sp.simplify((Pres - 1) * Wfail_match)

print("Success-side band width (C-form) =", success_band_widthA)
print("Failure-side band width (C-form) =", failure_band_widthA)
print("Success-side band width (P-form) =", success_band_widthB)
print("Failure-side band width (P-form) =", failure_band_widthB)

# Cross-check: under the relation Pres = 1/Cres^2, the two forms must agree.
expect_zero("success band C-form vs P-form (under Pres = 1/Cres^2)",
            sp.simplify((success_band_widthA - success_band_widthB).subs(Pres, 1/Cres**2)))
expect_zero("failure band C-form vs P-form (under Pres = 1/Cres^2)",
            sp.simplify((failure_band_widthA - failure_band_widthB).subs(Pres, 1/Cres**2)))
```

### Mathematica edits (file: `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`)

Replace the band-width block (after F2's Series-based "Way B" has already been written). The substitution to make is the same conceptual one: Way A becomes a `Cres`-form, Way B becomes a `Pres`-form, and the cross-check substitutes `Pres -> 1/Cres^2` before asserting equality.

Before (after F2 has been applied; the original lines 68-84 will look like):
```mathematica
banner["PROFILE-SENSITIVE BANDS"];
(* Way A: difference of profile and matched thresholds at C2 = 1/Pres. *)
successBandA = FullSimplify[(WsuffRes - WsuffMatch) /. C2 -> 1/Pres, Assumptions -> $Assumptions];
failureBandA = FullSimplify[(WfailRes - WfailMatch) /. C2 -> 1/Pres, Assumptions -> $Assumptions];

(* Way B: first-order Series expansion of Pres*Wmatch around Pres = 1. *)
...
Print["Success-side band width (A) = ", fmt[successBandA]];
Print["Failure-side band width (A) = ", fmt[failureBandA]];
Print["Success-side band width (B) = ", fmt[successBandB]];
Print["Failure-side band width (B) = ", fmt[failureBandB]];

expectZero["success band A vs B", successBandA - successBandB];
expectZero["failure band A vs B", failureBandA - failureBandB];
```

After:
```mathematica
banner["PROFILE-SENSITIVE BANDS"];
(* C-form: evaluate WsuffRes at C2 -> Cres^2 (NOT via Pres substitution). *)
WsuffResC = PeReq/(Cres^2 * Delta0);
WfailResC = PeReq/(Cres^2 * Deltainf);
successBandA = FullSimplify[WsuffResC - WsuffMatch, Assumptions -> $Assumptions];
failureBandA = FullSimplify[WfailResC - WfailMatch, Assumptions -> $Assumptions];

(* P-form: (Pres - 1) * Wmatch directly. *)
successBandB = FullSimplify[(Pres - 1)*WsuffMatch, Assumptions -> $Assumptions];
failureBandB = FullSimplify[(Pres - 1)*WfailMatch, Assumptions -> $Assumptions];

Print["Success-side band width (C-form) = ", fmt[successBandA]];
Print["Failure-side band width (C-form) = ", fmt[failureBandA]];
Print["Success-side band width (P-form) = ", fmt[successBandB]];
Print["Failure-side band width (P-form) = ", fmt[failureBandB]];

expectZero["success band C-form vs P-form (under Pres = 1/Cres^2)",
  FullSimplify[(successBandA - successBandB) /. Pres -> 1/Cres^2, Assumptions -> $Assumptions]];
expectZero["failure band C-form vs P-form (under Pres = 1/Cres^2)",
  FullSimplify[(failureBandA - failureBandB) /. Pres -> 1/Cres^2, Assumptions -> $Assumptions]];
```

(Note: the F2 Series-based "Way B" is superseded by this F3 reformation. The Series approach was an interim independence improvement; F3's `Pres`-form vs `Cres`-form is the stronger cross-check. Apply F3's replacement directly over whatever F2 left.)

**Verification:**
After Codex applies, the verifier confirms:
1. Both scripts still exit 0.
2. The printed expressions for `Success-side band width (C-form)` and `Success-side band width (P-form)` both equal `Pe_req*(P_res - 1)/Delta_0` after the substitution `Pres -> 1/Cres**2` (or equivalently, `Pe_req*(1 - Cres**2)/(Cres**2 * Delta_0)`).
3. The final assertion line numbers reference "C-form vs P-form", not "A vs B".
4. The `Wsuff_res` / `Wfail_res` symbols defined in section 2 are still printed unchanged earlier in the transcript; section 3 simply uses an alternate symbolic route.

---

After all three findings are applied, the verifier will run both engines to fresh transcripts and re-check the report's assertion inventory against the new line numbers.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`
- summary: Replaced the band-width comparison with C-form and P-form routes tied together only by the P_res = 1/C_res^2 relation.
- deviation: none
