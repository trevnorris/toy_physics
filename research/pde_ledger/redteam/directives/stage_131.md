---
unit_id: 131
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 131

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — missing_verification_script (subtype: script_doesnt_cover_claim)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py`

**Issue:**
The SymPy "audit" script has zero `assert` statements. It computes `Pi_*`, `V1_*`, `g'(Pi_*)`, and a "Parent bias mismatch formula", but only prints them — nothing can fail. The output transcript shows `Status: PASS` purely because nothing is tested.

**Required change:**

Append the following block at the end of the file (after the existing line 26 `print` of the threshold residual), keeping all current lines intact:

```python
# --- Anchored assertions vs notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md ---

# Anchor (1): Pi_* matches the notes-quoted value 1.50882951349316 (Sec. 1).
Pi_star_paper = sp.Float("1.50882951349316", 50)
Pi_star_num = sp.N(Pi_star, 50)
assert abs(Pi_star_num - Pi_star_paper) < sp.Float("1e-14", 50), (
    f"Pi_* mismatch vs notes Sec. 1: computed {Pi_star_num}, paper {Pi_star_paper}"
)
print("PASS: Pi_* matches notes Sec. 1 value 1.50882951349316")

# Anchor (2): g'(Pi_*) matches the notes-quoted slope 0.0714453558083195 (Sec. 3).
gprime_paper = sp.Float("0.0714453558083195", 50)
gprime_num = sp.N(gprime_star, 50)
assert abs(gprime_num - gprime_paper) < sp.Float("1e-14", 50), (
    f"g'(Pi_*) mismatch vs notes Sec. 3: computed {gprime_num}, paper {gprime_paper}"
)
print("PASS: g'(Pi_*) matches notes Sec. 3 value 0.0714453558083195")

# Anchor (3): parent threshold identity at Pi = Pi_*.
# notes Sec. 2:  T_m - q_* A_0' = Pi_* * Theta_sigma / L
threshold_at_star = threshold_residual.subs(Pi, sp.N(Pi_star, 50))
expected_form = (Tm - qstar*A0p) - sp.N(Pi_star, 50)*Theta_sigma/L
assert sp.simplify(threshold_at_star - expected_form) == 0, (
    "parent threshold identity at Pi_* does not match (T_m - q_* A_0') - Pi_* Theta_sigma/L"
)
print("PASS: parent threshold identity at Pi = Pi_* matches notes Sec. 2")

# Anchor (4): lower-branch discrimination — Pi_* sits on the g_- branch, NOT on a singular point.
# A point clearly away from Pi_* (e.g. 2*Pi_*) must give a residual visibly far from zero.
gPi_offstar = gPi.subs(Pi, 2*sp.N(Pi_star, 30))
offstar_residual = abs(sp.N(gPi_offstar - g_minus, 30))
assert offstar_residual > sp.Float("1e-3", 30), (
    f"counter-example failed: gPi(2*Pi_*) residual vs g_- = {offstar_residual}, "
    "expected clearly nonzero (lower-branch discrimination, paper Checks item 3)"
)
print(f"PASS: lower-branch discrimination — gPi(2*Pi_*) - g_- = {offstar_residual}")
```

**Claim manifest:**
- M1: `Pi_* = 1.50882951349316` (from notes line 8 / Sec. 1).
- M2: `g'(Pi_*) = 0.0714453558083195` (from notes line 93 / Sec. 3).
- M3: `T_m - q_* A_0' = Pi_* Theta_sigma / L` (parent threshold identity, paper card body + notes Sec. 2).
- M4: `Pi_*` lies on the lower (`g_-^{F1}`) branch, not on a singular branch (paper Checks item 3).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 131` and confirm: (a) script exits 0; (b) transcript contains four new `PASS:` lines matching M1–M4.

---

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:26,46`

**Issue:**
Line 46's `expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20]` is tautological: `piStar` was just defined at line 35 as the root of `FindRoot[gPi == gMinus, ...]`, so the residual is forced near zero by construction. The check verifies FindRoot's convergence, not any physics claim. Separately, line 26's banner string contains `"STAGE 114"` but this is stage 131.

**Required change:**

1. Replace line 26:
   ```
   banner["STAGE 114 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH COMPENSATION"];
   ```
   with
   ```
   banner["STAGE 131 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH COMPENSATION"];
   ```

2. Replace line 46:
   ```
   expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20];
   ```
   with the following block:
   ```mathematica
   (* Anchored checks vs notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md *)

   (* Anchor (1): Pi_* matches notes Sec. 1 quoted value. *)
   expectApprox["Pi_* notes Sec. 1 value",
     N[piStar, 50], 1.50882951349316`50, 10^-14];

   (* Anchor (2): g'(Pi_*) matches notes Sec. 3 quoted slope. *)
   expectApprox["g'(Pi_*) notes Sec. 3 value",
     N[D[gPi, piM] /. piM -> piStar, 50], 0.0714453558083195`50, 10^-14];

   (* Anchor (3): parent threshold identity at piM = piStar — notes Sec. 2. *)
   thresholdAtStar = FullSimplify[
     thresholdResidual /. piM -> piStar,
     Assumptions -> $Assumptions
   ];
   expectedForm = (tM - qStar*a0Prime) - piStar*thetaSigma/lM;
   If[TrueQ[Simplify[thresholdAtStar - expectedForm] === 0],
     pass["parent threshold identity at Pi = Pi_* (notes Sec. 2)"],
     fail["parent threshold identity at Pi = Pi_* (notes Sec. 2)",
          Simplify[thresholdAtStar - expectedForm]]
   ];

   (* Anchor (4): lower-branch discrimination — gPi(2*piStar) is clearly far from gMinus. *)
   offStarResidual = Abs[N[(gPi /. piM -> 2*piStar) - gMinus, 30]];
   If[TrueQ[offStarResidual > 10^-3],
     pass["lower-branch discrimination (paper Checks item 3)"],
     fail["lower-branch discrimination (paper Checks item 3)", offStarResidual]
   ];
   ```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 131` and confirm: (a) script exits 0; (b) banner reads `STAGE 131`; (c) transcript contains four `PASS:` lines for the four anchors; (d) NO `PASS: Pi_* compensation point` line (it was the tautological one).

---

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py` (anchored to F1 block above)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl` (anchored to F2 block above)

**Issue:**
Neither script substitutes `Pi -> Pi_*` into the "Parent bias mismatch formula" or distinguishes the lower branch from the singular branch (paper Checks item 3). Both gaps are closed by Anchors (3) and (4) prescribed in F1 and F2.

**Required change:**
The F1 and F2 blocks above already contain the Anchor (3) and Anchor (4) checks that resolve this finding. No additional edits required beyond F1 and F2.

**Verification command:**
Same as F1 and F2 — the verifier confirms the Anchor (3) and Anchor (4) PASS lines in both transcripts.

---

## F4 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:12`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:33`

**Issue:**
`g_minus` is imported as a bare numeric literal (SymPy) or as an exact closed form (Mathematica) with no in-script citation to the upstream stage that derives the F1 lower-branch value. The two engines also use different representations (float vs exact), weakening independence.

**Required change:**

1. SymPy line 12. Replace:
   ```python
   g_minus = sp.Float("0.758035078944663")
   ```
   with:
   ```python
   # F1 lower-branch value g_-^{F1}: closed form is (2*sqrt(4107 - 100*pi^2) - 37*sqrt(3)) / (20*pi).
   # Carried forward from the F1 family lower-branch derivation in the upstream
   # parent-mouth stages (see notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md Sec. 3).
   g_minus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) - 37*sp.sqrt(3)) / (20*sp.pi)
   g_minus_literal = sp.Float("0.758035078944663", 50)
   assert abs(sp.N(g_minus_exact, 50) - g_minus_literal) < sp.Float("1e-14", 50), (
       f"g_-^{{F1}} closed form vs literal: {sp.N(g_minus_exact, 50)} != {g_minus_literal}"
   )
   g_minus = sp.N(g_minus_exact, 50)
   print(f"PASS: g_-^F1 closed form matches literal 0.758035078944663")
   ```

2. Mathematica line 33. Replace:
   ```
   gMinus = N[(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi), 80];
   ```
   with:
   ```mathematica
   (* F1 lower-branch value g_-^{F1}: closed form below; see notes/stages/
      moving_throat_pde_stage131_parent_mouth_threshold.md Sec. 3 for context. *)
   gMinusExact = (2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi);
   gMinusLiteral = 0.758035078944663`50;
   expectApprox["g_-^F1 closed form vs literal",
     N[gMinusExact, 50], gMinusLiteral, 10^-14];
   gMinus = N[gMinusExact, 80];
   ```

**Verification command:**
After Codex applies, the verifier confirms: (a) both transcripts contain a PASS line tying the closed form `(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi)` to the numeric `0.758035078944663`; (b) a grep for `4107` in both `.py` and `.wl` returns a match; (c) both scripts exit 0.
