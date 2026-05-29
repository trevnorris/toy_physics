---
unit_id: 139
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 139

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:53`

**Issue:**
The scripts compute six headline numerical values (`r_F1`, `R_q^nat`, `M_s^nat,*`, `M_q^nat,*`, `M_s^comp,*`, `M_q^comp,*`) but only assert one identity (`R_q^comp = 1/4`), which holds algebraically for any `rF` and so does not exercise Family-1 specific content. None of the paper card's three `Checks` items are tested.

**Required change:**

In the SymPy script (`scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`), REPLACE the existing line 30 with the following assertion block (insert immediately AFTER the existing `print('mixed gain magnitude ratio = ...')` on line 28; remove the lone existing `assert` on line 30 and replace with the full block):

```python
# Numerical-deliverable assertions against the boxed values in
# notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md.
tol_literal = sp.Float('1e-12')   # notes give ~15 digits
tol_algebraic = sp.Float('1e-25') # algebraic identities

# r_F1 closed form vs notes literal
assert abs(rF - sp.Float('1.77799353547498', 20)) < tol_literal, (rF,)

# R_q^nat closed form vs notes literal
assert abs(Rq_nat - sp.Float('0.145454452260421', 20)) < tol_literal, (Rq_nat,)

# Natural-branch mouth gains vs notes literals
assert abs(Ms_nat - sp.Float('1.66854252965624', 20)) < tol_literal, (Ms_nat,)
assert abs(Mq_nat - sp.Float('-0.242696939724365', 20)) < tol_literal, (Mq_nat,)

# Compensated-branch mouth gains vs notes literals
assert abs(Ms_comp - sp.Float('1.80594111095636', 20)) < tol_literal, (Ms_comp,)
assert abs(Mq_comp - sp.Float('-0.451485277739090', 20)) < tol_literal, (Mq_comp,)

# Outlet consistency Pi_* = M_s + M_q * S_q(Pi_*), both branches.
# (Algebraically built in by M_s = Pi_*/(1 - R_q S_q), M_q = -R_q M_s, but
# verifies that the imported Pi_* and S_q literals satisfy the form.)
assert abs(Pi_star - (Ms_nat + Mq_nat * Sq_star)) < tol_algebraic
assert abs(Pi_star - (Ms_comp + Mq_comp * Sq_star)) < tol_algebraic

# Compensated R_q closed form
assert abs(Rq_comp - sp.Rational(1, 4)) < tol_algebraic

print('all assertions passed')
```

In the Mathematica script (`mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`), REPLACE line 53 (`expectApprox["R_q^comp - 1/4", rQComp, 1/4, 10^-25];`) with the following block (keeping the `Print[""]`, `Print["Stage 139 Mathematica audit passed."]`, `Exit[0]` lines that follow):

```mathematica
(* Numerical-deliverable assertions against the boxed values in
   notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md. *)
tolLit = 10^-12;
tolAlg = 10^-25;

expectApprox["r_F1", rF, 1.77799353547498, tolLit];
expectApprox["R_q^nat", rQNat, 0.145454452260421, tolLit];
expectApprox["M_s^nat,*", mSNat, 1.66854252965624, tolLit];
expectApprox["M_q^nat,*", mQNat, -0.242696939724365, tolLit];
expectApprox["M_s^comp,*", mSComp, 1.80594111095636, tolLit];
expectApprox["M_q^comp,*", mQComp, -0.451485277739090, tolLit];
expectApprox["outlet consistency nat", piStar, mSNat + mQNat*sQStar, tolAlg];
expectApprox["outlet consistency comp", piStar, mSComp + mQComp*sQStar, tolAlg];
expectApprox["R_q^comp - 1/4", rQComp, 1/4, tolAlg];
```

**Variable-independence self-check:** All assertions are on numerical scalars (no `diff` / `D[]` involved), so no zero-derivative trap. The outlet-consistency residuals reduce to `0` exactly by construction (the `M_s` and `M_q` definitions encode the same algebraic relation), giving `tol_algebraic` headroom. All other residuals reduce to `(computed_value - notes_literal)`, expected to be of order `10^-15` given that `Pi_star` and `Sq_star` are themselves only 15-digit literals.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:36`

**Issue:**
The Mathematica script assigns `gMinus = N[rF - Sqrt[1 + rF^2]/2, 30]` identically to the SymPy script's `g_minus = sp.N(rF - sp.sqrt(1 + rF**2)/2, 30)`. It does not derive `g_minus` from the compensated-branch defining condition; it just rewrites the SymPy algebra.

**Required change:**

REPLACE line 36 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl` (`gMinus = N[rF - Sqrt[1 + rF^2]/2, 30];`) with:

```mathematica
(* Derive g_c on the compensated branch by solving the defining condition
   (g_c - r_F1)^2 / (1 + r_F1^2) == 1/4, branch g_c < r_F1 (lower
   compensated branch, see notes/stages/moving_throat_pde_stage139*.md
   section "Exact compensated branch"). *)
gMinusSolutions = gc /. Solve[(gc - rF)^2 == (1 + rF^2)/4 && gc < rF, gc, Reals];
If[Length[gMinusSolutions] =!= 1,
  Print["FAIL: g_minus branch selection ambiguous, solutions = ", gMinusSolutions];
  Exit[1];
];
gMinus = N[First[gMinusSolutions], 30];
(* Cross-check the derived value against the closed form quoted in the notes. *)
expectApprox["g_minus closed form", gMinus, rF - Sqrt[1 + rF^2]/2, 10^-25];
```

This makes the Mathematica side an independent re-derivation (via `Solve`) that recovers the same `g_minus` and then cross-checks against the closed form. The SymPy script is unchanged for this finding.

**Variable-independence self-check:** `Solve[(gc - rF)^2 == (1 + rF^2)/4, gc, Reals]` returns two real roots `gc = rF ± Sqrt[1 + rF^2]/2`. The constraint `gc < rF` selects the minus branch, matching the notes. No symbolic-versus-numeric trap because `rF` is already numeric at this point (line 28 evaluated `rF` to 30 digits).

## F3 — hardcoded_result

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:6-7`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:29-30`

**Issue:**
`Pi_star` and `Sq_star` (sympy lines 6-7) and `piStar`, `sQStar` (mathematica lines 29-30) are bare numeric literals with no provenance comments. The notes attribute them to Stage 236.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`, REPLACE lines 5-7:
```python
rF = sp.N(sp.sqrt(sp.Rational(12,1)/sp.pi**2 * (sp.Rational(37,20))**2 - 1), 30)
Pi_star = sp.N('1.50882951349316', 30)
Sq_star = sp.N('0.658075937605429', 30)
```
with:
```python
# r_F1 closed form from Stage 223.
rF = sp.N(sp.sqrt(sp.Rational(12,1)/sp.pi**2 * (sp.Rational(37,20))**2 - 1), 30)
# Pi_* and S_q(Pi_*) imported from Stage 236.
Pi_star = sp.N('1.50882951349316', 30)
Sq_star = sp.N('0.658075937605429', 30)
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`, REPLACE lines 28-30:
```mathematica
rF = N[Sqrt[(12/Pi^2)*(37/20)^2 - 1], 30];
piStar = SetPrecision[1.50882951349316, 30];
sQStar = SetPrecision[0.658075937605429, 30];
```
with:
```mathematica
(* r_F1 closed form from Stage 223. *)
rF = N[Sqrt[(12/Pi^2)*(37/20)^2 - 1], 30];
(* Pi_* and S_q(Pi_*) imported from Stage 236. *)
piStar = SetPrecision[1.50882951349316, 30];
sQStar = SetPrecision[0.658075937605429, 30];
```

No assertion change here — F3 is purely a provenance comment fix.

## F4 — stale_output (banner typo subtype)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:26`

**Issue:**
The Mathematica banner reads `"STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS"` but this is the Stage 139 audit. The output transcript propagates the misleading label.

**Required change:**

REPLACE line 26:
```mathematica
banner["STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS"];
```
with:
```mathematica
banner["STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS"];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 139` and `redteam exec-mathematica 139` (one at a time per the single-seat policy) and confirm:
- both scripts exit 0
- the sympy transcript contains the line `all assertions passed`
- the mathematica transcript banner line reads `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`
- the mathematica transcript contains at least 9 `PASS:` lines (the new `expectApprox` block from F1+F2)
