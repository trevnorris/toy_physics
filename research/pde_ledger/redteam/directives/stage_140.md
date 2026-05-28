---
unit_id: 140
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 140

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:26`

**Issue:** The Mathematica banner is a copy-paste artifact from another stage. Line 26 currently reads `banner["STAGE 123 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"];` but this is stage 140. The underlying math is correct; only the printed banner string is wrong.

**Required change:**
Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl` line 26.

Before:
```
banner["STAGE 123 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"];
```

After:
```
banner["STAGE 140 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"];
```

Do not modify any other line in this file.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 140` and confirm the printed banner reads "STAGE 140 — ..." and the `expectZero["Sigma_0_hat - 20 That^2/9", ...]` residual remains 0 with the script exiting 0.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py:18-24`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:46-53`

**Issue:** The Section-3 numerics from the notes (`That_nat ≈ 0.866512630228382`, `That_comp ≈ 0.901484054174206`, fractional enhancement ≈ 0.0403588161624) are reproduced by `print` only. There is no `assert` / `expectZero` guarding them, so a regression on the input literals or the `sqrt(9*Ms/20)` formula would silently pass.

**Required change:**

(1) In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py`, after the existing line 24 (which prints `fractional traction enhancement = ...`), append three assertions that compare the computed values to the notes' boxed numerics with a 1e-12 tolerance.

Before (lines 18-24, unchanged):
```
Ms_nat = sp.N('1.6685425296562397', 30)
Ms_comp = sp.N('1.80594111095636', 30)
That_nat = sp.N(sp.sqrt(9*Ms_nat/20), 30)
That_comp = sp.N(sp.sqrt(9*Ms_comp/20), 30)
print('That_nat =', That_nat)
print('That_comp =', That_comp)
print('fractional traction enhancement =', sp.N(That_comp/That_nat - 1, 20))
```

After (append immediately after line 24, with no other changes to lines 1-24):
```
assert sp.Abs(That_nat - sp.Float('0.866512630228382', 30)) < sp.Float('1e-12', 30)
assert sp.Abs(That_comp - sp.Float('0.901484054174206', 30)) < sp.Float('1e-12', 30)
assert sp.Abs(sp.N(That_comp/That_nat - 1, 30) - sp.Float('0.0403588161624', 30)) < sp.Float('1e-12', 30)
print('numeric fixed-point checks PASS')
```

(2) In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl`, after the existing `Print["fractional traction enhancement = ", ...]` (line 53), and before the `Print[""]` at line 55, insert three numeric checks using the existing `expectZero` helper. Because `expectZero` calls `FullSimplify === 0`, wrap each diff in `If[Abs[...] < tol, 0, ...]` so a tiny floating-point residual still resolves to literal 0.

Before (lines 53-55, unchanged):
```
Print["fractional traction enhancement = ", fmt[N[tHatComp/tHatNat - 1, 20]]];

Print[""];
```

After (insert between line 53 and line 55):
```
Print["fractional traction enhancement = ", fmt[N[tHatComp/tHatNat - 1, 20]]];

Module[{diff1, diff2, diff3, tol},
  tol = 10^-12;
  diff1 = N[tHatNat - SetPrecision[0.866512630228382, 30], 30];
  diff2 = N[tHatComp - SetPrecision[0.901484054174206, 30], 30];
  diff3 = N[(tHatComp/tHatNat - 1) - SetPrecision[0.0403588161624, 30], 30];
  If[Abs[diff1] < tol, pass["That_nat matches notes"], fail["That_nat matches notes", diff1]];
  If[Abs[diff2] < tol, pass["That_comp matches notes"], fail["That_comp matches notes", diff2]];
  If[Abs[diff3] < tol, pass["fractional enhancement matches notes"], fail["fractional enhancement matches notes", diff3]];
];

Print[""];
```

Do not modify the existing `Print["Stage 140 Mathematica audit passed."]` or `Exit[0]` lines.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 140` and `redteam exec-mathematica 140`. Both scripts must still exit 0. The sympy transcript must contain the new line `numeric fixed-point checks PASS`. The Mathematica transcript must contain three new `PASS:` lines for "That_nat matches notes", "That_comp matches notes", and "fractional enhancement matches notes".
