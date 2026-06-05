---
unit_id: 043
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T13:58:02Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 043

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`

**Issue:** the SymPy source carries stale pre-renumber SELF-labels (pre-renumber self-number `26`; canonical = `043`) in the docstring header + the `26.k` subbanner indices, matching the file's canonical banner `STAGE 43`. The main banner was already fixed; these self-labels were missed and print INTO the transcript, so a plain re-run leaves them stale. (The committed `.txt` are also stale — the orchestrator re-run refreshes them. The `.wl` source labels are already canonical.)

**Required change (label-only — change ONLY the indicated numeric token, preserve the rest of each string verbatim):**
- line 3: `Moving-throat PDE — Stage 26 SymPy audit.` → `Moving-throat PDE — Stage 43 SymPy audit.`
- line 57: subbanner `26.1` → `43.1`
- line 83: subbanner `26.2` → `43.2`
- line 117: subbanner `26.3` → `43.3`
- line 150: subbanner `26.4` → `43.4`
- line 164: subbanner `26.5` → `43.5`

**DO NOT TOUCH:** the already-canonical `STAGE 43` banners (lines 55, 172) — LEAVE. No cross-refs in this `.py`. No `.wl` edit for F1 (its labels are already canonical).

**Verification command:**
The verifier confirms the `.py` docstring + subbanners read `43` / `43.1`–`43.5`, the regenerated sympy `.txt` line 3 reads `STAGE 43 — ...`, the mathematica `.txt` line 3 reads `STAGE 043 — ...`, every PASS line remains, and the strip-the-number diff is byte-identical.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
- summary: Updated the stale SymPy Stage 26 self-labels and subbanner indices to canonical Stage 43 labels.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:143-148`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:125-129`

**Issue:** The "baseline value identification" block substitutes the SAME literal `B = 8/pi^2` into both sides of an identity already proven at sympy line 141 / wl line 123 (`M_supp structural form (free baseline)` == 0). Subtracting two halves of a proven identity after an identical substitution cannot fail — the assertion (sympy 148, wl 129) is redundant and exercises nothing new. It also does NOT verify the stage's actual numeric claim that the support baseline `kappa0^2` equals `8/pi^2`; the literal is assumed, not derived. The notes derive it from the stage's own frozen constants: `kappa0^2 = (9/11) sigma` with `sigma = 88/(9 pi^2)`, giving `kappa0^2 = (9/11)*88/(9 pi^2) = 8/pi^2`.

**Required change:**
Keep the structural-form check (sympy 141 / wl 123) exactly as-is. Replace the redundant numeric-baseline assertion (sympy 144-148 / wl 126-129) with a check that actually derives the literal `8/pi^2` from the stage's overlap constants, so the number is exercised rather than assumed.

SymPy (replace lines 143-148):
```python
# F2: baseline value identification — derive B = kappa0^2 = (9/11) sigma from
# the frozen overlap sigma = 88/(9 pi^2), so the literal 8/pi^2 is PRODUCED, not assumed.
sigma_value = sp.Rational(88, 1) / (sp.Integer(9) * sp.pi**2)
B_value = sp.simplify(sp.Rational(9, 11) * sigma_value)
print(f"baseline B = kappa0^2 = (9/11) sigma = {B_value}")
expect_zero("baseline B = 8/pi^2 from frozen sigma", B_value - sp.Rational(8, 1) / sp.pi**2)
Msupp_cont_eval = sp.simplify(Msupp_cont_in_B.subs(B, B_value))
Msupp_expected = sp.simplify(Msupp_struct_expected.subs(B, B_value))
print("M_supp at baseline B = 8/pi^2 =")
sp.pprint(sp.factor(Msupp_cont_eval))
expect_zero("M_supp at baseline B = 8/pi^2", Msupp_cont_eval - Msupp_expected)
```
(The trailing `Msupp_cont_eval`/`Msupp_expected` evaluation may be kept as a print/sanity line; the load-bearing NEW check is `baseline B = 8/pi^2 from frozen sigma`.)

Mathematica (replace lines 125-129):
```wolfram
(* F2: baseline value identification — derive B = kappa0^2 = (9/11) sigma from
   the frozen overlap sigma = 88/(9 Pi^2), so the literal 8/Pi^2 is PRODUCED, not assumed. *)
sigmaValue = 88/(9 Pi^2);
bValue = FullSimplify[(9/11) sigmaValue, Assumptions -> $Assumptions];
Print["baseline B = kappa0^2 = (9/11) sigma = ", fmt[bValue]];
expectZero["baseline B = 8/Pi^2 from frozen sigma", bValue - 8/Pi^2];
mSuppContEval = FullSimplify[mSuppContInB /. bBaseline -> bValue, Assumptions -> $Assumptions];
mSuppExpected = FullSimplify[mSuppStructExpected /. bBaseline -> bValue, Assumptions -> $Assumptions];
Print["M_supp at baseline B = 8/Pi^2 = ", fmt[mSuppContEval]];
expectZero["M_supp at baseline B = 8/Pi^2", mSuppContEval - mSuppExpected];
```

**Self-test (done by auditor):** `(9/11)*88/(9 pi^2) = 88/(11 pi^2) = 8/pi^2`, so `baseline B = 8/pi^2 from frozen sigma` evaluates to a literal `0` (real, non-vacuous — it subtracts two concrete rationals-times-1/pi^2). This matches notes:149+151 exactly, so the fix introduces no new paper_misalignment. The structural-form check at line 141/123 is untouched.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 043` and `redteam exec-mathematica 043` and confirms a new check labeled `baseline B = 8/pi^2 from frozen sigma` (resp. `... 8/Pi^2 ...`) appears, prints `= 0`, the scripts exit 0, and all prior PASS lines remain.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- summary: Replaced the redundant baseline substitution check with sigma-derived derivations of `B = 8/pi^2` and `B = 8/Pi^2` before the existing baseline sanity checks.
- deviation: none
