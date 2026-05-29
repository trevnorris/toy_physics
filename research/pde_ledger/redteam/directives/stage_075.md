---
unit_id: 075
batch: III.4
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 075

The original F2 `paper_misalignment` (Upsilon_w conversion factor) was resolved by the user as direction (a) on 2026-05-27: paper-side and notes-side text update `117/168 -> 100`, no script change to the `alpha_r = 10` value. The orchestrator has already applied the paper-side and notes-side edits to `paper/stages/stage_075.tex`, `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md`. Codex now applies F1, F3, and the new F4 (script-side `alpha_r**2 == 100` lock per the resolved direction).

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F2 — paper_misalignment

**Subtype:** value_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_075.tex:7` quote: `\stagefield{Inputs}{\(\kappa=12321/5\), \(\eta=37\), and \(\Upsilon_w=117\Theta_w\).}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_075.tex:27-29` quote (boxed Output): `\Theta_w\leq3.62605617972939\times10^{-4}\,\mathrm{Pe}_{\rm req}\Rightarrow\text{fail}` — this number requires `Upsilon_w / Theta_w = 100`, not 117.
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:92` quote: `alpha_r = 10.`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:108` quote: `Upsilon_w = 168 Theta_w.` — inconsistent with both `alpha_r = 10` and with the notes' own subsequent Theta_w numerics on lines 124-128.
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:124-128` quote: `Theta_fail = Upsilon_fail / 168 ≈ 3.62605617972939e-4 * Pe_req` — but `0.0362605617972939 / 168 = 2.158e-4`, not `3.626e-4`. The text says "divide by 168" but the printed result is `divide by 100`.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:24` quote: `alpha_r = sp.Integer(10)` → `Upsilon_w = alpha_r^2 Theta_w = 100 Theta_w`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:37` quote: `alphaR = 10;` → same.
- Script's final numeric output (sympy line 35 / mathematica line 17): `Theta_fail / Pe_req = 0.00036260561797293886969` — matches the paper's boxed value `3.62605617972939e-4` exactly, confirming the script and the paper's BOXED output both agree on `Upsilon_w = 100 Theta_w`.

## Approved by user (2026-05-27)

- direction: (a) — paper Inputs `117` and notes section 3 `168` are both stale text drift; the script's `alpha_r = 10 -> 100` is the load-bearing value.
- applied_by: orchestrator on 2026-05-27 (paper and notes side):
  - `paper/stages/stage_075.tex:7` Inputs line `117 Theta_w` -> `100 Theta_w`
  - `paper/stages/stage_075.tex:24` body line `117 Theta_w` -> `100 Theta_w`
  - `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:108` `168 Theta_w` -> `100 Theta_w`
  - `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:116` `168 Theta_w` -> `100 Theta_w`
  - `notes/stages/moving_throat_pde_stage075_family1_threshold_window.md:124-128` `Upsilon_fail / 168` and `Upsilon_suff / 168` -> `/100`
- Codex now applies F4 (script-side assertion locking `alpha_r**2 == 100`) below, in addition to F1 and F3.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:40-61,99-104`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:75-96`

**Issue:**
The v1 audit's F1 and F2 patches both shifted the tautology rather than discharging it. Concretely:

1. The "round-trip" `Upsilon_fail - alpha_r^2 * Theta_fail == 0` (sympy line 103; mathematica line 95) is tautological because `Theta_fail` is defined as `Upsilon_fail / alpha_r^2` (sympy line 78; mathematica line 65). Substituting gives `Upsilon_fail - alpha_r^2 * (Upsilon_fail / alpha_r^2) = Upsilon_fail - Upsilon_fail = 0` regardless of any physical content.

2. The "free-symbol identity" check (sympy lines 50-61; mathematica lines 87-93) `(alpha*sinh+eta*cosh) * Delta_0 - eta*(cosh-1)/alpha^2 == 0` is the *definition* of `Delta_0` rearranged: substituting `Delta_0 := eta*(cosh-1) / (alpha^2 * (alpha*sinh+eta*cosh))` reduces the LHS to `eta*(cosh-1)/alpha^2 - eta*(cosh-1)/alpha^2 = 0` by canceling the `(alpha*sinh+eta*cosh)` factor. This is a CAS-level algebraic cancellation, not a non-trivial property of `Delta_0`.

**Required change:**

Add a non-trivial independent check for both `Delta_0` and `Delta_inf` that exercises the closed form against an asymptotic limit. Concretely:

In the SymPy script, after line 61 (after the existing identity assertions), insert:
```
# Asymptotic check: large-alpha limit of Delta_inf is 1/alpha,
# and Delta_0 is exponentially suppressed (specifically alpha * Delta_inf → 1).
# This is a non-trivial consequence of the closed form: if the numerator or
# denominator of Delta_inf had a wrong factor, the limit would not be 1.
large_alpha_check = sp.limit(alpha_sym * Deltainf_sym, alpha_sym, sp.oo)
print("alpha * Delta_inf large-alpha limit =", large_alpha_check)
assert large_alpha_check == 1

# Small-alpha leading order: Delta_0 → 1/2 as alpha → 0
# (since cosh(alpha)-1 ~ alpha^2/2 and alpha*sinh(alpha)+eta*cosh(alpha) ~ eta).
small_alpha_check_delta0 = sp.limit(Delta0_sym, alpha_sym, 0)
print("Delta_0 small-alpha limit =", small_alpha_check_delta0)
assert small_alpha_check_delta0 == sp.Rational(1, 2)
```

In the Mathematica script, after line 93 (closing of the existing identity `Module`), insert:
```
(* Asymptotic check: large-alpha limit of Delta_inf is 1/alpha,
   so alpha*Delta_inf → 1. This is a non-trivial consequence of
   the closed form and is computed by Mathematica's Limit
   independently from SymPy's sp.limit. *)
Module[{aSym, eSym, delta0Sym, deltaInfSym, largeAlphaLimit, smallAlphaLimit},
  ClearAll[aSym, eSym];
  delta0Sym = eSym*(Cosh[aSym] - 1)/(aSym^2*(aSym*Sinh[aSym] + eSym*Cosh[aSym]));
  deltaInfSym = (Cosh[aSym] + (eSym/aSym)*Sinh[aSym] - 1)/(aSym*Sinh[aSym] + eSym*Cosh[aSym]);
  largeAlphaLimit = Limit[aSym*deltaInfSym, aSym -> Infinity, Assumptions -> eSym > 0];
  Print["alpha * Delta_inf large-alpha limit = ", fmt[largeAlphaLimit]];
  If[TrueQ[largeAlphaLimit === 1], pass["alpha * Delta_inf → 1"], fail["alpha * Delta_inf → 1", largeAlphaLimit]];
  smallAlphaLimit = Limit[delta0Sym, aSym -> 0, Assumptions -> eSym > 0];
  Print["Delta_0 small-alpha limit = ", fmt[smallAlphaLimit]];
  If[TrueQ[smallAlphaLimit === 1/2], pass["Delta_0 → 1/2"], fail["Delta_0 → 1/2", smallAlphaLimit]];
];
```

Also rewrite the misleading existing comment block in the Mathematica script (lines 78-82) to refer accurately to the new asymptotic checks rather than describing the tautological identity as "the independent-derivation leg".

Replace lines 78-82 (the existing `(* This identity check is the independent-derivation leg... *)` comment) with:
```
(* Note: the algebraic identity below is structurally tautological (it follows
   from the definition of Delta_0 / Delta_inf by canceling a common factor).
   The genuine independent check is the asymptotic-limit block further below,
   which exercises a non-trivial property of the closed forms via Mathematica's
   Limit operator. *)
```

The existing round-trip checks (sympy 99-104; mathematica 95-96) may remain in place as informational consistency checks — they are tautological but harmless, and removing them is out of scope.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 075` and `redteam exec-mathematica 075` and confirm both transcripts now print:
- `alpha * Delta_inf large-alpha limit = 1` (or its Mathematica equivalent) followed by PASS,
- `Delta_0 small-alpha limit = 1/2` followed by PASS,

and that the limits were computed from the free-symbol form (`alpha_sym`/`aSym` not numerically substituted).

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:78-82`

**Issue:**
The v1 audit's F3 fix declared the free-symbol identity check (M1/M2) to be the "independent-derivation leg" via a comment at lines 78-82. As F1 above shows, that identity is itself tautological, so the comment is inaccurate and the transliteration concern is unresolved.

**Required change:**
F3 is structurally subsumed by F1: once F1 inserts the asymptotic-limit checks (which Mathematica's `Limit` evaluates independently from SymPy's `sp.limit`), the Mathematica script genuinely exercises an independent algorithm path. The comment edit specified in F1 (replacing lines 78-82) is the F3 fix.

No additional file edits beyond F1's comment replacement are required for F3.

**Verification command:**
After Codex applies F1, the verifier reads the Mathematica script and confirms (i) the asymptotic-limit `Module` block from F1 is present with `aSym`, `eSym` as free symbols, (ii) the comment at the former lines 78-82 no longer claims a tautological check is the independent leg, and (iii) the new comment accurately describes the asymptotic block.

## F4 — script-side assertion lock for `alpha_r^2 == 100`

This is a follow-up item from the F2 user resolution (direction (a)). Locking the value via a script-side assertion prevents future drift between the paper Inputs line and the script's `alpha_r` value.

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:24` (immediately after the `alpha_r = sp.Integer(10)` line)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:37` (immediately after the `alphaR = 10;` line)

**Required change (SymPy):**

After the existing `alpha_r = sp.Integer(10)` at line 24, insert:

```python
# Paper Inputs line and notes section 3 state Upsilon_w = alpha_r^2 Theta_w with alpha_r^2 = 100.
# Lock the value so any future drift in the paper Inputs line surfaces here.
expect_zero("alpha_r^2 - 100 (paper Inputs line lock)", alpha_r**2 - 100)
```

**Required change (Mathematica):**

After the existing `alphaR = 10;` at line 37, insert:

```mathematica
(* Paper Inputs line and notes section 3 state Upsilon_w = alpha_r^2 Theta_w with alpha_r^2 = 100.
   Lock the value so any future drift in the paper Inputs line surfaces here. *)
expectZero["alpha_r^2 - 100 (paper Inputs line lock)", alphaR^2 - 100];
```

**Verification:**
After Codex applies, both transcripts should print PASS for `alpha_r^2 - 100 (paper Inputs line lock)` with residual `0`. Both scripts still exit 0.
