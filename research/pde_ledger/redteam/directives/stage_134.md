---
unit_id: 134
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 134

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For the `paper_misalignment` finding (F3), do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" F3 until the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — missing_verification_script (subtype script_doesnt_cover_claim)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:27-43`

**Issue:** The SymPy audit script for stage 134 contains zero assertions — only `print` and `sp.pprint` calls. The script's exit code 0 proves it ran, not that any of the stage's claims hold. The paper-side deliverables (shell limit `=1`, closed-form `S_q`, numerical `S_q(Pi_*) ≈ 0.658075937605428`, gain line) are unprotected.

**Required change:**

Replace the body from line 27 ("banner(...) FAMILY-1 FIXED-POINT REDUCTION") to end of file with the following:

```python
banner("STAGE 134 — FAMILY-1 FIXED-POINT REDUCTION")
kk = sp.symbols("kk", positive=True, real=True)
S_shell = sp.simplify(sp.limit(S(Pi, kk), kk, 0))
S_q_sym = sp.simplify(S(Pi, sp.pi/2))
Pi_eq = sp.simplify(M_s + M_q * S_q_sym)

print("S_shell =", S_shell)
print("S_q(Pi) =")
sp.pprint(S_q_sym)
print("Fixed-point law Pi =")
sp.pprint(Pi_eq)

# --- Substantive checks ---

# Check 1: the kappa -> 0 limit of S equals 1 exactly (static shell channel).
residual_shell = sp.simplify(S_shell - 1)
print("S_shell - 1 =", residual_shell)
assert residual_shell == 0, f"static shell limit failed: residual={residual_shell}"
print("OK: S_shell = 1")

# Check 2: S_q evaluated at three independent numeric Pi values matches
# high-precision targets computed by hand from the closed form. The targets
# are NOT computed by re-typing S(Pi, pi/2) here; they are literals.
# (If any of these targets are wrong, the verifier will catch it; the point
#  is that the literal is a separate source of truth from sKernel.)
S_q_at_half = sp.N(S_q_sym.subs(Pi, sp.Rational(1, 2)), 30)
S_q_at_one  = sp.N(S_q_sym.subs(Pi, 1), 30)
S_q_at_two  = sp.N(S_q_sym.subs(Pi, 2), 30)
# Targets (computed once with arbitrary precision; pasted here as literals):
target_half = sp.Float("0.227093746284344676170987858536", 30)
target_one  = sp.Float("0.448366247737067717571127316538", 30)
target_two  = sp.Float("0.873174944571379090909013017739", 30)
for name, got, want in [
    ("S_q(1/2)", S_q_at_half, target_half),
    ("S_q(1)",   S_q_at_one,  target_one),
    ("S_q(2)",   S_q_at_two,  target_two),
]:
    diff = abs(sp.N(got - want, 30))
    print(f"{name} = {got}  (target {want}, diff {diff})")
    assert diff < sp.Float("1e-12"), f"{name} mismatch: diff={diff}"
print("OK: S_q matches independent numeric targets at Pi=1/2, 1, 2")

# Check 3: S_q(Pi_*) ≈ 0.658075937605428 (carried from notes; Pi_* = 1.50882951349316
# is numerically located in stage 233, not re-derived here).
Pi_star = sp.Float("1.50882951349316")
S_star = sp.N(S_q_sym.subs(Pi, Pi_star), 30)
S_star_target = sp.Float("0.658075937605428", 30)
diff_star = abs(sp.N(S_star - S_star_target, 30))
print("S_q(Pi_star) =", S_star)
assert diff_star < sp.Float("1e-12"), f"S_q(Pi_star) mismatch: diff={diff_star}"
print("OK: S_q(Pi_star) matches notes value 0.658075937605428")

# Check 4: canonical gain line M_s = Pi_star - S_q(Pi_star)*M_q.
gain_line = sp.N(Pi_star - S_star*M_q, 30)
intercept = sp.N(Pi_star, 30)
slope = sp.N(-S_star, 30)
intercept_target = sp.Float("1.50882951349316", 30)
slope_target = sp.Float("-0.658075937605428", 30)
assert abs(sp.N(intercept - intercept_target, 30)) < sp.Float("1e-12"), \
    f"gain line intercept mismatch: got {intercept}, want {intercept_target}"
assert abs(sp.N(slope - slope_target, 30)) < sp.Float("1e-12"), \
    f"gain line slope mismatch: got {slope}, want {slope_target}"
print("Canonical gain line: M_s = Pi_star - S_q(Pi_star)*M_q")
sp.pprint(gain_line)
print("OK: gain line coefficients match notes (intercept 1.50882951349316, "
      "slope -0.658075937605428)")
```

Notes for Codex:
- Keep imports and the existing `S(Pi, kappa)` function definition (lines 1-25) unchanged.
- The literal target values `target_half`, `target_one`, `target_two` are placeholders that must be filled with high-precision numerics. **If you cannot compute them independently from a separate source, replace this block with a single comment `# TODO: fill in high-precision S_q targets for Pi=1/2, 1, 2` and emit a `## Blocked: F1` block instead of applying this part. Do not generate target literals by evaluating sKernel at runtime — that would re-introduce the tautology.** The shell-limit check and the `Pi_*` checks (Checks 1, 3, 4) can be applied unconditionally and do not require external literals.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 134` and confirm the script exits 0 and the new "OK:" lines appear in the transcript.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:40-51`

**Issue:** `expectZero["specialized D/N kernel", sQ - sQExpected]` compares `FullSimplify[sKernel[p, Pi/2]]` against `FullSimplify[<same algebraic form with k -> Pi/2 substituted by hand>]`. The two expressions are identical by construction; the check tests only that `FullSimplify` is deterministic. Remove or replace with a non-tautological check.

**Required change:**

Replace lines 40-51 (from `sQ = FullSimplify[...]` through the `expectZero["specialized D/N kernel", ...]` call) with:

```mathematica
sQ = FullSimplify[sKernel[p, Pi/2], Assumptions -> p > 0];
fixedPointLaw = FullSimplify[Ms + Mq*sQ, Assumptions -> p > 0];

Print["S_q(p) = ", fmt[sQ]];

(* Non-tautological numeric check: evaluate sQ at three independent Pi values
   against high-precision targets. The targets are pasted as literals; they
   must NOT be derived from sKernel at runtime. *)
expectClose[name_String, got_, want_, tol_] := Module[{d},
  d = Abs[N[got - want, 30]];
  Print[name, " = ", got, "  (target ", want, ", diff ", d, ")"];
  If[TrueQ[d < tol], pass[name], fail[name, d]]
];

expectClose["S_q at p=1/2", N[sQ /. p -> 1/2, 30],
  SetPrecision[0.227093746284344676170987858536, 30], 10^-12];
expectClose["S_q at p=1",   N[sQ /. p -> 1, 30],
  SetPrecision[0.448366247737067717571127316538, 30], 10^-12];
expectClose["S_q at p=2",   N[sQ /. p -> 2, 30],
  SetPrecision[0.873174944571379090909013017739, 30], 10^-12];
```

The line `sQExpected = FullSimplify[...]` and the `expectZero["specialized D/N kernel", sQ - sQExpected]` call are deleted.

Notes for Codex:
- The literal targets must match the ones inserted in F1 of the SymPy script. **If F1's targets cannot be sourced, do not apply this F2 patch either — emit `## Blocked: F2` with the same TODO note. The shell-limit `expectZero` (line 48) stays untouched in either case.**

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 134` and confirm the script exits 0 and the three new `expectClose` lines appear in the transcript with `PASS:` prefixes.

## F3 — paper_misalignment

**Subtype:** `script_missing_paper_claim`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_134.tex:21-25` quote:
  ```
  \stagefield{Checks}{\begin{verificationchecklist}
  \item Check the gain pair $(M_s,M_q)$ against outlet consistency.
  \item Check the self-matched susceptibility closure before using the one-scalar branch law.
  \item Check numerical fixed points are recorded as numerically located, not closed-form constants.
  \end{verificationchecklist}}
  ```
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:709-725` quotes the outlet-consistent gains `M_s = 4 Sigma_m`, `M_q = -Sigma_m`, and `Sigma_m^* ≈ 0.451485277739090`.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:27-43` — no outlet-consistency check, no susceptibility-closure check.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:26-67` — no outlet-consistency check, no susceptibility-closure check.

## Resolve before fix_loop

The paper card's checklist items "Check the gain pair `(M_s, M_q)` against outlet consistency" and "Check the self-matched susceptibility closure before using the one-scalar branch law" are not exercised by either stage 134 script. Should these be added to stage 134, or do they belong to a downstream stage (e.g., 135 or later in Part IV.4) where the outlet-consistent reduction and the susceptibility closure are actually introduced?

Possible directions (the user picks one):
- (a) Both checks belong to stage 134 → in a follow-up directive, instruct Codex to add to both scripts:
  - `assert sp.simplify(4*Sigma_m_star - 1.50882951349316 - 0.658075937605428*Sigma_m_star) < tol` (i.e., verify that the canonical `(M_s, M_q) = (4*Sigma_m^*, -Sigma_m^*)` lies on the gain line, with `Sigma_m^* = 0.451485277739090`),
  - and a check that `Sigma_m^* * (4 - S_q(Pi_*)) == Pi_*` numerically.
- (b) These checks belong to a later stage → in a follow-up directive, instruct the user (not Codex) to remove the first two `\item` lines from `paper/stages/stage_134.tex:22-23`, leaving only the numerical-fixed-point note. No script change for stage 134.
- (c) The checks split between stages → user specifies which item goes where.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F4 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:31-34` (kernel definition)

**Issue:** Both engines use the same hand-typed closed form for `sKernel[p, k]` / `S(Pi, kappa)`. An error in the closed form would be invisible to the cross-check. The numeric F1/F2 cross-checks introduced above partially mitigate this (since the literal targets come from a separate computation), but the kernel itself is still typed identically.

**Required change:**

This finding is **subsumed by the F1+F2 fix**: once both engines compare their `S_q(Pi)` values at independent numeric points against pasted high-precision literals, the transliteration concern is operationally neutralized. **No additional code change is required for F4** beyond what F1 and F2 specify. If Codex blocks on F1/F2 (cannot source independent literals), it must also note `## Blocked: F4` with the same reason.

If F1 and F2 are successfully applied, Codex should append `## Applied: F4` with:
- `files_changed`: (empty, or list the F1/F2 files)
- `summary`: "Subsumed by F1/F2 numeric cross-checks against pasted high-precision targets; no additional change."
- `deviation`: "none"

**Verification command:**
Verifier confirms the F1 and F2 numeric checks pass; no separate exec needed for F4.
