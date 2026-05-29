---
unit_id: 123
batch: IV.3
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 123

Apply each non-`paper_misalignment` finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper, notes, or scripts to "fix" a `paper_misalignment` unless the user has explicitly chosen a direction in a follow-up directive.

If a non-`paper_misalignment` finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** value_mismatch

**Paper side (notes):**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md:25-31` quote:
  > `\Xi_v = -\frac{3\sqrt{30}\,\pi^{3/2}}{228}\,\mathfrak r.`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage123_parent_normalized_branch_values.md:37-41` quote:
  > `\Xi_v^{F1} = -\frac{3\sqrt{30}\,\pi^{3/2}}{228}\,\mathfrak r_{F1} \approx -1.01675633282526.`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:46` quote:
  > `expect_zero("Xi_v law", Xi_v_expr + 3*sp.sqrt(30)*sp.pi**sp.Rational(3,2)*r/160)`
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.txt:13` quote:
  > `Xi_v(r) = -3*sqrt(30)*pi**(3/2)*r/160`

## Resolve before fix_loop

The notes' boxed symbolic constant is `228` but the same notes' boxed numeric `Xi_v^{F1} ≈ -1.01675633282526` is only consistent with denominator `160`. The script's independent re-derivation from `Ks`, `Kq`, `lam` yields `160` and the numeric `-1.01675633282526...` to 20 digits. Which is correct?

- Manual check: with `r_{F1} = sqrt(12*(37/20)^2/pi^2 - 1) ≈ 1.7779`,
  - `-3*sqrt(30)*pi^(3/2)*r_{F1}/160 ≈ -1.0168` (matches the boxed numeric)
  - `-3*sqrt(30)*pi^(3/2)*r_{F1}/228 ≈ -0.713`  (does NOT match the boxed numeric)

So the notes are internally inconsistent. Almost certainly the `228` is a typo for `160`.

Possible directions (the user picks one):
- (a) `160` is correct (recommended; matches both the derivation and the notes' own numeric). Update the notes file (lines 25-31 and 37-41) replacing `228` with `160`. No script change. Re-render the paper if needed.
- (b) `228` is correct. Then BOTH the script (line 46) AND the notes' numerical box `-1.01675633282526` are wrong; the upstream `lam` definition or one of `Ks`, `Kq` would need to be re-derived from stages 220-221 to produce `228`. This would also cascade into the numeric Xi_v^F1 quoted by the notes and likely into downstream stages 125-139.
- (c) Both are derived from a third source that contradicts both → flag for deeper review.

The orchestrator will not invoke Codex on this unit (for F1) until the user has chosen a direction. F2 below is independent and can be applied immediately.

## F2 — paper_misalignment (low; cosmetic banner)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:16`

**Issue:** The script banner reads `"STAGE 106 — PARENT-NORMALIZED BRANCH VALUES"`, but the file name, the paper-card `\label{stage:123}`, and the script's role in the ledger all identify this unit as stage 123 (paper section header "Stage~140"). The `106` label is dead — it does not match any current convention. Cosmetic, but the captured output file mirrors the wrong banner, which will confuse later searches.

**Required change:**

Edit line 16 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py`.

Before:

```python
banner("STAGE 106 — PARENT-NORMALIZED BRANCH VALUES")
```

After:

```python
banner("STAGE 123 — PARENT-NORMALIZED BRANCH VALUES")
```

No other lines change. Do not adjust assertions, do not alter symbolic content. This is a string-only edit.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 123` and confirm:
- the new sympy output's banner line reads `STAGE 123 — PARENT-NORMALIZED BRANCH VALUES`,
- exit code remains 0,
- the four printed numeric values (`Xi_v(F1)`, `Xi_T(nat)`, `Xi_T(-)`, `Xi_T(+)`) are bit-identical to the previous output,
- both assertions (`Xi_v law = 0`, `Xi_T law = 0`) still hold.
