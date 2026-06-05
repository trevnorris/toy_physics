---
unit_id: 060
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T12:05:58-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 060

Apply each finding below in order. After applying, append an `## Applied: F<n>` block with `files_changed`, `summary` (one sentence), `deviation` (or "none"). If a change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead and continue with the rest.

Edit only what each finding names. No unrelated refactors.
After editing, RUN both engines (`python3 <path>`, `math -script <path>`) and iterate until they exit 0 with all in-file checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-label; number-only, format preserved)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py:3`

**Issue:** This is stage 060 but the docstring carries the pre-renumber "Stage 43" (43+17=60). Committed `.txt` outputs are stale and refresh on re-run. The `.wl` self-labels are already canonical. Label only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3`: `Stage 43 SymPy audit — entropic source microclosure and microscopic Xi.` → `Stage 60 SymPy audit — entropic source microclosure and microscopic Xi.`
2. Re-run both engines to refresh outputs.

**DO NOT TOUCH (deferred to the dedicated SCRIPT/OUTPUT-band numbering plan — CROSS-reference to another stage, not a self-label):**
- `py:159` `the zero-flux branch is exactly the Stage-39 exponential family.` — "Stage-39" is a cross-reference to stage 056 (39+17=56), not a self-label. Leave it exactly as-is.

Do NOT pad the already-correct `STAGE 60` banners (py:22, py:156).

**Verification:** `redteam exec-sympy 060` exits 0; `py:3` reads `Stage 60 ...`; `py:159` still reads `Stage-39`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.txt`
- summary: Updated the SymPy self-label from Stage 43 to Stage 60 and refreshed both saved engine outputs.
- deviation: none

## F2 — tautological_check (Mathematica side; genuine fix)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl` — the assertion currently at `wl:140`, with `xiMicro` defined at `wl:132`.

**Issue:** The headline microscopic gain is *defined* at `wl:132`:
`xiMicro = lambdaPhi^2*ell^2/(theta*tX);`
and then "checked" at `wl:140` against a verbatim copy of that same definition:
`expectZero["Xi_micro - Lambda^2 L^2/(Theta T_X)", xiMicro - lambdaPhi^2*ell^2/(theta*tX)];`
The residual is `lambdaPhi^2*ell^2/(theta*tX) - lambdaPhi^2*ell^2/(theta*tX) === 0` by construction — it cannot fail and exercises nothing. (The adjacent `wl:134-138` chi/D-M consistency checks and `wl:141-142` susceptibility/phenomenological-form checks DO exercise the deliverable; line 140 is dead scaffolding.)

**Requirement (you design the route; do not just restate the definition):**
Make the closed form `Xi_micro = lambdaPhi^2*ell^2/(theta*tX)` exercised by a NON-tautological check — i.e. one whose left side is built from an object the script independently constructs (e.g. the susceptibility/Peclet/diffusion chain already present in this file or assembled from its derived parameters), not from a literal re-statement of `xiMicro`'s definition. If no genuinely independent construction is available beyond what `wl:134-142` already cover, the acceptable alternative is to DELETE `wl:140` (the deliverable stays covered by the surrounding checks). Either way, the line `xiMicro - lambdaPhi^2*ell^2/(theta*tX)` must not survive.

**Acceptance criteria:**
- No assertion in the file subtracts `xiMicro` from a verbatim copy of its own `wl:132` definition.
- The `Xi_micro` deliverable remains verified by at least one non-tautological assertion (whether a new independent-route check or the retained `wl:134-142` set).
- `math -script` exits 0; all checks PASS; the SymPy engine is untouched and the two engines still agree on `Xi_micro`.
- `material_change` stays false (no verified RESULT value moves — this only removes a tautology / adds an independent assertion over the already-correct value).

**Verification:** the orchestrator re-runs `redteam exec-mathematica 060` and inspects the diff: `wl` no longer contains `xiMicro - lambdaPhi^2*ell^2/(theta*tX)`; exits 0; `Xi_micro` value unchanged.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.txt`
- summary: Replaced the tautological Xi_micro literal-definition assertion with a chi-route versus D/M-route consistency check.
- deviation: none
