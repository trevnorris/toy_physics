---
unit_id: 081
batch: III.4
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T21:38:03Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 081

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:38-42`

**Issue:** The SymPy script verifies the stage's primary deliverable — the exact inversion `Q = [1 + (1-2 eps_blk) zeta]/[1 - eps_blk zeta]` — only at two interpolation points (`Q(0)-1`, `Q(1)-2`, lines 40-41). Two points do not pin a rational function; the full form is only guaranteed by `sp.solve`, never asserted. The Mathematica engine already asserts the full closed-form identity (`mathematica/...stage081...wl:54`). Add the matching full-identity assertion to SymPy so both engines load-bearingly verify the conversion formula.

**Required change:**
Insert a new `expect_zero` call immediately after the `Q(zeta;eps_blk)` print (current line 39), before the existing `Q(0)-1` anchor check (line 40). Insert exactly:

```python
expect_zero("Q-closedform",
    Q - (1 + zeta - 2*eps_blk*zeta)/(1 - eps_blk*zeta))
```

The symbols `Q`, `zeta`, and `eps_blk` are already in scope (declared line 30, `Q` formed line 36). Do not change any other line. The RHS is the verbatim notes/paper closed form (do NOT alter the constants `1`, `2`, or the `(1-2 eps_blk)` coefficient).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 081` and confirm a new output line `Q-closedform = 0` appears AND the script exits 0 with all checks passing.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py`
- summary: Added the SymPy full closed-form residual check for `Q(zeta;eps_blk)`.
- deviation: none

## F2 — DEFERRED by orchestrator (do NOT apply)

The audit flagged two comment-only labels at `:32` (`Stage-35`) and `:44` (`Stage 63`). These are **CROSS-stage references** to OTHER stages (052 and 080), not stage-081 self-labels. Per the settled Reading-2 in-loop policy, the loop fixes only a stage's own UNAMBIGUOUS self-labels; cross-references are owned by the dedicated `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed, never offset-swept) and are LEFT UNTOUCHED here. **Codex: do nothing for F2.** Only F1 (the `expect_zero` insertion) is in scope. The docstring line 3 (`SymPy audit for Stage 081.`) self-label is already canonical.
