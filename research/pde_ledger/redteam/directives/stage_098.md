---
unit_id: 098
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 098

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F1), do nothing — the orchestrator is holding for user resolution. Do not edit `paper/stages/stage_098.tex`, `notes/stages/...md`, or the script banners/docstrings until the user has picked the canonical stage number.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch `paper/`, `notes/`, or any prose documents.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script (+ paper card internal inconsistency)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex:1` quote: `\section[Stage~115]{Stage~115: The Explicit Family-1 Support Test Is Automatic on the Actual Isotropic Branch}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex:2` quote: `\label{stage:098}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex:7` quote: `Stage~098 is a geometry-lane firewall ledger step.`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md:105` quote: `After Stages 80–81, the reduced theorem split is now fully sharp:`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:3` quote: `Stage 81 SymPy audit: actual isotropic branch support demand is automatic for any explicit family with zeta_max > 1 on the admissible blocked interval.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:35` quote: `print('\nSTAGE 81 AUDIT PASSED')`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:38` quote: `banner["STAGE 081 — FAMILY-1 SUPPORT IS AUTOMATIC"];`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:71` quote: `Print["Stage 098 Mathematica audit passed."];`

## Resolve before fix_loop

The stage number is inconsistent across four labels: the paper card section title says `Stage 115`, the LaTeX label and body say `Stage 098`, both scripts' banners/docstrings say `Stage 81 / 081` (and the notes section 4 references "Stages 80–81"), and the Mathematica final print says `Stage 098`. The math content is the same everywhere; only the numbering disagrees. Which is canonical?

Possible directions (the user picks one):
- (a) `098` is canonical (matches `\label`, file path, audit unit id, manifest entry) → update paper card section title from `\section[Stage~115]{Stage~115: ...}` to `\section[Stage~098]{Stage~098: ...}`; update sympy script docstring and final `print` to say `STAGE 098`; update Mathematica banner to `STAGE 098 — FAMILY-1 SUPPORT IS AUTOMATIC`; update notes section 4's "Stages 80–81" reference to reflect the new numbering (or leave it as historical if the renumbering is internal-only). The paper-side edits are blocked from Codex by policy; only after user confirmation will a follow-up directive authorize them.
- (b) `115` is canonical (the section title is correct and the rest is stale) → update label, body, file path, manifest, scripts, and notes to `115`. This is a much larger change touching the audit pipeline's unit id.
- (c) `081` is canonical (the original numbering before renumbering) → revert label, body, file path, scripts, notes back to `081`.
- (d) The numbers refer to different things (e.g., `115` is a paper-ordering index, `098` is the manifest unit id, `081` is the source-narrative stage) → add a clarifying note in the card that distinguishes the three numbering schemes, no script change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. F2 and F3 below are independent script-side fixes; they can be applied regardless of F1's resolution.

## F2 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:23-24` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:60-61`

**Issue:** The numeric value `zeta_max^(F1) = 2.46752922945601` is hardcoded as a literal in both scripts with no in-script derivation and no inline comment naming the upstream source. The notes file states the value once without citing a precursor stage. If the value originated from a typo or a stale precursor calculation, both scripts would still pass — `zeta_edge_F1` and `gap_F1` are derived from the same literal. Add a single-line inline comment at each hardcoded site identifying the value as an external carry-forward.

**Required change:**

Edit `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py`:

Before line 23 (`# Family-1 specialization.`) leave the existing comment; replace line 24 from:

```
zmax_F1 = sp.N('2.46752922945601')
```

to:

```
# zeta_max^(F1) is an external carry-forward; not derived in this stage.
# Source: notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md (Family-1 specialization).
zmax_F1 = sp.N('2.46752922945601')
```

Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl`:

Replace line 61 from:

```
zMaxF1 = SetPrecision[2.46752922945601, 20];
```

to:

```
(* zeta_max^(F1) is an external carry-forward; not derived in this stage.
   Source: notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md (Family-1 specialization). *)
zMaxF1 = SetPrecision[2.46752922945601, 20];
```

If you can locate the upstream stage that derives `zeta_max^(F1) = 2.46752922945601` from a Family-1 closed form, name it in the comment instead of pointing to the notes file. Do NOT search broadly — only use information already visible. If you cannot, leave the comment pointing to the notes file.

Do NOT change the numeric literal `2.46752922945601`. Do NOT change the `expectApprox` targets `0.456730991107963169017835980412` and `2.01079823834804688464927835412`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 098` and `redteam exec-mathematica 098` and confirm both scripts still exit 0 and produce the same numeric output. The new comment lines should appear at the indicated sites.

## F3 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:23-27`

**Issue:** The SymPy script's Family-1 check is just `assert gap_F1 > 0`. This is implied by line 21's symbolic gap factorization (`3 zmax (zmax-1)/(3 zmax-2)`) combined with `zmax_F1 > 1`, so the assertion adds nothing and cannot catch a typo in the hardcoded `2.46752922945601`. The Mathematica counterpart pins both `zetaEdgeF1` and `gapF1` to specific literals with `1e-15` tolerance. Mirror that on the SymPy side.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py`, replace the block at lines 23–27:

Before:
```
# Family-1 specialization.
zmax_F1 = sp.N('2.46752922945601')
zeta_edge_F1 = sp.N(zeta_edge.subs(zmax, zmax_F1), 30)
gap_F1 = sp.N(zmax_F1 - zeta_edge_F1, 30)
assert gap_F1 > 0
```

After (assuming F2's comment is added):
```
# Family-1 specialization.
# zeta_max^(F1) is an external carry-forward; not derived in this stage.
# Source: notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md (Family-1 specialization).
zmax_F1 = sp.N('2.46752922945601')
zeta_edge_F1 = sp.N(zeta_edge.subs(zmax, zmax_F1), 30)
gap_F1 = sp.N(zmax_F1 - zeta_edge_F1, 30)
assert gap_F1 > 0
# Numeric pin (matches Mathematica's expectApprox targets):
zeta_edge_F1_target = sp.Float('0.456730991107963169017835980412', 30)
gap_F1_target = sp.Float('2.01079823834804688464927835412', 30)
assert abs(zeta_edge_F1 - zeta_edge_F1_target) < sp.Float('1e-15', 30)
assert abs(gap_F1 - gap_F1_target) < sp.Float('1e-15', 30)
```

If F2 has not been applied (e.g., user blocked it), still add the four new lines after `assert gap_F1 > 0`; the comment lines are independent.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 098` and confirm the script exits 0 and the printed output still shows `Family-1 zeta_edge = 0.456730991107963169017835980412` and `Family-1 margin = 2.01079823834804688464927835412` (sympy's `sp.N(..., 30)` produces those exact 30-digit strings — same as the existing output transcript).

---
## Applied: F1 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster C direction (a))

- files_changed: scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py, mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl
- summary: Script-side banners and docstrings updated to STAGE 098 to match the audit-unit number (option (a) in the resolve block). Paper card section title "Stage 115" stays for now and is flagged to PAPER_CLEANUP_TRACKER for a future paper-side pass per Cluster C direction (a). Note section 4 reference to "Stages 80-81" left as historical context.
- deviation: paper-side edits explicitly deferred; only script-side relabeling applied.

## Applied: F2
- files_changed: scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py, mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl
- summary: Added provenance comment naming notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md as the source of zeta_max^(F1) = 2.46752922945601.
- deviation: none

## Applied: F3
- files_changed: scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py
- summary: Mirrored the Mathematica numeric pins on the SymPy side (zeta_edge_F1_target = 0.456..., gap_F1_target = 2.01..., abs() < 1e-15).
- deviation: none
