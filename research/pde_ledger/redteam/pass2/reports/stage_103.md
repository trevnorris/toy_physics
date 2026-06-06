---
unit_id: 103
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage103_reduced_25pn_conditional_closure.md]
  paper_appendix: present
---

# Audit unit 103 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_103.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage103_reduced_25pn_conditional_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (the stage's appendix row is the `\input{stages/stage_103}` at line 1240 — the card itself is the row; no separate summary paragraph)
- sympy: `(missing)` — by design (status-only stage)
- mathematica: `(missing)` — by design (status-only stage)
- sympy output: `(missing)`
- mathematica output: `(missing)`

(Observation, not a finding: `notes/stages/review/stage_103_review.md` exists but is a stale/mis-numbered first-pass review artifact whose body actually reviews "Stage 035 dimensionless normalization locus," not the Stage 103 conditional-closure content. It is outside this audit's script scope and is not an authority for Stage 103's claim; it has no bearing on the card↔notes consistency below. Flagged here for orchestrator awareness only.)

## What the paper claims

Stage 103 is a **status/consolidation ledger step** in the retarded 2.5PN factorization block. Per the card, `\stagefield{Verification}` reads verbatim: "SymPy audit: none yet. Mathematica audit: none yet." — i.e. the stage explicitly carries **no executable script of its own** and consolidates results carried forward from earlier units. Its bottom-line claim (`\stagefield{Derivation ledger}` + the boxed quote) is the **conditional closure statement**: the reduced product `\(\widehat m_0^{\,2}\chi_Q N_Q\)` separates the source, conservative, and outgoing factors, and "Reduced theorem closes iff the actual passive/outgoing branch has \(\chi_Q=1\)." The notes elaborate the conditional structure: inside the reduced hierarchy the conservative isotropic branch is geometry-clean through `O(omega^4)`, the conservative quadrupole split is exactly `3/4 + 1/4`, the selected support/source ratio is `rho_alpha = 4/3`, the source-map gives `mhat_0 = 1 + O(a^2/r^2)`, higher odd retarded data begin at `O(omega^7)` (irrelevant to point-particle 2.5PN), and therefore the only remaining branch datum is `chi_Q`. The failure measure if `chi_Q != 1` is defined as `Delta_Q := chi_Q - 1`. The claim status is explicitly conditional (`\StatusExactClosure{} / \StatusOpen{}`), and the downstream note instructs that the conditional status tag be carried when cited. No deliverable is an unconditional numeric/symbolic theorem; the deliverable IS the closure condition `chi_Q = 1` and the failure-measure definition `Delta_Q := chi_Q - 1`.

## What the script claims to verify

There is no script. Both engines are absent **by design**, as the card's `Verification` line states directly. Per the prompt's status-only handling (`is_status_only_candidate: True`), the absence of executable scripts is legitimate for a consolidation/carry-forward stage and is **not** a `missing_verification_script` finding unless the carry-forward references a result no upstream unit's scripts actually verify. The constituent results the card/notes carry forward are produced and tested in earlier units of this block: the `O(omega^4)` geometry-clean conservative branch, the `3/4 + 1/4` conservative quadrupole split, `rho_alpha = 4/3`, the Family-1 sufficiency, `mhat_0 = 1 + O(a^2/r^2)` (source map), and the `O(omega^7)` higher-odd irrelevance — the last of which is exactly the immediately-preceding unit 102 (audited clean, `higher_odd_irrelevance`). Stage 103 introduces no new computed quantity of its own; it states the conditional closure logic assembled from these. There is therefore nothing on the script side to attack.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Conditional closure: `mhat_0^2 chi_Q N_Q` separates factors; closes iff `chi_Q = 1` | none (status-only, by design) | n/a — no script required |
| Failure measure `Delta_Q := chi_Q - 1` | none (definitional) | n/a — definitional, no computation |
| Carry-forward: `O(omega^4)` geometry-clean / `3/4 + 1/4` split / `rho_alpha = 4/3` / `mhat_0 = 1 + O(a^2/r^2)` / `O(omega^7)` higher-odd irrelevance | verified in upstream units (e.g. 102 for higher-odd irrelevance) | carry-forward, not re-verified here |

No row is `mismatch` or `missing` in the finding sense. The card and notes agree: the card's `mhat_0^2 chi_Q N_Q = 1` is consistent with the notes' `mhat_0 = 1 + O(a^2/r^2)` (with `mhat_0 -> 1` the product reduces to `chi_Q N_Q`, and the conservative normalization `N_Q` carries the `3/4 + 1/4` split), so the closure condition collapses to `chi_Q = 1` exactly as both documents state. `paper_alignment: aligned`.

## Assertion inventory

No scripts exist, so there are no assertions to inventory.

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | (none) | — | (no executable script — status-only by design) | n/a | n/a |

## Findings

None. The absence of both engines is by design for this status-only consolidation stage and is explicitly declared in the card (`SymPy audit: none yet. Mathematica audit: none yet.`). No script-side defect can exist where no script exists, and the carry-forward chain it consolidates is verified in upstream units. The card and notes are mutually consistent on every stated deliverable. Per the audit instructions, a missing engine is NOT manufactured into a finding here.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists.

## Engine cross-check

Not applicable — neither engine is present; nothing to cross-check.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix context, and confirmed Stage 103 is a status-only consolidation step that, by design and by its own explicit `Verification` line, has no SymPy or Mathematica script. The conditional-closure claim (`chi_Q = 1` closes the reduced 2.5PN theorem; `Delta_Q := chi_Q - 1` measures any residual failure) is internally consistent between the card's `Derivation ledger`/boxed quote and the notes' statement, and the carried-forward ingredients trace to upstream units (notably the `O(omega^7)` higher-odd irrelevance to unit 102). Attacks attempted and failed to land: (1) treating the two absent engines as `missing_sympy`/`missing_mathematica` — rejected, because the status-only carve-out applies and the carry-forward results are verified upstream, not orphaned; (2) hunting a card↔notes contradiction in the `mhat_0^2 chi_Q N_Q = 1` vs `mhat_0 = 1 + O(...)` forms — they reconcile, since `mhat_0 -> 1` collapses the product to `chi_Q = 1`; (3) checking for an unconditional theorem overclaim — the card correctly tags the status as conditional (`\StatusExactClosure / \StatusOpen`) and instructs downstream callers to carry the conditional tag. There are no script-emitted values to reconcile (see below).

## Value Reconciliation (pass-2 augmentation)

This is a **status-only stage with no scripts at all** (no `.py`, no `.wl`, no committed `.txt` outputs). There is therefore **no script output to reconcile against** — the standard procedure (enumerate every value the scripts emit, then locate it in the docs) yields an **empty deliverable set on the script side**. The card and notes state symbolic conditions/definitions only (`chi_Q = 1`, `Delta_Q := chi_Q - 1`, `mhat_0^2 chi_Q N_Q = 1`, `mhat_0 = 1 + O(a^2/r^2)`, `rho_alpha = 4/3`, the `3/4 + 1/4` split, `O(omega^4)`/`O(omega^7)` orders); none of these is produced by an engine in this unit — they are carried forward from upstream units and consolidated here. Per the augmentation's status-only guard ("Status-only stages (no scripts ...): reconcile whatever the present engine emits; note the single-engine status"), there is no present engine and nothing to reconcile.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| (none — no script emits any value) | (no scripts exist) | — | n/a |

Internal-scaffolding name-only list: (none — no scripts).

reconciliation: complete; 0 script-emitted values checked, 0 misaligned (status-only stage, no engines by design — nothing for the engines to emit).

## Self-test notes

Traps checked: (1) status-only carve-out — confirmed `is_status_only_candidate: True` applies and the card explicitly declares both audits "none yet," so absent engines are by design, not findings; (2) carry-forward orphan check — the consolidated results trace to upstream units (higher-odd irrelevance to unit 102), so no `missing_verification_script` is warranted; (3) card↔notes consistency — the `mhat_0^2 chi_Q N_Q = 1` product reconciles with `mhat_0 = 1 + O(a^2/r^2)` and collapses to `chi_Q = 1`, matching both documents, with the conditional status tag correctly carried. No `sp.diff`/parity/trivial-case traps apply since there are no assertions. Conclusion: zero findings, verdict clean; no directive written.
