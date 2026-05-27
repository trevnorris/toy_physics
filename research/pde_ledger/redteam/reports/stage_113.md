---
unit_id: 113
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage113_outlet_model_status.md
  paper_appendix: present
---

# Audit unit 113 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_113.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage113_outlet_model_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows touching 113: line 27 audit-path map, line 86 chain summary, lines 399-456 outlet narrative containing the explicit Robin / mixed / hybrid claims that stage 113 consolidates, line 1260 `\input{stages/stage_113}`)
- sympy: (missing — manifest marks `is_status_only_candidate: true`, `files.sympy.path: null`)
- mathematica: (missing — manifest marks `files.mathematica.path: null`)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 113 is explicitly framed in its card as a **status update** — the `\stagefield{Verification}` line literally reads "SymPy audit: none yet. Mathematica audit: none yet." The card's `Purpose` says the stage is an "outlet deformation and compensation ledger step" and its `Inputs` import "the canonical outgoing DtN expansion, a general isotropic deformation, Robin mouth loading, and a hidden mixed side-channel pole." The card body is terse; the operative content is the boxed quote: "Low-frequency outlet audit leaves pure scale/argument classes and the compensated Robin--mixed class as the viable routes." The companion notes (`moving_throat_pde_stage113_outlet_model_status.md`) make the consolidation explicit by listing the three outlet classes already checked upstream: pure Robin core (`\chi_Q^R = 3/(3-\rho_R)`), standalone mixed pole (must vanish to preserve the even branch), and the compensated hybrid with `(\rho_R,\kappa_W)=(4\sigma_W, 1/3)`, `\chi_Q^{hyb} = (1-9\sigma_W\gamma_W)/(1-\sigma_W)`, canonical iff `\gamma_W = 1/9`. The appendix re-states the same three identities at lines 402-456 with equation labels. The `\stagefield{Checks}` items are statements about what was checked in the consolidation chain, not new claims this stage produces independently. `\stagefield{Downstream use}` flags that stages 114-141 cite this status with its tag attached.

## What the script claims to verify

No script exists for this unit. There is no `.py` and no `.wl` file under the conventional names and the manifest confirms both file paths are `null` with `exists: false`. There is therefore no script-side assertion surface to attack.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `\chi_Q^R = 3/(3-\rho_R)` (pure Robin) | none in this unit; relies on upstream stage in the 107-112 outlet block | n/a (status carry-forward) |
| Standalone mixed pole forces `\sigma_W=0` (no-go) | none in this unit; relies on upstream | n/a (status carry-forward) |
| Compensated branch `(\rho_R,\kappa_W)=(4\sigma_W,1/3)`, `\chi_Q^{hyb}=(1-9\sigma_W\gamma_W)/(1-\sigma_W)`, canonical iff `\gamma_W=1/9` | none in this unit; relies on upstream | n/a (status carry-forward) |
| Audit-path summary "pure scale/argument or compensated Robin–mixed are the only viable routes" | none — this is the consolidation claim | n/a (status carry-forward) |

The dominant pattern is "no script-side check because this unit is a status ledger entry that summarizes what the 107-112 chain has already established." The paper card itself transparently states `Verification: SymPy audit: none yet. Mathematica audit: none yet.`, so the paper and the (absent) script are mutually consistent: both sides agree there is no independent verification at this unit. `paper_alignment` is therefore `aligned`.

## Assertion inventory

(No assertions. No script.)

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | — | — | — | — | — |

## Findings

(None.)

The manifest entry `is_status_only_candidate: true` legitimizes the absence of both engines per the prompt's status-only carve-out. The carve-out reads: "A `missing_sympy` or `missing_mathematica` finding is only valid if the unit's scripts (or comments referencing source units) reference a result that no upstream unit's scripts actually verify."

I cannot file a `missing_sympy`/`missing_mathematica` finding under this rule because:

1. There are no scripts here at all — no comments to inspect for upstream references.
2. The reading list for this audit explicitly excludes other units' scripts and the `notes/` trackers, so I am not permitted to confirm or refute whether stages 107-112 actually verified the three carry-forward claims. The carve-out's "if no upstream unit's scripts actually verify" condition is unverifiable from inside this unit's reading scope, and the prompt forbids me from invoking it speculatively.
3. The paper card is transparent — it states "SymPy audit: none yet. Mathematica audit: none yet." rather than claiming a verification that does not exist. There is no paper↔script mismatch because both sides disclose the same null state.

If a sweep across all of Part IV reveals that no upstream unit verifies the compensated branch identities (`(\rho_R,\kappa_W)=(4\sigma_W,1/3)`, `\chi_Q^{hyb}`, `\gamma_W=1/9`), that is a cross-unit coverage problem properly raised at the batch-tracker level, not a stage-113 finding.

## Independent-derivation check (Mathematica)

Not applicable — no Mathematica script exists.

## Engine cross-check

Not applicable — neither engine present.

## Verdict justification

`clean`. The unit is status-only; the manifest, the paper card's `Verification` line, and the absence of script files are mutually consistent. I attempted three lines of attack: (a) check whether the card or notes assert a claim that should have an in-unit assertion — they do not; the card explicitly disclaims having a SymPy/Mathematica audit and the notes only list classes "already checked" elsewhere; (b) check whether the body equations in the appendix's `Robin and mixed-pole tests` / `Compensated Robin--mixed outlet` subsections (lines 399-456) introduce a numeric constant the card secretly relies on without provenance — every constant (`3/(3-\rho_R)`, `4\sigma_W`, `1/3`, `1/9`, `(1-9\sigma_W\gamma_W)/(1-\sigma_W)`) is derived in the appendix itself, not hard-coded in this stage; (c) check whether the card's `Downstream use` language overpromises (e.g., asserting that the compensated branch IS the actual branch) — it does not; the language reads "The card is a derivation ledger entry, not an unconditional actual-branch theorem," which is exactly the status-only framing the carve-out permits. All three attacks failed. I read paper card, notes, and the appendix rows that touch this stage before reaching this verdict.

## Self-test notes

- Variable independence / symmetry / trivial-case traps: not applicable — no derivatives, integrals, or assertions to mentally simulate, because no script exists.
- Path-specification trap: not applicable — no `missing_verification_script` finding is being filed, so no target path needs to be named.
- Paper round-trip: re-checked the appendix lines 402-456 against the notes; the equations match (Robin `\chi_Q^R = 3/(3-\rho_R)`; hybrid solutions `(\sigma_W,0)` and `(4\sigma_W,1/3)`; `\chi_Q^{hyb} = (1-9\sigma_W\gamma_W)/(1-\sigma_W)`; canonical at `\gamma_W=1/9`). No internal contradiction between card, notes, and appendix.
- Carve-out gate: confirmed `is_status_only_candidate: true` in MANIFEST.yaml line 3894 and `is_checkpoint: false` in line 3893, so the higher checkpoint bar does not apply.
