---
unit_id: 113
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage113_outlet_model_status.md]
  paper_appendix: present
---

# Audit unit 113 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_113.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage113_outlet_model_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows/sections referencing stage 113: lines 17, 27, 78, 86, 345–456, 1177, 1189, 1260)
- sympy: (missing — by design, status-only stage)
- mathematica: (missing — by design, status-only stage)
- sympy output: (missing — by design)
- mathematica output: (missing — by design)

## Status-only / no-engine confirmation

Stage 113 is a **status-only consolidation stage with NO scripts by design**. There is no
`scripts/*113*.py`, no `mathematica/*113*.wl`, and no committed output `.txt`. This is
explicitly stated in the card:

> `\stagefield{Verification}{SymPy audit: none yet.  Mathematica audit: none yet.}`

and the manifest marks `is_status_only_candidate: True`, `is_checkpoint: False`. The card
is self-described as "a outlet deformation and compensation ledger step" whose "audit
target is the verification output quoted below," and `\stagefield{Downstream use}` reiterates
"the card is a derivation ledger entry, not an unconditional actual-branch theorem." The
absent engines are therefore **not** a finding (per the status-only carve-out and the audit
prompt's explicit instruction for this unit). The substantive math (deformation algebra,
Robin/mixed/hybrid outlet results) is fully carried and derived in the Part IV appendix
section `Outlet DtN and Robin outlet tests` (lines 345–456), which the card points to via
`\stagefield{Verification note}` ("collected in the corresponding block narrative above and
in the source stage").

## What the paper claims

Stage 113 consolidates the **explicit outlet-deformation audit** for the moving-throat
canonical outgoing branch. The card's bottom-line claim (the quoted result) is:

> "Low-frequency outlet audit leaves pure scale/argument classes and the compensated
> Robin--mixed class as the viable routes."

The notes and the Part IV appendix make this concrete via three explicit outlet classes:
(1) a **pure geometric Robin core** `Λ₂^R = Λ₂^out + ρ_R`, shifting normalization to
`χ_Q^R = 3/(3−ρ_R)`; (2) a **standalone mixed `A_w/F_{μw}` hidden pole**
`Λ₂^mix = Λ₂^out − σ_W/(1−κ_W z²−iγ_W z⁵)`, which cannot preserve the fixed even `l=2`
branch unless it vanishes (appendix: forces `κ_W=−1/9` then `σ_W=0`); and (3) a **hybrid
Robin–mixed outlet** `Λ₂^hyb = Λ₂^out + ρ_R − σ_W/(1−κ_W z²−iγ_W z⁵)`, which admits one
nontrivial compensated even branch `(ρ_R, κ_W) = (4σ_W, 1/3)`, on which
`χ_Q^hyb = (1−9σ_W γ_W)/(1−σ_W)`, preserving canonical outgoing normalization exactly when
`γ_W = 1/9`. The PDE-facing deliverable: the remaining branch must be either a harmless
pure scale/argument deformation or this compensated Robin–mixed outlet; pure Robin loading
alone and a naive standalone hidden mixed pole are insufficient. The appendix summary rows
(lines 27, 86, 1189) and result anchor `MTDC-T8(.3)` corroborate this verbatim.

## What the script claims to verify

No script exists for this unit (status-only by design), so there is no script-side claim to
audit. The verification burden for the underlying algebra lives in the Part IV appendix
narrative and in the upstream/adjacent outlet stages (107–124 per appendix line 86), not in
a stage-113 script. Per the status-only carve-out, a `missing_sympy`/`missing_mathematica`
finding would only be valid if stage 113 referenced a result that no upstream unit verifies;
the card explicitly imports prior results ("imports the canonical outgoing DtN expansion, a
general isotropic deformation, Robin mouth loading, and a hidden mixed side-channel pole")
and the appendix locates the supporting derivations in the 107–124 band. No orphaned,
unsupported carry-forward is present.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| Robin outlet `χ_Q^R = 3/(3−ρ_R)` | none (status-only) | n/a — supported in appendix eq:app-part04-chi-robin |
| Mixed pole vanishes (`κ_W=−1/9 → σ_W=0`) | none (status-only) | n/a — supported in appendix line 423 |
| Hybrid compensated branch `(ρ_R,κ_W)=(4σ_W,1/3)` | none (status-only) | n/a — supported in appendix eq:app-part04-hybrid-even-solutions |
| `χ_Q^hyb = (1−9σ_W γ_W)/(1−σ_W)` | none (status-only) | n/a — supported in appendix eq:app-part04-hybrid-chiQ |
| `γ_W = 1/9` preserves canonical normalization | none (status-only) | n/a — supported in appendix eq:app-part04-gammaW-canonical |

No script-side checks exist; all deliverables are carried by the appendix narrative. The
card↔notes↔appendix triad is internally consistent (see Value Reconciliation). Dominant
pattern: `aligned` (no misalignment, no extra, no mismatch).

## Assertion inventory

No assertions — no scripts exist for this unit. (Status-only stage; assertion inventory
not applicable.)

## Findings

None. (Zero findings; the absence of scripts is by design and consistent with the card,
notes, and appendix. No internal inconsistency found in the notes; all stated deliverable
values reconcile against the supporting appendix.)

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` exists (status-only by design).

## Engine cross-check

Not applicable — no engines present.

## Verdict justification

`clean`. I read the paper card, the notes file, and the Part IV appendix outlet sections in
full before concluding. Stage 113 is a status-only consolidation/ledger stage with no
scripts by design; the card explicitly states "SymPy audit: none yet. Mathematica audit:
none yet," the manifest flags `is_status_only_candidate: True` / `is_checkpoint: False`, and
the supporting algebra is carried in the appendix. Attacks tried: (a) checked whether the
notes assert any deliverable value that conflicts with the appendix derivation — every value
(`χ_Q^R=3/(3−ρ_R)`, hybrid even branch `(4σ_W, 1/3)`, `χ_Q^hyb=(1−9σ_W γ_W)/(1−σ_W)`,
`γ_W=1/9`) matches the appendix verbatim; (b) checked whether the mixed-pole "must vanish"
claim in the notes contradicts the appendix's `κ_W=−1/9 → σ_W=0` — consistent (both say the
standalone pole cannot survive); (c) checked whether any carry-forward references a result
no upstream unit supports — none; the card imports prior results and the appendix locates
their derivations in the 107–124 band; (d) checked the appendix summary rows (27, 86, 1189)
for a status overclaim relative to the card's hedged "derivation ledger entry, not an
unconditional actual-branch theorem" — consistent. No internal inconsistency, no
paper_misalignment, no missing-script finding warranted. The absent engines are by design.

## Self-test notes

Variable-independence and symmetry/parity traps are not applicable (no scripts, no
derivatives or integrals to check). I focused the self-test on the status-only carve-out
(confirmed the missing engines are by-design and the carry-forward is supported upstream, so
no `missing_verification_script` finding) and on the value round-trip: every deliverable the
notes state matches the Part IV appendix verbatim, so no `paper_misalignment` is introduced
or present.

## Value Reconciliation (pass-2 augmentation)

No scripts emit values for this stage (status-only, no engine by design), so there is no
script output to reconcile against. Per the augmentation's guard for status-only stages, I
reconcile the deliverable values the **notes STATE** against the `.tex` card and the Part IV
appendix (the natural carrier, since the terse card defers to "the corresponding block
narrative above"). All values are symbolic closed-forms; none are numeric constants pinned
by a script.

| value (stated in notes) | source (notes) | .tex / appendix location | status |
|---|---|---|---|
| `Λ₂^R = Λ₂^out + ρ_R` | notes §1 | appendix eq:app-part04-robin-outlet (line 405) | MATCH |
| `χ_Q^R = 3/(3−ρ_R)` | notes §1 | appendix eq:app-part04-chi-robin (line 410) | MATCH |
| `Λ₂^mix = Λ₂^out − σ_W/(1−κ_W z²−iγ_W z⁵)` | notes §2 | appendix eq:app-part04-mixed-pole-outlet (lines 417–421) | MATCH |
| mixed pole cannot preserve even branch unless it vanishes | notes §2 | appendix line 423 ("κ_W=−1/9, and then σ_W=0") | MATCH (consistent) |
| `Λ₂^hyb = Λ₂^out + ρ_R − σ_W/(1−κ_W z²−iγ_W z⁵)` | notes §3 | appendix eq:app-part04-hybrid-outlet (lines 431–435) | MATCH |
| compensated even branch `(ρ_R, κ_W) = (4σ_W, 1/3)` | notes §3 | appendix eq:app-part04-hybrid-even-solutions (line 442) | MATCH |
| `χ_Q^hyb = (1−9σ_W γ_W)/(1−σ_W)` | notes §3 | appendix eq:app-part04-hybrid-chiQ (line 449) | MATCH |
| `γ_W = 1/9` (preserves canonical χ_Q) | notes §3 | appendix eq:app-part04-gammaW-canonical (line 454) | MATCH |

INTERNAL (accounted, no finding): none — every stated value is a deliverable and reconciles.

reconciliation: complete; 8 values checked, 0 misaligned. Status-only/no-engine stage —
reconciled the notes' stated deliverables against the `.tex` card and Part IV appendix; no
script output exists to reconcile against (noted as the no-engine status, not a finding).
