---
unit_id: 120
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T16:15:38Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage120_core_parameter_status.md]
  paper_appendix: present
---

# Audit unit 120 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_120.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage120_core_parameter_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (relevant rows/blocks read: line 19 anchor description; subsec "Parent-action overlap extraction" lines 575-601; core-coupling-balance + r/g/compensation-family block lines 520-554; result-anchor citation index lines 1175-1179; `\input{stages/stage_120}` at line 1274)
- sympy: (missing — status-only by design)
- mathematica: (missing — status-only by design)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 120 is a **status-only ledger entry** (manifest `is_status_only_candidate: true`; checkpoint: false). Its `\stagefield{Verification}` field states verbatim "SymPy audit: none yet. Mathematica audit: none yet." — so by design no executable script is expected. The card's bottom-line claim, in the quote block, is: "Records the surviving finite parent-level target after overlap extraction." The notes (the authoritative carrier here) state the deliverable concretely: after explicit GNLS + localized-Maxwell reduction the throat-core outlet is controlled by the microscopic set `{K_s, K_q, λ, g_q, g_s, L_W}` reduced to the two normalized ratios `r = λ/√(K_s K_q)` and `g = g_q√K_s/(g_s√K_q)`, and the surviving theorem gate is `1 + r² = 4(g − r)²` together with `L_W = (πa/2)√((1+r²)/3)`. The deliverable is explicitly a *status reframing*: the open question becomes "what branch values of (r, g) does the actual GNLS + localized-Maxwell throat core select?" rather than a closed numeric result. The card carries `\StatusExactClosure / \StatusOpen` and is anchored to result `MTDC-T8`.

## What the script claims to verify

No script exists for this unit (none in `scripts/`, `mathematica/`, or either `output/` directory — confirmed by directory listing). Per the prompt and manifest this is legitimate status-only: the unit consolidates/reframes results derived in the surrounding Part IV appendix narrative rather than running its own algebra. There is therefore nothing to attack at the assertion level; the audit reduces to (a) confirming the absent-engine status is by-design (it is), and (b) reconciling every value the card and notes *assert* against the appendix block they cite and against each other.

## Paper ↔ script cross-check

No script side exists, so the cross-check is card/notes ↔ appendix narrative (the carry-forward source the card points to via "the corresponding block narrative above").

| Deliverable (card/notes) | Source in appendix | Status |
|---|---|---|
| `r := λ/√(K_s K_q)` | `stage_appendix_part04.tex:533` | match |
| `g := g_q√K_s/(g_s√K_q)` | `stage_appendix_part04.tex:535` | match |
| Gate `1 + r² = 4(g − r)²` | `stage_appendix_part04.tex:540` (`eq:app-part04-parent-compensation-family`) | match |
| `L_W = (πa/2)√((1+r²)/3)` | `stage_appendix_part04.tex:548` (`eq:app-part04-LW-compensation`) | match |
| Microscopic controls `{K_s, K_q, λ, g_q, g_s, L_W}` | parent definitions `stage_appendix_part04.tex:581-599` | match |
| "surviving finite parent-level target after overlap extraction" (status reframing) | subsec parent-overlap-extraction lines 575-601, esp. line 601 ("future full branch calculation must decide whether [the compensation family] is actually realized") | match |

`paper_alignment: aligned` — every asserted value/identity in the card and notes is faithfully present in, and agrees with, the Part IV appendix block the card cites. No value mismatch, no missing deliverable.

## Assertion inventory

No scripts → no assertions. (Table omitted; status-only by design.)

## Findings

### F1 — stale self-label (paper-side prose; numbering drift)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_120.tex:7`

**What's wrong:**
The card's `\section`, `\label`, and filename all correctly identify the unit as **120** (`\section[Stage~120]{Stage~120: Core-Parameter Extraction Status}`, `\label{stage:120}`). But the `\stagefield{Purpose}` field self-labels with the wrong number:

> `\stagefield{Purpose}{Stage~137 is a core outlet realization ledger step.  Its audit target is the verification output quoted below.}`

`137 = 120 + 17`, the known EM-extension renumber-drift signature. I confirmed this is a systematic stale-label band, not a real cross-reference: the Purpose self-labels in this neighborhood are correct through 116, then drift by exactly +17 from 117 onward (117→134, 118→135, 119→136, 120→137, 121→138, 122→139, 123→140, 124→141). Stage 137's *own* card has a genuinely different Purpose ("Stage~137 is a coupled mouth fixed point and gain selection ledger step"), so the "137" in stage 120's Purpose is not a legitimate forward reference to stage 137 — it is a stale self-reference that should read "Stage~120". The template phrase "core outlet realization ledger step" is shared across the 114-124 cards, so the *content* of the Purpose line is correct for 120; only the embedded number is stale.

**Why this matters:**
A reader (or downstream citation) taking the Purpose field literally would mis-attribute this status entry to stage 137. It is cosmetic/label-only (no value or identity is affected, and no script depends on it), but it is a genuine stale self-label in a Part IV card.

**Required change:**
This is a **paper-side prose edit**, NOT a script edit. Codex must not auto-apply it under the red-team script contract. Change `stage_120.tex:7` `Stage~137` → `Stage~120`. Route to the paper/notes owner (per the file-ownership policy, Codex applies published-paper edits only after user/Claude review of the paper-side change; the orchestrator should fold this into the dedicated numbering script/output-band + label cleanup track rather than the script fix_loop). No directive for Codex script work is written for this unit.

**Verification:**
After the fix, `grep "Stage~" stage_120.tex` should show only `Stage~120` (no `Stage~137`); the neighbor band 117-124 should likewise be de-drifted if handled as a batch.

## Independent-derivation check (Mathematica)

N/A — no `.wl` exists. Status-only by design; not a finding.

## Engine cross-check

N/A — neither engine present.

## Verdict justification

`verdict: findings`. The unit is **genuinely notes-only**, not scripts-missing: the card explicitly declares "SymPy audit: none yet / Mathematica audit: none yet," the manifest marks it status-only (checkpoint false), and it carries forward results (the `r`/`g` ratios, the compensation-family gate, and the `L_W` law) that ARE derived and stated in the Part IV appendix block it cites — so the absent engines are correct, not a coverage gap. I attacked the value content: every identity and definition asserted by the card and notes matches the appendix (`r`, `g`, the gate `1+r²=4(g−r)²`, and `L_W=(πa/2)√((1+r²)/3)` all reconcile verbatim, lines 533/535/540/548), and the notes are internally consistent with the card's status framing. The one defect found is a stale self-label in the Purpose field ("Stage~137" for unit 120, a +17 renumber-drift artifact), which is paper-side prose, label-only, low severity, and routed to the paper owner — not a script finding and not stop-cold. No `paper_misalignment` value finding, no missing-script finding, no UNFIXABLE/CRITICAL_DOWNSTREAM condition.

## Self-test notes

No scripts means no derivative/parity/trivial-case traps to walk (those apply to proposed assertions; none are proposed here). I checked: (1) absent engines are by-design status-only, not a missing-script finding — confirmed via card "none yet" + manifest flag + directory listing; (2) value round-trip — all four asserted identities/definitions reconcile to the cited appendix equations with no constant or sign discrepancy; (3) the "137" is a stale label, not a real cross-ref — confirmed by the systematic +17 drift across the 117-124 Purpose band and by stage 137's distinct own-Purpose. The only finding is paper-side prose, so no Codex script directive is written.

## Value Reconciliation (pass-2 augmentation)

No scripts emit values for this unit (status-only, both engines null). Per the augmentation's status-only guard, I reconcile the values the card/notes *assert* for internal/cross-stage consistency against the Part IV appendix block the card cites. There are no script-emitted numeric constants; the deliverables are symbolic definitions/identities plus the status reframing.

| value | source (card/notes — no scripts) | .tex/.md location | status |
|---|---|---|---|
| `r := λ/√(K_s K_q)` | notes lines 13-14; card "parent overlap ratios" | appendix `stage_appendix_part04.tex:533` | MATCH |
| `g := g_q√K_s/(g_s√K_q)` | notes lines 17-19 | appendix `stage_appendix_part04.tex:535` | MATCH |
| gate `1 + r² = 4(g − r)²` | notes lines 25-27 | appendix `stage_appendix_part04.tex:540` (`eq:app-part04-parent-compensation-family`) | MATCH |
| `L_W = (πa/2)√((1+r²)/3)` | notes lines 29-31 | appendix `stage_appendix_part04.tex:548` (`eq:app-part04-LW-compensation`) | MATCH |
| microscopic control set `{K_s, K_q, λ, g_q, g_s, L_W}` | notes lines 9-20 | appendix parent-action defs `stage_appendix_part04.tex:581-599` | MATCH |
| status reframing ("surviving finite parent-level target after overlap extraction") | card quote block (lines 15-17); notes lines 33-41 | appendix subsec lines 575-601 (esp. 601) | MATCH |

INTERNAL items (scaffolding not expected in prose): none — there are no scripts, hence no verification scaffolding, residuals, flags, or tolerances to account for.

The single defect (F1) is a stale *self-label number* in the Purpose field, not a deliverable-value mismatch, so it does not appear as a `value_mismatch` row above; it is recorded as the standalone F1 finding instead.

reconciliation: complete; 6 asserted values checked, 0 misaligned
