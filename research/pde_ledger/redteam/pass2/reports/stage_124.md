---
unit_id: 124
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage124_core_branch_status.md]
  paper_appendix: present
---

# Audit unit 124 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_124.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage124_core_branch_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows at lines 28 and 86)
- sympy: (missing — by design, status-only)
- mathematica: (missing — by design, status-only)
- sympy output: (missing)
- mathematica output: (missing)

## Status-only determination (genuinely notes-only vs scripts-missing)

Stage 124 is **genuinely notes-only**. Verified:

- Manifest (`redteam/pass2/MANIFEST.yaml:4252-4274`): `is_status_only_candidate: true`, `is_checkpoint: false`, both `sympy.path` and `mathematica.path` are `null` with `exists: false`; both output paths null.
- No files matching `*124*` exist under `scripts/`, `mathematica/`, or their `output/` subdirs.
- The card's `\stagefield{Verification}` reads: "SymPy audit: none yet. Mathematica audit: none yet." — the card itself declares status-only.

Critically, the carve-out is **legitimate** because every concrete value the stage-124 card/notes carry forward is independently verified by upstream scripts (stages 121–123 + 125–127), so 124 does not assert any new result that lacks engine support. The carry-forward chain:

| Stage-124 deliverable (notes) | Upstream script that verifies it |
|---|---|
| `r_F1 ≈ 1.77799353547498` (geometric hybridization at L/a=37/20) | stage 121 sympy out:9 `r_F1 numeric = 1.7779935354749781185`; stage 121 mathematica out:11 same to 20 digits. Closed form `r_geom = sqrt(12L²−π²a²)/(πa)` (stage121 out:5) is algebraically identical to the note's `sqrt((12/π²)(L/a)²−1)`. |
| `g_-^F1 ≈ 0.758035078944663`, `g_+^F1 ≈ 2.79795199200529` | stage 122 mathematica out:10-11 (both, 20 digits); stage 125 sympy out:19-20; stage 126/127 corroborate. |
| `g_nat = 1` (natural equal-normalized branch) | stage 122 mathematica out:27-30: traction residual vanishes exactly at `gNat -> 1`. |
| `31.9%` traction enhancement on lower `(-)` branch | stage 122 sympy out:14 / mathematica out:31: `numeric T ratio (-) = 1.3192001633911203307` ⇒ +31.9% over natural. |

Because no stage-124 deliverable references a result that lacks an upstream engine check, **no `missing_sympy`/`missing_mathematica` finding is warranted** (per the status-only rule in the prompt). The absent scripts are by-design, not a defect.

## What the paper claims

Per `\stagefield{Purpose}` and `\stagefield{Derivation ledger}`, stage 124 is a core-outlet-realization ledger/status step: after the geometric branch selection of stages 121–123, it records that the explicit outlet-core branch is no longer a free `(r,g)` pair. The card's bottom-line claim (the boxed quote) is a status statement, not an identity: "Remaining question is whether the real core chooses the naive or lower compensated mouth branch." The notes sharpen this: they pin the four concrete values above and frame the residual question as whether the true GNLS + localized-Maxwell core picks the natural equal-normalized branch (`g_nat=1`) or shifts onto the lower compensated branch `g_-^F1` (a 31.9% traction enhancement). The appendix rows (lines 28, 86) place 124 as the closing stage of the "two-channel core realization, parent extraction, geometric core selection" / "classify admissible outlet/core deformations" block. There is no `\stagefield{Output}` field with a numeric closed form — the Output is the status verdict. `\stagefield{Checks}` lists three review reminders (Schur signs, D/N length normalization vs chosen L/a, parent overlap ratios) — none introduces a new asserted value.

## What the script claims to verify

No scripts exist (status-only by design). There is nothing to transliterate, no assertions to inventory, no engine cross-check. The verdict therefore applies to (a) internal consistency of the card+notes and (b) the value-reconciliation of every value the notes assert against the upstream scripts that actually produce them.

## Paper ↔ script cross-check

| Paper-side deliverable | Coverage |
|---|---|
| Status verdict: explicit branch is fixed to `r_F1`, choice is now naive `g_nat` vs lower `g_-^F1` | `match` (carry-forward) — all numeric inputs to this verdict verified upstream (121/122/125/126/127); the verdict itself is a status statement requiring no new check |
| `r_F1 ≈ 1.778` | `match` — stage 121 (both engines) |
| `g_-^F1, g_+^F1` | `match` — stage 122/125 (both engines) |
| `g_nat = 1` | `match` — stage 122 (both engines) |
| `31.9%` enhancement | `match` — stage 122 T-ratio(−) = 1.3192 (both engines) |

`paper_alignment: partial` is set not because of any value mismatch (all values reconcile) but because the card carries an internal stale self-label (F1 below): its Purpose field names the wrong stage number.

## Assertion inventory

No scripts → no assertions. Not applicable.

## Findings

### F1 — paper_misalignment (subtype: notes_contradicts_script / stale self-label)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_124.tex:7`

**What's wrong:**
The card is unambiguously the stage-124 card — title, `\label{stage:124}`, and anchor are all 124 — but the `\stagefield{Purpose}` body opens with the wrong stage number:

> `\stagefield{Purpose}{Stage~141 is a core outlet realization ledger step.  Its audit target is the verification output quoted below.}`

It should read "Stage~124", not "Stage~141". This is the known `+17` numbering-drift artifact (124 + 17 = 141) from the EM-extension realignment. It is content-keyed, not ambiguous: the descriptor "core outlet realization ledger step" matches stage 124's subject (the appendix audit-path map at `stage_appendix_part04.tex:28` assigns "core realization, parent extraction, and geometric core selection" to stages 114–124), whereas the genuinely-distinct stage 141 (`paper/stages/stage_141.tex:1,7`) is titled "Mouth-Gain Status Update" and its own Purpose correctly reads "Stage~141 is a coupled mouth fixed point and gain selection ledger step." So 124's "Stage~141" is a self-label that drifted, not a real cross-reference to 141.

**Why this matters:**
A reader following the card believes stage 124 is some other stage's step, and the stale number could mask or be masked by the larger script/output-band numbering work. It is cosmetic-to-load-bearing: it does not change any math, but it is a genuine inconsistency in a published-paper card.

**Required change (routes to USER, not Codex):**
This is a paper-side (`.tex`) prose edit. Per the red-team contract, Codex does not edit `paper/`. Surface to the user: change `paper/stages/stage_124.tex:7` "Stage~141" → "Stage~124". The fix is content-determined (do not offset-sweep), so the direction is clear, but the edit is in a file the red-team scripting loop does not own.

**Verification:**
After the user authorizes, `grep -n "141" paper/stages/stage_124.tex` returns no match in the Purpose field; the only stage numbers in the card are 124 (title/label) and the genuine downstream range 125–139.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` for this unit (status-only by design).

## Engine cross-check

Not applicable — no scripts for this unit. The values it carries forward were cross-checked between engines at their source stages (121/122/125): e.g. `r_F1` agrees to 20 digits between stage 121 sympy and mathematica; `g_-` agrees to 20 digits between stage 122 sympy and mathematica; `T ratio (−) = 1.3192...` agrees between stage 122 sympy out:14 and mathematica out:31.

## Verdict justification

Stage 124 is genuinely notes-only and the status-only carve-out is legitimate: the manifest pins both engines null, no script files exist, and — the binding test — every value the card/notes carry forward is independently verified by upstream engine scripts (121/122/125/126/127), so 124 asserts no unsupported new result. I attacked the carry-forward by checking each of the four asserted values against the committed upstream outputs and confirmed all four match to full precision in both engines, including re-deriving `r_F1` from `L/a=37/20` (≈1.778) and confirming the 31.9% figure is the stage-122 `(-)`-branch T-ratio 1.3192. The single finding is a low-severity stale self-label in the card's Purpose field ("Stage~141" should be "Stage~124"), a `+17` numbering-drift artifact that routes to the user since it is a paper-prose edit Codex does not own. Verdict is `findings` (count 1, paper-side, user-resolution) rather than `clean` solely because of that self-label; no value mismatch, no missing-deliverable, no internal numeric inconsistency exists.

## Value Reconciliation (pass-2 augmentation)

No scripts exist for stage 124, so there are no script-emitted values of its own to reconcile. Per the augmentation's status-only guidance, I instead reconcile every value the stage-124 **card/notes assert** against the upstream scripts that produce them (the values are deliverables 124 carries forward) and against the prose carriers (the notes `.md`; the card `.tex` is terse by design and legitimately omits the numerics).

| value | source (upstream py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_F1 ≈ 1.77799353547498` | stage121 sympy out:9 `1.7779935354749781185`; stage121 wl out:11 (20 digits) | notes `...stage124_core_branch_status.md:16`; card omits (terse) | MATCH |
| closed form `r_F1 = sqrt((12/π²)(L/a)²−1)` | stage121 sympy out:5 `r_geom = sqrt(12*L**2 - pi**2*a**2)/(pi*a)` (algebraically identical) | notes md:11-16 | MATCH |
| `L/a = 37/20` | stage121 (selection point) | notes md:16 `L/a=37/20` | MATCH |
| `g_nat = 1` | stage122 wl out:27-30 (residual vanishes at gNat→1) | notes md:21 | MATCH |
| `g_-^F1 ≈ 0.758035078944663` | stage122 wl out:10; stage125 sympy out:19 | notes md:26 | MATCH |
| `g_+^F1 ≈ 2.79795199200529` | stage122 wl out:11; stage125 sympy out:20 | notes md:27 | MATCH |
| `31.9%` traction enhancement (lower `(-)` branch) | stage122 sympy out:14 / wl out:31 `T ratio (-) = 1.3192...` | notes md:31 `31.9%` | MATCH |

Internal scaffolding (no finding): none — stage 124 emits nothing of its own; all listed values are deliverables, not scaffolding.

The one card-side discrepancy (the "Stage~141" self-label, F1) is a stale **label**, not a **value**, so it is folded into the standard `## Findings` as a `paper_misalignment` self-label rather than a `value_mismatch`.

reconciliation: complete; 7 carry-forward values checked, 0 value-misaligned (1 stale stage-number label flagged separately as F1)

## Self-test notes

- **Variable independence / parity / trivial-case traps:** not applicable — no scripts, no new assertions or derivatives to construct, so the derivative-of-independent-symbol and parity traps cannot arise here.
- **Carry-forward support:** I confirmed each of the four asserted numbers traces to a concrete upstream output line in BOTH engines (not just one), so the status-only carve-out is not hiding an unverified value.
- **Numbering-drift trap:** I confirmed the "Stage~141" is content-keyed stale (descriptor matches 124's block, stage 141 is a genuinely different mouth-gain card) and did NOT offset-sweep — I read 141's actual card to disprove a real cross-reference. No directive is written (the only finding is a paper-side self-label routed to the user; Codex owns no script here).
