---
unit_id: 141
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage141_mouth_gain_status.md]
  paper_appendix: present
---

# Audit unit 141 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_141.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage141_mouth_gain_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_141}` line at :1316 references this unit — there is no separate summary row)
- sympy: `(missing)` — none exists; confirmed by `find scripts -iname '*141*'` (empty) and grep for `stage141`/`stage_141` in `scripts/` (no hits)
- mathematica: `(missing)` — none exists; confirmed by `find mathematica -iname '*141*'` (empty) and grep in `mathematica/` (no hits)
- sympy output: `(missing)`
- mathematica output: `(missing)`

## What the paper claims

Stage 141 is explicitly a **status / consolidation card** in the "Coupled mouth fixed point and gain selection" block. The card's `\stagefield{Verification}` line states verbatim: "SymPy audit: none yet.  Mathematica audit: none yet." There is no `\stagefield{Output}` boxed equation on the card; the bottom-line claim is the *status* itself — that after the self-matched susceptibility closure (Stage 140), the mouth-gain ambiguity has narrowed to a finite-bias regular canonical mouth branch plus a definite gain pair (M_s, M_q), no longer a free profile or free gain pair. The notes file carries the substantive content as a carry-forward summary of Stages 137–140: the symbolic mouth gains M_s = L g_s²/(K_s Θ_σ) and M_q = −L(K_s g_q − λ g_s)²/(K_s(K_s K_q + λ²)Θ_σ); the normalized relation M_q = −M_s(𝔤_c − 𝔯)²/(1+𝔯²); the core-balance family M_q = −M_s/4, M_s = 4Σ_m, M_q = −Σ_m; the Family-1 numerics (natural M_s ≈ 1.66854, M_q ≈ −0.24270; canonical M_s ≈ 1.80594, M_q ≈ −0.45149); and the susceptibility closure M_s = (20/9) T̂_m² with a ~4% natural-vs-canonical traction difference. The notes are explicit that these are *carried forward* from upstream and that two microscopic questions (which branch the real mouth source picks; the actual T̂_m) remain open. The card's `\stagefield{Checks}` even instructs that numerical fixed points "are recorded as numerically located, not closed-form constants," consistent with a status card that owns no new derivation.

## What the script claims to verify

There is no script of any kind for this unit (no `.py`, no `.wl`, no saved output). The card itself declares no audit was run. Therefore the unit performs no in-stage verification; all of its substantive content is carried forward from upstream stages 133/137/138/139/140, each of which has both a SymPy and a Mathematica audit plus committed outputs.

## Paper ↔ script cross-check

There are no script-side checks to map, so every paper-side deliverable maps to an upstream carry-forward. The manifest marks this unit `is_status_only_candidate: True` and `is_checkpoint: False`. Per the status-only handling in the prompt, a `missing_sympy`/`missing_mathematica` finding is only valid if the card/notes reference a result that no upstream unit's scripts actually verify. I confirmed the opposite: every carried value is anchored in an upstream committed output (table below). The card is also honest about its own status ("none yet"), so there is no paper↔script disagreement to flag.

| Paper-side deliverable (notes/card) | Status | Upstream anchor |
|---|---|---|
| Symbolic M_s, M_q closed forms (notes item 1) | match | stage137 sympy output :4–5 |
| M_q = −M_s(𝔤_c−𝔯)²/(1+𝔯²) (notes item 2) | match | stage138 sympy output :10 |
| M_q = −M_s/4 core-balance family (notes item 3) | match | stage138 sympy output :12–13 (R_q = 1/4 both branches) |
| Family-1 natural/canonical numerics (notes item 4) | match | stage139 sympy output :4–5, :8–9 |
| M_s = (20/9) T̂_m²; ~4% difference (notes item 5) | match | stage140 sympy output :2 and :5 |

Set `paper_alignment: aligned`.

## Assertion inventory

No scripts exist; there are no assertions to inventory.

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | (none) | — | (no script for this unit) | n/a | n/a |

## Findings

None. This is a legitimate status-only / no-engine consolidation card. The card's `\stagefield{Verification}` truthfully reports "SymPy audit: none yet. Mathematica audit: none yet," and every value it consolidates is verified by an upstream stage's committed scripts (137/138/139/140), each of which carries both engines. Per the prompt's status-only rule, the absent engines are not findings, and there is no paper↔carry-forward mismatch.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit, so there is no transliteration question. (The upstream `.wl` scripts that actually verify the carried values — stage137–140 — are out of scope for this unit's audit.)

## Engine cross-check

Not applicable — neither engine is present for this unit.

## Verdict justification

`clean`. I read the paper card, the notes, and the part-04 appendix line for this unit before looking for scripts. The unit is a status/consolidation card that, by its own explicit statement, runs no SymPy or Mathematica audit. The attacks I tried and that failed to produce a finding: (1) I searched for any script naming or referencing 141 (`find` + grep across `scripts/` and `mathematica/`) — none exists, so this is genuinely no-engine, not a mislaid file; (2) I checked whether the card/notes carry forward any value that the upstream chain does NOT actually verify (the only thing that would make a missing-engine finding valid here) — every one of the five notes deliverables traces to a committed upstream output to printed precision; (3) I checked the card for an unsupported `\stagefield{Output}` claim or a status tag that overstates verification — it instead honestly tags "none yet" and labels the numerics as "numerically located, not closed-form constants." Nothing overstates what upstream proves, so there is no `paper_misalignment`.

## Self-test notes

No new check is being prescribed, so the variable-independence, symmetry/parity, and trivial-case traps do not apply (no `diff`/integral/assertion authored). I did exercise the status-only-specific trap: I verified the carry-forward chain is real by confirming each notes value against the upstream committed outputs (137 :4–5, 138 :10/:12–13, 139 :4–5/:8–9, 140 :2/:5), so the "missing engine" is the legitimate no-engine status rather than a gap in coverage. Path check is moot — no directive is written.

## Value Reconciliation (pass-2 augmentation)

This unit emits **no values of its own** (no `.py`/`.wl`/`.txt`). Per the augmentation's status-only guidance, I reconcile the values the card/notes *report* (which are carried forward from upstream) against (a) the stage-141 docs and (b) the upstream committed outputs that are their authoritative source. No upstream value is run; all are read from committed `.txt`.

| value | source (upstream py output line) | .tex/.md location | status |
|---|---|---|---|
| M_s = L g_s²/(K_s Θ_σ) | stage137 sympy output :4 (`M_s = L*g_s**2/(K_s*Theta)`) | notes :8–10 (item 1) | MATCH |
| M_q = −L(K_s g_q − λ g_s)²/(K_s(K_s K_q+λ²)Θ_σ) | stage137 sympy output :5 | notes :8–13 (item 1) | MATCH |
| M_q = −M_s(𝔤_c−𝔯)²/(1+𝔯²) | stage138 sympy output :10 (`-Sigma_0*(g_c - r)**2/(r**2 + 1)`) | notes :17–20 (item 2) | MATCH |
| M_q = −M_s/4 (R_q = 1/4 both branches) | stage138 sympy output :12–13 | notes :23–29 (item 3) | MATCH |
| M_s = 4Σ_m, M_q = −Σ_m | stage138 :9 (`M_s normalized = Sigma_0`) + :10 (consistent rename Σ_m = Σ_0/4) | notes :23–29 (item 3) | MATCH |
| M_s^nat ≈ 1.66854 | stage139 sympy output :4 (`1.66854252965…`) | notes :34 | MATCH |
| M_q^nat ≈ −0.24270 | stage139 sympy output :5 (`-0.242696939…`) | notes :35 | MATCH |
| M_s^comp ≈ 1.80594 | stage139 sympy output :8 (`1.80594111…`) | notes :40 | MATCH |
| M_q^comp ≈ −0.45149 | stage139 sympy output :9 (`-0.451485277…`) | notes :41 | MATCH |
| M_s = (20/9) T̂_m² | stage140 sympy output :2 (`Sigma_0 in terms of That = 20*That**2/9`) | notes :47 (item 5) | MATCH |
| ~4% natural-vs-canonical traction difference | stage140 sympy output :5 (`fractional traction enhancement = 0.040358816…`) | notes :48–51 (item 5) | MATCH |

INTERNAL items (no finding): none — this unit emits no scaffolding because it has no script. The upstream pass/fail flags, residuals, and That_nat/That_comp intermediates belong to the upstream stages' own audits, not to unit 141.

reconciliation: complete; 11 values checked, 0 misaligned
