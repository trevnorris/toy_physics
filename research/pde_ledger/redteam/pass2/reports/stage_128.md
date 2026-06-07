---
unit_id: 128
batch: IV.4
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage128_mouth_source_bias_status.md]
  paper_appendix: present
---

# Audit unit 128 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_128.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage128_mouth_source_bias_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1290: `\input{stages/stage_128}`; this part-level file is a pure aggregator — no per-stage prose row beyond the card itself)
- sympy: (missing) — confirmed: no `*128*` file under `scripts/`
- mathematica: (missing) — confirmed: no `*128*` file under `mathematica/`
- sympy output: (missing)
- mathematica output: (missing)

This unit is `is_status_only_candidate: true`, `is_checkpoint: false` per `redteam/pass2/MANIFEST.yaml:4426-4450` (all four file paths `null`, `exists: false`). Status-only handling applies.

## What the paper claims

Stage 128 is a **status / consolidation card**, not a fresh computation. Its `\stagefield{Purpose}` says it is "a positive mouth-source and boundary layer ledger step" whose "audit target is the verification output quoted below," and `\stagefield{Verification}` states verbatim: **"SymPy audit: none yet.  Mathematica audit: none yet."** The card's bottom line (boxed in the notes) is qualitative: *"Under positive localized mouth sourcing, the lower compensated branch is uniquely selected and is easily reachable."* The notes pin the supporting numeric carry-forwards: Family-1 geometry `r_F1 ≈ 1.77799353547498`; the cosine-moment bound `0 ≤ g ≤ 1` for any positive normalized localized mouth source; the two compensated branches `g_+^{F1} ≈ 2.79795199200529` (ruled out by positivity) and `g_-^{F1} ≈ 0.758035078944663` (unique admissible); and the explicit reach points `g(point)=1`, `g(self-matched D/N)=π/4≈0.785398` (within `3.61%` traction of exact compensation), `ξ_*≈0.183918` (convex family), `x_*^slab≈0.797839`, `x_*^exp≈0.662765`. The deliverable is the *status assertion* (branch sign fixed; open datum = actual positive mouth-source shape), and `\stagefield{Downstream use}` requires the status tag to be carried with the result into Stages 133–145.

## What the script claims to verify

There is no SymPy or Mathematica script for unit 128, and none is expected: the card itself declares "none yet" for both engines, and the manifest marks the unit status-only. Stage 128 emits no values of its own; it consolidates results derived and verified in upstream units. Therefore there is no script-side bottom-line assertion to attack — the audit reduces to (a) confirming the missing engines are legitimate per the status-only carve-out, and (b) confirming every value the card/notes carry forward is actually verified by an upstream unit's scripts and is faithfully transcribed.

## Paper ↔ script cross-check

The card's deliverable is a status statement backed by carry-forward values. Each carry-forward maps to an upstream unit whose scripts (both engines, where shown) verify it:

| Paper-side carry-forward (notes) | Upstream verifier | Status |
|---|---|---|
| `g_-^{F1} ≈ 0.758035078944663` (unique lower branch) | stage 125 sympy out L19 `g_- (numeric) = 0.75803507894466282692`; stage 144 (both engines) | match |
| `g_+^{F1} ≈ 2.79795199200529` (ruled out by positivity) | stage 125 sympy out L20 `g_+ (numeric) = 2.7979519920052934101` (+ `g_+ > 1 = True`); stage 144 | match |
| `0 ≤ g ≤ 1` cosine-moment bound | stage 125 (positive source theorem); `g_- < 1 = True`, `g_+ > 1 = True` | match |
| `r_F1 ≈ 1.77799353547498` | stage 121 (geometric r-selection, both engines); re-asserted as literal in stage 139 sympy L44 | match |
| `g(self-matched D/N) = π/4 ≈ 0.785398` | stage 126 (both engines) | match |
| `ξ_* ≈ 0.183918` (convex family) | stage 126 sympy out L22 `xi_* numeric = 0.18391840551153962831` (+ L23 `g_xi(xi_*) - g_- = 0`); both engines | match |
| `x_*^slab ≈ 0.797839` | stage 127 sympy out L9 `0.79783936090456...`; both engines | match |
| `x_*^exp ≈ 0.662765` | stage 127 sympy out L10 `0.66276540262316...`; both engines | match |
| `3.61%` traction gap (self-matched vs exact) | stage 134/140 outputs | match |

`paper_alignment: aligned` — every paper-side carry-forward has a faithful upstream verifier, and the card honestly declares its own no-script status.

## Assertion inventory

No script exists for this unit, so there are no assertions to inventory. The carry-forward chain is anchored as follows (these assertions live in the *named upstream units*, not in 128, and were not re-audited here beyond confirming the emitted values match the 128 notes):

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| (upstream) | sympy stage125 | out L19/L20 | `g_-`, `g_+` closed form + numeric, `g_+>1`, `g_-<1` | branch uniqueness (128 claim) | yes |
| (upstream) | sympy stage126 | out L22/L23 | `xi_* numeric`, `g_xi(xi_*) - g_- = 0` | convex-family reach (128 claim) | yes |
| (upstream) | sympy stage127 | out L9/L10/L11/L12 | `x_*^slab`, `x_*^exp` with residual ~1e-83 | penetration-family reach (128 claim) | yes |
| (upstream) | sympy stage139 | L44 | `assert abs(rF - 1.77799353547498) < tol` | r_F1 geometry (128 claim) | yes |
| (this unit, 128) | — | — | (no script) | — | n/a |

## Findings

None.

The missing SymPy and Mathematica engines do **not** constitute a `missing_verification_script` finding under the status-only rule: that finding is valid only "if the unit's scripts (or comments referencing source units) reference a result that no upstream unit's scripts actually verify." Every result the card/notes reference is in fact verified upstream (stages 121/125/126/127/144 etc., with both engines present for the load-bearing branch and penetration values). The card is also transparent — it states "SymPy audit: none yet. Mathematica audit: none yet." rather than implying a verification it does not have.

No `paper_misalignment`: the card and notes agree with each other, and every carried numeric matches the upstream output to all printed digits (see cross-check table and Value Reconciliation below).

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` for this unit. (Upstream load-bearing verifiers for the carried values, e.g. stages 125/126/127/144, do have Mathematica outputs present, so the carry-forward chain is dual-engine backed; auditing those scripts' independence is the job of their own audit units, not this one.)

## Engine cross-check

Not applicable — neither engine is present for unit 128, by design.

## Verdict justification

`clean`. I read the stage card, the notes, and the appendix line, then confirmed there is genuinely no script for this unit (full `find ... -name '*128*'` over `scripts/` and `mathematica/` returns nothing; manifest entry has all paths null). The attack I tried was to break the status-only carve-out by finding a carried result with no upstream verifier — every one of the eight distinct numeric carry-forwards (and the `0≤g≤1` bound and the qualitative branch-uniqueness claim) traces to an upstream unit whose committed output reproduces the value, with both engines present for the load-bearing branch (125/144), convex-family (126), and penetration-family (127) results. The card honestly reports its own unverified status and tags the result as a conditional ledger entry, matching the notes. Nothing to fix; no directive written.

## Value Reconciliation (pass-2 augmentation)

Unit 128 has **no script of its own** (status-only consolidation card), so it emits no values directly. Per the augmentation's status-only clause, I reconcile the values the card/notes *carry forward* against (a) the docs that report them and (b) the upstream committed outputs that actually produce them. No script was run; reconciliation is based on script source + committed `output/*.txt` + reasoning.

| value | source (upstream py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_F1 ≈ 1.77799353547498` | stage139 sympy L44 (literal assert) ← stage121 (both engines) | notes L11 `\mathfrak r_{F1}\approx 1.77799353547498` | MATCH |
| `0 ≤ g ≤ 1` (cosine-moment bound) | stage125 sympy out L11/L21–23 | notes L23 `0\le \mathfrak g\le 1`; card L22 checklist | MATCH |
| `g_+^{F1} ≈ 2.79795199200529` | stage125 sympy out L20 `2.7979519920052934101` | notes L27 `\mathfrak g_+^{F1}\approx 2.79795199200529` | MATCH |
| `g_-^{F1} ≈ 0.758035078944663` | stage125 sympy out L19 `0.75803507894466282692` | notes L33 `\mathfrak g_-^{F1}\approx 0.758035078944663` | MATCH |
| `g(point) = 1` | stage125 sympy out L11 (equal-normalized branch) | notes L41 `\mathfrak g=1` | MATCH |
| `g(self-matched D/N) = π/4 ≈ 0.785398` | stage126 (both engines) | notes L46–47 `\frac{\pi}{4}\approx 0.785398` | MATCH |
| `3.61%` traction gap | stage134/140 sympy out | notes L49 `3.61\%` | MATCH |
| `ξ_* ≈ 0.183918` (convex family) | stage126 sympy out L22 `0.18391840551153962831` (+ L23 residual = 0) | notes L57 `\xi_*\approx 0.183918` | MATCH |
| `x_*^slab ≈ 0.797839` | stage127 sympy out L9 `0.79783936090456...` | notes L62 `x_*^{\rm slab}\approx 0.797839` | MATCH |
| `x_*^exp ≈ 0.662765` | stage127 sympy out L10 `0.66276540262316...` | notes L64 `x_*^{\exp}\approx 0.662765` | MATCH |
| boxed status claim (lower branch uniquely selected & reachable) | stage125 sympy out L34–35 conclusion | notes L74–76 boxed; card L16 quote | MATCH |

Internal scaffolding (not expected in this card / no finding): none originate in unit 128 (no script). The upstream residual/pass-fail flags (e.g. stage127 `residual slab ≈ -2.67e-83`, stage126 `g_xi(xi_*) - g_- = 0`) are scaffolding of *their* units and reconcile within those audits, not here.

**reconciliation: complete; 11 carried deliverable values checked, 0 misaligned** (single-engine status: no engine present for unit 128 by design; carry-forward chain is dual-engine backed at the upstream sources for the load-bearing branch/penetration/convex-family values).

## Self-test notes

No directive is written (zero findings), so the directive self-test is vacuous, but I exercised the relevant traps anyway: (1) no `sp.diff`/`D[]` derivatives exist in this unit (no script), so the variable-independence trap is N/A; (2) the carried `g` is a cosine moment integral over `[0,L]` whose `0≤g≤1` bound is verified upstream at stage125, not re-asserted here; (3) trivial-case check applied as a value-match: every carried numeric was confirmed digit-for-digit against the upstream committed output (e.g. `g_- = 0.75803507894466...`, `ξ_* = 0.18391840551153...`, `x_*^slab = 0.79783936...`), all reproducing. The only live trap for a status-only unit — a carried result with no upstream verifier — was checked directly and did not fire.
