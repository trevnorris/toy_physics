---
unit_id: 145
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00-06:00
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage145_mouth_branch_selection_status.md]
  paper_appendix: present
---

# Audit unit 145 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_145.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage145_mouth_branch_selection_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read rows 17, 25-34, 86, 640-869, 1173-1180, 1324 referencing this unit / the Family-1 mouth branch)
- sympy: (missing) — confirmed absent via `find scripts -iname '*stage145*'` (no match) and MANIFEST `'145'.files.sympy.path: null, exists: false`
- mathematica: (missing) — confirmed absent via `find mathematica -iname '*stage145*'` (no match) and MANIFEST `'145'.files.mathematica.path: null, exists: false`
- sympy output: (missing) — MANIFEST `'145'.files.sympy_output.exists: false`
- mathematica output: (missing) — MANIFEST `'145'.files.mathematica_output.exists: false`

## What the paper claims

Stage 145 is a **status / derivation-ledger card** for "Mouth-Branch Selection." It consolidates the coupled-mouth fixed-point and gain-selection block, isolating the finite-bias **regular canonical mouth branch** and the gain pair `(M_s, M_q)` that replaces the naive equal-normalized limit. The card's `Verification` line states verbatim: "SymPy audit: none yet. Mathematica audit: none yet." (`stage_145.tex:11`), and its `Verification note` defers detail to the block narrative / source stage (`stage_145.tex:19`). There is **no `\stagefield{Output}` line**; the closest bottom-line statement is the quoted closure: "Branch choice is settled inside the explicit closure; finite corrections around `(\Pi_*,\widehat T_{m,*})` remain." (`stage_145.tex:16`). The card is explicitly labeled "a derivation ledger entry, not an unconditional actual-branch theorem" with a status tag that must be carried downstream (`stage_145.tex:27`).

The source notes enumerate five concrete deliverables this status summarizes: (1) the exact bias map `g_Π = 2Π(2Π e^Π+π)/((4Π²+π²)(e^Π−1))` with range `2/π < g_Π < 1` (`...stage145...md:8-20`); (2) the self-consistent Family-1 gain law `Π = Σ₀[1 − R_q(Π)S_q(Π)]` with `R_q(Π) = (g_Π − r_F1)²/(1 + r_F1²)` and self-matched closure `Σ₀ = (20/9) T̂ₘ²` (`md:22-32`); (3) the upper compensated branch is impossible because `g₊^F1 > 1` (`md:33-36`); (4) the naive equal-normalized branch `g_c = 1` is reached only as `Π→∞, T̂ₘ→∞` (`md:38-44`); (5) the lower compensated branch is the unique finite regular branch at `Π* ≈ 1.50882951349316`, `T̂ₘ(Π*) ≈ 0.901484054174205` (`md:46-52`). The part-04 appendix carries these as the block narrative (anchor MTDC-T8.4): `Π* ≈ 1.50882951349316` (`stage_appendix_part04.tex:663`), `g_Π` map (`:658`), the self-consistent law (`:763-766`), `Σ₀(Π*) ≈ 1.80594111095636` and `T̂ₘ(Π*) ≈ 0.901484054174205` (`:773-775`), and the "Unique regular Family-1 canonical mouth branch" theorem (`:787-790`).

## What the script claims to verify

Nothing — there is **no SymPy script and no Mathematica script** for this unit. The MANIFEST entry (`MANIFEST.yaml:4947-4971`) declares `is_status_only_candidate: true`, `is_checkpoint: false`, and all four file slots (`sympy`, `mathematica`, `sympy_output`, `mathematica_output`) as `path: null, exists: false`. A filesystem `find` over both `scripts/` and `mathematica/` for `*stage145*` returns no match. This is consistent with the card's own "none yet" verification declaration. The unit is therefore a no-engine consolidation/status card, the same shape as the notes-only units 103/113/120/124.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Bias map `g_Π` + range `2/π<g_Π<1` (md:8-20; appx:658,660,783) | none (no script) | n/a — status-only, no engine declared |
| Family-1 gain law + `Σ₀=(20/9)T̂ₘ²` (md:22-32; appx:751,763) | none | n/a |
| Upper branch impossible `g₊^F1>1` (md:33-36) | none | n/a |
| Naive `g_c=1` only as Π→∞ (md:38-44; appx:785) | none | n/a |
| Unique regular lower branch `Π*≈1.50882951…`, `T̂ₘ≈0.901484…` (md:46-52; appx:663,775) | none | n/a |

Every deliverable is a carry-forward consolidation of the part-04 mouth-branch block (anchor MTDC-T8.4, `stage_appendix_part04.tex:1178`), not a freshly-asserted result of this stage. The card declares no audit script ("none yet"), so no script-side check is required by the card itself, and no carry-forward references a result that lacks an upstream derivation in the visible block narrative. `paper_alignment: aligned` — paper card, notes, and appendix narrative agree, and the absence of scripts matches the card's own declaration.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | (none) | — | no scripts exist for this unit | — | — |

No assertions to inventory: this is a no-engine status-only card.

## Findings

None. (Per the status-only handling in the auditor contract: `is_status_only_candidate: true`, and the card explicitly declares "SymPy audit: none yet. Mathematica audit: none yet." A `missing_sympy`/`missing_mathematica` finding is only valid if a carry-forward references an upstream result no upstream script verifies — that is not the case here. The carried numerics `Π*`, `T̂ₘ`, `Σ₀`, and the `g_Π` map all live in the part-04 block narrative (MTDC-T8.4) that the upstream coupled-mouth-fixed-point stages (125–145) derive. No `paper_misalignment`: the card, notes, and appendix agree, and the value-reconciliation below confirms there is no script-emitted value left unreflected.)

## Independent-derivation check (Mathematica)

N/A — no `.wl` exists.

## Engine cross-check

N/A — neither engine is present; `engines_agree: n/a`.

## Verdict justification

**Clean.** I read the paper card, the source notes, and the part-04 appendix block (the natural carrier of this stage's narrative and numerics) before concluding. The card itself declares no audit script for this stage ("SymPy audit: none yet. Mathematica audit: none yet.", `stage_145.tex:11`), the MANIFEST flags the unit `is_status_only_candidate: true / is_checkpoint: false` with all four file slots null, and a filesystem search confirms no `stage145` script exists. I attacked the "missing-script is a finding" angle: it fails because the unit is a legitimate status-only consolidation whose every carried deliverable (bias map, gain law, branch-selection results, `Π*`/`T̂ₘ`/`Σ₀`) is derived upstream and reflected in the part-04 appendix (MTDC-T8.4); no carry-forward references an unverified upstream result. I attacked the paper-alignment angle: card, notes, and appendix narrative are mutually consistent (Π*, T̂ₘ, Σ₀, and the g_Π map all agree across notes and appendix to full quoted precision; the notes' tighter `2/π<g_Π<1` range and the appendix regularity theorem's `0<g_Π<1` are consistent statements, not a contradiction — the family's infimum is `2/π`, stated at `stage_appendix_part04.tex:660`). With no engine present, no assertion can be tautological, hardcoded, or mis-anchored. Verdict: clean, no directive.

## Self-test notes

I checked the four standard traps insofar as they apply: (1) variable-independence / derivative-zero traps — N/A, no derivatives in any script (no script). (2) Symmetry/parity of integrals — N/A, no integrals asserted in-stage. (3) Trivial-case pre-check on assertions — N/A, no assertions. (4) Path specifications for missing-script findings — N/A, no missing-script finding raised (status-only carve-out applies and is satisfied). I additionally cross-checked the doc-to-doc numeric consistency (Π*, T̂ₘ, Σ₀, g_Π range) between the notes and the appendix block and found exact agreement to quoted precision, so no `value_mismatch` arises even at the prose level.

## Value Reconciliation (pass-2 augmentation)

This is a **no-engine status-only stage**: no SymPy or Mathematica script, no saved output `.txt`. Therefore there are **zero script-emitted result values** to reconcile (no computed constants, no boxed symbolic results derived in-script, no labeled-result prints). Per the augmentation's status-only guard, I reconcile "whatever the present engine emits" — which is nothing — and note the no-engine status. The stage's numeric content is a carry-forward consolidation derived upstream and carried in the part-04 appendix block (anchor MTDC-T8.4); I confirm those carried numerics reconcile between the two prose carriers (notes `.md` and appendix `.tex`) for completeness, but none of them originates from a script for this unit, so none is a script→doc reconciliation target.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| (no script-emitted values) | — none (no engine) — | — | n/a — no-engine status-only stage |

Doc-carrier consistency confirmation (informational, not a script→doc check — no script exists):

| carried numeric | notes `.md` | appendix `.tex` | agree? |
|---|---|---|---|
| `Π* ≈ 1.50882951349316` | `...stage145...md:51` | `stage_appendix_part04.tex:663` | yes |
| `T̂ₘ(Π*) ≈ 0.901484054174205` | `...stage145...md:51` | `stage_appendix_part04.tex:775` | yes |
| `Σ₀ = (20/9) T̂ₘ²` (and `Σ₀(Π*) ≈ 1.80594111095636`) | `...stage145...md:31` | `stage_appendix_part04.tex:751, 773` | yes |
| `g_Π` exact bias map | `...stage145...md:12-16` | `stage_appendix_part04.tex:658` | yes |
| `g_Π` range `2/π < g_Π < 1` | `...stage145...md:18-20` | `stage_appendix_part04.tex:660, 783` | yes (appendix `0<g_Π<1` with infimum `2/π` at :660; consistent) |

INTERNAL items (genuine scaffolding, no finding): none (no script ⇒ no scaffolding to list).

reconciliation: complete; 0 script-emitted values checked, 0 misaligned (no-engine status-only stage; doc-carrier numerics independently confirmed consistent across notes and appendix).
