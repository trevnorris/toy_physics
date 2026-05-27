---
unit_id: 093
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: present
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage093_grouped_p2_status_update.md"]
  paper_appendix: present
---

# Audit unit 093 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_093.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage093_grouped_p2_status_update.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the row `\input{stages/stage_093}` at line 1220; appendix is long but the surrounding "Stage-by-stage verification cards" intro at lines 1207-1214 confirms the role)
- sympy: (missing — intentionally absent, status-only candidate per MANIFEST)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl`
- sympy output: (missing)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.txt`

## What the paper claims

Stage 093 is explicitly framed as a geometry-lane firewall ledger / status-update step. Per `\stagefield{Purpose}`: "Stage~093 is a geometry-lane firewall ledger step.  Its audit target is the verification output quoted below." The body identifies the forced conservative carrier `\widehat Y_Q^{\rm cons}=3/4+(1/4)/(1-\omega^2/\Omega_Q^2)` and the obstruction variables `(\epsilon_2,\epsilon_4)`. The card's quoted Output declares "the remaining question: whether actual isotropic geometry is dynamically inert or supplies (\epsilon_2,\epsilon_4)." The notes give the carry-forward static-limit derivations: `rho_alpha = 4/3`, `zeta_req = 1/3`, and the obstruction formula `c_pole = (1+eps_4)/[4(1+eps_2)^2]` (notes lines 28–31, 52–53). The verification checklist (Checks) names three items, of which only the first ("static limit eps_2=eps_4=0 returns c_pole=1/4") is a directly scriptable numeric check; items 2 and 3 are hygiene/orthogonality status statements that reference earlier units rather than introduce new content to verify here.

## What the script claims to verify

The Mathematica script sets `eps2 = 0`, `eps4 = 0`, computes `cPole = (1+eps4)/(4*(1+eps2)^2)`, then derives `cGeom = 1 - cPole`, `rhoAlpha = 1/cGeom`, `zetaReq = cPole/cGeom`, and asserts via `expectZero` that the four quantities equal `1/4, 3/4, 4/3, 1/3` respectively (lines 28–43). The script's role is to confirm that the obstruction formula imported from upstream collapses to the four canonical static-limit constants quoted by the paper card and notes. This is a recap / carry-forward check, consistent with `is_status_only_candidate: true`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Static limit `eps_2 = eps_4 = 0 ⇒ c_pole = 1/4` (Checks item 1; notes §1) | `expectZero["c_pole - 1/4", cPole - 1/4]` (wl:40) | match |
| `c_geom = 3/4` (implicit; module decomposition in notes §1) | `expectZero["c_geom - 3/4", cGeom - 3/4]` (wl:41) | match |
| `rho_alpha = 4/3` (notes §1) | `expectZero["rho_alpha - 4/3", rhoAlpha - 4/3]` (wl:42) | match |
| `zeta_req = 1/3` (notes §1) | `expectZero["zeta_req - 1/3", zetaReq - 1/3]` (wl:43) | match |
| Obstruction formula `c_pole = (1+eps_4)/[4(1+eps_2)^2]` (notes §2) | Encoded literally at wl:30 but only evaluated at `eps2=eps4=0` | partial (carry-forward — full formula derivation belongs to upstream Stage 75 per notes line 53) |
| `l=0` and `l=2` orthogonality (Checks item 2) | none | missing (status-only carry; checked upstream) |
| Minimal-module hypothesis still carries (Checks item 3) | none | missing (meta-hygiene statement, not a numeric check) |
| Full omega-dependent module `Yhat_Q^cons(omega) = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` (notes §1, card body) | none — only the omega-independent reduction is exercised | partial (omega-dependent pole structure is a recap from upstream; this unit's scope is the static-limit consolidation) |

Dominant pattern: scriptable Checks items align; remaining mismatches reflect the status-only role rather than misalignment.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | mathematica | 40 | `expectZero["c_pole - 1/4", cPole - 1/4]` | static-limit Checks item 1 | yes |
| A2 | mathematica | 41 | `expectZero["c_geom - 3/4", cGeom - 3/4]` | notes §1 (3/4 + 1/4 split) | yes |
| A3 | mathematica | 42 | `expectZero["rho_alpha - 4/3", rhoAlpha - 4/3]` | notes §1 (rho_alpha = 4/3) | yes |
| A4 | mathematica | 43 | `expectZero["zeta_req - 1/3", zetaReq - 1/3]` | notes §1 (zeta_req = 1/3) | yes |

Each assertion is non-tautological in the sense that a typo in the obstruction formula at wl:30 (e.g. missing the square on `(1+eps2)`, wrong placement of `eps4`) would cause one or more of A1–A4 to fail. The four checks are not independent — A2/A3/A4 are algebraic consequences of A1 — but A1 alone provides the discriminating signal, and A2–A4 serve as redundancy on the downstream symbolic algebra in the script.

## Findings

(no findings)

The reading order I followed: paper card → notes → appendix row → Mathematica script → output. Three concerns were considered and rejected:

1. **Potential `hardcoded_result` / `tautological_check` on A1.** The obstruction formula at wl:30 is imported literally from upstream (Stage 75 per notes line 53). Evaluating at `eps2 = eps4 = 0` reduces by hand to `1/4`, so the assertion's only failure mode is a typo in the formula. I considered flagging this as weak verification, but the stage card explicitly requests precisely this check ("Check the static limit eps_2=eps_4=0 returns c_pole=1/4"), and the unit's MANIFEST flag `is_status_only_candidate: true` makes this carry-forward role legitimate. The verification depth matches the paper's stated scope for this stage.

2. **Potential `script_missing_paper_claim` on Checks items 2 and 3 and on the omega-dependent module.** Items 2 ("l=0 / l=2 orthogonality") and 3 ("minimal-module hypothesis still carries") in the card's verificationchecklist are not numeric tests but status carry-forwards from upstream units; flagging them would conflict with the status-only carve-out. The omega-dependent module `Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` likewise originates upstream (notes §1 marks it as recovered from the grouped-`P2`+static-geometry split) and is not re-derived here.

3. **Banner mislabel.** Line 26 prints `banner["STAGE 076 — GROUPED-P2 STATUS UPDATE"]` while the stage number is 093 (and the paper section header is "Stage~110"). This is a cosmetic numbering legacy (notes line 4 references "Stages 74-75"; paper card header says Stage~110; file/MANIFEST id is 093). The terminal line at script line 46 correctly prints "Stage 093 Mathematica audit passed." No mathematical content depends on this label, and it does not fall into any of the ten finding categories. Noted for informational purposes only.

## Upstream carry-forward dependencies

Per the status-only carve-out rules, I note the result this unit imports without re-verifying:

- The obstruction formula `c_pole = (1+eps_4)/[4(1+eps_2)^2]` (used at wl:30) is sourced upstream from Stage 75 per notes line 53 ("Stage-75 obstruction formula"). I did not read other units' scripts (per scope), so I cannot here confirm that Stage 75's audit verified this formula's derivation; this is a known and intended carry-forward, not a finding under the carve-out.
- The full omega-dependent module form (notes §1) is similarly imported from the grouped-`P2` + static-geometry split done upstream.

## Independent-derivation check (Mathematica)

Not applicable — no SymPy script exists for this unit (status-only carve-out). The Mathematica side stands alone.

## Engine cross-check

Not applicable — single-engine status unit. The Mathematica output at lines 13–24 of the saved transcript shows all four `expectZero` checks passing with residual `0`.

## Verdict justification

The Mathematica script faithfully exercises the only directly-scriptable item on the stage card's verificationchecklist (static limit `eps_2=eps_4=0 ⇒ c_pole = 1/4`) and the four derived constants the notes record (`c_pole=1/4`, `c_geom=3/4`, `rho_alpha=4/3`, `zeta_req=1/3`). Items 2 and 3 of the Checks list are hygiene/carry-forward statements that do not require a new numeric check under the status-only carve-out (`is_status_only_candidate: true`). The obstruction formula at wl:30 is imported from upstream (Stage 75), which is legitimate carry-forward; the audit cannot probe upstream from within this scope. Output file is fresh (mtime 2026-05-11 13:03 > script mtime 2026-05-11 11:56). Verdict: clean.

Attacks attempted that failed: (a) checked whether `cPole - 1/4 == 0` could be algebraically forced by the formula's structure (it can only collapse to 1/4 if both factors `1+eps4` and `4(1+eps2)^2` are typed correctly — typos would flip the assertion); (b) checked whether `cGeom`, `rhoAlpha`, `zetaReq` are redundant (they are derived from cPole, so failure of A1 would propagate, but no independent signal — accepted because the stage is a recap); (c) checked that paper card constants (`c_pole=1/4`, `rho_alpha=4/3`, `zeta_req=1/3`) match the script literals (they do, exactly); (d) compared the obstruction formula's wl:30 form `(1+eps4)/(4*(1+eps2)^2)` against notes line 53 `(1+eps_4)/[4(1+eps_2)^2]` (match, including the square placement).

## Self-test notes

Checked: (1) variable independence — N/A, no derivatives in the script; (2) symmetry/parity — N/A, no integrals; (3) trivial-case pre-check — confirmed by hand that with `eps2=eps4=0`, the formula `(1+eps4)/(4*(1+eps2)^2)` reduces to `1/4`, then `cGeom = 3/4`, `rhoAlpha = 4/3`, `zeta_req = 1/3` — all four `expectZero` arguments become `0` as expected; (4) path specifications — N/A, no missing-script finding; (5) paper round-trip — re-verified the four numeric constants (`1/4`, `3/4`, `4/3`, `1/3`) and the obstruction formula appear in the notes and card exactly as the script encodes them. No corrections required.
