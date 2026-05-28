---
unit_id: 128
batch: IV.4
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
  notes_stage_files: [moving_throat_pde_stage128_mouth_source_bias_status.md]
  paper_appendix: present
---

# Audit unit 128 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_128.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage128_mouth_source_bias_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 128 is an explicitly status-only ledger entry in the "Positive mouth-source and boundary layer" block (Part IV, anchor `MTDC-T8`). The card declares `\claimstatus{\StatusExactClosure{} / \StatusOpen{}}` and `\stagefield{Verification}{SymPy audit: none yet.  Mathematica audit: none yet.}`, so by its own statement the stage does not itself run a verification script. The body claim, in narrative form, is: branch sign for the Family-1 outlet/core compensation is fixed by positivity of the localized mouth source — the upper compensated branch \(\mathfrak g_+^{F1}\approx 2.798\) is ruled out, and the lower compensated branch \(\mathfrak g_-^{F1}\approx 0.758\) is the unique physically admissible canonical branch; the open datum is the actual positive mouth-source profile \(\sigma(z)\). The card's `\stagefield{Checks}` list three explicit deliverables: (i) positivity of the mouth source (no signed fitting to the upper branch), (ii) zero-flux and boundary-layer normalizations in the GNLS/localized-Maxwell reduction, and (iii) the Family-1 compensation point against the lower branch (not the singular equal-normalized branch). The notes restate the same bound \(0\le\mathfrak g[\sigma]\le1\) using \(\mathfrak g[\sigma]=\int_0^L\sigma(z)\cos(\pi z/(2L))\,dz\), and quote the explicit-source values \(\mathfrak g=1\) (point), \(\mathfrak g=\pi/4\) (self-matched D/N), and the convex/slab/exponential broadening parameters \(\xi_*\approx 0.18392\), \(x_*^{\rm slab}\approx 0.7978\), \(x_*^{\exp}\approx 0.6628\). The boxed bottom line: "Under positive localized mouth sourcing, the lower compensated branch is uniquely selected and is easily reachable." The appendix's Section "Core-to-mouth compensation law" (around lines 603-672) is where these formulas live formally (`eq:app-part04-gsigma-def`, `eq:app-part04-positive-source-bound`, `eq:app-part04-g-match`, `eq:app-part04-gPi`, `eq:app-part04-Pi-star`).

## What the script claims to verify

There is no SymPy or Mathematica script for unit 128. The card itself states "SymPy audit: none yet. Mathematica audit: none yet." The unit is a status / consolidation card whose verifications are carried by the surrounding stages 125-132 (which constitute the "positive mouth-source theorem and explicit GNLS/localized-Maxwell mouth boundary layer" block per appendix line 29). No script-side claim exists to evaluate.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Positivity of mouth source / branch selection (lower vs upper) | (none for unit 128 — relies on upstream block) | carry-forward |
| Zero-flux and BL normalizations (GNLS/localized-Maxwell) | (none for unit 128 — relies on upstream block) | carry-forward |
| Family-1 compensation point vs lower branch | (none for unit 128 — relies on upstream block) | carry-forward |

Per the manifest, this unit is `is_status_only_candidate: True`. The card transparently declares its own verification status as "none yet" and is by intent a positioning / ledger card, not a stand-alone derivation unit. Setting `paper_alignment: aligned` because no script-side check contradicts the card; the carry-forward block (stages 125-132) is the locus of the underlying derivations and is audited under those units' own reports, not here.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | n/a | n/a | (no script) | n/a | n/a |

No assertions exist; both SymPy and Mathematica files are absent and the card declares the verification status as "none yet."

## Findings

(None.)

Under the status-only carve-out, a `missing_sympy` / `missing_mathematica` finding would only be valid if the carry-forward chain did not support the card's claims. The claims in this card (positivity bound on \(\mathfrak g[\sigma]\), exclusion of the upper Family-1 branch, admissibility of the lower branch, GNLS/localized-Maxwell BL normalizations, Family-1 compensation point) are exactly the equations expressed in the part-IV appendix sections "Positive mouth-source selection" and "Electrochemical boundary layer" (lines 606-672 of `stage_appendix_part04.tex`) and in the matching notes file. The block declares those derivations belong to stages 125-132 collectively; auditing whether the corresponding upstream units have scripts that exercise those identities is the job of those units' own audits, not this unit's. No script-side `paper_misalignment` is possible because there is no script-side claim to misalign.

## Independent-derivation check (Mathematica)

Not applicable: no Mathematica script exists.

## Engine cross-check

Not applicable: no scripts exist.

## Verdict justification

`clean`. I read the paper card, notes, and the relevant appendix section in full. The card is an explicit status / consolidation entry — it self-declares `Verification: SymPy audit: none yet. Mathematica audit: none yet.` and tags itself `StatusExactClosure / StatusOpen`. Its body claims are reproduced in the appendix's "Core-to-mouth compensation law" section with formal equations, and in the notes file with the same numerical values (\(\mathfrak g_-^{F1}\approx 0.758035\), \(\mathfrak g_+^{F1}\approx 2.79795\), \(\mathfrak g_{\rm match}=\pi/4\), \(\xi_*\approx 0.18392\), \(\Pi_*\approx 1.5088\)). Under the unit's `is_status_only_candidate: True` carve-out, the absence of stage-128-specific scripts is not a finding so long as the carry-forward block (stages 125-132) supports the claim, and the appendix structure indicates it does. I attacked: (a) is there a script-side numerical value disagreeing with the card? — no script exists; (b) does the card claim a value that the appendix/notes do not anchor? — no, every quoted value (\(\mathfrak g_{\rm match}=\pi/4\), the two branch values, \(\Pi_*\), the broadening fractions) has a corresponding appendix equation or notes entry; (c) does the card silently claim something the upstream block cannot deliver? — no, the listed `Checks` (positivity, zero-flux BL normalization, Family-1 compensation point) all map onto explicit appendix equations in lines 609-672. No directive is written.

## Self-test notes

Traps checked: (1) status-only carve-out applicability — confirmed by both manifest flag `is_status_only_candidate: True` and the card's own "Verification: none yet" tag; (2) paper-side numeric anchoring — every constant in the card/notes (\(\pi/4\), \(\xi_*\), \(\Pi_*\), \(\mathfrak g_\pm^{F1}\)) is reproduced in the appendix equations; (3) silent paper claim outside upstream support — the card's three `Checks` items each correspond to appendix equations within the same Part IV block, so the carry-forward chain is intact. No script self-test required (no script exists).
