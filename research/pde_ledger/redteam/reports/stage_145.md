---
unit_id: 145
batch: IV.5
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage145_mouth_branch_selection_status.md
  paper_appendix: present
---

# Audit unit 145 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_145.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage145_mouth_branch_selection_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 145 is a **status-only ledger entry** documenting the current state of the
mouth-branch selection problem inside the explicit positive-exponential
mouth-layer reduction. The card's purpose line states it is "a coupled mouth
fixed point and gain selection ledger step. Its audit target is the verification
output quoted below," where the quoted output is verbatim:

> Branch choice is settled inside the explicit closure; finite corrections
> around \((\Pi_*,\widehat T_{m,*})\) remain.

The notes elaborate five "now explicit" items already derived in upstream
stages 125–144 — the bias-factor formula \(\mathfrak g_\Pi\), the Family-1 gain
law and self-matched closure \(\Sigma_0 = (20/9)\widehat T_m^2\), the
impossibility of the upper compensated branch \(\mathfrak g_+^{F1}>1\), the fact
that the naive equal-normalized branch \(\mathfrak g_c=1\) is reached only as
\(\Pi\to\infty\), and the numerically-located unique finite regular branch at
\(\Pi_*\approx 1.50882951349316,\; \widehat T_m(\Pi_*)\approx 0.901484054174205\).
Two open items are explicitly enumerated: (1) accuracy of the explicit family as
a proxy for the actual moving-throat layer, and (2) the size of finite
corrections around \((\Pi_*,\widehat T_{m,*})\). The card's `Verification` field
declares plainly: "SymPy audit: none yet. Mathematica audit: none yet." The
part-04 appendix row (line 30) corroborates the block role: "Stages 133--145:
coupled mouth fixed point, core-to-mouth gain map, and canonical regular branch
selection."

## What the script claims to verify

No scripts exist for this stage (`scripts/.../stage145_..._sympy_audit.py` and
`mathematica/.../stage145_..._mathematica_audit.wl` are both absent, as are their
output files). The paper card itself acknowledges this state under
`\stagefield{Verification}`. This is a `is_status_only_candidate: True` unit per
the audit prompt manifest entry; the carve-out for status-only units applies.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| "Branch choice is settled inside the explicit closure" (declarative status) | (no script — carry-forward from upstream stages 133–144) | n/a (status-only) |
| "finite corrections around \((\Pi_*,\widehat T_{m,*})\) remain" (declarative open item) | (no script — this is an open-items note, not a verification target) | n/a (status-only) |
| Checks list: gain pair \((M_s,M_q)\), self-matched susceptibility closure, numerically-located fixed points carried as numerical (procedural reminders) | (no script — reminders to downstream consumers) | n/a (status-only) |
| Notes items 1–5 (\(\mathfrak g_\Pi\) formula, gain law, upper-branch impossibility, naive-branch limit, numerical \(\Pi_*\)) | (no script in this stage — these are carry-forward from earlier stages) | carry-forward (status-only) |

The paper card and notes both consistently frame this stage as a **status
record**, not a fresh derivation. The card's explicit "SymPy audit: none yet.
Mathematica audit: none yet." is in agreement with the absence of scripts on
disk. There is no on-stage assertion the scripts would be expected to discharge.

`paper_alignment` is set to `aligned`: the paper card, notes, and part-04
appendix row are mutually consistent in describing this as a coupled
fixed-point / gain selection status ledger; the absence of scripts is a
declared state, not a contradiction.

## Assertion inventory

No assertions exist (no scripts on disk).

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | — | — | — | — | — |

## Findings

None.

Per the prompt's status-only carve-out: a `missing_sympy` / `missing_mathematica`
finding would only be valid if the carry-forward chain references results no
upstream stage's scripts verify. The notes' five "now explicit" items are all
attributable to the upstream block (stages 125–144 per the part-04 row):
- Bias-factor \(\mathfrak g_\Pi\) closed form (positive-exponential mouth-layer
  family) — block stages 125–132 per appendix row, line 29.
- Family-1 gain law \(\Pi = \Sigma_0[1 - R_q(\Pi)\mathcal S_q(\Pi)]\) and
  \(\Sigma_0 = (20/9)\widehat T_m^2\) — block stages 133–144.
- Upper-branch impossibility \(\mathfrak g_+^{F1}>1\) — same block.
- Naive equal-normalized limit \(\Pi\to\infty\) — same block.
- Numerical \(\Pi_*,\widehat T_{m,*}\) — declared `\StatusNumerical{}` by the
  part-04 claim-status banner (line 37).

I did not read upstream stages' scripts (the prompt forbids reading other
stages' files), so I cannot independently confirm those upstream scripts each
discharge their share. What I can confirm is that **this stage's own card
declares the result as carry-forward and labels the unresolved part
`\StatusOpen{}`**, which is an honest declaration consistent with status-only
status. No script-side claim is silently introduced beyond the paper card.

`paper_misalignment` does not apply: the paper card, notes, and appendix
agree on framing.

## Independent-derivation check (Mathematica)

n/a — no Mathematica script exists.

## Engine cross-check

n/a — no scripts exist for either engine.

## Verdict justification

`clean`. The stage is a declared status-only ledger entry. The paper card
explicitly states "SymPy audit: none yet. Mathematica audit: none yet." and
flags `\StatusOpen{}` for the realization question. The notes file enumerates
which results are carry-forward and which remain open; both are consistent with
the part-04 appendix row describing stages 133–145 as the coupled mouth fixed
point / gain selection block. Under the prompt's status-only carve-out, the
absence of scripts is not a finding because (a) the carry-forward results
referenced (g_Pi formula, gain law, Sigma_0 closure, branch ordering, numerical
fixed point) are all attributed to upstream block stages 125–144 by the
appendix's own enumeration, and (b) there is no script-side assertion to
audit against — the stage introduces no new verifiable identity of its own.
The attacks I attempted: (1) check whether the notes smuggle in a fresh
identity the .tex doesn't cover — no, all five items are explicitly framed as
already-derived block results; (2) check whether the part-04 appendix row
attributes any new derivation to stage 145 alone — no, it bundles 133–145; (3)
check whether the quoted Output ("Branch choice is settled... finite
corrections remain") implies a numerical bound the scripts would have to
discharge — no, it is a categorical claim plus a declared open item; (4) check
whether the numerical values in the notes (\(\Pi_*\approx 1.508...\),
\(\widehat T_m\approx 0.901...\)) are introduced fresh here vs. carried — the
notes' framing ("What is now explicit") and the `\StatusNumerical{}` banner on
the part-04 claim-status line confirm carry-forward. No paper↔notes
disagreement, no script-introduced claims, no stale anything.

## Self-test notes

- Status-only carve-out: verified the missing scripts are declared in the
  paper card itself (`Verification: SymPy audit: none yet. Mathematica audit:
  none yet.`) and that the carry-forward items map to the upstream block per
  appendix row line 30.
- Paper round-trip: re-checked card body, notes, and appendix row — all three
  describe stage 145 as a coupled-fixed-point / gain-selection status ledger
  with the same Output sentence and the same open-items framing. No
  `paper_misalignment` introduced or latent.
- No directive needed (`findings_count: 0`, `verdict: clean`).
