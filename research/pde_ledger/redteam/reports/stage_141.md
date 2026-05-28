---
unit_id: 141
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage141_mouth_gain_status.md
  paper_appendix: present
---

# Audit unit 141 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_141.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage141_mouth_gain_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_141}` line at L1316; no inline appendix prose for this stage)
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 141 is explicitly framed as a **status update / derivation-ledger entry**, not a freestanding verification. The card's `Purpose` line says it is "a coupled mouth fixed point and gain selection ledger step. Its audit target is the verification output quoted below." The `Verification` field states verbatim: "SymPy audit: none yet.  Mathematica audit: none yet." The quoted block records "the narrowed gain-scale problem after self-matched susceptibility." Three checklist items are listed (`Checks`): (i) check the gain pair `(M_s,M_q)` against outlet consistency; (ii) check the self-matched susceptibility closure before using the one-scalar branch law; (iii) check numerical fixed points are recorded as numerically located, not closed-form constants. The notes flesh out the substance: explicit throat-core formulas `M_s = L g_s^2/(K_s Theta_sigma)` and `M_q = -L (K_s g_q - lambda g_s)^2/(K_s(K_sK_q+lambda^2)Theta_sigma)`; the normalized-core form `M_q = -M_s (g_c - r)^2/(1+r^2)`; the core-balance compensation family `M_q = -M_s/4`, `M_s = 4 Sigma_m`, `M_q = -Sigma_m`; Family-1 numerical pairs (natural equal-normalized: `M_s ≈ 1.66854`, `M_q ≈ -0.24270`; canonical: `M_s ≈ 1.80594`, `M_q ≈ -0.45149`); and the self-matched closure relation `M_s = (20/9) hatT_m^2`. The notes' "What remains open" section reiterates this is a status carry-forward, not a new theorem. The `Downstream use` line states that the status tag must be carried with the result and that this card is "a derivation ledger entry, not an unconditional actual-branch theorem."

## What the script claims to verify

No SymPy or Mathematica audit script exists for unit 141. The paper card itself declares this state ("SymPy audit: none yet.  Mathematica audit: none yet."). There is therefore no script-side claim to audit.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Mouth-gain pair `(M_s, M_q)` formulas (notes #1) | — | upstream carry-forward (Stages 188–191 per notes) |
| Normalized-core form `M_q = -M_s (g_c-r)^2/(1+r^2)` (notes #2) | — | upstream carry-forward |
| Core-balance compensation family `M_q = -M_s/4`, `M_s = 4 Sigma_m`, `M_q = -Sigma_m` (notes #3) | — | upstream carry-forward |
| Family-1 numerical pairs (notes #4) | — | upstream carry-forward |
| Self-matched closure `M_s = (20/9) hatT_m^2` (notes #5) | — | upstream carry-forward |
| Status: ledger-entry summary of "narrowed gain-scale problem" | — | this stage is a prose ledger; no new identity to verify |

`paper_alignment = aligned`: the paper card explicitly says no verification scripts have been written (and none exist), and the stage card's purpose is to record status, not to prove a new identity. The notes confirm "Stage 141: Mouth-Gain Status Update" framing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | — | — | (no assertions; no scripts) | — | — |

## Findings

None.

### Status-only handling

Per the audit prompt's status-only carve-out: this unit's manifest entry has `is_status_only_candidate: True`. The paper card's `Verification` field is upfront that no SymPy or Mathematica audit has been written ("none yet"); the card's `Purpose` and `Downstream use` lines explicitly frame this as a derivation-ledger entry that consolidates results from upstream stages (Stages 188–191 per the notes, plus the self-matched susceptibility closure block). The carry-forward chain references results that are produced and verified upstream (mouth-gain formulas, core-balance compensation family, Family-1 numerical pairs, susceptibility closure relation). I did not (and per the audit's hard rules must not) read upstream stages' scripts to confirm those verifications, but the references named in the notes correspond to identified upstream stages, and the paper card does not claim a new identity that has no upstream anchor.

Therefore `missing_sympy` and `missing_mathematica` are not findings here. No `paper_misalignment` arises: the paper says "none yet" and there are indeed none — paper and script side are mutually consistent in their (shared) absence.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists.

## Engine cross-check

Not applicable — neither engine is present.

## Verdict justification

Clean. Stage 141 is a status-only ledger entry that openly declares it has no verification scripts ("SymPy audit: none yet.  Mathematica audit: none yet.") and the absence of both scripts and saved outputs is consistent with that declaration. The notes file matches the card in framing this stage as a status update consolidating upstream results (Stages 188–191, plus the self-matched susceptibility closure). The status-only manifest flag carves out the otherwise mandatory `missing_sympy` / `missing_mathematica` findings, and there is no `paper_misalignment` because the paper does not claim a new identity beyond what upstream stages produce. Attacks attempted: (i) check whether the paper card asserts a freestanding identity that would demand its own verification — it does not, the `Checks` items are pointer-style ledger checks against upstream content; (ii) check whether the notes claim a derivation original to this stage — they explicitly attribute the gain formulas to "the explicit throat-core plus mouth-layer closure of Stages 188–191"; (iii) check whether the card's `Inputs` or `Downstream use` lines hide a new constant the script would need to introduce — they list only carry-forwards and a downstream use tag.

## Self-test notes

I checked: (a) the paper card does not introduce a new identity requiring its own script; (b) the notes attribute the operational formulas to upstream stages (188–191) and to the self-matched susceptibility closure block; (c) the part-04 appendix only `\input`s the stage card and adds no extra prose claim. No trap applies because there is no script to mentally execute. The verdict is `clean` under the status-only carve-out.
