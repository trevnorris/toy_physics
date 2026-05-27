---
unit_id: 103
batch: IV.2
auditor_model: claude-opus-4-7-1m
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage103_reduced_25pn_conditional_closure.md
  paper_appendix: present
---

# Audit unit 103 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_103.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage103_reduced_25pn_conditional_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read intro, retarded-factorization section, canonical DtN section, theorem block, and the `\input{stages/stage_103}` line; this stage has no standalone narrative beyond the imported card)
- sympy: (missing — manifest declares `is_status_only_candidate: true`, no path registered)
- mathematica: (missing — manifest declares `is_status_only_candidate: true`, no path registered)
- sympy output: (missing)
- mathematica output: (missing)

Manifest entry (`/var/projects/toy_physics/research/pde_ledger/redteam/MANIFEST.yaml:3600`): `is_checkpoint: false`, `is_status_only_candidate: true`, both `sympy.path` and `mathematica.path` are `null`.

## What the paper claims

The stage card (`stage_103.tex:7`) calls Stage 103 "a retarded \(2.5\)PN factorization ledger step" whose Derivation-ledger line is: "The computation isolates the reduced product \(\widehat m_0^{\,2}\chi_QN_Q=1\) and the canonical condition \(\chi_Q=1\)." The card's bottom-line quote (`stage_103.tex:16`) reads: "Reduced theorem closes iff the actual passive/outgoing branch has \(\chi_Q=1\)." The card explicitly states (`stage_103.tex:11`) "SymPy audit: none yet. Mathematica audit: none yet." — i.e., this card carries no verification of its own. The notes file (`moving_throat_pde_stage103_reduced_25pn_conditional_closure.md`) makes the same point: the conditional closure rests on six previously fixed items (geometry-clean conservative branch, 3/4+1/4 split, `rho_alpha=4/3`, Family-1 support sufficiency, `mhat_0=1+O(a^2/r^2)`, irrelevance of higher odd data) and reduces the remaining PDE-facing problem to a single scalar question about `\chi_Q`. The Part-IV theorem block (`stage_appendix_part04.tex:340-343`) restates the closure as Theorem `thm:app-part04-conditional-25pn-closure`, and the part-IV intro's compression chain (`stage_appendix_part04.tex:44-62`) lists `\chi_Q \to \Delta_Q` as one arrow in a larger reduction sequence.

Distinct deliverables enumerated by the card + notes:

1. **D1 — Factorization**: the point-particle 2.5PN normalization factorizes as `\widehat m_0^{\,2}\chi_Q N_Q = 1`, keeping source/conservative/outgoing factors separate (card Checks item 1; appendix eq. `eq:app-part04-main-factorization`).
2. **D2 — Source-map limit**: on the natural source-map branch, `\widehat m_0 = 1 + O(a^2/r^2)`, so the remaining obstruction is `\Delta_Q := \chi_Q - 1`.
3. **D3 — Higher-odd irrelevance**: higher odd terms begin beyond the point-particle 2.5PN coefficient (card Checks item 2; appendix line 295).
4. **D4 — Canonical outgoing DtN fingerprint**: the canonical outgoing `l=2` DtN branch has `\chi_Q = 1` (card Checks item 3; appendix eq. `eq:app-part04-chiQ-equals-one`).
5. **D5 — Conditional closure**: combining D1-D4 yields the conditional reduced 2.5PN closure.

Each of D1-D5 is itself proved in a neighbouring stage (see "Verdict justification" below).

## What the script claims to verify

No script exists for unit 103. The unit is registered as a status-only consolidation card; both `sympy.path` and `mathematica.path` are `null` in the manifest, and the stage card explicitly disclaims any standalone audit ("SymPy audit: none yet. Mathematica audit: none yet."). The card therefore carries no script-side claim of its own — the claim it consolidates is the conjunction of the deliverables proved in stages 100, 101, 102, 104, 105, and 106.

## Paper ↔ script cross-check

Because unit 103 is a status-only card with no script of its own, the cross-check is whether each paper-side deliverable D1-D5 is verified somewhere in the carry-forward chain. Visible upstream/neighbouring scripts named in the `scripts/` directory:

- `moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py` — by filename targets D1 (factorization `\widehat m_0^{\,2}\chi_Q N_Q = 1`).
- `moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py` — by filename targets D2 (`\widehat m_0 \to 1`).
- `moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py` — by filename targets D3 (higher-odd irrelevance).
- `moving_throat_pde_stage104_outgoing_dtn_fingerprint_sympy_audit.py` — by filename targets the DtN fingerprint needed for D4.
- `moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py` — by filename targets D4 (`\chi_Q = 1`).
- `moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py` — by filename targets D5 (canonical outgoing reduced closure).

| Deliverable | Paper anchor | Carry-forward verifier (filename evidence) | Status |
|---|---|---|---|
| D1 — `\widehat m_0^{\,2}\chi_Q N_Q = 1` | `stage_103.tex:13`, appendix eq. `eq:app-part04-main-factorization` | stage 100 sympy/mathematica scripts present | `match` (consolidation) |
| D2 — `\widehat m_0 = 1 + O(a^2/r^2)` | notes Statement bullet 5; appendix line 290 | stage 101 scripts present | `match` (consolidation) |
| D3 — higher-odd irrelevance | `stage_103.tex:23` Checks item 2; appendix line 295 | stage 102 scripts present | `match` (consolidation) |
| D4 — canonical `\chi_Q = 1` | `stage_103.tex:16`, appendix eq. `eq:app-part04-chiQ-equals-one` | stage 104 (DtN fingerprint) + stage 105 (chiQ fix) scripts present | `match` (consolidation) |
| D5 — conditional reduced closure | `stage_103.tex:16`, appendix Thm. `thm:app-part04-conditional-25pn-closure` | stage 106 scripts present | `match` (consolidation) |

Paper alignment: `aligned`. The card's status-only stance is honest (it explicitly states "none yet" for SymPy/Mathematica), and every load-bearing deliverable maps to a neighbouring unit whose scripts exist on disk. No carry-forward target is orphaned; no upstream script is missing. Per the audit instructions for `is_status_only_candidate: true` units, this is the configuration in which "missing_sympy / missing_mathematica is **not** a finding."

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | (none — status-only card) | — | — | — | — |

No assertions exist for unit 103. Inventory is empty by design; the v2 audit instructions explicitly permit this for status-only consolidations as long as the carry-forward chain is intact (it is).

## Findings

(none)

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists for unit 103, and none is required (status-only card; manifest `is_status_only_candidate: true`; no carry-forward references an unverified upstream result).

## Engine cross-check

Not applicable — neither engine is present for this unit. The deliverables it consolidates are each independently exercised in their own neighbouring units (stages 100, 101, 102, 104, 105, 106), which are audited separately.

## Verdict justification

Stage 103 is a status-only consolidation card whose paper claim is explicitly conditional ("Reduced theorem closes iff … `\chi_Q = 1`") and whose card body openly declares no standalone SymPy or Mathematica audit. The manifest confirms `is_status_only_candidate: true` with both script paths `null`. Per the v2 instructions, the only way this can be a finding is if the carry-forward references a result no upstream unit verifies — but every load-bearing piece (factorization, source-map limit, higher-odd irrelevance, canonical `\chi_Q = 1`, conditional closure) maps to a neighbouring unit whose scripts are present on disk (100, 101, 102, 104, 105, 106). The paper card, the notes file, and the Part-IV appendix narrative agree on the deliverables and on the consolidation role. No paper↔script mismatch is possible because there is no script side; the consolidation is bookkeeping. Adversarial attacks attempted: (a) check whether the card silently introduces a *new* numeric constant absent from upstream — no, every symbol (`\widehat m_0`, `\chi_Q`, `N_Q`, `\Delta_Q`, `\rho_\alpha`, the 3/4+1/4 split) is sourced upstream; (b) check whether the card's "iff" claim asserts more than the upstream chain proves — it does not, the upstream chain proves conditional closure modulo `\chi_Q = 1`, which is exactly what the card states; (c) check whether the appendix theorem on lines 340-343 quietly adds the natural source-map limit as a hypothesis the card omitted — the card includes it implicitly via Inputs (line 9: "the source-map factor `\widehat m_0`") and the notes file enumerates it as one of the six fixed items. No finding. Verdict: `clean`.

## Self-test notes

Traps checked: (1) verified the `is_status_only_candidate: true` flag in MANIFEST.yaml:3608 to confirm absent scripts are not by themselves a finding; (2) cross-walked every distinct deliverable in the card/notes/appendix theorem to a neighbouring unit whose script files exist on disk (stages 100, 101, 102, 104, 105, 106) so the carry-forward chain has no orphans; (3) verified the card's "iff" wording matches the appendix theorem's "if … closes; if not, the failure is `\Delta_Q = \chi_Q - 1`" framing — both sides are consistent.
