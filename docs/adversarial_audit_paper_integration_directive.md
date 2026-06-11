# Adversarial Audit — Paper Integration Directive (for Codex)

**Status:** Ready to execute. Written 2026-06-10 after the Codex review of `docs/adversarial_audit_deployment.md`.

**Scope:** The two paper-side edits that should land *now* — before the adversarial-audit skill is built — plus one definition-only item. Everything else from the Codex review (fit-insertion-point index appendix, parameter-value provenance ledger appendix, post-audit prose patching) is **deliberately deferred** until Phases A/B produce real artifacts; those appendices will be *generated* from `redteam_adversarial/` YAML, never hand-authored, to avoid the stale-label drift this project already paid for once (591 stale references in the numbering remediation).

**Ownership:** Codex applies these edits (paper files are Codex-owned per the calibrated convention); Claude reviews afterward. This directive states requirements and acceptance criteria only — wording is Codex's call.

**Path convention:** every `paper/...`, `stages/...`, `parts/...`, `appendices/...`, and `notes/...` path below is relative to `research/pde_ledger/` — e.g. `paper/README.md` means `research/pde_ledger/paper/README.md`, and `notes/stages/` means `research/pde_ledger/notes/stages/`. Only `docs/...` paths are repo-root-relative. `redteam_adversarial/...` (future audit artifacts) is also under `research/pde_ledger/`.

---

## Item 1 — Resolve the canonical-source contradiction (gates Phase B)

### Verified facts (do not re-derive; spot-check if desired)

- `paper/README.md:43` says `stages/`: "generated or archival stage templates; not the active narrative source for either build."
- `paper/README.md:46-51` ("Current source rule") says the canonical stage narrative lives in `appendices/stage_appendix_part*.tex` and `stages/` is "template inventory."
- `paper/appendices/reader_provenance_summary.tex:38-42` ("Current stage-source rule") says the opposite: the canonical narrative source is the per-stage TeX tree under `paper/stages/`, with the archive appendices merely assembling those files.
- Ground truth (verified 2026-06-10): the eight `appendices/stage_appendix_part0*.tex` files each carry substantial **inline theorem-block narrative** (e.g. part 06: 18 `\input`s but 1,360 lines), then a closing "Original-stage verification cards" section that **`\input`s the stage cards from `stages/`** — 253 inputs total across the eight parts, covering every stage.
- `stages/*.tex` is live, load-bearing content, not template inventory: both red-team passes audited it, and fixes land there (e.g. commit `b052471` corrected `stages/stage_100.tex`/`stage_101.tex`).

### Required outcome

Both documents describe the same, true, **two-layer structure**, in Codex's wording:

1. **Theorem-block synthesis layer** — inline narrative in `parts/` and in the archive stage appendices.
2. **Stage-card layer** — `paper/stages/*.tex`, the canonical per-stage provenance and audit anchors, assembled into the appendices via `\input`.

Additionally, the source-hierarchy description must name `notes/stages/` (and the other per-stage notes trees) as the **decisive audit layer beneath both** for provenance questions — the red-team record shows origin disputes are resolved by the notes, not by cards or scripts.

### Acceptance criteria

- [ ] `paper/README.md` no longer calls `stages/` "template inventory" or "not the active narrative source"; it states the two-layer structure and that all 253 cards are `\input` into the archive appendices.
- [ ] `reader_provenance_summary.tex` "Current stage-source rule" no longer implies the appendices are mere assembly; it names the inline theorem-block layer.
- [ ] The two documents agree with each other and with the ground truth above.
- [ ] The escalation path in `reader_provenance_summary.tex` includes the notes layer as the decisive provenance authority.
- [ ] No change to any stage card, appendix narrative, or `\input` list — this item edits *descriptions of* the structure only.

---

## Item 2 — Falsification-stack disclosure in the reader verification summary

### Requirement

`paper/appendices/reader_verification_summary.tex` currently reports audit *counts* (§"Current baseline") and repo *locations* (§"Verification layers") but does not distinguish the three falsification layers. A reader can mistake "verified" (internal consistency, complete) for "externally falsified" (fit-vs-derive, not yet run). Add a compact section — natural position: adjacent to "Verification layers" — presenting the falsification stack with honest, current status:

| Layer | Status (as of 2026-06-10) |
|---|---|
| Internal-consistency red-team (SymPy + Mathematica dual-engine, script audit) | Complete — 253/253, two full independent passes |
| Script→doc value reconciliation (second-pass augmentation) | Complete |
| Adversarial fit-vs-derive audit | Designed (`docs/adversarial_audit_deployment.md`); **not yet run** |
| Stage-1 branch realization (target-blind numerical PDE) | Gated on the adversarial layer (`docs/branch_realization_execution_plan.md`) |

Row content above is the substance; table shape and TeX idiom are Codex's call (match the existing `longtable` style).

### Acceptance criteria

- [ ] The section states plainly that internal-consistency verification **cannot** detect a fit dressed as a derivation, and that the adversarial layer exists for exactly that — one or two sentences, in the register of the existing claim-status firewall.
- [ ] Nothing in the new section overstates: the adversarial layer is *designed, not run*; branch realization is *gated, not started*.
- [ ] Reader build stays lean: status + escalation pointer only; no adversarial machinery, no per-stage tables.
- [ ] The existing "Current baseline" and "Verification layers" sections are not removed (light edits to avoid redundancy are fine).

---

## Item 3 — Audit-status vocabulary: define only, never hand-tag

### Requirement

Adjacent to the claim-status firewall (`paper/frontmatter/02_claim_status_firewall.tex`), add a short **definition-only** passage introducing the audit-status vocabulary, orthogonal to the mathematical claim-status tags:

`CAS-checked` · `dual-engine-checked` · `value-reconciled` · `value-provenance-checked` · `external-benchmark-sourced` · `fit-vs-derive-audited` · `branch-realization-tested`

It must state explicitly that per-stage audit-status *values* are **generated from machine-readable audit state** (red-team manifest, trackers, and the future `redteam_adversarial/` artifacts) and are **never hand-written into individual stage cards**. Hand-maintained inline tags are the exact stale-label failure mode the numbering remediation cleaned up.

### Acceptance criteria

- [ ] Vocabulary defined once, next to the firewall; presented as a second axis, not a replacement for the six mathematical status tags.
- [ ] The generated-never-hand-tagged rule is stated in the text.
- [ ] **No stage card receives an audit-status tag in this pass.** Zero edits under `stages/`.

---

## Out of scope (do not do in this pass)

- No fit-insertion-point index appendix, not even a stub table (deferred until Phase A emits `fit_insertion_points.yaml`).
- No parameter-value provenance appendix (deferred until Phase B emits the ledger).
- No edits to `docs/adversarial_audit_directive.md` (calibration confirmed; frozen).
- No per-stage prose patching ("derived" → "target condition" etc.) — that happens only after Phase C verdicts exist.
- No changes to scripts, notes, or the atlas graph.

## Verification

After applying: both paper builds (reader + archive) must compile clean. Run the standard build check and iterate until exit 0, per the usual contract.
