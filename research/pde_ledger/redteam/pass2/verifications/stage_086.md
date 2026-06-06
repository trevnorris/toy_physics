---
unit_id: 086
batch: III.5
verdict: verified
overall: verified
material_change: false
verified_at: 2026-06-05
---

# Verification — unit 086 (pass 2)

**Verdict: verified — clean (0 findings).** Audit verdict was `clean`; no directive, no Codex.
Stage 086 (`family1_loading_ratio_window`), not a checkpoint, both engines present.

The load-bearing physics is genuinely exercised: the symbolic unblocked reduction
`Q(zeta;0)=1+zeta` and the two derivative closed forms `dQ/dzeta=(1-eps)/(1-eps·zeta)^2` and
`d rho_max/d eps = zeta_max(zeta_max-1)/(1-eps·zeta_max)^2` (both re-derived by hand, signs
confirmed). The Mathematica side is NOT a transliteration — it adds independent symbolic
assertions (the two derivative closed forms, plus zeta-vs-paper input checks) that SymPy only
prints/consumes. The numeric `rho` threshold checks are weak-but-valid corroboration (round-trip
of `1+zeta`) — the thresholds are derived upstream (Stage 080 transcendental solves) and
legitimately carried; not tautological. Self-labels canonical ("STAGE 086").

**Value reconciliation:** 13 deliverable values checked, 0 misaligned (the `rho_suff^(J)`
envelope and the `eps_blk` cap live in the `.md` notes, not the terse `.tex` card → MATCH per the
augmentation guard).

**Deferred (not a finding):** SymPy py:37 comment "carried from Stages 63-64" is a stale
CROSS-reference (notes attribute these to Stages 080/081; +17 pre-renumber drift) → routed to
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed, never offset-sweep).

Orchestrator independent re-run (reliability gate): SymPy exit 0, Mathematica exit 0; both
committed outputs **byte-identical to HEAD**. Engines agree to 30 digits. No paper/notes edits.
