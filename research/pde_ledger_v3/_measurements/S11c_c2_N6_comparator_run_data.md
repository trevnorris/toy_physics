# S11c-c2 N6 cross-engine comparator — RUN DATA (raw; ⛔ disposition NOT yet adjudicated)

⚠ This is the RAW comparator output picture. Per [[feedback_reconcile_representational_bridge]] a NONZERO residual is
⛔ NOT yet a disagreement; the disposition (reconcilable representational difference vs genuine cross-engine difference)
is the reconcile/step-record stage. ⛔ Nothing here is adjudicated.

## Reproduce
`python3 research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py` (defaults: the committed WL `.out`
`ae73b884` + the 3 committed SymPy `.out`s `7c0790ab`). Ran ~2.86 h, peak RSS ~330 MB, output ~1.6 GB — ⛔ NOT
committed, REPRODUCIBLE by re-running the command. The mechanical tally is COMMITTED alongside this record at
`_measurements/S11c_c2_N6_comparator_run_tally.txt`. ⚠ session `/tmp` scratchpad artifacts (the raw `.out`, the
`n6cmp_tally.py` helper) are EPHEMERAL — ⛔ do not expect them across a compact; re-run to regenerate.

## RUN_ACCOUNTING (clean)
`{families: 34, families_with_join: 15, families_with_unpaired: 34, parse_failed: 0, deferred_oversize: 0,
zero_extract_failures: 0, peak_rss_kib: 330412, runtime_s: 10294}`. ⭐ 0 deferrals — ALL N6 operands were tractable
on this box (contrary to the ≥64 GB worry; N6's largest ~80 MB); nothing forced, nothing operationally failed.

## SYMBOLIC channel (CASE A_minus_B) — matched keys only; UNDECIDED = unmatched typed key (pairing table)
⭐ **DIRECT cross-engine AGREEMENT (matched keys → residual 0, 0 nonzero):**
- `N6COV_R_COV`, `N6COV_R_COV_BASELINE`, `N6COV_R_COV_CONTROL_DELTA` — 160 ZERO / 0 NONZERO each. ⇒ the per-engine
  **R_cov = 0 (Reading B, source-naturality) is CROSS-ENGINE CONFIRMED** on matched keys.
- `N6COV_SOURCE_CONTROL_DELTA` — 160 ZERO / 0.
- `N6RC_CARRIER_BRIDGE_RESIDUAL` — 320 ZERO / 0. ⇒ the **geometric carrier reconciles cross-engine (C_E=C_M)**.
- `N6RC_ADVECTION_ABSENCE` 6/0; `N6RC_FROZEN_RELATIONS` 8/0.

⚠ **SURFACED computed residuals (⛔ disposition deferred to reconcile — rep-difference vs genuine?):**
- `N6RC_CARRIER_EULERIAN` / `CARRIER_MATERIAL` — 280 ZERO / **40 NONZERO** each. The blind-WL graph-geometry carrier vs
  the imported-slab SymPy carrier — the DO-NOT-FOLD representational residual. ⚠ The carrier BRIDGE residual vanishes on
  matched keys, but that does ⛔ NOT dispose of these 40 cross-engine operand residuals; whether they reconcile
  representationally or genuinely differ is UNADJUDICATED (the reconcile stage's job).
- `N6COV_SOURCE_ACTUAL` / `SOURCE_BASELINE` / `SOURCE_PREDICTED` — 16 ZERO / **76 NONZERO** each; `N6COV_FROZEN_PHI`
  42 ZERO / 18 NONZERO. The constitutive SOURCE channel + the field map Φ — consistent with the per-engine finding
  that the residual localizes to the constitutive source (source bridge nonzero in 3 cases). ⚠ Whether the two
  engines' sources AGREE (rep-difference) or genuinely differ is THE open cross-engine question for the reconcile.

**ENTIRELY UNMATCHED cross-engine (different WL block/kernel/component vocabularies):**
`N6RC_R_N6`, `SPLIT_CHECK`, `SPLIT_SUM`, `EULERIAN/MATERIAL_OPERAND`, `*_CHANNEL`, `CROSS_CHANNEL`, `DIMENSIONS`, all
guards — all UNDECIDED (unmatched typed keys → pairing table), with **NO matched symbolic sibling AND no structural
support line** (the support-tally 400/0 is entirely SOURCE+CARRIER, above; these families contribute none). ⇒ there
is **NO direct cross-engine comparison of R_N6 itself** — the N6 cross-engine evidence rests on R_cov + the carrier
bridge + the SOURCE/CARRIER support, ⛔ not on R_N6.

## STRUCTURAL channel — support AGREES where matched
- support: 400 ZERO / **0 NONZERO** / 1680 UNDECIDED ⇒ the nonzero-SUPPORT agrees on every matched key (0 support
  disagreements); the rest are one-sided/unmatched (honest UNDECIDED). SOURCE_ACTUAL/BASELINE/PREDICTED 96 ZERO
  support each; CARRIER_EULERIAN/MATERIAL 56 ZERO each.
- dimension_vectors: honestly surfaced (zero-vectors agree; `STRUCTURE_DISAGREE` where shapes differ). ⚠ the tally's
  "NONZERO" bucket here miscounts zero-VECTORS (e.g. `[0,0,0]`) as nonzero — the review legs verified this channel is
  agree-or-honestly-surfaced; `dimension_record` is always KEY_DISAGREE (redundant) and `addition_consistency` always
  BooleanNotResidualable (required bool rejection) — both NON-load-bearing (comparator review record).

## NEXT — the reconcile DISPOSITION (⛔ not done here)
Frame + adjudicate whether the SURFACED residuals collapse under justified representational identities (deeper
agreement) or are genuine cross-engine differences: (1) the blind-vs-imported CARRIER rep (40); (2) the constitutive
SOURCE channel (76) + Φ (18). ⛔ Per the CAS-authorship rule any collapse-testing instrument is CODEX-written +
G1-reviewed, ⛔ never orchestrator-authored; and the disposition gets legs (the c1 reconcile's correction-verify scoped
the orchestrator's verdict TWICE — [[feedback_reconcile_representational_bridge]]). Then the c2 step record (surface,
⛔ not pre-adjudicate, these + the carried 2 S11c-b signs / 6 §3d / c1 ENERGY; carry the 3 N6 premise caveats; fix
`I_{M→E}`; NO per-substep card).
