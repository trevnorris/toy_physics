# S11c-b cross-engine comparator — build review, two legs, consolidated (INSTRUMENT SOUND)

**Artifact:** `scripts/S11c_b_cross_engine_comparator.py` (Codex-built, `bqj3ftjqc`, exit 0) + synthetic
fixtures `scripts/test_S11c_b_cross_engine_comparator.py`. **Brief:**
`directives/S11c_b_comparator_build_directive.md`. **Legs:** fresh Claude agent + Grok (Codex-written ⇒
fresh-Claude + Grok), rendered prompt `_legs/S11c_b_comparator_build_review.md` (byte-identical to both).
**Raw:** Grok `~/.s11_build/S11c_b_comparator_review_grok.txt` + `/tmp/S11c_b_comparator_review/` (ablation
`.py`/`.stdout`, `REVIEW_REPORT.md`); Claude-agent
`/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/`
(`baseline_tests.out`, `ablation_a..d.out`, `fold4_probe.out`, `faceflux_probe{,2}.py`,
`real_run_full.out`, `real_run_status.txt`). Verdicts are the orchestrator's after adjudication (rule 13).

## Deliverable verification (before legs)
16/16 synthetic fixtures pass (`python3 test_S11c_b_cross_engine_comparator.py` → `OK`). Build task exit 0;
builder report plausible (28 families, all join, 1079s run). Memory healthy at launch (20Gi avail).

## Load-bearing folds — ALL SOUND (both legs, FORM-ablated, literal diffs saved)
- **(a) μ_θ per-branch registry** — forcing the forbidden global collapse `mu_theta_L/M → mu_theta` FAILs
  `test_mu_registry_is_branch_specific_and_argument_preserving` +
  `test_mu_heads_are_never_globally_collapsed` (LAB−MATERIAL residual → **0** under the collapse). The real
  path keeps `mu_theta_L`(LAB)/`mu_theta_M`(MATERIAL), args preserved. Branch-consulted. SOUND.
- **(b) ENERGY_BASIS compared RAW** — `reduce_divergence=True` appears ONLY inside `extract_coupling`
  (comparator L872/880/909/921); `extract_energy` (L736–757) never sets it. Injecting a divergence reducer
  on a first-background-jet difference (`W_bg·θ_d1 + W_bg_d1·θ` vs `0`) masks it to `0`; the real path
  retains the `W_bg_d1` jet term. No representative-identity fold. SOUND (⛔ never blanket-collapse honored).
- **(c) SECTOR rename load-bearing** — PY uses `*_THICKNESS` spellings, WL uses `*_THETA_EW_UL` (disjoint
  sets, 8 vs 36 each); only the rename joins them. Baseline COUPLING join=24; breaking the rename →
  join=8, py_only=16, wl_only=16 (all 16 sector-bearing cases unpair). SOUND.
- **(d) Residual is an independent typed recursion** — one-sided PY corruption moves `operand_A` + `A−B`,
  leaves `operand_B` fixed (symmetric for WL). Not zero-by-construction. No `assert` in the comparator
  (only token is in a docstring, L560). SOUND.
- **Fold-4 adjointness reducer** — (i) a true total-divergence density reduces to `0`; (ii) a non-adjoint
  `Λ_X`-face density survives nonzero. The reducer zeros only genuine total divergences; it cannot mask the
  dissipative `Λ_X` non-self-adjointness nor misread an adjoint density. SOUND.
- **Hand re-key (NOT via `make_key`)** — independent keying of MU_THETA(4)/ENERGY_BASIS(2)/
  ADMISSIBILITY(20 = 12 body-force + 8 faces)/COUPLING sectors matched the comparator joins; no
  hand-found collision the comparator failed to flag as `duplicate_key`.

## Real run — two independent runs AGREE TO THE TOKEN
```
RUN_ACCOUNTING families=28 families_with_join=28 families_with_unpaired=7 parse_failed=0 duplicate_key=0
```
(Grok `runtime_seconds=1202.5`; Claude-agent `1192.6`.) 555 CASE lines, each `operand_A`+`operand_B`+`A−B`
emitted before any guard; **no `PASS`/`FAIL`/`VERDICT`/`FINAL_STATUS`/`EXPECTED` token** (output or source);
`return 2` only on operational error, else `return 0`. No family `join=0`. The 7 unpaired families (108
unpaired cases, each an emitted `UNDEFINED_UNJOINED` + `py_only/wl_only/axis_set_mismatch` note, 0 duplicate,
0 parse_failed): `CONTROL_INDEPENDENCE_{BASE,CORRUPTED,RESIDUAL}` wl_only=8; `COUPLING_KERNEL_TERM_ORIGINS`
py_only=4/wl_only=16; `DIMENSIONS` wl_only=4; `SLAB_OPERATOR` py_only=28/wl_only=16;
`SLAB_OPERATOR_TERM_ORIGINS` py_only=4/wl_only=12. These are genuine engine structural asymmetries (real
residuals) surfaced for reconcile, **not** comparator bugs.

## The one finding — FACE_FLUX (Grok flagged; verified + downgraded, rule 13)
- **Structural fact (confirmed by orchestrator + both legs):** PY `COUPLING_KERNEL` VALUE has 4 subtrees
  `{SECTOR_LABELS, BULK_BLOCKS, FACE_FLUX_BLOCKS, VARIATIONAL_ADJOINTNESS}`; `extract_coupling` reads ONLY
  `VARIATIONAL_ADJOINTNESS` (comparator L856). So `FACE_FLUX_BLOCKS` (sub-keys `SOURCE_SUBSTRATE`,
  `TRANSVERSE_TO_FACE_FLUX`, `MU_THETA_TO_TRANSVERSE`, `ADVECTIVE_MASS_BLOCK`) produces **zero cases in the
  main COUPLING_KERNEL family**, whose `py_only=0` does not reflect it. WL emits **no** face-flux operand
  (0 occurrences), so no cross-engine partner exists.
- **Why it is NOT a silent drop:** the identical `FACE_FLUX` object is re-emitted by the engine into
  `COUPLING_KERNEL_TERM_ORIGINS` (engine L2064, key `FACE_FLUX`), and `extract_term_origins` (comparator
  L809) iterates the full origins map → it surfaces as **`py_only`, 4 cases** (LAB/MATERIAL × RHO4/RHOBR),
  `operand_B=<MISSING>`, `A−B=UNDEFINED_UNJOINED`, `case_note=py_only`. **Counted, visible.** It
  manufactures no false agreement (never subtracted to 0) and hides no cross-engine disagreement (WL has no
  counterpart; any boundary-representation difference in WL would surface as a nonzero adjointness residual
  in-family). Under the physics filter this is **not** a comparison-integrity finding. Both legs converge on
  this reading (Grok initially called it a silent drop; the term-origins visibility downgrades it).

## OWED → carried to the reconcile + step record (legibility + physics, NON-blocking for the instrument)
1. **Legibility:** the main `COUPLING_KERNEL` `py_only=0` line does not advertise that 3 of its 4 PY
   subtrees (`SECTOR_LABELS`/`BULK_BLOCKS`/`FACE_FLUX_BLOCKS`) go unread in-family, and that its bulk is
   compared via the `VARIATIONAL_ADJOINTNESS` sector-pairings rather than the top-level `BULK_BLOCKS`. Not
   an integrity defect (WL has no `BULK_BLOCKS`; the sector pairings ARE compared cross-engine, 24 joins).
2. **Reconcile obligation (hand-code):** (a) PY `FACE_FLUX_BLOCKS` is py_only — determine whether WL is
   genuinely missing the boundary/interface coupling or folds it into its weak-form adjointness (§3c weak
   variational restriction, PY-explicit vs WL-weak — a never-blanket-collapse hand-check); (b) confirm PY
   top-level `BULK_BLOCKS` ≡ the `VARIATIONAL_ADJOINTNESS` forward sector-pairings across ALL cases (Grok
   verified equality only for one sample), so nothing bulk goes uncompared.
3. Inherited: S11c-a's owed control-family keying; SymPy dead code (`paired_kernel_from_density`/
   `mixed_variation`). Admissibility §5 control-coverage (from the SymPy repair review).

## Verdict
Comparator PASSES both build legs; the instrument is SOUND (folds ablation-verified, run reproducible and
token-identical across two independent runs, no verdict/target/false-agreement). No comparator repair
warranted. NEXT: commit the reviewed artifact → run the committed comparator → hand-coded reconcile
(`S11c_b_handcoded_comparison.py`, matching the S11c-a pattern) honoring the OWED reconcile obligations →
step record.
