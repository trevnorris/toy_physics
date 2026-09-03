# Deferred heavy runs — re-run these when a bigger box is available

Some in-band engine self-checks OOM the current **30 GB** dev box (the un-reduced variable-coefficient
operators are ~50k+ LeafCount and several are held/expanded at once). Where that happens we run a memory-fitting
variant in-band and prove the skipped/scoped checks **out-of-band** (independent build-leg ablations + numerics),
recording the provenance. This file lists exactly what to re-run, and how, on a machine with more RAM to get the
full-fidelity in-band run. Everything here is **already verified out-of-band** — this is belt-and-suspenders
in-band self-certification, not open physics.

Guidance for the re-run box: the operators peak ~15–16 GB per engine on the current controls; give it **≥64 GB**
of headroom to be safe. Watch RSS; a Mathematica kernel can balloon and orphan — `kill -9` by exact pid if so.

> ⚠ **SCOPE — what the ≥64 GB requirement DOES and DOES NOT bound (measured 2026-09-03,
> `directives/_measurements/S11c_b_step0_residual_scope.md`).** The ≥64 GB requirement is for the **full 4-case
> in-band `.out` regen with the tower-depth control variants** (`operatorTruncated/Extended`,
> `kernelTruncated/Extended`, built unconditionally at `…mathematica_audit.wl:2204-2231`) — those variants are the
> ~16 GB/case peak. It does **NOT** bound the cross-engine **residual**, which needs only the PRIMARY single-case
> operator + `FINAL_KERNEL`: a guarded single-case WL run of exactly those (no tower variants) peaked **7.95 GB RSS
> / 0.99 GB in-kernel, min MemAvailable 14.94 GB, ~26 min** on this 30 GB box. ⇒ the cross-engine single-case
> residual is **doable on the 30 GB box**; only a comparator-PARSEABLE single-case `.out` needs a residual-mode
> emit (engine's own `emitShared` on the primary objects, tower variants skipped).

---

## PY — S11c-b #89 (SymPy engine): the PRIMARIES_ONLY-skipped in-band controls

- **Engine:** `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
- **What's skipped in the committed `.out`:** it was generated with `S11CB_PRIMARIES_ONLY=1`, which skips the
  build-heavy in-band controls — `PROJECTION_EQUIVALENCE`, `HESSIAN_FREEZE`, and the `FORM` ablation (see the
  guard at `…sympy_audit.py:4334`, `if not os.environ.get("S11CB_PRIMARIES_ONLY")`).
- **Re-run command (full in-band):**
  ```
  cd /var/projects/toy_physics
  python3 research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py \
    > research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out 2>&1
  # (do NOT set S11CB_PRIMARIES_ONLY)
  ```
  ⚠ On the 30 GB box this OOMs ~3.5 h in, still inside the first control (the MATERIAL coupling-kernel
  `restrict_expression` builds dominate).
- **What it would add / verify in-band:** the basis 40→26 form ablation, the lever-C projection equivalence,
  and the Hessian-freeze control — all currently proven **out-of-band** by both build legs' ablations/numerics
  (`_measurements/S11c_b_89_sympy_buildleg_clearance.md`). A green in-band run is confirmation, not new physics.

---

## PY — S11c-b #90 (SymPy engine): enlarged face-coupling all-case output and heavy controls

- **Engine:** `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
- **What's deferred:** regeneration of the full four-case `.out` and the build-heavy §5a/§5b, tower-depth,
  uniform-limit, and homogeneity control payloads after folding the supplied face generalized-force rows into
  the operator. The #90 build directive explicitly excludes the full package loop on this box.
- **Completed targeted computation:** one requested-depth coupling case, a reduced-depth origin-schema smoke,
  and independent source-feed removal/form-ablation probes. These print the folded operator rows, the complete
  weak blocks, construction origins, adjointness operands, and literal source-ablation differences without
  using a residual value as an acceptance condition.
- **Re-run command (full in-band, on a ≥64 GB box):**
  ```
  cd /var/projects/toy_physics
  python3 research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py \
    > research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out 2>&1
  # (do NOT set S11CB_PRIMARIES_ONLY)
  ```
- **What it would add:** the all-case computed primary objects and every affected control's operand A, operand
  B, and `A−B`. Those residual values are emitted findings, not build acceptance gates.

---

## WL — S11c-b #89b (Wolfram engine): the ENTIRE in-band `.out` regen (operator emit + the two controls)

- **Engine:** `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
- **Status:** DEFERRED — and the deferral is BROADER than first thought. Measured 2026-09-02 (`wolframscript`
  run with `S11CB_SKIP_HEAVY_CONTROLS=1`, memory-watched): even with **both heavy equivalence controls gated
  off**, the run reached **16.1 GB RSS while still in the per-case operator loop at only 13 (pre-operator-emit)
  tags** and was reaped under memory pressure (NOT an OS OOM — dmesg clean; the orphaned kernel was killed by
  exact pid). Cause is intrinsic to the un-freeze physics itself, ⛔ not the double-checking: the correct
  un-frozen operator must hold its **full background jet tower un-reduced (~16 GB/case)** right up to the
  **final** `applyBackgroundProfileRetainingJets` reduction — reducing any earlier is the very re-freeze #89b
  removes (rule 17). Per-case streaming (`directives/S11c_b_89b_wl_percase_streaming.md`; the streaming behavior — `KeyDrop` of `KERNEL_SOURCE_*` + `spoolAssociationValue` — is in the committed engine) frees
  between cases, but a single case's un-reduced operator is ~16 GB, so producing all four cases' operator emits
  does not fit 30 GB. ⇒ the **whole** in-band `.out` regen (the operator/kernel/μ_θ emits AND the two deferred
  controls below) is bigger-box work; on 30 GB the committed engine ships WITHOUT a fresh full `.out` and the
  #89b build legs verify by source-fidelity + reduced-scale ablation + the out-of-band records.
- **The two controls additionally gated** behind `S11CB_SKIP_HEAVY_CONTROLS` (directive
  `directives/S11c_b_89b_wl_defer_heavy_controls.md`), each emitting a `DEFERRED_HEAVY_CONTROL` marker:
  - `CONTROL_SURFACE_BACKGROUND_JET_DIFFERENCE_ATOMS` (per-surface Hessian witness)
  - `CONTROL_TRACTABLE_ACTIVATION_EQUIVALENCE` (full-operator equality vs per-row top-down reference)
- **Re-run command (full in-band, on a ≥64 GB box):**
  ```
  cd /var/projects/toy_physics
  wolframscript -file research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl \
    > research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out 2> /tmp/wl.stderr
  # (do NOT set S11CB_SKIP_HEAVY_CONTROLS; no timeout; .out is git-annex — datalad unlock or it's a plain file)
  ```
  ⚠ The operator loop alone peaks ~16 GB/case on 30 GB (reaped there); needs ≥64 GB for the full 4-case run.
- **What it verifies in-band (on the bigger box):** (a) the emitted operator/kernel for all four cases retains
  the retained-grade curvature (the whole un-freeze), and (b) the tractability speed-ups (innermost-first `Div`
  activation + per-summand `Series`) drop no jet term — full-operator equality vs the per-row top-down reference
  on **every** case. Both are already proven OUT-OF-BAND: the un-freeze/tractability approach by the diagnosis
  agent + 2 tractability decision legs (`PossibleZeroQ`), and the engine's implementation by the #89b build legs
  (source-fidelity + reduced-scale ablation). A green in-band all-cases run is confirmation, not new physics.
