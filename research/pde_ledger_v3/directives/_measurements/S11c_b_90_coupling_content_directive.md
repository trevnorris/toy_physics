# S11c-b #90 — §3c coupling content: decision + build review record

Directive `directives/S11c_b_90_coupling_content_directive.md`. Target: `scripts/S11c_b_brane_operator_sympy_audit.py`
(SymPy engine) + a §0 clarity pin to `directives/S11c_b_SHARED_PHYSICS.md`. Repairs PY's coupling extraction to the
settled §3c INCLUDE/INCLUDE content: the reversible tilted-face geometry + the irreversible face response (Λ
symbolic), computed as face generalized-force rows INTO the §3b operator, then weak-restricted.

## Settled verdict (do not re-litigate)
§3c mandates INCLUDE/INCLUDE — WL spec-correct, PY under-extracts (bulk-only). Adversarially confirmed Codex+Grok
×2 rounds, user-endorsed. Records `directives/_measurements/S11c_b_coupling_84_{diagnosis,consult,basis_verification}.md`.

## Decision-leg review (rule 7 — orchestrator-written → Codex + Grok)
Prompt `directives/_legs/S11c_b_90_coupling_decision_review.md`; logs `~/.s11_build/S11c_b_90_coupling/decision_{codex,grok}.log`.
**v1 REJECTED by both legs, convergent + computation-backed** (orchestrator-verified, rule 13):
1. **Architectural (both):** `FACE_FLUX_BOUNDARY_OPERANDS` is the RAW T-substrate bundle, not operator rows.
   Weak-restricting it directly is the §3c-forbidden parallel route. Fix: compute the face generalized-force rows
   (coefficients of `δ_vu`/`δ_ve_W` in the consumed virtual work, as WL's `faceGeneralizedRows`), ADD to the
   reduced operator rows (origin FACE, not through the θ-fold), then weak-restrict the full operator.
2. **`closure_shape_deriv` skipped** — the boundary code selected only `traction`/`virtual_work_shape_deriv`/
   `evolution_mass_balance`, skipping the T-i `closure_shape_deriv` that carries Λ_A/Λ_V (verified: `:389` key,
   `:3535` selection; Λ inventory = traction/virtual_work→Λ_X, closure_shape_deriv→Λ_A/Λ_V, all else none).
3. **`A_T` token collision** — my "A_T = tilted-face geometry" clashes with PY's `A_T_s11cb` = the solenoidal test
   potential. Named the T-objects instead.
4. **Live `μ_θ` unbound**, double-count / over-reach (re-contract `evolution_mass_balance`, re-pair the already-weak
   virtual work, dump `virtual_constraint`/projection/`face_shift`, a `ζ_c` channel), §0 "every" overstatement, and
   rule-5 leaks in my acceptance/build-leg guidance.
v2 folded both legs (rule 7, one pass, no re-leg). Leak-gated clean.

## Build (Codex, `--sandbox danger-full-access`, 631k tok)
Log `~/.s11_build/S11c_b_90_coupling/build_codex.log`. Deliverable verified vs the DIFF (Codex-corrected numstat):
commit `7677aa18` = engine **+489/−116**, §0 +6, `DEFERRED_HEAVY_RUNS.md` +22 (the earlier "+517" conflated
engine+spec+deferred); parses OK; **`operator_from_density`/`committed_strong_rows` BYTE-IDENTICAL to HEAD**
(#88 refs preserved, both legs re-hashed them); the dead `bulk_kernel_from_density`/`paired_kernel_from_density` have
NO callers. ⚠ Codex ran its OWN "review legs" with two Codex variants (Claude/Grok "unavailable") — DISCARDED as
self-review (invalid per the authorship rule); the real legs are below.

## Build-leg review (rule 7 — Codex-written → fresh Claude agent + Grok). BOTH VERDICT CLEAR (8/8 probes)
Prompt `directives/_legs/S11c_b_90_coupling_build_review.md`; evidence
`~/.s11_build/S11c_b_90_coupling/{claude,grok}_buildleg_evidence/` (each leg's ablation scripts + literal stdout).
Both legs independently ablated a COPY; working tree untouched. The mandatory FORM ablations BITE on both legs:

- **Provenance (probe 1):** the face rows are computed INTO the operator (`U/E_EXPANDED = bulk + face`, verified
  `True` both legs; Grok's independent face-row reconstruction matches the emit, `RECON_*_MINUS_EMITTED = 0`); the
  kernel has NO `FACE_FLUX_BLOCKS`/`SOURCE_SUBSTRATE`, and `build_kernel` never reads `FACE_FLUX_BOUNDARY_OPERANDS`
  — NOT the parallel raw-bundle route.
- **FORM ablation (probe 2):** structural remaps of the virtual-work / closure(Λ) sources move DISJOINT kernel
  channels (virtual-work → E_W/DIV_U/reciprocal; Λ-closure → θ-forward, and `Lambda_A_0→Lambda_V_0` collapse changes
  the kernel's Λ); zeroing a NON-face key (projection) moves NOTHING. A typed/frozen face tag would be byte-identical.
- **Delete face feed (probe 3):** face-origin terms leave, bulk `γ·profile-jet` (count 26) stays.
- **Λ symbolic (probe 4):** `Λ_I⁰/(1−iωτ_I)`; no `Z`/DtN/impedance bulk-response solve in the build_kernel/operator/
  face slices (⚠ "no `.solve`" is scoped to the RESPONSE path — the UNRELATED constraint fold does use `sp.solve` at
  ~`:2459`/`:2475`; the point is no bulk-response/DtN solve, not that the engine is `sp.solve`-free). Grok independently validated
  the closure-fold recipe: the T-i identity `TRUE_AREA − RELFLUX_SUM = 0`, so `mass − closure` = replacing kinematic
  `J` by `Λ𝒜+ΛV` (`RECIPE_MINUS_TRUEAREA_REPLACE_ZERO True`).
- **Exact-once (probe 5):** `evolution_mass_balance` not re-contracted (`Σ origins − θ EXPANDED = 0`); no `δ_v` /
  test-field survives in the face rows (single weak pairing).
- **μ_θ bound (probe 6):** reserved `mu_theta_L/M` absent from kernel/blocks/operands; live operand bound
  (`MU_THETA_FACE_BINDING` live termcount 170).
- **No ζ_c / adjointness (probe 7):** 6 semantic operands, no center channel; adjointness `NO_INDEPENDENT_SECOND_ROUTE`
  over the complete (bulk+face) blocks.
- **Controls/#88/leak (probe 8):** #88 refs byte-identical; 0 asserts; §0 pin consistent with §1c/T-i; controls
  PRINT A/B/A−B (no assert); no expected-value leak in the diff.

## Two CROSS-ENGINE / step-record flags (NOT build defects — both legs CLEAR the construction)
A single engine cannot settle the CONTENT of the face response; only the PY↔WL residual can:
1. **Closure-fold sign/magnitude** (the Λ response into the θ-row via `mass_balance − Σ_faces closure_shape_deriv`,
   unit coefficient) — Grok's T-i identity check corroborates the RECIPE, but the sign/magnitude is a cross-engine
   check. "The one place a subtle physics error would hide" (Claude leg).
2. **Uniform-limit residual is nonzero and Λ-bearing:** `UNIFORM_RESIDUAL_TERMCOUNTS (2,4,0,0,0,0)` carrying
   `Lambda_A_0/Lambda_X_0`, γ-count 0 (Grok). §5c is a smoke test (no residual value supplied), so this is NOT a
   build gate — but whether the irreversible face-response coupling legitimately survives at `(η,σ_W)=0` or would
   violate §1d's uniform decoupling is a cross-engine + step-record adjudication (does WL's uniform coupling carry
   the same Λ survivors?). Also `Λ_V_0` appears in the θ OPERATOR row but not the weak transverse kernel for the
   LAB_HELD case — a per-case content question, not a presence quota.

## Status
Build CLEAR (both legs, 8/8, biting form ablations). No must-fix. The response-content questions (closure-fold
sign/magnitude, uniform-limit Λ survivors) are the DEFERRED cross-engine residual's job — and that residual is
doable HERE (single-case WL operator ~0.9 GB; not ≥64 GB-gated). Follow-on: the cross-engine integration pass
(comparator transliteration + `extract_slab` + `row_residual`, reviewed) → #88 KINETIC+θ re-adjudication + 2 owed
control-hardenings → step record + S11c card + close.
