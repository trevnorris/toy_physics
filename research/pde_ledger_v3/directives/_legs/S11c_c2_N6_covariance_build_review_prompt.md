# Independent BUILD review — S11c-c2 N6 covariance instrument (a SCRIPT)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py` (astra-written, 277 lines).
Working dir `/var/projects/toy_physics`. SCRIPT review — derive independently, then ABLATE.

## Role + the ONE boundary
The instrument computes the source-naturality residual `R_cov = ms − ms_pred` that closes the OPEN per-engine N6
question: is the constitutive residual the covariant transformation of the source under the declared field map `Φ`, or
a real non-covariance / `material_pullback` defect? ⛔⛔ **Verify the instrument COMPUTES `R_cov` CORRECTLY — do NOT
adjudicate whether `R_cov` vanishes.** There is NO expected value for `R_cov`; ⛔ never report "the residual came out
nonzero/zero" as a finding. The strict-vs-covariance interpretation is the orchestrator's.

## The governing construction (implement THIS — derive from it, ⛔ not from the output)
Directive `directives/S11c_c2_N6_covariance_directive.md` (CLEAR TO BUILD); framing
`_measurements/S11c_c2_N6_sufficient_test_vet_adjudication.md`. Settled points to verify the script implements:
- **`ms` = ACTUAL material source** (pull-back-then-vary): `source_terms(..., μ_M, V_M)`, `μ_M` via
  `b.material_pullback` then vary (the diagnostic `ms`).
- **`ms_pred` = PREDICTION** (vary-then-substitute-Φ): `source_terms(..., μ_E.subs(Φ), V_E)`, where `μ_E` is the
  IMPORTED Eulerian amplitude `es` uses (`inputs.mu/ε`) and `Φ` is the declared field map.
- ⛔⛔ **NON-CIRCULARITY (the crux).** `ms_pred` must be built from the IMPORTED `μ_E` + the SUPPLIED `Φ` — ⛔ NOT by
  calling `b.material_pullback` (the builder under test), ⛔ NOT by subtracting material from Eulerian, ⛔ NOT reusing
  `ms`/`ΔS`. If `ms_pred` re-imports `material_pullback` on the prediction side, the test is vacuous — REPORT it.
- **F1 — the Φ prolongation must be COMPLETE.** `μ_E = EL(E)` carries **second** jets (`theta_didi`, `e_W_didi`;
  confirmed at `S11c_b_exports.py:5629`). `μ_E.subs(Φ)` must prolong the field map through EVERY θ/`e_W` jet in `μ_E`
  (rank-2: `θ_I↦D_I(θ+a_ρ)`, `e_{W,I}↦D_I(e_W+h_α)`), generated via the LIVE
  `b.total_derivative(...,background_depth=3)`+`DERIVATIVE_MAP` chain, `simultaneous=True`, with `PHI_DOMAIN_CENSUS`
  showing **zero uncovered** atoms. ⛔ Reusing the energy-level 0+1 dict (fields + first jets) would leave second jets
  unshifted (omitting `D_iD_j a_ρ`, `O(σ_W)`, live at retained order) and mislabel a covariant construction as a
  defect.
- **`R_COV_INCREMENT`** = reconcile `closed_response(m_coeff, R_cov)` (signatures 6/9/12 only), ⛔ NOT
  `build_increment`/`I(C_M,R_cov)` (which re-adds `−C_M·p` → nonzero at `R_cov=0`).
- **PIT:** one in-process joint PIT per case; nonzero = one-sided certificate; all-zero = "no nonzero found" at
  conditional δ; residual-zero ⛔ never an exit.

## Required method — SCRIPT branch (derive, then ablate)
1. **Derive independently.** Write your OWN probe: confirm on the real objects that `ms_pred` built from
   prolonged-`Φ`·`μ_E` differs from a naive 0+1-jet substitution EXACTLY by the second-jet prolongation (so F1 is
   load-bearing), and that `ms` = `material_pullback`-then-vary. Save script + literal stdout to named /tmp paths;
   report them. ⛔ Without them your claims are discarded.
2. ⛔⛔ **FORM ABLATION MANDATORY** (on /tmp COPIES; ⛔ never the working tree). E.g.: (a) truncate the Φ prolongation to
   0+1 jets (drop the rank-2 map) → `R_cov` must MOVE (proves the F1 prolongation is load-bearing, not decorative);
   (b) make `ms_pred` call `material_pullback` (the circular version) → `R_cov` collapses to ≈0 for ALL inputs (proves
   non-circularity is doing work); (c) make `R_COV_INCREMENT` use `build_increment` → a spurious `−C_M·p` appears at
   `R_cov=0`. Report the LITERAL diffs.
3. **Verify the knives bite (the instrument ships uncorrupted: `κ_a=1`, `κ_j=0`).** Run with the Φ-coefficient knife
   (`κ_a=2` on the ACTUAL path, prediction fixed at `κ=1`) → `R_cov`/`R_COV_CONTROL_DELTA` must move even though the
   `a_ρ+h_α` truth table is unchanged. Run the junk knife (`κ_j=1`, `J_μ·e_W`) → `R_cov` nonzero at
   `MATERIAL_ADVECTED.RHO4` (where `R_N6=0`). If a knife does not bite, the control is vacuous — report it.
4. **Classic defects + WHICH LINE COMPUTED THIS** for every emitted object; report any `assert` before an emit;
   confirm `PHI_DOMAIN_CENSUS` has zero uncovered and is a real computed coverage check (not a hardcoded pass).
5. **PIT soundness:** shared samples, primes, joint singular rejection, on-shell `q` from `k`, honest FN bound,
   residual-zero never an exit.

## Ablation sandbox / ops
⛔ Copy to /tmp, ablate the COPY. Pure SymPy (no Mathematica). Carrier-class so light, but wrap each run in
`timeout 900`, ONE case at a time. astra's baseline `.out` are at `/tmp/S11c_c2_N6_covariance_sympy.<...>.out`
(a_knife/junk_knife are astra's own knife runs — you generate your own from copies). Save every ablation + stdout to
named /tmp paths and report them.

## Physics filter
Report a finding only if it catches a way the covariance INSTRUMENT could be circular, wrong, vacuous, intractable, or
answer-leaking — ⛔ not "wrong on a different input", ⛔ not the `R_cov` disposition.

## Output
Findings each with the line, the ablation + literal diff, why it matters, minimal fix. If nothing outstanding changes
what the instrument computes or may be claimed, say **BUILD CLEAR**. Evidence-first, brief.
