# Independent review — S11c-c1 step record + its cross-engine reconcile verdict

You are an independent, adversarial reviewer. The artifact is an **orchestrator-written step record** that claims
a physics conclusion — that the two blind engines (SymPy + Wolfram) **AGREE** on the S11c-c1 curved-interface bulk
closure. Your job: decide whether that verdict is EARNED, whether the reconcile's central computational claims are
real (⛔ not tuned to zero), whether the rule-17 caveat is correct, and whether the record faithfully reports its
sources without overstatement. ⛔ A prose judgement is worth nothing — where a claim rests on computation, RE-RUN
it (or write your own) and report the LITERAL stdout with the command. Save every script + output to named
absolute paths and report them.

Working dir `/var/projects/toy_physics` (branch `ledger-v3-rebuild`). All paths below are under
`research/pde_ledger_v3/`.

## The artifact
`steps/S11c_c1_curved_bulk_closure.md`

## Sources of truth (read these and form your OWN view BEFORE trusting the record)
- The reconcile record + its committed bridge scripts: `directives/_measurements/S11c_c1_comparator_reconcile.md`
  and `directives/_measurements/S11c_c1_reconcile_*.py` / `*.txt`.
- The comparator (SOUND, do NOT re-review it): `scripts/S11c_c1_cross_engine_comparator.py`.
- The two engine transcripts (git-annex — `datalad get` first; ⚠ large ~82/91 MB, sample with grep/awk, NEVER
  full-CAS-parse): `scripts/out/S11c_c1_bulk_closure_sympy_audit.out`,
  `mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out`.
- The two engine sources: `scripts/S11c_c1_bulk_closure_sympy_audit.py`,
  `mathematica/S11c_c1_bulk_closure_mathematica_audit.wl`.
- The frozen contract: `directives/S11c_c1_SHARED_PHYSICS.md:580-587` (the comparator computes and PRINTS, deciding
  nothing; representational differences are adjudicated AFTER the run, ⛔ never pre-folded).

## Checks (CLEAR / DEFECT, each with the command + literal stdout)

1. **The central claim is REAL, not tuned to zero (the most important check).** Re-run the committed bridge scripts
   `S11c_c1_reconcile_dtn_kernel_bridge.py` and `S11c_c1_reconcile_response_coeff_bridge.py` and confirm they print
   the collapse (kernel residual → 0 on-shell; `ε·A≡B` worst_rel ~0). Then ADVERSARIALLY test that the collapse is
   MEANINGFUL by one-sided corruption of an identification, e.g.: (a) use a WRONG jet relation
   `jet_hat_i → i·L_W·(k_i + k′_i)·hat_transfer` (a sign flip) — the residual MUST become nonzero; (b) DROP the
   on-shell substitution — the Stage-2 remainder MUST be the dispersion polynomial `Q_IN²c_s0²+|k′|²c_s0²−ω²`, NOT
   an arbitrary expression; (c) map WL `qOut(ω,{k})` and `qOut(ω,{k′})` to the SAME opaque symbol (freeze the two
   legs into one) — the kernel residual MUST become nonzero (this proves seal 3's two-momentum structure is really
   carried by both engines, not manufactured by the map). If any corruption FAILS to move the residual, the
   collapse was vacuous → DEFECT. Report each literal result.

2. **The identifications are physically justified, read from the engine sources — not invented to force a match.**
   For each map entry in the reconcile record §2, confirm it against the two engine sources: the FT-of-derivative
   jet identity (`.wl:174-179` slopeHat built from profileHat; `.py:200-207,606` jet opaque), the
   `σ_W=η_bg·W_0/L_W` binding (`.wl` SIGMA_BINDING), the ε-placement (`.py:840` `mu_theta=ε·…` vs `.wl:748`
   `ε·muTheta`), on-shell dispersion (WL encodes it in `qOut`, PY writes `κ²=ω²/c_s0²`). Any identification with NO
   source basis, or that hides a real content difference, is a DEFECT.

3. **Rule 17 (seal 5, background density).** Independently verify the record's claim that NEITHER engine
   differentiates the density field: grep both engines for any `D[...rhoBrBgRho4Constant...]` /
   `Derivative`/`diff(...rho_br_bg_rho4...)`. Confirm WL binds `rhoBrBgRho4Constant(x)=(ρbr/W_0)·W_bg(x)` (a LIVE
   varying field) yet uses it only pointwise (`μ_s=μ_θ/ρ`). Is the AGREE verdict + the c2 re-adjudication caveat
   correct, or is a frozen varying field being waved through (a rule-17 violation)? This is the one seal where a
   wrong "AGREE" would be most costly — probe it hardest.

3b. **The comparator really SURFACES the seals (no pre-registered fold).** Confirm from the run that the seal
   families are UNRECONCILED/nonzero as printed (DTN_KERNEL joins nonzero; regime family unjoined by
   `OUTPUT_/INPUT_` naming; density appears bare-PY vs applied-WL) — i.e. the record's "AGREE" is OUR post-run
   adjudication, not something the comparator baked in (contract N8, rule 5/6).

4. **Deferrals are honest (⛔ not narrowing to fit).** Check that the 4 giant families + the full symbolic residual
   are deferred by SIZE/RESOURCE (30 GB box), not to dodge a hard comparison, and that the two UNDECIDED coverage
   gaps (raw DTN_OPERATOR; ENERGY closed-form-vs-integral) are honestly flagged as OPEN, not glossed as AGREE.
   Confirm `DEFERRED_HEAVY_RUNS.md` records the re-run command.

5. **No overstatement / faithful to sources.** Does the record claim more than the per-engine records + the
   reconcile support? In particular: is "per-engine SOUND (2-leg each)" accurate (incl. the WL repair-directive
   rule-7 gap + its remediation), and is "cross-engine AGREE" correctly scoped to what actually ran (46/50 families,
   first order for the Hermitian parts, giants+full-residual deferred)? Flag any claim the sources don't back.

## Method + output
⛔ Write your own verification scripts; save them + literal stdout to named `/tmp` absolute paths and report those
paths. ⛔ Copy any file you ablate to `/tmp` and ablate the COPY; never modify the working tree. Per check:
CLEAR / DEFECT with the command, literal stdout, and your script/output paths. Rank real defects most-severe first;
a clean pass with citations is equally useful. ⛔ Do NOT propose an expected measured residual for a deferred family
— a physics-bearing difference stays a SURFACED residual for the ≥64 GB adjudication.
