# Independent review — S11c-c1 T7 cross-engine comparator (SCRIPT; ablate it)

You are an independent, adversarial reviewer of a **Codex-written measurement instrument**: the S11c-c1 T7
cross-engine comparator. It reads the two committed c1 tag streams (the reviewed SymPy "PY" engine and the
reviewed blind Wolfram "WL" engine), joins them by object name on an axis-typed key, and prints `operand_A`
(PY), `operand_B` (WL), and the typed `A − B` residual per case, plus per-family accounting. **It must compute
and print, deciding nothing (rule 2).** Your job: find every way it manufactures a false AGREEMENT (a printed
zero where the engines disagree) or hides a real disagreement, and every way it violates the frozen contract.

⛔ **A prose judgement is worth nothing. ABLATE the instrument and report the LITERAL stdout diff.** Code-reading
alone has repeatedly missed real defects here. For any physics/algebra claim, write your own check as a script
and save both the script and its literal stdout to named absolute paths, or the claim is discarded.

Working dir: `/var/projects/toy_physics` (git repo, branch `ledger-v3-rebuild`).

## Artifacts under review
- `research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py`
- `research/pde_ledger_v3/scripts/test_S11c_c1_cross_engine_comparator.py`

## What you are handed (read these; ⛔ do NOT read the build directive — an artifact can satisfy its directive and still be wrong)
- ⭐ The FROZEN T7 contract (the authority): `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:580-587`
  and the "three-valued" source `research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md:42,53-60`.
- The reconciliation schema: `research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md:233-253`
  ("never blanket-collapse"; AppliedUndef→Symbol arg-strip is the named defect; μ fold is reviewed registry
  data, never comparator operand surgery).
- The sibling instrument + its tests: `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`,
  `test_S11c_b_cross_engine_comparator.py`; the transliteration `S11b_cross_engine_comparator.py:147`.
- The two committed inputs (⛔ `datalad get` both first; git-annex pointers):
  `research/pde_ledger_v3/scripts/out/S11c_c1_bulk_closure_sympy_audit.out` (PY, 63 tags),
  `research/pde_ledger_v3/mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out` (WL, 51 tags).
- The c1 PY engine construction `research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py` (the
  authority for what each PY key slot MEANS — its axis constants at `:119-124`, its `key=(...)` lines).

## OPERATIONAL (both legs — identical constraints)
- The `.out` payloads are LARGE (~82 MB WL / ~91 MB PY); several tags hold multi-megabyte srepr. ⛔ Do NOT
  full-CAS-parse them by hand; sample with `grep`/`awk`/`cut`.
- Running the comparator on the FULL streams may be memory-heavy (per S11c-b the full cross-engine residual can
  need ≥64 GB; this box is 30 GB). ⭐ Wrap every comparator run in `timeout 900` and bound memory by running a
  few families at a time with `--family <NAME>` (repeatable) for your ablations. If a FULL run OOMs, that is an
  operational deferral to REPORT, ⛔ not a defect of the instrument and ⛔ not a reason to skip the family-scoped
  ablations. ⛔ Do NOT run two heavy comparator processes at once.
- ⛔ Copy the script to `/tmp` and ablate the COPY; ⛔ never modify the working tree. Save every ablation script
  AND its literal stdout to named absolute paths and report those paths.
- First run `python3 -m pytest research/pde_ledger_v3/scripts/test_S11c_c1_cross_engine_comparator.py -q` (or
  `python3 .../test_....py`) and report the literal result.

## What to check — ABLATE, don't assert
1. **No verdict; three-valued preserved (contract `:583-587`).** Grep the script and a real run for any
   per-case `PASS`/`FAIL`/`VERDICT`/`AGREE`/`DISAGREE`/`STATUS`/`FINAL_STATUS` token — there must be none. But
   the inherited three-valued residual OBJECTS (`BooleanNotResidualable`, `UndecidedResidual`/`ResidualFailure`,
   `Mismatch`) MUST still be printed as `A_minus_B` where they arise — confirm the instrument did not strip them
   to force a scalar. A native boolean must be REJECTED as a residual operand, not silently subtracted.
2. **False-agreement ablations (the core).** For a family that JOINS with a small/zero printed residual, prove
   the residual is LOAD-BEARING:
   - **One-sided corruption:** perturb ONE engine's operand for one case (e.g. flip a sign / add a term to the
     PY `.out` copy leaf) and re-run that family. The residual MUST move. If breaking PY also leaves the residual
     zero, the comparator is not actually subtracting the operands — a defect.
   - **Repoint ablation:** substitute a DIFFERENT object's payload under a previously-paired NAME (model on
     `S10_cross_engine_comparator_repoint_ablation.py`); the residual MUST move. ⛔ A symbol rename is not a
     repoint.
   - **Form ablation of the axis typing:** make FACE and DIRECTION merge (or make `Integer(-1)` type as
     DIRECTION) and re-run; confirm the change is VISIBLE (joins/accounting move). If merging FACE↔DIRECTION
     silently increases joins with zero residuals, that is the manufacture-agreement path.
3. **The SEALS must remain unreconciled (rule 6 / contract `:585`).** Confirm the comparator does NOT collapse
   any of these — each must SURFACE as a residual/accounting row, not a printed zero:
   - the two-momentum legs: PY opaque `s11cc1_q_out_output/_input` vs WL live `qOut[omega,{kOne..}]` /
     `qOut[omega,{kPrime..}]` — ⛔ `qOut[...]` must NOT be arg-stripped to a bare `q_out`; `kOne`/`kPrimeOne`
     must NOT be mapped onto the PY `k_output`/`k_input` legs;
   - μ_θ: PY face-specific composites (`s11cc1_mu_theta_lab_held_plus` …) vs WL opaque `muTheta` — no registry
     fold;
   - the ω real-assumption: `Symbol('omega')` vs `Symbol('omega', real=True)` not normalized;
   - the background density: PY bare `Symbol('rho_br_bg_rho4_constant')` vs WL applied field
     `rhoBrBgRho4Constant[x…]` — ⛔ not arg-stripped (rule 17);
   - regime/parity: WL `OUTPUT_PROPAGATING`/`PARITY_DELTA_W` vs PY `PROPAGATING`/`THICKNESS` not declared equal.
   ⭐ Verify each by a name-map ablation: re-point ONE seal's name in the comparator's map (on your /tmp copy)
   and confirm the previously-surfaced residual MOVES — i.e. the seal is real, the map is load-bearing.
4. **`raw_control_case` scope.** The outer-signature (`ASSOCIATION_KEYS`/`TUPLE_ARITY`) shortcut is legitimate
   ONLY for `DTN_OPERATOR` and `NONINVERTIBILITY_CONDITION`'s `OPERATOR` leaf (different operator algebras). For
   every other family — especially `CONTROL_FORM_*`, `ENERGY_*`, `FACE_RESPONSE*` — the comparator MUST descend
   to the paired scalar leaf and residual it. Find any family where a reachable scalar leaf was replaced by an
   outer-signature "agreement". Show the leaf you can reach and the signature the comparator emitted instead.
5. **No silent 0-extract; coverage is real.** Every one of the 50 shared non-`LOCAL_` families must print an
   `ACCOUNTING` line and reach real leaves — a family that "joins" one generic/`raw_control_case` case while its
   real content goes unmeasured is under-measured. Check the per-family extracted-leaf count; name any family
   that silently extracts 0 or measures only an outer shell. Confirm the `_LOCAL_` inventory is emitted and the
   13 PY / 1 WL locals are excluded from the join (not compared).
6. **Blanket-collapse / arg-strip audit.** Grep the name/CAS map for any `AppliedUndef→Symbol` that strips
   arguments, any blanket `X(args)→X`, any spelling entry that is NOT `mechanical_lower_camel(<py>)` of a real
   PY `Symbol('<py>')`. Report each with the line and why it could manufacture agreement.
7. **Injectivity + `Inactive[Greater]`.** Confirm the spelling map is injective (no two reserved names collapse)
   and that `Inactive[Greater]`/`Inactive[FourierTransform]` are kept UNEVALUATED (not turned into a native
   boolean or a false 0). WL `DTN_BY_REGIME_PAIR` carries `Inactive[Greater][…,0]`.

## Physics filter
Report a finding only if it catches a way the MEASUREMENT could be wrong — a false agreement, a hidden
disagreement, a dropped family, a broken seal, a non-load-bearing residual. ⛔ Do not report "it would be wrong
on a different input" for inputs these two committed streams do not contain.

## Output
For each of the 7 items: CLEAN / DEFECT, with the ablation command, its LITERAL stdout, and the absolute path of
your ablation script + output. Rank real defects most-severe first. A clean pass with citations is equally
useful. ⛔ Do NOT propose an expected residual value or a physics reconciliation to bake in — a physics-bearing
difference stays a SURFACED residual for post-run adjudication.
