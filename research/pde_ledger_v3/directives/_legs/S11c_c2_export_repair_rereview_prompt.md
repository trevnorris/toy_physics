# Independent re-review — the S11c-c2 EXPORT REPAIR (publication-only diff)

You are ONE of two independent legs reviewing a **publication-only** change to an already-physics-reviewed SymPy
script. You will not see the other leg's output. The physics of this script was independently reviewed clear
(the fold is sound); this change touched **only** the export serialization to stop over-exporting. Your job: catch
any way this change (a) altered the physics, (b) exports the wrong set or loses information, or (c) installed a
guard that does not actually bite. A build is reviewed **until clear**; a real defect is a finding, not a failure.

## Artifact under review
The diff of `research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py` from commit `8f3a017f`
(the reviewed baseline) to the current working tree — plus the regenerated `research/pde_ledger_v3/scripts/S11c_c2_exports.py`.
See the diff with: `git -C /var/projects/toy_physics diff 8f3a017f -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py`.

**What the change was supposed to do** (the requirements): `research/pde_ledger_v3/directives/S11c_c2_export_repair_directive.md`.
Read it — but your job is to verify the RESULT is correct, ⛔ not merely that it matches the directive's wording.

## What you are handed (read any; nothing you may use is withheld)
- the diff + the new `scripts/S11c_c2_exports.py` (60 MB → 21.4 MB claimed);
- the requirements directive (above) and the physics adjudication for context
  `research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md` (the physics is sound; the
  self-energy increment is a comparator/EMIT-only object, not a downstream binding);
- the fold machinery `scripts/ledger_fold.py` (`load_model`@102, `check_consumer`, `assert_delta_is_minimal`) and
  the parent exports `scripts/S11c_b_exports.py`, `scripts/S11c_c1_exports.py`;
- the emitted output of the build's final run: `/tmp/S11c_c2_export_repair.out` (~499 MB; the physics tags are the
  `PY_S11CC2_*` lines — the operator VALUEs are emitted there in full/expanded srepr form).

## The checks (derive/verify independently — ⛔ do not trust the script's own printed EXPORT_SEMANTIC results)
1. **Physics unchanged.** Confirm the diff touches ONLY the export path (`EXPORT_ROOTS`, the `export_key` map in
   `run()`, `publish`, and a new publication helper) and **no construction** (`build_case`/`build_face`/`extract`/
   `kernel_bridge`/`retained_shape`/`emit`/`traction_pairing`/`control`/`grade_object` or the physics loop). Then
   confirm the emitted physics is unchanged: re-run one case in TRIAGE mode (`S11CC2_PACKAGE=TRIAGE`, one case, the
   3 core objects — fast) on BOTH `git show 8f3a017f:…/S11c_c2_selfenergy_fold_sympy_audit.py` and the current
   script, and confirm the emitted `CLOSED_SLAB_OPERATOR` / `CLOSED_COUPLING_KERNEL` VALUEs are identical.
2. **Membership.** In the new `exports.py`: `s11cc2SelfEnergyIncrement` is **absent**; `s11cc2ClosedSlabOperator`
   and `s11cc2ClosedCouplingKernel` are present, each a 4-case tree; the delta is **minimal** (exported keys ==
   required closure, no extra rows); `check_consumer` closure resolves and a `load_model` over the two parents +
   this delta binds them. ⛔ No term_origins / parity / §3d / control rows leaked into the delta.
3. **Semantic equivalence — RECOMPUTE it yourself.** For at least one case and both operators, independently
   compare the **compact exported VALUE** (from `exports.py`, via its `_restore`) against the **emitted VALUE**
   (from your TRIAGE re-run, or extracted from `/tmp/S11c_c2_export_repair.out`): leafwise `sp.expand(compact −
   emitted) == 0`, **Integral-aware** (top-level `expand` leaves `Integral(0,…)` unevaluated — `doit` the integrand
   or protect it). ⚠ Verify the container shapes match exactly (case keys, matrix shapes, mapping keys) — a
   dropped trailing channel could compare as zero under a naive `zip`.
4. **Transparency + no information loss.** The stored value must be a real evaluable `sp.Basic` a downstream step
   can import and **differentiate** — ⛔ no `UnevaluatedExpr`, no surviving `cse`/`Dummy(` temporaries, no
   string-only wrapper. Grep the new `exports.py` for `Dummy(` / `UnevaluatedExpr`. And confirm the compaction did
   not change the **singular locus** (the denominator/pole set of each leaf is unchanged — a `cancel` that removed
   a pole would be lost physics for the downstream resonance step even if the value "equals" elsewhere).
5. **⛔⛔ FORM ABLATION (MANDATORY) — does the guard actually bite?** Copy the script to /tmp and:
   (a) corrupt the **compaction** so it silently drops or alters a real term of one operator VALUE, re-run, and
   confirm the added semantic check **hard-stops** (nonzero difference). If a corrupted compact form passes, the
   check is decorative — a finding.
   (b) corrupt the **membership** (re-add the increment to `EXPORT_ROOTS`, or drop one operator), and confirm the
   minimality/closure guard catches the over/under-export.
   Report the literal outcomes.
6. **No overstep.** Confirm the working tree carries no builder-created review sprawl (self-review /
   derived-or-declared / output-checker / finalize artifacts) from this build — only the script, the export, the
   builder report, and the expected `_measurements/S11c_c2_sympy_guard_evidence.json` checkpoint.

## Physics / correctness filter
Report a finding only if it catches a real defect: physics changed, information lost in compaction, wrong export
membership, a broken/decorative guard, or an un-importable/un-differentiable stored form. ⛔ Not style nits.

## Sandbox + limits
Copy anything you run to `/tmp` and work there (`ROOT` resolves from the script's own location, so a copy under
`/tmp/<dir>/scripts/` writes there). ⛔ Never modify the working tree. TRIAGE mode runs one case fast; a full run
is ~19 min / ~1.5 GB — budget it, one heavy job at a time, watch `free -h`. Save every script + its literal stdout
to named absolute paths and report them; a prose-only derivation is discarded.

## Output
For each of 1–6: finding + evidence (your script path + literal stdout, the diff `file:line`). Separate CONFIRMED
defects from unsettled questions. End: is the export repair correct to commit, or what must change?
