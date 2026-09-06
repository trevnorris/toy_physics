# Decision review — the F/G diagnostic BUILD DIRECTIVE (before the builder runs). STATIC; ⛔ run nothing.

You are ONE of two decision legs reviewing an **orchestrator-written build directive** before the builder (astra)
runs. Check the requested F/G diagnostic **decisions + physics** are correct, complete, buildable, and leak no
answer. ⛔ **No fictional-script ablation** — the script does not exist yet; defer executable-control tests to the
build review. A real problem is a finding; ⛔ do not rubber-stamp.

## ⚠ Paths
Working dir = `/var/projects/toy_physics` (repo root). **Every** `directives/…`, `scripts/…`, `_measurements/…`,
and `scratchpad/…` path in this prompt is under `research/pde_ledger_v3/` — resolve each there. E.g. the directive
under review is `research/pde_ledger_v3/directives/S11c_c2_FG_diagnostic_build_directive.md`; the audit script is
`research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py`; the vet transcript is
`/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/codex_FG_question_vet.txt`.
⛔ Review THAT F/G diagnostic directive — ⛔ not any WL-engine / T7 / export-chain directive.

## Artifact
`research/pde_ledger_v3/directives/S11c_c2_FG_diagnostic_build_directive.md` (read it in full).

## Context (read for grounding)
This re-grounds F (uniform-limit decoupling) and G (directionality / adjointness) of the c2 self-energy increment,
which were resolved with **retired biased orchestrator scripts** (`verify_F.py`, `verify_EG.py`). A Codex-sol
question-vet sharpened both questions: `scratchpad/codex_FG_question_vet.txt` is the transcript — you may read it,
but ⛔ judge the directive on its own correctness, not merely that it matches the vet. Physics sources: c2 §3c
(increment), §3b (adjointness clause), §5a (`ORDERING_COMMUTATOR`), §5e (uniform limit = secondary smoke test),
`S11c_decisions.md` N4/N6; the audit machinery `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`.

## Checks (derive independently; quote lines)
1. **F predicate.** Is the **weak-kernel-identity** test — `I_closed|_uniform` vanishes as a weak operator after
   **re-specializing at the imported inputs**, with integral kernels inspected but Fourier integrals NOT evaluated
   — the right test, and does it actually exclude the retired `.doit()==0` proxy? Are `I_closed`, `I_bare`, and
   `ORDERING_COMMUTATOR` defined **distinctly and correctly** (the slot-carrier split `S_P = SLAB − SLAB|_{P=0}` is
   sound)? Is any co-specializing background jet / density binding missing from "the complete uniform
   specialization"? Is the live-sentinel + FORM ablation enough to prove the zero-checker bites? Any residual way
   the directive's F still admits a proxy zero (special field, parity cancellation, k=k′, dropping `∂_w δp_s`)?
2. **G object + adjoint audit.** Are the six directional blocks + the **generic-nonzero-witness** requirement +
   **per-face/assembled** reporting correctly specified? Is the **adjointness independence audit** (independent
   reverse-vs-adjoint routes, one-sided source corruption, ⛔ no `A−A`) the right way to settle the open
   "no-residual" sub-question rather than presuming it? Is the **representation scope** correct (Eulerian-only now;
   cross-representation directionality deferred to the not-yet-built N6 material route)? Does it correctly forbid
   the **dissipativity/passivity/energy-sign** overclaim?
3. **Buildability.** Can astra build this by importing/reusing the audit machinery (`bind_inputs`, `build_case`,
   `extract`, the §3a `close`, `retained_shape`, `shape_coefficients`) without a missing or nonexistent object? Is
   `S_P = SLAB − SLAB|_{P=0}` and the closed-pressure map `χ` reachable from those functions? Flag any object the
   directive names that the machinery cannot produce, or any under-specification that would let two builders
   compute different objects.
4. **No leaked acceptance criterion / no manufactured path.** Does the directive **withhold** the expected values
   (whether F weak-decouples, which G blocks are live, whether an independent adjoint route exists)? Any place it
   states the answer, supplies a target, or implies an iterate-to-zero exit? Are the objects named without
   dictating an algebraic path that pre-answers the result?
5. **Physics soundness + completeness.** Is the F/G physics correct? Any missing premise, wrong object, or control
   that would breed a wrong-but-consistent result (a frozen varying field, a check that can't fail, a residual
   structurally zero by construction)?

## Output
For each of 1–5: the requested-decision issues, CONFIRMED vs unsettled, with quoted lines. End: is the directive
**ready to build**, or the exact decisions/wording that must change first. ⛔ Do not ablate a nonexistent script.
