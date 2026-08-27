# Independent physics review — S11c-a cross-engine comparator (Codex-written script)

## Artifact
`research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py` (+ any synthetic self-test file it wrote).
It reads two committed CAS engine transcripts and prints, per tag family, `operand_A` (PY), `operand_B` (WL),
and the residual `A − B`, plus per-family accounting. It is a MEASUREMENT instrument: it must compute and
print, decide nothing, and — above all — never fold a genuine cross-engine DISAGREEMENT into false agreement.

## What to check (the physics the instrument bears)
This comparator will judge whether two independent engines agree on 39 tag families of the S11c-a interface
shape derivative. A defect here corrupts that judgment. The three highest-stakes failure modes:
1. **A smuggling fold that manufactures false cross-engine agreement.** Especially the perturbation-current
   fold: it MUST be a strict arity-preserving AST-head/symbol-base rename (`currentWPerturbation(args)→
   delta_j_bulk_4(args)`, etc.), and MUST NOT (a) map an `AppliedUndef` to a bare `Symbol` (arg-strip — this
   is exactly what hid finding #3), nor (b) reintroduce any "held-context"/`J(w; held=x,t)` reduction that
   drops or annotates away WL's `x,t` arguments. PY `delta_j_bulk_4(w)` and WL `delta_j_bulk_4(x1,x2,x3,{w,t})`
   must remain DISTINCT after the fold. Also check MUMAP is a per-branch registry (never a global
   `mu_theta_L/M→mu_theta` collapse) and CONORMAL is compared raw (no §3c Taylor fold).
2. **An integral canonicalization that manufactures agreement OR a false disagreement.** It must retain the
   binder and `(lo,hi)`; combine two integrals ONLY over identical limits; pull a factor out of an integral
   ONLY when that factor is free of the integration variable. A canon that pulls a w-dependent factor out, or
   that combines mismatched-limit integrals, is a defect.
3. **Silent under-coverage.** Every one of the 39 families must print an accounting line and either `join>0`
   or a documented `axis_set_mismatch`/`py_only`/`wl_only`. A family that silently extracts 0 cases (e.g. the
   object-nested controls under a `VALUE`-only extractor) is a defect — it reads as "nothing to disagree
   about" when in truth nothing was compared.

## What you are handed
The artifact above; both committed transcripts (PY
`research/pde_ledger_v3/scripts/out/S11c_a_interface_geometry_sympy_audit.out`, WL
`research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`); the spec
`research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (§1a/§1b/§3a/§3b/§3c/§5–8). Derive your own view
of what the comparator should do from these; do NOT assume its self-report is correct.

## Required method (SCRIPT)
- ⚠ The full comparator run is ~37 minutes. You do NOT need to complete a full run. Your PRIMARY method is
  **fold ablation on small synthetic fixtures** (fast) plus **reading the extractor/keying/accounting logic**.
  The builder's reported per-family accounting (for coverage context, NOT physics results) is:
  simple families join cleanly; FACE_SHIFT `join0/py160/wl80/axis_mismatch240` (WL missing DENSITY, not
  deduped); ADMISSIBILITY `join0` (PY missing BRANCH, not broadcast); CONTROL_FORM PY528/WL960 granularity
  surfaced; object-nested controls join (not 0-extracted). **Verify by reading the code that the object-nested
  control families (REP_INVARIANCE/CONTROL_INDEPENDENCE/UNIFORM_LIMIT/HOMOGENEITY) actually flatten and
  extract — that they do NOT silently extract 0 under a `VALUE`-only path.** If you want to confirm coverage
  empirically, run the comparator but you may kill it once the per-family accounting has printed; do not wait
  for all 6,714 residual triplets.
- **Ablate every load-bearing fold and report the literal diff:**
  - Current fold: construct a synthetic PY bare `Symbol('delta_j_bulk_1')` vs WL `currentXPerturbation1[x1,w]`
    field; confirm the comparator's fold leaves them a NONZERO residual (arity preserved). Then try to make the
    fold strip args (the failure mode) and confirm the real code does NOT do that. Grep the code for any
    `Symbol(` construction on a current head, and for any `held`/spectator reduction.
  - Integral canon: feed two integrals with different `(lo,hi)` and a shared integrand; confirm they are NOT
    combined. Feed `c(w)*f(w)` inside an integral and confirm `c(w)` is NOT pulled out.
  - Keying: perturb one axis token and confirm the join count changes (the key uses that axis).
- ⛔⛔ **A FORM ablation is mandatory:** change the STRUCTURE of a load-bearing object (flip the current fold to
  arg-stripping; drop bound-equality from the integral canon; collapse two axes in the key) and report the
  literal change in the per-family residual/accounting. A COEFFICIENT change tests arithmetic; only a FORM
  change tests whether the check has teeth.
- Report any `assert` that runs on the measured payloads (it must not — assertions belong in a separate
  synthetic test file), and any place the comparator emits a conclusion/verdict rather than operands+residual.

## Physics filter
Report a finding only if it catches a way the comparator could smuggle a false cross-engine agreement, hide a
real disagreement, or silently fail to compare a family. Do not report "would be wrong on a different input."

## Ablation sandbox
Copy the artifact to `/tmp` and ablate the COPY; ⛔ never modify the working tree. Save every ablation script
AND its literal stdout to named absolute paths and report those paths. A prose re-derivation without a script
and its stdout is discarded.
