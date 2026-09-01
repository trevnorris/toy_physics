# Build review — S11c-b #89 PY §3a repair engine (+ lever-C projection)

You are one of two independent legs reviewing a **built engine**. Derive independently; ablate load-bearing
objects and report the LITERAL diff (code-reading alone has repeatedly missed real defects here). Report a
numbered list (file:line — problem — how you found it), then a one-line verdict.

## Artifact
`research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (working tree = committed) and its output
`research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out`.
⚠ The `.out` is a **git-annex pointer** (symlink); its content is present locally — read/grep it directly (if it ever
resolves empty, run `datalad get <path>`). Treat the `.out` as read-only INPUT: do NOT modify it. The ablation
target is the ENGINE `.py` — copy IT to /tmp and ablate the copy (never touch the working tree). This is a
SymPy/Python engine (no Mathematica kernels), so run python freely — no kernel-seat serialization needed.

## What it is (context, all readable)
Two changes to review, distinct in purpose:
- **A · the §3a physics repair** — the engine was under-counting the variable-coefficient energy basis by
  FREEZING the background profile's spatial-jet tower (treating `∂W_bg`, `∂μ_R,bg` as constant so the Hessian
  `∂²W_bg` a divergence must generate was dropped). The repair un-freezes it across four consumers: the §3a
  basis quotient (`basis_euler_signatures`/`quotient_independent_indices`), the §3b strong rows
  (`operator_from_density`), the MATERIAL pullback (`material_pullback`), and the coupling differentiation
  (`operator_dx`/`operator_divergence`/`operator_curl`). The corrected basis is the object under review.
- **B · lever-C projection** — a pure performance change to the retained-grade projector. The original
  `first_shape_series_reference` uses `sp.series(...,eta_bg,0,2)`; the new `first_shape_series_fast` computes
  the first-order Taylor `f(0)+eta_bg*f'(0)` (identical to `series(...,0,2).removeO()` for these analytic-at-0
  expressions), falling back to the reference for `Integral`/`Derivative`-bearing scalars. It is meant to
  change SPEED, not physics — but it emits results in a **reduced rational form** (e.g. `X/D` where the
  reference leaves `2X/(2D)`); confirm that difference is representational, never a value change.

Supplied references (you may use them — you are the reviewer): the corrected basis is **40** = 10 uniform +
15 `∂W_bg` + 15 `∂μ_R,bg` (per-source 15; `_measurements/S11c_b_86_reference_result.md`), reducing to the
committed frozen **26** when the Hessian is zeroed. The spec is `directives/S11c_b_SHARED_PHYSICS.md`
(§1c/§1d/§3a/§3b) — §3a lines 251–254 require retaining a second spatial derivative of `W_bg` at background
order.

## Required method
⛔ Wrap every kernel/engine run appropriately and copy the artifact to `/tmp` before ablating — never modify
the working tree. Save every ablation script AND its literal stdout to named absolute paths and report them.
A prose derivation with no script+stdout is discarded.

## Checks — part A (the physics; this is the load-bearing review)

1. **Is the corrected §3a basis genuinely complete, and does it reduce to 26 frozen?** Independently
   reconstruct (your own script) the variable-coefficient new-invariant count per source with the background
   first-jet spurion carried WITH its own first jet (the Hessian) in the divergence/quotient, and confirm it
   is 15/source (nullity 0), total basis 40. Then FORM-ablate: zero the Hessian map (the second-jet entries)
   in a `/tmp` copy and re-run — the basis MUST drop to the committed 26 and the strong rows to their frozen
   forms. Report the literal counts both ways. A basis that does not MOVE under the Hessian ablation means the
   repair did nothing.
2. **Spurion enters the divergence map, not the variation-field set.** Confirm the quotient differentiates the
   spurion via a derivative/Hessian table but does NOT take an Euler–Lagrange variation `δ/δ(spurion)` (that
   would wrongly treat fixed background data as a dynamical field). Check `basis_euler_signatures` and its
   `fields` argument.
3. **All four frozen paths are actually repaired.** For each of the strong rows (`operator_from_density`),
   the MATERIAL pullback (`material_pullback`), and the coupling cascade (`operator_dx` depth), construct a
   test expression carrying a background first jet, apply the engine's divergence, and confirm the Hessian
   (and, for the coupling, the third jet) term is generated — not dropped. Report the literal generated term.
4. **Retained grade / bookkeeping.** Confirm the generated higher background jets are σ_W¹ (retained), and
   that `first_shape_series` keeps them (η≤1 ∧ σ_W≤1) rather than truncating them away.

## Checks — part B (lever C; physics must be unchanged)

5. **Physics-identity of the fast projector.** On a broad set of real operator AND kernel scalars, compare
   `first_shape_series_reference(x)` vs `first_shape_series_fast(x)` NUMERICALLY (substitute random values;
   symbolic `simplify` crashes on these) and confirm every difference evaluates to ~0. Any genuinely nonzero
   difference is a physics bug — report it with the scalar.
6. **Integral fallback + determinism.** Confirm `Integral`/`Derivative`-bearing scalars route to the exact
   reference (byte-identical there). Confirm the parallel path (`retained_grade_parallel`) reassembles in
   order (serial vs parallel object residual = 0).
7. **The reduced-form change is acceptable downstream, or it is not.** Judge whether emitting reduced rational
   coefficients (vs the reference's un-reduced form) could shift the PY-vs-WL comparison. State your reasoning
   (the comparator canonicalizes rationals — does that cover this?).

## Checks — part C (soundness)

8. No hard-coded counts/operators; no tautological residuals; the emitted `PROJECTION_EQUIVALENCE_RESIDUAL`
   rows are representational-only (numerically zero) not hidden value changes; the frozen-limit switch yields
   26. Report anything that would let a wrong result pass.

## Physics filter
Report a finding only if it catches a way the physics is wrong (part A), the projection changed a value
(part B), or a wrong result could pass unnoticed (part C). Not style.

## Output
Numbered discrepancies (file:line — problem — how found, with script+stdout paths), then `VERDICT: CLEAR` or
`VERDICT: N defects`. Nothing else.
