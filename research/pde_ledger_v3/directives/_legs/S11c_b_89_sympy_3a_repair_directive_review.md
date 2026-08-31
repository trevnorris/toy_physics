# Decision-list review — S11c-b #89 PY §3a repair directive (before any builder)

You are one of two independent legs reviewing an **orchestrator-written build directive** before it is
handed to a builder. Your job is to find defects in the DIRECTIVE — ambiguities, wrong premises, missed
code paths, answer leaks, and non-decisive controls — not to build anything. The directive governs an
in-place repair of a SymPy engine; a defect here costs a full build round plus its two review legs, so this
gate exists to catch it now. Report a numbered list (file:line — problem — correction), then a one-line
verdict.

## Artifact under review
- Directive: `research/pde_ledger_v3/directives/S11c_b_89_sympy_3a_repair_directive.md`
- Its premises: `research/pde_ledger_v3/directives/_measurements/S11c_b_89_sympy_3a_repair_directive.md`

## What you are handed (all readable — you are a reviewer, not the blind builder)
- The spec: `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (esp. §1c/§1d/§3a/§3b/§3d).
- The engine to be repaired: `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`.
- The committed output: `research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out`.
- The prior results this repair rests on (these DO state the withheld corrected answer — you may use them to
  check the directive is aimed right and does not leak; the builder never sees them):
  `research/pde_ledger_v3/directives/_measurements/S11c_b_86_reference_result.md` and
  `.../S11c_b_88_blast_radius_result.md`.

## Required method
Open the spec and the engine and verify every premise against the actual bytes. For any claim you contest —
that a line says something else, that a fix is wrong or incomplete, that a control does not bite — show the
exact command (`sed`/`grep`/a short `python3` snippet) and its **literal** output. A prose assertion with no
command behind it is discarded (this binds the review exactly as it binds the builder). Copy anything you
execute to `/tmp` and never modify the working tree.

## Checks (report any failure; each is a way the repair or the directive could be wrong)

1. **Spec faithfulness and the spec-COMPLIANCE claim.** Read §3a lines 251–256, §1d lines 163–168, §3b lines
   276–278. Do they actually require retaining the generated second background jet (Hessian) at `σ_W¹` and
   forbid freezing a coefficient before differentiation, as the directive §1 and premises P1–P3 claim? Or is
   the directive over-reaching — is #89 actually a spec *change* dressed as a compliance repair? If the spec
   does not mandate this, that is the finding.

2. **Code-grounding accuracy (P4–P6).** Verify the cited functions and line numbers: that `basis_fields`
   (~L1025) carries no spurion and is what `basis_euler_signatures` is called with (~L1272); that
   `operator_from_density` (~L1459) uses the global `dx` (~L616) whose `DERIVATIVE_MAP` (~L612–613) has no
   `grad_W→second-jet` entry so `dx(grad_W[i])=0`; and that `operator_dx` (~L1850) / `background_dx`
   (~L2122) already add the Hessian entries. Flag any citation that does not match the current file.

3. **Completeness of the fix scope — the highest-value check.** Are Path A (basis quotient) and Path B
   (strong rows) the ONLY frozen total-divergence paths that bear on the §3a energy basis and the §3b
   strong operator? Search the engine for every use of the global `dx`/`divergence`/`euler_derivative` and
   every constant-coefficient total-divergence, and decide whether any *other* path that feeds the emitted
   basis or the strong rows also freezes the spurion and is silently omitted from the directive (e.g. a
   material-pullback route, the θ `mu_theta_amplitude`, a uniform-vs-new interaction). Conversely, flag any
   path the directive tells the builder to change that should be left alone (the directive says the UNIFORM
   family and the already-Hessian-aware `operator_dx` routes stay as they are — is that correct, or would
   the repair double-correct the coupling/admissibility routes that already retain the Hessian?).

4. **Answer leak (rule 5).** The corrected basis count and per-source count and the corrected operator
   structure must NOT appear as a target in the directive or its premises; the only target may be the frozen
   public `26`. Grep both files for any live count or corrected-operator value and confirm none leaks. Also
   judge: could the frozen-switch `26` regression, or any wording, let a builder infer and tune toward the
   withheld value? (It should not — the corrected count is not derivable from `26` alone.)

5. **Are the controls decisive and falsifiable?** The frozen-switch regression tests only the *frozen limit*
   — it does NOT verify the live path is correct, and the directive §8 says so. Confirm the directive does
   not present the frozen regression as verification of the corrected physics. Confirm the FORM ablation is
   decisive (zeroing the second-jet map must MOVE the basis count; a count byte-identical with and without
   the live Hessian means the repair did nothing), and the COEFFICIENT control is genuinely distinct (rank/
   count invariant under a nonzero rescale). Flag any control that would pass with the repair still frozen,
   or any that is tautological (zero for any input).

6. **Bookkeeping / retained-grade consistency.** Verify the directive's claim that the generated Hessian is
   `σ_W¹` and is KEPT by `first_shape_series` (~L713–725, keeps `η≤1 ∧ σ_W≤1`). Could the prescribed fix
   change the background-amplitude order, double-count a factor, or generate a `σ_W²` term that the
   truncation would silently drop (making the "fix" invisible)? Check the `operator_dx` atoms
   (`σ_W·w1_profile_dij/L_W`) against `PROFILE_GRADE_SUBS` (~L662–671) for a consistent `σ_W` count.

7. **Under/over-specification (rules 3, 12).** Does the directive name the object (the variable-coefficient
   divergence that differentiates the spurion) without prescribing a brittle exact edit that could bias the
   builder toward a wrong implementation? Is anything ambiguous enough that a builder could satisfy the
   letter while freezing a different quantity, or narrow the anchorings/representatives to force
   convergence? The directive forbids narrowing scope — is that enforceable from what it says?

8. **Does the repair, as scoped, actually restore what #88 measured?** #88 found the correction disturbs the
   strong rows non-absorbably via the background Hessian. Confirm the directive's Path B change is exactly
   the operation whose absence #88 witnessed (retain `∂²W_bg` in the strong EL rows), so the repaired engine
   will produce the structure #88 predicted — not a different or partial correction. If the directive would
   fix the basis count but leave the strong rows frozen (or vice versa), that split is the finding.

## Physics filter
Report a finding only if it catches a way the repaired engine would compute the wrong physics, or a way the
directive would mislead the builder into doing so. Do not report style, or "the script would be wrong on a
different input." Prefer findings backed by a command and its literal output.

## Output
Numbered discrepancies (file:line — problem — correction), then `VERDICT: CLEAR` or
`VERDICT: N defects`. Nothing else.
