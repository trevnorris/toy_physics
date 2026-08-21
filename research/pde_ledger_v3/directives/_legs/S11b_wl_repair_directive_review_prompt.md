# Directive-design review — S11b WL engine repair directive

You are one of two independent legs reviewing an ORCHESTRATOR-written repair directive BEFORE the builder
runs. Find defects: a fix that won't close its defect, a fix that OPENS a new one, a rule-5 leak, or an
acceptance that a bad fix survives. A leg that finds nothing is weak — show the checks you ran.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_wl_engine_repair_directive.md`

## Sources of truth — read first, then the directive
- The engine being repaired: `research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl`
  (F-WL-1 `ZPERM_SLICE_MAP` ~L865 + `equationZeroForm` L367; F-WL-2 `PRESSURE_WORK_SIGN_CHECK` ~L1176 +
  `energyEquationRules` ~L1118; F-WL-3 `KERNEL_ORIENTATION_IDENTITIES` ~L1102, `CAUSALITY_CHECK` ~L1105,
  `GRAZING_MODE_CLASSIFICATION` ~L1368).
- Physics authority: `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` (§4 face closure / affinity;
  §6 energy accounting; §8 no-verdict + locus protocol; B2c).
- The adjudicated findings: `research/pde_ledger_v3/steps/S11b_wl_engine_review_disposition.md`.

## Checks
1. **F-WL-1 correctness of the fix.** Does "extract `Λ_p⁰` as the coefficient of `facePressure` in the FLUX
   RHS (not the residual zero-form), as the static (ω-independent) map" actually yield the coefficient §4's
   affinity forces on the raw-pressure channel? Derive it yourself from §4 (write+run a script; paste stdout)
   and confirm the prescription is right WITHOUT the directive stating the value. Is the static-vs-dynamic
   split unambiguous?
2. **F-WL-2 / F-WL-3 fixes.** For each, will the prescription make the check GENUINE (able to fail under a
   form change), and is the "or emit honestly (NOT_ESTABLISHED)" escape correctly bounded (never a tautology
   dressed as a check)? Could any fix open a NEW tautology or hide a real defect? Is the grazing
   Laurent/nullspace prescription implementable, or under-specified?
3. **Rule-5 leak.** Does the directive state or let a builder derive any withheld value — the F-WL-1 sign/
   relation, the basis count, an expected residual? Quote any leak. (⚠ The F-WL-1 correct value is withheld.)
4. **Acceptance able-to-fail.** Is the FORM-ablation acceptance decisive (a bad fix fails it)? Does it avoid
   demanding the withheld value?
5. **Scope / no-regression.** Does it adequately bar behavioural change to the correct derived physics
   (impedance, added mass, ZPERM_SLICE, transverse, breathing, roots)?

## Method
Document review; you may run small SymPy/WL scripts to check a claim (save script+stdout to named absolute
paths, report them). ⛔ Do not edit the working tree.

## Report — under ~20 lines
Numbered findings (directive line, the failure it permits, the concrete fix), then a one-line verdict: safe
to build from, or fold first?
