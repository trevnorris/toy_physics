# v2 CONFIRM-PASS. Your prior review gave 8 fixes (central: the sign depends on the OPEN throat-shear coupling zeta -> compute the sign as a function of zeta, do not assume generic point-vortex dynamics; also: O(V^2) bilinear; radial-force vs trajectory are DIFFERENT observables; drop the imported static-charge-repulsion fixture; fix ghost transforms so reversing both circulations leaves k1*k2 invariant; anti-import as a provenance dependency-graph; expand the verdict set). v2 folds all 8. CONFIRM they landed and flag anything still wrong or under-specified. If clean say DIRECTIVE_CLEAN.

# Design-review the directive (do NOT execute it) — native two-throat magnetism SIGN gate

You are design-reviewing a directive, not running it. Toy 4D superfluid analog program. Falsification is the goal; a computed departure (native magnetism repels/co-rotates instead of Ampere-attract) is a first-class, welcome result.

## Read
- software/em_charge_attribute/directive_native_vortex_sign.md  (THE directive under review)
- docs/conceptual_foundation.md sec 3-4 (magnetism = moving throat's 4D-body swirl, +/-w; gravitomagnetism vs EM-magnetism distinction)
- docs/em_u1_body_definition.md (E4/E5 roll-vs-slip boundary-operator endpoints — supersedes the retired em_gravity_native_ontology.md §9)

## The question
The force_visualizer sim shows two parallel currents (like-circulation swirls) attracting = correct EM. But pathA_39 got that sign by IMPORTING an EM current source (j=sV) + Maxwell exchange, not from native moving-throat/vortex dynamics. This gate computes the sign natively (Magnus/Berry/inertia/derived-flow), with the honest possibility it gives repel/co-rotate (a falsifier). Two like superfluid vortices co-rotate; parallel currents attract; these are opposite — that is the crux.

## Review for (concise, able-to-fail)
1. Is the native vortex-dynamics setup physically correct and well-posed (Magnus force F = rho*kappa x V_rel; collective-coordinate two-body interaction; core inertia; the derived-flow-not-inserted-current requirement)? Anything mis-stated?
2. Are the four controls' KNOWN answers actually correct as stated — vortex-antivortex co-translating dipole; two like vortices co-rotate with no radial attraction (ideal 2D); single vortex Magnus deflection in uniform flow; static like charges repel? Any wrong, and is a needed control missing?
3. Is the anti-import guard (must FAIL if j=sV or an asserted Maxwell sign is used as the source) genuinely enforceable as a computed check, not a source-grep?
4. Can the pipeline genuinely return ATTRACT / REPEL / CO-ROTATE (three qualitatively different outcomes), or is it accidentally biased to one? Is CO-ROTATE (neither attract nor repel) properly handled as a distinct verdict?
5. Tautology / false-pass / false-fail risks; is the verdict set exhaustive and unambiguous?
6. Scope honesty (collective-coordinate, not full nonlinear sim; EM-magnetism vs gravitomagnetism kept distinct).

Do NOT design the scripts. Critique the directive.

## Output
- VERDICT: DIRECTIVE_CLEAN or DIRECTIVE_NEEDS_FIXES.
- If fixes: numbered, each with the specific change + why (physics error / under-specification / vacuous control / missing control / scope). Then <=5 bullets reasoning. Concrete over polite.
