# Authoring task (Codex-sol): produce v3 of the N6 reconcile-collapse build directive

You are AUTHORING the next version of an orchestrator build directive (rule-15 author change: the orchestrator's
hand-written v1→v2 dictionary drew two convergent NOT-SOUND decision-gate rounds; you author v3). This is a build
directive for a CAS collapse-test instrument — NOT the instrument itself. ⛔ Do not write the instrument code;
write the DIRECTIVE. Ground every frozen table/scale in the actual engine sites (cite file:line). Produce a
directive that a fresh-Claude + Grok review pair will find SOUND.

## Write your output to this file (overwrite it)
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md`
Keep the same overall structure/sections as the current v2; change the version marker to v3 and add a one-line
"v3 — Codex-authored (rule 15); folds round-2 MUST findings" note. Then print a short summary of exactly what
you changed and the engine sites you grounded each frozen table in. ⛔ Do NOT edit any other file.

## Read first (form the frozen tables from the engines yourself)
- Current directive (your base): `directives/S11c_c2_N6_reconcile_collapse_build_directive.md` (v2).
- The two round-2 decision-leg reports you are folding:
  `directives/_legs/S11c_c2_N6_reconcile_collapse_directive_codex_r2.md`,
  `directives/_legs/S11c_c2_N6_reconcile_collapse_directive_grok_r2.md`.
- The framing (the QUESTION this instrument implements): `_measurements/S11c_c2_N6_reconcile_question.md`.
- Engines (ground the maps ONLY here): `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py`,
  `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`, `scripts/S11c_b_brane_operator_sympy_audit.py`,
  `mathematica/S11c_c2_N6_mathematica_audit.wl`. Extraction reuse: `scripts/S11c_c2_N6_cross_engine_comparator.py`.

## The round-2 MUST findings you MUST fold (verify each against the engines, then freeze it)
- **M-a: FREEZE map 4 (energy-basis coefficient pairing) as an explicit table.** The pairing carries physics
  normalization/scale, not mechanical spelling. Supply the COMPLETE contraction-ID table: for each retained
  energy contraction, the WL coefficient (e.g. `bRho`, `cCoupling`, `energyCoefficient<i>` — `.wl:211,217,222-226`)
  ↔ the SymPy coefficient (e.g. `B_rho_3`, `C`, `kappa_theta`, first-jet `gamma_s11cb_*` — `brane:1584-1594,357,1820,1727-1872`),
  WITH the exact leftover scale/orientation (e.g. the `W_0`/`WBg` and `1/2` factors, `B_rho_3 = bRho·W_0`, etc.).
  Derive the scales from the two energy-density constructors yourself and cite them. If a fully static table is
  not derivable, instead specify a DETERMINISTIC basis-change algorithm (with its expected pair table) that a
  reviewer can check — but a frozen table is preferred. ⛔ It must be a coefficient-name+scale correspondence,
  ⛔ never an equality of the resulting μ objects, ⛔ never a runtime match/solve against residuals.
- **M-b: REMOVE map 6 (leftover density) from the production rewrite dictionary.** Both engines rebind live
  density INSIDE source construction (PY `diagnostic:378,382`; WL `sourceBind` `rhoFace→density3`, `.wl:380,839`),
  so no leftover density atom NAME survives in the emitted operands — the map is either inert or a forbidden
  constant↔live-expression rewrite (c1 rule-17 hazard). Keep the live-density agreement ONLY as a separately
  emitted premise/census comparison, ⛔ not an operand rewrite.
- **M-c: give every map a TYPED family/stage SCOPE.** The source-wave bare↔applied-at-`Y` map (`jet[f,I] ↔
  ∂_I s11cc2Field_f(Y,t)`) applies ONLY to `SOURCE_*` (source extraction inserts applied jets at Y —
  `diagnostic:392`; WL `sourceMap` bare jets `.wl:707`). Φ is an abstract jet map with NO `Y` evaluation
  (`cov:63`; WL `:236`) — Φ receives only domain-name/multi-index/derivative-syntax spelling + the pre-grade
  profile substitution, ⛔ never the source-point-evaluation map. Carriers receive the pressure-slot/jet name
  maps, ⛔ not the source-wave map. Make the scope of each map explicit (which families it may touch).
- **M-d: complete the CONTROLS (locality + dropped-map).** For each ACTIVE map: corrupt it one-sided → require a
  nonzero control-delta on ≥1 declared-use leaf AND require ZERO control-delta OUTSIDE its declared use-set
  (locality — catches an over-broad map). Add a DROPPED-map control: delete each active primitive → require a
  residual on a predeclared defining-relation operand. Run the blanket-collapse ablation through the SAME
  extraction/rewrite/witness path. Keep the reachability census (honest inert-map emit) + the grade-combining
  tripwire. ⛔ Drop any "shipped run must/‑not be all-collapsed" statement.
- **M-e: FREEZE the census crosswalk exactly (no "e.g.").** Enumerate the exact pairings: `kappa_a↔ADVECTION`,
  `kappa_j↔JUNK`, and the domain-census `DOMAIN`/`COVERAGE`/`UNCOVERED`/`MAX_RANK` ↔ their SymPy counterparts
  (`cov:105-115,146-153`; WL `:977-980`). ⛔ `MATERIAL_NORMAL` and `THICKNESS` are WL-ONE-SIDED — retain them as
  surfaced one-sided fields; ⛔ do NOT synthesize a SymPy zero counterpart.

## Keep (both legs found these SOUND — do not regress them)
Retained-order coefficientwise over the 4 independent grades `(η^i σ_W^j)`; the exclusion list (whole-object /
construction-stage / cross-anchoring / σ_W-binding / default on-shell-Fourier); Φ graded via
`extract_meta`→profile-expand→4-coefficient extraction; the three-valued witness (residual printed first; fresh
PIT with primes/seeds disjoint from both engines; denominator rejection; degree bound; max attempts/time/RSS;
undecided-on-exhaustion; provenance emitted; no verdict/assert); BASELINE as nominal-control duplicate; the
import-comparator-primitives + group-by-key extraction (no comparator edit; STOP+report if a helper is needed);
the three clauses (PRINT never assert; residual printed not asserted; interpretation in the record); the builder
fence (build→verify→run→report→STOP; no self-review; no outside edits); the DoD + builder-report shape.

## ⛔ Value-free (hard constraint)
The directive states NO expected collapse outcome, NO expected count, NO "the residual is zero/nonzero", NO pass
condition, NO statement about whether the shipped run collapses. Prior-run nonzero counts are historical facts,
⛔ not targets. The map-4 SCALE factors (`W_0`, `1/2`, etc.) are engine-construction facts (coefficient
definitions), ⛔ not collapse outcomes — supplying them is correct and is not a leak.

## Output
1. Overwrite the directive file with v3.
2. Print: the list of changes; the exact map-4 frozen table with the engine sites you derived each scale from;
   confirmation map 6 is removed from rewrites; the typed scopes; the completed controls; the frozen census
   crosswalk. ⛔ No physics conclusions about whether anything collapses.
