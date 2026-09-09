# Build directive — S11c-c2 N6 reconcile COLLAPSE instrument (v2, folded from 2 decision legs)

⚠ **Orchestrator-written build directive.** It governs the CODE astra writes; astra does NOT decide physics.
This instrument implements the reconcile collapse TEST framed + vetted in
`_measurements/S11c_c2_N6_reconcile_question.md`. ⛔ The instrument PRINTS; it decides NOTHING.

⭐ **v2 — folded from 2 CONVERGENT decision legs** (Codex-sol + Grok, both DIRECTIVE NOT-SOUND; record
`_measurements/S11c_c2_N6_reconcile_collapse_directive_gate.md`; leg reports in `_legs/`). The v1 dictionary
was NOT executable/frozen and carried three computation-changing hazards (a degree-2 construction projector, a
profile η-re-expansion, and a whole velocity/source-factor equality) plus gaps (ungraded Φ, a missing
source-wave map, tautological controls, no census crosswalk, no operand-pair API, no PIT bounds). This v2 fixes
all of it. ⭐ **The frozen dictionary is now specified as ALLOWED MAP TYPES + exact ENGINE SITES + a per-map
defining-relation control; astra CONSTRUCTS the explicit rewrite tables FROM those sites and FREEZES + EMITS
each as an object** — ⛔ astra does not invent maps, ⛔ the orchestrator does not hand-type the symbol tables
(the route-2 failure mode). The re-review legs verify each constructed map is a mechanical name/structure fact
that bites.

## Deliverable
`research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_collapse.py` (+ a `test_…` alongside).

## Object
Per surfaced-nonzero cross-engine family, per matched key, per retained grade `(η^i σ_W^j), i,j∈{0,1}`, per
fixed `(anchoring,density)` case: apply ONE frozen bridge dictionary of **mechanical name/structure
correspondences** (§ Dictionary) UNCHANGED to the two engines' emitted operands, then PRINT the bridged residual
`bridge(operand_SymPy) − bridge(operand_WL)` and a three-valued collapse witness. ⛔ Never a verdict, ⛔ never an
assertion.

Surfaced families tested:
- `N6RC_CARRIER_EULERIAN`, `N6RC_CARRIER_MATERIAL` (channel a);
- `N6COV_SOURCE_ACTUAL`, `N6COV_SOURCE_PREDICTED`, `N6COV_FROZEN_PHI` (channel b — Φ IS a collapse operand,
  see § Φ);
- `N6COV_SOURCE_BASELINE` — computed but LABELLED the nominal-control DUPLICATE of `SOURCE_ACTUAL`
  (`kappa_a=1, kappa_j=0` ⇒ μ_actual=μ_baseline; `cov:34-35,140-154`; WL `:18-20,845-849`; the 160/0
  `SOURCE_CONTROL_DELTA` confirms), ⛔ NOT a third independent discriminator.
Plus a control/premise census (`N6COV_ACTUAL_CONTROL_PARAMETERS`, `N6COV_PHI_DOMAIN_CENSUS` — § Census).

## Inputs (all committed; content present on disk)
- `DEFAULT_PY` = the 3 SymPy `.out` (`scripts/out/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.out`,
  `7c0790ab`); `DEFAULT_WL` = `mathematica/out/S11c_c2_N6_mathematica_audit.out` (`ae73b884`).
- The two engines + shared module (the maps are grounded ONLY in these):
  `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py`,
  `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`, `scripts/S11c_b_brane_operator_sympy_audit.py`,
  `mathematica/S11c_c2_N6_mathematica_audit.wl`.
- `scripts/S11c_c2_N6_cross_engine_comparator.py` — reuse its extraction PRIMITIVES (below).

## Operand-pair extraction (⛔ do NOT re-implement the join/schema bridge; ⛔ do NOT edit the comparator)
`compare_family` materializes operands, computes the residual, emits `CASE`, then RELEASES them
(`comparator:794-825`); `object_work` returns only accounting (`:862-898`). So you CANNOT get operand pairs by
"reusing the CASE-line path." Instead **import the comparator's extraction primitives** — `load_py_jsonl`,
`load_wl`, `verified_spelling_maps`/`ACTIVE_NAMES`/`POINT_NAMES`, the leaf keying + `materialize` — and, SCOPED
to the surfaced families, group matched leaves by their typed key to obtain, per key, the SymPy operand and the
WL operand as CAS objects BEFORE any residual. The group-by-key uses the comparator's OWN keying + mechanical
spelling maps (so keying/schema-bridge stay identical); the ONLY new content is the physics name/structure
dictionary + the witness. ⛔ If a needed primitive cannot be imported without editing the comparator, STOP and
report — a comparator helper is a SEPARATE reviewed patch (its own gate), ⛔ not an in-band edit.

## THE FROZEN DICTIONARY — ALLOWED MAP TYPES + SITES (astra constructs, freezes, and EMITS each table)
⛔⛔ ONLY **mechanical name/structure correspondences between the two engines' emitted atoms**, each grounded in
a named engine site below. ⛔ NO whole-object equality (carrier/μ/source/Φ), ⛔ NO construction-stage operation
(the engines already performed those; applying them to the surfaced post-construction leaves is either vacuous
or reintroduces η), ⛔ NO map that identifies independently-computed subobjects. astra builds each table from its
site, freezes it, and EMITS it as an object `{map_id, domain→image pairs, engine_site}` so the re-review legs
verify it is a naming/structure fact. Applied UNCHANGED to `SOURCE_ACTUAL`, `SOURCE_PREDICTED`, `C_E`, `C_M`,
**and Φ**.

1. **grad-θ jet duals** — PY geometry `theta_d{i}` (`S11c_a…:140`) ↔ operator `grad_theta_{i}` (`S11c_b…:227-228`)
   via `wave_jet` (`selfenergy:146-147`); WL `thetaD_i`/`gradTheta_i` → `jet["theta",{i}]` (`.wl:946-948`). The
   comparator already maps WL `thetaJet{i}`→PY `theta_d{i}` (`comparator:290-311`); add the PY dual name only.
2. **source-wave point-preserving map (MUST — was missing)** — WL `sourceMap` retains BARE registered jet atoms
   (`.wl:707-711`); SymPy `source_value` inserts APPLIED wave jets at `Y` (`diagnostic:392-396`; `wave_jet`
   builds applied fields+derivs `selfenergy:139-169`). The comparator deliberately does NOT identify bare with
   applied fields (`comparator:287-289,990-993`). Construct the exact POINT-PRESERVING correspondence
   `jet[f,I] ↔ ∂_I s11cc2Field_f(Y,t)` (bare atom ↔ applied field derivative at the SAME point Y), for the
   source jets that appear.
3. **profile-jet name/scaling** — leftover profile-jet NAMES + scaling factors only: `w1ProfileJet{ij}` ↔
   `w1_profile_d{ij}`, the `m1`-profile analogue, and the `σ_W/L_W`, `μ_R/W_0` leftover scale factors
   (`.wl:127-132`; `selfenergy:281-301`; `brane:767-821`). ⛔ NO `WBg→W0(1+η·w1)` / `muRBg→muR(1+η·m1)`
   re-expansion — both engines expand profiles BEFORE grade extraction (WL `:136`→`:141`; PY
   `diagnostic:247-248`); re-expanding after grading reintroduces η into a single-grade leaf (forbidden — §
   Retained order).
4. **energy-basis coefficient-name pairing TABLE** — an INJECTIVE table pairing WL quotient-basis coefficient
   names (`.wl:208,222-226`) with SymPy termwise-EL monomial coefficients (`diagnostic:318-349`;
   `brane:1727-1872`), by CONTRACTION IDENTITY. ⛔ NOT a runtime match/solve against μ/source residuals; ⛔ NOT
   equating the resulting μ objects. Emit the pair table + an injectivity/coverage check as an object.
5. **Φ domain-name / multi-index / derivative-syntax spelling** — spelling of the map's DOMAIN atom names +
   sorted-multi-index + derivative syntax only (`cov:63-129`; `.wl:236-247`). ⛔ The Φ VALUES stay in the
   residual (their equality is what channel b TESTS) — ⛔ never a whole-Φ equality.
6. **leftover density name** — after BOTH engines rebind live density, the leftover-name identification PY
   `rho_br_bg_rho4_constant`/`inputs.density[(ρ,)][1]` (`diagnostic:382`) ↔ WL `density3` (=`density4·WBg`,
   `.wl:839-840`, used in `sourceBind` `:381`). ⛔ A leftover-NAME identification only; ⛔ NOT a constant↔live-
   field equality (c1 made re-adjudication mandatory, `S11c_c1_comparator_reconcile.md:137-163`; it did NOT
   authorize an unconditional constant↔field equality).
7. **Jacobian spelling** — `1 + tr(∇u)` ↔ `1 + Sum jet[u_i,{i}]` (`reconcile:65-66`; `.wl:246`), applied ONLY as
   a leftover-FACTOR identification where that factor still appears in a leaf. ⛔ NO wave-projection "degree 2"
   projector — that is a construction stage inside μ_M (PY `brane:1970-1981`; WL `waveScale[…,2]` `:246`);
   emitted sources are already wave-LINEAR (PY `wave_terms` one wave `diagnostic:368-374`; WL `waveScale[…,1]`
   `:381`), so a degree-2 substitution would project them to 0 (vacuous collapse).
8. **c1 conditional atom maps** — ε-placement + `omega` realness ONLY on a leaf whose expression actually
   contains those atoms after maps 1-7 (split into exact conditional rules; ⛔ omit any rule with no occurrence).
   ⛔ NO default on-shell / Fourier-of-derivative — those are c1 KERNEL identities absent from N6 source
   construction (PY strips DtN/resolvent before source extraction `diagnostic:383`; WL kernel begins after
   `sourceBind` `:380` vs `:385-418`).

⛔ **Excluded (each would make a collapse vacuous or leave the retained rectangle):** whole-carrier/μ/source/Φ
equality; any construction-stage operation as an operand rewrite (degree-2 projection; profile η-re-expansion;
whole velocity or source-solve factor equality — SymPy builds material velocity + source factors independently
`reconcile:82-92`, `diagnostic:378-389`, and WL derives its own `.wl:274,304,366-381,862-863`, and RESOLVED.md
keeps face-velocity correctness OPEN `:41-43`, so these stay in the residual); `LAB_HELD ↔ MATERIAL_ADVECTED`
(anchoring is a retained comparator key `comparator:86`); `σ_W↔η` / `σ_W→0`; default on-shell/Fourier.

## Φ as a graded collapse operand (MUST — Φ is ungraded metadata)
`FROZEN_PHI` is emitted as an UNGRADED substitution map (PY `cov:116-124`; WL `:975-976`) and the comparator
materializes it via `extract_meta` as `MAP_VARIABLE`/`FIELD_PATH` keys with NO grade axes (`comparator:688-765`).
So for Φ the instrument must: extract the map images via `extract_meta`, apply the profile expansion, then
EXPLICITLY extract all four `(η^i σ_W^j)` coefficients, and only then apply the dictionary + witness
coefficientwise. ⛔ Do NOT compare ungraded Φ images; ⛔ do NOT reuse the numeric CASE-line grade path for Φ.

## Retained order (⛔ coefficientwise; ⛔ no proxy; ⛔ no η re-introduction)
Operate per retained grade `(η^i σ_W^j), i,j∈{0,1}` independently (numeric leaves are already grade-keyed;
Φ is graded per § Φ). ⛔ No grade-combining, ⛔ no cross-grade cancellation, ⛔ no `η·σ_W→η²`, ⛔ no `σ_W→0` or
`σ_W↔η` binding, ⛔ no map that reintroduces `η` into a graded leaf. A physical σ/η relation may be reported as
a SEPARATE object after the coefficientwise comparison; ⛔ never applied during it. Enforce this structurally (a
grade-combining tripwire).

## Collapse witness (three-valued; PRINT residual first)
Per leaf, after the dictionary: PRINT the bridged residual, THEN a bounded witness ∈
{`collapsed` (no nonzero found), `residual_remains` (computed nonzero, printed full or as its arithmetic-DAG),
`undecided` (budget/parse exhaustion)}:
- attempt a bounded symbolic zero-test (`cancel`/`together`/`factor`, per-leaf budget);
- else a FRESH joint finite-field PIT with the instrument's OWN disjoint primes + seeds — ⛔ NOT the engines'
  primes (PY `(1000000009,998244353,1000000033)` `diagnostic:669`; WL `(1000000009,998244353,1004535809)`
  `.wl:733`); with denominator rejection, a derived degree/exclusion bound, and enforced max attempts / time /
  RSS; `undecided` on exhaustion. Emit the PIT provenance (primes, seeds, degree bound, draws, conditional δ) as
  an object.
⛔ The witness is NOT a pass/fail token and carries NO interpretation.

## MANDATORY controls (the re-review legs will ablate these; build them so they BITE, and demonstrate by running
the shipped controls — ⛔ NOT by a self-review)
1. **Per-map defining-relation control** — for EACH constructed map, corrupting its defining relation (in a
   /tmp copy) must move a nonzero **CONTROL DELTA** (baseline bridged residual − corrupted bridged residual) on
   the leaves that use it. ⛔ NOT merely "the corrupted residual is nonzero" (a leaf already nonzero satisfies
   that vacuously). This is the per-map load-bearing + independence test.
2. **Production-reachability census** — for EACH map, emit (as an object) whether it actually occurs in any
   surfaced leaf. Some maps may have NO surviving occurrence (the operation was done inside the engine) — that
   is HONEST and emitted, ⛔ NOT forced; a map with no occurrence is simply inert on these leaves.
3. **Blanket-collapse ablation** — adding a whole-object equality (identify the two engines' operands outright)
   in a /tmp copy must force `collapsed` on every leaf. This proves the shipped dictionary is not doing that.
   ⛔ Do NOT state or imply any expectation about whether the SHIPPED run is all-collapsed or not (no such
   target; the honest shipped run's collapse pattern is whatever it is).
4. **Grade-combining tripwire** — an attempt to reintroduce `η` into a graded leaf or bind `σ_W` must be
   structurally impossible or hard-flagged.
⭐ For every control, emit BOTH operands + the residual + the witness/delta; ⛔ never assert.

## Census audit (channel d — ⛔ LOAD-BEARING control/premise, not dismissible metadata)
Predeclare a **census-only semantic crosswalk** grounded in both engines (e.g. `kappa_a↔ADVECTION`,
`kappa_j↔JUNK`, `max_present_jet_rank↔MAX_RANK`, material-normal↔`MATERIAL_NORMAL`; PY `cov:146-153`; WL
`:977-980`; the comparator adds NO such aliases `comparator:683-685`). Emit, per matched census key: SymPy value,
WL value, and a structural production-control / domain-coverage equivalence object. RETAIN unmatched/one-sided
fields (e.g. WL `THICKNESS` that PY does not emit — PY `actual_amplitudes` retags only theta-advection+junk
`cov:137-154`) as SURFACED objects. ⛔ Do NOT dismiss as metadata, ⛔ do NOT pre-decide; a differing μ-jet
coverage routes to channel b.

## ⭐⭐⭐ THE THREE CLAUSES (verbatim — non-negotiable)
> **1. The script may PRINT computed objects. It may NOT state conclusions.** An emit payload must be a CAS
> object (an expression, a residual, a symbolic zero-witness), ⛔ never prose describing a result.
> **2. PRINT the residual; do NOT assert it.** ⛔ No `assert residual == 0`, ⛔ no residual-zero exit.
> **3. Interpretation belongs to the reconcile record.** ⛔ The script does not editorialise.

⭐ The ONLY hand-combination of physical symbols allowed is the dictionary's name/structure correspondences
(each grounded in a named engine site + verified by its defining-relation control); every bridged residual +
witness is REACHED BY COMPUTATION from the emitted operands + the frozen tables.

## Value-free / leak discipline
⛔ NO expected collapse outcome, NO expected count, NO "the residual is zero/nonzero", NO pass condition, NO
statement about whether the shipped run collapses. The prior-run nonzero counts are facts of the committed
comparator run, ⛔ NOT targets — do not reproduce them. The builder must see only "apply the frozen name/
structure dictionary and PRINT what results."

## Builder fence (astra) — build → verify → run → report → STOP
⛔ Build the instrument, verify the deliverable exists + is non-empty + its `test_…` passes, run it SCOPED to the
surfaced families over the 4 cases (bounded per-object budgets like the comparator's
`--object-seconds/--object-rss-mib`; watch RSS), produce the `.out`, then REPORT. ⛔ Do NOT author or run any
review of your own work; ⛔ do NOT call other AI; ⛔ do NOT edit any file outside the deliverable + its test +
its `.out`. The independent review is the orchestrator's to launch.

## Definition of Done
- `scripts/S11c_c2_N6_reconcile_collapse.py` exists; imports the comparator's extraction PRIMITIVES + groups
  matched operand pairs by key (no comparator edit); applies ONLY the frozen name/structure dictionary (each map
  EMITTED as `{map_id, pairs, site}`; ⛔ no whole-object / construction / η-reintroducing map reachable);
  operates coefficientwise per grade at retained order; Φ graded via `extract_meta`+profile-expand+4-coefficient
  extraction; PRINTS per leaf both operands, applied-maps, bridged residual, and the three-valued witness.
- Each map's defining-relation control moves a nonzero control-delta on its leaves; the production-reachability
  census is emitted; the blanket-collapse ablation forces all-collapsed; the grade-combining tripwire holds.
- The census crosswalk + structural equivalence objects for the 8 control/premise leaves are emitted; unmatched
  fields retained.
- PIT bounds (disjoint primes/seeds, denominator rejection, degree bound, max attempts/time/RSS, undecided-on-
  exhaustion) enforced; PIT provenance emitted.
- `test_…` covers: dictionary is name/structure-only (no whole-object/construction map reachable), coefficientwise
  enforcement + grade-combining tripwire, witness three-valuedness + native-boolean rejection, Φ graded path,
  extraction reuse (no comparator edit), no-assert / no-verdict.
- A scoped `.out` is produced. ⛔ No verdict token, ⛔ no assert-before-emit, ⛔ no expected value anywhere.

## Builder report (return this)
Deliverable + test paths; token usage; the emitted per-map site table + production-reachability census; the run
accounting (families, keys, grades, NEUTRAL per-witness tally collapsed/residual_remains/undecided, ⛔ not an
interpretation); the per-map control-delta results; peak RSS + runtime; confirmation the fence held. ⛔ No
physics conclusions.
