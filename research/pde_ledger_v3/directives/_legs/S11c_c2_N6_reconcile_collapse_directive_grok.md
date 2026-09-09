I will review the directive as an independent decision leg: engines and framing first, then the frozen bridge dictionary and collapse-witness design against what those engines actually compute.I’ll form the maps from the engines and comparator first, then read the directive. Starting with the production scripts and the shared diagnostic.The question framing is in hand. Next I’ll read the engine construction sites those maps are supposed to rest on.I have the engine construction paths. Next I’ll pin the comparator extraction, the c1 bridge precedent, and the retained-order/profile conventions.The directive asks the right *kind* of question (primitive bridge, retained rectangle, print-don’t-decide), but as written it would make the builder apply two construction-time identities as post-extraction substitutions, implement map 6 from the wrong WL objects, and mishandle Φ’s ungraded metadata leaves. Those change what the instrument computes.

Engines were read first. Maps below are derived from those sources, then checked against the directive.

---

## Independent map table (directive list vs engines)

| # | Directive map | Primitive? | Grounded at cited lines? | Collapse-safe? |
|---|---|---|---|---|
| 1 | jet-vocabulary `a.grad_theta[i]↔b.grad_theta[i]` via `wave_jet` | Yes — symbol convention. PY geometry `theta_d{i}` (`S11c_a…:140`) vs operator `grad_theta_{i}` (`S11c_b…:227-228`); `wave_jet` rewrites `grad_theta_*` → `theta_d*` (`selfenergy_fold:146-147`). WL `thetaD_i`/`gradTheta_i` → `jet["theta",{i}]` (`:946-948`). Comparator spelling already sends WL `thetaJet{i}` → PY `theta_d{i}` (`comparator:290-311`). Remaining work is the PY dual name. | Yes (`reconcile:59-69`; WL `:946-948`) | Yes |
| 2 | Jacobian `1+tr(∇u)` **and wave-projection degree 2** | Jacobian spelling is primitive. Degree-2 projection is **not** a map on the surfaced leaves. | Citations exist (`reconcile:65-66`; WL `:246`) but name a **construction** already baked into μ_M | **No — vacuous on SOURCE/Φ** |
| 3 | live-density rebind | Yes — c1 seal-5 field-vs-constant, mandatory in c2. PY already rebinds `rho_br_bg_rho4_constant` → `inputs.density[(ρ,)][1]` inside `source_terms` (`diagnostic:382`, function starts `:378`). WL `density4` / `density3=density4·WBg` (`:839-840`); `sourceBind` uses `rhoFace→density3` (`:381`). `density3` matches PY `density_pair` 3-density (`S11c_b:1886-1890`). | Function-level cite is slightly early (rebind is `:382`) but the object is right | Yes, as leftover-name identification, not a second rebind |
| 4 | `WBg→W0(1+η·w1)`, `muRBg→muR(1+η·m1)` + σ-scaled jets | Primitive **before** grade extraction. Both engines already expand then grade: WL `finish` applies `profileRules[]` (`:127-132,136`) then `gradePart` (`:141`); PY `Compiler.grades` does `xreplace(profiles)` then `shape_coefficients` (`diagnostic:247-248`; profiles from `selfenergy_fold:204-207`, jets `:281-301`). | Yes | **No on pre-graded leaves — reintroduces η** |
| 5 | energy-basis by contraction identity; ⛔ not μ-equality | Primitive **if** it is coefficient-name pairing (`B_rho_3` on `W θ²` ↔ `bRho`, etc.). WL quotient `BASIS` is `:208`, coefficient assignment `:222-226`; `:211` is only `constructEnergy`’s start. PY is termwise EL of S11c-b monomials (`diagnostic:318-349`). Explicit μ-ban is the right whole-object fence. | Weakly (function entry, not the pairing table) | Yes **only** as name-pairing |
| 6 | face-velocity + source-solve factors | Primitive **objects** exist. **Cited WL lines are the wrong objects.** | PY `cov:189` is only Eulerian `V_E`. WL `:359-362` is `massSubstrate`’s virtual-constraint / Jacobian return; `:364` is a comment. Real WL `VELOCITY` is `:274`/`:304`, bound `:862-863`. Real source-solve is `sourceConstruction`/`sourceBind` `:366-381`. PY factors are `diagnostic:388-389` (`mu_factor`, `v_factor`, `velocity/eps`). | Builder following the cites bridges mass-constraint, not V/source-solve |
| 7 | Φ spelling; equality tested not assumed | Yes — name/structure of `θ↦θ+a_ρ`, `e_W↦e_W+h_α`, prolongation via `td`/`total_derivative`. PY `prolonged_phi` (`cov:63-129`); WL `phiMap` (`:236-243`). Whole-object Φ equality is correctly forbidden. | Yes | Yes as spelling; **grade handling is unspecified** (below) |
| 8 | c1 ε / ω-real / FT-of-derivative, only if atoms remain after 1–7 | Yes; same pattern as c1 staged bridge. Default on-shell correctly excluded. | Conditional — no N6-specific lines, by design | Yes |

**Excluded list is right:** whole-carrier/μ/source/Φ equality; `LAB_HELD↔MATERIAL_ADVECTED` (anchoring is a comparator key, `comparator:86-91,103`); `σ_W↔η` / `σ_W→0`; default on-shell/Fourier.

**Retained rectangle is real in both engines:** PY `GRADES = product((0,1), repeat=2)` (`diagnostic:43`), kept in `grades` (`:248`) and `gmul` (`:181-182`); WL `gradeIndices = Tuples[{Range[0,1],Range[0,1]}]` + `Series[...,{etaBg,0,1},{sigmaW,0,1}]` (`:123-124`); comparator rejects `ETA,SIGMA∉{0,1}` (`:96-97`). No engine binds `σ_W=η·W0/L_W` at this order (that was c1’s numeric test, and is correctly forbidden here).

---

## MUST findings

### M1. Map 2’s wave-projection degree 2 would vacuous-zero surfaced SOURCE (and Φ) leaves

Both engines already take the quadratic-wave piece **while building μ_M**, not while emitting sources:

- PY `material_pullback`: `(jacobian * pulled)` then `diff(scaled, scale, 2)/2` (`S11c_b:1971-1981`).
- WL `materialAmplitude`: `waveScale[(1+Sum u_{i,i}) (density/.MAP), 2]` (`:246`).

Emitted sources are already **wave-linear**: PY `wave_terms` requires exactly one wave symbol (`diagnostic:368-374`); WL `sourceBind` uses `waveScale[..., 1]` (`:381`). Φ images are `θ+a`, `e_W+h` and their jets (`cov:75-76`; WL `:238-239`), not quadratic densities. Carriers are pressure-slot derivatives of rows (`reconcile:410-412`; WL `carrier` `:336`).

The dictionary is applied **unchanged** to `SOURCE_*` and `C_E`/`C_M` (directive:68). Map 2 does **not** have map 8’s “only if the atom remains” restriction. A builder who applies “wave-projection degree 2” as a substitution on those leaves projects degree-1 sources / non-energy Φ images to 0, so `bridge(A)-bridge(B)` is `0-0`. That is a blanket collapse of channel (b) by a construction filter, i.e. the §3 question answered by deleting the operands.

The Jacobian **spelling** `1+tr(∇u)` ↔ `1+Sum jet[u_i,{i}]` is a legitimate leftover-factor identification **if that factor still appears**. It must not be bundled with a degree-2 projector.

### M2. Map 4 as `WBg→W0(1+η·w1)` on already-graded leaves reintroduces η (forbidden grade-combining)

Directive retained-order clause: operands are **already graded** by the comparator key; no `η·σ_W→η²`, no `σ_W↔η` (directive:71-75).

Both engines expand profiles **before** taking `(η^i σ_W^j)`:

- WL `:136` then `:141`
- PY `diagnostic:247-248`

Substituting `WBg→W0(1+η·w1)` into a grade-`(0,*)` coefficient puts a new `η` into a leaf that is supposed to be the coefficient of a single grade. That leaves the independent rectangle — the same class of defect as the excluded `σ_W` binding.

The primitive that **is** needed is leftover **name/scaling** of profile jets (`w1_profile_d{ij}` ↔ `w1ProfileJet{ij}`, `σ_W/L_W` and `μ_R/W_0` factors in `selfenergy_fold:281-301` / WL `:127-132`), applied **without** re-expanding `η` after grade extraction.

### M3. Map 6 is not grounded at the cited WL (or PY source-solve) lines — the builder would map the wrong objects

Directive: “face-velocity + source-solve factors (`cov:189`; WL `:359,364`).”

What those lines actually are:

- `cov:189` — PY **Eulerian** `inputs.geometry['face_velocity']` only. Material `V` is `:190` (`build_material_velocity`). Source-solve **factors** are `mu_factor=∂src/∂μ_slot`, `v_factor=∂src/∂V_slot`, and `velocity/eps` at `diagnostic:388-389`.
- WL `:359-362` — `massSubstrate` return: evolution bases, **virtual constraints**, `DENSITY_JACOBIAN`. `:364` is the comment *before* `sourceConstruction`.

Actual WL face velocity: `graphGeometry`/`materialGeometry` `"VELOCITY"` (`:274`, `:304`), stored `:862-863`. Actual source-solve: `sourceConstruction` driving at `pFace→0` (`:366-375`) and `sourceBind` `{muTheta→μ, vFace→V, rhoFace→density3}` (`:380-381`).

A builder implementing the **cited** WL map pulls virtual-constraint / Jacobian identities (Jacobian already in map 2) and misses the `μ`/`V`/`ρ` factors that actually sit in `SOURCE_*`. That changes the bridged residual for every source leaf.

### M4. Φ is not a pre-graded comparator leaf; the extraction cite omits `extract_meta`; line 68 drops Φ from the unchanged-dictionary sentence

Question §3(b) requires Φ **coefficientwise per grade, same dictionary as SOURCE**.

What the engines/comparator actually emit:

- PY `FROZEN_PHI` is an ungraded `substitution_map` (`cov:116-124`).
- WL `FROZEN_PHI` is ungraded `"MAP"` (`:975-976`).
- Comparator materializes it via `extract_meta` as `MAP_VARIABLE` + `FIELD_PATH` (`comparator:688-765`, PY `:704-706`, WL `:758-759`), **not** via `(ETA,SIGMA)`.
- Numeric `extract_py_numeric` / `decode_*_key` grade path is `:143-201`, `:631+`. `object_work` sends META through `extract_meta` (`:877-878`).

Directive Object section demands “each retained grade” for every surfaced family, including `N6COV_FROZEN_PHI` (lines 12-21). Retained-order says “operands are already graded by the comparator key” (line 72). Line 31-35 tells the builder to reuse “the `CASE`-line path around `comparator.py:794-825`” — that is pairing of **already-built** leaves, not META walk. Line 68 applies the dictionary unchanged to ACTUAL/PREDICTED/`C_E`/`C_M` and **does not mention Φ**.

A builder can therefore (i) skip Φ, (ii) compare ungraded map images (not coefficientwise), or (iii) re-extract grades from Φ images and then apply map 4’s `η` expansion (M2). Any of those changes channel (b).

Census/control families are correctly split out; this finding is only about `FROZEN_PHI` as a collapse operand.

### M5. Control 1 leaks a shipped witness pattern (value-free violation)

Directive:107: “the shipped run must NOT be all-collapsed by construction.”

That is an expected outcome of the **shipped** witness, not merely “don’t put whole-object equality in the dictionary.” It contradicts the value-free clause (lines 117-121). A builder iterating toward visible targets can withhold a legitimate primitive to keep a `residual_remains` so Control 1 “bites.” That changes the dictionary, hence the residuals.

Control 1’s actual bite is: adding a whole-object equality in a `/tmp` copy forces `collapsed` on every leaf. That remains valid even if the honest shipped run is all-`collapsed`. Drop the parenthetical.

---

## Checks that do **not** fail the gate

**Witness.** Three-valued `collapsed` / `residual_remains` / `undecided`; residual printed first; fresh PIT with own primes/points; bounded; no assert-before-emit; no AGREE/DISAGREE token. Tokens are objects, not verdicts, if the three clauses hold. (Builder must not copy PY primes `(1000000009, 998244353, 1000000033)` at `diagnostic:669` or WL `{1000000009, 998244353, 1004535809}` at `:733` — the “OWN primes” sentence is enough if followed.)

**Controls 2–4, with nits not rising to MUST.** One-sided corruption is the right independence test. Dropped-primitive is conditional on some leaf actually collapsing — unable to fail if nothing collapses, but that is honest given value-free, not a leaked count. Grade-combining tripwire matches the rectangle; M2 is the hole that would still let η back in through a “primitive.”

**BASELINE + census.** `SOURCE_BASELINE` as nominal-control duplicate matches shipped `kappa_a=1, kappa_j=0` (PY `cov:34-35,140-145`; WL `:18-20,845-848`) and the 160/0 `SOURCE_CONTROL_DELTA` tally. Census of the 8 nonzero `ACTUAL_CONTROL_PARAMETERS` + `PHI_DOMAIN_CENSUS` leaves as production-control / domain-coverage, with `n6cov_J_mu`↔`junkMu` and WL-thickness / PY-no-thickness (`cov:137-154` vs WL `:977-980`) surfaced as objects, is right. Not dismissed, not pre-decided.

**Comparator nonzero counts.** Treated as scope (which families/leaves to walk), not collapse targets. Honored except for M5.

**Fence + DoD.** `build→verify→run→report→STOP`, no self-review, no outside edits: explicit. Reusing comparator extraction (not re-joining schema) is the right architecture **if** META goes through `extract_meta` (M4). DoD can detect a missing deliverable and an obvious whole-object dictionary; it cannot catch M1–M3 as currently written.

**Maps 1, 3, 5, 7, 8** are primitive in intent. Map 5 must stay coefficient-name pairing (contraction identity ≠ replacing PY termwise μ by WL’s quotient representative). Map 7’s “tested, not assumed” ban is the correct anti-blanket for Φ.

**No default on-shell:** PY strips DtN/resolvent before source extraction (`diagnostic:383`); WL kernel begins after `sourceBind` (`:380` vs `:385-418`). Exclusion is right.

---

**DIRECTIVE NOT-SOUND** — M1 (degree-2 projector vacuums SOURCE/Φ), M2 (profile expansion reintroduces η into graded leaves), M3 (map 6 cites mass-constraint, not V/source-solve), M4 (Φ ungraded / extract_meta omitted / dropped from the unchanged-dictionary sentence), M5 (Control 1 leaks “shipped is not all-collapsed”).
