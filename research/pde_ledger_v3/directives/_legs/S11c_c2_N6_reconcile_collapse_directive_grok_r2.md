I'll review the directive as an independent decision leg: engines first, then each admitted map, then whether the collapse design would actually compute or claim the right thing.Engines first, then the framing packet. I'll pull the four CAS sources, the shared diagnostic, and the comparator together so the maps are derived from code rather than the directive.The c1 reconcile file isn’t at the cited path; I’ll locate it and then ground each map in the engine sources.The c1 precedent lives under `directives/_measurements`. Next I’ll check the cited engine sites for each of the eight maps.Engines first, then the eight maps against those sites. One map is not an executable primitive: the energy-basis table as specified is a name pairing, but the two constructors split \(W_0\)/\(W_{\mathrm{bg}}\) differently.

## Method

Read both engines and the comparator before treating the directive’s dictionary as given. Citations below are to those sources. This is a design review (the instrument does not exist); no script ablation.

---

## 1. The eight admitted maps

### Map 1 — grad-θ jet duals — **primitive, grounded**

PY geometry atoms are `theta_d{i}` (`S11c_a_interface_geometry_sympy_audit.py:140`). Operator atoms are `grad_theta_{i}` (`S11c_b_brane_operator_sympy_audit.py:227–228`). `wave_jet` remaps `grad_theta_*` → `theta_d*` (`S11c_c2_selfenergy_fold_sympy_audit.py:146–147`). WL `JET_VOCABULARY_BRIDGE` sends `thetaD_i` / `gradTheta_i` → `jet["theta",{i}]` (`.wl:946–948`); `jet[]` itself builds `thetaJet{i}` (`.wl:84–85`). The comparator already aliases `thetaJet{i}` → `theta_d{i}` (`comparator.py:290–311`). Adding the PY dual name is a jet-identity, not a carrier equality.

### Map 2 — source-wave point-preserving map — **primitive, grounded, not whole-object**

WL `sourceMap` keeps `D[f, atom] * atom` on registered wave jets (`.wl:707–711`) after `atPoint[..., "Y"]`. PY `source_value` multiplies the coefficient by `wave_jet(wave, Y)` (`diagnostic.py:392–396`), and `wave_jet` builds applied `s11cc2Field_*(Y,t)` plus spatial/time derivatives (`selfenergy:139–169`). The comparator explicitly refuses field-to-bare collapse (`comparator.py:287–289`, `applied_heads='arguments preserved'` at `:990–993`). Atomwise `jet[f,I] ↔ ∂_I Field_f(Y,t)` is a representation convention for the source evaluation point; it does not identify whole sources.

The domain clause “source jets that appear” is load-bearing: carriers compile wave *symbols* as globals (`diagnostic.py:229–236`, `grades` default `point='X'` at `reconcile.py:95–101`; WL `carrierMap` uses `atPoint[..., "X"]` at `.wl:713–714`). A rewrite that sent every bare jet to `Field(Y)` on carriers would be a different map. The written domain keeps that from being the shipped rule.

### Map 3 — leftover profile-jet names/scales — **primitive; η-re-expansion correctly excluded**

WL `profileRules[]` (`.wl:127–132`) and PY `dx` / `background_jet_expression` (`selfenergy:281–301`; `brane.py:767–778`) both leave `w1ProfileJet{ij}` / `w1_profile_d{ij}` (and the `m1` analogue) with leftover factors `σ_W/L_W` and `μ_R/W_0`. Both engines expand profiles *before* grade extraction (WL `finish` = `retain[(f /. harmonicRules) /. profileRules]` at `:136` then `gradePart` at `:141`; PY `grades` does `xreplace(self.inputs.profiles)` then `shape_coefficients` at `diagnostic.py:247–248`). Re-applying `WBg → W0(1+η·w1)` to a single-grade leaf would reintroduce `η`. Exclusion is required.

### Map 4 — energy-basis coefficient-name table — **MUST: not an executable primitive as written**

The two energy constructors do **not** emit matching leftover *names* for the same contraction.

WL (`constructEnergy`, `.wl:211–226`):

- rewrites `theta^2 → WBg theta^2`, `theta eW → WBg theta eW` (`.wl:217`);
- known rule `WBg theta^2 -> bRho/2` (`.wl:222`);
- so the density term is `(bRho/2) * WBg * theta^2`.

PY (`uniform_coefficient` / `construct_energy`, `brane.py:1584–1594`, `:1763–1771`):

- catalogue entry `(btheta**2, B_rho_3 * W_bg / (2 * W0))` (`:1589`);
- so the density term is `(B_rho_3 * W_bg / (2*W0)) * theta^2`.

After the comparator’s `WBg ↔ W_bg` spelling those terms match iff **`bRho = B_rho_3 / W0`**, not iff `bRho ↔ B_rho_3`. That scale is the S11c-b Bridge A convention (`S11c_b_SHARED_PHYSICS.md:102–108`: \(B_\rho^{(3)}\equiv B_\rho W_0\); the same split is recorded as `BRIDGE_A_RULE = (BRHO, B_RHO_3 / W_0)` in `S11c_b_adjudicated_comparison.py:85–94`). It is a leftover **normalization**, of the same kind map 3 already admits for `σ_W/L_W`.

Two further defects in the written map type:

1. **“Injective coefficient-name pairing”** does not carry that \(W_0\) (nor the WBg-in-monomial vs WBg-in-coefficient split). A name table `bRho ↔ B_rho_3` computes the wrong residual on every θ² channel in μ/source. Other theta-dependent terms are not uniform either: `WBg theta eW -> cCoupling` (`.wl:222`) vs `C * W_bg` (`brane.py:1590`) *is* name-only after spelling — so the map has to be **term-by-term (coeff × monomial) identities including leftover scales**, not a uniform name bijection.
2. **The cited PY site `diagnostic.py:318–349` is constitutive EL of μ**, which is channel (b) — the question. Pairing monomials *in μ* is a disguised whole-μ identification. The only legitimate construction sites are the energy-density constructors (`.wl:211–226`; `brane.py:1584–1594` and `:1727–1872`). The directive already says “not a runtime match against μ/source residuals”; the constitutive citation contradicts that.

This is not the §3 question in disguise *if* the table is built from energy-density contractions with scales. As specified (name table + μ-EL site), astra can either omit \(W_0\) (wrong residual) or match μ monomials (vacuous). That changes what is computed.

### Map 5 — Φ domain/multi-index/derivative spelling — **primitive**

`prolonged_phi` (`covariance.py:63–129`) vs `phiMap` (`.wl:236–243`). Values stay in the residual. Comparator already name-tokens Φ keys (`extract_meta`, `:688–765`, `:757–759`). Leftover derivative-syntax / sorted multi-index (WL `Sort[indices]` at `.wl:84`) is spelling, not whole-Φ.

### Map 6 — leftover density name — **primitive leftover-name, not the c1 freeze**

PY `source_terms` rebinds `rho_br_bg_rho4_constant` → `inputs.density[(rho,)][1]` (`diagnostic.py:382`). That `[1]` is the live 3-density: `BACKGROUND_DENSITY_MAP` stores `Tuple(profile_context, finalize(rho4*W_bg), gradient)` (`S11c_a:1004–1010`); `cases` + `named(...,'VALUE')` then `[1]` is `σ_e = ρ₄ W_bg`. WL `density3 = density4*WBg` (`.wl:839–840`) used as `rhoFace` in `sourceBind` (`:381`). Same constructed 3-density, two leftover names, after both engines have already rebound. Matches c1’s “field-vs-field, not constant↔field” (`S11c_c1_comparator_reconcile.md:137–163`). Not whole-source.

### Map 7 — Jacobian leftover factor — **primitive; degree-2 correctly excluded**

Spelling `1+tr(∇u)` (`reconcile.py:65`) ↔ `1+Sum jet[u_i,{i}]` (`.wl:246`), only where that factor still appears. Degree-2 is a μ_M construction stage (`brane.py:1970–1981`; WL `waveScale[...,2]` at `:246`). Emitted sources are already wave-linear (`wave_terms` one wave, `diagnostic.py:368–374`; WL `waveScale[...,1]` at `:381`). A degree-2 rewrite on those leaves would project them to 0.

### Map 8 — conditional ε / ω — **primitive if occurrence-gated**

Default on-shell / Fourier-of-derivative are c1 kernel identities (PY strips DtN/resolvent at `diagnostic.py:383`; WL `constructKernel` starts at `.wl:385` after `sourceBind` at `:380`). Restricting to atoms that survive maps 1–7 is right. The exact ε identity, if the atom occurs, is the c1 placement `PY·ε = WL` (`S11c_c1_comparator_reconcile.md:100–102`); N6 PY already divides by `eps` in `source_terms` (`diagnostic.py:383,389`). Conditional-omit-if-absent keeps this from becoming a kernel import.

### Excluded list — **right**

Whole carrier/μ/source/Φ; construction-stage degree-2 and profile η-re-expansion; whole velocity / source-solve factor equality (independent builders: PY `reconcile.py:82–92`, `diagnostic.py:378–389`; WL `.wl:274,304,366–381,862–863`; face-velocity correctness still OPEN, `RESOLVED.md:41–43`); `LAB_HELD ↔ MATERIAL_ADVECTED` (anchoring is a retained key, `comparator.py:67–68,91–92`); `σ_W↔η` / `σ_W→0`; default on-shell/Fourier. None of those belong in the dictionary.

---

## 2. Completeness vs over-reach

Needed and missing: the energy-coefficient **scale** in map 4 (above). Map 3 already has leftover scales for profile jets; map 4 needs the same treatment for \(W_0\) / WBg-splitting.

Not missing: face-velocity / source-solve *whole-factor* equality (excluded for the right reason). Jacobian degree-2 (vacuous). Profile η-re-expansion (grade-illegal). Comparator spelling (`verified_spelling_maps`) already covers lowerCamel (`etaBg↔eta_bg`, `W0↔W_0`, …); the physics dictionary should not duplicate that.

Over-broad risk is map 4 if built from constitutive μ rather than energy density. Maps 1–3 and 5–8 are not over-broad as written.

---

## 3. Retained order — **SOUND**

Both engines keep the independent rectangle `(η^i σ_W^j)`, `i,j∈{0,1}`: PY `GRADES = product((0,1), repeat=2)` (`diagnostic.py:43`), extracted independently (`:245–248`); WL `gradeIndices = Tuples[{Range[0,1],Range[0,1]}]` and nested `Series` (`:123–124`). Comparator joins `ETA, SIGMA ∈ {0,1}` as independent axes (`:96–97`). Directive forbids grade-combining, `η·σ_W→η²`, `σ_W→0`, `σ_W↔η`, and η-reintroduction into a graded leaf. Φ is ungraded metadata (`covariance.py:116–124`; WL `:975–976`; `extract_meta` keys `MAP_VARIABLE`/`FIELD_PATH` with no grade axes, `comparator.py:688–765`); grading via extract_meta → profile-expand → four coefficients **before** the dictionary is the numeric path, not the forbidden post-grade re-expansion.

---

## 4. Collapse witness — **SOUND**

Three-valued `{collapsed, residual_remains, undecided}`; residual printed first; no pass/fail interpretation. Bounded `cancel`/`together`/`factor`, then a fresh PIT with primes/seeds disjoint from PY `(1000000009, 998244353, 1000000033)` (`diagnostic.py:669`) and WL `(1000000009, 998244353, 1004535809)` (`.wl:733`); denominator rejection; degree/exclusion bound; max attempts/time/RSS; `undecided` on exhaustion; provenance emitted. No designed path that asserts or emits a verdict token. Clause 1’s “symbolic zero-witness” covers the enum.

---

## 5. Controls — **SOUND as redesigned (not tautological)**

| Control | Can it fail? |
|---|---|
| Per-map defining-relation → **control-delta** (baseline bridged residual − corrupted) | Yes, on leaves that use the map. “Corrupted residual nonzero” would be vacuous on already-nonzero leaves; control-delta is the able-to-fail form. |
| Production-reachability census | Honest emit of inert maps; not a bite-test, and should not be forced. |
| Blanket whole-object equality → all `collapsed` | Yes, if the witness cannot report `collapsed`. Does **not** by itself prove the shipped dictionary isn’t whole-object when the honest run is all-collapsed; control-delta plus “no whole-object map reachable” in `test_…` cover that. The “do not expect shipped all-collapsed” sentence is right. |
| Grade-combining tripwire | Yes if a map reintroduces `η` into a **graded** leaf or binds `σ_W`. Φ’s pre-grade profile expansion must not trip it (the text already says “graded leaf”). |

“Dropped primitive leaves a residual” on already-nonzero leaves is the v1 tautology. For maps that occur, dropping them is a defining-relation corruption (control 1). Completeness of the *set* of maps is this review’s job, not a fourth tautological ablation.

---

## 6. BASELINE + census — **SOUND**

Shipped `ACTUAL_A_RHO=1`, `ACTUAL_JUNK=0` (`covariance.py:34–35`); `actual_amplitudes` vs baseline `(1,0)` (`:140–154`). WL `actualAdvectionCoefficient=1`, `actualJunkCoefficient=0` (`:18–20`); `materialAmplitude(..., actual, ...)` vs `baseline = materialAmplitude(..., 1, 1, 0)` (`:845–848`). At those settings μ_actual = μ_baseline, so `SOURCE_BASELINE` is a nominal-control duplicate, not a third discriminator. The `160/0 SOURCE_CONTROL_DELTA` line is a **comparator** fact used only to justify that label, not a collapse target.

Census: `ACTUAL_CONTROL_PARAMETERS` selects material amplitude + junk (`covariance.py:146–153`; WL `:977–980`); `PHI_DOMAIN_CENSUS` abort-on-uncovered (`covariance.py:105–115,125`; WL `:850–852`). Comparator adds no `kappa_a↔ADVECTION` aliases (`comparator.py:683–685`). Crosswalk + retain unmatched (WL `THICKNESS`; PY `actual_amplitudes` retags only theta-advection + junk, `:137–154`) is production-control / domain-coverage, not metadata and not a pre-decided collapse. The “8 leaves” are the committed matched nonzero census leaves (`run_tally.txt`: 4+4); unmatched stay surfaced.

---

## 7. Value-free / leak — **SOUND**

No expected collapse outcome, count, or pass condition for the **instrument**. Prior-run nonzero counts are not reproduced as targets. `160/0` is confined to labeling BASELINE from the comparator. Builder report asks for a **neutral** witness tally, not an interpretation.

---

## 8. Fence, extraction reuse, DoD — **SOUND** (once map 4 is fixed)

`compare_family` materializes, emits `CASE`, then `release_case` (`comparator.py:794–851`); `object_work` returns only accounting (`:862–898`). Reusing `load_py_jsonl` / `load_wl` / `verified_spelling_maps` / leaf keying / `materialize` / `extract_meta` without editing the comparator is the right API. Fence is explicit (build→verify→run→report→STOP; no self-review; no outside edits). DoD can detect a missing deliverable, a comparator edit, a reachable whole-object/construction map, ungraded Φ, asserts, and a witness that cannot report `collapsed` (blanket ablation). It cannot detect a name-only map 4 that omits \(W_0\) — that is the dictionary defect above.

---

## MUST finding

**M1 — Map 4 is not a mechanical name/structure correspondence as specified.** The energy-density constructors split the θ² term as `(bRho/2) WBg θ²` (`.wl:217,222`) vs `(B_rho_3 W_bg)/(2 W0) θ²` (`brane.py:1589`). The primitive is the contraction-term identity including leftover scale \(B_\rho^{(3)} = bRho\cdot W_0\), not an injective leftover-name table. Citing constitutive EL (`diagnostic.py:318–349`) as a pairing site points at μ (channel b). Until map 4 is specified as termwise (coeff × monomial) identities built only from `.wl:211–226` and `brane.py:1584–1594,1727–1872`, with leftover \(W_0\)/WBg factors first-class (as map 3 already does for `σ_W/L_W`), the instrument will either compute the wrong residual or equate μ.

---

**DIRECTIVE NOT-SOUND** — M1 (map 4 name-table / μ-EL site vs contraction-term identity with \(W_0\)).
