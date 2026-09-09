I'll form an independent view from the engines' actual outputs and the governing spec first, then check the directive against that. Gathering the sources and the N6 tag/schema evidence.I'll start from the engines, SHARED_PHYSICS, and tag schemas, then check the directive against that view.Next I’ll inventory actual emitted tags, PIT primes/sample points, and the `_NODES`/ARITHMETIC schemas from the engine outputs.I’ll inventory the actual JSON object names, PIT primes, and `_NODES` schemas from the engine transcripts with mechanical lookups only.Independent view of the engines is in place. Next I’ll read the directive and check it against that evidence.**Verdict: NOT-SOUND.** Sealing PIT is the right call, and the three-valued / no-verdict / no-target frame is mostly intact. The join set, inner-axis bridge, and `_NODES` reconstruction are not specified well enough for a blind builder to produce a correct comparator.

Empirically, the two engines cannot share a numeric PIT channel. WL `LOCAL_PROBE` primes are `{1000000009, 998244353, 1004535809}` with 17 cells (16 face×kinematic plus `{0, EVANESCENT³}`) and 408 samples (`3×17×8`). SymPy primes are `{1000000009, 998244353, 1000000033}` with 16 cells `[±1,0|1,0|1,0|1]`, different seeds (`110602` / `110603`), and 432–576 samples. Digests are WL integer fingerprints vs SymPy sha256. Symbolic operands do exist on both sides (WL `ARITHMETIC` `Inactive[Plus/Times/Power]`; SymPy `_NODES` `{op,args,degree}`). So channels (A) structural and (B) symbolic are the real ones; differencing PIT residues would be a false measurement.

---

### F1 — JOIN SET drops shared covariance families (missed channel)

**Decision.** Covariance join is only `R_COV`, `R_COV_BASELINE`, `R_COV_INCREMENT`, `SOURCE_ACTUAL`, `SOURCE_PREDICTED`, `SOURCE_BASELINE`, `FROZEN_PHI` (directive lines 84–86).

**Evidence.** Both engines emit the control-delta pair:

- WL `mathematica/S11c_c2_N6_mathematica_audit.wl:30-33` `covNames` includes `SOURCE_CONTROL_DELTA`, `R_COV_CONTROL_DELTA`, `PHI_DOMAIN_CENSUS`, `ACTUAL_CONTROL_PARAMETERS`.
- WL tags 25–26, 29–30 in `mathematica/out/S11c_c2_N6_mathematica_audit.out`.
- SymPy `S11CC2_N6COV_SOURCE_CONTROL_DELTA`, `_NODES`, `S11CC2_N6COV_R_COV_CONTROL_DELTA`, `_NODES`, `PHI_DOMAIN_CENSUS`, `ACTUAL_CONTROL_PARAMETERS` in `scripts/out/S11c_c2_N6_covariance_sympy.out`.

Those are the covariance instrument’s one-sided control movement, jet-domain census, and live knobs — not bookkeeping.

**Why wrong.** A builder handed only this join never measures whether the two engines agree on the control-delta of the commuting square. That is a missed cross-engine channel. `PROVENANCE` is correctly excluded (digest seal).

**Minimal fix.** Add `SOURCE_CONTROL_DELTA` and `R_COV_CONTROL_DELTA` to the covariance join (symbolic + support, PIT sealed). Add `PHI_DOMAIN_CENSUS` and `ACTUAL_CONTROL_PARAMETERS` as structural joins (seal digest subfields). Keep `PROVENANCE` surfaced-not-joined.

---

### F2 — Inner-axis “bridge” is under-specified; support comparison will be guessed or counted

**Decision.** Lines 57–62, 123–126: column bases differ; the column↔slot bridge is “a computed residual adjudicated after the run, not a pre-registered collapse”; if ambiguous, emit `axis_set_mismatch`.

**Evidence.** Keys are not a permutation of one vocabulary:

| axis | SymPy | WL |
|---|---|---|
| increment block | `BLOCKS = product((TRANSVERSE_TO_THICKNESS, THICKNESS_TO_TRANSVERSE), (THETA, E_W, DIV_U))` (`diagnostic_sympy.py:45-46`) | `pair[[1]]<>"_FROM_"<>pair[[2]]` with `{THETA,E_W,LONGITUDINAL,TRANSVERSE}` (`audit.wl:483-484,719-721`) |
| kernel | integer `0/6/9/12` | `LOCAL_BARE` / `FOURIER_KOUT_Y` / `FOURIER_KOUT_KIN_Y` / `FOURIER_KOUT_KIN_Y_MIDDLE` (`audit.wl:485`) |
| grade | `[η,σ]` 2-tuple | `{1,η,σ}` (`LOCAL_SETUP`; `audit.wl:713`) |
| face-sum | integer `0` (`diagnostic_sympy.py:512`) | `"SUM"` (`audit.wl:697-699`) |
| SOURCE | `(face, grade)` only (`reconcile_sympy.py:105-111`) | `(face, wave, grade)` (`audit.wl:707-711`) |
| SLOT_GUARD route | JSON `probe` label, not a column (`diagnostic.out` `probe: EULERIAN`; columns `['E_W','linear',[0,0]]`) | inside the slot key `{"EULERIAN", 1, 1, {1,0,0}}` (WL tag 32) |
| components | one DAG per column | `COMPONENT_AXES -> {1,2,3}` on transverse families (WL tags 10–11) |

**Why wrong.** T7’s “don’t collapse a disagreement” is right. It does not replace a typed inner-axis decoder. c1 required `decode_py_key`/`decode_wl_key` from the constructors (`S11c_c1_comparator_build_directive.md:86-90`). Here a builder will either (i) string-match and get a silent 0-join, (ii) compare nonzero *counts* (false structural agreement), or (iii) invent `THETA_FROM_TRANSVERSE ↔ (THICKNESS_TO_TRANSVERSE, THETA)` and maybe swap test/trial. The support channel the directive names cannot be computed without that decoder.

**Minimal fix.** Require typed inner axes from the two constructors, as distinct axes, with an emitted pairing table `{matched, sympy_only_column, wl_only_slot}`. Restrict closed spelling folds to name/CAS identities (grade `(1,η,σ)↔(η,σ)`, face `SUM↔0`). Do not pre-register block or kernel-family equality. Compare support/symbolic only on matched keys; unmatched keys are the residual, not a count.

---

### F3 — `_NODES` is not a self-contained DAG (missing trap → false-zero symbolic channel)

**Decision.** Lines 82–83, 116–118: pair `<TAG>_NODES` to WL `ARITHMETIC`; reconstruct from `op/args/degree/literal`.

**Evidence.** `reconcile_sympy.py:151-181` emits `ARITHMETIC_DAG = {nodes: [...]}` and then `_NODES = {columns, root_ids, root_nodes: [nodes[index] for ...]}`. Non-leaf args are `{'ref': i}` into **that shared `nodes` list**, not into `root_nodes`. Sample `S11CC2_N6RC_R_N6_NODES`: 228 roots are `op=number` with `literal 0`; 60 are `op=add` with `ref`s. `S11CC2_N6RC_ARITHMETIC_DAG` has 8477 nodes (`number/variable/add/mul/pow`). `ARITHMETIC_DAG` is not in the join set and is not named as the ref table. Diagnostic guards have **no** `_NODES` at all.

**Why wrong.** A reader of `_NODES` alone cannot reconstruct a non-leaf. The likely failure is comparing only the zero `number` roots and reporting symbolic agreement on the empty circuit.

**Minimal fix.** State that reconstruction is `ARITHMETIC_DAG.nodes[root_id]` with `ref` indices into that table; `_NODES.root_nodes` is a copy of the root record only. Add a trap: guards have PIT tables and no `_NODES` (symbolic is WL-only). Add a synthetic test that a `ref` DAG round-trips and that dropping `ARITHMETIC_DAG` is `parse_failed`, not a zero residual.

---

### F4 — `load_py_jsonl` specifies a `CASE` field the streams do not have

**Decision.** Line 113: new loader is `object`/`CASE`/`data`. Lines 96–97 (correct): case is top-level `anchoring`+`density`, **not** a `CASE` array.

**Evidence.** Every SymPy line is `{"object","anchoring","density","data","dimension","N12",...}` (e.g. `S11CC2_N6RC_R_N6` keys). There is no `CASE` field.

**Why wrong.** A builder following the mechanical-base signature keys nothing and silently 0-joins every family.

**Minimal fix.** Delete `CASE` from the loader contract. Key on `anchoring`×`density`. For diagnostic guards, also type the `probe` label as an inner axis (F2).

---

### F5 — “SEAL” is overloaded; carrier symbolic comparison can be dropped

**Decision.** Lines 44–51: PIT/digests are sealed = surface, **never difference**. Lines 137–138, 166: the blind-carrier representation is “SEALED (seal 2)” while seal 2 is the digest seal.

**Evidence.** `CARRIER_EULERIAN` has WL `ARITHMETIC` (tag 3) and SymPy `_NODES` (`S11CC2_N6RC_CARRIER_EULERIAN_NODES`). Lines 56–68 still list symbolic Inactive↔DAG as channel (B).

**Why wrong.** In c1, “do not reconcile” means don’t *fold* a difference; still print A−B. Here “SEAL” for PIT means don’t subtract. Calling carrier-rep “seal 2” lets a builder skip the one symbolic channel that exists for `C_E`. That is a missed channel, not T7 discipline.

**Minimal fix.** Reserve SEAL for PIT residues, sample points, and digests. For carrier-rep: surface both, compute the three-valued residual, do not fold it to identity.

---

### F6 — Expected-value leak on a sympy-only family

**Decision.** Line 71: ``REP_INVARIANCE_RESIDUAL` (⚠ the §5c-**superseded** wrong-object — **nonzero is EXPECTED**, ⛔ never read as cross-engine disagreement)``.

**Why wrong.** Rule 5 / N8: no residual target, no “should be zero/covariant.” “Nonzero is EXPECTED” is a target. It is also unnecessary: the family has no WL sibling and is already excluded. A builder can bake a nonzero-is-OK guard into accounting.

**Minimal fix.** Keep it as `sympy_only` with reason “no WL sibling.” Drop the expected-nonzero / wrong-object sentence. Enumerate the exclusion as: every diagnostic `S11CC2_N6_*` except the six guards, plus `REP_INVARIANCE_*` and `CONTROL_INDEPENDENCE_*` (the listed three residuals + `PREMISES` are incomplete).

---

### F7 — `load_wl` reuse will not parse this stream (silent 0 WL extract)

**Decision.** Lines 111–112: reuse c1 `load_wl` (multiline `<|…|>` through the next tag line).

**Evidence.** c1 `S11c_c1_cross_engine_comparator.py:54-55,189-212`: regex `^(WL_S11CC1_[A-Z0-9_]+):`, split on `:`. N6 (`audit.wl:48-53` and the `.out`): `WL_S11CC2_<NAME> = <|…|>` , 40 **single-line** tags. Inputs line 92 does state the `=` form; the reuse bullet does not.

**Why wrong.** Unchanged `load_wl` yields an empty WL map. Every family becomes `sympy_only` and the comparator never measures.

**Minimal fix.** Require a new N6 loader: prefix `WL_S11CC2_`, delimiter ` = `, one tag per line, lazy line-addressing. Reuse `parse_wl_value` / `canonical_basic` / `residual` / `materialize`/`release_case` only.

---

### F8 — Per-family `[L,T,M]` on SymPy is not inside `data` or `N6RC_DIMENSIONS`

**Decision.** Lines 57–58, 156–159: dimensions come from `DIMENSIONS`; the named trap is only the `N6RC_DIMENSIONS` container mismatch.

**Evidence.** Increment dimensions are the top-level JSON `dimension` array parallel to `data.columns` (`emit(..., dimension=[schemas...])` in `diagnostic_sympy.py:769-772`; `R_N6` top-level keys include `dimension`, not inside `data`). `N6RC_DIMENSIONS` is only `{eulerian, material}` keyed by `(row, slot)` from `build_increment`’s carrier-coeff record (`diagnostic_sympy.py:483-513`, `reconcile_sympy.py:274`). WL per-slot `DIMENSIONS` lives *inside each family* (`audit.wl:804`) and `N6RC_DIMENSIONS` dumps all objects (`audit.wl:987-997`).

**Why wrong.** Comparing only the `DIMENSIONS` family misses the increment/source/guard `[L,T,M]` channel. Looking in `data` for it is a silent 0-extract.

**Minimal fix.** Trap: SymPy per-column units are the top-level `dimension` list; WL’s are the per-slot `DIMENSIONS` field; `N6RC_DIMENSIONS` is a third, differently keyed object.

---

## Checks that pass

**PIT seal.** Correct. Different third primes, cell sets, draw counts, and sample generators. Re-sampling would be a third engine, not a join of emitted objects. Nonzero-support as structure is the legitimate numeric remnant; residue values are not.

**Three-valued / no-verdict.** Lines 31–40 preserve `residual()` returns, ban `PASS`/`FAIL`/`AGREE`/`COVARIANT`, exit 0 on disagreement. Aside from F6, no residual target, no prime literals, no “should be zero.”

**Authority.** Pointing at §N8/T7 (~477), §5c (~303), §1b (~76) is the right trio. N8’s generic contract (join by name, three-valued, no boolean, no blanket collapse, name-don’t-adjudicate what this box cannot close) applies to N6. N8’s *load-bearing* residual is the §3c self-energy increment — this instrument correctly does not hunt `S11CC2_SELF_ENERGY_INCREMENT`. PIT blindness is not in N8; the extra seal clause is required, not a spec-currency error. “Giants deferred” as per-object `deferred_oversize` on a 30 GB box is a faithful reading of “names, does not pre-adjudicate, whatever it cannot close on this box” (`SHARED_PHYSICS.md:477-486`, `418-437`), not the ≥64 GB S11c-b bulk residual. Do not copy that 64 GB figure onto N6 tags (largest WL line is `CLOSURE_GUARD_RESIDUAL` ~80 MB).

**Reuse + `_NODES` reader.** Accurate that c1 provides `load_wl` (format-incompatible, F7), `parse_wl_value`/`canonical_basic` as Inactive→SymPy, three-valued `residual` from `S11b_cross_engine_comparator.py:682`, and lazy `materialize`/`release_case` (`c1:532-550`). A new `_NODES` reader is required: emitted SymPy is JSON, not srepr (`srepr` is only an internal sha, `diagnostic_sympy.py:92`). DAG ops are `number|variable|add|mul|pow` — reconstructable **if** F3 is fixed.

**Named container traps.** `DIMENSIONS` shape, `ARITHMETIC`+PIT dual form, `FROZEN_PHI` Wolfram vs expr-string, dedicated vs embedded `PRIMES`/`PIT_PROVENANCE` are real. WL `N6_` is thin: `LOCAL_*` + six guards only (`N6RC_DIMENSIONS` is `N6RC_`, already in the reconcile join). Still missing: F2 inner axes, F3 `ARITHMETIC_DAG` refs, F8 top-level `dimension`, reconcile’s `S11CC2_N6RC_N6_PRIMES` double prefix vs covariance `N6COV_PRIMES`.

**DoD.** Can fail a silent 0-join (`extracted-leaf count`, `axis_set_mismatch`) and a PIT-seal violation (synthetic test). It cannot fail a manufactured *column* pairing or an `_NODES` false-zero unless F2–F3 tests are added. `join>0` is correctly not an exit criterion.

---

**Ranked:** F1 missed covariance control-delta join → F2 inner-axis under-spec → F3 `_NODES`/DAG refs → F4 `CASE` vs `anchoring`/`density` → F5 seal overload on carriers → F7 `load_wl` grammar → F8 dimension location → F6 expected-nonzero leak.

Fix F1–F5 before a builder launch. F6–F8 are cheaper and still load-bearing.
