I'll start from the source engines and committed evidence, then judge the reconcile framing against those—not the other way around.I'll read the engine sources first and form a view from the constructions, then look at the framing.Next I’ll pull the engine constructions that actually define those objects, then the committed run evidence.I have a source-grounded picture of the residuals and order. Next I’ll read the framing and check it against that, not the reverse.The open gap is the **operands**, not the within-engine differences. §3 asks that question. §4 would let the collapse test treat the question as an identity, which is the change that has to be made before any instrument is written.

---

## Independent view (engines + tally, before the framing)

**Within each engine the differences vanish; the operands need not.**

SymPy computes `CARRIER_BRIDGE_RESIDUAL = residual(C_E, C_M)` at `reconcile_sympy.py:254` and `R_COV = actual − predicted` at `covariance_sympy.py:194`. WL computes the same differences at `mathematica_audit.wl:874` (`bridge = mapCombine[eMap, mMap, sub]`) and `:878` (`covDelta = mapCombine[sMMap, pMap, sub]`). Vanishing is not on those lines; it is the committed PIT/tally:

- Per-engine SymPy: `CARRIER_BRIDGE_RESIDUAL` 0/80 in all 4 cases (`S11c_c2_N6_reconcile_adjudication.md`); `R_COV` 0 in all 4 (`S11c_c2_N6_covariance_build_clearance.md`).
- Cross-engine matched keys (`run_tally.txt`): `N6RC_CARRIER_BRIDGE_RESIDUAL` 320 ZERO / 0 NONZERO; `N6COV_R_COV` 160 ZERO / 0 NONZERO.

Because SymPy’s bridge is already 0, a 320-zero cross-engine bridge residual means WL’s bridge is 0 on those keys too. Same for `R_cov`. That is `0 − 0 = 0`. It does **not** imply `C_E^WL = C_E^PY` or `ACTUAL^WL = ACTUAL^PY`.

The tally on the **operands** is the actual gap:

| family | ZERO | NONZERO |
|---|---:|---:|
| `N6RC_CARRIER_EULERIAN` / `MATERIAL` | 280 | **40** each |
| `N6COV_SOURCE_{ACTUAL,PREDICTED,BASELINE}` | 16 | **76** each |
| `N6COV_FROZEN_PHI` | 42 | **18** |
| `N6COV_ACTUAL_CONTROL_PARAMETERS` | 0 | **4** |
| `N6COV_PHI_DOMAIN_CENSUS` | 0 | **4** (+44 BOOL) |

On the same 40 matched carrier keys, both engines have `C_E = C_M` internally and both differ from the other engine — `CARRIER_MATERIAL` is also 280/40. Same cancellation pattern is available for the 76 source keys: `R_cov = 0` in each engine is compatible with `ACTUAL^WL − ACTUAL^PY = PREDICTED^WL − PREDICTED^PY ≠ 0`.

`R_N6`, channels, operands, RC `SOURCE_{EULERIAN,MATERIAL,BRIDGE}` are entirely unmatched (`run_tally.txt`: 0/0/8640 and 0/0/640). PY RC sources are `(face, grade)` (`reconcile_sympy.py:105–111`); WL `sourceMap` is `{face, wave, grade}` (`mathematica_audit.wl:708–711`). That is a schema non-join, not a measured residual. No cross-engine comparison of `R_N6` exists.

Retention is the independent rectangle `(η,σ_W) ∈ {0,1}²`: `GRADES = product((0,1), repeat=2)` (`diagnostic_sympy.py:43`), `Compiler.grades` keeps `g in GRADES` (`:248`), emit `declared_retained_rectangle=[[1,*g] for g in GRADES]` (`:772`). WL: `gradeIndices = Tuples[{Range[0,1], Range[0,1]}]` and `Series[..., {etaBg,0,1}, {sigmaW,0,1}]` (`mathematica_audit.wl:123–124`). Comparator joins `ETA,SIGMA ∈ {0,1}` as independent axes (`comparator.py:96–97`). Per-engine N6 residual lives in retained `σ_W^1` (`reconcile_adjudication.md`).

---

## 1. Right question, not a proxy — §2 structural argument

**The structural argument is correct.** Quote both sides:

Framing §2: “A difference of two objects that are each 0 in each engine is 0 cross-engine **trivially** — it does NOT establish that the operands agree.”

Engines + tally: as above. Line-level construction:

- SymPy `C_E`: `pressure_coefficients(expanded_rows(inputs.slab[α,ρ]))` (`reconcile_sympy.py:196–197`) — imported slab.
- SymPy `C_M`: `pressure_coefficients(face_factory(…MATERIAL…))` (`:88–92, 205`).
- WL `C_E`/`C_M`: `carrier[eulerianSlabFace[…]]` / `carrier[materialFaceFold[…]]` with `carrier[rows_] := Table[D[rows[[r]], p] /. slotZero, …]` (`mathematica_audit.wl:327–336, 866–868`).
- SymPy sources: `source_terms(..., mu_pred, V_E)` vs `source_terms(..., mu_actual, V_M)` (`covariance_sympy.py:191–193`); `R_cov = source_difference(actual, predicted)` (`:194`).
- WL: `SOURCE_ACTUAL = sMMap`, `SOURCE_PREDICTED = pMap`, `R_COV = sM − p` (`:871–878, 882–883`).

Cite nit (does not change the question): vanishing is the tally / per-engine PIT, not `reconcile:254` or `cov:204` (`:204` is `R_COV_CONTROL_DELTA`, not `R_cov`). The argument still holds.

§3(a)–(c) target exactly those operand residuals. That is the actual open gap. `SOURCE_CONTROL_DELTA` 160/0 in §1.i is even more trivial: shipped `ACTUAL_A_RHO = 1`, `ACTUAL_JUNK = 0` (`covariance_sympy.py:34–35, 140–145`) makes actual = baseline by construction.

---

## 2. Channel (b) — two-errors-cancel / single identity

**The cancellation risk is the right discriminator; “a single identity” is slightly the wrong object.**

Right test: **one uniform, independently justified identification chain**, applied to `{Φ, ACTUAL, PREDICTED, BASELINE}` together, at the retained rectangle. That is the c1 staged bridge (`S11c_c1_comparator_reconcile.md` §2), not a different map per operand.

- Weaker (wrong): collapse each operand with its own map. That is exactly `[[feedback_matching_number_is_not_evidence]]` — two representational stories that keep `R_cov = 0`.
- Stronger and better: include `Φ` in the same chain. `PREDICTED = μ_E.subs(Φ)` (`covariance_sympy.py:132–134, 185–191`; WL `:849, 853, 871`). Collapsing sources without collapsing `Φ` (or vice versa) reintroduces cancellation between the map and its image. Framing splits (b) and (c); the instrument must not.

Framing §3(b) consequent “If yes ⇒ the two engines' constitutive sources genuinely AGREE” is too strong in the c1 sense: collapse-to-zero is AGREE **on matched keys under that chain**. Unmatched RC sources / `R_N6` stay UNDECIDED. c1 over-claimed closure-wide AGREE and had to walk it back.

---

## 3. Retained order (§5)

**The asked order is right. A `σ_W→0` (or grade-mixing) reduction would repeat L-CAS.**

§5 correctly forbids `σ_W→0`. The engines actually retain the independent rectangle (cites above). `frozen_relations` also emits `'retained_grades': n.GRADES` (`reconcile_sympy.py:77`) and FROZEN_PHI `'retained_grades': n.GRADES` (`covariance_sympy.py:123`). The metadata string `'retained_eta_sigma_rectangle'` at `reconcile:74` is a label; the load-bearing object is `GRADES` / `gradeIndices` / `declared_retained_rectangle`.

**Finding (changes the test):** §4’s bullet “the **σ_W binding** and the retained-order rectangle” conflates a domain with a substitution. c1’s kinematic `σ_W = η_bg·W_0/L_W` **identifies the two grade axes**. N6 emits them independently (`ETA`/`SIGMA` each in `{0,1}`). The per-engine residual is in retained `σ_W^1`. Binding `σ_W` to `η` leaves that rectangle; it is a cousin of the L-CAS `σ_W→0` proxy, not a representational identity. Retained order is **where the test runs**, not a map to apply.

---

## 4. Candidate identities (§4)

**§3 asks the right physics question. §4 names that question as if it were an identity — that can make a collapse vacuous.**

c1 listed *specific* independently justified maps (`q_out ↔ qOut`, profile FT-of-derivative, `ε` placement). It did **not** list “the kernel identity.”

§4’s second and last bullets:

- “the **blind-graph-Eulerian-slab-face ↔ imported-slab identity**”
- “the **constitutive-source blind-derivation ↔ imported-μ identity**”

Those are §3(a) and §3(b) **questions**, not maps. Applying “the two slabs are the same representation” is a blanket collapse (`[[feedback_handcode_comparison_never_blanket_collapse]]` / rule 5). The comparator already refused to pre-register block/kernel/component equality (`S11c_c2_N6_comparator_directive_review.md` finding B). Putting the desired equality in the identity list undoes that.

**Wrong axis:** “S11c-a anchoring/coordinate map relating LAB_HELD and MATERIAL_ADVECTED charts” maps **between cases**. The comparator already joins at fixed `(ANCHORING, DENSITY)`. A cross-case chart map is not a WL−SymPy bridge.

**Missing, and actually specified in both engines:**

- `jet_vocabulary_bridge` / `jet_physical_field_images`: `(a.grad_theta[i], b.grad_theta[i])` with common `wave_jet` (`reconcile_sympy.py:59–69`); WL `JET_VOCABULARY_BRIDGE` (`mathematica_audit.wl:946–948`).
- `quadratic_density_jacobian = 1 + tr(∇u)` and `density_wave_projection_degree = 2` (`reconcile_sympy.py:65–66`); WL `materialAmplitude` uses `(1 + Sum[jet[u_i,{i}]]) (density /. MAP)` at waveScale degree 2 (`:246`).
- Φ as constructed, to be *tested* not assumed: PY `prolonged_phi` along `DERIVATIVE_MAP` (`covariance_sympy.py:63–129`); WL `phiMap` via `tdMany` on `theta+a` / `eW+h` (`:236–243, 849`).
- c1 maps **only where those symbols actually occur** in these N6 operands (`ε`, `omega` assumption, Fourier-of-derivative). Do not import the whole c1 table unread.

On-shell dispersion is a c1 kernel identity. These objects are graded operator coefficients at the formal rectangle, not the DtN kernel. Use on-shell only if a matched residual still carries an explicit `q² − (ω²/c² − |k|²)` factor after the name maps — do not fold it in as a default.

---

## 5. Pre-adjudication / leak (§6)

§6 and the header do **not** pre-decide that a channel reconciles or give an expected collapse value. Conditionals in §3 are questions.

Leak-ish, not MUST: §1.i titles the trivial zeros “DIRECT cross-engine AGREEMENT” and inherits `run_data.md`’s “`R_cov = 0` is CROSS-ENGINE CONFIRMED” / “geometric carrier reconciles cross-engine (`C_E=C_M`)”. §2 then correctly retracts that. Relabel §1.i as “matched residual-of-residuals is 0 (trivial given within-engine vanishing)” so §1 cannot be quoted as the disposition.

`R_N6(=18/288)` and `R_cov(=0)` in §7.6 are per-engine SymPy facts (`RESOLVED.md`); say so. Caveats stay open — consistent with `RESOLVED.md`.

---

## 6. Scope fork (§8) — **A**

A faithful retained-order collapse on the **matched nonzero leaves** is in the same resource class as the comparator that already finished.

Evidence (`run_data.md`, `run_tally.txt`):

- `peak_rss_kib: 330412`, `runtime_s: 10294`, `deferred_oversize: 0`, largest operand ~80 MB.
- Those 334 NONZERO leaves already survived `residual(..., leaf_budget_seconds=2)` and returned a computed nonzero, not `ResidualFailure` / deferral.
- N6 carries none of the c1 ≥64 GB families (directive review H).
- Unmatched `R_N6` / channels **cannot** be collapsed without the forbidden block/kernel fold. They stay UNDECIDED under both A and B (c1 UNDECIDED-surfacing, not a heaviness deferral).

Path B would be right if the test required expanding the full ~80 MB objects rather than the already-materialized matched residuals. The comparator already materialized `operand_A` / `operand_B` per matched key (`comparator.py:816–825`). Default A for tractable channels is correct.

**(d) is part of A:** confirming 8 metadata leaves is a fact-lookup-sized walk, not a heavy CAS.

---

## 7. Census / minor (§3d)

**Correctly provisional — not assumed non-load-bearing.** Confirm, or route into (a)–(c).

`ACTUAL_CONTROL_PARAMETERS` (`covariance_sympy.py:146–153`): `kappa_a`, `kappa_j`, `baseline_parameters`, `material_tag`, `junk_symbol` (`n6cov_J_mu`), junk dimension/wave/amplitude, sampler. WL (`:977–980`): `ADVECTION`, `THICKNESS`, `JUNK`, `JUNK_CASE`, `JUNK_SYMBOL` (`junkMu`), `MATERIAL_NORMAL`. Comparator explicitly does **not** alias `kappa_a = ADVECTION` (`comparator.py:684–685`; review record). Tally: 0 ZERO / **4 NONZERO** / 80 UNDEC — some keys **matched and differed**. Likely `junk_symbol` (`n6cov_J_mu` vs `junkMu`) = bookkeeping. WL also emits a **thickness** actual control PY does not (`actual_amplitudes` retags only theta-advection + junk, `:137–154`). That is not obviously metadata; confirm.

`PHI_DOMAIN_CENSUS` is PY `:105–115` (the framing’s `105–124` bleeds into `FROZEN_PHI` at `:116`). It is a computed domain: coverage tuples, uncovered, max rank, before substitution. Clearance: `uncovered=[]`, covered 11/11, max rank 2 — load-bearing for Φ construction. 4 NONZERO + 44 BOOL: the bools are the required `BooleanNotResidualable` rejection; the 4 nonzero could be a matched coverage/rank leaf. If coverage of imported μ jets differs, that is channel (c), not bookkeeping.

---

## 8. Architecture (§7)

**Rule-compliant.** Question and G4 adjudication stay with the orchestrator; collapse instrument is Codex-written + G1 (fresh Claude + Grok); its directive gets G2 decision legs; disposition gets legs (c1 correction-verify pattern); mechanical fact-lookup = regenerate comparator output and retrieve residual expressions verbatim; the zero-test-after-substitution is the instrument. Matches E1 / G1 / G4 / L-CAS.

---

## Re-framing required (§4 + the simultaneous set)

Keep §2–§3, §5–§8 (with the §1.i relabel and §3(b) “matched keys” scope). Replace §4 with independently specified maps, and couple (b)+(c):

**Test (uniform chain, independently justified, not tuned to zero):**

1. Jet vocabulary already emitted: `a.grad_theta[i] ↔ b.grad_theta[i]` via common `wave_jet` (`reconcile_sympy.py:59–69`; WL `:946–948`).
2. Density Jacobian / wave-projection degree as emitted (`reconcile_sympy.py:65–66`; WL `:246`).
3. c1 name/assumption maps **only** where those atoms appear in these residuals (`ε`, `omega` realness, Fourier-of-derivative).
4. Φ spelling as both engines construct it (`prolonged_phi` / `phiMap`) — **test** equality of the maps; do not assume it.

Apply that **same** chain to `{FROZEN_PHI, SOURCE_ACTUAL, SOURCE_PREDICTED, SOURCE_BASELINE}` and, separately, to `{CARRIER_EULERIAN, CARRIER_MATERIAL}`. Print collapsed residuals at the independent `(η^{≤1}, σ_W^{≤1})` rectangle. Decide nothing.

**Do not treat as identities:** blind-slab ↔ imported-slab; blind-source ↔ imported-μ; LAB_HELD ↔ MATERIAL_ADVECTED; `σ_W → η W_0/L_W`.

Whatever does not collapse, and every unmatched family (`R_N6`, channels, RC sources), stay UNDECIDED — not disagreement, not closed.

---

**QUESTION NEEDS-WORK** — §3 is the right operand question at retained order, but §4 would authorize a vacuous/off-rectangle collapse. **A** — matched nonzero leaves already ran at ~330 MB peak / ~80 MB largest with 0 deferrals; unmatched `R_N6` stays UNDECIDED by schema, not by heaviness.**QUESTION NEEDS-WORK.** **A** — matched nonzero leaves already computed at ~330 MB peak / ~80 MB largest with 0 deferrals; unmatched `R_N6` stays UNDECIDED by schema, not heaviness.
