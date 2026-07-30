# Census rows — stage016 (ℓ=2 SO(3) covariance), both engines

**Status:** builder output for ONE stage, per `CENSUS_SCHEMA.md`. ⛔ Nothing here is adjudicated
identity (§8.3 reserves merges for the physics review leg); the QIDs below are *minted*, not merged.
⛔ No dimensional-correctness verdict is reused as an oracle (§11): the stage's `21/0` and its
"12 typed / 9 computed" split are a different axis and appear below only where a locus is quoted.

**Artifacts censused (§7.2 — two artifacts):**

- `PY` := `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py`
- `WL` := `/var/projects/toy_physics/research/pde_ledger_v2/mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl`

Every `PY:NNN` / `WL:NNN` below expands to the full path above; no locus is a bare `:NNN`.
Other loci are written in full.

**Result set used (§7.1.1 step 1, fixed before classification):**
`/var/projects/toy_physics/research/pde_ledger_v2/notes/census/REPORTED_RESULTS.md` §1 —
`RES:016:L2-IRREP` (`REPORTED_RESULTS.md:64-80`) and `RES:016:K2-FORM` (`REPORTED_RESULTS.md:82-106`).
⛔ Not amended.

**Every locus cited below was opened and read** (§9.1 rule 1), including the parameter-register and
stage043 loci. Where a document's attribution did not check out it is recorded as an attribution
defect, not as a pass.

---

## 0. Builder decisions and selection rules (stated so coverage is checkable)

1. **Binding site = §7.2 as written**: the line of the assignment/assertion expression, head line for a
   multi-line expression. So `PY:218` (`theta, phi = sp.symbols(...)`) is ONE binding site carrying TWO
   quantities ⇒ two occurrences; `WL:21` (`$Assumptions`, a system global the artifact sets — §7.2
   rule 4) is ONE binding site carrying FIVE quantities ⇒ five occurrences; `PY:269` / `WL:165` (the
   harmonic dict/association heads) are ONE binding site each carrying FIVE quantities.
2. **Channel-indexed collections produced by one uniform loop at one binding site are ONE quantity**
   with `channels = 5` recorded (e.g. `lambdas` at `PY:289`). ⚠ The alternative reading — five
   quantities per collection — multiplies six row-groups by five and would take the in-universe count
   from **56 to 104**. Stated so the choice is visible; the tier *distribution* is unchanged by it.
3. **Out-of-scope granularity**: counted as *listed objects*, not binding sites (a 12-entry dimension
   table counts 12), so no exclusion hides inside a container. Groups carry exact counts and loci.
4. **`print` / `Print` emissions bind no value** and are not occurrences — with one exception: a prose
   statement of a *model premise* (route 2) is recorded as an occurrence, since route-2 propositions
   characteristically have no other binding site. See gap **S5**.
5. **Scope**: the two engines. The stage note is *not* censused as an artifact in this pass; its
   in-universe binding sites are enumerated in §7 below (14 occurrences), not dropped.

---

## 1. The universe and its reachability witnesses (§7.1.1 step 2)

Both admission routes were applied.

- **Route 1 (valued quantities)** reached: the harmonic literals, the channel set, the Gram matrix,
  `−Δ_S²Y`, `λ_m`, the eigenfunction residual, the degeneracy boolean, the typed `6` target, `K₂`, the
  extracted `T̃_Ω` coefficient, the K₂-coefficient residual, `M₂`, and the three free radial scalars
  and two coordinates reached as closure **terminals**.
- **Route 2 (constitutive propositions)**, admitted through the `premise-dependence` hop of §3.3.1(4):
  the S² domain/measure, the `−Δ_S²` angular operator, the wall 3-volume measure `dV = a²·dw·dΩ`, the
  `M₂` density integral, the `K₂` density integral, and the frozen-radial-calibration premise. Each was
  checked against **both** of §7.1 route 2's bounds; the checks are recorded on the rows.

**In-universe occurrences: 56** — **28 in `PY`, 28 in `WL`**. The two engines carry the *same* 28
occurrences (both engines bind every one of the 27 QIDs, and `QID:lambda_m` twice: computed and as a
typed assert target). ⚠ They differ in *how* two of them are bound — `WL:21` binds at one site what
`PY` binds at four (`PY:218-221`), and `WL:164` **types** the channel set that `PY:279` **computes** —
and the second of those is the stage's only cross-engine classification divergence (R11).

**Hop kinds used:** `claim-content` (the result's own value-bearing locus), `expression` (input
closure), `premise-dependence` (§3.3.1(4)).

---

## 2. IN-UNIVERSE ROWS

Field order per block: quantity · occurrences (OID) · witness · axis A + §9 evidence · axis B + §9
evidence · axis C (leaves → `PHYSICS-FED`, `SELF-REFERENTIAL`) · `CONVENTION-LADEN` + basis ·
`is_tier` · `should_be_tier` · `should_be_basis` · `delta` · LIVE/RETIRED.

### R01 — θ (polar angle)

- **Occurrences (2):** `QID:theta@PY#218`, `QID:theta@WL#21`.
- **Witness:** `RES:016:L2-IRREP`, `expression` — `PY:289 → PY:287 → PY:233-235/PY:241-242`;
  `WL:210 → WL:180-186 → WL:109-111/WL:116-118`.
- **Axis A:** `A-INDEPENDENT-VARIABLE`. §9 evidence (loci where the ledger *uses* it as one):
  integrated over at `PY:233-235` (`sp.integrate(..., (theta, 0, sp.pi))`) and `WL:109-111`;
  differentiated with respect to at `PY:241-242` (`sp.diff(sp.sin(theta)*sp.diff(expr, theta), theta)`)
  and `WL:116-118`.
- **Axis B:** `B-DECLARED-UNASSIGNED`. Evidence: declaration `PY:218`; no assignment of `theta`
  anywhere in `PY` (only the declaration and uses). `WL:21` declares `Element[{theta, phi}, Reals]`;
  no `Set` of `theta` in `WL`.
- **Axis C:** closure terminates immediately at the occurrence's own declaration; leaf `C-FREE` ⇒
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined (§3.3). ⚠ See gap **S1** — the leaf table
  has no tag for a coordinate.
- **CONVENTION-LADEN:** `false` — *not a candidate*: no transformation group is claimed for this
  occurrence anywhere in the corpus (a chart-change group would make it one; none is documented).
- **is_tier:** `no-tier:independent-variable` (§5.7; the axis-C state does not move it, §3.3).
- **should_be_tier:** `no-tier:independent-variable` · **basis:** `none` · **delta:** no.
- **LIVE.**

### R02 — φ (azimuth)

Identical to **R01** in every field. **Occurrences (2):** `QID:phi@PY#218`, `QID:phi@WL#21`.
Evidence differs only in the loci: integrated over `(phi, 0, 2*sp.pi)` at `PY:233` / `WL:109`;
differentiated at `PY:242` (`sp.diff(expr, phi, 2)`) / `WL:117`.

### R03 — M̃ (ℓ=2 radial mass scalar)

- **Occurrences (2):** `QID:M_tilde@PY#219`, `QID:M_tilde@WL#21`.
- **Witness:** `RES:016:K2-FORM`, `expression` — `PY:507 (m2_core = Mtilde) → PY:219`;
  `WL:228 → WL:21`.
- **Axis A:** `A-REDUCIBLE-UNDERIVED`. §9 evidence — **named route**: R35, recorded at
  `/var/projects/toy_physics/research/pde_ledger_v2/notes/parameter_register.md:302`
  ("the grouped-lane scalars **= ∫ density·β₂² dV** — `M̃=∫μ_η β₂²` …") and at `:184`.
  **Executable within this framework (§3.1.2):** the integrand is already assembled in-artifact at
  `PY:336` / `WL:232`, and stage043 records the route as *un-run*, not *un-buildable* —
  `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:390`
  carries `moment_integral_executed=False` with `source_status="DERIVED-IN-FORM-UNEXECUTED"` (`:384`)
  for `REG:C1:Mtilde` (`:371`).
- **Axis B:** `B-DECLARED-UNASSIGNED`. Evidence: declaration `PY:219`
  (`Mtilde = sp.Symbol("Mtilde", positive=True, real=True)`); the artifact contains **no assignment** of
  it — every other `Mtilde` line in `PY` is a use (`:323`, `:475`, `:507`, `:722`, `:835`). In `WL` the
  only binding is the positivity assertion at `WL:21`; no `Set`.
- **Axis C:** leaf `C-FREE` (the artifact never assigns it, nothing it stands for is fixed here) ⇒
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
  ⛔ **No `C-PEER` leaf.** `PY:835` / `WL:646` claim "CONSUMED-from-011/012/013 … Mtilde/Ktilde/
  TomegaTilde … cited as provenance" — that names *stages*, which §3.3 says is not a citation.
- **⚠ Attribution defect (§9.1 rule 2), measured:** the attributed source does not carry the symbols at
  all. `grep` over `notes/stages/ledger_stage011*.md`, `…012*.md`, `…013*.md` and their `scripts/`
  counterparts returns **no** occurrence of `Mtilde|Ktilde|TomegaTilde|M̃|K̃|T̃`. Both loci recorded:
  attributed `PY:835` / `WL:646`; computing locus — **none exists in the ledger** (see the value note
  below).
- **⚠ Second attribution defect:** `parameter_register.md:184` cites the dimension of `M̃, K̃, T̃_Ω` as
  "computed in stage016's SymPy dimension rules, `scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:355-366`".
  Opened: `PY:355-366` is the tail of `dimension_eval`'s `ok = bool(...)` conjunction and the head of
  its return dict. The twelve rule literals are at **`PY:314-325`**. The stage note repeats the stale
  range at `notes/stages/ledger_stage016_l2_so3_covariance.md:194` (`sympy:355-366`) and at `:210`
  (`sympy:364`). Values agree; pointers do not.
- **CONVENTION-LADEN:** `false` — not a candidate; no transformation group is claimed for the scalar.
- **is_tier:** `tier1-debt` (§9.0: an unresolved closure does not touch an evidenced axis A).
- **should_be_tier:** `tier3-emergent` · **basis:** `named-route` (R35, recorded at
  `parameter_register.md:302`) · **delta:** YES.
- **⚑ Conflict (§10.3, intra-occurrence)** — the substrate carries incompatible classes for this
  quantity, recorded with both/all loci and ⛔ not resolved here:
  `parameter_register.md:184` and `:302` class it **`DERIVED`**; `parameter_register.md:353` (R87) lists
  it in `C1` as "frozen as `calibration_inputs` in the stage-017 build with the moment-integral never
  evaluated to bind them"; `:357` (R91) as "written-as-DERIVED-form but un-executed"; this stage's own
  note `:311-314` and `PY:828` / `WL:639` call the radial scalars "FROZEN calibration inputs". stage043
  carries the same quantity under **two** IDs in two categories: `REG:derived:Mtilde_definition`
  (`…stage043…py:298`, `CAT_DERIVED` = *derived-not-counted*) and `REG:C1:Mtilde` (`:371`, `CAT_OPEN` =
  *extension-convention-open*). The stage note's H1 (`:210-213`) additionally records R35's
  "dual-engine dim-verified" as **overstated**.
- **Value note (not Route C):** the QID is *not* `unvalued-in-universe` — stage017 assigns it numerically
  at `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage017_grouped_p2_lane_isotropy_sympy_audit.py:41`
  (`"Mtilde": 3.0`) and `mathematica/ledger_stage017_grouped_p2_lane_isotropy_mathematica_audit.wl:273`.
  Those are stage017 occurrences, out of this stage's scope.
- **LIVE** (no retirement marker; contrast the marked retired rows `parameter_register.md:139`, `:159`).

### R04 — K̃ (ℓ=2 radial stiffness scalar)

**Occurrences (2):** `QID:K_tilde@PY#220`, `QID:K_tilde@WL#21`. All fields as **R03**, with:
witness `RES:016:K2-FORM`, `expression` — `PY:499 → PY:305 → PY:220`; `WL:219 → WL:216 → WL:21`.
Route-R35 evidence line for `K̃`: `parameter_register.md:302` (`K̃=∫[T_w β₂'²+K_η β₂²]`);
stage043 `REG:C1:Ktilde` (`…stage043…py:371`) and `REG:derived:Ktilde_definition` (`:299`).
**is_tier** `tier1-debt` · **should_be** `tier3-emergent` / `named-route` · **delta** YES · **LIVE**.

### R05 — T̃_Ω (ℓ=2 reduced angular-stiffness scalar)

**Occurrences (2):** `QID:T_Omega_tilde@PY#221`, `QID:T_Omega_tilde@WL#21`. All fields as **R03**, with:
witness `RES:016:K2-FORM`, `expression` — `PY:499 → PY:305 → PY:221`; `WL:219 → WL:216 → WL:21`.
Route-R35 evidence: `parameter_register.md:302` (`T̃_Ω=∫T_Ω β₂²`); stage043 `REG:C1:Ttilde_Omega`
(`:371`) and `REG:derived:Ttilde_definition` (`:299`).
⚠ Additionally carried, not adjudicated here: the stage note's identity hazard `T̃_Ω`(016) vs `T_Ω`(023)
— "different quantity until R42 exists" (`notes/stages/ledger_stage016_l2_so3_covariance.md:187-192`).
⛔ Left as two QIDs; §8.3 reserves the merge for the physics review leg.
**is_tier** `tier1-debt` · **should_be** `tier3-emergent` / `named-route` · **delta** YES · **LIVE**.

### R06–R10 — the five real ℓ=2 harmonics `Y20, Y21c, Y21s, Y22c, Y22s`

- **Occurrences (10):** `QID:Y_20@PY#269` … `QID:Y_22s@PY#269` (five quantities, one binding site; value
  loci `PY:270-274`) and `QID:Y_20@WL#165` … `QID:Y_22s@WL#165` (value loci `WL:166-170`).
- **Witness:** `RES:016:L2-IRREP`, `expression` — `PY:289 → PY:287 → PY:286 → PY:269`;
  `WL:210 → WL:180-186 → WL:165`. Also reached from `RES:016:K2-FORM` via `PY:501`/`WL:221`.
- **Axis A:** `A-UNADJUDICATED` (§3.1.1 third branch). Recorded reason: **no route is recorded anywhere
  in the ledger** for the basis (searched `notes/stages/ledger_stage011*|012*|013*`, `docs/model_map.md`,
  `notes/parameter_register.md`; `parameter_register.md:301` (R34) *states* the harmonics form an
  orthonormal SO(3) irrep and derives nothing), and it is **not stated as a postulate about the medium**
  at any locus, so neither `A-REDUCIBLE-UNDERIVED` (§9: named route + where recorded) nor
  `A-IRREDUCIBLE-POSTULATE` (§9: where posited) nor `A-IRREDUCIBLE-STRUCTURAL` (§9: the foreclosing
  framework property — nothing forecloses writing the harmonics down) can be evidenced. ⚠ See gap **S2**.
- **Axis B:** `B-DECLARED-LITERAL` — typed closed forms at `PY:270-274` / `WL:166-170`.
- **Axis C:** leaf `C-SELF` ⇒ `PHYSICS-FED = false`; `SELF-REFERENTIAL = false` (plain literal, §3.3).
  `C-SELF`-over-`C-MATH` per §3.3's test: the artifact *could* have declared a different orthonormal
  ℓ=2 basis (any SO(3)-rotated one) and stayed internally consistent — `λ_m` and `Gram=I₅` both survive.
- **CONVENTION-LADEN:** `UNADJUDICATED`. A normalisation claim **is** implied — the artifact's own
  orthonormality gate (`PY:676-677`, `WL:462-463`) and `REPORTED_RESULTS.md:111` ("a property of the
  sphere and of the basis this artifact writes down") — but §3.4 clause (a)'s transformation group is
  documented **nowhere** for this occurrence, and clause (b)'s invariance is demonstrated nowhere. ⛔ Not
  set `true` by intuition (§3.4). ⇒ counted in `convention-unadjudicated`; buys **no** exclusion.
- **is_tier:** `no-tier:unadjudicated` (flag `UNADJUDICATED` ⇒ the axis-A projection runs, §5.7).
- **should_be_tier:** `no-tier:convention` · **basis:** `convention-candidate` · **delta:** YES.
  (⛔ Sets nothing and excludes nothing — §6.1's guard.)
- **LIVE.**

### R11 — the ℓ=2 channel set (the "5-dimensional" half of `RES:016:L2-IRREP`)

- **Occurrences (2):** `QID:l2_channel_set@PY#279`, `QID:l2_channel_set@WL#164`.
- **Witness:** `RES:016:L2-IRREP`, `claim-content` (the claim's "single 5-dimensional" clause) +
  `expression` (`PY:279 = list(harmonics)`).
- **⭐ The two engines differ, and this is the only such divergence in the stage.**
  - `PY:279` `order = list(harmonics)` — **`B-EXECUTED`** (evidence: the computation `PY:279`; input
    leaves `PY:269`/`:270-274`). **Axis A `A-REDUCED`** — §9: the reduction is performed at `PY:279`
    and reduces to the basis declaration at `PY:269`. **Axis C:** leaves `{C-SELF ×5}` ⇒
    `PHYSICS-FED = false`, `SELF-REFERENTIAL = false`. ⇒ **is_tier `no-tier:unclassified-nonfed`**
    (§5.1) and §4 near-miss **executed-but-not-physics-fed**.
  - `WL:164` `order = {"20","21c","21s","22c","22s"}` — **`B-DECLARED-LITERAL`** (typed).
    **Axis A `A-UNADJUDICATED`** (no recorded route to `2ℓ+1`; same evidence as R06–R10).
    **Axis C:** leaf `C-SELF` ⇒ `PHYSICS-FED = false`. ⇒ **is_tier `no-tier:unadjudicated`**.
- **CONVENTION-LADEN:** `false` (not a candidate).
- **should_be_tier:** `tier3-emergent` (both) · **basis:** `physical-picture-expectation` (the sector's
  dimension ought to follow from the model's angular sector; no route is recorded) · **delta:** YES both.
- **§10.2.2 rule 5:** the QID mixes `no-tier:unadjudicated` with `no-tier:unclassified-nonfed` ⇒
  reported in **both** buckets, flagged `mixed-adjudication`, and ⛔ **not** entered in the conflict set.
- **LIVE.**

### R12 — the S² domain and measure `dΩ = sinθ dφ dθ`, `θ∈[0,π]`, `φ∈[0,2π]` (route 2)

- **Occurrences (2):** `QID:S2_domain_measure@PY#230` (body `PY:233-235`),
  `QID:S2_domain_measure@WL#107` (body `WL:109-111`).
- **Witness:** `RES:016:L2-IRREP`, `premise-dependence` (§3.3.1(4)) — the Rayleigh quotient and the
  eigenvalue are defined *with respect to* this measure (`PY:287`, `WL:183`).
- **Route-2 bounds, both checked (§9):** (i) depended on by `RES:016:L2-IRREP` and `RES:016:K2-FORM`;
  (ii) **not derived elsewhere in the ledger** — stages 011/012/013 carry no angular sector at all
  (`notes/stages/ledger_stage011_frozen_reduction_certificate.md`,
  `…stage012_dtn_pole_ladder_robin.md`, `…stage013_breathing_harmonic_mk_projection.md`: searched for
  `S²|angular|Laplace|harmonic`; the hits are radial/ℓ=0 only), and `docs/model_map.md:94` merely
  *restates* this stage's own `−Δ_S²Y_2m=6Y_2m`.
- **Axis A:** `A-UNADJUDICATED` (§3.1.1 third branch). The `STRUCTURAL`/`POSTULATE` call **cannot be
  decided from evidence**: no locus states the throat's spherical symmetry as a defining property of the
  medium (`docs/model_map.md:23-31` describes the throat without asserting it), and no route to derive
  the wall's angular geometry is recorded. ⛔ Not resolved by picking the better-reading branch.
  ⚠ **Sensitivity, stated:** a `STRUCTURAL` reading is available with the extension named — the deferred
  throat-interior solve (`parameter_register.md:49` "all siblings of the ONE deferred throat-interior
  solve"; `notes/stages/ledger_stage016_l2_so3_covariance.md:315-316` "`G = GENUINE_BLOCKED`") plus a
  wall action over the throat's angular geometry, which the ledger does not posit (`docs/model_map.md:70`
  records the two existing wall actions as un-reconciled and neither is applied to this sector).
  Adopting it would move **4 occurrences** (this row + R13) from the tier-1 upper span to its lower bound.
- **Axis B:** `B-DECLARED-LITERAL` — the measure and both limits are typed at `PY:233-235` / `WL:109-111`.
  (Route-2's `B-POSTULATED` reading changes no §4 bucket — see gap **S6**.)
- **Axis C:** §3.3.1(4) — the closure terminates where the artifact itself declares it ⇒ leaf `C-SELF` ⇒
  `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN:** `false` — no transformation group is claimed for the *S²* measure. (The **wall**
  measure is a different occurrence and IS convention-claimed — see R18.)
- **is_tier:** `no-tier:unadjudicated` · **should_be_tier:** `tier3-emergent` ·
  **basis:** `physical-picture-expectation` · **delta:** YES · **LIVE.**

### R13 — the angular operator `−Δ_S²` (route 2)

- **Occurrences (2):** `QID:angular_operator_neg_Delta_S2@PY#239` (body `PY:241-242`),
  `QID:angular_operator_neg_Delta_S2@WL#115` (body `WL:116-118`).
- **Witness:** `RES:016:L2-IRREP`, `premise-dependence` — the reported spectrum is the spectrum *of this
  operator* (`REPORTED_RESULTS.md:66-68`).
- **Route-2 bounds:** as R12 — depended on by the result; not derived elsewhere in the ledger.
- **⭐⭐ The stage's decisive axis-C determination, recorded here once and referenced by every row
  downstream.** Under §3.3.1(1) the operator is a `C-FIELDEQ` leaf **iff the artifact cites the locus of
  the equation being solved**. It does not: `PY:2-10` names "pathA_32 II-G3a" (a *source name*, which
  §3.3 says is not a citation), `WL:1-8` the same; the stage note's provenance
  (`notes/stages/ledger_stage016_l2_so3_covariance.md:16-18`) cites the predecessor **implementation**
  files and a report, and states "the derivation below is **self-contained**". ⇒ **the operator
  contributes no leaf**, and every closure it appears in is `{C-SELF, …}`. ⚠ This is what makes the
  whole angular chain (R14–R17, R11-PY) `PHYSICS-FED = false` rather than `true`.
- **Axis A:** `A-UNADJUDICATED`, same evidence and same stated sensitivity as R12.
- **Axis B:** `B-DECLARED-LITERAL` — the operator is typed out in coordinates at `PY:241-242` /
  `WL:116-118`.
- **Axis C:** leaf `C-SELF` ⇒ `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN:** `false` (not a candidate).
- **is_tier:** `no-tier:unadjudicated` · **should_be_tier:** `tier3-emergent` ·
  **basis:** `physical-picture-expectation` · **delta:** YES · **LIVE.**

### R14 — the Gram matrix (orthonormality of the five channels)

- **Occurrences (2):** `QID:gram_l2@PY#281`, `QID:gram_l2@WL#208` (computed `WL:178`).
- **Witness:** `RES:016:L2-IRREP`, `premise-dependence` — linear independence underwrites the
  "5-dimensional" half of the claim. ⚠ **Contestable, flagged:** `REPORTED_RESULTS.md:111` demotes
  `Gram = I₅` from the *result* set; this row admits it as a *premise* of a result, not as a result.
  If a reviewer rejects the hop, 2 occurrences move to `reached-by-no-reported-result`.
- **Axis A:** `A-REDUCED`. §9: the reduction is performed at `PY:281` /`WL:178` and reduces to the
  harmonics (`PY:269`/`:270-274`; `WL:165`/`:166-170`) and the S² measure (`PY:233-235`; `WL:109-111`).
- **Axis B:** `B-EXECUTED`. §9: computation `PY:281` (`integrate_s2(ys[i]*ys[j])`) / `WL:178`; input-leaf
  loci as above. Genuine symbolic integration (`sp.integrate` / `Integrate`), not a lookup.
- **Axis C:** leaves `{C-SELF (harmonics), C-SELF (measure)}` ⇒ `PHYSICS-FED = false`;
  `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN:** `false` — the convention question attaches to the basis (R06–R10), not to the
  computed matrix.
- **is_tier:** `no-tier:unclassified-nonfed` (§5.1) — §4 near-miss **executed-but-not-physics-fed**.
- **should_be_tier:** `no-tier:unclassified-nonfed` · **basis:** `none` · **delta:** no.
  (Orthonormality on S² is a fact about the sphere and the artifact's basis; the model's picture makes no
  claim that it should reduce to the medium.)
- **LIVE.**

### R15 — `−Δ_S² Y_m` (per channel)

- **Occurrences (2):** `QID:neg_laplacian_Y@PY#288` (computed `PY:286`), `QID:neg_laplacian_Y@WL#209`
  (computed `WL:179`). `channels = 5`.
- **Witness:** `RES:016:L2-IRREP`, `expression` (`PY:287 → PY:288`; `WL:183 → WL:209`).
- **Axis A:** `A-REDUCED`. §9: reduction at `PY:286` / `WL:179`; reduces to the harmonics (`PY:269`;
  `WL:165`) and the operator (`PY:239-243`; `WL:115-119`).
- **Axis B:** `B-EXECUTED`. §9: computation `PY:286` (`compact(-laplacian_s2(y_expr))`, native
  `sp.diff`) / `WL:179` (native `D`); input leaves as above.
- **Axis C:** leaves `{C-SELF}` — the operator contributes **no** leaf (R13) ⇒ `PHYSICS-FED = false`;
  `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unclassified-nonfed` — near-miss 1.
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` (it becomes
  physics-fed exactly when the operator is a cited model equation — R13) · **delta:** YES · **LIVE.**

### R16 — `λ_m` (the `−Δ_S²` eigenvalue), computed

- **Occurrences (2 of this QID's 4):** `QID:lambda_m@PY#289` (computed `PY:287`),
  `QID:lambda_m@WL#210` (computed `WL:180-186`). `channels = 5`.
- **Witness:** `RES:016:L2-IRREP`, `claim-content` (`REPORTED_RESULTS.md:73` cites `PY:287`).
- **Axis A:** `A-REDUCED`. §9: the reduction is the Rayleigh quotient performed at `PY:287`
  (`integrate_s2(y_expr*neg_lap)/integrate_s2(y_expr**2)`) / `WL:183`; it reduces to `PY:288` (`−ΔY`),
  the harmonics `PY:269` and the measure `PY:233-235` (`WL:209`, `WL:165`, `WL:109-111`).
- **Axis B:** `B-EXECUTED`. §9: computation loci as above; input-leaf loci as above.
- **Axis C:** leaves `{C-SELF (harmonics), C-SELF (measure)}` ⇒ `PHYSICS-FED = false`;
  `SELF-REFERENTIAL = false` (the walk leaves the occurrence and terminates at the basis; it does not
  return to `PY:289`).
- **CONVENTION-LADEN:** `false` (not a candidate).
- **is_tier:** `no-tier:unclassified-nonfed` — §4 near-miss **executed-but-not-physics-fed**. ⚠ Code
  genuinely ran and can genuinely fail; what it ran over is the artifact's own declarations (§3.5).
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` · **delta:** YES.
- **Substrate note (not a conflict):** stage043 carries `REG:derived:lambda_m_SO3` in `CAT_DERIVED`
  (`…stage043…py:301`, `:139`). ⛔ Not inherited (§2); recorded because the census's axis-C finding is
  precisely what that label cannot express.
- **LIVE.**

### R17 — `λ_m` asserted against the typed literal `6`

- **Occurrences (2 of this QID's 4):** `QID:lambda_m@PY#679`, `QID:lambda_m@WL#466`.
- **Witness:** `RES:016:L2-IRREP`, `claim-content` (`REPORTED_RESULTS.md:73` cites the assert block
  `PY:679-682`).
- **Axis A:** `A-REDUCED` — the *quantity* is obtained, in this same artifact, at `PY:287` / `WL:183`
  (§9 loci as R16). Axis A is about the quantity, not about this occurrence's execution (§3).
- **Axis B:** `B-ASSERTED-TARGET` — `PY:679` `expect_zero(..., lambdas[name] - 6)` / `WL:466`
  `expectZero[..., lambdas[name] - 6]`: the artifact checks *against* the typed `6`; it does not produce
  it here.
- **Axis C:** leaf `C-SELF` (the typed `6`) ⇒ `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
  ⚠ `C-SELF` over `C-MATH` by §3.3's tie-break: neither engine anywhere writes `ℓ(ℓ+1)`; the identity is
  stated only in prose (`notes/stages/ledger_stage016_l2_so3_covariance.md:52`,
  `REPORTED_RESULTS.md:78`), so no cited-or-executed mathematical identity forces the literal.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unclassified-nonfed`.
  ⛔ In **no** §4 near-miss (near-miss 1 requires `B-EXECUTED`; near-miss 3 requires `PHYSICS-FED`).
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` · **delta:** YES.
- **LIVE.**

### R18 — the m-degeneracy (`all five λ = 6`)

- **Occurrences (2):** `QID:l2_m_degeneracy@PY#299`, `QID:l2_m_degeneracy@WL#213` (computed `WL:202`).
- **Witness:** `RES:016:L2-IRREP`, `claim-content` (`REPORTED_RESULTS.md:73` cites `PY:299`
  `lambda_all_six`). This boolean **is** the degeneracy that the result claims.
- **Axis A:** `A-REDUCED`. §9: reduction at `PY:299` / `WL:202`; reduces to the computed λ set
  (`PY:289`, `WL:210`).
- **Axis B:** `B-EXECUTED` (the five comparisons run). ⚠ Recorded on the row: the comparison **target
  `6` is typed at the same binding site** (`PY:299` `compact(value - 6)`; `WL:202` `clean[# - 6]`); per
  §3.3 a sub-expression within one binding site is not a separate step, so no second occurrence.
- **Axis C:** leaves `{C-SELF}` through R16 ⇒ `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unclassified-nonfed` — near-miss 1.
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` · **delta:** YES.
- **LIVE.**

### R19 — the eigenfunction residual `(−Δ_S²)Y_m − λ_m Y_m`

- **Occurrences (2):** `QID:eigenfunction_residual@PY#290`, `QID:eigenfunction_residual@WL#211`
  (computed `WL:187-193`). `channels = 5`.
- **Witness:** `RES:016:L2-IRREP`, `claim-content` (`REPORTED_RESULTS.md:73` cites `PY:290`).
- **Axis A:** `A-REDUCED`. §9: reduction at `PY:290` / `WL:190`; reduces to `PY:288`, `PY:289`,
  `PY:269` (`WL:209`, `WL:210`, `WL:165`).
- **Axis B:** `B-EXECUTED`. §9: computation `PY:290` / `WL:190`; input leaves as above.
- **Axis C:** leaves `{C-SELF}` ⇒ `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unclassified-nonfed` — near-miss 1.
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` · **delta:** YES.
- **LIVE.**

### R20 — the K₂ assembly form `K₂ = K̃ + coeff·T̃_Ω`

- **Occurrences (2):** `QID:K2_assembly_form@PY#305`, `QID:K2_assembly_form@WL#216`.
- **Witness:** `RES:016:K2-FORM`, `claim-content` (`REPORTED_RESULTS.md:91` cites `PY:304-305`).
- **Axis A:** `A-REDUCIBLE-UNDERIVED`. §9 — **named route**: the density-integral decomposition of the
  ℓ=2 stiffness, stated at `notes/stages/ledger_stage016_l2_so3_covariance.md:76-79` and carried as edge
  R35 at `parameter_register.md:302`; the integrand is already assembled in-artifact at `PY:337-340` /
  `WL:233-236`. **Executable within this framework:**
  `…stage043…py:384`/`:390` record the moment integral as `DERIVED-IN-FORM-UNEXECUTED` /
  `moment_integral_executed=False` — machinery present, run absent (§3.1.2's debt branch).
- **Axis B:** `B-DECLARED-LITERAL` — the expression tree `Ktilde + coeff*TomegaTilde` is typed at
  `PY:305` / `WL:216`. ⛔ **Not** `B-DERIVED-IN-FORM-UNEXECUTED`: no symbolic derivation of the form
  exists in either engine — the stage's own H1 records "**Neither engine anywhere writes**
  `M̃=∫μ_ηβ₂²dV`; it lives only in print strings" (`notes/stages/ledger_stage016_l2_so3_covariance.md:210-213`).
- **Axis C:** leaves `{C-FREE (K̃), C-FREE (T̃_Ω)}`; no leaf establishes physics-feeding ⇒
  `PHYSICS-FED = C-UNRESOLVED` (§3.3); `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN:** `false` (not a candidate).
- **is_tier:** `tier1-debt` — §9.0: the unresolved closure gates tier 3 and `DERIVED`, and does **not**
  touch the evidenced axis A.
- **should_be_tier:** `tier3-emergent` · **basis:** `named-route` (R35, `parameter_register.md:302`) ·
  **delta:** YES.
- **⚠ Carried, not folded into any axis:** the stage's adversarial finding **H8**
  (`…ledger_stage016_l2_so3_covariance.md:259-270`) — dimension-preserving rewrites of the walked
  expressions, *including dropping `λ_m` from the angular term*, leave the audit green. That bears on how
  well the artifact protects this row, not on its provenance.
- **LIVE.**

### R21 — `K₂` (assembled value, per channel)

- **Occurrences (2):** `QID:K2_l2@PY#499`, `QID:K2_l2@WL#219`. `channels = 5`.
- **Witness:** `RES:016:K2-FORM`, `claim-content` (`REPORTED_RESULTS.md:91` cites `PY:499`).
- **Axis A:** `A-REDUCED`. §9: the reduction is performed at `PY:499`
  (`build_K2(lambdas[name])`) / `WL:219`, and reduces to `PY:220` (K̃), `PY:221` (T̃_Ω) and `PY:289`
  (the live computed λ) — `WL:21`, `WL:21`, `WL:210`.
- **Axis B:** `B-EXECUTED`. §9: computation `PY:499` → `PY:305`; input leaves as above. The λ consumed
  is the live computed one, not a literal (`PY:499` passes `lambdas[name]`).
- **Axis C:** leaves `{C-FREE (K̃), C-FREE (T̃_Ω), C-SELF (the λ chain)}` ⇒
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN:** `false`.
- **is_tier:** `no-tier:unadjudicated` — §5.7's one case where an unresolved closure decides the line:
  `A-REDUCED` projects only to tier 3, that tier's `PHYSICS-FED` conjunct fails, and there is no other
  tier to fall to. ⛔ Not `no-tier:unclassified-nonfed` (that is a positive finding).
  ⇒ contributes to the **tier-1 upper span**.
- **should_be_tier:** `tier3-emergent` · **basis:** `named-route` (R35) · **delta:** YES.
- ⛔ **Not** in near-miss 1: §3.3 forbids asserting a `C-UNRESOLVED` row into
  `executed-but-not-physics-fed`.
- **LIVE.**

### R22 — the extracted `T̃_Ω` coefficient of the assembled `K₂`

- **Occurrences (2):** `QID:K2_TomegaTilde_coefficient@PY#500`, `QID:K2_TomegaTilde_coefficient@WL#220`.
- **Witness:** `RES:016:K2-FORM`, `claim-content` (`REPORTED_RESULTS.md:91` cites `PY:500-504`).
- **Axis A:** `A-REDUCED`. §9: reduction at `PY:500` → `PY:308-309` (`sp.diff(k2_expr, TomegaTilde)`) /
  `WL:220` → `WL:217`; reduces to `PY:499` (`WL:219`).
- **Axis B:** `B-EXECUTED`. §9: computation `PY:500`/`WL:220`; input leaf `PY:499`/`WL:219`.
- **Axis C:** the producing expression consumes `k2_core`, whose closure carries `C-FREE` leaves ⇒
  `PHYSICS-FED = C-UNRESOLVED`. ⚠ **Gap S3:** the differentiation *annihilates* K̃, so the C-FREE leaf is
  in the input closure but contributes nothing to the value (which is exactly λ). The schema has no rule
  for that; the letter of §3.3 ("walk back the expression that produced its value") was followed.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unadjudicated` (as R21).
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` · **delta:** YES.
- **LIVE.**

### R23 — the K₂-coefficient residual `(−Δ)Y − coeff·Y`

- **Occurrences (2):** `QID:K2_coefficient_residual@PY#501` (expression `PY:501-503`),
  `QID:K2_coefficient_residual@WL#221` (expression `WL:221-224`). `channels = 5`.
- **Witness:** `RES:016:K2-FORM`, `claim-content` (`REPORTED_RESULTS.md:91` cites `PY:500-504`).
- **Axis A:** `A-REDUCED`. §9: reduction at `PY:502` / `WL:223`; reduces to `PY:288`, `PY:500`,
  `PY:269` (`WL:209`, `WL:220`, `WL:165`).
- **Axis B:** `B-EXECUTED`. §9: computation `PY:502` / `WL:223`; input leaves as above.
- **Axis C:** inherits the `C-FREE` leaf through `PY:500` ⇒ `PHYSICS-FED = C-UNRESOLVED`
  (same gap **S3**); `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unadjudicated`.
- **should_be_tier:** `tier3-emergent` · **basis:** `physical-picture-expectation` · **delta:** YES.
- **LIVE.**

### R24 — `M₂ = M̃`

- **Occurrences (2):** `QID:M2_l2@PY#507`, `QID:M2_l2@WL#228`.
- **Witness:** `RES:016:K2-FORM`, `claim-content` — the claim states `M₂ = M̃`
  (`REPORTED_RESULTS.md:84-85`). ⚠ **Recorded seam (gap S9):** the result set's *script*-locus list
  (`REPORTED_RESULTS.md:91`) does **not** cite `PY:507` / `WL:228`, although the claim contains the
  relation. ⛔ The result set was not amended; the row is admitted on the claim's content.
- **Axis A:** `A-REDUCED`. §9: the reduction `M₂ → M̃` is performed at `PY:507` / `WL:228` and is
  recorded at `notes/stages/ledger_stage016_l2_so3_covariance.md:57-58`; it reduces to `PY:219` /
  `WL:21`.
- **Axis B:** `B-DECLARED-LITERAL` — `m2_core = Mtilde` is typed; no computation produces it.
- **Axis C:** leaf `{C-FREE (M̃)}` ⇒ `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN:** `false`. **is_tier:** `no-tier:unadjudicated` (as R21).
- **should_be_tier:** `tier3-emergent` · **basis:** `named-route` (R35's `M̃=∫μ_η β₂²`,
  `parameter_register.md:302`) · **delta:** YES.
- ⛔ **Not** near-miss 3: that requires `PHYSICS-FED = true`.
- **LIVE.**

### R25 — the wall 3-volume measure `dV = a²·dw·dΩ` (route 2)

- **Occurrences (2):** `QID:wall_volume_measure@PY#335`, `QID:wall_volume_measure@WL#231`.
- **Witness:** `RES:016:K2-FORM`, `premise-dependence` — the measure the ℓ=2 radial reduction integrates
  over; without it `K₂ = K̃ + 6T̃_Ω` is algebra on two free symbols.
- **Route-2 bounds:** (i) depended on by `RES:016:K2-FORM` (through R26/R27); (ii) not derived elsewhere
  — the register carries it as pathA_32's *convention*, not as a derivation
  (`notes/stages/ledger_stage016_l2_so3_covariance.md:74`, `:322-324`).
- **Axis A:** `A-UNADJUDICATED`, same reasoning and same stated sensitivity as R12 (the wall's geometry
  is contingent on the deferred throat-interior solve; nothing states it as a postulate).
- **Axis B:** `B-DECLARED-LITERAL` — typed at `PY:335` (`a_dim**2 * dw_dim * dOmega_dim`) / `WL:231`.
- **Axis C:** §3.3.1(4) terminal, artifact-declared ⇒ leaf `C-SELF` ⇒ `PHYSICS-FED = false`.
- **CONVENTION-LADEN:** **`UNADJUDICATED`** — a convention claim **is** made: the stage calls it
  "pathA_32's own convention — VOLUME densities on the wall measure `dV = a²·dw·dΩ`"
  (`…ledger_stage016_l2_so3_covariance.md:74`) and contrasts it with stage013's LINE-density convention
  at `:322-324`. §3.4 clause (a): no transformation group is documented — the same locus states the two
  conventions are related by "an `∫a²dΩ ≈ L²` bridge, **not equal**". Clause (b): the corpus records
  **non**-invariance ("the stage-013 relation `K_η = T_w β²` does **NOT** transfer"), so the invariance is
  not merely undemonstrated. ⇒ `convention-unadjudicated` bucket; buys **no** exclusion; the row keeps
  its axis-A tier (§3.4, §9.0).
- **is_tier:** `no-tier:unadjudicated` · **should_be_tier:** `no-tier:convention` ·
  **basis:** `convention-candidate` · **delta:** YES · **LIVE.**

### R26 — the `M₂` density integral `M₂ = ∫ μ_η β₂² dV` (route 2)

- **Occurrences (2):** `QID:M2_density_integral_premise@PY#336`,
  `QID:M2_density_integral_premise@WL#232`.
- **Witness:** `RES:016:K2-FORM`, `premise-dependence` (§3.3.1(4)) — this is what makes `M₂ = M̃` a claim
  about the wall rather than a renaming. ⛔ The hop confers no `PHYSICS-FED` on anything it passes
  through.
- **Route-2 bounds:** (i) depended on by `RES:016:K2-FORM`; (ii) not derived elsewhere in the ledger —
  no artifact derives the relation from an action; `parameter_register.md:302` *states* it as edge R35
  and `…stage043…py:390` records it un-executed.
- **Axis A:** `A-UNADJUDICATED` — the relation would follow from a wall action over the throat's angular
  geometry; the ledger posits none (evidence as R12). ⛔ Not `A-REDUCED` on the strength of the
  register's `DERIVED` label (§2: no row inherits a substrate class).
- **Axis B:** `B-DECLARED-LITERAL` — the integrand is typed at `PY:336` / `WL:232` and is **walked only
  for dimensions**; no code path integrates it (H1, `…ledger_stage016_l2_so3_covariance.md:210-213`).
  ⛔ Not `B-DERIVED-IN-FORM-UNEXECUTED`: no symbolic derivation is present in either engine, only the
  integrand product.
- **Axis C:** §3.3.1(4) terminal, artifact-declared ⇒ leaf `C-SELF` ⇒ `PHYSICS-FED = false`.
  ⚠ Gap **S8**: the terminal rule means the nine density/measure symbols named inside this expression
  never become leaves.
- **CONVENTION-LADEN:** `false` (the convention claim attaches to the measure — R25).
- **is_tier:** `no-tier:unadjudicated` · **should_be_tier:** `tier3-emergent` ·
  **basis:** `physical-picture-expectation` · **delta:** YES.
- **⚑ Substrate conflict carried (same as R03's):** `parameter_register.md:302` labels R35 `DERIVED`;
  `:357` labels it "DERIVED-in-form but un-executed". Cross-referenced to the R03/R04/R05 conflict
  entries; ⛔ not counted twice in the conflict set.
- **LIVE.**

### R27 — the `K₂` density integral `K₂ = ∫ (T_w β₂'² + K_η β₂² + λ·T_Ω β₂²) dV` (route 2)

- **Occurrences (2):** `QID:K2_density_integral_premise@PY#340` (terms `PY:337-339`),
  `QID:K2_density_integral_premise@WL#236` (terms `WL:233-235`).
- All fields as **R26**, witness `RES:016:K2-FORM`, `premise-dependence`.
- ⭐ Recorded on the row: this expression is the **only** statement in either engine of *where the K₂
  form comes from*, and it is consumed by nothing but the dimension walk — which is precisely the H8
  finding (`…ledger_stage016_l2_so3_covariance.md:259-270`) that dropping `λ_m` from the angular term
  leaves the audit green.
- **is_tier:** `no-tier:unadjudicated` · **should_be_tier:** `tier3-emergent` ·
  **basis:** `physical-picture-expectation` · **delta:** YES · **LIVE.**

### R28 — the frozen-radial-calibration premise (route 2)

- **Occurrences (2):** `QID:frozen_radial_calibration_premise@PY#828`,
  `QID:frozen_radial_calibration_premise@WL#639`.
  Statement: "angular structure is earned; **radial profile/scalars are frozen calibration inputs**, so
  the joint is CALIBRATED not PASS."
- **Witness:** both results, `premise-dependence` — it is the scope the result set carries with
  `RES:016:K2-FORM` (`REPORTED_RESULTS.md:96-99`, quoting the stage's §3 at
  `…ledger_stage016_l2_so3_covariance.md:311-314`).
- **Route-2 bounds:** (i) depended on by `RES:016:K2-FORM`; (ii) not derived elsewhere — it is a status
  premise, and the register records the *route out of it* as PENDING, not as executed.
- **Axis A:** `A-REDUCIBLE-UNDERIVED`. §9 — **named route**: R36, recorded at
  `parameter_register.md:303`: "a Gate-1 `R0`-support-equation derivation of the ℓ=2 frozen calibration
  `{β₂(w), T_Ω}` from the straight-reference throat `R0(w)` … `PENDING`". **Executable within this
  framework:** the same locus scopes R36 explicitly to "the ℓ=2 support-equation level (the frozen wall
  profile + radial/support scalars, **not the nonlinear throat**)", and the Gate-1 `R0` machinery exists
  and is cited as consumed by this stage (`…ledger_stage016_l2_so3_covariance.md:326`).
- **Axis B:** `B-POSTULATED` — asserted as a premise; ⛔ no code locus is owed (§3.2). §9 evidence
  (where the model posits it): `PY:828`, `WL:639`, `…ledger_stage016_l2_so3_covariance.md:311-314`.
- **Axis C:** §3.3.1(4) — the artifact itself declares it ⇒ leaf `C-SELF` ⇒ `PHYSICS-FED = false`.
- **CONVENTION-LADEN:** `false` (not a candidate).
- **is_tier:** `tier1-debt` · **should_be_tier:** `tier3-emergent` · **basis:** `named-route`
  (R36, `parameter_register.md:303`) · **delta:** YES.
- ⛔ In **none** of the three §4 near-misses (§4: a premise is not a derivation that fell short).
- ⚠ Gap **S5**: the binding site is a `print` / `Print`. §7.2 does not say whether a prose emission that
  states a model premise is a binding site; it was recorded as one.
- **LIVE.**

---

## 3. THE REQUIRED OUTPUT SET (§5.8) — stage016, both engines

⛔ These do not partition anything and ⛔ must not be summed (§5.8).
⛔ Every count is an `is_tier` count (§6.1); `should_be_tier` appears in **no** numerator or denominator.

| # | output | occurrence level | QID level |
|---|---|---:|---:|
| 1 | **TIER 1, as a RANGE** (⛔ never a scalar) | **[10, 39]** | **[5, 20]** |
| 1a | — `tier1-debt` | 10 | 5 |
| 1b | — `tier1-structural` | 0 | 0 |
| 1c | — `tier1-postulate` | 0 | 0 |
| 2 | **TIER 2 — calibrated** (§10.2.1: counted regardless of any other reduction) | 0 | 0 |
| 3 | **TIER 3**, split by axis B **and** by calibration-propagation | 0 | 0 |
| 3a | — by axis B: `B-EXECUTED` / `B-DERIVED-IN-FORM-UNEXECUTED` | 0 / 0 | 0 / 0 |
| 3b | — `tier3-calibration-propagated` / `tier3-held-out` | 0 / 0 | 0 / 0 |
| 4 | **`DERIVED`** (`B-EXECUTED ∧ PHYSICS-FED ∧ A-REDUCED`) | **0** | **0** |
| 5 | near-miss **`executed-but-not-physics-fed`** | 11 | 6 |
| 6 | near-miss **`derived-in-form-but-unexecuted`** | 0 | 0 |
| 7 | near-miss **`physics-fed-but-declared-literal`** | 0 | 0 |
| 8 | **`unclassified-nonfed`** | 13 | 6 |
| 9 | **convention** bucket (flag `true`) | 0 | 0 |
| 10 | **`unadjudicated`** | 29 | 15 |
| 11 | **`convention-unadjudicated`** | 12 | 6 |
| 12 | **`self-referential`** | 0 | 0 |
| 13 | **conflict set** (all intra-occurrence) | 6 | 3 |
| 14 | **stage043 delta** — see §5 | — | — |
| 15 | **out-of-scope list** — see §4 | 197 objects | — |
| 16 | **`reached-by-no-reported-result`** | 136 objects | — |
| 17 | **reported-result set** — `REPORTED_RESULTS.md:64-106`, 2 results with loci | 2 | — |
| 18 | **`unvalued-in-universe`** (Route C; QID-level only) | 0 by construction | 0 |
| 19 | **`C-PEER` populations, reported SEPARATELY** — `peer-cited-in-artifact` / `peer-cited-in-stage-note` | **0 / 0** | 0 / 0 |
| 20 | **`no-tier:independent-variable`** | 4 | 2 |

**Occurrence totals check:** `tier1-debt` 10 + `unadjudicated` 29 + `unclassified-nonfed` 13 +
`independent-variable` 4 = **56** ✓ (the four `is_tier` values are disjoint; buckets 5, 11, 13 overlap
them by construction).

**Fractions (⛔ never bare — both denominators, §10.1):**
`DERIVED` = **0 / 56 occurrences** and **0 / 27 QIDs**.
⚠ The rank denominator of §10.1 is `OPEN-METHOD` and is not computed here.

**Why tier 3 and `DERIVED` are empty, stated as a mechanism not a verdict:** every closure in this stage
terminates in `C-SELF` literals or `C-FREE` symbols. The one candidate physics leaf — the `−Δ_S²`
operator — contributes none under §3.3.1(1) because no locus of a model equation is cited (R13).

**Tier-1 range, both levels, per §5.6:**
- occurrence level: lower `10` = the evidenced `tier1-debt` rows — R03, R04, R05 (the three radial
  scalars, 2 engines each = 6) + R20 (2) + R28 (2); upper `39` = `10 + 29 unadjudicated`.
- QID level with §10.2.2 rule 4 applied: lower `5` = {M̃, K̃, T̃_Ω, K₂-assembly-form,
  frozen-radial-premise}; upper `20` = `5 + 15` QIDs that have **no** tiered occurrence and ≥1
  `no-tier:unadjudicated` occurrence. `QID:l2_channel_set` widens the span once (rule 4/5), not twice.

**⚠ Reporting obligation discharged, not interpreted:** the span is wide at both levels (29 of 56
occurrences; 15 of 27 QIDs). Its causes are itemised in §6.

---

## 4. OUT-OF-SCOPE LIST (§7.3, §5.8 output 15) — 197 objects, each with its reason

Counted as *listed objects* (§0 rule 3). ⛔ Nothing is silently dropped; every group carries its loci.

### 4.1 `PY` — 101 objects

| loci | objects | reason |
|---|---:|---|
| `PY:27-28` (`PASS_COUNT`, `FAIL_COUNT`) | 2 | test scaffolding |
| `PY:30-34` (five verdict-token strings) | 5 | display and formatting constants |
| `PY:158-175` (`DIMENSION_BASIS`, `Dim`, `ZERO_DIM`, `EXPECTED_M/K/RATIO`, `DIMENSIONLESS_FUNCTIONS`) | 7 | `reached-by-no-reported-result` |
| `PY:223-227` (`a_dim, dw_dim, dOmega_dim, beta2_dim, beta2_prime_dim, mu_eta_density, T_w_density, K_eta_density, T_Omega_density`) | 9 | `reached-by-no-reported-result` ⭐ see gap **S8** |
| `PY:246` (`expr_hash`) | 1 | test scaffolding |
| `PY:250` (`scoped_verdict`) | 1 | `reached-by-no-reported-result` (gate machinery) |
| `PY:295`, `PY:300` (`gram_is_identity`, `residuals_zero`) | 2 | `reached-by-no-reported-result` (gate flags) ⚠ see gap **S9** |
| `PY:308` (`extract_k2_coeff` rule) | 1 | `reached-by-no-reported-result` (its *result* is in-universe at `PY:500`) |
| `PY:314-325` (the twelve `dim_rules.*` literals) | 12 | `reached-by-no-reported-result` (dimension exponent vectors, §7.1.1) |
| `PY:342-350` (the nine walked dimensions) | 9 | `reached-by-no-reported-result` (dimension exponent vectors) |
| `PY:352` (`ok`) | 1 | `reached-by-no-reported-result` (gate flag) |
| `PY:406-422` (`corrupt_rules_for` + the four corrupt maps) | 1 | controls and deliberate negatives |
| `PY:427`, `:431`, `:441` (baseline, `density_corruptions`, `FAIL_DIMENSIONAL_probe`) | 3 | controls and deliberate negatives |
| `PY:466-488` (the 21 emitted `dimension_records`) | 21 | `reached-by-no-reported-result` (dimension exponent vectors) |
| `PY:504` (`k2_coeff_residuals_zero`) | 1 | `reached-by-no-reported-result` (gate flag) |
| `PY:506`, `PY:508` (`lambda_ref`, `k2_ref`) | 2 | `reached-by-no-reported-result` ⭐ see gap **S7** — second binding sites of in-universe QIDs, consumed only by the dimension block |
| `PY:509` (`dimensional`) | 1 | `reached-by-no-reported-result` |
| `PY:511-512` (`input_hashes`, `distinct_hashes`) | 2 | test scaffolding (tautology guard) |
| `PY:513-514` (`self_overlaps`, `tautology_clear`) | 2 | `reached-by-no-reported-result` (duplicates the in-universe Gram diagonal; feeds only the gate) |
| `PY:516-627` (probe builders, four probe instances, `probes`, the `able_to_fail` block) | 15 | controls and deliberate negatives |
| `PY:628`, `:634`, `:640` (`covariant_ok`, `gates`, `verdict`) | 3 | `reached-by-no-reported-result` (gate flags) |
| `PY:665-884` (the assertion + emission block) | 0 new value bindings | 82 recorded checks and ~60 prints; the only occurrence admitted from it is `PY:679` (R17). The rest assert out-of-universe quantities (dimension records, probe verdicts, gate flags) or re-assert already-counted values. |

### 4.2 `WL` — 96 objects

| loci | objects | reason |
|---|---:|---|
| `WL:12-13` (`passCount`, `failCount`) | 2 | test scaffolding |
| `WL:15-19` (five verdict tokens) | 5 | display and formatting constants |
| `WL:23-26` (`zeroDim`, `expectedM/K/Ratio`) | 4 | `reached-by-no-reported-result` |
| `WL:28` (`failureMessage`) | 1 | test scaffolding |
| `WL:30-105` (raise/heading/clean/fmt/assert/expect*/residual/verdict helpers) | 16 | `reached-by-no-reported-result` (assert + verdict machinery) |
| `WL:121-162` (`dimensionAxisSlots` … `dimOf`) | 6 | `reached-by-no-reported-result` (dimension machinery) |
| `WL:172` (`ys`) | 1 | `reached-by-no-reported-result` (a re-listing of the already-counted harmonics) |
| `WL:212`, `WL:214` (`gramIsIdentity`, `residualsZero`) | 2 | `reached-by-no-reported-result` (gate flags) |
| `WL:217` (`extractK2Coeff`) | 1 | `reached-by-no-reported-result` |
| `WL:225` (`k2CoeffResidualsZero`) | 1 | `reached-by-no-reported-result` (gate flag) |
| `WL:227`, `WL:229` (`lambdaRef`, `k2Ref`) | 2 | `reached-by-no-reported-result` ⭐ gap **S7** |
| `WL:238-251` (the twelve `makeDimRules` entries) | 12 | `reached-by-no-reported-result` (dimension exponent vectors) |
| `WL:253-298` (nine walked dimensions + `Ok`) | 10 | `reached-by-no-reported-result` |
| `WL:300` (`corruptRulesFor`) | 1 | controls and deliberate negatives |
| `WL:315`, `:316`, `:322` (`dimRules`, `baselineDim`, `dimensionalOk`) | 3 | `reached-by-no-reported-result` |
| `WL:318`, `:323` (`densityCorruptions`, `dimProbe`) | 2 | controls and deliberate negatives |
| `WL:342`, `:343`, `:344`, `:346` (hash guard) | 4 | test scaffolding |
| `WL:345` (`selfOverlaps`) | 1 | `reached-by-no-reported-result` |
| `WL:348-439` (probe builders, instances, `probes`, expected-verdict maps, able-to-fail block) | 17 | controls and deliberate negatives |
| `WL:441`, `:442`, `:448` (`covariantOk`, `gateBooleans`, `verdict`) | 3 | `reached-by-no-reported-result` (gate flags) |
| `WL:450-451` (`matrixSquaredResidual`, `symbolNames`) | 2 | `reached-by-no-reported-result` |
| `WL:455-702` (run-functions, arity self-check, tallies) | 0 new value bindings | 91 recorded checks and the print block; only `WL:466` (R17) is admitted. `WL:614-628`'s arity self-check is test scaffolding. |

### 4.3 Out-of-scope totals by reason

| reason | objects |
|---|---:|
| `reached-by-no-reported-result` | **136** |
| controls and deliberate negatives | 39 |
| test scaffolding | 12 |
| display and formatting constants | 10 |
| `retired-or-excluded` | 0 |
| **`paused-or-pending` (provisional sub-count, re-checked every pass)** | **0** |
| `derived-proposition-not-constitutive` | 0 |
| **total** | **197** |

⛔ **The 136 is a finding, not housekeeping** (§7.1.1): most of both engines — the entire dimension
block, the whole probe battery, every gate flag — lies outside the closure of the stage's own two
reported results.

---

## 5. §8.4 reconciliation with stage043 — partial, and declared partial

⛔ Full reconciliation of all 152 `REG:` IDs is corpus-level and is **not** attempted here. What this
stage can settle, with loci opened:

| census QID | `REG:` ID(s) | note |
|---|---|---|
| `QID:M_tilde` | `REG:C1:Mtilde` (`…stage043…py:371`, `CAT_OPEN`) **and** `REG:derived:Mtilde_definition` (`:298`, `CAT_DERIVED`) | one quantity, **two** IDs in two categories → conflict entry C1 |
| `QID:K_tilde` | `REG:C1:Ktilde` (`:371`) **and** `REG:derived:Ktilde_definition` (`:299`) | same shape → C2 |
| `QID:T_Omega_tilde` | `REG:C1:Ttilde_Omega` (`:371`) **and** `REG:derived:Ttilde_definition` (`:299`) | same shape → C3 |
| `QID:lambda_m` | `REG:derived:lambda_m_SO3` (`:301`, `CAT_DERIVED`) | substrate label recorded, ⛔ not inherited (§2) |
| the other **23** census QIDs | none | **census extension** — §5.8 output 14, direction 1 |

Directions 2 and 3 of output 14 (`REG:` IDs reconciled out-of-scope; Route-C `unvalued-in-universe`)
are **not** answerable from one stage and are left open; ⛔ recorded as owed, not as zero.
For this stage's QIDs, Route C is empty: `M̃/K̃/T̃_Ω` are valued at
`scripts/ledger_stage017_grouped_p2_lane_isotropy_sympy_audit.py:41-43` and
`mathematica/ledger_stage017_grouped_p2_lane_isotropy_mathematica_audit.wl:273-275`.

---

## 6. CONFLICT SET (§10.3) — 6 occurrences / 3 QIDs, all **intra-occurrence**

| id | QID | occurrences | the two (or more) incompatible substrate claims, with loci |
|---|---|---|---|
| C1 | `QID:M_tilde` | `@PY#219`, `@WL#21` | **`DERIVED`** at `notes/parameter_register.md:184` and `:302` (R35) · **frozen calibration input / moment-integral never evaluated** at `notes/parameter_register.md:353` (R87), `:357` (R91), `scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:384`,`:390`, and this stage's own `notes/stages/ledger_stage016_l2_so3_covariance.md:311-314` / `PY:828` / `WL:639` |
| C2 | `QID:K_tilde` | `@PY#220`, `@WL#21` | as C1 |
| C3 | `QID:T_Omega_tilde` | `@PY#221`, `@WL#21` | as C1 |

⛔ Not resolved by picking a claim (§10.3); ⛔ no axis was demoted on account of the substrate's
disagreement. R26/R27 carry the same disagreement about the R35 *relation* and are cross-referenced,
⛔ not counted again.

**Cross-occurrence conflicts: 0.** The only cross-engine divergence is `QID:l2_channel_set`
(`no-tier:unclassified-nonfed` in `PY`, `no-tier:unadjudicated` in `WL`), which §10.2.2 rule 5
explicitly rules **not** a conflict; it is reported in both buckets and flagged `mixed-adjudication`.

---

## 7. WHERE THE SCHEMA COULD NOT DECIDE — 10 items, none of them silent

⛔ **Zero rows are blocked**: all 56 carry all three axes, the flag, both tiers and their evidence. What
follows are rules the schema does not supply; each records the reading adopted and its sensitivity.

- **S1 — no axis-C leaf tag for a coordinate / bound integration variable.** §3.1 gives axis A
  `A-INDEPENDENT-VARIABLE`; §3.3's leaf table has no counterpart, and `C-FREE` ("nothing it stands for
  has been fixed") is substantively wrong for θ, φ. **Adopted:** a bound integration/differentiation
  variable is **not a leaf** of the closure it is quantified away in; the coordinate's own declaration
  row takes `C-FREE`. **Sensitivity:** if coordinates were leaves everywhere they appear, all 13
  `no-tier:unclassified-nonfed` rows become `C-UNRESOLVED` ⇒ `no-tier:unadjudicated`; occurrence
  TIER 1 goes **[10, 39] → [10, 52]** and near-miss 1 goes **11 → 0**. This single missing rule moves
  more of this stage's result than any other item here.
- **S2 — no axis-A value for a quantity fixed by mathematics.** The five harmonics and the channel count
  are mathematical objects typed in; no route is recorded, no framework property forecloses one, and
  they are not postulates about the medium ⇒ they fall to `A-UNADJUDICATED` and inflate the tier-1
  **upper span** with 11 occurrences carrying no physical input. Axis C has `C-MATH`; axis A has no
  counterpart.
- **S3 — no rule for an input the producing expression consumes but that does not survive into the
  value.** `sp.diff(K̃ + λT̃_Ω, T̃_Ω)` annihilates the `C-FREE` leaf K̃. **Adopted:** the letter of §3.3
  (input closure) ⇒ `C-UNRESOLVED`. **Sensitivity:** 4 occurrences (R22, R23) would move from
  `no-tier:unadjudicated` to `no-tier:unclassified-nonfed`, narrowing occurrence TIER 1 to `[10, 35]`.
- **S4 — route 1 says "a named quantity to which the ledger **assigns a value**", but §3.2's
  `B-DECLARED-UNASSIGNED` and §3.3's `C-FREE` require rows for symbols with no value.** Admission was
  taken from §7.1.1 (closure membership), which is stated to be the only test. 10 occurrences
  (θ, φ, M̃, K̃, T̃_Ω across both engines) depend on this reading.
- **S5 — §7.2 does not say whether a prose `print`/`Print` that states a model premise is a binding
  site.** Route-2 propositions frequently have no other site. **Adopted:** it is one. 2 occurrences (R28).
- **S6 — route 1 vs route 2 for a function definition that types an expression *form*** (`build_K2`,
  `integrate_s2`, `laplacian_s2`): axis B differs (`B-DECLARED-LITERAL` vs `B-POSTULATED`). **Adopted:**
  `B-DECLARED-LITERAL` where the artifact types an expression, `B-POSTULATED` where it states a premise
  in prose. 6 occurrences; ⛔ **no §4 bucket changes either way**, since none is `PHYSICS-FED`.
- **S7 — the universe is occurrence-level by closure, but aggregation (§10.2) is QID-level.** A second
  binding site of an in-universe QID that the walk does not reach (`PY:506`, `PY:508`, `WL:227`,
  `WL:229`) leaves the occurrence count with no rule saying so. **Adopted:** out-of-scope,
  `reached-by-no-reported-result`, recorded. 4 objects.
- **S8 — §3.3.1(4) terminates a proposition's closure at the positing locus**, so the nine
  density/measure symbols named inside R26/R27's expressions never become leaves. **Adopted:**
  terminate. **Sensitivity:** the alternative admits **18** more occurrences (9 per engine:
  `μ_η, T_w, K_η, T_Ω, β₂, β₂', a, dw, dΩ`), each of which would land in tier 1 or its upper span.
- **S9 — the fixed result set admits and excludes objects of the same kind asymmetrically.**
  `REPORTED_RESULTS.md:73` cites `PY:299` (`lambda_all_six`) as a result locus but not its sibling
  `PY:300` (`residuals_zero`); `:91` does not cite `PY:507` / `WL:228` although the claim states
  `M₂ = M̃`. ⛔ The result set was **not** amended (§7.1.1). Effect: `PY:300`/`WL:214` are out-of-scope
  while `PY:299`/`WL:213` are in; `M₂ = M̃` was admitted on the claim's content with the gap recorded.
- **S10 — §3.4 clause (b) is arguably *satisfiable* for the harmonic normalisation while clause (a) is
  undocumented.** `λ_m` is a model observable and is invariant under rescaling or SO(3)-rotating the
  basis, so (b) has a real test; but no transformation group is documented for the occurrence, and §3.4
  forbids convention-by-intuition ⇒ `UNADJUDICATED`. 10 occurrences whose flag turns on the absence of
  documentation rather than on the physics. Had the flag been set `true`, those 10 would leave every
  tier total for `no-tier:convention` — the flattering direction, which §1 forbids taking on an
  undemonstrated claim.

---

## 8. COVERAGE STATEMENT

**Classified:** the **full** binding-site universe of both engines. Every value-binding site in
`PY` (906 lines) and `WL` (702 lines) was enumerated and assigned exactly one disposition:
**56 in-universe occurrences** (§2) + **197 out-of-scope objects** (§4) = **253 objects**.
⛔ Nothing was truncated, and no group was summarised without a count.

**Not covered, with its exact count — the stage note as an artifact.** §7.2 makes a stage note an
artifact; this assignment scoped the census to the two engines. The note's in-universe binding sites are
therefore **uncensused, and are enumerated here rather than dropped — 14 occurrences**, all in
`/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage016_l2_so3_covariance.md`:

1–5. the five harmonics, `:29-35` · 6. the S² measure, `:36` · 7. `Gram = I₅`, `:38` ·
8. the `−Δ_S²` operator, `:45-47` · 9. `λ_m = 6` (Rayleigh), `:49` · 10. the m-degeneracy, `:52-53` ·
11. `K₂ = K̃ + λ_m·T̃_Ω`, `:58` · 12. `M₂ = M̃`, `:58` · 13. `M₂ = μ_η β₂² dV`, `:78` ·
14. `K₂ = (T_w β₂'² + K_η β₂² + λ·T_Ω β₂²) dV`, `:79`.

The note's remaining value-bearing rows are dimension exponent vectors and coverage/enumeration tables
(`:73-83`, `:85-149`, `:175-207`) — out of universe by the same rule that excludes the engines'
dimension blocks. ⚠ The note's own citations `sympy:355-366` (`:194`) and `sympy:364` (`:210`) are
**stale** (the rules are at `PY:314-325`); recorded as an attribution defect on R03.

**Other stages, and cross-stage identity:** out of scope for this pass. ⛔ No merge was adjudicated
(§8.3) — in particular `T̃_Ω`(016) vs `T_Ω`(023) is left as **two** QIDs and flagged as an open identity
question, per the stage's own hazard note (`…ledger_stage016_l2_so3_covariance.md:187-192`).
