# stage043 — census reconciliation frame

**Status:** the §8.4 reconciliation frame for stage043, taken as the census pilot. It answers two
questions: whether stage043 contributes any census occurrence, and which §8.4 route each of the
stage's ten `ratified_category` blocks takes.

⛔ **This file classifies nothing.** No census row, no tier, no axis value, no `is_tier` /
`should_be_tier`, no derived/declared judgement, no count of anything as derived or declared. Every
number below is either a count of manifest IDs or a count of routes.

**Method.** Every locus cited here was opened and read (`CENSUS_SCHEMA.md` §9.1 rule 1). Loci are full
paths (§7.2). The manifest total was measured from the file, not taken from any document's assertion —
see §3.

**Sources of record used.**

- `/var/projects/toy_physics/research/pde_ledger_v2/notes/census/CENSUS_SCHEMA.md`
- `/var/projects/toy_physics/research/pde_ledger_v2/notes/census/REPORTED_RESULTS.md`
- `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage043_irreducible_count_range_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage043_irreducible_count_range.md`
- `/var/projects/toy_physics/research/pde_ledger_v2/notes/parameter_register.md`
- `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`

---

## 1. TASK 1 — does stage043 contribute any census occurrence?

> ⭐ **No. Stage043 contributes ZERO census occurrences.** Every binding site in
> `ledger_stage043_irreducible_count_range_sympy_audit.py` and
> `notes/stages/ledger_stage043_irreducible_count_range.md` is out of the census universe.

Two independent rules reach that answer. The first is the one that decides the manifest rows the task
asks about; the second decides the stage as a whole.

### 1.1 The deciding rule for a `REGISTER_TO_COUNT_MANIFEST` row

> **§7.2 read against §7.1: a manifest row asserts a CATEGORY for a knob, and a category is not a
> value of that knob. So the row is not a binding site at which the artifact "assigns or asserts a
> value", and it is not an occurrence.**

§7.2 defines an occurrence as one `(artifact, binding-site, quantity)` triple "at which an artifact
assigns or asserts a value" (`CENSUS_SCHEMA.md:680-681`). §7.1 fixes what "a value" means for the
objects the census ranges over: "any named quantity to which the ledger **assigns a value — numeric or
closed-form symbolic**" (`:583-584`).

A manifest row is a `RegisterFact` — the dataclass at
`scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:179-191`, instantiated at `:213-408`.
Take the first one, `:216-218`:

```python
fact_group(
    ("REG:a:hbar", "REG:a:m_GNLS", "REG:a:K_EOS", "REG:a:rho0"),
    "ACTION", CAT_CONTINUOUS, part="I-II", bucket="a",
)
```

For the quantity `ħ` this row asserts `source_status = "ACTION"` and
`ratified_category = "continuous-counted-irreducible"` (the constant at `:137`). Neither is a numeric
value of `ħ` nor a closed-form symbolic expression for `ħ`. What the row asserts is a **provenance
label about how some other artifact treats `ħ`** — a statement about the ledger's bookkeeping, not an
assignment to the quantity. ⇒ no occurrence, for any of the 152 rows.

⚠ **Reinforcement, not the deciding rule.** Even on a reading that admitted the row, §2 bars what it
carries: "⛔ **No census row inherits a substrate class**", naming "stage043's `ratified_category`
constants" as one of the three fused taxonomies (`CENSUS_SCHEMA.md:46-54`). The thing a manifest row
asserts is precisely the thing §2 forbids entering a census row. The category is evidence *about* the
substrate; §8.4 tells the census to **reconcile against** it, never to ingest it (`:802`).

### 1.2 The second rule, which decides the whole stage

§7.1.1: "A named quantity is IN THE CENSUS UNIVERSE **iff** it lies in the transitive input closure
(§3.3) of at least one REPORTED PHYSICAL RESULT of the ledger" (`CENSUS_SCHEMA.md:617-618`).
`REPORTED_RESULTS.md:225-228` records stage043's reported-physical-result count as **ZERO**.

⚠ That alone does not settle it — §7.1.1 is corpus-wide, and "⛔ Membership is NOT decidable by reading
one script in isolation" (`:664-668`). A stage043 binding site would be in universe if some *other*
stage's reported result reached it. It does not, and the stage says why: it "mints NO reduction/codim
edge", "FEEDS Part VII's unified count" (note `:315-319`), and its consumers are stages 044/045/046/047/049
(note `:392-397`). Nothing physical is computed **from** stage043's outputs; it is terminally
downstream of the register it counts. Under the fixed result set (`REPORTED_RESULTS.md`, five results,
all at stages 016 and 023), no closure reaches this stage.

⇒ Every stage043 binding site is recorded out of scope with the reason **`reached-by-no-reported-result`**
(§7.1.1 `:670-672`; §7.3 `:705-707`; §5.8 output 16 `:510`).

⚠ **Scope, stated honestly.** The result set is enumerated for 3 of 43 stages. If a corpus-wide pass
admits a result whose closure reaches stage043, this must be revisited — but such a result would have
to be a physical claim *downstream of the count*, which the stage's own scope class forecloses:
"Scope class: **STRUCTURAL / dimensionless** … a **bookkeeping object**" (note `:30-34`).

### 1.3 The other candidate binding sites in stage043, checked rather than assumed

Two families of binding site in the script are value-shaped and were opened individually.

1. ⭐ **The `Δr` witness points — the only genuine value assignments to model quantities in the
   stage.** The `AlgebraBlock` witness maps at `script:617-622` and `:638-643` assign literals to named
   model symbols: `hbar: 5`, `mass: 5`, `big_k: 1`, `rho0: 1`, `cs0: 1`, `xi_h: sqrt(2)`,
   `h0: Rational(5,4)`, `scale_a: 1`, and in the wall block `delta: sqrt(2)/2`,
   `sigma_wall: sqrt(2)/6`. These are **not** claims that `ħ = 5`. The script's own docstring names
   their role: "exact positive smooth witnesses **guard the real locus**" (`:9-12`) — they exist so the
   Groebner/Krull computation is taken over a nonempty positive real locus. That is **test
   scaffolding** and a **control**, both named in §7.3's exclusion list (`CENSUS_SCHEMA.md:708-709`),
   and the rows are out of universe by §1.2 regardless.
   ⚠ The relations they sit under (`mass*cs0**2 - 5*big_k*rho0**4`, `scale_a*mass*cs0 - hbar`,
   `xi_h**2*mass**2*cs0**2 - 2*hbar**2`, `4*h0 - mass*cs0**2`, `lambda_gamma*cs0 - c_gamma` at
   `:610-614`; `k_eta - t_w*beta**2`, `2*a_b*delta**2 - kappa_b`, `36*sigma_wall**2 - 2*a_b*kappa_b`
   at `:633-635`) are **transcribed** edges earned in other stages, carried with no in-script source
   locus. Where those relations are census occurrences, they are occurrences **of the stages that
   perform them**, not of stage043.
2. **The count integers** — `[40, 49]`, `[27, 36]`, `13`, `9`, `3`, `6`, `11`, `34`, `43`, the ten
   `EXPECTED_CATEGORY_COUNTS` (`:165-176`) and the `Δr` payloads `(10,5,5)` / `(10,7,3)` (`:646-653`).
   None is a value of a model quantity; all are tallies over the ledger's own rows, which §7.1.1
   excludes by name — "coverage figures and row counts, … category tallies" (`:638-640`).

---

## 2. TASK 2 — the §8.4 route table

§8.4 (`CENSUS_SCHEMA.md:789-806`) requires every `REG:` ID to be reconciled in exactly one of two ways:

- **route A** — "it **maps to ≥ 1 census occurrence**";
- **route B** — "it is **recorded as out of scope**, carrying **its stage043 category** … **and a
  reason** referred to §7.1.1 or §7.3".

⚠ Because §1 establishes that stage043 contributes no occurrence of its own, **every route-A mapping
lands on an occurrence elsewhere in the corpus.** Route A is therefore always an *expectation* here,
discharged by the corpus-wide pass, never a fact this pilot can close. Route B assignments are
decidable now, because they turn on what the ID *is*, not on who consumes it.

### 2.1 The table

| # | `ratified_category` (script `:137-146`) | n | route | what decides it |
|---|---|---:|---|---|
| 1 | `continuous-counted-irreducible` | 22 | **A** (expected, conditional) | the model's counted knobs. Route A holds iff each is (i) reached by a result closure and (ii) actually assigned a value somewhere — see ⭐F1 |
| 2 | `reduction-debt-counted-once` | 18 | **A** (expected) — ⭐ **measured counterexample** | in-universe knobs, but ⭐F1 fires here with a worked case (`REG:b:K0c`) |
| 3 | `derived-not-counted` | 40 | **A** — firmest | each is a value a stage computes in closed form; anchor verified below |
| 4 | `calibrated-external-bridge-not-counted` | 10 | **A** | a calibration stays counted (§5.2 `:351-368`, §10.2.1 `:939-953`) and is a census row; the label "not counted" is stage043's rule, not the census's |
| 5 | `departure-no-knob` | 4 | **B** | characterized departures assign no value (note `:80`) — ⭐F4 on the *reason* |
| 6 | `convention` | 2 | **A** | a convention-laden occurrence is a census row in its own reported bucket (§5.2, §5.8 output 9) — ⭐F6 |
| 7 | `discrete-structural-postulate` | 11 | ⭐⭐ **SPLIT — genuinely unclear** | members differ in kind; see ⭐F2 |
| 8 | `extension-convention-open` | 9 | **A** (C1 firm, C2 conditional) | C1 anchors verified below; C2 subject to ⭐F1 |
| 9 | `parallel-track-out-of-scope` | 8 | **B** | register class `GAP` = a missing derivation, not a knob — ⭐F5 on the *reason* |
| 10 | `structural-no-knob` | 28 | **B** (26) / **unclear** (2) | register edges and controls; ⭐F3 and ⭐F4 |

**Route totals.**

| route | categories | IDs |
|---|---|---:|
| A (expected to map to ≥1 census occurrence elsewhere) | 1, 2, 3, 4, 6, 8 | **101** |
| B (recorded out of scope) | 5, 9, 10 | **40** |
| unclear as a block — routes per member, not per category | 7 | **11** |
| **total** | | **152** |

`22 + 18 + 40 + 10 + 2 + 9 = 101`; `4 + 8 + 28 = 40`; `101 + 40 + 11 = 152`.

⚠ Within route B, 2 of category 10's 28 are individually unclear (⭐F3); they are left inside the
route-B count and flagged rather than moved, because the block's default is B.

### 2.2 Route-A anchors that were opened and verified

Route A is an expectation, but three of its blocks have anchors that can be checked now, and were:

- **category 3** — `REG:derived:lambda_m_SO3` is stage016's `λ_m`, computed as a Rayleigh quotient at
  `scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:287` and asserted across all five channels
  at `:299` (`lambda_all_six`). That is inside the closure of `RES:016:L2-IRREP`
  (`REPORTED_RESULTS.md:64-76`). A value is computed at a binding site ⇒ route A is available.
- **category 8, the C1 half** — `REG:C1:Ktilde` and `REG:C1:Ttilde_Omega` are live symbols in
  stage016's `build_K2`, `scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:304-305`
  (`compact(Ktilde + coeff * TomegaTilde)`), inside the closure of `RES:016:K2-FORM`
  (`REPORTED_RESULTS.md:82-99`). The stage043 rows carry `frozen_calibration_input=True` and
  `moment_integral_executed=False` (`script:390-391`) — a frozen calibration input is a value that
  *was* assigned, so route A holds for C1.
- **category 2** — `REG:b:K0c` and `REG:b:T_Omega` enter stage023's rank computation directly:
  `base_constraints = [port["P0_raw"], K0c, Keta + 2 * TOmega]` at
  `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:384`, and `K0c` again in the
  transfer at `:284`/`:302`. Both are in the closure of `RES:023:RETURN-UNDERDETERMINED`
  (`REPORTED_RESULTS.md:132-147`). ⚠ **In universe, and yet no occurrence** — this is ⭐F1.

---

## 3. The manifest total — MEASURED, not taken from any document

⚠ §13 already carries the stage043 note-vs-script disagreement as `OPEN-SUBSTRATE`, with the ruling
that "the census cites the **script** as authoritative for counts" (`CENSUS_SCHEMA.md:1027`). The
measurement below was taken independently of both.

**What was measured.** Every double-quoted identifier literal beginning `REG:` inside the
`MANIFEST_FACTS` definition region, `script:213-408` — the region that ends where `MANIFEST_FACTS` is
last extended (`:396-408`) and before the first consumer function (`category_of`, `:411`). Identifier
literals appearing elsewhere in the file (`:437`, `:444`, `:707`, `:712`, and the mutation blocks at
`:809`–`:1160`) were excluded; they reference manifest rows rather than declaring them.

| measurement | result |
|---|---:|
| identifier literals in `script:213-408` | **152** |
| distinct identifiers among them | **152** (no duplicate) |

**Per-category, measured by prefix and cross-checked against the group declarations:**

| `ratified_category` | measured | `EXPECTED_CATEGORY_COUNTS` (`script:165-176`) | agrees |
|---|---:|---:|---|
| `continuous-counted-irreducible` | 22 | 22 | ✅ |
| `reduction-debt-counted-once` | 18 | 18 | ✅ |
| `derived-not-counted` | 40 | 40 | ✅ |
| `calibrated-external-bridge-not-counted` | 10 | 10 | ✅ |
| `departure-no-knob` | 4 | 4 | ✅ |
| `convention` | 2 | 2 | ✅ |
| `discrete-structural-postulate` | 11 | 11 | ✅ |
| `extension-convention-open` | 9 | 9 | ✅ |
| `parallel-track-out-of-scope` | 8 | 8 | ✅ |
| `structural-no-knob` | 28 | 28 | ✅ |
| **total** | **152** | **152** | ✅ |

Prefix tally as measured: `derived` 40, `struct` 28, `b` 15, `E` 13, `disc` 11, `bridge` 10, `R49` 8,
`c` 7, `C2` 6, `departure` 4, `a` 4, `C1` 3, `conv` 2, `force` 1 — summing to 152. The `E` block splits
`III` 5 / `IV` 7 / `V` 1; the `disc` block splits 7 Parts-I–II / 4 `D1`–`D4`. Mapping those onto the
categories: continuous = `a`(4) + `c`(7) + `force`(1) + `E:III`(5) + the 5 `ACTION` rows of `E:IV`
(`:263-269`) = 22; debt = `b`(15) + `E:IV` `{c_E, Q_E}`(2) + `E:V`(1) = 18; open = `C1`(3) + `C2`(6) = 9.

⇒ **The total measured is 152, and it reconciles against every one of the manifest's own asserted
category counts and against the uniqueness assertion at `script:436`** (`len(set(identifiers)) == 152
and len(identifiers) == 152`).

### 3.1 Two further reconciliations that fell out of the measurement

Recorded because they were measured, not because they were sought.

1. ⭐ **The script's disjoint peers reproduce the note's headline range exactly.** The note presents
   `reduction-debt-counted-once` and `extension-convention-open` as **roles that roll up** into the
   continuous count, not as disjoint peers (note `:212-217`); the script makes them disjoint
   categories. The two presentations agree numerically:
   `continuous(22) + debt(18) = 40` = the note's LOW endpoint, and `+ open(9) = 49` = its HIGH
   endpoint (note `:140-150`). The `OPEN-SUBSTRATE` item at `CENSUS_SCHEMA.md:1027` is a presentation
   difference, not an arithmetic one.
2. ⚠ **The note's own §10 table totals 152 only at the HIGH endpoint.** Reading its
   `continuous counted-irreducible | [40, 49]` row (note `:271`) at HIGH:
   `49 + 11 + 8 + 4 + 2 + 40 + 10 + 28 = 152`. At LOW it gives 143, short by exactly the spread 9 —
   the `C1`/`C2` rows, which the note lists as a sub-row (`:273`) rather than a peer. The script's
   flat 152 is the unambiguous form.

---

## 4. ⭐ Flagged routes — where §8.4 does not decide, recorded rather than forced

These are findings about §8.4 and §7.1, which is what this pilot is for. ⛔ None is resolved here by
picking the reading that reads better.

### ⭐⭐ F1 — an in-universe knob that satisfies NEITHER route (categories 1, 2, 8-C2 — up to 46 IDs)

§8.4 offers exactly two routes, and its own note says "⛔ Only an ID that is **neither** mapped nor
recorded is a defect in the census" (`CENSUS_SCHEMA.md:799`). **A knob can be neither, without any
defect.** Worked case, every locus opened:

- `REG:b:K0c` is declared in stage023 as a bare free symbol:
  `K0c, Keta, TOmega = sp.symbols("K0c K_eta T_Omega", positive=True, real=True)`
  (`scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py:170`). Nothing assigns it a
  value there or anywhere.
- It **is** in the closure of a reported physical result — `base_constraints` at `:384` and the
  transfer `t0 = compact(k0 / (k0 + z0))` at `:302` with `k0 = K0c` at `:284`, feeding
  `RES:023:RETURN-UNDERDETERMINED`. So it is **in universe** by §7.1.1, and route B ("out of scope") is
  false.
- Its register row, `notes/parameter_register.md:170`, carries a dimension (`M T⁻²`), a provenance
  class (`FREE-UNREDUCED`, PENDING scalar-reduction) and prose — **no value**.
- The only other appearance in a source of record is `CANONICAL_DIMENSIONS.md:148`, a dimension
  exponent vector `{0, 1, -2}` — and §7.1 says outright "⛔ **The universe is NOT the dimension-record
  universe**" (`CENSUS_SCHEMA.md:590`).

⇒ `REG:b:K0c` is in the census universe and has **zero** binding sites at which any artifact assigns or
asserts a value. It maps to no occurrence, and it is not out of scope. §8.4's binary is **not
exhaustive**, and the missing third route is *"in universe, value never assigned — a `C-FREE` symbol
(§3.3 `:137`) all the way down"*.

⚠ **Scale.** This is structural, not incidental: `FREE-UNREDUCED` means exactly "no value has been
established". It therefore threatens route A across `reduction-debt-counted-once` (18),
`continuous-counted-irreducible` (22) and the C2 half of `extension-convention-open` (6) — **up to 46
of 152**. How many actually fire is a corpus-wide question this pilot cannot close; one is measured.

⚠ §3.3 already anticipates the *leaf*: "`C-FREE` — a **free symbol carrying no value** … ⛔ it is not
`C-SELF` — nothing was declared" (`:137`, `:150-156`). What has no home is the *ID*: §8.4 has no route
for an object that is a leaf of somebody's closure and the subject of nobody's assignment.

### ⭐⭐ F2 — `discrete-structural-postulate` (11) cannot be routed as a block

The eleven members differ in kind, and §7.1's "numeric or closed-form symbolic" (`:583-585`) splits
them:

| member | what it is | route |
|---|---|---|
| `REG:disc:EOS_exponent_5` | a **number** — the EOS exponent 5, which stage043's own block writes as `mass*cs0**2 - 5*big_k*rho0**4` (`script:610`) | **A** |
| `REG:disc:D1:H_existence_and_Poschl_Teller_law` | an **existence claim + a law-form** (`O_⊥ = A†A ≥ 0`, note `:176-178`) | **split** — a Pöschl–Teller form *is* closed-form symbolic (⇒ A); an existence claim is not (⇒ B) |
| `REG:disc:D2:s_i_Z2_topology` | `s = ±1` orientation topology; note `:180-182` says outright "the `Q_χ = s` VALUE is DERIVED (R57); the postulate is the ±w-orientation topology — **no magnitude**" | **B** — the value lives at a *different* ID (`REG:derived:charge_value_R57`) |
| `REG:disc:D3:R63_mouth_BC_class` | a **class selection** from `{V, M, J, MIXED}` (note `:183-185`); its per-class `{φ, q, j, λ}` fillers are separate continuous rows | **unclear** |
| `REG:disc:D4:R68_tau_d_time_arrow` | a **T-parity basis** / arrow (note `:186-187`) | **unclear** |
| `chi_B_field`, `Gamma_B_law`, `chi_B_gating`, `G0_postulate_1/2/6` | field-existence, law-forms and boundary postulates (note `:113-114`) | **mostly B**, same split as D1 |

⭐⭐ **The finding.** §7.1 admits only values that are "numeric or closed-form symbolic". A **Z₂ label**,
a **BC class**, and a **time-arrow parity** are none of those — yet §3.1 provides
`A-IRREDUCIBLE-POSTULATE` for exactly this kind of object, and §9 requires "where the postulate is
**stated** as a defining property" (`:822`) as its evidence. The schema wants these classified and its
universe rule cannot admit them. ⇒ either §7.1's value clause needs a **discrete/label-valued** case,
or §8.4 needs a third route for a postulate that fixes a *choice* rather than a *magnitude*.

⚠ Stage043 itself records the D1 split as a live sensitivity: "were D1 split into field-existence +
operator-law-form, `D_III_VI → 5`" (note `:312-313`). The seam is already known to the substrate.

### ⭐ F3 — `structural-no-knob`: 2 of 28 unclear, 24 of 28 with no reason token

Measured composition of the 28 (`script:347-366`): 24 register edges (`REG:struct:R23…R84`) plus
`z_g_opaque`, `z_b_opaque`, `control_inventory`, `census_bookkeeping`.

- `REG:struct:control_inventory` → route B with a **clean** reason: §7.3 "controls and deliberate
  negatives" (`:708`).
- `REG:struct:census_bookkeeping` → route B with a **clean** reason: §7.1.1's exclusion of "the
  census's own outputs" (`:640`).
- The 24 R-edges → route B, but see ⭐F4: they are obligations, firewalls and validity records that
  assign no value, and no §7.3 token says that.
- ⭐ `REG:struct:z_g_opaque` / `REG:struct:z_b_opaque` → **unclear**. `notes/parameter_register.md:220`
  records `z_g` as an "**opaque** escape factor; strict positivity `0<z_g≤1` (**witness `z_g*=1`**)
  POSTULATED — 031 verifies NO formula for `z_g`". No formula ⇒ route B; but a *witness value* `z_g* = 1`
  is an assertion of a value at a binding site ⇒ possibly route A, then excluded at row level as a
  control. §8.4 gives no rule for an ID whose only value is a witness.

### ⭐ F4 — `departure-no-knob` (4) and the 24 R-edges: route B is clear, the REASON is not

`REG:departure:R66_native_P_no_emergent_Gauss`, `:R73_bT_T_even`, `:R81_scalar_sim_gated`,
`:light_stray_longitudinal` (`script:325-333`) are **characterized departures** — statements that the
model's behaviour differs from a target. Note `:80`: "A characterized DEPARTURE adds NO knob." They
assign no value, so §7.1's universe never contains them and no closure could reach them.

⛔ **But §8.4 requires the reason to be "referred to §7.1.1 or §7.3", and neither has a token for
"this is not a value-bearing object at all".** §7.1.1 supplies `reached-by-no-reported-result`, which
describes a value-bearing object that nobody's result happens to reach; §7.3's list is controls,
test scaffolding, display constants, loop bounds/tolerances, and retired approaches (`:708-712`). Using
the §7.1.1 token here **misdescribes the row and corrupts a reported count**: §5.8 output 16 makes
`reached-by-no-reported-result` a *finding about the ledger* — "⛔ **If that set is large, that is
itself a census finding**" (`:674-676`) — and padding it with 28 objects that were never candidates
inflates precisely the finding it exists to report.

⇒ §7.3 needs a `not-a-value-bearing-object` reason, or §8.4 needs to say which existing token covers
it.

### ⭐ F5 — `parallel-track-out-of-scope` (8): route B clear, reason token mismatched

The eight `REG:R49:profile_E0…E3` / `profile_B0…B3` (`script:338-346`) are the register's
`GAP_OPEN_FIELD_PROFILE` rows. `notes/parameter_register.md:37` defines class `GAP` as "a *missing*
derivation (**not a knob**) — a deferred obligation | tracked, not counted". A missing derivation
assigns no value ⇒ route B, and ⭐F4's reason gap applies here too.

⚠ Additionally, the nearest §7.3 token — "values belonging to a **retired or excluded** approach"
(`:712`) — is **factually wrong** for these. The track is **paused and pending**, not retired:
`parameter_register.md:887` reads "**Edge R49 (PENDING):** all eight ↦ discharged by the native tilted-sleeve
BVP solve (U2-gated)", and note `:219-226` calls it "the **PAUSED** EM/U1 Phase-C parallel track".
§7.3's own rule bites here: "⚠ an exclusion reason that cannot be checked against the row is not a
reason" (`:716`). Recording these as retired would fail that check.

### ⭐ F6 — `convention` (2): §8.4 does not say whether a quotiented row is "mapped"

`REG:conv:a_pin` and `REG:conv:lambda_gamma_ratio` (`script:334-337`) are value-bearing — a pin
assigns a value. Read as route A: §5.2 puts convention-laden occurrences **outside all three tiers**
but keeps them as census rows in their own reported bucket (§5.2 `:351-353`; §5.8 output 9 `:503`), and
§5.7's enum reserves a value for them, `no-tier:convention` (`:470`). A classified row is a mapped row.

⚠ **But §10.1 removes them from the denominator anyway** — "the convention-laden occurrences of §3.4,
quotiented out" (`:915-916`) — sitting in the same sentence as "the out-of-scope rows of §7.3". Two
different mechanisms with the same effect on the denominator, and §8.4 distinguishes routes by
*scope*, not by *denominator*. The reading taken here is **route A** (the row is classified, not
excluded); recorded so a reviewer can contest it.

⛔ Two traps carried with this block, neither resolved here: **(i)** §2 and §3.4 forbid inheriting the
`CONV` label — the census must run the operational (a)/(b)/(c) test (`:226-240`) and the author's label
is not evidence; **(ii)** §3.4's vacuous-invariance rule means "⛔ **A quantity from which no
dimensionless observable of the model is currently computable is `UNADJUDICATED`. It is NEVER
`CONVENTION`**" (`:246-248`). Whether an `a`-pin clears clause (b) is a classification question, not a
routing one.

### ⚠ F7 — recorded, not a route ambiguity

- **Category 4's name is a trap.** `calibrated-external-bridge-**not-counted**` states stage043's
  counting rule (note `:69`, `count = nominal − DERIVED-and-EXECUTED − CONV − external-bridge`). The
  census counts these: §5.2 and §10.2.1 both quote "**IMPOSED calibrations STAY COUNTED** — a tuning is
  NOT a reduction" from note `:78` (`CENSUS_SCHEMA.md:363-365`, `:950-951`). Route A stands; the label
  must not be inherited (§2).
- **`REG:E:IV:ell_over_a` is the §3.4.1 laundering-hazard row.** `ℓ/a = ε_r = 1/20` is the one tagged
  `[CALIBRATED]` two rows above a `[CONVENTION]` twin in one six-row table
  (`CENSUS_SCHEMA.md:261-272`). Its **route** is A either way — a convention row is still a census row
  (⭐F6) — so this is a classification hazard, not a routing one. Recorded so the classifier meets it
  already flagged.
- **Category 8's name fuses two source statuses.** `extension-convention-open` covers both
  `DERIVED-IN-FORM-UNEXECUTED` (C1) and `DOWNSTREAM-DEFERRED-OPEN` (C2) (`script:159-160`,
  `:384`, `:399`) — an instance of the three-question fusion §2 describes (`:46-54`). Only the C1 half
  corresponds to §4's near-miss 2 (`derived-in-form-but-unexecuted`).

---

## 5. What this frame does NOT deliver

1. ⛔ **No route-A mapping is closed.** Route A names the *expectation* that an ID meets ≥1 occurrence
   elsewhere; discharging it needs the corpus-wide closure pass §7.1.1 requires and §13 carries as
   `OPEN-METHOD` (`CENSUS_SCHEMA.md:1024`). Three anchors were verified (§2.2); 98 of the 101 were not.
2. ⛔ **No classification of any kind appears above** — no census row, no tier, no axis value, no
   derived/declared count.
3. ⚠ **Route B is decidable now and was decided; its *reasons* are not all available** — ⭐F4 and ⭐F5
   are gaps in §7.3's reason list, not judgements withheld.
4. ⚠ **The §8.4 delta this feeds** (§5.8 output 14, `:508`) is the *option-2* half only. The other
   direction — census QIDs mapping to no `REG:` ID — cannot be computed until rows exist.
