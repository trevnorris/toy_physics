# Census rows — stage023, both engines (PILOT)

**Status:** builder deliverable. Rows classified under `CENSUS_SCHEMA.md` against the three reported
physical results of `REPORTED_RESULTS.md` §2.

⛔ **This file is not a verdict on the stage.** It records provenance-axis classifications only.
⛔ **The stage's `24/0/10` and `27/0/7` tallies are DIMENSIONAL-CORRECTNESS verdicts on a different
axis** (schema §11). Their row universes were *not* used as a target and are not reproduced here; the
dimension block is in fact **out of this census's universe** (§F.1) under the closure reading adopted
at §B.2.

⛔ Every locus below was opened and read. No attribution was taken from the stage note, the parameter
register, or `REPORTED_RESULTS.md` without opening the file it points at.

**Artifacts in scope (§7.2 — two artifacts, so a dual-engine quantity has ≥ 2 occurrences):**

- `PY` = `/var/projects/toy_physics/research/pde_ledger_v2/scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py`
- `WL` = `/var/projects/toy_physics/research/pde_ledger_v2/mathematica/ledger_stage023_nullspace_underdetermination_mathematica_audit.wl`

---

## A. The reported-result set (schema §5.8 output 17) — the fixed starting set

Taken **unchanged** from `REPORTED_RESULTS.md` §2.1. ⛔ Not amended after classification began.

| id | claim (abbreviated) | loci opened for this pass |
|---|---|---|
| `RES:023:RETURN-UNDERDETERMINED` | over 11 generator dofs the 3 collected named constraints have rank 3; `{Z0_ret, Z1_ret}` are untouched by all of them yet move `T0`/`T1`; `return_moving_nullity = 2` | `PY:205-217`, `:384`, `:387-393`, `:395-413`, `:807-812`, `:819`; `WL:123-124`, `:143-148`, `:150-151`, `:560-563`; note `:40-61` |
| `RES:023:A0-FORM` | `A0 = i·v₀·(aω/c_s)·M0·(1−T0)` coincides with the pathA_29 form and is realized at `ω¹` | `PY:337`, `:342`, `:345`, `:349`, `:351`, `:882`, `:884`; `WL:197`, `:199`, `:201`, `:203`, `:452`, `:458`, `:638`, `:640`; note `:77-94` |
| `RES:023:A1-FORM` | `A1 = i·v₁·(aω/c_s)³·D1·(1−T1)` coincides with the pathA_29 dipole form and is realized at `ω³` | `PY:338-339`, `:343`, `:346-348`, `:350`, `:352`, `:883`, `:885`; `WL:198`, `:200`, `:202`, `:204`, `:453`, `:459`, `:639`, `:641`; note `:77-94` |

Reachability-witness codes used in the row tables: **`W-RU`**, **`W-A0`**, **`W-A1`** = reached from
the correspondingly-named result; the hop is given per row.

---

## B. Method, and the interpretive decisions this pass had to make

Every item in §B.2–§B.9 is a place the schema **did not decide the case**. Each states the rule I
adopted, and the rule is applied uniformly. These are reported as findings of the pilot, per the
instruction to record blocked/undecided cases rather than resolve them silently.

### B.1 Order of work

Result set (§A, fixed) → closure walk marking every reached binding site with its witness (§C, §D) →
classification. The closure walk and the axis-C determination were the same walk, done once (§7.1.1
step 2 note).

### B.2 ⛔ GAP-1 — *does the closure of a result run through the verdict token?* (largest single effect)

`REPORTED_RESULTS.md` §2.1 lists `PY:1068` (`PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE`)
among `RES:023:RETURN-UNDERDETERMINED`'s loci. The verdict token is produced by `base_verdict`
(`PY:672-686`) from **seven** conditions (`PY:703-717`), of which `dimensional`, `tautological`,
`decoupled`, `overcancel`, `no_consistent_return` and `epsilon_mismatch` are guards. Under strict
data-dependence the token's value depends on all seven, so the **whole dimensional gate, the whole
provenance-class table and the control gains would enter the universe**. The stage's own step-(h)
evidence confirms the dependence is real (note `:898-902`: corrupting one dimension declaration flips
the verdict before `run_dimensional_gate()` is reached).

- **Reading (A), ADOPTED:** the closure is walked from the result's **value-bearing content** — the
  rank/nullity/witness quantities and the `A0`/`A1` forms — not from the gate token that reports it.
  §3.3 defines the closure as "the expression that produced its value", and the guard conditions do
  not produce those values.
- **Reading (B):** the token is the result's locus, so its data-dependence is the closure.
- **Consequence, measured:** reading (B) would add **at minimum 122 further occurrences** — `PY`
  22 `SOURCED_DIMS` (`:460-483`) + 7 `EXPECTED_DIMS` (`:486-494`) + 29 emitted records (`:539-571`),
  `WL` 22 `baseDims` (`:247-255`) + 7 `expectedDims` (`:257-262`) + 29 `printDimRecord` calls
  (`:298-332`) — plus 9 provenance items × 2 class strings per engine (`PY:575-616`, `WL:341-361`)
  and the 4 control symbols per engine. ⛔ **The schema does not decide this and it is the single
  biggest lever on this stage's universe size.**

### B.3 GAP-2 — precedence between `no-tier:unadjudicated` and a determined axis-A value

§5.7 assigns `no-tier:unadjudicated` when `PHYSICS-FED = C-UNRESOLVED`, and §3.3 says "`C-UNRESOLVED`
rows go to `unadjudicated`". But §5.7 *also* assigns `tier1-debt` on `A-REDUCIBLE-UNDERIVED`, and §5
states the tiers are "a projection of axis A **with one stated axis-C guard on tier 3**". A row that
is `A-REDUCIBLE-UNDERIVED` **and** `C-UNRESOLVED` matches both lines; §5.7's precedence note covers
only the convention flag.

**Adopted:** the axis-C `C-UNRESOLVED` demotion applies only where the tier assignment *depended* on
`PHYSICS-FED` — i.e. to `A-REDUCED` rows (tier 3 / `unclassified-nonfed`). Tier 1 and tier 2 are pure
axis-A projections and are not drained by an unresolved closure. This is also what §1's second half
requires: the other reading would move well-evidenced tier-1 rows out of the **lower** bound.
**Effect: 15 occurrences** (all `tier1-debt` rows) sit in the lower bound rather than the span.

### B.4 GAP-3 — `A-REDUCIBLE-UNDERIVED` vs `A-IRREDUCIBLE-STRUCTURAL` when the named route needs an
extension of the framework

`Z0_ret`/`Z1_ret` have a route **named and recorded** (`notes/parameter_register.md:169`: "named
route: the track-3/Gate-6 nonlinear return would derive the transmissions from the medium … —
`PENDING`"), which is §3.1's `A-REDUCIBLE-UNDERIVED`. But that route requires a **nonlinear closure
the framework does not currently contain** (sim-deferred, note `:774-777`), and §3.1.1's test — "would
an extension of the framework open a route? **Yes** ⇒ `A-IRREDUCIBLE-STRUCTURAL`" — points the other
way. Both are tier 1; the **sub-bucket** differs, and the sub-bucket is the deliverable (§5).

**Adopted:** §3.1.1 is titled and framed as the *STRUCTURAL-versus-POSTULATE* test, so it is applied
only once **no** route is named. Where a route is named and recorded, `A-REDUCIBLE-UNDERIVED` wins.
**Effect: all 15 tier-1 occurrences are `tier1-debt`; under the other reading the 7 `Z_ret`-family
occurrences would be `tier1-structural`, and `tier1-debt` would be 8.**

### B.5 GAP-4 — axis A for a composition over equally-unreduced symbols

`T0 = K0c/(K0c + Z0_ret)`, `K1 = K_eta + 2·T_Omega`, `P0_raw`, `A0`, `A1` are all *computed* from
symbols that are themselves tier-1. §3.1's `A-REDUCED` asks whether the value can be obtained "from
something **more fundamental**".

**Adopted:** a composition is `A-REDUCED` when the ledger's own accounting treats the targets as the
counted degrees of freedom and the composite as not a knob — which the register does explicitly
(`:169` "the residual prediction is parameterized by `ε_ℓ`"; `:170` "the `GENERATOR_DOFS` treat
`K0c`, `K_eta`, `T_Omega` as separate directions"). A **pure identification** `X ≡ Y` where `Y` is a
single free symbol (`K0 = K0c` at `PY:284`; `Z0 = Z0_ret` at `PY:287`; `P0port` at `WL:121`;
`ou = 1·OmegaU` at `WL:106`) is **not** a reduction and **inherits the target's axis-A value**.

### B.6 GAP-5 — whose citation establishes a `C-PEER` leaf

§3.3 requires "a **cited source locus**" for `C-PEER`; §3.3.1(1) requires the artifact to "cite the
locus of the equation being solved". **Neither engine cites a locus anywhere.** `PY:1049-1054` and
`WL:765-770` print stage *names* (`stage022`, `stage017`, `stage008`, `pathA_29`); `PY:336`, `:344`
and `WL:196` are comments naming a stage. The stage **note** does carry section loci (`:42-44`
"handoff §10.2–10.3", "§8.2 / Gate-2", "§9.4") — but the note is a different artifact (§7.2).

**Adopted:** the citation must be in the artifact carrying the occurrence. A stage name is not a
locus. ⇒ **no `C-PEER` leaf is established anywhere in either engine**, and by §3.3's rule an uncited
match is `C-SELF` hand transcription.

**⭐ Consequence, stated plainly and not softened:** no closure in either artifact reaches a
`C-FIELDEQ` or `C-EXTERNAL` leaf. **`PHYSICS-FED` is `false` or `C-UNRESOLVED` for all 136 rows; no
row reaches TIER 3; `DERIVED = 0`.** This is the §3.5 mechanism the axis exists to detect, measured
here rather than assumed.

### B.7 GAP-6 — the rank/`NullSpace` call under §3.3.1(1)

`sp.Matrix(rows).rank()` (`PY:375`) and `NullSpace`/`MatrixRank` (`WL:132`, `:136`, `:146`) are the
"operator" of §3.3.1(1). The operator contributes a `C-FIELDEQ` leaf **iff the artifact cites the
locus of the equation being solved**; it does not, and the constraint rows are assembled from the
artifact's own free symbols. §3.3.1(1)'s own able-to-fail test applies verbatim: "a solve assembled
entirely from the artifact's own constants, with no cited equation, is **not** `PHYSICS-FED`".

### B.8 GAP-7 — schema values that do not exist for cases present here

| # | case | missing rule | what I recorded |
|---|---|---|---|
| a | a **free symbol** declaration (`PY:163-177`; `WL:26-36`) | §7.2 defines a binding site as where a value is *assigned or asserted*; these assert a **domain**, not a value. §3.2 has no axis-B value for "declared with no value" | binding site kept; axis B = `B-DECLARED-LITERAL` (the declaration is typed in). ⛔ `B-UNADJUDICATED` would be wrong — the execution state is perfectly well established (nothing computes it) |
| b | Wolfram globals are never declared | the `WL` free symbols have **no assignment line at all**; their only assertion site is `$Assumptions` (`WL:26-36`) | `WL:26` used as the binding site for all 16 |
| c | `SELF-REFERENTIAL` for an unresolved closure | §3.3 says it is "undefined rather than false" | recorded `n/d` (not determined), not `false` |
| d | `should_be_basis` for "should be **convention**" | the enum (§6) offers only `named-route` / `physical-picture-expectation` / `none`, all phrased about *reduction* | `physical-picture-expectation` used for the 4 convention rows, flagged here |
| e | the **time convention** as a closure member | §3.3's walk is over expressions; the `exp(−iωt)` convention fixes the *sign* the result claims but is not an expression input to `A0` | included via a **semantic** hop (flagged in the row), not a data-dependence hop |
| f | an independent **variable** (`ω`) | axis A asks about reducibility; a frequency argument has no reducibility status | `A-UNADJUDICATED`, flagged |

### B.9 Identity / merges (§8.1, §8.3)

⛔ **The builder does not adjudicate merges.** Every QID below is minted by the builder; every merge
is **PROPOSED** and left for the physics review leg.

- **Cross-engine merges** (`PY` ↔ `WL`, e.g. `K0c`↔`K0c`, `T0`↔`T0dc`, `A0`↔`A0lead`): proposed on
  the strength of the stage's own dual-engine statement (note `:839-841`), **not** on name match.
  65 QIDs merged / **136 unmerged** — both denominators reported at §G.
- **`c_s`@023 vs `c_s0`@register `:268`** — different spellings, and the register's R1 row is about
  `c_s0`. ⛔ **Left as two QIDs, flagged as an open identity question.** This is why `c_s` is *not*
  recorded `A-REDUCED` here (§C row PY-03).
- **`a`@023 vs `a` (pin)@register `:132`** — same open question, same treatment.
- `QID_REGISTRY.md` (§8.3) does not exist (schema §13 `NOT-YET-CREATED`); the merges above are the
  registry's input, not a substitute for it.

### B.10 Legend

Axis A: `RED`=`A-REDUCED` · `RU`=`A-REDUCIBLE-UNDERIVED` · `IS`=`A-IRREDUCIBLE-STRUCTURAL` ·
`IP`=`A-IRREDUCIBLE-POSTULATE` · `CAL`=`A-CALIBRATED` · `AU`=`A-UNADJUDICATED`.
Axis B: `EX`=`B-EXECUTED` · `DIF`=`B-DERIVED-IN-FORM-UNEXECUTED` · `DL`=`B-DECLARED-LITERAL` ·
`AT`=`B-ASSERTED-TARGET` · `BU`=`B-UNADJUDICATED`.
`PF` = `PHYSICS-FED`: `F`=false · `U`=`C-UNRESOLVED`. `SR` = `SELF-REFERENTIAL` (`n/d` = undefined
under an unresolved closure, §3.3). `CONV` = `CONVENTION-LADEN`. `is`/`sb` = `is_tier`/
`should_be_tier`, values from §5.7 abbreviated: `T1d`=`tier1-debt` · `T3`=`tier3-emergent` ·
`NTu`=`no-tier:unadjudicated` · `NTn`=`no-tier:unclassified-nonfed` · `NTc`=`no-tier:convention`.
`bas` = `should_be_basis`: `nr`=`named-route` · `ppe`=`physical-picture-expectation` · `–`=`none`.
`Δ` = delta flagged. **Every row is LIVE** — no retired quantity appears in either engine (§7.4;
checked against the register's two retired rows `λ_Pu` `:139` and `α_aniso` `:159`, neither of which
occurs in either artifact).

Occurrence ID = `QID:<name>@<artifact>#<line>`; artifact is `PY` or `WL` as defined above. `EV` cites
the evidence block at §E.

---

## C. Rows — artifact `PY` (66 occurrences)

### C.1 Free-symbol declarations (16)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-01 | `omega` | `PY:163` | W-A0: `A0`(`:342`)←`ω` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E9 |
| PY-02 | `a` | `PY:164` | W-A0: `A0`(`:342`)←`a` | AU | DL | `{C-FREE}` | U | n/d | **UNADJ** | NTu | NTc | ppe | ✔ | E1,E5 |
| PY-03 | `c_s` | `PY:165` | W-A0: `A0`(`:342`)←`c_s` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | T3 | nr | ✔ | E1,E4 |
| PY-04 | `M0` | `PY:166` | W-A0: `A0`(`:342`)←`M0` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E10 |
| PY-05 | `D1` | `PY:167` | W-A1: `A1`(`:343`)←`D1` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E10 |
| PY-06 | `D0` | `PY:169` | W-RU: `rank0`(`:388`)←`rows`(`:387`)←`base_constraints`(`:384`)←`P0_raw`(`:225`)←`D0` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| PY-07 | `K0c` | `PY:170` | W-RU: `base_constraints`(`:384`)←`K0c`; also `T0`(`:302`) | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E2 |
| PY-08 | `K_eta` | `PY:170` | W-RU: `base_constraints`(`:384`)←`K_eta` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E2 |
| PY-09 | `T_Omega` | `PY:170` | W-RU: `base_constraints`(`:384`)←`T_Omega` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E2 |
| PY-10 | `Z0_ret` | `PY:171` | W-RU: `grad_T0`(`:390`)←`T0`(`:302`)←`Z0_ret`; witness `:397` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E3 |
| PY-11 | `Z1_ret` | `PY:171` | W-RU: `grad_T1`(`:391`)←`T1`(`:303`)←`Z1_ret` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E3 |
| PY-12 | `Omega_U` | `PY:175` | W-RU: …←`P0_raw`(`:225`)←`P_port`(`:222`)←`Omega_U` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| PY-13 | `Omega_W` | `PY:175` | W-RU: …←`Delta_port`(`:223`)←`Omega_W` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| PY-14 | `R_mix` | `PY:175` | W-RU: …←`P_port`(`:222`)←`R_mix` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| PY-15 | `g_U` | `PY:175` | W-RU: …←`P_port`(`:222`)←`g_U` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| PY-16 | `g_W` | `PY:175` | W-RU: …←`P_port`(`:222`)←`g_W` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |

### C.2 The ℓ=2 port kernel (5)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-17 | `omega_u` | `PY:221` | W-RU via `P_port`(`:222`) | AU | EX | `{C-FREE(Omega_U), C-SELF(branch flag)}` | U | n/d | false | NTu | NTu | – | – | E6,E11 |
| PY-18 | `P_port` | `PY:222` | W-RU via `P0_raw`(`:225`) | RED | EX | `{C-FREE ×4, C-MATH(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| PY-19 | `Delta_port` | `PY:223` | W-RU via `P0_raw`(`:225`) | RED | EX | `{C-FREE ×3, C-MATH}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| PY-20 | `N0_from_port` | `PY:224` | W-RU via `P0_raw`(`:225`) | RED | EX | `{P_port, Delta_port ⇒ C-FREE ×5, C-MATH}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| PY-21 | `P0_raw` | `PY:225` | W-RU: `base_constraints`(`:384`)←`P0_raw` — **the ℓ=2 named constraint row** | RED | EX | `{C-FREE ×6, C-MATH}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |

### C.3 Transfers, admittances and the branch condition (9)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-22 | `K0` | `PY:284` | W-RU: `T0`(`:302`)←`K0` | **RU** | EX | `{C-FREE(K0c)}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E2,E7 |
| PY-23 | `K1` | `PY:285` | W-RU: `T1`(`:303`)←`K1`; also re-typed inline at `:384` | RED | EX | `{C-FREE(K_eta,T_Omega), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| PY-24 | `Z0` | `PY:287` | W-RU: `T0`(`:302`)←`Z0` | **RU** | EX | `{C-FREE(Z0_ret)}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E3,E7 |
| PY-25 | `Z1` | `PY:288` | W-RU: `T1`(`:303`)←`Z1` | **RU** | EX | `{C-FREE(Z1_ret)}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E3,E7 |
| PY-26 | `T0` | `PY:302` | W-RU: `grad_T0`(`:390`)←`T0`; W-A0: `A0`(`:342`)←`T0` | RED | EX | `{C-FREE(K0c,Z0_ret)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| PY-27 | `T1` | `PY:303` | W-RU: `grad_T1`(`:391`)←`T1`; W-A1: `A1`(`:343`)←`T1` | RED | EX | `{C-FREE(K_eta,T_Omega,Z1_ret), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| PY-28 | `epsilon0` | `PY:304` | W-A0: `expected_A0`(`:345`)←`eps0` | RED | EX | `{C-FREE(Z0_ret,K0c)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| PY-29 | `epsilon1` | `PY:305` | W-A1: `expected_A1`(`:346`)←`eps1` | RED | EX | `{C-FREE(Z1_ret,K_eta,T_Omega), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| PY-30 | `positive_bounded` | `PY:316` | W-A0/W-A1: the branch condition the residual class is earned **conditional on** (`REPORTED_RESULTS` §2.1 cites `:316`) | RED | EX | `{T0,T1,eps0,eps1 ⇒ C-FREE ×5, C-SELF}` | U | n/d | false | NTu | NTu | – | – | E7 |

### C.4 The `A0`/`A1` residual family (11)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-31 | `v0` | `PY:337` | W-A0: `A0`(`:342`)←`v0` | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| PY-32 | `v1` | `PY:338` | W-A1: `A1`(`:343`)←`v1` | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| PY-33 | `power1` | `PY:339` | W-A1: `A1`(`:343`)←`power1` (the ω-power typed in) | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| PY-34 | `A0` | `PY:342` | W-A0 hop 1 (the result's own quantity) | RED | EX | `{C-MATH(i), C-SELF(v0), C-FREE(a,ω,c_s,M0), T0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| PY-35 | `A1` | `PY:343` | W-A1 hop 1 | RED | EX | `{C-MATH(i), C-SELF(v1,power1), C-FREE(a,ω,c_s,D1), T1 ⇒ C-FREE ×3}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| PY-36 | `expected_A0` | `PY:345` | W-A0: `A0_residual`(`:349`)←`expected_A0` | RED | EX | `{C-MATH(i), C-FREE(a,ω,M0,c_s), eps0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| PY-37 | `expected_A1` | `PY:346` | W-A1: `A1_residual`(`:350`)←`expected_A1` | RED | EX | `{C-MATH(i), C-SELF(2,3), C-FREE(a,ω,D1,c_s), eps1 ⇒ C-FREE ×3}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| PY-38 | `A0_residual` | `PY:349` | W-A0 hop 0 — the reported coincidence | RED | EX | `{A0, expected_A0 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-39 | `A1_residual` | `PY:350` | W-A1 hop 0 | RED | EX | `{A1, expected_A1 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-40 | `A0_order` | `PY:351` | W-A0 hop 0 — the reported `ω¹` | RED | EX | `{A0 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-41 | `A1_order` | `PY:352` | W-A1 hop 0 — the reported `ω³` | RED | EX | `{A1 leaves, incl. C-SELF(power1)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |

### C.5 The rank / nullspace audit (15)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-42 | `GENERATOR_DOFS` | `PY:205` | W-RU hop 0 — the "11 genuine generator dofs" of the claim | AU | DL | `{C-FREE ×11}` | U | n/d | false | NTu | NTu | – | – | E15 |
| PY-43 | `base_constraints` | `PY:384` | W-RU hop 0 — "the three collected named constraints" | AU | DL | `{P0_raw, C-FREE(K0c,K_eta,T_Omega), C-SELF(2)}` | U | n/d | false | NTu | NTu | – | – | E16 |
| PY-44 | `constraint_jacobian` | `PY:387` | W-RU: `rank0`(`:388`)←`rows` | RED | EX | `{base_constraints, GENERATOR_DOFS}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-45 | `rank0` | `PY:388` | W-RU hop 0 — the reported rank 3 | RED | EX | `{constraint_jacobian leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E17 |
| PY-46 | `native_nullspace_dimension` | `PY:389` | W-RU hop 0 — the reported 8 | RED | EX | `{rank0, len(dofs) ⇒ C-SELF}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-47 | `grad_T0` | `PY:390` | W-RU: `return_augmented_rank`(`:392`)←`grad_T0` | RED | EX | `{T0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-48 | `grad_T1` | `PY:391` | W-RU: `return_augmented_rank`(`:392`)←`grad_T1` | RED | EX | `{T1 ⇒ C-FREE ×3, C-SELF}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-49 | `return_augmented_rank` | `PY:392` | W-RU hop 0 — the reported 5 | RED | EX | `{constraint_jacobian, grad_T0, grad_T1}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-50 | `return_moving_nullity` | `PY:393` | W-RU hop 0 — ⭐ **the verdict-bearing quantity, 2** | RED | EX | `{return_augmented_rank, rank0}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E17 |
| PY-51 | `witness_unit_vector` | `PY:397` | W-RU hop 0 — the constructive witnesses `e_{Z_ℓ,ret}` | RED | EX | `{C-SELF}` (the value is 0/1 by dof position only) | **F** | false | false | **NTn** | NTn | – | – | E12,E18 |
| PY-52 | `constraint_deltas` | `PY:399` | W-RU hop 0 — "preserves every constraint row" | RED | EX | `{constraint_jacobian, witness_unit_vector}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-53 | `delta_T0` | `PY:402` | W-RU hop 0 — the reported `ΔT0 ≠ 0` | RED | EX | `{grad_T0, witness ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-54 | `delta_T1` | `PY:403` | W-RU hop 0 — the reported `ΔT1 ≠ 0` | RED | EX | `{grad_T1, witness ⇒ C-FREE ×3, C-SELF}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-55 | `moves_return` | `PY:411` | W-RU hop 0 | RED | EX | `{delta_T0, delta_T1}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| PY-56 | `native_null_moves_return` | `PY:422` | W-RU hop 0 — the boolean the verdict reads | RED | EX | `{return_moving_nullity}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |

### C.6 Asserted targets — values the artifact checks *against* (9)

⚠ These are the **same quantities** as their computed twins above, at a second binding site, typed in.
That is a §7.2 second occurrence and a §10.3 conflict; see §G.13.

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-57 | `n_constraint_rows` | `PY:807` | W-RU (target 3) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-58 | `n_generator_columns` | `PY:808` | W-RU (target 11) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-59 | `rank0` | `PY:809` | W-RU (target 3) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-60 | `native_nullspace_dimension` | `PY:810` | W-RU (target 8) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-61 | `return_augmented_rank` | `PY:811` | W-RU (target 5) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-62 | `return_moving_nullity` | `PY:812` | W-RU (target 2) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-63 | `n_witnesses` | `PY:819` | W-RU (target 2) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-64 | `A0_order` | `PY:884` | W-A0 (target 1) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| PY-65 | `A1_order` | `PY:885` | W-A1 (target 3) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |

### C.7 The time convention (1)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PY-66 | `time_convention` | `PY:612` | W-A0/W-A1 — **semantic hop** (GAP-7e): the results claim a **sign**, and `exp(−iωt)` is what fixes the sign of the `i` in `A0`(`:342`)/`A1`(`:343`) | AU | DL | `{C-SELF}` | F | false | **UNADJ** | NTu | NTc | ppe | ✔ | E20 |

---

## D. Rows — artifact `WL` (70 occurrences)

### D.1 Free symbols — asserted only in `$Assumptions` (16)

All at binding site **`WL:26`** (the `$Assumptions` head; the `Element` list runs `:28-31`, the
inequality block `:33-36`). GAP-7b. Witnesses mirror `PY` hop-for-hop through the `WL` loci given.

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WL-01 | `omega` | `WL:26` | W-A0: `A0lead`(`:199`)←`omega` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E9 |
| WL-02 | `a` | `WL:26` | W-A0: `A0lead`(`:199`)←`a` | AU | DL | `{C-FREE}` | U | n/d | **UNADJ** | NTu | NTc | ppe | ✔ | E1,E5 |
| WL-03 | `cs` | `WL:26` | W-A0: `A0lead`(`:199`)←`cs` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | T3 | nr | ✔ | E1,E4 |
| WL-04 | `M0` | `WL:26` | W-A0: `A0lead`(`:199`)←`M0` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E10 |
| WL-05 | `D1` | `WL:26` | W-A1: `A1lead`(`:200`)←`D1` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E10 |
| WL-06 | `D0` | `WL:26` | W-RU: `rank0`(`:146`)←`Jbase`(`:143`)←`baseConstraints`(`:124`)←`P0port`(`:121`)←`p0`(`:110`)←`D0` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| WL-07 | `K0c` | `WL:26` | W-RU: `baseConstraints`(`:124`)←`K0c` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E2 |
| WL-08 | `Keta` | `WL:26` | W-RU: `K1dc`(`:122`)←`Keta` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E2 |
| WL-09 | `TOmega` | `WL:26` | W-RU: `K1dc`(`:122`)←`TOmega` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E2 |
| WL-10 | `Z0ret` | `WL:26` | W-RU: `Greturn`(`:145`)←`T0dc`(`:138`)←`Z0ret` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E3 |
| WL-11 | `Z1ret` | `WL:26` | W-RU: `Greturn`(`:145`)←`T1dc`(`:139`)←`Z1ret` | **RU** | DL | `{C-FREE}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E1,E3 |
| WL-12 | `OmegaU` | `WL:26` | W-RU: …←`pport`(`:107`)←`OmegaU` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| WL-13 | `OmegaW` | `WL:26` | W-RU: …←`dport`(`:108`)←`OmegaW` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| WL-14 | `Rmix` | `WL:26` | W-RU: …←`pport`(`:107`)←`Rmix` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| WL-15 | `gU` | `WL:26` | W-RU: …←`pport`(`:107`)←`gU` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |
| WL-16 | `gW` | `WL:26` | W-RU: …←`pport`(`:107`)←`gW` | AU | DL | `{C-FREE}` | U | n/d | false | NTu | NTu | – | – | E1,E11 |

### D.2 Port kernel and the top-level constraint objects (10)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WL-17 | `omega_u` | `WL:106` | W-RU via `pport`(`:107`); `ou = 1·OmegaU` at `baselinePort`(`:120`) | AU | EX | `{C-FREE(OmegaU), C-SELF(factor)}` | U | n/d | false | NTu | NTu | – | – | E6,E11 |
| WL-18 | `P_port` | `WL:107` | W-RU via `p0`(`:110`) | RED | EX | `{C-FREE ×4, C-MATH(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| WL-19 | `Delta_port` | `WL:108` | W-RU via `p0`(`:110`) | RED | EX | `{C-FREE ×3, C-MATH}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| WL-20 | `N0_from_port` | `WL:109` | W-RU via `p0`(`:110`) | RED | EX | `{C-FREE ×5, C-MATH}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| WL-21 | `P0_raw` | `WL:110` | W-RU: `baseConstraints`(`:124`)←`P0port`(`:121`)←`p0` | RED | EX | `{C-FREE ×6, C-MATH}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6,E12 |
| WL-22 | `P0_raw` (2nd site) | `WL:121` | W-RU: re-export `P0port = baselinePort["P0Raw"]` | RED | EX | `{P0_raw@:110}` | U | n/d | false | NTu | T3 | ppe | ✔ | E6 |
| WL-23 | `K1` | `WL:122` | W-RU: `baseConstraints`(`:124`)←`K1dc` | RED | EX | `{C-FREE(Keta,TOmega), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-24 | `GENERATOR_DOFS` | `WL:123` | W-RU hop 0 — the 11 dofs | AU | DL | `{C-FREE ×11}` | U | n/d | false | NTu | NTu | – | – | E15 |
| WL-25 | `base_constraints` | `WL:124` | W-RU hop 0 — the three collected constraints | AU | DL | `{P0_raw, C-FREE(K0c), K1}` | U | n/d | false | NTu | NTu | – | – | E16 |
| WL-26 | `positive_bounded` | `WL:471` | W-A0/W-A1 branch condition (`caseFor` baseline) | RED | EX | `{t0,t1,e0,e1 ⇒ C-FREE ×5, C-SELF}` | U | n/d | false | NTu | NTu | – | – | E7 |

### D.3 Top-level transfers and rank objects (12)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WL-27 | `T0` | `WL:138` | W-RU: `Greturn`(`:145`)←`T0dc`; W-A0: `A0lead`(`:199`)←`T0dc` | RED | EX | `{C-FREE(K0c,Z0ret)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-28 | `T1` | `WL:139` | W-RU: `Greturn`(`:145`)←`T1dc`; W-A1: `A1lead`(`:200`)←`T1dc` | RED | EX | `{C-FREE(Keta,TOmega,Z1ret), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-29 | `epsilon0` | `WL:140` | W-A0: `expectedA0`(`:201`)←`eps0` | RED | EX | `{C-FREE(Z0ret,K0c)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-30 | `epsilon1` | `WL:141` | W-A1: `expectedA1`(`:202`)←`eps1` | RED | EX | `{C-FREE(Z1ret,Keta,TOmega), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-31 | `constraint_jacobian` | `WL:143` | W-RU: `rank0`(`:146`)←`Jbase` | RED | EX | `{base_constraints, GENERATOR_DOFS}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-32 | `null_basis` | `WL:144` | W-RU hop 0 — ⭐ the **constructive** route (`WL`-only object) | RED | EX | `{constraint_jacobian leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E21 |
| WL-33 | `return_gradients` | `WL:145` | W-RU: `returnMovingNullity`(`:148`)←`Greturn` | RED | EX | `{T0,T1 ⇒ C-FREE ×5, C-SELF}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-34 | `rank0` | `WL:146` | W-RU hop 0 — the reported rank 3 | RED | EX | `{constraint_jacobian leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E17 |
| WL-35 | `native_nullspace_dimension` | `WL:147` | W-RU hop 0 — the reported 8, as `Length[basis]` | RED | EX | `{null_basis}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E21 |
| WL-36 | `return_moving_nullity` | `WL:148` | W-RU hop 0 — ⭐ the verdict-bearing 2, as an **image rank** | RED | EX | `{null_basis, return_gradients}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E21 |
| WL-37 | `witness_unit_vector` | `WL:150` | W-RU hop 0 — `e_{Z0_ret}` | RED | EX | `{C-SELF}` | **F** | false | false | **NTn** | NTn | – | – | E12,E18 |
| WL-38 | `witness_unit_vector` (2nd) | `WL:151` | W-RU hop 0 — `e_{Z1_ret}` | RED | EX | `{C-SELF}` | **F** | false | false | **NTn** | NTn | – | – | E12,E18 |

### D.4 Top-level `A0`/`A1` family (8)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WL-39 | `v0` | `WL:197` | W-A0: `A0lead`(`:199`)←`v0` | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| WL-40 | `v1` | `WL:198` | W-A1: `A1lead`(`:200`)←`v1` | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| WL-41 | `A0` | `WL:199` | W-A0 hop 1 | RED | EX | `{C-MATH(I), C-SELF(v0), C-FREE(a,omega,cs,M0), T0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-42 | `A1` | `WL:200` | W-A1 hop 1 (ω-power **typed inline** as `^3`) | RED | EX | `{C-MATH(I), C-SELF(v1,3), C-FREE(a,omega,cs,D1), T1 ⇒ C-FREE ×3}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-43 | `expected_A0` | `WL:201` | W-A0: `resA0`(`:203`)←`expectedA0` | RED | EX | `{C-MATH(I), C-FREE(a,omega,M0,cs), eps0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-44 | `expected_A1` | `WL:202` | W-A1: `resA1`(`:204`)←`expectedA1` | RED | EX | `{C-MATH(I), C-SELF(2,3), C-FREE(a,omega,D1,cs), eps1 ⇒ C-FREE ×3}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-45 | `A0_residual` | `WL:203` | W-A0 hop 0 | RED | EX | `{A0, expected_A0 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-46 | `A1_residual` | `WL:204` | W-A1 hop 0 | RED | EX | `{A1, expected_A1 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |

### D.5 The `caseFor["baseline"]` recomputation — second binding sites in the same artifact (18)

⚠ §7.2: the same quantity at a second binding site in one artifact is a **second occurrence**. The
`caseFor` memoized block rebuilds the whole baseline chain; the baseline verdict of
`RES:023:RETURN-UNDERDETERMINED` is read from **this** copy (`WL:613`, `:762`), not from the
top-level copy, so it is in the closure and cannot be dropped.

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WL-47 | `Z0` | `WL:430` | W-RU via `t0`(`:447`) | **RU** | EX | `{C-FREE(Z0ret)}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E3,E7 |
| WL-48 | `Z1` | `WL:431` | W-RU via `t1`(`:448`) | **RU** | EX | `{C-FREE(Z1ret)}` | U | n/d | false | **T1d** | T3 | nr | ✔ | E3,E7 |
| WL-49 | `v1` (2nd) | `WL:432` | W-A1 via `a1`(`:453`) | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| WL-50 | `power1` | `WL:433` | W-A1 via `a1`(`:453`) | RED | **DL** | `{C-SELF}` | **F** | false | false | **NTn** | T3 | nr | ✔ | E8,E13 |
| WL-51 | `T0` (2nd) | `WL:447` | W-RU via `moving`(`:460`) | RED | EX | `{C-FREE(K0c,Z0ret)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-52 | `T1` (2nd) | `WL:448` | W-RU via `moving`(`:460`) | RED | EX | `{C-FREE(Keta,TOmega,Z1ret), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-53 | `epsilon0` (2nd) | `WL:449` | W-A0 via `expected0`(`:454`) | RED | EX | `{C-FREE(Z0ret,K0c)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-54 | `epsilon1` (2nd) | `WL:450` | W-A1 via `expected1`(`:455`) | RED | EX | `{C-FREE(Z1ret,Keta,TOmega), C-SELF(2)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E7,E12 |
| WL-55 | `A0` (2nd) | `WL:452` | W-A0 hop 1 (⚠ `v0` is **inlined as 1** here, not read from `WL:197`) | RED | EX | `{C-MATH(I), C-SELF(1), C-FREE ×4, T0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-56 | `A1` (2nd) | `WL:453` | W-A1 hop 1 | RED | EX | `{C-MATH(I), C-SELF(coeff1,power1), C-FREE ×4, T1 ⇒ C-FREE ×3}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-57 | `expected_A0` (2nd) | `WL:454` | W-A0 via `form0`(`:456`) | RED | EX | `{C-MATH(I), C-FREE ×4, e0 ⇒ C-FREE ×2}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-58 | `expected_A1` (2nd) | `WL:455` | W-A1 via `form1`(`:457`) | RED | EX | `{C-MATH(I), C-SELF(2,3), C-FREE ×4, e1 ⇒ C-FREE ×3}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E14 |
| WL-59 | `A0_form_flag` | `WL:456` | W-A0 hop 0 — the coincidence, as a boolean | RED | EX | `{A0, expected_A0 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-60 | `A1_form_flag` | `WL:457` | W-A1 hop 0 | RED | EX | `{A1, expected_A1 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-61 | `A0_order` | `WL:458` | W-A0 hop 0 — the reported `ω¹` | RED | EX | `{A0 leaves}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-62 | `A1_order` | `WL:459` | W-A1 hop 0 — the reported `ω³` | RED | EX | `{A1 leaves, incl. C-SELF(power1)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12 |
| WL-63 | `return_moving_nullity` (2nd) | `WL:460` | W-RU hop 0 — the copy the baseline **verdict** reads (`:481`, `:613`) | RED | EX | `{null_basis, return_gradients(t0,t1)}` | U | n/d | false | NTu | T3 | ppe | ✔ | E12,E21 |
| WL-64 | `time_convention` | `WL:360` | W-A0/W-A1 — **semantic hop** (GAP-7e), as `PY-66` | AU | DL | `{C-SELF}` | F | false | **UNADJ** | NTu | NTc | ppe | ✔ | E20 |

### D.6 Asserted targets (6)

| # | QID | locus | witness | A | B | closure leaves | PF | SR | CONV | is | sb | bas | Δ | EV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WL-65 | `jacobian_shape` | `WL:560` | W-RU (target `{3, 11}`) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| WL-66 | `rank0` | `WL:561` | W-RU (target 3) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| WL-67 | `native_nullspace_dimension` | `WL:562` | W-RU (target 8) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| WL-68 | `return_moving_nullity` | `WL:563` | W-RU (target 2) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| WL-69 | `A0_order` | `WL:640` | W-A0 (target 1) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |
| WL-70 | `A1_order` | `WL:641` | W-A1 (target 3) | RED | **AT** | `{C-SELF}` | F | false | false | **NTn** | NTn | – | – | E19 |

---

## E. Evidence blocks

Every locus in this section was opened.

- **E1 — free-symbol declarations carry no value.** `PY:163-177` declares 22 symbols via
  `sp.Symbol`/`sp.symbols` with `real`/`positive`/`nonzero` assumptions and **no value**; no later
  line assigns one (the only re-bindings are the local identifications at `PY:284`, `:287-288`,
  which are recorded as their own occurrences). `WL` never declares its symbols at all; `WL:26-36`
  asserts only domains and inequalities. ⇒ leaf tag `C-FREE` (§3.3). Axis-B value per GAP-7a.
- **E2 — `K0c` / `K_eta` / `T_Omega`: named route, recorded.** `notes/parameter_register.md:170`,
  opened: class `FREE-UNREDUCED` **PENDING scalar-reduction — COUNTED as reduction debt**; the route
  is stated as "an explicit scalar-reduction (profile+measure) of the ℓ=0 Gate-2 collective (§8.2,
  stage013) + the ℓ=1 §9.4 harmonic (stage017) sectors to the calibrated wall packet would DISCHARGE
  this — `PENDING` (no equation yet; do NOT assert DERIVED)". Corroborated at the same file `:309`
  (edge **R42**) and note `:808-815`. ⇒ `A-REDUCIBLE-UNDERIVED` with §9's "named route + where
  recorded" satisfied. Also `notes/stages/ledger_stage023_nullspace_underdetermination.md:812-814`
  records the **non**-identification with stage013/017 (dims `M T⁻²` vs `M L⁻¹T⁻²` / `M L⁻³T⁻²`) —
  i.e. the route has *not* been executed.
- **E3 — `Z0_ret` / `Z1_ret`: named route, recorded.** `notes/parameter_register.md:169`, opened:
  the `ε0, ε1` row, `FREE-UNREDUCED`, "named route: the track-3/Gate-6 nonlinear return would derive
  the transmissions from the medium (the same wall that discharges R23; R42 = the sharpened form) —
  `PENDING`", and "stage023's `Z0_ret/Z1_ret` are the COORDINATE ALIASES of these SAME 2 dofs". Edge
  R42 at `:309`. ⇒ `A-REDUCIBLE-UNDERIVED`. ⚠ See GAP-3/GAP-4 for the `tier1-structural` alternative.
- **E4 — `c_s`, why NOT `A-REDUCED`.** `notes/parameter_register.md:268` (R1) reads
  `c_s0 = √(5K ρ0⁴/m_GNLS)` | `DERIVED` | I-1/I-2 (004/005) | "collapses `c_s0` into
  `{K, ρ0, m_GNLS}`". The register's symbol is **`c_s0`**; this artifact's is `c_s`. §8.1 forbids
  joining by name and §8.3 forbids the builder adjudicating the merge ⇒ axis A stays
  `A-UNADJUDICATED`; the reduction is carried as the `should_be` **named route**, conditional on the
  merge. Stage note `:803` asserts the identity; per §9.1 rule 1 an attribution is not evidence.
- **E5 — `a`, the convention claim.** `notes/parameter_register.md:132`: `` `a` (pin) | `L` | I-1
  (004) | `CONV` | `= ħ/(m_GNLS c_s0)` | pin choice — never free ``. `PY:1054` / `WL:770` print "a
  from CONV". ⇒ a convention claim **is** made (§3.4), so `false` is not available. Clause (a) is
  plausible but no transformation group is exhibited in either artifact; clause (b) is **not
  demonstrated** — no dimensionless observable of the model is computed from `a` here (`aω/c_s`
  appears but `ω`, `c_s` carry no values), which §3.4's ⭐⭐ rule makes `UNADJUDICATED`, never
  `CONVENTION`. ⇒ flag `UNADJUDICATED`, no exclusion bought, row keeps its axis-A tier.
- **E6 — port-kernel computation.** `PY:220-232` (`p_port` `:222`, `delta_port` `:223`, `n0_port`
  `:224`, `p0_raw` `:225`); `WL:105-118` (`pport` `:107`, `dport` `:108`, `nport` `:109`, `p0`
  `:110`), instantiated `WL:120` and re-exported `WL:121`. Input leaves are the five port generators
  plus `D0` (`PY:169`, `:175`; `WL:26`). No locus is cited in either artifact for the kernel's form;
  the only provenance marker is the tag string `"stage017_grouped_port_kernel"` (`PY:589`,
  `WL:351`) and the print at `PY:1052` / `WL:768` — stage names, not loci (GAP-5).
- **E7 — transfer computation.** `PY:284-305` (`k0`, `k1`, `z0`, `z1`, `t0`, `t1`, `eps0`, `eps1`),
  branch condition `PY:273-280` + `:316`; `WL:122`, `:138-141`, `:420-422`, `:447-450`, `:471`.
  Input leaves: `K0c`, `K_eta`, `T_Omega`, `Z0_ret`, `Z1_ret` (`PY:170-171`; `WL:26`) and the
  literal `2` (C-SELF).
- **E8 — `v0`, `v1`, `power1` are typed literals here, and the reduction lives at stage022.**
  Declared `PY:337-339` (`sp.Integer(1)`, `sp.Rational(1,2)`, `3`) and `WL:197-198`, `:432-433`;
  the `WL:200` copy inlines the power as `^3`. The comment at `PY:336` ("CHECKABLE cited input from
  stage022") and `WL:196` name a stage, not a locus ⇒ `C-SELF`, not `C-PEER` (§3.3). The reduction
  itself was opened: `scripts/ledger_stage022_cross_l_fingerprints_sympy_audit.py:209-210` computes
  `radiative_power = 2*lval + 1` and `radiative_coeff = y_series.coeff(z, radiative_power)/sp.I`
  from the spherical-Hankel series, checked against the typed table at `:136-140`
  (`{0: 1, 1: 1/2, 2: 1/27}`) at `:274-275`; recorded as edge **R41** at
  `notes/parameter_register.md:308` ("series-expands to the DERIVED radiative fingerprint"). ⇒ axis A
  `A-REDUCED` (reduction exists and is recorded, elsewhere); axis B `B-DECLARED-LITERAL` here. This
  is §4 near-miss shape 3 in structure but **not** in fact: `PHYSICS-FED` is `false`, not `true`, so
  the row is `unclassified-nonfed`, not `physics-fed-but-declared-literal`.
- **E9 — `ω` is an independent variable** (GAP-7f): it is the drive frequency the response is a
  function of, not a model constant; the schema has no axis-A value for it.
- **E10 — `M0`, `D1`.** No register row assigns them a class; the `ε0/ε1` row
  (`notes/parameter_register.md:169`) instead carries two open **work pointers** on the `M0` used
  here (`WORK-023-STAGE009-MOMENT0`, `WORK-023-MOMENT-CONVENTION`), explicitly "**no
  adjudication**". No route, no benchmark, no stated postulate ⇒ `A-UNADJUDICATED`. This is §5.3's
  named case (a bare knob satisfying no §9 evidence row).
- **E11 — `D0` and the four port generators.** Grepped `notes/parameter_register.md` for
  `Omega_U`/`OmegaU`/`g_U`/`g_W`/`R_mix`: the only hits are inside the R42 narrative `:309`, which
  *lists* them as generator dofs and assigns no class or route. `D0` appears in R40 `:307` as "the
  carried reduced static conservative denominator `D₀ = K − B₀ − Z₀`" — a candidate named route, but
  under a different stage's symbol and therefore an unadjudicated merge (§8.1) ⇒ not usable as axis-A
  evidence here. All five ⇒ `A-UNADJUDICATED`.
- **E12 — execution loci.** Every `EX` row's computation locus is the row's own locus; the exactness
  wrappers are `PY:59-66` (`compact`) / `WL:53` (`clean`), the rank kernel `PY:369-375` and
  `WL:126-136`, `:146`, and the whole baseline chain is driven from `PY:690-731` (`run_gate`) and
  `WL:424-494` (`caseFor`). Input-leaf loci are as listed in each row's closure column, resolving to
  `PY:163-177` / `WL:26` for every `C-FREE` leaf.
- **E13 — `should_be` named route for `v0`/`v1`/`power1`** = E8's stage022 loci
  (`…stage022_cross_l_fingerprints_sympy_audit.py:209-210`, `:136-141`) and register `:308`.
- **E14 — the `expected_A*` forms.** `PY:344` comments "Independently cited pathA_29 forward residual
  form"; `PY:345-348`, `WL:201-202`, `:454-455`. No locus is cited in either artifact. The nearest
  recorded relation is register R24 `:291` (`Z = −M0·ε0/(1+ε0)`, `DERIVED` at II-B2/009) — a
  different expression under a different symbol, so it is **not** used as evidence here.
- **E15 — the dof set.** `PY:205-217`, `WL:123`. The *membership* of the 11-dof list is typed; the
  register's R42 row `:309` restates the same list but records no derivation of it ⇒
  `A-UNADJUDICATED` for the set as an object.
- **E16 — the constraint set.** `PY:384`, `WL:124`. The claim that these three are *the* collected
  Gate-5 named constraints is asserted, not computed; the note's citations for them (`:42-44`:
  "handoff §10.2–10.3", "§8.2 / Gate-2", "§9.4") live in a **different artifact** (GAP-5) ⇒
  `A-UNADJUDICATED` for the set. ⚠ This is the single assertion the whole nullity result rests on.
- **E17 — the reported rank/nullity.** Computed `PY:388-393` from `PY:387`; `WL:146-148` from
  `WL:143-145`. Both engines' asserted targets are §C.6 / §D.6.
- **E18 — the witnesses are position-only.** `PY:397` builds `[1 if dof == return_dof else 0 …]`;
  `WL:150-151` build `UnitVector[11, Position[rankDofs, Z0ret]]`. The value depends on the dof
  **ordering** (`PY:205`, `WL:123`), not on any symbol's value ⇒ closure `{C-SELF}`, determined,
  `PHYSICS-FED = false`. These are the only `B-EXECUTED ∧ PHYSICS-FED=false` rows in the stage.
- **E19 — asserted targets.** `PY:807-812`, `:819`, `:884-885`; `WL:560-563`, `:640-641`. Each is a
  literal the artifact checks the computed value **against**; the reduction is performed at the
  adjacent computing locus (E17), so axis A is `A-REDUCED` while axis B is `B-ASSERTED-TARGET`.
- **E20 — the time convention.** `PY:612-615` and `WL:360` declare the item
  `time_convention` with tag `exp_minus_i_omega_t` and `computed_class = "convention"`. Clause (a)'s
  group (conjugation `e^{-iωt} ↔ e^{+iωt}`) is not exhibited and clause (b) is not demonstrated ⇒
  flag `UNADJUDICATED`.
- **E21 — the `WL` engine is a different algorithm, so its objects are not `PY`'s under another
  name.** `WL:132` (`NullSpace`), `:136` (`MatrixRank[basis . Transpose[gradients]]`), `:144`,
  `:147-148`. There is **no** `return_augmented_rank` in `WL` — `PY:392`'s object has no `WL`
  counterpart. Recorded at §G.14.

---

## F. Out of scope (§7.3) — recorded with reasons and counts, never dropped

### F.1 `reached-by-no-reported-result` (§5.8 output 16) — **the count that is itself a finding**

Under the adopted closure reading (§B.2 reading A):

| set | loci | occurrences |
|---|---|---|
| `PY` `SOURCED_DIMS` exponent triples | `PY:460-483` | 22 |
| `PY` `EXPECTED_DIMS` | `PY:486-494` | 7 |
| `PY` emitted `dimension_records` join keys | `PY:539-571` | 29 |
| `PY` `P0_physical` (used only by the dimensional gate, `PY:517`) | `PY:231` | 1 |
| `PY` `forward_relations_ok` (a labelled SELF-CONSISTENCY identity, `PY:310-315`, `:886-889`) | `PY:310` | 1 |
| `PY` provenance items + their computed/emitted class strings | `PY:575-616` | 9 |
| `PY` `pathA29_premise` (`Z_is_premise`, `boundary_dof`) — ⚠ see F.4 | `PY:425` | 2 |
| `WL` `baseDims` | `WL:247-255` | 22 |
| `WL` `expectedDims` | `WL:257-262` | 7 |
| `WL` `printDimRecord` emissions | `WL:298-332` | 29 |
| `WL` `P0Physical` | `WL:116` | 1 |
| `WL` provenance items | `WL:341-361` | 9 |
| `WL` `pathA29Premise` | `WL:183` | 2 |
| **total** | | **141** |

⛔ **141 value-bearing occurrences in these two artifacts reach no reported physical result** under
reading (A) — more than the 136 that do. ⚠ Under reading (B) (§B.2) most of this set moves **into**
the universe. The number is reported because §7.1.1 requires it, not as a judgement.

### F.2 Controls and deliberate negatives (§7.3)

`q_free`, `eta_null`, `gain0`, `gain1` (`PY:172-174`; `WL:26`) — 4 per engine, **8**; the `Mutation`
flag defaults (`PY:184-202`) — 13; the counterfactual selector `{Z0_ret=K0c, Z1_ret=K1}`
(`PY:236`, `:240-242`; `WL:153`, `:437-438`) — 4; the mutation-mode branches in `WL:436-445` — 9;
the corrupted-copy objects (`PY:507`, `:509`; `WL:205-207`, `:289`, `:293`) — 7. **Total 41.**
Each is a construction built to make a check fail (`REPORTED_RESULTS` §2.2 items 1, 12).

### F.3 Test scaffolding, display and framework constants

The `expect_zero`/`expectZero` implicit zero target (`PY:124`; `WL:71`), the verdict token strings
(`PY:27-35`; `WL:16-24`), the failure-priority ladder ordering (`PY:674-682`; `WL:401-409`), the
banner/format helpers, `PASS_COUNT`/`FAIL_COUNT`. Not enumerated individually; **not** value-bearing
about the medium.

### F.4 ⚠ One exclusion a reviewer should attack

`pathA29_premise = {"Z_is_premise": True, "boundary_dof": "none"}` (`PY:425`; `WL:183`) is called
"the keystone premise" by the stage (note `:59-61`) and is the stated **reason** the
underdetermination is earned — yet **no expression consumes it**; it is printed only (`PY:814-818`,
`WL:623-627`). Under §3.3's expression-walk it is therefore out of the closure. If a reviewer holds
that the *justification* of a claim belongs to its closure, these 2 occurrences move in, and they are
the one place in the stage where `A-IRREDUCIBLE-POSTULATE` would be a live candidate (the premise is
stated as a defining property of the flat-slab family). ⛔ **Not adjudicated here.**

---

## G. The complete required output set (§5.8)

All counts are `is_tier` counts (§6.1). Occurrence denominator **136**; QID denominator **65**
(cross-engine merges PROPOSED, §B.9); unmerged QID denominator **136**. No convention-laden
occurrence exists, so nothing is quotiented out of either denominator (§10.1).

| # | output | occurrences | QIDs |
|---|---|---:|---:|
| 1 | **TIER 1, as a RANGE** — ⛔ never a scalar | **[15, 111]** | **[8, 57]** |
| — | `tier1-debt` (lower bound is entirely this bucket) | 15 | 8 |
| — | `tier1-structural` | 0 | 0 |
| — | `tier1-postulate` | 0 | 0 |
| 2 | **TIER 2 — calibrated** | 0 | 0 |
| 3 | **TIER 3 — emergent** | **0** | **0** |
| — | tier-3 split by axis B (§5.4): executed / unexecuted | 0 / 0 | — |
| — | tier-3 split by propagation (§5.5): `tier3-calibration-propagated` / `tier3-held-out` | 0 / 0 | — |
| 4 | **`DERIVED`** (`B-EXECUTED ∧ PHYSICS-FED ∧ A-REDUCED`) | **0** | **0** |
| 5 | near-miss `executed-but-not-physics-fed` | 3 | 1 |
| 6 | near-miss `derived-in-form-but-unexecuted` | 0 | 0 |
| 7 | near-miss `physics-fed-but-declared-literal` | 0 | 0 |
| 8 | `unclassified-nonfed` | 25 | 8 wholly + 6 mixed |
| 9 | convention bucket (`CONVENTION-LADEN = true`) | 0 | 0 |
| 10 | `unadjudicated` | 96 | 43 wholly + 6 mixed |
| 11 | `convention-unadjudicated` | 4 | 2 |
| 12 | `self-referential` | 0 | 0 |
| 13 | conflict set (§10.3) | — | **6** |
| 14 | stage043 delta, both directions | see §G.14 | |
| 15 | out-of-scope list with reasons | §F — 141 + 41 + scaffolding | |
| 16 | `reached-by-no-reported-result` | **141** | |
| 17 | reported-result set with loci | §A — 3 | |

⛔ **Do not sum this table** (§5.8): buckets 5–7 overlap the tiers by construction.

### G.1 Tier-1 range, computed

- **Occurrence level.** lower = 15 (`tier1-debt`, evidence E2/E3); upper = 15 + 96 (`unadjudicated`)
  = 111. ⇒ **TIER 1 = [15, 111]** of 136 occurrences.
- **QID level** (§5.6 level-consistency + §10.2.2 rule 4). lower = 8 QIDs with ≥ 1 tier-1 occurrence
  (`K0c`, `K_eta`, `T_Omega`, `Z0_ret`, `Z1_ret`, `K0`, `Z0`, `Z1`); upper = 8 + 43 wholly
  `no-tier:unadjudicated` QIDs + 6 mixed tier-less QIDs = 57. ⚠ §10.2.2 rule 4 says only QIDs
  **wholly** `no-tier:unadjudicated` widen the span, and does not say where a QID that mixes
  `unadjudicated` with `unclassified-nonfed` goes; I placed the 6 in the upper bound per §1's
  asymmetry. **Under the narrow reading the QID range is [8, 51].**
- **The 65 QIDs.** 16 imported free symbols · 5 port-kernel objects · 9 transfer/admittance objects ·
  11 `A0`/`A1`-family objects · 15 rank-audit objects · 4 `WL`-only objects (`null_basis`,
  `return_gradients`, `A0_form_flag`, `A1_form_flag`) · 4 assertion-only objects
  (`n_constraint_rows`, `n_generator_columns`, `n_witnesses`, `jacobian_shape`) · `time_convention`.

### G.13 The conflict set (6 QIDs) — computed from `is_tier` only

Each is a quantity the ledger **computes** at one binding site and **types in** at another, giving
two different `is_tier` values (`no-tier:unadjudicated` vs `no-tier:unclassified-nonfed`):

| QID | computed occurrence | asserted occurrence |
|---|---|---|
| `rank0` | `PY:388`, `WL:146` | `PY:809`, `WL:561` |
| `native_nullspace_dimension` | `PY:389`, `WL:147` | `PY:810`, `WL:562` |
| `return_augmented_rank` | `PY:392` | `PY:811` |
| `return_moving_nullity` | `PY:393`, `WL:148`, `WL:460` | `PY:812`, `WL:563` |
| `A0_order` | `PY:351`, `WL:458` | `PY:884`, `WL:640` |
| `A1_order` | `PY:352`, `WL:459` | `PY:885`, `WL:641` |

⛔ Not resolved by precedence (§10.3). ⚠ A reviewer may hold that an asserted target is not an
occurrence of the *same* quantity; that is an identity adjudication and belongs to the physics leg.

### G.14 stage043 reconciliation (§8.4) — **partial, by construction**

⛔ The 152-ID manifest is a corpus-wide object; one stage cannot reconcile it. What this pass can
map, from `scripts/ledger_stage043_irreducible_count_range_sympy_audit.py` (opened):

| `REG:` ID | locus | census QID(s) |
|---|---|---|
| `REG:b:K0c` | `:225` | `QID:K0c` |
| `REG:b:K1_eta_direction` | `:225` | `QID:K_eta`, `QID:T_Omega` (⚠ 1→2, an open identity question) |
| `REG:b:epsilon0` | `:224` | `QID:epsilon0`, and `QID:Z0_ret` as its coordinate alias (register `:169`) |
| `REG:b:epsilon1` | `:224` | `QID:epsilon1`, `QID:Z1_ret` |
| `REG:derived:Z0_ret_alias` | `:297` | `QID:Z0_ret` |
| `REG:derived:Z1_ret_alias` | `:297` | `QID:Z1_ret` |

- **Direction 1 — census QIDs with no `REG:` ID: 59 of 65.** Everything except the six above. Most
  are intermediate computed objects (`P0_raw`, `T0`, `A0`, `rank0`, …) that the register does not
  carry as knob rows; `M0`, `D1`, `D0` and the four port generators are the notable knob-shaped ones.
- **Direction 2 — `REG:` IDs recorded out of scope: not attempted.** The remaining ~146 IDs belong
  to other stages; recording each with its `ratified_category` is a census-wide task. ⛔ **Reported
  as not done, not as zero.**

---

## H. Coverage statement

**Covered: 136 of 136 in-universe occurrences** in the two artifacts, under the closure reading
adopted at §B.2 (reading A) and the occurrence rule at §B.5. Every row carries all three axes, the
closure with leaf tags, the two derived axis-C fields, the convention flag, `is_tier`,
`should_be_tier`, `should_be_basis`, `delta`, LIVE/RETIRED, IDs, a reachability witness, and §9
evidence for every value that requires it. Nothing was truncated.

**Not covered, with counts:**

1. **141 occurrences** excluded as `reached-by-no-reported-result` (§F.1) — enumerated by set with
   loci, not classified. ⚠ **Reading-sensitive**: under §B.2 reading (B) at least 122 of them are
   in-universe and would need classifying.
2. **41 occurrences** excluded as controls / deliberate negatives (§F.2) — enumerated, not classified.
3. **Test scaffolding and display constants (§F.3)** — not individually enumerated. This is the one
   place I did not produce a count; the class is excluded by §7.3 by name and none of it is
   value-bearing about the medium.
4. **The stage note as an artifact.** `notes/stages/ledger_stage023_nullspace_underdetermination.md`
   is an artifact under §7.2 and its tables are binding sites, so the stage has a **third** artifact's
   worth of occurrences (the §1.1–§1.3 display blocks restate `P0_raw`, `K1`, `A0`, `A1`,
   `expected_A0/A1`, `rank0`, `8`, `5`, `2` as prose/display values). The task scoped me to the two
   engines; **not covered, count not taken.**
5. **Direction 2 of the stage043 reconciliation** (§G.14): ~146 `REG:` IDs unreconciled.
6. **All merges are PROPOSED** (§B.9). If the physics leg rejects the cross-engine merges, the QID
   denominator is 136, not 65, and the QID-level tier-1 range changes accordingly.

**Blocked / unadjudicated rows:** 0 rows are *blocked* (every row received a full classification).
**96 of 136 rows carry `no-tier:unadjudicated`**, and that is a classification, not a blockage. The
96 split by *why*:

- **30** because **axis A itself** could not be established from evidence — the imported free symbols
  (`ω`, `a`, `c_s`, `M0`, `D1`, `D0`, `Omega_U`, `Omega_W`, `R_mix`, `g_U`, `g_W`, `omega_u`), the
  11-dof set, the 3-constraint set and the time convention, 15 per engine (E4, E9–E11, E15, E16, E20);
- **66** because axis A **is** determined (`A-REDUCED`) but the closure is `C-UNRESOLVED` — every
  branch terminates in a free symbol carrying no value, so no physics leaf is reachable (§B.6). ⭐
  These are the rows a single cited field-equation locus would move; they are not rows anyone failed
  to look at.

⚠ The `convention-unadjudicated` flag value on 4 rows (§G row 11) buys no exclusion (§3.4) and is
reported separately.

**Where the schema did not decide a case:** GAP-1 (§B.2, verdict-token closure — largest effect),
GAP-2 (§B.3, `C-UNRESOLVED` vs a determined axis A), GAP-3/4 (§B.4, named route that needs a
framework extension), GAP-5 (§B.6, whose citation makes a `C-PEER` leaf), GAP-6 (§B.7, the rank
kernel as an operator), GAP-7a–f (§B.8, six missing enum values), §G.1 (§10.2.2 rule 4 silent on
mixed tier-less QIDs), §F.4 (justification-vs-closure for the keystone premise). Each is stated with
the rule I adopted and the number of rows it moves.
