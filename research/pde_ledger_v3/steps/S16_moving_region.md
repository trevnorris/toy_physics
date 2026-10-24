# S16 — a live region around one moving throat

**Status, 2026-10-24:** record and register pass written for source-first review. The governing
S16 decision list is revision 4, cleared at `0464b80a`; the constructions are accepted at
`1d990ee6`, and the comparator at `f08db34a`, with its scoped citation-containment claim.
Acceptance of those instruments is not acceptance of this record. No engine or comparator was
rerun. The round-0 record is preserved at `4ae7e6cf`; this repair addresses Y1–Y4 and makes
no new commit. The edit account remains against `4d4d22cb`. (Retrieve M01–M03/M16–M17.)

**Deliverable:** record what the supplied balances give for a moving region outside one throat's
excluded core, preserve the cross-engine findings, file the specified future requirements and
S22 debts, and correct the named status/provenance lines. The conditional regional objects do
not deliver the plan's compact-body response approximation. No mass, inertial factor, passive
factor, surroundings field, potential, inter-body force law or throat solve is computed.
(Retrieve M01/M04/M07.)

Paths below are relative to `research/pde_ledger_v3/`, except where marked repository-root.
**Retrieve Mnn** means the exact command
`python3 research/pde_ledger_v3/scripts/S16_record_measurements.py --lookup Mnn`, run from the
repository root. Its complete literal stdout is filed in
[`_measurements/S16_record_measurements.md`](_measurements/S16_record_measurements.md).
M07 retrieves native printed payloads, M08 the filed comparator's tables/summary, M09 every
binding's own object text and comparison fields, M10 literal counts and exhaustive binding
membership, M11 sample coverage, M12 every unjoined tag and unpaired component, and M13 every
binding's parsed/compared counts. These are retrievals, not a new comparison instrument.
M17 retrieves the preserved round-0 dispositions and repair scope; M18 retrieves the supplied,
integral-definition and geometry comparisons added to the status account for Y4.

## 1. Supplied inputs, choices and their standing

The equations are supplied and **"unfalsifiable in this build"**. On S14a part 1's projected
three-coordinate domain, with the slab already integrated between its live faces, they are

```
∂t Σ + ∂i F^i − n₁ a₁,₋ + n₃ a₃,₊ = 0
∂t G^i + ∂j Π^{ij} = f^i .
```

The first is accepted part 1's projected material balance, where its committed material law
holds. The second is the user's general momentum balance. `Σ`, `F`, `G`, `Π`, `f` and the
face quantities stay live. The list says **"G^i is not related to F^i or Σ by any supplied
law, and Π^{ij} is not assumed symmetric."** It says `f` includes every face/bulk momentum
exchange but **"stays undivided"**. Each construction prints these supplied equations and its
solved local rates; this is no derivation of the physical inputs. The inherited real-fraction
projected-material standing/A13 remains as the accepted projection and register state it.
(Retrieve M01/M04/M06–M09/M14.)

The region is a shell about live `X_c(t)`, between live `r_in(t)` and `R(t)`, outside the
excluded core. **"That is a choice of control surface, not a symmetry of the fields"**;
**"r_in is not identified with the mouth radius."** The fields are **"the total fields on
the region"**. Nothing uses core values or a body/surroundings split. The native integral
definitions print the shell volume and each spherical area measure. The inner/outer surface
indices are Py 0/1 and WL 1/2 respectively, as the comparator's explicit surface table states.
(Retrieve M01/M06–M08.)

This table gives provenance, the list's literal standing (or expressly identified inherited
part-1 standing), and dependence. Its date is the date
the list records the decision/input, not a newly inferred date of physical validity. No O2
material/drive/exchange/bulk postulate is silently imported into these general equations.
(Retrieve M01/M04/M06/M14.)

| Input or decision; source/date | Standing quoted from the step list | Printed dependence or absence of one |
|---|---|---|
| Framing, user 2026-10-19, carried in the 2026-10-21 list | "This step states what the supplied balances give. It picks no material law, drain law or flow."; "Each force states the flow it needs, and S21 compares them" | Scope of the whole record and S21 provenance; no numerical/algebraic input is added. |
| "Keep S12", user 2026-10-20, same list | "the conversion law's forms stay S12's and stay OPEN" | The projected material input is supplied, not a completed conversion law; no chosen form enters an object. |
| Boundary, user 2026-10-21 | "Keep the flux live"; "Any condition under which these terms drop out is printed as a condition, not imposed." | Every count/dipole/momentum rate, each surface piece and the differentiated count/centre terms. The model point states the condition is that the corresponding object vanishes; no zero-flux premise is applied. |
| Momentum, user 2026-10-21 | "General balance"; "Supplied as an equation and labelled supplied." | Momentum definition, rates, rest/motion pieces, force companions and received momentum, independently of the count-centre objects. |
| Outside field, user 2026-10-21 | "No potential" | No potential or separately named outside field enters. The later two decisions expressly narrow the earlier promises of a gradient condition and two masses: those become requirements, not computed objects. |
| Passive-factor decision, user 2026-10-21 | "Record what it needs"; "It computes no passive factor."; "The missing owner is stated as missing." | Received momentum is printed apart from surface motion. Field/separation/normalization enter no computed object; `R-S16-01` records their absence. |
| Mass decision, user 2026-10-21 | "Count material; mass open"; "No mass is defined." | `N`, its moments and its centre are count objects. The relation of momentum to motion is OPEN, with the required momentum–material relation filed for S8. |
| Coordinate domain, S16 list 2026-10-21, inherited from part 1 | "On S14a part 1's projected domain: in-plane coordinates x^i (three) and time t" | The three-component balances, shell geometry/measures, moments and centre relations use this domain. It consumes the existing `R-S1-01` dimension obligation at S6 without deriving a dimension. |
| U1 support identification, user 2026-10-19, inherited through the accepted part-1 balance supplied by the 2026-10-21 list | Part 1 calls the support encoding "an identification, flagged" and the inputs "adopted substrate input to a conditional model" | The projected-material object's interpretation as the material in the shear-holding region inherits that identification. It consumes `R-S1-02` at S1 for this interpretation only; no shear law, confinement or order profile is derived. U1's re-ordering clause is not enforced. |
| Projected material equation, accepted part 1, carried 2026-10-21 | "Supplied equations (unfalsifiable in this build)"; "It holds where the committed material law holds" | Count/dipole rates, governed `Σ` time derivatives, count normalization and count-centre kinematics; moments use the supplied live density. |
| General projected momentum equation, user/list 2026-10-21 | "All three are general live fields of x and t."; `f` "stays undivided" | Momentum rates, integrated force companions, received momentum and surface-motion objects. No relation to material motion follows. |
| Reference point, core and total fields, list 2026-10-21 | "The throat's centre is a live point X_c(t)."; "The core's content stays with S22."; "No object of the build depends on values inside the core or on a split of a field into parts" | Shell measures, moving-reference terms and count moments/centres; no core solve or surroundings decomposition enters. |
| Control shape and retained moments, list/model point 2026-10-21 | "Labelled choices"; "the control surfaces' shape"; "moments printed through second order" | Both constructions choose concentric spheres and print count, dipole and all second-tensor components. No isotropic field or Gaussian width is adopted. |
| Regime, list/model point 2026-10-21 | "none chosen"; "no steady, quasi-static or slow limit is taken" | All time derivatives, including higher geometry jets, remain in the rates and differentiated rate terms. |
| Freezes, list/model point 2026-10-21 | "none intended"; "Nothing is imposed: no compactness, vanishing dipole, boundary-flux suppression, potential, mass-to-count ratio, or symmetry of the fields." | No such premise enters a printed object. The surface terms remain separate objects. |
| Transfer, list/model point 2026-10-21 | "wherever part 1's projected balance holds, outside the excluded core" | Domain of the conditional construction, without any claim about the core or a substrate solution. |
| Normalization, list 2026-10-21 | "Each object normalized by N is stated on the domain where N is nonzero."; "The unnormalized balances and moments keep the full domain" | Relative/absolute count-centre definition, rate and acceleration require `N ≠ 0`; unnormalized objects do not. |

These are supplied balances, user-selected scope, definitions and labelled control choices.
There is no newly adopted substrate response form or derived material identification. In
particular, the flagged missing identifications are count-to-mass, momentum-to-material motion,
control radius-to-mouth/profile width, and total-to-surroundings field. None is supplied by
notation or by dimensional homogeneity. O2's `ℐ_br^live` and `ℛ_br` are OPEN interface objects;
its designation of S16 as an orbital-response matching interface is explicitly an **inference**.
(Retrieve M01/M04/M06/M14.)

## 2. What items 1–5 print, part by part

The following is an index of the actual printed objects, with the compact integrand descriptions
read from their native payloads. `V` here denotes the printed shell-volume operator (`Volume`
in Py); `B_s` denotes its printed surface operator (`Surface` in Py). `y_i = x_i − X_ci`
is the reference displacement already used by the constructions. All arguments, orientations,
measures, derivatives, alias definitions and unsummed pieces remain verbatim in M07. This
notation is an index to that output, not a new derivation. (Retrieve M06–M09.)

### Item 1: regional balances

Both engines print `N = V[Σ]`, the three `D_i = V[y_i Σ]`, and the three `P_i = V[G_i]`.
The native `Eq` versus rule definitions are joined by explicitly named lhs/rhs components.
For each rate they print the following pieces and its total. (Retrieve M07–M09.)

| Rate family | Bulk pieces, kept separate | Each inner and outer surface |
|---|---|---|
| Count: Py `count`, WL `count` | Face-minus `V[a₋ n₁]`, face-plus `V[−a₊ n₃]`; reference centre, reference radius and weight are each printed as `V[0]` | Rest flux `−F` contracted with the shell's outward normal; centre-motion transport of `Σ`; radius-motion transport of `Σ` |
| Dipole: Py `dipole.0/.1/.2`, WL `dipole.1/.2/.3` | Weighted face-minus and face-plus; moving reference `V[−Σ Ẋ_ci]` apart from weight `V[F_i]`; reference-radius piece `V[0]` | Each piece carries `y_i`: rest material flux, centre motion and radius motion separately |
| Momentum: Py `momentum.0/.1/.2`, WL `momentum.1/.2/.3` | Integrated force `V[f_i]`; reference-centre, reference-radius and weight pieces `V[0]` | Rest total momentum flux `−Π^{ij}` contracted with the shell normal, transport of `G_i` by centre motion, and transport by radius motion separately |

Every rate also prints a three-component `integratedForceCompanion`, including the count and
dipole families. It is a companion object, not a force added to their material balance. The
rate tags use Py `bulk.faceMinus/facePlus/force`, `bulk.reference.*`, `bulk.weight` and
`surface0/1.{rest,centerMotion,radiusMotion}`; WL uses `bulkSource1/2`, `bulkReference.*`,
`bulkWeight` and `surface1/2.*`. The printed source/force pieces and totals are each compared,
not inferred from agreement of a total. `Π` is an undivided supplied total flux: separately
identified stress and advected-momentum constituents are **not produced** by either engine,
because no split is supplied. Face momentum remains in `f`, not in a computed face-by-face
momentum account. (Retrieve M01–M03/M06–M09.)

### Item 2: the count centre and its derivatives

`center.relative.definition` prints the three `D_i/N` displacements. The absolute definition
prints `X_ci + D_i/N`. Thus the scripts' relative object named `X` is not the reference point
`X_c`, nor a mass centre. The relative/absolute rate and acceleration vectors retain the
change of `N`, every count/dipole rate piece and the relation to `X_c`; the absolute vectors
also retain its velocity/acceleration. `count.rateForNormalization` is the sum of the eleven
named count pieces; `count.accelerationForNormalization` is the sum of their printed time
derivatives. The complete expanded centre expressions are in M07; no count is frozen.
(Retrieve M01/M06–M09.)

For the eleven count pieces and the eleven pieces of each dipole component, both engines print
term aliases, their time-derivative definitions, every derivative-origin piece, and the direct
derivative, sum of pieces, residual integrands and route names. These are **44 rate-term
derivatives**. A volume term's derivative has `localEvolution`, bulk centre/radius motion,
and inner/outer centre/radius transport; a surface term's derivative has `localEvolution`
and that surface's centre/radius motion. The derivatives of reference-motion pieces retain
their parent names, including the separately printed moving-reference contribution.
`coverage.ungovernedRates` prints `F1,F2,F3,n1,n3,a_minus,a_plus` (WL `aminus,aplus`): their
ungoverned time rates remain terms. Only the supplied material balance governs `Σ`'s partial
time derivative in this construction. (Retrieve M02/M06–M10.)

The shared normalized-domain line applies to every normalized object, even though each line
does not repeat it. Py prints `0 < r_in < R` and `N ≠ 0`; WL also prints its angular/reality
chart restrictions. The unnormalized rates and moments retain `N = 0`. The printed domain
difference, common-domain sampling and exclusions are preserved in §3. (Retrieve M02/M07–M11.)

### Item 3: momentum received, distinct from the rate of stored momentum

`momentumReceived` is a three-vector made from the integrated `f_i` and the two **rest**
momentum-flux surface pieces. `momentumRate` prints the full regional rate, and
`momentumSurfaceMotion` prints the motion pieces separately. The shared twelve surface-motion
components are bound one by one. Py's latter payload additionally contains two `Volume(0)`
bulk-reference entries per component, which WL's payload omits; those values/keys remain
unpaired, as do their six extra Py dimension tags. The corresponding separately emitted
rate pieces are compared elsewhere. No equality of received momentum and `dP/dt` is imposed,
and no moving-frame reinterpretation replaces the supplied `x`-coordinate object. Momentum
received is not a passive factor or the force a drain exerts on another body. (Retrieve
M01–M03/M06–M09/M12–M13.)

### Item 4: count moments

There are thirteen requested moment components: the count, three dipoles and all nine second
moments about live `X_c`. Py encodes indices in `moment.()`, `moment.(i,)` and `moment.(i,j)`;
WL prints `moment.0` through `moment.12` and a separate `moment.indices`. The explicit
binding table names each object. No scalar width is defined and no Gaussian or vanishing
dipole is imposed. These count moments are not the outstanding multipoles of `M`.
(Retrieve M01/M06–M09/M12.)

### Item 5: dimensions

The engines print per-object dimensions and native sum-homogeneity operands/residuals, using
independent length, time, count and momentum units. No mass unit or count–momentum law is
introduced. Examples are `count.rate: count/time`, `dipole.*.rate: count·length/time`,
`momentum.*.rate: momentum/time`, and count-centre rate/acceleration `length/time` and
`length/time²`. The comparator prints `EmptySet` for all 662 dimension differences; the
35 common field-dimension bindings are separate exact-value comparisons. This is dimensional
homogeneity under the printed assignments, not physical delivery of an inertial law.
The native dimension-sum trees are not paired by ordinal. WL prints `Missing["zeroTerm"]`
on 42 lines; that empty zero-term audit is not a homogeneity failure or a zero-valued physical
check. Use the per-object dimensions, including `count.rate.dimensions` and
`dipole.*.rate.dimensions`, as the build acceptance requires. (Retrieve M02/M06–M10/M12.)

### Cross-engine status for supplied, geometric and requested objects

| Printed object | Filed comparator status and limit |
|---|---|
| Supplied equations: `supplied_material` (1 binding) and `supplied_momentum` (3) | All four as-computed exact residuals and every sample residual are `0`, on the unnormalized contexts stated below. This compares the printed supplied equations, not their truth for arbitrary fields. |
| Integral definitions: `integral_definitions` (3 bindings: volume and two surfaces) | All three as-computed exact residuals and every sample residual are `0`, with the stored integrands, regions and measures and the unnormalized contexts below. |
| Geometry: `geometry_surfaces` (6), `geometry_normal` (6), `geometry_areaMeasure` (2), `geometry_volumeMeasure` (1) | All fifteen as-computed exact residuals and every sample residual are `0`, on the unnormalized contexts below. Differing free-chart/split-versus-total geometry objects in §3 remain unjoined. |
| Count, dipole and momentum definitions, all rate totals, both material-face pieces, bulk force/reference/weight pieces, every surface's rest/centre/radius pieces, force companions and solved local rates | Each shared requested value has an explicit binding and exact residual `0` as computed. Domain conditions remain separate. |
| Relative/absolute count-centre definitions and rates | All six definition and six rate bindings have as-computed exact residual `0`; `N ≠ 0` applies. |
| Relative/absolute count-centre accelerations | Six nonzero-text as-computed residuals, each with resolved exact residual `0`; `N ≠ 0` applies. |
| Normalization rate/acceleration, aliases, derivative-origin pieces and partition operands/residuals | Every value is exact `0` as computed or resolves to `0`, with the per-class counts below. The domain difference is not erased, and derivative-attribution limits apply. |
| Received momentum, momentum rate and the twelve shared surface-motion components | As-computed exact residuals `0`; six Py bulk-reference values in the grouped motion payload and their dimension tags remain unpaired. |
| All thirteen moments | As-computed exact residuals `0`, with the WL-only index-label tag unjoined. |
| Per-object dimensional homogeneity and common field dimensions | All 662 dimension differences `EmptySet`; common field dimensions exact `0`. Native sum audits and additional Py alias-dimension content are outside those comparisons. |

The status in this table is retrieved per binding from M08–M11, not inferred from a successful
exit or from agreement of totals. §3 and exhaustive M09–M13 are qualifications on every row.
For the 22 supplied-equation/integral-definition/geometry bindings, M18 retrieves each binding,
operand, residual, domain context and all 45 samples. Their contexts are Py `0 < r_in < R`
and WL the same with `0 < theta < Pi` and real `theta,phi,t`; no nonzero-count restriction
applies. Samples include the count-zero stratum. The domain-predicate difference remains OPEN
as §3 states, despite zero value residuals on the comparison domain. (Retrieve M18.)

For item 2's finer printed objects, literal counts of **as-computed `0` / additional resolved
`0`** are: term-alias sides **88 / 0**; term-rate-alias sides **56 / 32**; local-evolution
pieces **3 / 41**; centre-motion derivative pieces **40 / 44**; radius-motion derivative
pieces **40 / 44**; direct partition operands **28 / 16**; partition sums **0 / 44**;
partition residual components **176 / 0**; count-normalization rate/acceleration bindings
**1 / 1** (the rate is as-computed `0`, acceleration resolves). Each member and its own
printed object text is retained, including zero reference terms. (Retrieve M07/M09–M10.)

## 3. What the accepted comparison establishes and leaves apart

### Scope and instrument limits

The filed output's printed scope is its **47-rule, 2028-binding, 91-name, two-surface**
configuration, including each binding's own object text. It prints **"No residual target was
given."** Its summary states **"zero for findings; nonzero only for operational failure"**.
Its printed resolution method is **"evaluate substitution wrappers preserving partial
derivative slots, orders and points; expand and group integrals; no supplied equation or
domain assumption"**. Its domain-residual note says **"symbolic restriction; evaluated
values are printed per sample without replacing the as-computed residual"**. These bound
the comparison; it is not a proof of whole-stream identity or a test of supplied physics.
(Retrieve M03/M08–M09.)

The printed non-comparison statements include:

- **"No source-cited value selection covers this component. Constructor and key leaves are
  not separately compared values; any unselected content remains printed. The binding table
  supplies the selected component identities and native layouts."**
- **"Missing constructor and reason-code metadata (WL 15,246) has no separate Py text
  component (Py 23-26,437). The message text itself is compared; metadata remains unpaired."**
- Geometry totals **"cannot bind individual parts"**; native dimension audits have **"no
  term/reference correspondence by tree ordinal"**.

All other literal unjoined/unpaired reasons are in M08/M12 and indexed below. The output
contains 885 `exact_comparison` rows beside 2028 `comparison` rows; the other row classes
are configuration (1), test functions (45), progress markers (66), operands (2028), binding
accounting (2028), tag accounting (63606), resources (1) and summary (1). Thus exact rows
are not 885 additional bindings. (Retrieve M08–M10/M12–M13.)

The build acceptance's limits **C1–C5** travel unchanged in substance:

- `localEvolution` is differentiation at fixed geometry. Under a motion parent it retains
  that parent's **undifferentiated velocity**. It is free of motion only under rest, flux,
  face, source and weight parents. Do not add it across parents and call that sum motion free.
- Each engine's 44 partition residuals test that the pieces sum to the derivative, not which
  piece owns a term: they share the splitter. Attribution rests on the review legs'
  independent splits, not those self-residuals alone.
- Both engines use the same jet-slot partition construction, written by the same repair
  author. Cross-engine agreement cannot detect a shared convention. The accepted disposition
  reports the independent legs' splits agreeing with both engines.
- Zero reference terms create empty WL top-level sum audits; homogeneity is still printed
  per object, with the 42 `Missing["zeroTerm"]` lines kept distinct from a failure.
- The build ran at `gpt-6-astra`; its repair at `gpt-6.1-sol`, in the same session.

The comparator acceptance also carries: each citation **contains** the lines building or
printing its object, but ranges may span neighbours; bindings rest on round-1 independent
measurements, form ablations and one-sided corruption, not citation strings; six rules use
their first row's object text, so this record uses each **binding's own** text. Its tests
are synthetic fixture controls, not locks on production tables: production-table ablations
could leave fixture tests green while moving production output. Round 2 was the user's
**citation-only repair with a scoped check**, followed by acceptance with containment, not
a new full comparison-method review. The filed output is that reviewed run itself. The
build's and comparator's historic findings/dispositions remain literal in M02/M03; none is
re-adjudicated here. (Retrieve M02–M03/M08–M09.)

In particular, a volume integrand's partial time derivative is at fixed Cartesian `x`,
then evaluated at the chart point; it is not replaced by a total derivative along the
moving chart (the comparator's repaired K1). The build's round-0 `G = 2F` probe that
stopped at `Solve::ivar` was not counted as a mathematical control; W4b replaced it.
Neither operational failure nor the existence of two engines supplies a physics check.
(Retrieve M02–M03/M06.)

### Residuals, domains and evaluation

Of 885 exact residuals, **655 print `0` as computed**, and **228 more have a nonzero-text
as-computed residual and print resolved exact `0`**. Their complete membership and both
forms remain in M09/M10. The accepted K3 resolution settles these as representation/wrapper
differences within that instrument, not permission to discard the original residuals.
The remaining two are `domains:0` and `domains:1`; their resolved residuals remain the
printed predicate symmetric differences. (Retrieve M03/M09–M10.)

Py's unnormalized geometric condition is `0 < r_in(t) < R(t)`; WL adds `0 < theta < Pi`
and real `theta,phi,t`. Normalized conditions additionally require nonzero count. The two
domain comparisons retain the nonzero as-computed/resolved predicates, their exclusive
domains and three angular probes: at `theta = −1/3` and `7/2`, Py is `True` and WL `False`;
at `theta = 1/3`, both are `True` (probe `phi = 1/4`). The common-domain condition samples
print `False` residuals. Value-binding `condition_residual` fields also retain the differing
geometric/normalized contexts even when their value residual is `0`. This is an **open
cross-engine domain difference, no owner named**; the sources do not settle whole-domain
equivalence. No angular/reality condition is silently dropped. (Retrieve M07–M11.)

The printed test forms are mixed quadratic, cubic and quartic live profiles, at
`(1/5,2/7,−1/4,1/6)`, `(-2/9,3/8,4/11,−1/7)` and `(3/10,−1/3,2/5,2/9)`.
The five printed strata are generic, centre velocity zero, both radius velocities zero,
all those velocities zero, and count zero; higher geometry derivatives stay live. The
count-zero construction subtracts the exact spatial shell mean at the probe time from
the generic live density, rather than replacing the density by a constant. The supplied
balances are not imposed on the profiles. There are **39573 evaluated sample residuals
`0`**, **90 evaluated condition residuals `False`**, and **162 `outside normalized domain`
samples**, retained per binding in M11. The accepted run record calls the evaluated total
**39663**; it records no numerical integral fallback. The summary prints
`unresolved_coverage: {}`, `numerical_integrals: []`, and `not_formed: []`; no binding has
a printed `not_compared` result. These coverage statements do not close the domain
difference or extend sampled values into excluded domains. (Retrieve M03/M08–M11.)

### Every textual difference

M09/M10 retrieve all of these, including every member of repeated classes. They are not
hidden by translating a tag/name or by replacing a difference with its total.

| Binding class | Literal differences and classification |
|---|---|
| 44 `partition_routes` labels | Py `['direct product derivative', 'partial geometric jets and normal-speed transport']` versus WL `['direct product derivative', 'chain-rule origins and window-only motion']`. Documented construction-route difference, grounded by engine sources/run record; no equality of methods is claimed. No owner named for a further route-equivalence question. |
| `schema:1556/1557` | Volume/surface descriptions: `shell integral in x` versus `volume integral; x remains a bound Cartesian argument`; `angular integral of chart-pulled integrand times area` versus `surface integral; integrand already pulled to chart`. Representational wording, settled by each printed integral definition, not by matching strings. |
| `schema:1558/1560` | `dipole/N, displacement relative to Xc` versus `dipole/N: displacement from Xc`; `concentric spheres; no field symmetry` versus `concentric spheres, not field symmetry`. Representational wording, settled by printed definitions and labelled control choice. |
| `schema:1562/1563` | `undivided supplied force density, including slab-face exchange` versus `undivided supplied f`; `total supplied flux, without a constitutive split` versus `undivided supplied total momentum flux`. Representational scope wording, settled by the supplied equations/list. No new decomposition follows. |
| `coverage:1600` | The full native `coverage.objects` lists differ. Every list member is retrieved in M07/M09 and accounted for in M12; this is a native inventory difference, not evidence that whole inventories are equal. Any further full-inventory equivalence remains open, no owner named. |
| `coverage:1601` | Py `a_minus,a_plus` versus WL `aminus,aplus` within the same seven ungoverned-rate names. Representational naming, settled by the explicit name table and sources. |
| `coverage:1602` | `no supplied separate stress/advection constituents of Pi` versus `Pi has no separately supplied stress or advection constituents`. Representational wording of the same unproduced split, settled by the list and accepted build run record. WL `Missing` metadata is still unpaired. |
| 425 `comparison_notes:1603`–`:2027` | `nothing to compare it with; no second route` versus `nothing to compare it with`. Representational engine-note wording, settled by the two `emit` definitions. These are not claims of a second derivation or cross-engine value checks inside either engine. WL constructor/reason metadata remains unpaired. |

The other text comparisons—schema physical-premise text, moment order and dimension-basis
labels—print empty differences. The 35 common field dimensions are exact-value bindings,
not text matches. (Retrieve M06–M10.)

### Every unjoined tag and unpaired component

The summary accounts for **26028 Py / 37578 WL tags**, **238488 / 580543 parsed leaves**,
and **161466 / 212319 actually compared leaves**. Bound-end parsed leaves are
**161466 / 213171**. Per-tag bound counts are **1706 / 1557**; unjoined counts are
**24322 / 36021**. Counts do not imply that a bound tag's every component was compared.
(Retrieve M08/M12–M13.)

| Unjoined class, using its stored reason | Py / WL tags | Standing |
|---|---:|---|
| Native dimension-sum audits enumerate different expression trees | 24264 / 36014 | Uncompared native audits; no ordinal term correspondence. |
| Free-radius chart and split centre/radius velocities/speeds versus surface totals | 26 / 6 | Different objects. Named unjoined geometry is Py `geometry.chart`, `geometry.centerVelocity`, `geometry.radiusVelocity`, `geometry.centerNormalSpeed`, `geometry.radiusNormalSpeed`, and WL `surfaces.velocity`, `surfaces.normalSpeed`, including their listed companion/dimension tags. Surface charts themselves are bound. |
| Extra Py bulk-reference dimensions in `momentumSurfaceMotion` | 6 / 0 | No unused WL counterpart; the grouped zero reference values themselves also remain unpaired. |
| Py unwindowed IBP intermediates | 21 / 0 | `weightedLocalRate`, `afterParts.bulk`, `afterParts.current` differ from WL's window-weighted intermediates. |
| Py engine-only Leibniz/generic pullback diagnostics | 5 / 0 | No matching WL object; WL's window differentiation is a different intermediate. |
| WL explicit moment-index labels | 0 / 1 | Py indices live in tag names. The moment values are bound, the extra label tag is not. |

These classifications are the filed comparator's source-cited different-object reasons,
not newly reconciled expressions. **Every member**, its source transcript line, occurrence
and stored leaf counts is retrieved in M12. No shared requested balance, moment, rate or
printed derivative piece is left unjoined, as the accepted run record states. This does
not cover the different named geometry or engine-local audit objects. No owner is named
for additional comparisons of those objects. (Retrieve M03/M08/M12.)

M12 also retrieves **every unpaired bound component**, with its native path, parse, leaf
count and literal reason. The exhaustive reason classes are:

- Generic constructors/keys and other unselected content (283 Py / 208 WL component records).
- Py's grouped motion payload, including extra bulk-reference values and dictionary keys
  (18 Py records). The shared four surface pieces per momentum component are compared.
- Alias-dimension and field-name-key packaging: Py prints aliases with the field table;
  WL prints its table before adding aliases (79 Py / 35 WL records).
- WL solved-law lhs, window-weighted `windowDerivative`/`afterParts` intermediates and keys
  versus Py's solved rhs/unwindowed intermediates (28 WL records).
- WL `Missing` constructor/reason-code metadata (852 WL records), with the message compared.

These are counts of stored **component records**, not asserted leaf totals. The 426 bindings
below parsed count are exactly `coverage:1602` and `comparison_notes:1603`–`:2027`; their
message is compared and WL metadata remains unpaired. M13 preserves both counts for every
binding, and M08 prints the full below-count list. No missing component is credited by
another component's agreement. This coverage remains uncompared; no owner is named for
extending it. (Retrieve M03/M08/M12–M13.)

## 4. Separate inventories and the plan boundary

The inventories required by the step list are source-row inventories, kept apart below;
their contents are retrieved literally in M04 rather than restated as S16 results.

| Inventory | Source rows (retrieve M04 for each) |
|---|---|
| Plan's S16 section | `V3_STEP_PLAN.md:670–746`, including statement/transcription note `:700–703`, gravity guardrail `:705–730`, A11 gate `:732–737` |
| Source assumptions | `../4d_1pn_full/paper/4d_1pn_full.tex:695–720` |
| Source isotropic-Gaussian width | Same source `:761–765` |
| Source reduced equation | Same source `:885–891` |
| 2.5PN statement as written | `../4d_2_5pn/paper/4d_2_5pn.tex:605–614` |
| Charter scope boundary | `CHARTER.md:74–80` and `:109–119` |
| S22 width/bridge/mass-multipole debt rows | `V3_STEP_PLAN.md:1303–1305` |
| Drain–mass gap C3 | `DEFECT_REGISTER.md:171` |
| O2 S16 interface | `directives/O2_SHARED_PHYSICS.md:607` |
| O2 scoping | `directives/O2_steady_brane_balance_scoping.md:164` |
| Three-mass guardrail | Repository-root `notes/muonium_gravity_corrections.md:21–25` |

The plan's response-side approximation, supplied mass/multipoles, compactness, smooth
external field, calibrated potential, discarded boundary fluxes and reduced-force equation
are **not printed as an S16 reduction** here. The list implements the later user decisions
with live count, count-moment and boundary-flux objects. The list author's filing leaves
`M`, `a`, multipoles of `M`, compactness, source normalization and any flux suppression
as debts for S22's table.
`a_WT` remains distinct from `a_mouth`; neither control radius nor a count second moment
delivers their A11 bridge. S16 supplies no smooth surroundings field to expand, no
Gaussian profile width and no source normalization. (Retrieve M01/M04/M07.)

The plan's three-mass/species guardrail and its deferred branch remain constraints on
later work, not tasks run here. S19 owns the C10 tail question; S22 owns the core and the
width bridge handoff. S15, conditional on S14a, owns the source-side question of what force
a drain exerts. `R-S14a-01` remains OPEN for slow/far selection with its existing joint
deliverers. O2's orbital-`GM` matching interface remains an inferred assignment; it does
not assign the missing surroundings field's delivery or establish a coupling. Where these
sources do not assign an owner, this record says **no owner named**. No broader charter
scope or later-step result is earned by the regional identities. (Retrieve M01/M04/M14.)

Two plan handoffs remain beside this step's limited delivery. The plan says **"S16 must be
run BEFORE S19 is trusted"** and **"its result feeds back into the charter"**, adding
**"the scope boundary is materially narrower than the charter claims"** as a likely outcome
and **"Say so rather than preserving the boundary."** Whether this regional-object step
disposes of the S19 trust gate remains OPEN; **no owner named** for that disposition in the
cited rows. S19's separately named C10/tail ownership is retained. The charter-feedback
handoff also remains OPEN; **no owner named** for its disposition in those rows. No gate
verdict, charter amendment or scope-boundary decision is made here. (Retrieve M04/M17.)

## 5. Prior art as an oracle, and the second route

The accepted repository-root gravity section, headed **"ATTRIBUTION PASS — the gravity
sector and S16"**, attributes the region algebra to G1/G2 and places G3 beside the centre
relations. Its row labels are retained below. The sources file records what was opened;
some citations were available only as abstracts or bibliographic records, and secondary
claims remain secondary. Its round-2 acceptance requires checking computed objects only
where the row's conditions match S16's. No literature is newly searched or verified here,
and no row supplies a premise for a construction. (Retrieve M05.)

The table quotes each cited row's conditions beside the list's S16 conditions. A comparison
here means retrieval of named objects/standing only; it is not a new algebraic residual.
Where conditions differ, or comparison requires computation, the result comparison stays
OPEN. No owner is assigned to a new prior-art applicability/equivalence check.
(Retrieve M01/M05–M09.)

| Row and its filed label | Prior conditions quoted from the gravity row | S16 conditions quoted from the step list | Retrieval comparison / open comparison |
|---|---|---|---|
| G1 **STANDARD** | Reynolds: transport "relative to the surface's own velocity"; "The instantaneous rate of production within the surface is not altered by the continuous motion of the surface". Variable-mass extension: "sources inside a material volume, or mass flowing through the surface of a non-material volume". "Splitting a sphere's surface velocity into the centre's motion and the radius change was not looked up as a separate result." | "both surfaces move"; "the part present when that surface is at rest"; motion in "two pieces" | Named transport/rest-motion objects are present in M07 under the general transport conditions: attribution to standard algebra by retrieval. No full prior-formula subtraction or separate historical split theorem is computed. |
| G2 **STANDARD** | "Lamb and Irschik–Holl balance the momentum of matter, mass × velocity (Reynolds states the theorem for any entity), whereas S16 relates G^i to nothing" | "G^i is not related to F^i or Σ by any supplied law" | Any-entity transport object and stored/rest/motion distinction are identifiable by retrieval. The mass-velocity specializations are not identified with S16; a specialized equality would need extra identification/computation and stays OPEN. |
| G3 **ADJACENT** | "a mass density whose momentum is ∫ρ × velocity"; "§4 excludes sources inside the volume, and §3 treats them only for a material volume"; "one control surface moving at a velocity u, with no split into a centre's motion and a radius change". GR alternative: "post-1-Newtonian, vacuum field equations between the bodies, moments from the metric, no drain" | "No mass is defined"; "both surfaces move"; "G^i is not related to F^i or Σ by any supplied law" | Located centre-rate/acceleration precedent, not S16's result. Conditions differ; comparison OPEN. |
| G4 **ADJACENT** | "curved spacetime, Papapetrou's test particles in a given gravitational field, no drain"; the constraint summary is "secondary only" and Dixon's centre statement is from "an opened abstract" | "the material count's moments about X_c through second order"; "No mass is defined" | Reference-point/moment analogy only. No spin/centroid theorem is transferred; result comparison OPEN. |
| G5 **ADJACENT**, earlier question | "curved spacetime, and vacuum field equations outside the body"; "point singularities (EIH) or slowly changing moments of an isolated body (Thorne–Hartle); no drain" | "no field of the surroundings"; "no throat solve"; "Time dependence is kept" | Earlier problem-of-motion question, not a result for this regional object; conditions differ, comparison OPEN. |
| G6 **ADJACENT**, earlier question | "only the transverse, nondissipative part of the force; a line defect carrying circulation; an infinite, homogeneous, neutral superfluid with a complex scalar order parameter; a uniform background; no drain" | "a region around one throat"; "All three are general live fields of x and t" | Different force/defect and conditions; no core-independence result is transferred. Comparison OPEN. |
| G7 **KNOWN**, prior result for the open passive requirement | "inviscid and incompressible, and Lagally's flow is steady". Wang adds "irrotational, barotropic, a uniform constant density, no external body forces on the fluid", "an unbounded fluid at rest at infinity, a point sink, and a control volume between two spheres, fixed in space and then shrunk onto the sink. All 3D, no sheet". Landweber–Yih covers "Unsteady flow and moving bodies"; Landweber–Miloh "Multipoles and deformable bodies". | "It computes no passive factor"; "a named field of the surroundings, a rule separating it from the body's own field, and its normalization"; "no steady, quasi-static or slow limit is taken" | Name G7 as the oracle for `R-S16-01`, retaining Wang's no-external-body-force condition for Theorems 2/3 and Corollary 4. S16 has no passive object to compare. No identification with Lagally/Wang force; comparison OPEN. |
| G8 **ADJACENT** | "inviscid flow"; body radius "small compared with the scale over which the strain rate changes"; for Auton's sphere, "the changes in the ambient velocity field over the scale of the sphere are small compared with the velocity of the body relative to the flow"; "volume and added-mass coefficient, not by a drain" | "Nothing is imposed: no compactness"; "no inertial factor" | Non-draining-body oracle only; conditions differ, comparison OPEN. |
| G9 **KNOWN**, later-step oracle | "incompressible potential flow"; "Lamb's sources are the ends of narrow tubes, so the fluid is conserved overall"; "Bjerknes' forces are averages over synchronous pulsation". Wang's inverse-square-only force is for "a sink at rest"; the moving result retains "a velocity term ρQ₂v₂". | "no force law between bodies"; "no potential" | S16 supplies no pair-force object. No inverse-square result is claimed; comparison OPEN for later work. The Bjerknes charge analogy and its reversed sign are not gravity premises. |
| G10 **KNOWN**, an unexecuted dated sketch | "a 3D substance, with no sheet or throats; inertia-free; light comes from line-element resistance, not MacCullagh's rotation; no deflection computation was found in either fragment". Dated fragment has mass-proportional inflow in text and force-proportional inflow in a footnote; the undated fragment keeps amount "constant" and states no atom inflow. | "No mass is defined"; "The face terms are the material leaving through the faces" | Architectural ancestor as filed, not a reproduced S16 result. Distinct fragments are not fused or jointly dated; result comparison OPEN/inapplicable without additional structure. |
| G11 **KNOWN** | "Particles are sinks in an ideal-fluid Ω(0) substratum"; "A single-author line of work; no independent assessment was found". River/superfluid-inflow works in that row are labelled "ADJACENT". | "Each of the three stays a general live field, with no material law chosen" | Filed sink-gravity precedent, not a force-law result of S16. No prior model is selected; result comparison OPEN. |
| G12 **KNOWN as posits**, not derivations | Riemann "posits an inflow ∝ mass" (and force in the footnote); Wang's "Assumption 6 and dm/dt = −ρQ make mass the accumulated inflow, increasing linearly"; Yarkovsky is "secondary only" | "Count material; mass open"; "No mass is defined" | C3's mass-link gap is not closed by a historical posit. No mass derivation to compare; comparison OPEN. |
| G13 **KNOWN**, objections | "A moving drain feels a drag: Wang's formula (32) gives a universal damping force"; "Energy disappears without explanation, some drag must arise, and no process that creates matter is known (Wikipedia, secondary only)"; "every pulsation in the universe would have to be in phase (same source)" | "no passive factor"; "no mass"; "no material law, drain law or flow" chosen | Each objection is retained; S16 computes no energy/mass-creation/pulsation model or sink drag law. Applicability and whether the model answers any objection remain OPEN. |

Thus the retrieval comparisons made are identification of standard transport objects and
their stated domains for G1/the any-entity part of G2, and explicit checks of the remaining
rows' standing and conditions against S16's packet. No prior-art algebraic equality is
newly established. The gravity section's NOT FOUND claims remain search-failure statements,
not originality proofs; this record claims no originality. (Retrieve M01/M05–M09.)

For the register's **second route**, G13's moving-drain damping remains beside the passive
and inertial questions; unexplained energy disappearance remains beside S12 exchange and
the existing S22 energy-split requirement; matter creation remains beside C3/the missing
mass account; universal pulsation phase remains an unapplied-model question. S16 answers
none. The relevant future interfaces are not owners of a newly invented literature check:
none is named for that check. No energy-conservation, matter-creation or synchronization
premise is added to the register as if S16 had reproduced that model. The same applicability
question is retained for G3–G12 with their different objects/conditions. This examines every
identified gravity row under the register's broader rule; it does not apply only the
historical reproduced-result restriction. (Retrieve M01/M04–M05/M14.)

## 6. Requirements filed and the S16 register pass

The register's schema and rest-on method govern. A requirement names a future object; a
postulate is assumed here, and a defect is already wrong. The user's selected options
expressly file the passive-factor needs and inertial relation. The list author separately
files compactness (E3/F3), face momentum (E2), density response (O2 `:138`) and the S22
debt account. The current contract, item 7, directs this record to file the list's named
objects. Their admission grounds are the distinct sourced future-delivery conditions below,
not attribution of every filing to the user or a claim that the unsolved regional balances
require a selected material state. A future owner alone does not admit an unspecified
closure form. (Retrieve M01–M02/M04/M14/M17.)

| Entry; status | Object filed once | Target; source/failure boundary |
|---|---|---|
| `R-S16-01`; OPEN | Named surroundings field, separation from the body's own field, normalization | **No owner named**. User's 2026-10-21 "Record what it needs": "What one needs … goes to the register as a requirement" (list `:49–53`; pass-1 decision). Admitted on that express passive-factor requirement; without it received momentum cannot be given that interpretation. |
| `R-S16-02`; OPEN | Body compactness relative to both separation and surroundings-field variation scale | **No owner named** for physical delivery; S22 receives the debt. Orchestrator's list filing `:158–161`, "The body's compactness, separately", from E3/F3 and source assumption `:697–700`. Admitted as the named small-body condition of the later response approximation, not a user-selected or imposed compactness premise; count moments settle neither comparison. |
| `R-S16-03`; OPEN | Momentum exchanged at each slab face, split by the faces' motion | **No owner named**. Orchestrator's list filing `:162–163`, "For the faces", from E2's direction to keep "the face-resolved momentum accounting OPEN". Admitted as that explicitly retained missing accounting object; undivided `f` and regional control-surface motion do not deliver it. |
| `R-S16-04`; OPEN | Momentum–material relation needed for an inertial factor | **S8**, with S1.5 antecedents when available. User's 2026-10-21 "Count material; mass open": "The law it needs goes to the register, as O2 keeps it (S8)" (list `:54–56`, O2 `:133`). Admitted on this express future inertial-relation requirement; no inertia follows from printing `P` beside the count centre. |
| `R-S16-05`; OPEN | Slab density response | **S8** substrate reduction / **S22** completion. Orchestrator's list filing `:164–166`: "The slab's density response is also OPEN there", citing O2 `:138`. Admitted as the list's expressly carried future density-response delivery, required to give that interpretation to material/mass accounting; retaining live `Σ` does not supply it. It is not admitted merely because O2 names an owner. |
| `R-S22-02`; OPEN | S22's missing reduction-data inventory: `M`, `a`/`a_WT`, multipoles of `M`, source normalization, any boundary-flux suppression, and the separately identified `a_WT`/`a_mouth` bridge | **S22 for the debt-table account**, not an invented physical supplier. Orchestrator's filing `:28–30`, "What stays a debt for S22's table is what this step does not deliver", and plan A11/debt rows. Admitted as the named debt-account delivery; compactness is cross-referenced to `R-S16-02`, not filed twice. Printed count/moments/fluxes are not these debts. |

Each entry has the register's source, target, requirement, on-failure and OPEN fields.
None names a kinetic/material/drain law, candidate, premise or means of meeting the
requirement. The momentum–material **relation required** is not a prescribed law for
obtaining it. The missing physical compactness owner is not inferred from S22 being
the debt-table recipient. Nor is the passive-field requirement's normalization conflated
with source normalization. (Retrieve M01/M04/M14.)

The rest-on test also admits two inherited conditions of the present conditional object.
The three in-plane coordinates and spherical shell geometry consume the existing
`R-S1-01` dimension obligation at **S6**; without delivery/reaffirmation of that coordinate
setting the printed objects do not transfer unchanged. The projected-material object's
shear-holding-region interpretation inherits part 1's flagged, adopted U1 support
identification and consumes `R-S1-02` at **S1** for that interpretation only. Failure of
that identification removes the interpretation, not the algebra. No dimension, shear
response, confinement or constitutive law is derived; U1's re-ordering clause remains
unenforced. Both merge into existing entries, as the earlier passes did. (Retrieve M01/M04/M14/M17.)

Together with `R-S12-02` for the supplied projected-material object/"Keep S12" only and
`R-S21-01` for the framing decision only, these make **four** existing entries gaining
scoped S16 provenance. Their requirement, target and status do not change. The original
homogeneous quadratic/slab inertia in `R-S8-06` is not the new live momentum–material
relation; O2's earlier no-widening boundary stays. No density-response entry existed.
S22's existing `R-S22-01` is an energy split, not the reduction-debt account. No new entry
duplicates those objects. No source is added to `R-S14a-01`; its source-step id does not
name its deliverer, and a note now makes that explicit. (Retrieve M14/M16.)

No further requirement whose object is a printed S16 value is admitted. These are general
conditional balances, definitions and moments; they select no material state or vanishing
flux. The inherited coordinate/support dependencies are covered above; the supplied projected-material/S12 scope
is already covered by `R-S12-02`. No printed nonzero flux is converted into a requirement
to suppress it. The unowned general `G,Π,f` forms remain open handoffs, not extra entries
merely because a closed solve might choose them. The core is excluded, so its future
S22 contents are not premises of these integrals. Control-sphere shape, second-order
moments and independent dimension units are labelled choices, distinct from those inherited
coordinate/support conditions. The separately sourced future requirements above concern
their named interpretations/accounting deliveries, not defects in the conditional regional
objects. (Retrieve M01/M06–M07/M14/M17.)

The pass applies both population routes; §5 and the pass retain every prior-art condition
and historic objection without assuming its applicability. It adds six entries:
**19 → 25, all OPEN**. New by target: **three no physical owner named**, **one S8**,
**one S8/S22**, **one S22 debt-table account**. Exactly **four** existing entries gain
source/provenance: `R-S1-01`, `R-S1-02`, `R-S12-02` and `R-S21-01`. The Y1–Y4 repair is
**25 → 25**, with no new entry and two additional existing-entry S16 sources. No status
changes or requirement closes. The original before-count is the register at `4d4d22cb`,
not an estimate from the decision list. (Retrieve M14/M16–M17.)

## 7. What stays OPEN, and source disagreements

Every retained-OPEN item of the list remains open with its sourced owner:

| Item | Owner / standing (retrieve M01/M04/M14) |
|---|---|
| Forms of `G^i`, `Π^{ij}`, `f^i` | No owner assigned here; not inferred from the O2 matching interface. |
| Slab material law | No owner assigned by S16; accepted S12 gives material completion to S8/S22. |
| Faces' law | No owner named. |
| Inherited real-fraction/A13 qualification | S1's gate, propagated through S5/S12; no phase-carrying extension is computed. |
| Surroundings field/separation/normalization, compactness, face momentum/split | `R-S16-01/-02/-03`, no physical owner named; compactness also in S22's debt table. |
| Momentum–material relation and density response | O2's S8 and S8/S22 targets, `R-S16-04/-05`. |
| Core interior, count and momentum | S22. |
| `a_WT`/`a_mouth` bridge | S22, under A11; no radius substitution. |
| Slow limit and far behaviour | `R-S14a-01`, jointly S12 source/flux forms, S8/S22 material law, unowned face/interior laws. This step takes no slow limit. |
| Source side: what force a drain exerts | S15, conditional on S14a. |
| Drain-rate/mass link | C3, a gap. No delivery owner is named by that gap row; it is also retained in the S22 debt account without claiming delivery. |
| C10 tail question | S19, not run here. |
| Plan's S19 trust gate: "S16 must be run BEFORE S19 is trusted" | OPEN handoff; **no owner named** for disposition of the gate in the cited plan rows. No gate decision here (M04/M17). |
| Plan's charter feedback: "its result feeds back into the charter" | OPEN handoff; **no owner named** for disposition in the cited plan rows. No charter/scope-boundary decision here (M04/M17). |
| Species-dependent response branch | Plan's deferred branch, gated on the throat/moving-throat PDE and S11 medium chain; not run here. |
| Missing `M`, profile width, mass multipoles, source normalization, flux suppression | S22's debt table; physical suppliers remain unassigned where the sources do not name them. |

The following disagreements/qualifications are preserved with both sides, not repaired by
new calculation:

- **Who files the requirements.** The orchestrator's heading says **"Requirements the
  user's decisions file"**. The selected passive option instead specifies **"What one
  needs (a named field of the surroundings, a rule separating it from the body's own
  field, and its normalization) goes to the register as a requirement"**, and the
  selected mass option says **"The law it needs goes to the register, as O2 keeps it
  (S8)."** Those option texts file the two user-selected objects. Compactness, face
  momentum and density response are the orchestrator's separate list filings at
  `:158–166`, from E3/F3, E2 and O2 `:138`; the S22 debt account is its filing at
  `:28–30`. The heading's broader attribution disagrees with the option texts. Each
  entry retains its own filer and admission ground in §6/the register; none of those
  three additions is relabelled a user decision. (Retrieve M01–M02/M04/M14/M17.)
- **Plan versus the later user-selected object.** The plan says **"Rewrite S16 as a
  conditional response-side monopole approximation"** and calls `M`, `a`, multipoles,
  compactness, source normalization and boundary-flux suppression **"supplied inputs"**.
  The list says **"No mass is defined"**, **"No source normalization is used, because
  no potential is supplied"**, and **"The count, its moments and the boundary fluxes
  are printed objects, not those debts"**. This record follows the later user decisions
  and does not say the monopole approximation is delivered. (Retrieve M01/M04.)
- **O2's earlier interface versus current scope.** O2 `:607` says **"Its compact-body,
  supplied external-field/worldtube approximation retains calibrated potential, supplied
  mass/multipoles and discarded boundary-flux assumptions"**; its S16 designation is
  **"an inference"**. S16 says **"no mass"**, **"no field of the surroundings"**,
  **"no potential"**, with boundary terms live. Neither interface text supplies those
  missing objects here or names their deliverers. (Retrieve M01/M04.)
- **Historical mass-link filing.** Accepted S12's OPEN table says **"Drain-strength/mass
  link"** with **"S13"** as owner; the S16 list retains **"C3, a gap: the link between
  the drain rate and mass"**, and the C3 row names no delivery owner. S13 is folded into
  S14a by the user's later filing decision, as the accepted status/provenance states.
  These are different source assignments, not a delivered mass bridge or a new owner
  inferred here. The current S16 list's C3 gap and S22 debt-table recipient are retained
  with that historical qualification. (Retrieve M01/M04/M15.)
- **Revision heading.** The source still starts **"Revision 3."**, while its later
  paragraph says **"Revision 4 repairs pass 3's findings F1–F3"** and the acceptance
  clears revision 4 at `0464b80a`. The accepted unchanged content governs; the source
  heading is not edited. (Retrieve M01/M02.)
- **2.5PN quotation versus 1PN equation.** The former prints
  `+O((a²/r²) Φ_ext/r)` in a force equation; the latter prints
  `+O(M (a²/r²) Φ_N/r)`. The plan calls the missing mass factor a **"typo"** and
  instructs **"Do not present the repaired formula as a quotation."** The inventory
  retrieves both as written. No repaired equation is presented as S16's output and no
  dimensional calculation is newly performed. (Retrieve M04.)
- **Derivative-partition readings.** The build disposition records Grok's demand that
  **"localEvolution hold no velocity factor"**, versus Claude/the adopted reading that
  each derivative of a motion parent **"is kept apart and filed in its piece"**. The
  disposition confirms the fact, not a defect, and carries C1. This record retains
  undifferentiated velocities within those parents, not a motion-free sum. (Retrieve M02.)
- **Citation scope.** Round 2 records Grok's range objections as **"FACTS CONFIRMED;
  NOT A FINDING"** and permits only **"Each printed citation contains the lines that
  build or print its object"**. It forbids claiming **"That every line a citation
  names belongs to that object."** The user accepted containment. (Retrieve M03.)
- **Population route and no-commit wording.** The register's historical disagreement
  quotes **"Run it only on a prior result that a record in the reading list identifies
  as reproduced by the sector"**, whereas its method says **"Run this route for every
  sector with identified prior art"**. The broader register route is used. The current
  contract requests **"a committed generator script"** and also says **"No commit
  and no push"**. The reviewed generator and measurements are preserved at `4ae7e6cf`;
  this repair is supplied uncommitted under the no-commit instruction. No operative
  register/list conflict is left by that choice;
  the conflicting wording is reported, not silently satisfied. (Retrieve M01/M14–M17.)

The two unresolved domain predicates, native inventory comparison and all unpaired/unjoined
coverage in §3 remain as measured. No owner is named for extending their comparison. Whether
two different domains or inventories can be treated as equivalent beyond those filed
comparisons needs work beyond retrieval and stays OPEN. So do new algebraic comparisons
with prior art, matched-condition passive response, a field split/normalization, inertial
identification, density response, face momentum resolution, compactness, the slow/far
selection, source matching, any mass bridge or throat/species/tail solve. None is attempted
as a prerequisite to this record. (Retrieve M01/M04–M05/M08–M14.)

## 8. Account of edits and authoring STOP

Outside this record, its measurements and generator, the following are **all edits by this
work** against `4d4d22cb`; M16 retrieves their exact `git diff`. No physics or inventory
source, accepted engine, comparator or transcript changes.

| File | Changes |
|---|---|
| `SUBSTRATE_REQUIREMENTS.md` | Updates S14a part 2 header/pass heading and generator/measurement provenance to acceptance at `79511272`; no other part-2 pass content changes. Adds the source-step-id note for `R-S14a-01`. Adds S16 status/census summary and population pass with both routes, merge/rest-on applications and limitations. Adds six entries and six by-target rows for §6's objects, all OPEN; updates the two current total-count lines from nineteen to twenty-five. Each new entry names its filer, source quotation/line/finding and admission ground: selected user options for passive/inertial needs, list-author filings for compactness/face momentum/density response/debt account. Reports the heading/option-text disagreement. Adds inherited coordinate-dimension/support-identification S16 sources and scope/rest-on paragraphs to `R-S1-01`/`R-S1-02`, supplied-projection/Keep-S12 provenance to `R-S12-02`, and framing provenance to `R-S21-01`, including all four by-target rows. Records four existing-source gains overall and 25→25 for Y1–Y4, plus preserved-review/no-new-commit provenance and the OPEN S19-gate/charter-feedback handoffs; no existing requirement/target/status changes. |
| `steps/S14a_2_slow_limit.md` | Header and starting-commit/repair lines state acceptance at `79511272`. Generator/measurements/repair provenance states their subsequent commit and acceptance there. A provenance sentence makes the quoted historical pending-review/commit account refer to the earlier authoring pass. Closing authoring/review/STOP provenance is historical and states subsequent acceptance instead of review still to come. No other part-2 content changes. |
| Repository-root `STATUS.md` | Heading and part-2 paragraph/closing register paragraph state S14a part-2 acceptance at `79511272`. Adds `R-S14a-01` to the open-question table with its exact joint owners. Adds S16's accepted construction/comparator versus written-for-review record/pass status, conditional regional claim and limits, 19→25 OPEN census and no closure. Distinguishes the user-selected passive/inertial filings from the list-author filings; names the four existing entries gaining scoped S16 provenance and 25→25 repair census. States preservation at `4ae7e6cf` with repair awaiting review. Adds two open-question rows for S16 requirements, distinguishing missing physical owners from S22's debt-table recipient, and separately carries the plan's S19 trust gate and charter feedback in prose and two open-question rows, each OPEN with no disposition owner named and no decision made here. |

The round-0 disposition and lookups, `directives/_measurements/S16_record_r0_review_disposition.md`
and `S16_record_r0_review_lookups.md`, are read-only orchestrator additions preserved at
`4ae7e6cf`, not edits by this author. M17 retrieves them and their commit provenance;
their filenames also appear in M16's whole-tree diff inventory against `4d4d22cb`.
The pre-existing `.claude/skills/prep-for-compact/SKILL.md` modification is unrelated and
untouched; its filename is visible in M16's whole-tree edit inventory and is excluded from
this author's account above. No orchestrator memory or other checkout was read. Generator
working files are under `_scratch/s16/rec/work/`; no other job's scratch files are changed.
The final generator regenerates the measurements after the deliverables and ancillary
edits are written. No review is launched, no commit is made and nothing is pushed. The
record's source-first review remains with the orchestrator at this authoring STOP.
(Retrieve M01/M16.)
