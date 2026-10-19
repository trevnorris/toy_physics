# S9b — light bending and radar delay: Parts A–C

**Record author:** Codex, 2026-10-19; record review round 1 revision. **Standing:** revised for
the record STOP from `36f59ff1`; this revision awaits review. No new commit. Part D is held, not built or accepted.

The accepted constructions give conditional optical functionals and profile requirements, not a
substrate-derived light equation or a demonstrated match to GR. At the selected comparison grades,
M13 retrieves byte-identical paired condition payloads in PY. In WL, `GM` and `ProfileEquation`
match literally but `InverseWeight`, `ReferenceInverse` and velocity-reduction fields differ (§3).
Whether those differing fields state the same condition is OPEN. The record therefore establishes
no second independent test from the PY pairs and leaves the WL separate-test count OPEN; it does
not collapse both engines to one test or infer independence from differing text.
Their complete cross-engine conditions, implied exchange and effective γ objects remain uncompared.
The numerical comparisons that do form retain their printed residuals and unresolved coverage (§6).
Sources: `B_CONDITION_*`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_*`, measurements M3–M7/M12–M13.

## 1. Sources and retrieval contract

Paths below are relative to `research/pde_ledger_v3/`. The governing sources are
`directives/S9b_record_directive.md`; Parts A–C of `directives/S9b_SHARED_PHYSICS.md` at `34a994d5`;
D2/D5 of `directives/S9b_repair_decision_list.md` at `79055918`; and O2 §8 at `72866fcf`.
The build is accepted at `fb2d6d25`, with engine sources unchanged since `503736ff`:
`scripts/S9b_light_bending_sympy_audit.py` and
`mathematica/S9b_light_bending_mathematica_audit.wl`. Their production transcripts are the corresponding
`out/*_audit.out` files; `directives/_measurements/S9b_production_runs.md` identifies their runs.
The comparator is accepted at `88bb1f41`; its filed, reviewed production output is
`scripts/out/S9b_cross_engine_comparator.out`, with
`scripts/out/S9b_cross_engine_comparator_accounting.jsonl`. The run record is
`directives/_measurements/S9b_comparator_production_run.md`. Its r4 acceptance and the build's r11
acceptance, with the preceding dispositions, bound every claim here. No engine, harness, comparator
or oracle was rerun. Harness evidence concerns its own engine only. (M1/M8/M12.)

**Retrieval commands.** Every `M#` citation in this record means the command
`python3 research/pde_ledger_v3/scripts/S9b_record_measurements.py --section M#`, run from the
repository root; its literal stdout is in
[`_measurements/S9b_record_measurements.md`](_measurements/S9b_record_measurements.md).
M3/M4 retrieve original WL/PY field substrings with native tags, line numbers and syntax paths;
M2 preserves general gates; M11 preserves the entire WL stream, including restricted/response/forward
gates. M5 retrieves comparator configuration, comparison fields, literal counts and every numeric
sample whose `computed_18` or `computed_30` string differs from `(0.0 + 0.0j)`; no tolerance is applied.
M6 is the complete accounting; M7 is the accepted report, full coverage appendix and name table;
M8 preserves the acceptance history; M9 is O2; M10 is the register census/pass; M12 verifies literal
source-version matches and retains the PY γ-difference solution sets. M13 compares literal within-engine
condition payloads and fields; M14 retrieves accepted harness marks/coverage; M15 preserves the
repository's prior-art sources; M16 retrieves complete engine-specific mass and N-dependence objects.
Every source read is pinned with `git show REV:path`; annex targets are resolved from that revision
inside this checkout. Native/harness sources use `fb2d6d25`, comparator sources `88bb1f41`, and
repository context `36f59ff1`, except the explicit spec/decision/O2 pins above. M10 reads the
pre-S9b register explicitly at `caa0710d`; its after-pass text is labelled WORKTREE. Advancing HEAD
alone cannot change these measurements. Full native objects can also
be retrieved by `--native py TAG` or `--native wl TAG`. Compact listings never replace attached domains.
The final generator is `scripts/S9b_record_measurements.py`; regenerate with `--write`.
The current instruction to commit nothing governs the directive's description of a “committed generator”.

## 2. Supplied object, domain and standing inputs

The supplied local leading-eikonal dispersion is

```
(ω − Vⁱkᵢ)² = cγ² gⁱʲkᵢkⱼ,
gᵢⱼ = δᵢⱼ + ∂ᵢξw ∂ⱼξw,       cγ² = μ⊥/ρbr,       ξw = ℓh.
```

`u` is the material displacement of the material with density `ρbr` and velocity `V`; the same
material carries the optical branch and the supplied coordinate mass balance
`∇·(ρbr V)=−jn`. The polarization identification gives one isotropic `cγ` to both transverse
polarizations. `cγ` is the nonnegative root. `ℓ` is a fixed reduction scale, not a selected thickness.
`c₀=lim cγ` is identified with measured asymptotic light speed; bulk `cs` is separate. `GM` is the
independent orbital mass parameter; no amplitude, `Φ` or density profile is identified with it.
All these are supplied identifications, not findings. LAB_HELD supplies neither a physical holder
nor a material-reference evolution law. (Spec, M1; native supplied fields, M11.)

The domain is one isolated spherical mass at rest, Eulerian steady draining profiles, far field,
linear waves, lab time in the far-field rest frame. Profiles are general radial functions, not a
selected power family. The optical box retains every `δᵃ(V/c₀)ᵇ((∂ξw)²)ᶜ` with
`0≤a≤1`, `0≤b≤2`, `0≤c≤1`. The twelve grades are G000, G001, G010, G011, G020, G021,
G100, G101, G110, G111, G120, G121, with supplied ε powers
`0,1,1/2,3/2,1,2,1,2,3/2,5/2,2,3`. Here `ε=GM/(c₀²r)` and the supplied counting is
`δ=O(ε)`, `V/c₀=O(ε^(1/2))`, `(∂ξw)²=O(ε)`. Only G100/G010/G020/G001 are summed in
Part B. The other retained grades remain printed, not compared with a first-order reference.
The counting is owned by gravity/S12, not derived here. (M1, M3/M4 retained-grade fields.)

The native gates require real coefficients, nonzero material density, a real positive optical
ratio `μ⊥/ρbr`, its equality to `c₀²(1+δ)²` with the nonnegative speed, and the ray's own turning
and angular-closure construction. Traversal additionally requires the printed along-ray inequality
`μ⊥/ρbr > V²(1+ξw′²)`. The turning roots, outward monotonicity and finite-endpoint closure stay
attached; a pointwise inequality alone is not the full domain. WL uses its turning-radius/impact
conditions; PY uses its constructed ray-root sets. Their equivalence is not established by the
comparator. PY declares `rho_br_profile` positive; WL's gate states real, nonzero `rhoBr` and the
ratio condition. This literal assumption/domain difference remains OPEN, no adjudication owner named;
it is distinct from the bulk-density partition below. (M2–M4/M11.)

WL's condition domain additionally prints `r>bFar`, `rAnchor>bFar`, required derivatives, convergence
of the displayed improper integrals and existence of the Abel inverse on the exterior half-line,
with `δ=O(1/r)`, `V/c₀=O(1/√r)`, `ξw′²=O(1/r)`. PY attaches the exterior interval and its
ray gates to its radial identity. Branch-type rows retain growing/decaying, absent/nonpropagating and
failure-to-traverse alternatives; outside the observable gate the output is `NOT_ESTABLISHED`.
These are conditional classifications, not a proof that any real substrate supplies an admissible ray.
(M2, M3/M4 `BRANCH_TYPE`, M11 `B_CONDITION_*`.)

Part A prints full flyby turning, full finite-endpoint round-trip excess time, both one-way excess
times and half their difference, grade by grade. Part B compares only
`𝒮RT=lim ∂ΔtRT/∂ln(1/b²)` with `ZE,ZR` fixed during differentiation and then far endpoints.
Endpoint-dependent terms independent of `b` are not compared. The supplied references are
`Δθref=(1+γ)2GM/(bc₀²)` and logarithmic round-trip slope `2(1+γ)GM/c₀³`, with `γGR=1`.
They are comparison inputs, not a demanded construction outcome. The nonreciprocal path-dependence
predicate prints `False` on its domain in each engine: endpoint-only in this radial setting, not
zero one-way nonreciprocity in general; outside that domain it is `NOT_ESTABLISHED`.
(M1; M3/M4 `A_NONRECIPROCAL_PATH_DEPENDENCE`; M5/M6.) The build's K11 limit applies to
all these shared nonreciprocal objects, including conditional forward copies: “The SymPy harness
evaluates 3 tags and prints 376 `NOT_EVALUATED`. The Wolfram harness reaches 10 tags and
prints 364 `NOT_EVALUATED`. Neither reaches the shared `A_NONRECIPROCAL_*` tags. Claim limit:
K11 supports claims about the local one-form and closedness construction only, not item 8's
shared observables.” (Build r7, lines 72–74; unchanged through r11, M8/M14.)

The accepted harnesses print baseline/corrupted shared objects for K1–K10 (including K4a/b/c,
K9a/b), with literal structural marks. This is emitted-object reach, not a witness at every object.
For the 13 general `A_NONRECIPROCAL_*` tags, the marks are:

| Knives | PY | WL |
|---|---|---|
| K1, K2, K3, K4c, K10 | all 13 `not compared` | all 13 `not compared` |
| K4b | all 13 `not compared` | all 13 `identical` |
| K5 | G010/G011/G110/G111 `not compared`; other 9 `identical` | same named mark partition |
| K4a, K6, K7, K8, K9a, K9b | all 13 `identical` | all 13 `identical` |
| K11 | all 13 not evaluated; difference/residual pairs `not compared` | all 13 not evaluated; difference/residual pairs `not compared` |

These marks describe printed reach and literal structure, not a mathematical change verdict.
K10's WL witness is specifically the shared path-dependence **domain**, with corrupted/baseline
`{False, True}`; PY's K10 witness is instead at `A_DEFLECTION_G000`, not a nonreciprocal value.
M14 retains every mark separately, including the 24 conditional forward grade objects: K8 marks
those `not compared`, though its 13 general tags are `identical`; other engine-specific forward
partitions are retained there too. No local K11 witness validates a shared nonreciprocal observable.
The K11 quoted limit also governs §6's comparator rows; cross-engine samples do not enlarge
that knife's reach. (M8 r7–r11; M14.)

## 3. Each mechanism alone

This section precedes all combined conditions. The notation in the tables transcribes printed
fields: `r` is the integration/radial coordinate, primes are its derivatives, `ξ=ξw`, and
`D=√(r²−b²)` abbreviates the displayed denominator. These are readable quotations of emitted
integrands, not new reductions. All integrals run from `b` to infinity on §2's restricted ray and
functional domain. The restrictions are substitutions into the general objects and their gates;
they add no premise. In particular, flow-only here does **not** assume `jn=0`.

| Restriction | Printed deflection contribution | Printed radar-slope contribution | PY printed condition, for every exterior radius |
|---|---|---|---|
| `R_FLOW_ONLY_`: δ=ξ=0, V live | G020: `∫ b V(V−2rV′)/(c₀² r D) dr` | G020: `∫ b² V(V−2rV′)/(c₀³ r D) dr` | `GM = r(V²−2rVV′)/4` for each observable |
| `R_SPEED_ONLY_`: V=ξ=0, δ live | G100: `∫ 2bδ′/D dr` | G100: `∫ 2b²δ′/(c₀D) dr` | `GM = c₀²r²δ′/2` for each observable |
| `R_TILT_ONLY_`: δ=V=0, ξ live | G001: `∫ bξ′²/(rD) dr` | G001: `∫ b²ξ′²/(c₀rD) dr` | `GM = c₀²rξ′²/4` for each observable |

The displayed integrands are WL `OnDomain/0` fields, and the radial identities are PY condition
field `4`; M3/M4 retain both native constructions and M11 every domain. Each G010 deflection/radar
restriction prints zero. The other monomials remain individually printed; an inactive integral of
zero is left as printed, not privately evaluated. WL conditions instead give the `GM` rule at
`rAnchor`, together with a profile equation over `r`: respectively
`V²=2r(rV′²+V(V′+rV″))`, `2δ′+rδ″=0`, and `ξ′(ξ′+2rξ″)=0`.
M13 verifies literal identity of the WL `GM` and `ProfileEquation` fields just quoted, not of
whole conditions. The PY/WL condition representations have no formed value comparison; no
cross-engine equivalence proof is supplied here. (M3/M4 `R_*_B_CONDITION_*`, M6/M13.)

**Within-engine condition retrieval, all cases.** M13 reads all 14 deflection/radar condition
pairs in each accepted transcript: general B; the three single optical restrictions; constant,
fixed-ratio and power C; the two bulk-only C restrictions; and the five conditional forward
B/C pairs. All 14 PY paired payloads are byte-identical, including their attached domains and
both density-label wrappers where present; no field differs. None of the 14 WL paired payloads
is byte-identical. At the condition's `Cases/0` location, every WL pair has identical `GM`,
`ProfileEquation`, `ProfileEquationRadius` and `Domain`; outer `Domain` and `BranchDomain` also
match literally. The complete differing-field inventory is:

| WL pair family | Differing fields relative to `Cases/0` |
|---|---|
| General B, R_FLOW_ONLY B, combined C_CONSTANT/C_FIXED_RATIO/C_POWER | `InverseWeight`, `ReferenceInverse`, `VelocityReduction/SquaredVelocityEquation` |
| R_SPEED_ONLY B, R_TILT_ONLY B, R_BULK_ONLY C_FIXED_RATIO/C_POWER | `InverseWeight`, `ReferenceInverse`, `VelocityReduction/Cases/0/Domain` |
| Conditional F_FLOW B, F_LIVE B, F_LIVE C_CONSTANT/C_FIXED_RATIO/C_POWER | `InverseWeight`, `ReferenceInverse`, `VelocityReduction/Cases/0/Domain` |

Fixed-ratio/power locations are under each `DensityParts/0/Object` and `DensityParts/1/Object`;
forward locations have the additional outer `Object`. Their density labels/label conditions,
`SuppliedBulkDensityPremise`, and forward `Premise` match literally within each pair. M13 lists
all other identical subtrees and retrieves every differing field verbatim. For example, general
WL `ReferenceInverse` is `(2*GM)/(c0^2*r)` for deflection and `(2*GM)/(c0^3*r)` for radar.
No rescaling or equivalence calculation is made. Whether any differing fields state the same
condition is OPEN, no owner named. Consequently PY supplies no second independent profile test
from these identical pairs; WL's separate-test count remains OPEN. Matching `GM`/`ProfileEquation`
fields alone do not settle the full conditions, and differing serialization alone does not
establish independent physical tests. These limits apply to every condition comparison in
§§3–5, Appendix A and the register. (M13.)

**Implied exchange and domain.** Speed-only and tilt-only both print `jn=0`, with `ρbr` still
symbolic. Flow-only prints both velocity signs on each positive-squared-velocity component, the
zero intervals, differentiable zero-set joins, and excludes negative squared velocity. WL prints

```
Y = -r Inactive[Integrate][(-4 GM)/tailRadius^3, {tailRadius,r,Infinity}],
V = -√Y or +√Y,
```

with its general `r velocityIntegrationConstant` term excluded by supplied asymptotic counting,
and the retained constant domain `velocityIntegrationConstant==0`. Its two positive-component
exchange numerators are the negative and positive versions, respectively, of the literal field

```
(rhoBr[r]*(4*GM + 5*r^2*Inactive[Integrate][(-4*GM)/tailRadius^3,
 {tailRadius,r,Infinity}]) + 2*r^3*Inactive[Integrate][(-4*GM)/tailRadius^3,
 {tailRadius,r,Infinity}]*Derivative[1][rhoBr][r]) /
 (2*r^2*Sqrt[-r*Inactive[Integrate][(-4*GM)/tailRadius^3,
 {tailRadius,r,Infinity}]])
```

No integral is evaluated here. PY separately prints its signed-component and zero-locus exchange
objects. The common coordinate mass expression is printed engine-locally as
`jn=−ρbr V′−Vρbr′−2ρbr V/r`; each restriction retains its own computed substitutions, not a
constant-density shortcut. WL prints component-boundary limits and zero-interval `jn=0`; PY prints
one-sided differentiable-join conditions and the corresponding exchange. Constants with different
engine names are not identified. (M3/M4 `LOCAL_MASS_BALANCE`, `R_FLOW_ONLY_B_IMPLIED_JN_*`;
full WL restrictions M11; cross-engine bounds M6/M8.)

**Bulk density alone, before combined Part C.** `R_BULK_ONLY_` sets V=ξ=0, with `ρbr` symbolic.
On each separately printed `ρ>0` and `ρ≤0` part, PY's condition is
`GM=c₀²(n−1)r²f′/4` for fixed ratio and `GM=c₀²s r²f′/2` for power response, for each observable.
WL prints those anchor rules with profile equations `(n−1)(2f′+rf″)=0` and
`2sf′+rsf″=0`, respectively. Both restrictions print `jn=0`. The exponents `n,s` stay live,
including their degenerate values; no division by `n−1` or `s` and no chosen f-profile is made here.
The exact response domains below travel with both labelled parts. These rows have label/domain
comparisons only, no formed condition or implied-exchange value comparison. (M3/M4
`R_BULK_ONLY_C_{FIXED_RATIO,POWER}_{CONDITION,IMPLIED_JN}_*`; M5–M7/M11/M13.)

**Cross-engine status of the three optical restrictions.** All their deflection/radar/residual rows
have formed value samples; their full along-ray domains do not. Their conditions and implied-jn
rows have no formed comparison; their three γ rows compare only `NOT_ESTABLISHED` markers.
Each active contribution also has printed numeric residuals other than the literal zero string,
retained OPEN in M5/§6. Thus “formed” does not mean equality, and no restriction is reported as a
verified GR-matching substrate. (M5/M6, full coverage M7.)

## 4. Combined Part B and Part C conditions

**Part B.** The native PY `B_CONDITION_DEFLECTION` and `B_CONDITION_RADAR` both print the radial identity

```
GM = r [2c₀²rδ′ + c₀²ξ′² − 2rVV′ + V²]/4.
```

This is a condition on the combination of live speed gradient, embedding slope and radial velocity
and its gradient, relative to independent orbital `GM`, at every far-zone radius on the printed
domain. It does not determine all three profiles. WL prints its `GM` anchor rule and profile
equation quoted in Appendix A; M13 verifies just those paired WL fields as identical while
`InverseWeight`, `ReferenceInverse` and `VelocityReduction/SquaredVelocityEquation` differ.
The whole-condition equivalence is OPEN; the identity is not asserted to have been cross-engine compared. No profile is selected and no matching outcome is
claimed. (M3/M4 `B_CONDITION_*`; M6/M7.)

The printed effective-γ objects preserve the `GM=0` stratum. PY's difference solution set is `{0}`
when `GM≠0`; when `GM=0` and both printed observable functionals vanish it is the full real interval;
the other three combinations of vanishing/nonvanishing functionals print `EmptySet`. WL prints
`(gammaDifference==0 && GM!=0) || (GM==0 && radarFunctional==0 && deflectionFunctional==0)`.
No nonzero-mass division is extended across `GM=0`. These two emitted forms are recorded separately;
M6 counts marker-only comparisons for all 18 γ-family rows. The PY payload identity and the WL
field differences retrieved in M13 yield §3's bounded counting conclusion: no second independent
test established by PY, WL count OPEN. The printed γ relation does not settle equivalence of
the differing WL condition fields. Nothing extends this to mixed grades or full radar time.
(M3 `B_GAMMA_*`; M4/M12 PY solution sets; M6/M8/M13.)

The general implied exchange leaves δ, ξ and `ρbr` live. WL's squared-velocity field is

```
2*c0^2*deltaProfile[r] - r*Inactive[Integrate][
 (-4*GM + 2*c0^2*tailRadius*deltaProfile[tailRadius]
  + c0^2*tailRadius*Derivative[1][xiProfile][tailRadius]^2)/tailRadius^3,
 {tailRadius,r,Infinity}]
```

Both signs, reality/zero-component conditions, asymptotic constant exclusion and differentiable
joining remain attached, as in §3; each `B_IMPLIED_JN_*` prints its corresponding mass-law exchange
with the derivative of live `ρbr`. These objects and their complete domains are **uncompared**,
not an inferred drain law. Full signed expressions are M3/M4 and M11, under the named tags.

**Part C response domains.** The supplied bulk EOS is `P=Kρⁿ`, `cs²=nKρ^(n−1)/m`,
`f=ρ/ρ₀−1`, with `ρ₀>0` a user premise (2026-10-11). `ρ` is bulk number density;
it is not `ρbr`. Local `ρ>0` is a label, not a domain assumption. The outputs carry
`RHO_POSITIVE` / `RHO_NONPOSITIVE` in PY and `rho > 0` / `rho <= 0` in WL. Both parts retain the
same constructed object and their own label conditions; no nonpositive-density part is discarded.
The accepted name table settles these label spellings as names for the same parts (comparator r3 G1,
r4; M5/M7/M8). This is the scoped representational classification, not equality of the conditions.

| Supplied response | Printed first-order δ | Printed condition relative to GM (PY, each observable) | Domain retained |
|---|---|---|---|
| Constant cγ=c₀ | `0` | `GM=r[c₀²ξ′²−2rVV′+V²]/4` | positive real c₀, ρ₀>0 and the substituted ray/functional gates |
| Pointwise fixed cγ/cs=c₀/cs0 | `(n−1)f/2` | `GM=r[c₀²(n−1)rf′+c₀²ξ′²−2rVV′+V²]/4` | exact ratio of nonnegative sound-speed roots: local and asymptotic cs² real, local cs²≥0, asymptotic cs²>0, and exact cγ real and >0; all ray gates use that exact response |
| Power cγ/c₀=(ρ/ρ₀)ˢ | `sf` | `GM=r[2c₀²srf′+c₀²ξ′²−2rVV′+V²]/4` | exact `c₀ densityPower[1+f,s]` real and >0, ρ₀>0, and substituted ray/functional gates |

`densityPower` retains the explicitly printed zero-power value 1. The fixed-ratio exact response
is **not** replaced by `(1+f)^((n−1)/2)` on the nonpositive part. The build r5 disposition records
that this germ-equivalent replacement changes membership on `ρ≤0`; the first-order series alone
cannot supply its domain. Parts C and the bulk-only restrictions print both density labels even
when an individual exponent/parameter choice has no admissible point in a part. No expansion about
`f=0` is promoted to a finite-negative-density accuracy claim. (M3/M4 `LOCAL_C_*`, `C_*`, M8 r5,
M11 exact response domains.)

WL's constant and power N-dependence objects print zero `ResponseDerivative`, `KernelDerivatives`
and `ConditionGMderivatives`; its fixed-ratio fields retain f and its derivative. PY's
N-dependence objects instead print a pair of condition-GM derivative entries, with attached
domains: both zero for constant/power, both `c₀²r²f′/4` for fixed ratio. Those PY tags do not
print WL's response/kernel derivative fields. (M3/M4/M16.) This states where the EOS exponent n
enters the printed first-order conditions, not a selected n or a new stiffness/density relation.
Each combined response has its own signed-velocity reduction and implied `jn`, with live `ρbr`
and its derivative; no relation between bulk and brane densities is supplied. First order in f
supplies no f–ε relation. All these complete conditions, N-dependence and exchange objects are
uncompared values; fixed-ratio/power rows may compare labels and label domains only. (M3–M7/M11.)

## 5. Conditional forward case

**Every object in this section and every `F_` row below is conditional on the adopted no-far-zone-loss
premise:** material leaves the brane only through throats, `jn≡0` in the far zone, the same brane mass
crosses each sphere. `Φ=∮ρbr V·nout dA` uses the coordinate sphere and outward normal, with no supplied
value or sign and no supplied relation to `GM` (user 2026-10-06/07/08; M1).

Under that premise, both engines' `F_MASS_SOLUTION` prints the flux velocity `Φ/(4πr²ρbr)`.
PY's tuple prints the attached `jn=0` domain, `Φ=4π C1`, and the velocity equality; it prints
no mass-residual field. WL prints `Premise/{NormalExchange,Flux,Measure}` and
`Object/{GeneralSolution,FluxVelocity,Residual}`: the general `C[1]/(r²rhoBr[r])` solution,
the flux velocity, and `Residual -> 0`. The zero residual is attributed to WL only. (M16.) This is a mass-balance solution only. The full forward mass object has no
formed comparator value comparison. (M3/M4 `F_MASS_SOLUTION`, M6.)

Conditional `F_FLOW_` substitutes this velocity with δ=ξ=0 into Part A, Part B and branch/traversal
objects. WL's conditional `GM` anchor rule is
`Φ²(5ρbr+2rρbr′)/(64π²r³ρbr³)` at `rAnchor`, with its printed profile equation
`Φ[15ρbr+6r²ρbr′²/ρbr−2r(−7ρbr′+rρbr″)]=0`. PY retains its substituted radial identity.
WL's conditional `F_LIVE_` keeps δ and ξ live and prints the anchor rule

```
GM -> (5*Phi^2)/(64*Pi^2*rAnchor^3*rhoBr[rAnchor]^2)
 + (c0^2*rAnchor^2*Derivative[1][deltaProfile][rAnchor])/2
 + (Phi^2*Derivative[1][rhoBr][rAnchor])/(32*Pi^2*rAnchor^2*rhoBr[rAnchor]^3)
 + (c0^2*rAnchor*Derivative[1][xiProfile][rAnchor]^2)/4
```

with its profile equation and Part C rewrites, all conditional on the same premise. M13 verifies
byte-identical PY forward condition pairs. WL's forward `GM`/`ProfileEquation` fields match,
but `InverseWeight`, `ReferenceInverse` and `VelocityReduction/Cases/0/Domain` differ in each
pair (with the wrappers listed in §3); equivalence and the WL separate-test count remain OPEN.
Their cross-engine condition comparison is unformed. These are computed relations among Φ, free profiles and GM,
not supplied flux normalization and not a matching outcome. No implied `jn` is printed for this
case: `jn=0` is already its premise. Full conditional gates and objects: M3/M4/M11, `F_FLOW_*`,
`F_LIVE_*`; comparison status M5–M7. All conditional forward observable samples inherit §6's
numeric differences, domain gap and vanishing-brane-density coverage limitation. Their shared
nonreciprocal rows also carry §2's r7 K11 limit: K11 supports only the local one-form/closedness
construction and reaches none of the shared observables; the conditional forward marks from
K1–K10 are retrieved separately in M14.

## 6. Cross-engine evidence and limits

The production configuration prints this numerical scope (M5):

> exact substitutions; absolute convergence screening and tanh-sinh error, precision and panel checks;
> complete-functional derivatives checked at two rational step sizes and working precisions;
> 12 printed digits only after reliability checks; numerical evidence not a certified bound;
> unresolved strata are coverage findings

The configuration prints all profile formulas, exact points, name/function bindings and input hashes.
The accepted report says **753 native occurrences (379 PY, 374 WL), 330 joined, 0 unjoined,
49 PY/44 WL local**. Local tags are not joined or counted as value comparisons. The complete local
list and each reason are M6/M7, including engine-local optical/Abel/mass/response construction,
endpoint certificates, import/export metadata and progress names. No export delta or ablation
transcript enters the production comparison. (M5 completion, M6/M7, M8.)

| Printed row class | What the production comparison reaches | What it does not establish |
|---|---|---|
| A_DEFLECTION, A_ROUND_TRIP, A_ONE_WAY_ER/RE, A_NONRECIPROCAL, A_RADAR_LOG; conditional F_FLOW/F_LIVE copies | Formed value samples and markers; all twelve applicable grades retained, four radar grades | No full along-ray domain comparison; no all-profile equality; unresolved brane-zero samples and printed numeric differences remain. Shared nonreciprocal objects retain §2/M14: K11 supports local one-form/closedness only and reaches no shared observable. |
| A_NONRECIPROCAL_PATH_DEPENDENCE | 33 value Boolean samples print `boolean_difference:false`; native on-domain predicate False | No along-ray domain equality or nonradial path statement. K11 does not reach this shared object; only local one-form/closedness claims are supported by K11 (§2/M14). |
| B_RESIDUAL_DEFLECTION/RADAR and R_*/conditional F_* copies | Formed value samples, with printed differences retained | No selected profiles or GR match, and no domain equality |
| B_GAMMA_DEFLECTION/RADAR/DIFFERENCE and R_*/conditional F_* copies (18 rows) | Markers only | γ values/solution sets and their differences uncompared: container/scalar mismatch |
| B_CONDITION_*, B_IMPLIED_JN_*; R_FLOW/SPEED/TILT copies; conditional F_* conditions | No formed value comparison | Full relation/domain/velocity/constant/exchange equivalence OPEN |
| C_CONSTANT conditions/exchange/N-dependence; conditional F_LIVE_C_CONSTANT conditions | No formed value comparison | Complete response-rewritten objects uncompared |
| C_FIXED_RATIO/C_POWER conditions/exchange/N-dependence, R_BULK_ONLY, conditional F_LIVE copies | Label and label-domain comparisons on 22 rows | The condition, N-dependence or exchange value is not compared by a label match |
| BRANCH_EXISTENCE, PATH_TRAVERSAL, BRANCH_TYPE and conditional F_* copies | Stored objects and attempted predicates | No formed along-ray comparison; no equivalence of type/turning-ray constructions |
| Conditional F_MASS_SOLUTION | Raw operands/accounting | No formed comparison of the full mass solution |
| Engine-local rows | Native inventory and parse census only | No cross-engine comparison |

The exact per-row accounting (M6) controls this grouped table. Overlapping formed-row counts are
**value 253, domain 22, marker 271, label 22**; **77 rows have no formed value**, **37 no comparison
of any kind**, and **all 330 retain uncompared leaves**. Parsed / formed across all kinds /
attempted-unformed / unattempted leaf counts are PY **7249253 / 97316 / 2766776 / 4385161**, WL
**1486620 / 119863 / 803332 / 563425**; formed value leaves alone are PY **96781**, WL **119328**.
These are coverage counts, not equality counts. (M6/M7.)

The printed sample census is **5786 value residuals, 33 value Boolean comparisons, 319 unresolved
value samples; 1452 domain Boolean comparisons, 111280 unresolved domain samples; 118870 total
samples**. The domain Booleans are label-domain evidence, not an along-ray test. The literal
`computed_30` residual string `(0.0 + 0.0j)` occurs **4425** times; **1361** other strings occur over
**85 named rows**. These counts use string matches only, not rounding, magnitude, tolerance, or a
mathematical zero test. M5 retains every differing `computed_18` or `computed_30` sample with both
operands, exact assignments, profile formulas, printed residual/reliability fields and location.
The record classifies **all those printed numeric differences as OPEN cross-engine differences;
no adjudication owner named**. Their cause is not assigned to roundoff or to physics by this record.
The 33 value and 1452 domain Boolean comparisons all print `boolean_difference:false` (M5).
Separately, symbolic value leaves print `zero_class:"zero"` on **81 exact-scalar** leaves and
**50 sampled-functional** leaves, and `zero_class:"unresolved"` on **121 sampled-functional**
leaves. These are the comparator's printed classifications, not a new symbolic test. M5 retains
each symbolic operand pair and residual; every unresolved symbolic subtraction is an OPEN
cross-engine comparison, no adjudication owner named. Domain/predicate residual expressions are
attempted truth comparisons with the same unresolved membership/free-name qualifications; no
unresolved Boolean expression is promoted to True or False. Appendix C names every row with a
numeric residual-string exception at either printed precision.

**Every remaining coverage reason and every affected object is retained in M7**, verbatim from the
accepted coverage appendix, including the full no-value/local lists. The two remaining `not_formed`
reasons are **“container/scalar representation differs” (18)** and
**“nested sibling has no bound counterpart” (8722)**. The first is a diagnosed representation/coverage
mismatch for the γ containers (r3/r4 dispositions), not a proof of mathematical equivalence. The
second leaves the missing native field/condition/value content OPEN, no adjudication owner named.
No field is reconstructed from the other engine. M5 lists every false serialization-match location;
its source labels these **“representation only; not a mathematical residual”**. All unreconciled
content behind those diagnostic differences remains OPEN except the sourced density-label spelling
map. The two velocity integration constants remain engine-qualified (M7/M8).

Appendix B reproduces the remaining reason/count inventory compactly; M7 gives **every affected native
spec-object name**, with no omitted reason. Reason/sample pairs overlap: **2 not-formed reasons /
8740 occurrences; 31 sample-unresolved reasons / 197319 occurrences; 8 stratum-unresolved reasons /
56637 occurrences** are not counts of distinct samples or values. No absent sample is called agreement.

**Acceptance limits carried forward (M8).**

- No along-ray domain predicate forms because WL carries free turning radii; a one-sided domain
  corruption could move printed condition text without moving a formed comparison (r2/r4).
- The 18 γ rows compare markers only. Label matches do not compare Part C values. Every no-value row
  remains listed; no `comparison_formed` flag is promoted to a value result (r2–r4).
- Vanishing-brane-density samples retain interior-denominator-zero/integrability or pole-set reasons.
  The r4 acceptance counts 183 and 27 occurrences of its two pole-reason families and 105 unresolved
  value samples in that stratum. The broader printed value-unresolved census remains 319; these
  differently scoped counts are not substituted for one another.
- C2's repaired formation property is **amplitude-independence only**. Integrals cancelling near zero
  over an O(1) integrand can remain unresolved; r4's probe has that issue, no full-02 sample does.
  Precision/step/panel consistency is numerical evidence, not a certified bound.
- The tests do not cover every WL precedence/chained-comparison behavior. Production real-stream
  recomputations support this filed run; reuse on other transcripts inherits no such check. The label
  table orientation rests on its cited engine lines, not on the relative swap control alone (r4).
- Domain-stratum inventory order can differ across runs; the accepted full-02 itself is filed.
  No rerun or byte-identity promise is made (r2 observation O2, r4 c).
- Earlier comparator r0/r1 failures and coverage restrictions concern those versions; repairs r2–r4
  changed their scope. In particular the 36 mixed-grade rows now form value samples on generic
  inverse/shifted/exponential profiles, with the remaining brane-zero gap, rather than remaining
  wholly untested. The old `approx_30` caveat is not a claim of 30 reliable production digits.
- Build r5 retains exact fixed-ratio response-domain dependence on both density parts. Build r11
  accepts each final K10 gate witness in its own engine; it does not compare full engines. The
  `WholeEmittedObjectsEvaluated` flag means absence of unevaluated heads, not correct values. Superseded
  runs mentioned by the disposition source no result here. The r7 K11 claim limit remains unchanged
  through r11: PY evaluates 3 tags/376 `NOT_EVALUATED`; WL reaches 10/364 `NOT_EVALUATED`;
  neither reaches shared `A_NONRECIPROCAL_*`. K11 supports the local one-form and closedness
  construction only, not the shared observables. §2/M14 give the other knives' printed reach
  and structural marks. Unreached parts remain “not evaluated” and “not compared”; identity
  marks or witness absence cannot supply cross-engine equality.
  This record cites no software-check Boolean as physics evidence.

## 7. Conditional premises and remaining owners

Every adopted premise below is an **adopted substrate input to a conditional model**, conditional
and re-tested at the knit (user 2026-10-09). None is derived. The dependence column describes
Parts A–C actually reported here; held Part D supplies no computed object. (M1 D2; M9 O2.)

| Premise, source/date | Content and dependence |
|---|---|
| No far-zone loss, record directive/user 2026-10-06/07/08 | `jn≡0`, conserved coordinate-sphere flux; only `F_MASS_SOLUTION`, `F_FLOW_*`, `F_LIVE_*` depend on it. Ordinary R_FLOW_ONLY does not. |
| w-parity, D2, 2026-10-07 | Neutral far field symmetric under w→−w, ξw and material w-velocity vanish there. Scoped to held Part D; no A–C object uses it. Tilt stays live in the general A–C objects. |
| P5, D2, 2026-10-07; carriage 2026-10-08 | Momentum density ρbrV and carried current ρbrVⁱVʲ. Part D only; no A–C object. Additional stress-carried momentum/retained grade remains OPEN, S8. |
| P6, D2, revised 2026-10-08 | Isotropic pressure plus linear viscous Cauchy stress; live pbr(ρbr), ccomp²=dpbr/dρbr, η, ζ. Part D only; no A–C object. Not derived from P1; power-law creep excluded from that adopted form, relaxation and supplier/power identification OPEN. |
| Density link, D2, 2026-10-06 | `cγ²=μ⊥/ρbr`, with `μ⊥∝ρbr^α`, α live. The extra constitutive link is Part D only; no A–C object uses it. A–C keep speed change and brane density independent live profiles, subject to the supplied optical ratio; **no effect of the link is stated**. |
| O2 P1, 2026-10-06 | Elastic optical shear, relaxation under steady load; explicit relaxation power and any supplier/budget. A–C use only the separately supplied optical dispersion/ratio, not a relaxation law. No computed A–C relaxation/power object; compatibility is later. |
| O2 P2, 2026-10-06 | Dynamical order-conversion drive, no separate external body-force entry. No A–C momentum/drive object uses it; no drive(GM) supplied. |
| O2 P3, 2026-10-06 | Converted material carries local V, closing in-plane carried jnV only. A–C's coordinate mass input is distinct; no A–C momentum/exchanged-energy object uses P3. |
| O2 P4, 2026-10-06 | Postulated shear-free scalar bulk, native face-normal load, no independent tangential bulk stress. No A–C load object uses it. Part C's supplied bulk EOS/acoustic response is not a calculation of this loading. |

`ρ₀>0` is the separate supplied density premise, not an inferred sign. The supplied dispersion,
advection, metric, optical ratio, common polarization speed, reference speed, mass law, LAB_HELD,
counting and bulk responses remain conditional inputs even where they are not one of these adopted
brane premises. The uniform S11b anchor remains at V=ξ=0 with uniform coefficients; it supplies no
live flowing material law. (M1/M3/M4/M11.)

| Remaining dependence | Relation and source-named owner |
|---|---|
| S8 | Substrate material identity and stiffness; P5's possible additional stress-carried momentum remains its question. No physical constitutive realization is delivered. |
| S11c | Nonuniform brane operator unfinished; the supplied local eikonal equation is not derived from it. No subprincipal/amplitude/polarization-transport completion. |
| O1 ℳ⊥ | Stiffness response remains live/unselected in A–C; S8/substrate reduction/S22 completion (O2 §7). |
| O2 | Conditional live steady accounting is the held Part D handoff, not an optical closure; its returned interpretation remains below. |
| O3 Πn | Momentum/exchange partners not computed; S12 reactions and O6/normal reduction. A–C computes only coordinate mass exchange. |
| O4 ℰh | ξw=ℓh is supplied; live embedding/longitudinal relation and identity with O2 normal content remain OPEN. Q1/Q2 static ingredients, O4 live identification, S8/S22 completion. |
| O5 ℋcore | Holder/mouth data not selected by LAB_HELD; Q2/S22. |
| O6 𝒥map | No native-face/finite-thickness or induced-measure law is computed. S12 dynamical conversion/native data; S14a projected source bridge distinct. |
| O7 ℛbr | `ρbr` retained live, not eliminated with the bulk EOS. S8/substrate reduction/S22 completion; unsupplied derivative/grade contracts have no owner named. |
| GM, jn, f, counting | GM owned by gravity; jn, f and optical counting by gravity or S12. Relation of jn to S11b's normal drain vdr OPEN; Φ–GM relation OPEN, no separate adjudication owner named. |
| Source/response bridges | O2 §7 retains S13 GAP, S14 conditional on S14a, S15 source-side conditional on S14a, S14a projected-source/flux/return bridge and S16 response-side GM matching (S16 assignment an inference). None is closed here. |
| Knit | S21 owns integration/input-versus-consequence sorting; material that meets all sectors is chosen there, after the sectors, not in S9b. |

The following remain excluded/open, with **no owner named by the A–C source unless stated**:
radial/tangential anisotropic stiffness, coupling to thickness/bulk fields,
polarization-dependent propagation, mixed ωk other than advection, explicit subprincipal terms,
strong field, mouth/interior, moving/rotating mass and time-changing drain. The retained historical
zero-width sheet, continuum, amplitude→0, nondissipative/frequency-independent optical limits are not
lifted. In-plane V and the isotropic strain content in cγ are live; anisotropic strain remains outside.
S11b's normal vdr freeze and uncarried `O(v₀|qn|/ω)` correction remain distinct from both V and jn;
S12 owns native drain/return data, with the historical correction's C/nonlinear-program handoff
retained in `R-S12-02`. Optical compatibility with steady relaxation remains later, owner unassigned.
No induced-measure correction or order is claimed. (M1; M9 §7; M10.)

**O2 obligation, undischarged and owned by the orchestrator.** D5 carries the WL-only `ξw″` key
with first-derivative keys `V_r, delta, f, h, j_n, mu_perp, o2_rho_br_live`, all at the stored radial
argument `sqrt(x1²+x2²+x3²)`. O2 §4 locates all eight at each momentum-flux role in
`hold_inplane[0..2]` (three roles each), `hold_bulk` (three), `hold_normal` (all twelve), and
`energy_balance` (three `OPEN_MaterialEnergyFlux_i` roles). The interpretation is an
“undischarged sub-step-7 obligation returned to the orchestrator under M1”. S9b without Part D
cannot discharge it; neither a finite dependence inventory nor an adopted Part-D premise erases
it. O2's other OPEN action/orientation/nesting differences and force/power compatibility remain its
bounded handoff, not S9b results. An unresolved difference alone establishes no substrate requirement;
no new entry is created for this interpretation. (M1 D5, M9 §§4/8.)

## 8. Register pass, source tensions and open questions

The S9b pass in [`SUBSTRATE_REQUIREMENTS.md`](../SUBSTRATE_REQUIREMENTS.md) applies its rest-on test
and both routes. **16 (`caa0710d`, before S9b) → 17 entries**, all OPEN; round 1 itself keeps 17 → 17. One new entry, **R-S21-01**, target **S21**, records
the conditional optical requirement, both observable conditions and Part C rewrites, with symbolic
ρbr, the printed implied-jn branches and all free quantities retained. It names no means of meeting
it. Existing **R-S8-01** and **R-S8-06** gain S9b sources for the supplied transverse stiffness input
and material identity respectively, without widening their original inertia/stiffness domains into
a completed live constitutive law. M10 retrieves the explicitly pinned before text/headings and the full after-pass register.

The requirement is admitted because the reported conditional reference-matching statement rests
on these printed profile conditions; failure removes that match, while the conditional functionals
remain. It is not an assertion that a substrate meets them. Other unselected closure forms,
excluded effects and unresolved engine differences remain handoffs rather than additional entries.
Adopted premises retain their dates/status for knit re-testing; no postulate is renamed a derivation.
No ansatz or defect register is edited.

For the second route, the repository identifies more than GR/PPN references. The inventory below
names the prior-art objects relevant to the optical, flow, embedding and carried mechanical
questions. `survey` means repository-root `docs/medium_requirements_and_prior_art.md` at
`36f59ff1`; `scoping` means `directives/O2_steady_brane_balance_scoping.md` §5 at that revision.
M15 retrieves both sources and `directives/prior_art_4d_bulk_search_prompt.md`; M10 retrieves
the pre-S9b register's O2 second route. The access qualifications below are the sources' reports,
not fresh literature verification. **No oracle is applied in this record.**

| Identified object and repository source | Use in this pass |
|---|---|
| GR/PPN deflection and logarithmic radar references, A–C spec (M1) | Supplied comparison inputs in the accepted constructions; no new oracle check or GR solution reproduced here. |
| Cahill (2003) inflow; Hamilton–Lisle (2008), *The river model of black holes*; Sbitnev inflow; Unruh acoustic sinks — survey line 246, cautions 263–266; 4D search prompt | All identified, none applied. Survey scopes them to 3D with no off-sheet bulk and flags Cahill's status. No inflow profile or GM relation is imported. |
| Hamilton–Lisle (2008), AJP 76, 519–532; Visser (1998), *Acoustic black holes*; Unruh (1981), *Experimental Black-Hole Evaporation?* — scoping §5 | None applied. Scoping reports manuscripts read for Hamilton–Lisle/Visser, only abstract/metadata for Unruh. GR river and barotropic/inviscid/irrotational scalar acoustics require matching time, velocity, dimension, EOS/exchange and independent mass calibration; they supply no shear-sheet/off-sheet law. |
| Gurtin–Murdoch (1975), *A continuum theory of elastic material surfaces*; Scriven (1960), *Dynamics of a fluid interface*; Helfrich (1973), *Elastic Properties of Lipid Bilayers* — scoping §5; O2 register second route | None applied. Surface stress/traction, Newtonian surface-fluid and curvature-sector candidates; Gurtin–Murdoch abstract/bibliography only, Scriven/Helfrich full texts unverified. Ordinary 2D surfaces in 3D bulk do not establish this 3D-brane/4D-bulk law. |
| Abels–Garcke–Grün (2012), *Thermodynamically Consistent, Frame Indifferent Diffuse Interface Models…*; Rodriguez–Hoger–McCulloch (1994), *Stress-dependent finite growth in soft elastic tissues* — scoping §5; O2 register second route | None applied. Diffusive exchange-momentum and evolving-reference/growth candidates; manuscript read for the former, metadata only/full text unverified for the latter. No two-fluid, mobility, growth or reference law is adopted. |
| MacCullagh (1839) curl-only light / Whittaker Ch. V; Cauchy–Navier elasticity — survey attribution table; register R-S8-04/05 (M10) | None applied as an oracle. Existing light-carrier, longitudinal/compression, angular-momentum and objectivity obligations remain OPEN; the supplied optical ratio does not re-certify them. |
| Structural acoustics (fluid-loaded plate, radiation impedance, added mass); Fruchart–Scheibner–Vitelli (2023) passivity/Onsager–Casimir; Cembranos–Dobado–Maroto (2003) branon fluctuation — survey attribution table | None applied. Adjacent interface, driven-response and embedding context; no interface law, drive or curvature dynamics is computed in A–C. |
| Analog/emergent gravity (Unruh, BLV, BEC); Volovik/³He; Sakharov/Volovik induced gravity — survey prior-art table/unified findings | None applied. The survey's kinematics-versus-dynamics and undelivered-G qualifications remain; optical conditional requirements are not a derivation of gravity or G. |
| Rubakov–Shaposhnikov (1983), Akama (1982), RS2, DGP (2000); material-medium superfluid-vacuum class (Volovik, Sbitnev, Consoli, Huang, Fedi) — survey 4D-bulk table | None applied. Adjacent embedding/material-bulk context only; no localization or material model is supplied to A–C by these names. |

This route supplies no new entry: no additional sourced requirement follows from an oracle result
applied here. It does not restrict the register's route to reproduced results or erase identified
candidates. Existing MacCullagh obligations remain OPEN; the conditional optical functionals and
M13's bounded condition comparison do not discharge them. Applicability checks of these candidate
oracles remain unperformed, not negative results. The register's earlier pass-2 disagreement stays
unresolved (M10).

**Source tensions, with both sides.**

- O2 §8 says an induced-measure reading has a **“relative O(ε) correction to j_n”**. D4 says
  **“No source supplies that order”**; the accepted A–C spec says **“No induced-measure object,
  re-expression or correction is computed or printed.”** This record uses only coordinate d³x;
  the O2 qualification is preserved as a source tension, not inherited as an S9b result. (M1/M9.)
- The earlier decision list D3 says **“Part D re-expresses the Part B conditions with that balance
  in force.”** The later record directive says **“Part D is not built in S9b”** and its spec source's
  Part D is **“held and not accepted”**. The later user scope controls; D2 premises are recorded but
  Part D sources no computed claim. (M1.)
- O2 acceptance assigned **“Interpreting it belongs to the O2 record (sub-step 7)”**; O2 returned it
  as **“undischarged sub-step-7 obligation returned to the orchestrator under M1”**; D5 now explicitly
  says **“S9b does not discharge it.”** The obligation remains owned and undischarged. (M1/M9.)
- PY's printed density-function assumptions include **`'positive': True`**; WL's branch predicate
  states **`rhoBr[r] != 0`**, real coefficients and positive `muPerp[r]/rhoBr[r]`. No comparator
  domain equality settles this difference. The density in these quotations is brane density,
  not either bulk label. (M2–M4/M11.)
- The register's earlier pass-2 directive limited the second route to a result **“that a record in
  the reading list identifies as reproduced by the sector”**; the register says its MacCullagh
  entries were **“Neither is anywhere in S9's or S10's records”** and **“Run this route for every
  sector with identified prior art”**. S9b follows the register's method; the prior disagreement
  stays unresolved. There is no new operative S9b directive/register disagreement. (M10.)

Questions requiring work beyond retrieval remain OPEN: whether differing within-WL condition fields
state the same condition and whether the observables furnish separate WL tests; complete cross-engine equivalence or
physical interpretation of numeric/domain/condition/γ/exchange differences; the free turning-radius,
constant and uncovered stratum questions; a substrate derivation and physical holder for the local
wave equation; optical/relaxation compatibility; whether any joint live profiles and a physical
source satisfy the conditions relative to orbital GM; Φ–GM and jn–vdr relations; O1–O7 closure and
O2's returned derivative interpretation. Ownership remains as in §7; no computation to answer these
questions is part of this record.

## Appendix A. Native WL combined-condition fields

Retrieved by M3/M11/M13, `B_CONDITION_DEFLECTION` and `B_CONDITION_RADAR`, `Cases/0`; domains
are §2/M11. M13 verifies literal identity of the two fields quoted here, `GM` and `ProfileEquation`.
`InverseWeight`, `ReferenceInverse` and `VelocityReduction/SquaredVelocityEquation` differ;
whole-condition equivalence remains OPEN (§3). Cross-engine comparison remains unformed.

```
GM -> -1/4*(rAnchor*(-velocityProfile[rAnchor]^2
 - 2*c0^2*rAnchor*Derivative[1][deltaProfile][rAnchor]
 + 2*rAnchor*velocityProfile[rAnchor]*Derivative[1][velocityProfile][rAnchor]
 - c0^2*Derivative[1][xiProfile][rAnchor]^2))

velocityProfile[r]^2 + c0^2*(4*r*Derivative[1][deltaProfile][r]
 + Derivative[1][xiProfile][r]^2 + 2*r^2*Derivative[2][deltaProfile][r]
 + 2*r*Derivative[1][xiProfile][r]*Derivative[2][xiProfile][r])
 == 2*r*(r*Derivative[1][velocityProfile][r]^2
 + velocityProfile[r]*(Derivative[1][velocityProfile][r]
 + r*Derivative[2][velocityProfile][r]))
```

## Appendix B. Remaining coverage reasons

The table below retrieves the accepted coverage appendix's literal kind/reason/count columns.
The full affected-object column is M7; exact native field/sample paths remain in the filed comparator.
All entries are unresolved coverage, not evidence of equality; no adjudication owner is named.

| Kind | Reason | Count |
|---|---|---:|
| not_formed | container/scalar representation differs | 18 |
| not_formed | nested sibling has no bound counterpart | 8722 |
| sample_unresolved | domain translation: no bound counterpart | 92290 |
| sample_unresolved | domain translation: unapplied live function rho_br_profile | 2178 |
| sample_unresolved | domain translation: unapplied live function s9b_V | 2178 |
| sample_unresolved | domain translation: unapplied live function s9b_delta | 2178 |
| sample_unresolved | domain translation: unapplied live function s9b_f | 2178 |
| sample_unresolved | domain translation: unapplied live function s9b_xi | 2178 |
| sample_unresolved | domain translation: unsupported constructor Text | 19965 |
| sample_unresolved | free names py::s9b_velocity_integration_constant | 462 |
| sample_unresolved | free names py::s9b_velocity_zero_radius; set equality/emptiness could not be decided | 2032 |
| sample_unresolved | free names wl::impactRadar,wl::rTurnRadar; free names wl::rTurnFlyby; free names wl::rTurnRadar; quantified condition could not be decided over its full bound domain | 21708 |
| sample_unresolved | free names wl::rAnchor | 4304 |
| sample_unresolved | free names wl::rTurnFlyby; free names wl::rTurnRadar | 2376 |
| sample_unresolved | free names wl::velocityIntegrationConstant | 8592 |
| sample_unresolved | operand evaluations unavailable for subtraction | 105 |
| sample_unresolved | quadrature UnresolvedEvaluation: held derivative: complete functional unavailable at b=129979/70000: interior denominator zero at 23/7; local integrability not established | 9 |
| sample_unresolved | quadrature UnresolvedEvaluation: held derivative: complete functional unavailable at b=21997/10000: interior denominator zero at 17/5; local integrability not established | 9 |
| sample_unresolved | quadrature UnresolvedEvaluation: held derivative: complete functional unavailable at b=69991/30000: interior denominator zero at 11/3; local integrability not established | 9 |
| sample_unresolved | quadrature UnresolvedEvaluation: interior denominator zero at 11/3; local integrability not established | 26 |
| sample_unresolved | quadrature UnresolvedEvaluation: interior denominator zero at 17/5; local integrability not established | 26 |
| sample_unresolved | quadrature UnresolvedEvaluation: interior denominator zero at 23/7; local integrability not established | 26 |
| sample_unresolved | quadrature UnresolvedEvaluation: nonpolynomial denominator zero set not established for -11/3 + 7*sec(_bound_0)/3 | 9 |
| sample_unresolved | quadrature UnresolvedEvaluation: nonpolynomial denominator zero set not established for -17/5 + 11*sec(_bound_0)/5 | 9 |
| sample_unresolved | quadrature UnresolvedEvaluation: nonpolynomial denominator zero set not established for -23/7 + 13*sec(_bound_0)/7 | 9 |
| sample_unresolved | quantified/set truth requires root isolation with non-rational coefficients; not implemented; set equality/emptiness could not be decided | 6867 |
| sample_unresolved | sampled predicate: Invalid NaN comparison | 250 |
| sample_unresolved | sampled predicate: Invalid comparison of non-real -55*sqrt(221)*I/3094 | 122 |
| sample_unresolved | sampled predicate: Invalid comparison of non-real 19*I/8 | 144 |
| sample_unresolved | sampled predicate: Invalid comparison of non-real I | 100 |
| sample_unresolved | sampled predicate: Invalid comparison of non-real zoo | 247 |
| sample_unresolved | set equality/emptiness could not be decided | 24555 |
| sample_unresolved | set is not a scalar truth value | 2178 |
| stratum_unresolved | no bound counterpart | 2614 |
| stratum_unresolved | no sampled point established membership; piece not proved empty | 53088 |
| stratum_unresolved | unapplied live function rho_br_profile | 66 |
| stratum_unresolved | unapplied live function s9b_V | 66 |
| stratum_unresolved | unapplied live function s9b_delta | 66 |
| stratum_unresolved | unapplied live function s9b_f | 66 |
| stratum_unresolved | unapplied live function s9b_xi | 66 |
| stratum_unresolved | unsupported constructor Text | 605 |

## Appendix C. Rows with printed numeric residual-string differences

M5 `SAMPLE_LITERAL_EXCEPTION_COUNTS` retrieves this inventory. Each count is of stored samples
whose `computed_18` or `computed_30` differs literally from `(0.0 + 0.0j)` (or a stored Boolean
comparison reports a difference). Both precision fields, operands, profile/point, residual and
reliability text are preserved per sample in M5. Every listed numeric difference remains **OPEN,
no adjudication owner named**. No entry is labelled a physical disagreement or rounded away.
The F_FLOW/F_LIVE rows remain conditional on §5's no-far-zone-loss premise.

| Native shared object | Literal-exception samples |
|---|---:|
| `A_DEFLECTION_G001` | 9 |
| `A_DEFLECTION_G020` | 9 |
| `A_DEFLECTION_G021` | 8 |
| `A_DEFLECTION_G100` | 9 |
| `A_DEFLECTION_G101` | 8 |
| `A_DEFLECTION_G120` | 8 |
| `A_DEFLECTION_G121` | 7 |
| `A_ONE_WAY_ER_G001` | 10 |
| `A_ONE_WAY_ER_G020` | 9 |
| `A_ONE_WAY_ER_G021` | 24 |
| `A_ONE_WAY_ER_G100` | 10 |
| `A_ONE_WAY_ER_G101` | 24 |
| `A_ONE_WAY_ER_G120` | 24 |
| `A_ONE_WAY_ER_G121` | 21 |
| `A_ONE_WAY_RE_G001` | 10 |
| `A_ONE_WAY_RE_G020` | 9 |
| `A_ONE_WAY_RE_G021` | 24 |
| `A_ONE_WAY_RE_G100` | 10 |
| `A_ONE_WAY_RE_G101` | 24 |
| `A_ONE_WAY_RE_G120` | 24 |
| `A_ONE_WAY_RE_G121` | 21 |
| `A_RADAR_LOG_G001` | 26 |
| `A_RADAR_LOG_G020` | 19 |
| `A_RADAR_LOG_G100` | 20 |
| `A_ROUND_TRIP_G001` | 9 |
| `A_ROUND_TRIP_G020` | 8 |
| `A_ROUND_TRIP_G021` | 24 |
| `A_ROUND_TRIP_G100` | 3 |
| `A_ROUND_TRIP_G101` | 24 |
| `A_ROUND_TRIP_G120` | 24 |
| `A_ROUND_TRIP_G121` | 21 |
| `B_RESIDUAL_DEFLECTION` | 18 |
| `B_RESIDUAL_RADAR` | 26 |
| `F_FLOW_A_DEFLECTION_G020` | 10 |
| `F_FLOW_A_ONE_WAY_ER_G020` | 23 |
| `F_FLOW_A_ONE_WAY_RE_G020` | 23 |
| `F_FLOW_A_RADAR_LOG_G020` | 33 |
| `F_FLOW_A_ROUND_TRIP_G020` | 26 |
| `F_FLOW_B_RESIDUAL_DEFLECTION` | 10 |
| `F_FLOW_B_RESIDUAL_RADAR` | 13 |
| `F_LIVE_A_DEFLECTION_G001` | 9 |
| `F_LIVE_A_DEFLECTION_G020` | 10 |
| `F_LIVE_A_DEFLECTION_G021` | 8 |
| `F_LIVE_A_DEFLECTION_G100` | 9 |
| `F_LIVE_A_DEFLECTION_G101` | 8 |
| `F_LIVE_A_DEFLECTION_G120` | 8 |
| `F_LIVE_A_DEFLECTION_G121` | 7 |
| `F_LIVE_A_ONE_WAY_ER_G001` | 10 |
| `F_LIVE_A_ONE_WAY_ER_G020` | 23 |
| `F_LIVE_A_ONE_WAY_ER_G021` | 27 |
| `F_LIVE_A_ONE_WAY_ER_G100` | 10 |
| `F_LIVE_A_ONE_WAY_ER_G101` | 24 |
| `F_LIVE_A_ONE_WAY_ER_G120` | 27 |
| `F_LIVE_A_ONE_WAY_ER_G121` | 24 |
| `F_LIVE_A_ONE_WAY_RE_G001` | 10 |
| `F_LIVE_A_ONE_WAY_RE_G020` | 23 |
| `F_LIVE_A_ONE_WAY_RE_G021` | 27 |
| `F_LIVE_A_ONE_WAY_RE_G100` | 10 |
| `F_LIVE_A_ONE_WAY_RE_G101` | 24 |
| `F_LIVE_A_ONE_WAY_RE_G120` | 27 |
| `F_LIVE_A_ONE_WAY_RE_G121` | 24 |
| `F_LIVE_A_RADAR_LOG_G001` | 26 |
| `F_LIVE_A_RADAR_LOG_G020` | 33 |
| `F_LIVE_A_RADAR_LOG_G100` | 20 |
| `F_LIVE_A_ROUND_TRIP_G001` | 9 |
| `F_LIVE_A_ROUND_TRIP_G020` | 26 |
| `F_LIVE_A_ROUND_TRIP_G021` | 27 |
| `F_LIVE_A_ROUND_TRIP_G100` | 3 |
| `F_LIVE_A_ROUND_TRIP_G101` | 24 |
| `F_LIVE_A_ROUND_TRIP_G120` | 27 |
| `F_LIVE_A_ROUND_TRIP_G121` | 24 |
| `F_LIVE_B_RESIDUAL_DEFLECTION` | 13 |
| `F_LIVE_B_RESIDUAL_RADAR` | 29 |
| `R_FLOW_ONLY_A_DEFLECTION_G020` | 9 |
| `R_FLOW_ONLY_A_RADAR_LOG_G020` | 19 |
| `R_FLOW_ONLY_B_RESIDUAL_DEFLECTION` | 16 |
| `R_FLOW_ONLY_B_RESIDUAL_RADAR` | 24 |
| `R_SPEED_ONLY_A_DEFLECTION_G100` | 9 |
| `R_SPEED_ONLY_A_RADAR_LOG_G100` | 20 |
| `R_SPEED_ONLY_B_RESIDUAL_DEFLECTION` | 11 |
| `R_SPEED_ONLY_B_RESIDUAL_RADAR` | 23 |
| `R_TILT_ONLY_A_DEFLECTION_G001` | 9 |
| `R_TILT_ONLY_A_RADAR_LOG_G001` | 26 |
| `R_TILT_ONLY_B_RESIDUAL_DEFLECTION` | 9 |
| `R_TILT_ONLY_B_RESIDUAL_RADAR` | 18 |
