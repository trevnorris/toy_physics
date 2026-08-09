# S10 — the transverse mode count supplied by the brane dimension

## Record boundary

This record reports the S10 computation that was actually emitted by the live
SymPy and Wolfram engines and the narrower set of claims that survives their
live comparator. It does not turn a supplied action into an explanation of that
action, and it does not turn the symbolic sweep in D into a derivation that the
physical brane has D = 3.

The measured generic result is:

> For the supplied curl-only in-plane action, at nonzero wavevector and away
> from allowed exceptional strata, the nonzero root has D − 1 transverse null
> directions in every measured MAIN case D = 2, 3, 4, 5. Thus the D = 3 member
> has two transverse directions.

The form controls prohibit a stronger reading of that sentence. At both
dimensions where the form-changing packages and MAIN were emitted, D = 3 and
D = 4, FULLGRAD has the same nonzero root and the same D − 1 transverse
nullity as MAIN. The transverse count therefore cannot be attributed to the
curl-only form. What that form determines here is the disposition of the one
remaining, longitudinal direction: unlike FULLGRAD, it leaves that direction
at the zero root instead of putting it on the propagating root. DIVONLY, which
reverses the two sectors, confirms that the stiffness form can change the
count even though curl-only is not what uniquely produces D − 1.

A separate condition comes from the inertia structure. `XFORM_ANISO` is not a
stiffness-form control and does not multiply the inertia by one scalar: it
keeps `S_curl` and changes the kinetic quadratic form on one distinguished
axis to `s_ρ(∂_t u_1)² + Σ_{j=2..D}(∂_t u_j)²`, with `s_ρ > 0` and
`s_ρ ≠ 1`. Thus the isotropic inertia value is not admitted by this control.
Writing each root as `N2 / N3`, where N2 is total nullity and N3 is the
basis-independent exactly-transverse nullity, the committed readings are:

| D | case | roots in emitted order: N2 / N3 | total N2 / N3 over propagating roots |
|---:|---|---|---:|
| 3 | MAIN, generic | zero `1 / 0`; positive `2 / 2` | `2 / 2` |
| 3 | ANISO, generic | zero `1 / 0`; positive `1 / 1`; positive `1 / 0` | `2 / 1` |
| 3 | ANISO, allowed stratum | zero `1 / 0`; positive `2 / 2` | `2 / 2` |
| 4 | MAIN, generic | zero `1 / 0`; positive `3 / 3` | `3 / 3` |
| 4 | ANISO, generic | zero `1 / 0`; positive `2 / 2`; positive `1 / 0` | `3 / 2` |
| 4 | ANISO, allowed stratum | zero `1 / 0`; positive `3 / 3` | `3 / 3` |

Therefore the measured inertia-form change does not remove a propagating
mode: at both shared dimensions its total N2 nullity over the positive roots
remains D − 1, exactly as in MAIN. In the same breath, it does generically
reduce the number of those modes that are exactly transverse from D − 1 to
D − 2: ANISO puts N3 nullity D − 2 on one positive root and zero on its
additional positive root, instead of MAIN's N3 nullity D − 1 on one positive
root. MAIN itself is unchanged at every measured D; this is an added condition
on the headline's generality, not a retraction of the MAIN result.

The allowed ANISO stratum is `k_2 = ... = k_D = 0` with `k_1 ≠ 0`, so the
wavevector has zero obliquity to the distinguished inertia axis. There the two
generic positive branches coalesce, and the rerun restores N2 / N3 totals
`2 / 2` at D = 3 and `3 / 3` at D = 4. This is distinct from restoring
isotropic inertia: `s_ρ = 1` remains excluded. The scalar stiffness-coefficient
control XCOEF_SCALE and the stiffness-sign control SIGNFLIP leave N3 unchanged
in their measured cases, but ANISO shows that this control-specific result is
not a general rule that only a stiffness-form change can move the
exactly-transverse count. No control in this sweep varies the in-plane field
content, the rest-background premise, dissipation, or the quadratic
truncation; those remain supplied conditions.

XCOEF_SCALE nevertheless moves the nonzero root's coefficient, and SIGNFLIP
moves that root's sign and stability; their unchanged N3 counts do not mean
that nothing physical moved.

The N3 rerun also shows that modal obliquity vanishes on that stratum: every
propagating null direction there is exactly transverse.

SymPy emits the MAIN operands at :366, :397, :404, :462, :469; :677,
:708, :715, :773, and :780, and the ANISO generic and stratum operands at
:2953, :2984, :2991, :3049, :3056, :3114, :3121, :3269, :3290, :3305,
:3312, :3348, :3355; :3468, :3499, :3506, :3564, :3571, :3629, :3636,
:3775, :3796, :3811, :3818, :3854, and :3861 in
scripts/out/S10_brane_mode_spectrum_sympy_audit.out. Wolfram emits the
corresponding operands at :240, :257, :260, :270, :273; :462, :480, :483,
:493, :496; :2113, :2157, :2160, :2170, :2173, :2183, :2186, :2246,
:2259, :2276, :2279, :2289, :2292; :2453, :2515, :2518, :2528, :2531,
:2541, :2544, :2601, :2614, :2631, :2634, :2644, and :2647 in
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out.

The generic per-root N2 and N3 integer rows above have joined comparator PASS
verdicts; the stratum locus and rerun are paired raw readings only, because
their engine tag names do not join. The positivity of ANISO's third root is a
premise-backed inference from its emitted expression under `ρ_br > 0`, `μ_R > 0`,
`s_ρ > 0`, and nonzero real wavevector, not a joined sign verdict: its live
sign rows remain unresolved by the comparator. N3, as specified at
directives/S10_SHARED_PHYSICS.md:230-257, is the explanatory object for every
transverse count here; none is inferred from display-only N6.

The generic verdict blocks are at
scripts/out/S10_cross_engine_comparator.out:344-369, :417-442, :582-607,
:655-680, :1319-1344, :1386-1411, :1453-1478, :1603-1628, :1670-1695,
and :1737-1762. The corresponding stratum names appear only in the engine-only
sets at :4698-4749, :5136-5187, :8009-8033, and :8254-8278.

Both engines consume the same shared action and assumptions. Their agreement is
therefore protection against independent implementation errors, not independent
physical evidence for the curl-only premise. This limitation belongs to the
headline because the headline is otherwise independently quotable.

The physical selection D = 3 is not made in S10. The live S10 computation keeps
D symbolic for dimensions and evaluates an indexed sweep at D = 2, 3, 4, 5.
Consequently, S10 establishes the conditional map D ↦ D − 1 for the cases
measured; it does not establish which D nature selects.

## What was supplied and what was computed

The shared specification supplies:

- a D-component in-plane displacement u, with its separation from every
  out-of-plane field inherited rather than tested;
- the real cosine plane-wave ansatz;
- positive inertia and stiffness coefficients, nonzero real wavevector,
  unstrained rest background, no dissipation, and linear response;
- the curl-only stiffness density

      S_curl = (1/2) Σ_i Σ_j (∂_i u_j − ∂_j u_i)²

  in

      L = (ρ_br/2) Σ_j (∂_t u_j)² − (μ_R/2) S_curl.

These premises and the action are at
directives/S10_SHARED_PHYSICS.md:13-28, :30-47, and :82-107. The action is an
input, not a result of the mode count.

From that input each engine constructs the Euler–Lagrange equations; obtains the
dynamical matrix by an equation-of-motion route and a quadratic-form route;
uses the quadratic-form matrix downstream; solves its determinant for
omegaSquared; and, at every root, measures:

- N2: rank and nullity of the root matrix;
- N3: rank after stacking the wavevector row, hence the basis-independent
  dimension of null(M) intersect k-perp;
- N4: the difference between total and transverse nullity;
- N6: a displayed null-space basis and residuals;
- N7: the returned basis count.

The two matrix routes share one action and test coding consistency only
(directives/S10_SHARED_PHYSICS.md:190-208). The basis-independent definition and
the display-only status of N6 are specified at
directives/S10_SHARED_PHYSICS.md:230-261.

## Measured MAIN spectrum and count

The raw transcripts give the following generic table. Here nu is N2 nullity and
nu_T is N3 transverse nullity.

| D | distinct roots | zero root: nu / nu_T | positive root: nu / nu_T |
|---:|---:|---:|---:|
| 2 | 2 | 1 / 0 | 1 / 1 |
| 3 | 2 | 1 / 0 | 2 / 2 |
| 4 | 2 | 1 / 0 | 3 / 3 |
| 5 | 2 | 1 / 0 | 4 / 4 |

SymPy evidence is at
scripts/out/S10_brane_mode_spectrum_sympy_audit.out:73-74, :105, :112,
:170, :177; :365-366, :397, :404, :462, :469; :676-677, :708, :715,
:773, :780; and :978-979, :1010, :1017, :1075, :1082.
Wolfram evidence is at
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out:26-27, :43,
:46, :56, :59; :239-240, :257, :260, :270, :273; :461-462, :480,
:483, :493, :496; and :676-677, :696, :699, :709, :712.

The raw root formula is 0 and μ_R |k|² / ρ_br. Under the supplied
positive-coefficient and nonzero-wavevector assumptions, the latter is
positive. The zero-root null direction is longitudinal and has no transverse
nullity; the positive root has D − 1 total and D − 1 transverse null
directions. This paired raw reading is not a blanket comparator pass for every
Q3 expression; the joined integer families below are the cross-engine result
claimed here.

The zero root retains a degree of freedom; it does not remove one. What the
curl-only stiffness removes is the restoring stiffness for the longitudinal
direction. The per-root nullities make the distinction explicit:
`1 + (D − 1) = D` in every MAIN dimension measured, so all D amplitude
directions remain in the spectrum. Within the light sector, Maxwell has no
counterpart to this surviving zero-frequency longitudinal direction. That is
a characterised departure, and its onward disposition belongs to S11
(`stray_longitudinal`). S10 assigns it no further interpretation.

The live name join supplies stronger, but still delimited, support for the
integer count. Across the 26 shared root cases, each of N2 rank, N2 nullity, N3
stacked rank, N3 transverse nullity, N4 nullity difference, and N7 basis count
has 26 PASS and 0 FAIL rows. Q3 root count has 13 PASS and 0 FAIL rows. These
are bare-integer comparisons, not empty-container passes. N6 canonical
null-space spans have 26 PASS and 0 FAIL rows; one literal example is
scripts/out/S10_cross_engine_comparator.out:447-456.

Both engines consume the same shared action and assumptions, so those
cross-engine passes test independent implementations of one stipulated model;
they do not independently validate the model.

The N7 residual itself is not a joined family. The engines emit
PY_S10_MAIN_D3_ROOT2_N7_BASIS_COUNT_RESIDUAL = 0 and
WL_S10_MAIN_D3_ROOT2_N7_COUNT_RESIDUAL = 0 at the raw-transcript loci
scripts/out/S10_brane_mode_spectrum_sympy_audit.out:492 and
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out:280, but the
names land in engine-only sets at comparator lines 3727 and 7271. This record
therefore does not claim cross-engine agreement for that residual.

Every number in the MAIN table is generic-rank evidence. The exceptional
qualification is measured separately below and may not be dropped when quoting
the generic result.

## Controls measured from their actions

The six live packages and their actions are defined at
directives/S10_SHARED_PHYSICS.md:421-472. MAIN is the D sweep. FULLGRAD and
DIVONLY change the stiffness form. SIGNFLIP and ANISO also change coefficients
or the kinetic form despite their XFORM spelling. XCOEF_SCALE changes a
dimensionless coefficient while retaining the curl form.

At D = 3 the emitted generic spectra and basis-independent counts are:

| package | roots / signs | measured nu / nu_T |
|---|---|---|
| MAIN | 0; +μ_R k²/ρ_br | zero 1/0; positive 2/2 |
| XFORM_FULLGRAD | +μ_R k²/ρ_br | positive 3/2 |
| XFORM_DIVONLY | 0; +μ_R k²/ρ_br | zero 2/2; positive 1/0 |
| XFORM_SIGNFLIP | 0; −μ_R k²/ρ_br | zero 1/0; negative 2/2 |
| XFORM_ANISO | three roots | nonzero transverse nullities 1 and 0 generically |
| XCOEF_SCALE | 0; +s μ_R k²/ρ_br | zero 1/0; positive 2/2 |

The corresponding SymPy transcript loci are :1285-1324 and :1436 for
FULLGRAD; :1726-1830 and :1960 for DIVONLY; :2339-2443 and :2573 for
SIGNFLIP; :2952-3121 and :3269 for ANISO; and :3973-4077 for XCOEF_SCALE in
scripts/out/S10_brane_mode_spectrum_sympy_audit.out.
The Wolfram loci are :892-902 and :928; :1238-1272 and :1313;
:1675-1709 and :1750; :2112-2186 and :2246; and :2802-2836 respectively in
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out.

This is a paired reading of raw control rows, not a claim that the comparator
passed every root or sign expression. The live comparator support asserted for
the table is limited to the integer families stated next; its unresolved sign
and Q6 rows are enumerated in the comparator section.

The controls show that the mode count is action-sensitive: the two form
controls alter the nullity structure, the sign flip preserves the count while
making the propagating root unstable, anisotropic inertia creates an additional
branch and an exceptional stratum, and a dimensionless coefficient rescaling
moves the positive root without moving the count.

The comparison that fixes the headline's scope is visible at every shared
form-control dimension:

| D | package | nonzero root | transverse nullity there |
|---:|---|---|---:|
| 3 | MAIN | μ_R |k|² / ρ_br | 2 |
| 3 | XFORM_FULLGRAD | μ_R |k|² / ρ_br | 2 |
| 3 | XFORM_DIVONLY | μ_R |k|² / ρ_br | 0 |
| 4 | MAIN | μ_R |k|² / ρ_br | 3 |
| 4 | XFORM_FULLGRAD | μ_R |k|² / ρ_br | 3 |
| 4 | XFORM_DIVONLY | μ_R |k|² / ρ_br | 0 |

SymPy emits those root-ordering and N3 operands at :366, :469, :677,
:780, :1286, :1324, :1511, :1549, :1727, :1830, :2038, and :2141 in
scripts/out/S10_brane_mode_spectrum_sympy_audit.out. Wolfram emits the
corresponding operands at :240, :273, :462, :496, :893, :902, :1070,
:1079, :1239, :1272, :1461, and :1495 in
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out. Their shared
root-ordering rows and integer N3 rows pass the live comparator. FULLGRAD also
has total nullity D at its sole root, whereas MAIN partitions total nullity as
1 at zero and D − 1 at the nonzero root. Thus curl-only determines the
missing longitudinal restoring stiffness, not the transverse count itself.

Q7 has a separate, narrower blind spot. The shared specification defines MAIN
and SIGNFLIP with the same `S_curl` stiffness density but opposite signs for
that density in the action. At D = 3 each raw engine emits byte-identical Q7
stiffness payloads for the two packages: SymPy at :548 and :2522, and Wolfram
at :295 and :1731. Q7 compares `S_pkg` with ordinary curl-squared, so it can
distinguish stiffness forms and normalisations but cannot distinguish the sign
with which the density enters the action. The Q7 tag families themselves are
not joined: SymPy spells them `Q7_STIFFNESS`, `Q7_CURL_DOT`, and
`Q7_DIFFERENCE`, while Wolfram spells them
`Q7_PACKAGE_STIFFNESS_DENSITY`, `Q7_ORDINARY_CURL_NORM`, and
`Q7_PACKAGE_STIFFNESS_VS_ORDINARY_CURL_RESIDUAL`. Their raw equality is
therefore neither comparator agreement nor disagreement; the shared Q3 sign
rows, not Q7, are what distinguish MAIN from SIGNFLIP.

The D = 3 Q7 rows also contain substantially less independent content than
their row count suggests. In each raw engine transcript the six package
stiffness payloads collapse to three distinct source-mapped classes: one curl
payload repeated by MAIN, SIGNFLIP, ANISO, and XCOEF_SCALE, one FULLGRAD
payload, and one DIVONLY payload. Of the six difference payloads, four are
identically zero; the two nonzero payloads, FULLGRAD and DIVONLY, are distinct.
Thus the twelve stiffness-and-difference rows in either engine carry five
informative raw-payload classes: three stiffness classes plus two nonzero
difference classes. These are groupings of committed source payloads, not live
comparator joins; Q7 does not join. The rows are at
scripts/out/S10_brane_mode_spectrum_sympy_audit.out:548-558, :1403-1413,
:1909-1919, :2522-2532, :3200-3210, and :4156-4166, and at
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out:295-297,
:917-919, :1294-1296, :1731-1733, :2215-2217, and :2858-2860.

Nor does the live join supply the gradient-symbol correspondence. Its only
automatic symbol normalization is mechanical snake-case to lower-camel, with
no exceptions (scripts/S10_cross_engine_comparator.py:75-78, :719); it applies
no `g{r}{c}` to `g{r}x{c}` map. A separate committed declaration records the
intended relation

    WL g{r}x{c} = PY g{r}{c} = d(u_c)/d(x_r),  r,c in 1..3,

on the authority that both engine definitions index the auxiliary symbol as
`(coordinate, field)`
(directives/S10_gradient_map_declaration_build.md:9-22). That source-level
declaration fixes the intended correspondence; it is not a value inference by
the current comparator and does not turn the Q7 rows into joined evidence.

The load-bearing control counts are included in the same 26 PASS / 0 FAIL
integer families reported above. Both engines nevertheless share the action
definitions and assumption set, so this cross-engine agreement remains
common-mode with respect to the stipulated physics.

## Prospective pre-registration

`steps/S10_PREREGISTERED_PREDICTION.md` was committed as `bc276485` at
2026-08-02 01:35:40 −06:00, before the first Wolfram engine commit
`a942c950` at 01:42:06 and the first SymPy engine commit `90c30ed8` at
01:52:54. It prospectively records the MAIN `D − 1` pattern, the FULLGRAD
one-root/nullity-`D` form control, the coefficient-rescaled roots with
unchanged nullities, and two explicitly least-certain expectations: the D = 3
curl normalization and whether D = 2 behaves in kind.

Every emitted prediction holds in both live transcripts. The MAIN table above
covers the four dimensions. FULLGRAD D = 3 has one root and total nullity 3 at
SymPy :1285 and :1317 and Wolfram :892 and :899. XCOEF_SCALE D = 3 has the
rescaled root and nullities 1 and 2 at SymPy :3972, :4005, and :4070 and
Wolfram :2800, :2820, and :2833. For the least-certain items, MAIN D = 3 has
equal Q7 operands and zero difference at SymPy :548, :553, and :558 and
Wolfram :295-297, while MAIN D = 2 has total/transverse nullities 1/0 and 1/1
at SymPy :105, :112, :170, and :177 and Wolfram :43, :46, :56, and :59.

This prospective ordering guards against rewriting the recorded expectations
after seeing the outputs. It does not establish target blindness or that an
engine builder could not read the file. In particular, the `1/2` convention
on which the Q7 normalization prediction turns is supplied by the shared
specification that both builders read, so the pre-registration adds no
independent route for that convention.

## Exceptional strata

ANISO has an allowed exceptional stratum. It changes the generic answer rather
than merely relabelling it:

| case | generic nonzero transverse nullities | allowed-stratum spectrum | allowed-stratum positive-root transverse nullity |
|---|---|---|---:|
| ANISO, D = 3 | 1 and 0 | 2 roots | 2 |
| ANISO, D = 4 | 2 and 0 | 2 roots | 3 |

SymPy emits the D = 3 allowed locus and stratum rerun at
scripts/out/S10_brane_mode_spectrum_sympy_audit.out:3269, :3289-3290,
:3305, :3312, :3348, :3355, and the D = 4 result at :3775, :3795-3796,
:3811, :3818, :3854, :3861. Wolfram emits the corresponding D = 3 result at
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out:2246,
:2258-2259, :2276, :2279, :2289, :2292, and D = 4 at :2601, :2613-2614,
:2631, :2634, :2644, :2647.

There is no live cross-engine stratum verdict. Measurement of the live join
found 275 Python stratum-tag names, 102 Wolfram stratum-tag names, and 0 shared
stratum-tag names. The shared specification itself records that the as-built
stratum names are not aligned at directives/S10_SHARED_PHYSICS.md:484-502.
Thus the table above is a paired reading of raw artifacts, not a comparator
pass, and the comparator establishes no exceptional-stratum result.

Both raw engines report the same exceptional movement, but they consume the
same shared action and assumptions. That common-mode limitation applies to the
exceptional table as well. Completeness is also not established inside either
engine: the specification assigns completeness to an orchestrator and says the
engine does not assert it (directives/S10_SHARED_PHYSICS.md:385-408). No live
artifact supplies that completeness verdict.

## Dimensional result

With the supplied displacement dimension [u] = (1,0,0) in
(length,time,mass), the two raw transcripts each emit

    [ρ_br] = (−D, 0, 1)
    [μ_R]  = (2−D, −2, 1)
    [μ_R] − [ρ_br] = (2, −2, 0).

SymPy emits these values for every MAIN D at
scripts/out/S10_brane_mode_spectrum_sympy_audit.out:223-225, :515-517,
:826-828, and :1128-1130. Wolfram emits them at
mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out:100-127,
:317-344, :537-564, and :753-780. This is a homogeneity result under the
supplied action and supplied [u], not an independent measurement of either.
This is a paired raw observation, not a claimed live-comparator verdict for a
three-object dimension family. The two engines also share the action and [u]
premise, so matching values would remain common-mode at those inputs.

No bulk quantity, ambient dimension, embedding dimension, or codimension
occurs in either S10 engine or raw transcript. The `D − 1` here is therefore
the corank after adding the wavevector row inside a D-dimensional in-plane
calculation; it is not a codimension computed from bulk and brane dimensions.
That observation has a strict bound. `SUBSTRATE_REQUIREMENTS.md` `R-S1-01`
records that absence of codimension from this mode-count calculation says
nothing about whether `D_brane` can be derived as a codimension elsewhere; it
leaves that object open for S6 and expressly identifies a codimension-one wall
in a four-dimensional bulk as a natural possible route. S10 neither supplies
nor excludes that later derivation.

The action-homogeneity booleans have a narrower scope than the rest of Q6.
For all 13 package-and-dimension runs, SymPy's
`Q1_LAGRANGIAN_EXPANDED_Q6_SOLVED_HOMOGENEITY` and Wolfram's
`Q6_SOLVED_ACTION_TERM_HOMOGENEITY` are `True`, but these are the action-term
booleans evaluated under dimensions solved by requiring those same action
terms to have the energy-density dimension. The engines mark that circular
scope explicitly: every SymPy `Q6_SOLVED_HOMOGENEITY_VACUOUS` is `True` with
operands `(6, 6, 0)`, and every Wolfram
`Q6_SOLVED_ACTION_HOMOGENEITY_VACUITY` is
`<|"EquationCountDifference" -> 0, "Vacuous" -> True|>`. Those marker
families do not join because their names differ. The qualification applies to
the solved action-term homogeneity booleans only. It does not downgrade
homogeneity on roots, null-space objects, minors, or Q7, and it does not
downgrade either the dimension solution or the cross-step overwrite
residuals; the shared specification says an inhomogeneity outside the solved
action scope is a real finding.

### How the control scales enter the dimension solves

All three live unknown-count residual rows carry the same route difference:
XCOEF_SCALE D = 3, ANISO D = 3, and ANISO D = 4 each compare SymPy 6 with
Wolfram 9 and emit residual −3
(scripts/out/S10_cross_engine_comparator.out:1065-1070, :1308-1313,
:1592-1597). The repeated residual is not three independent dimensional
findings. It is one construction difference exercised by every package row
that contains a control scale.

The premise surface is explicit in both raw transcripts. For SymPy, `s_rho`
appears as a standalone symbol on 148 transcript lines (493 occurrences) and
`s` on 44 lines (215 occurrences). The premise-bearing surface is four tags
for `s_rho`—the joint-assumption and dimension-premise tags at D = 3 and
D = 4—and two for `s` at D = 3
(scripts/out/S10_brane_mode_spectrum_sympy_audit.out:2894-2895,
:3409-3410, :3915-3916). For Wolfram, `sRho` appears on 138 lines (1081
occurrences) and `coefficientScale` on 31 lines (150 occurrences); again there
are four premise-bearing tags for `sRho` and two for `coefficientScale`
(mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out:2092, :2303,
:2432, :2658, :2782, :2881). The joint tags supply positivity and reality
(and, for ANISO, non-unit distinctness); each dimension-premise tag supplies a
zero dimension. The XCOEF_SCALE joint tags do not contain the current
`s ≠ 1` distinctness condition.

The engines then obtain that zero dimension by different routes:

- SymPy removes both control-scale symbols before forming its solver unknowns
  and installs each directly as `ZERO_DIM`
  (scripts/S10_brane_mode_spectrum_sympy_audit.py:448-470). In all three rows
  it emits six equations and six solution slots, all belonging to `rho_br` and
  `mu_R`; the scale has no equation or solution slot. The literal equations,
  solutions, counts, and slot lists are at SymPy :3168-3177, :3683-3692, and
  :4124-4133.
- Wolfram selects every coefficient symbol present in the action, so the
  control scale contributes three unknown slots, and it appends three
  zero-dimension control equations before solving
  (mathematica/S10_brane_mode_spectrum_mathematica_audit.wl:1002-1023).
  It therefore emits nine unknown slots and solutions for
  `dimSRhoLength/Time/Mass` or `dimScaleLength/Time/Mass`, each zero. The
  literal premise, equations, counts, and solutions are at Wolfram :2303-2313,
  :2658-2668, and :2881-2891.

Those routes must be judged by their build-time texts. The Wolfram inclusion
construction dates to `3bf61a6b` (2026-08-04 17:10:50 −06:00), and its
unknown-count diagnostic to `933398d6` (23:19:33). The then-current spec said
to emit the number of unknown coefficient dimensions and declared both scales
dimensionless, but did not say whether a declared dimensionless symbol still
counted as an unknown. The SymPy exclusion construction dates to `d6248be2`
(2026-08-05 02:59:54); its governing spec, preserved at `s10-as-built`, added
that the count came from the package's own action but still declared both
scales dimensionless and still did not settle that counting convention (the
`s10-as-built` version of directives/S10_SHARED_PHYSICS.md:432-436). The three residual
rows are therefore build-time specification-convention divergences, not
engine shortfalls against the texts the engines were built to.

The current rule was added later, in `e5a2c695` on 2026-08-06 06:57:13
−06:00: a declared dimensionless symbol is excluded, while `s_rho` is no
longer declared dimensionless and is included as a solved unknown
(directives/S10_SHARED_PHYSICS.md:443-454). Against that later rule, the
frozen SymPy ANISO premise and six-slot route, the frozen Wolfram ANISO
dimensionless premise, and the frozen Wolfram XCOEF_SCALE nine-slot route are
current-spec divergences caused by a post-build spec change. SymPy's
XCOEF_SCALE six-slot route matches the later rule. None of the postdated
divergences is evidence of a build-time engine shortfall. The same commit added
the `s ≠ 1` premise after both engines were built, so its absence from both
XCOEF_SCALE joint tags has the same postdated-spec disposition.

The same provenance qualification applies to the two current-spec findings
listed below. At `s10-as-built`, Q7 supplied the three ordinary-curl
components explicitly and did not require a Levi-Civita construction
(the `s10-as-built` version of directives/S10_SHARED_PHYSICS.md:369-376); the Levi-Civita
requirement was added by `e5a2c695` after the hand-typed engine constructions.
The as-built tag grammar likewise had no `STRATUM<s>` scope token
(the `s10-as-built` version of directives/S10_SHARED_PHYSICS.md:465-474); the aligned grammar
was added in that same later commit. These frozen behaviors do not conform to
the current spec, but neither was an engine shortfall against its governing
build-time text.

One premise has lost its only live falsifier. The supplied field dimension
`[u] = (1,0,0)` was formerly policed by comparing coefficient dimensions
derived from it with independently declared registry dimensions. The shared
specification and export-chain directive both identify that exact comparison
as the one place the premise could fail. The live tree has no reduction
directory and neither S10 source, transcript, nor export contains an S10
registry-dimension comparison. Consequently `[u]` is unfalsifiable within
this build. The antecedent is this one supplied premise, not the dimensional
layer as a whole; the emitted dimension solve and the other non-vacuous Q6
objects retain their stated evidential status.

## Export chain and what its guard can detect

The committed export modules contain 165 S9 records and 617 S10 records. All
165 S9 keys are present in S10. Of those, 162 records are field-for-field
unchanged. Exactly three shared records are author-flagged corroborated
overwrites:

| record | S9 value | S10 value | exact value equality |
|---|---|---|---|
| inertia_coefficient_dimension | (−D,0,1) | (−D,0,1) | true |
| stiffness_coefficient_dimension | (2−D,−2,1) | (2−D,−2,1) | true |
| coefficient_dimension_difference | (2,−2,0) | (2,−2,0) | true |

The other 452 records are S10-only. The S10 class tally is
KNOB 2, STRUCTURAL 4, COORDINATE 27, CONTROL 9, PREMISE 78, DERIVED 497.
The live S10 transcript prints the three overwrite operands and zero residuals
at scripts/out/S10_brane_mode_spectrum_sympy_audit.out:4215, the live count at
:4216, and the class tally at :4220. The published records are at
scripts/S10_exports.py:405-427; the S9 operands are at
scripts/S9_exports.py:395-414. S10 imports and binds S9 at
scripts/S10_brane_mode_spectrum_sympy_audit.py:31-46 and constructs the
overwrite check at :2077-2111.

Two scratch ablations measured the reach of this guard. They were made in
/tmp, not in the repository, by copying the S9 and S10 engines, both export
modules, and the export extractor. The scratch S10 export target was removed
before the measured S10 run, and S9 was run before S10.

Form ablation:

    main_action = construct_flexural_action().subs(mu_F, mu_R)

This changed the S9 stiffness dimension to (4−D,−2,1). S10 printed residuals
0, 1, 1 for inertia, stiffness, and their difference, then raised
AssertionError. Its exit status was 1 and no S10 export was written.

Coefficient ablation:

    main_action = construct_curl_action(rho_br * identity3, 7 * mu_R)

This moved the S9 speed display to 7 μ_R/ρ_br while leaving all three compared
dimensions unchanged. S10 printed three zero residuals, reported live count
617, exited 0, and wrote an S10 export.

Therefore the export comparison catches the measured form mutation because it
moves two of the three compared dimensions. It does not catch the measured
dimensionless coefficient mutation even though that mutation moves the
spectrum. It compares three cross-step SymPy values; it is not a second
cross-engine check. The assertion can also be stripped by optimized Python,
the author flag does not prove fresh derivation, and premises imported by both
steps can move both operands consistently. Those are hard limits of the live
guard, not hypothetical guarantees.

The compared objects are three exact SymPy value records classified `DERIVED`:
the inertia-coefficient dimension, stiffness-coefficient dimension, and their
difference. They are not mutually independent. S10 freshly derives its two
coefficient dimensions, but it imports from S9 the dimension symbol, the
coefficient symbols, the unit objects, and the supplied field-dimension
premise; its third record is then constructed from the first two. The guard is
therefore blind to the whole class of mutations that leave the three dimension
records equal on both sides. That class includes every dimension-preserving
action or coefficient mutation, not only the coefficient ablation that was
run, and every common-mode premise mutation accompanied by a consistent move
of the upstream records. The form ablation and coefficient ablation are one
caught member and one missed member of that class, not an exhaustive boundary.

To reproduce either ablation, copy the five named files to a scratch directory,
apply the single action mutation shown above to the scratch S9 engine, run the
scratch S9 engine to regenerate its export, remove only the scratch S10 export
target, and then run the scratch S10 engine. No ablated artifact is part of the
committed ledger.

## Live comparator: measured scope, not a blanket verdict

The comparator joins only mechanically matched names, compares sequences in
order, and canonicalizes null-space bases as spans
(scripts/out/S10_cross_engine_comparator.out:1-5). Its measured summary at
:9288-9313 is:

| quantity | count |
|---|---:|
| Python names | 4233 |
| Wolfram names | 2983 |
| shared names | 562 |
| compared shared names | 552 |
| unparsed shared names | 10 |
| agreements | 388 |
| bare-integer agreements | 215 |
| empty-container agreements | 14 |
| symbolic or structured agreements | 159 |
| disagreements | 164 |
| naming-only disagreements | 23 |
| representational disagreements | 13 |
| content divergences | 128 |
| shape or type mismatches | 109 |
| route-token divergences | 13 |
| numeric or algebraic residuals | 6 |
| null-space bases compared | 26 |
| duplicate rows | 0 |
| format issues | 10 |

For a value-comparable width, this record counts only 388 PASS rows plus the 6
rows with a genuine numeric or algebraic residual: 394. The remaining 168
shared names are 23 naming-only, 13 representational, 109 shape/type, 13
route-token, and 10 unparsed rows. They do not establish value agreement.
Likewise, the thousands of engine-only names are not cross-engine evidence.

The six genuine residual rows are one MAIN sign row and five control rows:

- MAIN D = 5 root 2 sign:
  1 − sign(k1²+k2²+k3²+k4²+k5²), at comparator :938-943;
- XCOEF_SCALE D = 3 Q6 unknown-coefficient count: −3, at :1065-1070;
- ANISO D = 3 Q6 unknown-coefficient count: −3, at :1308-1313;
- ANISO D = 3 root 3 sign:
  undecidedUnderJointAssumptions − 1, at :1498-1503;
- ANISO D = 4 Q6 unknown-coefficient count: −3, at :1592-1597;
- ANISO D = 4 root 3 sign:
  undecidedUnderJointAssumptions − sign(k1² sRho+k2²+k3²+k4²), at
  :1782-1787.

These are unresolved comparator rows. In particular, the MAIN D = 5 positive
sign is not cross-engine established even though the raw premise
Σ k_i² > 0 makes the mathematical sign positive; this record does not convert
that inference into a comparator pass.

Both engines still share the action and assumptions, so even the 388 agreement
rows are common-mode with respect to the physical premise. The comparator's
FINAL_GUARD is FAIL. Its source makes any disagreement, unparsed row,
duplicate, or format issue sufficient for failure
(scripts/S10_cross_engine_comparator.py:927-934). FAIL means the live joined
surface is not clean; it is not itself a refutation of the MAIN mode count.

The D12 worklist has two unrefuted symbol pairs, D ↔ braneDimension and
s ↔ coefficientScale, at comparator :9284-9287. Mechanical name matching has
no exception map, but injective symbol mapping is not established.

The comparator states of the families used in this record are as follows.
These states are not interchangeable: an engine-only name is neither agreement
nor disagreement, a shared row that fails on an unrefuted inner-symbol mapping
is a compared failure, and a shared passing row is a comparison result.

| family used here | live comparator disposition |
|---|---|
| Q1 action / Euler–Lagrange system | The action payloads have no shared name: Python uses `Q1_LAGRANGIAN_EXPANDED` and Wolfram uses `Q1_LAGRANGIAN`, so the join never compares them. All 13 Euler–Lagrange names are shared and compared, but fail representationally because Python emits expressions and Wolfram emits equalities. |
| Q2 matrices and route | `M_A`, `M_B`, and their residual each have 12 PASS rows; the XCOEF_SCALE row is shared and FAILS on the unrefuted `s` ↔ `coefficientScale` mapping. All 13 scalar entry-ratio rows pass. Both sources send downstream work through the quadratic-form matrix, but all 13 shared route-account rows fail as route-token content because the emitted tokens are `M_B` and `quadraticFormRoute`; a passing value family must not be quoted as a passing route account. |
| Q3 root expressions, counts, and signs | Root ordering is shared: 12 package-and-dimension rows pass, while XCOEF_SCALE is a shared naming-only failure on `s` ↔ `coefficientScale`. All 13 root-count rows pass. Of 26 sign rows, 23 pass and the three residual rows listed above fail. |
| N2, N3, N4, and N7 count | N2 rank, N2 nullity, N3 stacked rank, N3 transverse nullity, N4 nullity difference, and N7 basis count each have 26 shared PASS rows and no failed row. These are the joined integer families supporting the generic count. |
| N6 basis and its displayed diagnostics | The 26 `N6_NULLSPACE_BASIS` names are shared and their canonical spans pass. The displayed dot products and vector residuals have no shared name because Python uses `N6_BASIS_DOT_K` / `N6_BASIS_VECTOR_RESIDUALS` and Wolfram uses `N6_BASIS_DOTS` / `N6_BASIS_RESIDUALS`; those diagnostics are paired raw readings only. |
| N7 residual | No shared name: Python emits `N7_BASIS_COUNT_RESIDUAL` and Wolfram emits `N7_COUNT_RESIDUAL`. This is neither agreement nor disagreement. |
| Q6 displayed coefficient dimensions | The tuples displayed in the dimensional-result section have no shared family name: Python emits indexed coefficient `...Q6_DIMENSIONS` and a singular `DIMENSION_SOLUTION`, while Wolfram emits aggregate `...COEFFICIENT_DIMENSIONS` and plural `DIMENSION_SOLUTIONS`. The shared energy-density-dimension rows are 13 naming-only failures on D ↔ `braneDimension`. The shared unknown-count family has 10 passes and the three numeric failures listed above. The action-vacuity markers likewise have no shared name, for the spellings stated in the dimensional section. |
| Q7 form comparison | No shared name, because the two three-name vocabularies are the spellings stated in the controls section. The paired raw operands support the limited form comparison; the comparator supplies neither agreement nor disagreement for it. |
| Q8 allowed strata and stratum reruns | Eleven shared empty-list rows pass. The two ANISO allowed-strata rows share a name but are unparsed and fail without a value comparison. The reruns have no shared names because Python inserts `Q8_STRATUM1` while Wolfram inserts `STRATUM1`; hence no exceptional-stratum result is comparator-established. |
| cross-step export records | Outside the cross-engine join: all operands are SymPy records from the S9→S10 export chain. Their guard state and lack of independence are reported in the export section, not converted into a cross-engine state. |

At the package level, the same evidence reads:

| package | comparator state of the spectrum/count claims quoted here |
|---|---|
| MAIN | Root ordering and every generic integer count family pass. Signs pass except the D = 5 nonzero-root sign, which is a shared numeric/algebraic failure. All four empty allowed-strata rows pass. |
| XFORM_FULLGRAD | Root ordering, signs, every generic integer count family, and both empty allowed-strata rows pass. |
| XFORM_DIVONLY | Root ordering, signs, every generic integer count family, and both empty allowed-strata rows pass. |
| XFORM_SIGNFLIP | Root ordering, signs, every generic integer count family, and both empty allowed-strata rows pass. Q7 remains raw-only because its names do not join. |
| XFORM_ANISO | Root ordering and every generic integer count family pass. The root-3 sign fails at D = 3 and D = 4; both unknown-dimension counts fail; both allowed-strata rows are shared but unparsed; and the stratum reruns have no shared names. |
| XCOEF_SCALE | Every generic integer count and sign row passes, and its empty allowed-strata row passes. Its root-ordering row is a shared naming-only failure on `s` ↔ `coefficientScale`, as are its Q2 matrix rows; its unknown-dimension count is a shared numeric failure. |

## Claims this step still does not establish

Each item below is a current disposition, not inherited prose from the deleted
instrument.

1. **The in-plane/out-of-plane split remains supplied.** S10 does not derive or
   dynamically test the decoupling of u from out-of-plane fields.

2. **The physical value D = 3 remains owed.** The symbol D is a STRUCTURAL S9
   export imported by S10; no live reduction registry selects D = 3.

3. **The curl-only action remains supplied.** Absence of a longitudinal
   propagating mode follows from that action and does not explain why other
   allowed elastic terms are absent.

4. **Basis-span comparison is now established on its joined generic surface.**
   N6 canonical spans have 26 PASS and 0 FAIL rows, while N3 independently
   supplies the headline count. Both results remain common-mode with respect to
   the shared action and assumptions.

5. **The second matrix route remains only an internal consistency test.** The
   live comparator classifies 13 route-token rows as content divergences, but it
   does not compare the underlying route constructions. The shared specification
   expressly denies that the two routes are independent physics derivations.

6. **Exceptional-stratum cross-engine agreement remains open.** There are zero
   shared stratum names and no live completeness verdict, despite matching raw
   tables.

7. **The historical before/after Q7 failure count is not reproducible from the
   live tree.** The harness that produced it was deleted. Re-establishing that
   historical statement requires retrieving the old harness, its pinned engine
   outputs, and its configuration from version history and rerunning both
   versions. No historical count is carried forward here.

8. **Declared-name coverage is superseded, not inherited.** The live comparator
   uses mechanical matching without an exception map. Its current D12 result is
   the two-pair unrefuted worklist above; injectivity remains open.

9. **Package spelling is not physics classification.** SIGNFLIP and ANISO are
   not treated as pure form controls merely because their names begin XFORM.

10. **Whole-layer dimensional agreement remains open.** The three-record export
    comparison and its measured ablations do not certify every Q6 row in both
    engines; the live comparator has the three Q6 count residuals listed above.

11. **A green overall cross-engine result is not available.** Ten shared rows
    are unparsed, 164 compared rows disagree, format issues number 10, and the
    final guard fails. Engine-only and shape/type rows are not evidence.

12. **Current-spec conformance remains open at two points.** The current Q7
    specification requires Levi-Civita construction, but the live SymPy and
    Wolfram sources hand-type the three curl components at
    scripts/S10_brane_mode_spectrum_sympy_audit.py:1616-1637 and
    mathematica/S10_brane_mode_spectrum_mathematica_audit.wl:1295-1326.
    The current stratum grammar also is not aligned. These are evidence/spec
    findings; they do not license changing the frozen physics result here.

13. **The stated regime remains narrow.** The result assumes an unstrained
    medium at rest, no convective or dissipative terms, and a quadratic action.
    It does not establish the count for moving, strained, dissipative, or
    nonlinear backgrounds.

14. **Three premise-sensitive signs remain unresolved by the live join.** They
    are the MAIN D = 5 sign and the ANISO D = 3 and D = 4 root-3 signs listed
    above. A reader may not quote them as comparator-established signs.

15. **Finite thickness, microscopic discreteness, and frequency-dependent
    moduli are absent from the model.** Exact-token searches of the live shared
    specification and both action builders found no thickness or width
    variable, no lattice/site/spacing or other discrete-medium structure, and
    no frequency-dependent modulus. The shared action instead declares
    `rho_br` and `mu_R` as real constants
    (directives/S10_SHARED_PHYSICS.md:30-43), and the built actions use only
    those coefficients, the control scales, fields, and first derivatives
    (scripts/S10_brane_mode_spectrum_sympy_audit.py:378-400;
    mathematica/S10_brane_mode_spectrum_mathematica_audit.wl:150-190). These
    are measured absences, not physical idealisations that S10 took. The result
    was never asked to establish what happens with finite-thickness structure,
    a discrete medium, or dispersive moduli.

## Open requirements on interpreting the count as light

`SUBSTRATE_REQUIREMENTS.md` names seven OPEN entries whose source includes
S10. None is tested by the supplied S10 action. Their `on failure` text makes
the interpretation boundary explicit:

| entry | untested condition | what the register says fails |
|---|---|---|
| `R-S1-01` | delivery of `D_brane` | “Without a delivered `D_brane` the sentence is an assumption restated, not a result.” |
| `R-S6-02` | the unstrained rest reference state is an equilibrium of the substrate action | “`ω² = (μ_R/ρ_br)k²`, the mode count, and the dimensions all inherit the defect.” |
| `R-S8-01` | the substrate delivers the curl-only stiffness functional | “If S8 delivers any other form, S9's and S10's transverse results may survive while the no-longitudinal claim — Maxwell's third demand, and the whole reason the sector exists — does not.” |
| `R-S8-02` | the joint quadratic operator on in-plane `u` and out-of-plane `h`, including whether `h` is degenerate with the transverse pair | “S10's headline number changes”; block mixing alone is not the named mechanism. |
| `R-S8-03` | the substrate delivers the positive stiffness sign | “the transverse sector is two exponentially growing modes rather than two waves, and every mode count in S9 and S10 is unchanged.” |
| `R-S8-04` | the substructure supplies internal angular momentum or couple stress | “the curl-only stiffness functional is not an admissible continuum mechanics”; the counts would come from “a functional no medium can have.” |
| `R-S8-05` | the frame against which rotational stiffness is measured is admissible | if the reference is external, “every result in the sector inherits a preferred-orientation signature.” |

The quoted failure scopes are at `SUBSTRATE_REQUIREMENTS.md:103-119,
:151-163, :165-183, :185-210, :212-228, :230-254, and :256-275`. S10 need
not discharge these substrate obligations, but until they are discharged its
result is only a conditional mode count for the supplied quadratic action. It
does not establish a viable physical light sector.

## Prior art and departure

The relevant prior-art note says that the MacCullagh curl-energy construction
and its transverse count are standard. Its own search did not find the
thin-phase confinement architecture
(../../docs/medium_requirements_and_prior_art.md:127-156). That is a
search-failure statement, not an originality proof.

The supplied in-plane/out-of-plane split is not ours. The prior-art register
classifies “the `h`-branon as the brane's own
transverse fluctuation” as KNOWN and attributes that object to Cembranos,
Dobado, and Maroto, PRL 90 241301 (2003)
(../../docs/medium_requirements_and_prior_art.md:136-143). That entry
establishes attribution for the excluded out-of-plane field; it does not
establish this ledger's `u`/`h` decoupling, its quadratic operator, any
degeneracy with the in-plane modes, or the D − 1 count. Prior art is used here
as an attribution oracle, never as a premise for those uncomputed properties.

Accordingly:

- agreement with the ordinary D = 3 curl is an oracle check inside the standard
  MacCullagh regime, not a novel result and not a premise for the D sweep;
- sweeping D and recovering D − 1 is evidence for this implementation, not a
  novelty claim.

The broader project may make architectural claims, but S10 does not establish
them. The step plan places S10 in PHASE 2 — light and gives the surviving
longitudinal direction its onward home at S11. S10's permanent result is the
conditional, measured light-sector mode count and the limitations recorded
here (V3_STEP_PLAN.md:354-380).

## Registry disposition

The former reduction registry is absent: no tracked reduction files exist and
the reduction directory is absent. The live replacement is distributed:

- D is defined as an S9 STRUCTURAL export in scripts/S9_exports.py:31-36 and
  imported by S10 at scripts/S10_brane_mode_spectrum_sympy_audit.py:31-35;
- the per-D spectra and mode counts live in the indexed records of
  scripts/S10_exports.py, beginning with D = 2 at :1911, D = 3 at :2744,
  D = 4 at :3476, and D = 5 at :4208;
- the three cross-step dimension operands live at the export loci cited above;
- the physical selection D = 3 is owed and has no live registry entry.

The old provenance defect was a defect against a file that no longer exists.
This record neither cites that dead file nor treats its absence as evidence that
the physical dimension has been selected.
