# S9b Wolfram round 10 — current scoped-repair 8 GiB measurements

This packet describes only the current round-10 final-code runs. The scoped repair moves K10 to the completed branch condition, reads the source-written Profiles list literally, and corrects the round-8 historical note. The engine is unchanged. The harness and checks match their pinned final sources. No superseded run supplies a result. The full harness, checks and both engine portability demonstrations finished with exit code 0. This is an execution record, not a review verdict or a universal equality claim.

## Current sources and literal artifacts

- [Engine](../../mathematica/S9b_light_bending_mathematica_audit.wl), [harness](../../mathematica/S9b_light_bending_mathematica_ablation.wl), [checks](../../mathematica/S9b_light_bending_mathematica_ablation_checks.wl).
- [Repair input source snapshots](S9b_repair_wolfram_round10_input_sources.json.gz) preserve the handed-in harness/checks and note for the scope comparison.
- [Source hashes](S9b_repair_wolfram_round10_source_hashes.json), [mechanical validator](S9b_repair_wolfram_round10_validate.py), [validation JSON](S9b_repair_wolfram_round10_validation.json), [literal validation output](S9b_repair_wolfram_round10_validation_stdout.txt).
- [Complete harness stdout, gzip](S9b_repair_wolfram_round10_harness.stdout.gz), [baseline stdout, gzip](S9b_repair_wolfram_round10_baseline.stdout.gz), [all current captures and single-mutation/compact source copies, gzip JSON](S9b_repair_wolfram_round10_captures.json.gz), [artifact hashes and sizes](S9b_repair_wolfram_round10_artifacts.json).
- [Every actual attempt and density-sign trial](S9b_repair_wolfram_round10_witness_searches.txt), [site-part unions, exclusions, readings and search plans, gzip](S9b_repair_wolfram_round10_site_inventory.txt.gz), [reached K11 ordered pairs](S9b_repair_wolfram_round10_k11_reached_pairs.txt), [literal checks output](S9b_repair_wolfram_round10_structural_checks.txt), [literal harness stderr](S9b_repair_wolfram_round10_harness.stderr.txt).

The validator also restores the named edits and checks every remaining witness-method byte against the round-9 input; the other knife manifest entries, literal comparison and K11 method are unchanged. It checks every raw pair operand against the current engine captures, vocabulary and local-name registry, unique single source replacements, K11 mechanical prefix extraction, source pinning, final exit/limit records, all printed witness check flags and sign trials, fixture records and engine portability. It launches no CAS and makes no mathematical comparison. Original payload syntax is preserved in the raw streams.

Structural marks use `SameQ` on integer byte lists from `BinarySerialize` of parsed `HoldComplete` emitted InputForm, before evaluation, simplification or canonicalization. The harness also prints approximate-number byte encodings. A mark is only **identical** or **not compared**. It does not decide mathematical equality or change, and the harness prints no verdict. Search decisions are separate.

## Current commands and literal guard output

All jobs use pool `s9b`, a **8 GiB** cap, zero swap and TasksMax **64**, one Wolfram kernel at a time. No wall-clock, CPU-time or inactivity deadline was imposed. The local driver polls until completion and writes durable scratch logs. No job was retried.

Checks script (cwd `/var/projects/toy_physics-s9b-r8-wl`):

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r10-checks-final -- math -script /var/projects/toy_physics-s9b-r8-wl/_scratch/s9b_wl_repair10/checks-pinned-final.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 1.5388288290123455,
  "unit": "s11c-guard-20340721fac9",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 1.468697472999338,
    "childPid": 1174575
  }
}
```

Full harness (cwd `/var/projects/toy_physics-s9b-r8-wl`):

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r10-full -- math -script /var/projects/toy_physics-s9b-r8-wl/_scratch/s9b_wl_repair10/harness-pinned-final.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 299.79964380303863,
  "unit": "s11c-guard-6b52d49ca478",
  "stderrBytes": 351,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 299.7361351440195,
    "childPid": 1175006
  }
}
```

Engine in repository (cwd `/var/projects/toy_physics-s9b-r8-wl`):

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r10-engine-repo -- math -script /var/projects/toy_physics-s9b-r8-wl/_scratch/s9b_wl_repair10/engine-pinned-final.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 15.74101212003734,
  "unit": "s11c-guard-104d1b141c1e",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.672229306015652,
    "childPid": 1193756
  }
}
```

Engine copied alone into empty scratch directory (cwd `/var/projects/toy_physics-s9b-r8-wl/_scratch/s9b_wl_repair10/engine-alone`):

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r10-engine-alone -- math -script /var/projects/toy_physics-s9b-r8-wl/_scratch/s9b_wl_repair10/engine-alone/S9b_light_bending_mathematica_audit.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 15.736352089967113,
  "unit": "s11c-guard-7d4a6dd0680f",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.672553798998706,
    "childPid": 1195212
  }
}
```

All four runs exited **0**; no guard refusal/kill, memory-exhaustion, task-limit or item-12 event occurred. Full-run peak cgroup memory **2867781632 bytes**, maximum sampled tasks **32**, sampled swap **0 bytes**. These are cgroup memory figures, not a native virtual-address-space peak. Effective cgroup and native address-space limits were both **8589934592 bytes**. Per-run resource summaries and final memory-event counters are in validation JSON.

The unchanged harness emitted **3 `Limit::alimv` warnings** (**351 stderr bytes**, plus a suppression notice): global assumptions mentioning the limit variable were ignored by those limit calls. Three warning messages are printed before Wolfram suppresses further messages of that type. The resulting computations and actual outputs are preserved. The checks and both standalone engine runs have zero stderr bytes. No warning is represented as a new source repair or a guard kill.

Mechanical validation and packaging:

```sh
python research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round10_validate.py s9b-r10-full
python _scratch/s9b_wl_repair10/package.py
```

## Tags and actual search results

**374 tags = 330 shared + 44 local:** [complete tag list](S9b_repair_wolfram_round10_tags.txt). No required shared name is missing; no extra shared name or duplicate is present. The local-name registry exactly matches local emissions. Every non-K11 knife copy emits all tags. The full harness prints all tags for every knife and the baseline-against-baseline control, including K11 not-evaluated records.

**271** baseline tags carry conditional `NOT_ESTABLISHED`: [baseline list](S9b_repair_wolfram_round10_not_established_tags.txt), [all structural harness occurrences](S9b_repair_wolfram_round10_harness_not_established_tags.txt), [all engine-capture occurrences](S9b_repair_wolfram_round10_capture_not_established_tags.txt). These count literal conditional branches, not numerical decisions about a chosen full ray.

Each knife's search first prints its direct/transitive site-part union and exclusions, including engine-local parts and both item-9 density signs. Unioning is independent of structural marks. K11 includes downstream full-baseline parts the compact method does not reach. Named mathematical values/domains are parts; associations and branch/wrapper structures are containers. Response rows, density signs, coordinate bases and defining branch meanings supply semantic location identities. Conditional values/domains and ordered first-match selections are explicitly read; branch lists are read as branch sets. Missing entries are printed as not emitted without inventing a value. The two additional exclusions are extracted symbol/function/row-name lists and replacement-source text differences. The complete actual inventories and readings are in the linked gzip file; counts are in validation JSON.

Searches stop at their first recorded witness; all other parts were enumerated but not tried. No mathematical no-change or cannot-be-tried decision is reported from a point. No global equality is inferred from agreeing values, samples, structural marks or the K11 method residual. The unused unsuccessful-search branch is not described as having run.

| Knife | Actual search outcome | Emitted evaluations (corrupted, baseline) |
|---|---|---|
| K1 | witness | `{-2249111/810000, -1279808/455625}` |
| K2 | witness | `{-21067/7500, -1279808/455625}` |
| K3 | witness | `{-5070551/1822500, -1279808/455625}` |
| K4a | witness | `{2917376911/2916000000, 323071723/324000000}` |
| K4b | witness | `{1937902843/1944000000, 323071723/324000000}` |
| K4c | witness | `{0, -1/90}` |
| K5 | witness | `{0, -1207049/145800000}` |
| K6 | witness | `{-1/450*1/Sqrt[5], -1/225*1/Sqrt[5]}` |
| K7 | witness | `{"velocity undetermined by one fixed-b condition", "velocity determined on signed components with smooth zero-set joining"}` |
| K8 | witness | `{2/135, 1/90}` |
| K9a | witness | `{2, 4}` |
| K9b | witness | `{2, 8}` |
| K10 | witness | `{False, True}` |
| K11 | witness | `{-5238707/1822500, -1279808/455625}` |

K7's witness uses the amended literal-text rule. Its site is `WL_S9B_B_IMPLIED_JN_DEFLECTION["VelocityStatus"]`, copied from the engine's computed solution status at line 575. Both emitted literal statuses, their source branches and the Wolfram `SameQ` decision are printed. It is not a value/domain proof for a reduced every-b mathematical condition. It is neither a list of extracted names nor the literal EVERY_B/FIXED_B replacement source.

K6 evaluates the computed grade-100 local radar integral kernel, its finite-endpoint source and the printed constructor conditions. This local container has no attached full-ray gate. Actual integral radii remain bound during profile/point substitution. Both frequency constructions, their joint solution sets/intersection and the common selected solution are printed.

K10 compares the complete `WL_S9B_BRANCH_EXISTENCE` condition. Its single source replacement wraps the entire `propagationGate = rayConstruction && ForAll[...]` result at engine 262–264 in `Not`; the local predicate at 166, the quantifier body and `PATH_TRAVERSAL` are unchanged. The K10 inventory is derived from this later entry: direct branch/gate-domain fields and inherited attached conditions; upstream local objects and branch-type/traversal objects are excluded as independent.

The K10 attempt chooses flat test profiles only for the witness: delta, velocity, xi and bulk-fraction zero; brane density one, stiffness four, light speed two, and exchange zero. It sets `GM=0`, `b=2`, `bFar=1`, `ZE=ZR=2`, with all remaining free inputs printed. Both copies' original ray constructors are independently read and solved. Algebraic turning candidates are filtered by their complete constructor domains. The endpoint integrals are evaluated ordinarily on the derived admissible interior; the closed upper endpoint is separately evaluated and excluded by the endpoint equation, so no boundary solution is discarded. Active quantifiers retain their whole traversed-radius ranges.

Both complete real solution sets give `rTurnFlyby=rTurnRadar=impactRadar=2`; their intersection is computed, and that common solution is used. The two full emitted branch conditions evaluate **False, True** (corrupted, baseline); source assumptions and the independent full traversal gate evaluate true in both. No turning/impact parameter is an extra independent assignment, and no local pre-quantifier predicate furnishes the witness. All fixing statements, candidate domains, ordinary integrals, active decisions, roles and complete object evaluations are printed.

`Profiles` is a literal list of five profile symbols written at engine 401. It is a descriptive value, not a condition, not text, and not a list extracted from an expression. Its elements receive no profile-function substitution. The generic name-list exclusion now applies only to actually extracted fields. The current union includes one-copy literal Profiles parts where reached; no absent value is invented. The checks print this actual literal reading and a synthetic whole-quantifier-versus-pointwise-negation demonstration.

K9a and K9b each perform a combined trial of the shared response object on both density signs, with actual object evaluations, not just label tests. Both sign computations are printed before returning the first witness of that combined trial. The negative profile is `f(r)=-9/r`, with `rho0=1`, `r=3`, `n=3`; the positive profile is `f(r)=1/(100r)` at the printed point. Both approach the supplied asymptotic bulk density. Each trial prints all profile functions and required derivatives, concrete input assignments, density in both copies, sign decision, exact fixing-solution intersection and condition/value evaluations.

Actual K9a sign records:

```json
[
  {
    "sign": "\"rho <= 0\"",
    "densities": "{-2, -2}",
    "evaluations": "{2, 4}",
    "result": "\"witness\""
  },
  {
    "sign": "\"rho > 0\"",
    "densities": "{301/300, 301/300}",
    "evaluations": "{2, Sqrt[301/3]/5}",
    "result": "\"witness\""
  }
]
```

Actual K9b sign records:

```json
[
  {
    "sign": "\"rho <= 0\"",
    "densities": "{-2, -2}",
    "evaluations": "{2, 8}",
    "result": "\"witness\""
  },
  {
    "sign": "\"rho > 0\"",
    "densities": "{301/300, 301/300}",
    "evaluations": "{2, 90601/45000}",
    "result": "\"witness\""
  }
]
```

The remaining numeric/domain witnesses print their semantic payload location, source entry, all assigned/fixed/bound roles, exact profiles and inputs, conditions, operations and evaluations. Frequency solutions are paired at common solutions whenever the computed intersection is nonempty; distinct solutions are used only when the intersection is decided empty. These test profiles enter no engine object and are not engine profile families.

## Item 13, bindings, restrictions and knife sites

General live radial functions and integral functionals represent the profiles. Every-b reduction uses exterior-half-line Abel inversion, a constant-GM anchor and signed squared-velocity branches with integration-constant domains. The mass balance uses Cartesian component differentiation before radial evaluation. No new fold/LEDGER binding, flux-to-GM relation or premise is introduced. No shared-vocabulary or item-10 object is omitted.

There is no engine profile-family restriction. Labelled item-10 substitutions retain their defined scope: flow-only sets delta/xi to zero; speed-only sets velocity/xi to zero; tilt-only sets delta/velocity to zero; bulk-only sets velocity/xi to zero. The forward case uses its labelled `jn=0` premise and flux-defined velocity, first with delta/xi zero and then live. These are additional labelled outputs, not replacements for general-profile objects. Part C retains its supplied response substitutions and first-order f counting, with exact shared response/domain objects at their response sites.

Supplied-input/spec-line bindings in `S9b_SHARED_PHYSICS.md`: dispersion/metric **45**; optical ratio **52**; polarizations **59**; advection **64**; mass balance **69** and coordinate measure **72–74**; embedding **87**; historical uniform anchor **95**; LAB_HELD/background anchoring **103–104**; reference speed **116**; delta/nonnegative root **124**; retained grades **128** and optical counting **135,138–145**; independent GM **202–204**; comparison references **208–210**; pressure, sound square and density fraction/asymptote **219–224**; Part-C responses **244–246**. Directive item 9 supplies root, zero-power and density-sign conventions; item 10 supplies its labelled forward premise. No Part D input is used.

| Knife | Actual construction site | Engine line |
|---|---|---:|
| K1 | advective kinetic shift | 60 |
| K2 | local two-dimensional metric | 58 |
| K3 | speed square in dispersion | 61 |
| K4a | differentiated delta in rayDerivative | 103 |
| K4b | differentiated velocity in rayDerivative | 103 |
| K4c | differentiated xi ansatz | 52 |
| K5 | return-leg direction | 54 |
| K6 | source whose radar slope is taken | 356 |
| K7 | every-b/fixed-b reduction selection | 395 |
| K8 | density in coordinate mass balance | 464 |
| K9a | local sound square in fixed-ratio response | 94 |
| K9b | local density in power response | 95 |
| K10 | completed propagationGate after whole-condition logical negation | 262–264 |
| K11 | azimuthal component of advected velocity | 59 |

All constructions are present. Every single mutation's original/replacement, source line and unique occurrence are printed in its SITE record; the validator checks the generated source contains exactly that one intended replacement.

## K11 and executable checks

K11 uses the unchanged exact symbolic mechanical prefix extraction through `COMPACT_K11_BOUNDARY`, applied with identical configuration to corrupted and unmutated copies. The mutation remains in force. Its coverage is untruncated local dispersion/Fermat and the supplied twelve-grade one-form ring, with no added truncation, samples, seeds or inexact arithmetic. Coverage metadata is not an extra attached domain or a universal no-change claim.

**10 reached tags:** [complete list](S9b_repair_wolfram_round10_k11_evaluated_tags.txt). Each prints the full baseline, ordered compact-corrupted/compact-unmutated difference pair, and ordered compact-unmutated/full-baseline method residual pair, each with its structural mark and coverage. The full baseline requires no additional configuration reduction beyond the same engine retained ring. The method residual cannot furnish a witness. The full corrupted K11 engine is not executed.

**364 not-evaluated tags:** [complete list](S9b_repair_wolfram_round10_k11_not_evaluated_tags.txt). Each explicitly prints not evaluated and not compared for both pairs, with no no-change decision. The site-part inventory still includes downstream full-baseline parts, with compact reachability recorded separately.

The checks script prints **23 output records**, including parsing and synthetic literal/reading fixtures, before the mechanical validator can inspect them. These are not expected physics values. The tests additionally print the whole-condition K10 replacement and literal Profiles reading. They retain literal approximate-number distinctions, different domains, unevaluated syntax, ordered case selection, common/disjoint fixed solution sets, decided/unresolved domain booleans, scoped bound variables through associations, semantic branch identities and density-sign profile/asymptote bookkeeping.

The engine was run in the repository and, byte-for-byte copied alone, in an initially empty scratch directory. Both streams are nonempty and byte-identical to the current full baseline. The standalone directory afterwards contains only the copied `.wl` file. The current full-run capture grew while its kernel was alive; actual observations are in validation JSON. Production-out git status is empty.

The round-8 note is corrected: round-9 amendment-12 runs completed. That historical result supplies no round-10 measurement. All required source and result files remain uncommitted. No commit, push, AI agent or reviewer was launched, and no other checkout or orchestrator memory was read. The requested execution and measurement work is complete; interpretation/review is left to the orchestrator.
