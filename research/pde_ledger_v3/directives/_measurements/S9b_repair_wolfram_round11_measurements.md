# S9b Wolfram round 11 — final-code measurements

This packet uses only the final round-11 sources and their current executions. The engine is unchanged. The only harness repair is K10's observable gate, site inventory and domain witness; the scratch output destination is isolated for this round. The checks exercise synthetic gate factors and conditional-object readings. Input snapshots are scope evidence, not a source of measured results. This is an execution record, not a review verdict or universal equality claim.

## Artifacts and reproduction

- [Complete harness stdout](S9b_repair_wolfram_round11_harness.stdout.gz); [baseline stdout](S9b_repair_wolfram_round11_baseline.stdout.gz); [all engine captures and generated mutation/compact sources](S9b_repair_wolfram_round11_captures.json.gz).
- [Every attempt and density-sign trial](S9b_repair_wolfram_round11_witness_searches.txt); [site inventories, exclusions and search plans](S9b_repair_wolfram_round11_site_inventory.txt.gz); [reached K11 ordered pairs](S9b_repair_wolfram_round11_k11_reached_pairs.txt).
- [Literal checks stdout](S9b_repair_wolfram_round11_structural_checks.txt); [literal harness stderr](S9b_repair_wolfram_round11_harness.stderr.txt); [mechanical validation](S9b_repair_wolfram_round11_validation.json); [literal validator stdout](S9b_repair_wolfram_round11_validation_stdout.txt).
- [Mechanical validator](S9b_repair_wolfram_round11_validate.py); [packager](S9b_repair_wolfram_round11_package.py); [source hashes](S9b_repair_wolfram_round11_source_hashes.json); [handed-in source snapshots](S9b_repair_wolfram_round11_input_sources.json.gz); [artifact hashes and sizes](S9b_repair_wolfram_round11_artifacts.json).

The validator checks every printed operand against its actual capture, the vocabulary and local-name registry, unique single mutations, unchanged engine/other-knife source, K11 mechanical prefix extraction, pinned sources, all witness check flags and sign trials, fixture records, resources and portability. It invokes no CAS. Literal structural comparison in the harness is `SameQ` on `BinarySerialize` integer bytes of frozen parsed `HoldComplete` InputForm; no evaluation, simplification or canonicalization precedes a mark. Approximate-number encodings are printed. Marks are only **identical** or **not compared** and decide no mathematical equality or change.

```sh
python research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_validate.py s9b-r11-full
python research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_package.py
```

## Guarded commands and literal outputs

All commands use 8 GiB, zero swap, TasksMax 64 and pool s9b. Before each, the local driver polls `pgrep -x WolframKernel` until it is empty; only one kernel is launched. No wall-clock, CPU-time or inactivity deadline is imposed. The final harness/checks executions were awaited by polling until completion; no detached driver, watcher, queued resumption or recurring schedule was launched for them. Logs and final source pins remain under `_scratch/`.

### s9b-r11-checks-final

Working directory: `/var/projects/toy_physics-s9b-r11-wl`.

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r11-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r11-wl/_scratch/s11c/s9b-r11-checks-final -- math -script /var/projects/toy_physics-s9b-r11-wl/_scratch/s9b_wl_repair11/checks-pinned-final.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 1.5975352139994357,
  "unit": "s11c-guard-f5d3ea0ca9cc",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 1.5207850469996629,
    "childPid": 393030
  }
}
```

### s9b-r11-full

Working directory: `/var/projects/toy_physics-s9b-r11-wl`.

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r11-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r11-wl/_scratch/s11c/s9b-r11-full -- math -script /var/projects/toy_physics-s9b-r11-wl/_scratch/s9b_wl_repair11/harness-pinned-final.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 316.7162472399996,
  "unit": "s11c-guard-8bbcce25edf3",
  "stderrBytes": 351,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 316.646517444,
    "childPid": 393267
  }
}
```

### s9b-r11-engine-repo

Working directory: `/var/projects/toy_physics-s9b-r11-wl`.

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r11-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r11-wl/_scratch/s11c/s9b-r11-engine-repo -- math -script /var/projects/toy_physics-s9b-r11-wl/_scratch/s9b_wl_repair11/engine-pinned-final.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 16.147273796999798,
  "unit": "s11c-guard-4d39d728ca3c",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 16.071089791001214,
    "childPid": 387538
  }
}
```

### s9b-r11-engine-alone

Working directory: `/var/projects/toy_physics-s9b-r11-wl/_scratch/s9b_wl_repair11/engine-alone`.

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r11-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 8 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r11-wl/_scratch/s11c/s9b-r11-engine-alone -- math -script /var/projects/toy_physics-s9b-r11-wl/_scratch/s9b_wl_repair11/engine-alone/S9b_light_bending_mathematica_audit.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 16.05604638099976,
  "unit": "s11c-guard-489323f085ab",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.983144043999346,
    "childPid": 389017
  }
}
```


Guard admission refusals: **10**, all solely because available memory did not cover outstanding commitments plus this job. Admission was retried after waits of about three minutes at the unchanged 8 GiB limit; no CAS workload started in those refused attempts. Exact refusal tracebacks are in validation JSON. Guard kills, Wolfram memory exhaustion, task-limit events and item-12 events: **none**. Full-run peak cgroup memory: **2867449856 bytes**; maximum sampled tasks: **32**; maximum sampled swap: **0 bytes**. The cgroup and native address-space caps are **8589934592 bytes**. Per-run outcomes and final memory-event counters are recorded in validation JSON.

Harness stderr is preserved verbatim. It prints **3** `Limit::alimv` warning messages and a suppression notice: those limit calls ignore global assumptions involving the limit variable. The checks and both standalone engine runs have empty stderr. No warning is reported as a guard event or a method failure.

## Tags and search outcomes

**374 tags = 330 shared + 44 local**: [complete list](S9b_repair_wolfram_round11_tags.txt). No required shared name, local registry entry or item-10 output is missing. Every non-K11 copy emits all tags; the harness prints every baseline/corrupted pair and baseline-against-baseline control. K11 prints all full-baseline tags with reached/unreached status.

**271** baseline tags carry conditional `NOT_ESTABLISHED`: [baseline occurrences](S9b_repair_wolfram_round11_not_established_tags.txt), [harness pair occurrences](S9b_repair_wolfram_round11_harness_not_established_tags.txt), [all engine capture occurrences](S9b_repair_wolfram_round11_capture_not_established_tags.txt). These are literal conditional branches, not numerical decisions for an independently chosen ray.

| Knife | Search outcome | Evaluations (corrupted, baseline) |
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

Each search prints the union of all direct/transitive site parts across copies and every exclusion, independently of structural marks. K11 includes compact-unreached parts. Mathematical values/domains are distinguished from containers; branch sets and ordered first-match selections retain their engine readings. Attached conditions are carried through containing levels. Density labels partition the shared response object without imposing an extra condition. Absent entries are not emitted and receive no invented value. Additional exclusions are extracted symbol/function/row-name lists and literal replacement-source text. All readings and computations are in the inventories and attempts.

All searches stop at their recorded witness. Other enumerated parts receive no no-change decision. No sample, unresolved evaluation, agreeing point or method residual supplies a universal conclusion. K7 uses the engine's literal status text at `B_IMPLIED_JN_DEFLECTION["VelocityStatus"]`, with a printed CAS literal comparison. K6 evaluates the emitted local radar integral kernel and finite-endpoint source with its constructor conditions. Other numeric witnesses print all roles, exact profiles, points, fixing equations/solution intersections, attached conditions and finite evaluations.

K9a and K9b each try the shared object at concrete inputs in both density signs and print both trials before returning the first witness. The negative input has `f(r)=-9/r`, `rho0=1`, `r=3`, `n=3`; the positive input has `f(r)=1/(100r)`. Both have the supplied asymptotic density. The exact density evaluations, condition/value computations and profile certificates are in the attempt file and validation JSON. Test profiles enter no engine object and are not profile-family restrictions.

## K10 construction and actual witness

The single replacement is at `rayGate = propagationGate && traversalGate;` (engine **266**). `With` inserts the computed traversed-radius domain into a held pattern. In the already computed `propagationGate`, that pattern replaces the complete along-ray branch-existence `Inactive[ForAll]` by `Not` of itself. The quantifier body is retained, not negated pointwise. The enclosing ray-construction conjuncts and the independently computed traversal gate are retained. `propagationGate` itself is not assigned a new value.

The generated source has exactly this one replacement. All local tags, ray-domain construction and all base/forward `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE` payloads are byte-identical to baseline. The synthetic gate-factor fixture applies the actual manifest replacement and prints unchanged diagnostic factors and the selected negated quantifier. These unchanged diagnostics are excluded from K10's site inventory. Observable/reduction ray-gate domains and values depending on their attached domains remain included.

The witness site is **`WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE`**, payload **`{"OnDomain","conditional domain"}`**, reached through engine **266 → 271 → 287**. Its domain evaluations are **`{False, True}`**. The full traversal evaluations are **`{True, True}`**. Whole emitted object evaluations, corrupted then baseline:

```wl
<|"OnDomain" -> Undefined, "OutsideDomain" -> "NOT_ESTABLISHED", "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>
<|"OnDomain" -> False, "OutsideDomain" -> Undefined, "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>
```

Whole-object evaluation first resolves each scoped conditional domain, then evaluates its value only on a true domain. `OnDomain`'s false condition supplies a decided absence (`Undefined`); its complementary `OutsideDomain` supplies `NOT_ESTABLISHED`. The association introduces no additional condition. This is a witness in an emitted observable domain, not in the unchanged branch diagnostic or a local pre-gate predicate.

Witness-only profiles are flat: delta, velocity, xi, bulk fraction and exchange zero; brane density one, stiffness four and light speed two. The printed independent point has `GM=0`, `b=2`, `bFar=1`, `ZE=ZR=2`, with every remaining input and role shown. Each copy's actual ray construction is independently solved. Turning candidates are filtered by their full domains; ordinary endpoint integrals are evaluated on the derived admissible interior and the closed upper endpoint is checked separately. Active quantifiers retain the full traversed ranges. Complete real solution sets and their intersection are computed; turning radii and radar impact are fixed by those equations, not independently assigned. The intersection is:

```wl
rTurnFlyby == 2 && rTurnRadar == 2 && impactRadar == 2
```

Both copies use that common solution. Every other attached condition and the traversal gate is computed true, and both emitted objects are evaluated. The attempt includes scoped bound variables, equations, solution selection and all CAS operations.

## Item 13, bindings and restrictions

General live radial functions and integral functionals represent the profiles. Every-b reduction uses exterior-half-line Abel inversion, a constant-GM anchor and signed squared-velocity branches carrying integration-constant domains. Coordinate Cartesian component differentiation computes mass divergence before radial evaluation. No fold/LEDGER binding, flux-to-GM relation or premise is introduced.

There is no engine profile-family restriction. Item-10 flow-only sets delta/xi to zero; speed-only sets velocity/xi to zero; tilt-only sets delta/velocity to zero; bulk-only sets velocity/xi to zero. The forward case uses its labelled `jn=0` premise and flux-defined velocity, first with delta/xi zero, then live. These additional restricted outputs do not replace general-profile objects. Part C applies its supplied responses and first-order f counting, retaining the exact shared response/domain at each response site.

Supplied-input/spec-line bindings in `S9b_SHARED_PHYSICS.md`: dispersion/metric **45**; optical ratio **52**; polarizations **59**; advection **64**; mass balance **69**, coordinate measure **72–74**; embedding **87**; historical uniform anchor **95**; LAB_HELD/background anchoring **103–104**; reference speed **116**; delta/nonnegative root **124**; retained grades **128**, optical counting **135,138–145**; independent GM **202–204**; references **208–210**; pressure, sound square and bulk fraction/asymptote **219–224**; Part-C responses **244–246**. Directive item 9 supplies root, zero-power and density-sign conventions; item 10 supplies the labelled forward premise. Part D is not built.

| Knife | Engine line | Original construction |
|---|---:|---|
| K1 | 60 | `kineticSymbol = (omega - vVector.kCov)^2;` |
| K2 | 58 | `metric2 = DiagonalMatrix[{aMetric, r^2}];` |
| K3 | 61 | `dispersionSpeed = cSquared;` |
| K4a | 103 | `rayDerivative[z_, x_] := D[z, x];` |
| K4b | 103 | `rayDerivative[z_, x_] := D[z, x];` |
| K4c | 52 | `xiSlope = Simplify[D[xiAnsatz, r], $Assumptions];` |
| K5 | 54 | `returnDirection = -1;` |
| K6 | 356 | `radarSource = roundTrip;` |
| K7 | 395 | `quantifierMode = "EVERY_B";` |
| K8 | 466 | `massDensity = rhoBr[r];` |
| K9a | 94 | `localSoundSquared = bulkSoundSquared /. rho -> rho0 (1 + fLocal);` |
| K9b | 95 | `localBulkDensity = rho0 (1 + fLocal);` |
| K10 | 266 | `rayGate = propagationGate && traversalGate;` |
| K11 | 59 | `vVector = {vLocal, 0};` |

Every construction is present. Exact replacements and unique occurrence counts are preserved in the full SITE records and validation JSON.

## K11 and executable checks

K11 uses exact mechanical extraction through `COMPACT_K11_BOUNDARY`, with identical configuration for corrupted and unmutated copies. It retains untruncated local dispersion/Fermat and the twelve-grade one-form ring, with no added truncation, samples, seeds or inexact arithmetic. Coverage is not an attached condition or a universal equality claim.

**10 reached tags**: [complete list](S9b_repair_wolfram_round11_k11_evaluated_tags.txt). Each prints the full baseline, compact-corrupted/compact-unmutated difference pair and compact-unmutated/full-baseline method residual pair, with structural marks and coverage. The method residual is not a witness. The full corrupted K11 engine is not run.

**364 not-evaluated tags**: [complete list](S9b_repair_wolfram_round11_k11_not_evaluated_tags.txt). Both pairs are explicitly not evaluated/not compared; no no-change decision is made. The inventory includes their downstream site parts.

The checks script prints **25** records before validation: syntax and structural fixtures, exact approximate-number encodings, domains, binders, ordered cases, common/disjoint solution sets, semantic branch identities, literal text/list readings, density-sign/asymptote bookkeeping, whole-quantifier versus pointwise negation, and the actual K10 replacement on synthetic gate factors. There are no added expected physics values.

The repository and copied-alone engine runs exit 0, produce nonempty byte-identical streams, and also match the full harness baseline byte for byte. Only the copied `.wl` remains in the standalone directory. Capture growth was observed while the kernel was alive; observations are in validation JSON. Production-out `git status --porcelain` is empty. No other checkout or orchestrator memory was read. No commit, push, other AI, reviewer or agent was launched. Build/report work stops here; interpretation/review remains with the orchestrator.
