# S9b Wolfram repair round 8 — measurements

Builder: Codex `gpt-6.1-sol`, xhigh effort, new harness author. Amendment 11 is the governing repair; the engine was not changed. No AI/reviewer/agent was launched and nothing was committed or pushed. The only repository accessed was `/var/projects/toy_physics-s9b-r8-wl`.

## Final code and literal output

- [Harness](../../mathematica/S9b_light_bending_mathematica_ablation.wl), [checks](../../mathematica/S9b_light_bending_mathematica_ablation_checks.wl), [unchanged engine](../../mathematica/S9b_light_bending_mathematica_audit.wl).
- [Source hashes](S9b_repair_wolfram_round8_source_hashes.json), [mechanical validation script](S9b_repair_wolfram_round8_validate.py), [validation](S9b_repair_wolfram_round8_validation.json), [literal validation stdout](S9b_repair_wolfram_round8_validation_stdout.txt).
- [Complete literal harness stdout, gzip](S9b_repair_wolfram_round8_harness.stdout.gz), [literal baseline stdout, gzip](S9b_repair_wolfram_round8_baseline.stdout.gz), [all engine captures and mutated/compact source copies](S9b_repair_wolfram_round8_captures.json.gz). These contain output from the final harness code. [Archive sizes, uncompressed stdout hashes and resource summary](S9b_repair_wolfram_round8_artifacts.json).
- [Every witness search, literal stdout](S9b_repair_wolfram_round8_witness_searches.txt), [all reached K11 operand pairs](S9b_repair_wolfram_round8_k11_reached_pairs.txt), [structural checks, literal stdout](S9b_repair_wolfram_round8_structural_checks.txt).

The per-tag mark is **identical** or **not compared**. Its named comparison is `SameQ` on integer bytes from `BinarySerialize` of `HoldComplete` parsed emitted InputForm. The parser's third argument freezes syntax before evaluation. There is no normalization, mathematical equality test or old difference classification in this path. Both full original emitted operands are copied verbatim into each ordinary row. Approximate operands have literal byte encodings printed as well; the executable examples distinguish adjacent machine numbers and different precisions. Polynomial rewrites, inert zero integrals and renamed binders remain not compared. The baseline-against-baseline control exercises all emitted tags.

A witness search is performed once per knife, separately from these marks. Successful local value witnesses print the payload/source locations, scalar point, function-valued profiles, consistent derivatives, active smoothness decisions, conditions and their exact evaluations, fixing equations, solution sets, selected roots and computed solution intersection. The two-copy ordering throughout witness records is **corrupted, baseline**. Ordered case readings include the default, the negation of every earlier condition, their exact evaluations and the selected effective condition. The joint solution sets include fixed auxiliary parameters as well as frequency; the response records include the fixed local density-fraction parameter. Frequency-derived parts compare at their common selected solution whenever the solution intersection is nonempty; different solutions are paired only when exact `Reduce` decides the intersection empty. Parts with no solved dependent quantity print the common empty assignment explicitly. All successful evaluations use exact real algebraic arithmetic, with error bound `0`. No successful candidate contains an unevaluated integral, limit or quantifier.

The successful witnesses concern the printed local construction parts, not the later gated ray observables. Only symbols and profiles read by the selected part, its conditions or its printed fixing equations are inputs to that witness; additional examples in the shared profile certificate are unused there. The only search configuration uses the printed rational radial functions on `q > 0`. It verifies derivatives through order four and the supplied optical tail counting by active computations. K9a/b witnesses lie in the printed `rho > 0` partition; density-sign labels are not treated as attached conditions. The local response parameter `fLocal` is fixed by the same `fProfile[r]` input in both copies. The full output retains both density partitions for every response copy.

K6, K7 and K10 print **none found**. Their candidate parts and the endpoint/turning construction after the test-profile substitution are printed, together with the unresolved dependent radii and radar impact. The search obtained no exact solution of that whole ray construction and its ranged conditions. It assigns no guessed turning radius, substitutes no bound variable, treats no inactive condition as decided, and makes no claim of mathematical equality or change for these searches. These finite unsuccessful searches do not change the engine or introduce a second construction for any physics object.

## Item 17 inventory

The full [tag list](S9b_repair_wolfram_round8_tags.txt) is **374 tags: 330 shared and 44 engine-local**. The local-name registry is complete. No shared-vocabulary name is missing, and no item 10 object is omitted. Every knife/subknife and the baseline-against-baseline control prints every tag. Marks and witness summaries per knife are in the validation JSON.

Every baseline occurrence of conditional `NOT_ESTABLISHED` is covered by the [271-tag list](S9b_repair_wolfram_round8_not_established_tags.txt); [all harness rows carrying it](S9b_repair_wolfram_round8_harness_not_established_tags.txt) are listed separately. These are conditional output branches, not claims that 271 observables were evaluated outside their domains.

Item 13 is unchanged: general radial functions and integral functionals represent the live profiles; every-`b` conditions use the engine's exterior-half-line Abel inversion, the constant-`GM` anchor and signed squared-velocity reduction with integration-constant domains; the coordinate mass law uses Cartesian component differentiation before radial evaluation. There is no engine profile-family restriction. Item 10's labelled restrictions and forward substitutions retain their stated scope; witness functions restrict only tests. No fold or LEDGER row is bound to a profile, no new binding is introduced, and no relation between flux and `GM` is imposed.

The unchanged supplied-input mapping to lines of `S9b_SHARED_PHYSICS.md` is: dispersion/metric **45**; optical ratio **52**; polarization identification **59**; advection **64**; coordinate mass balance **69, 72–74**; embedding **87**; historical uniform anchor **95**; LAB_HELD **103–104**; far reference speed **116**; delta/nonnegative root **124**; retained grades/counting **128, 135, 138–145**; reference objects **208–210**; pressure/sound/density **219–223**; responses **244–246**. These remain supplied inputs rather than newly derived bindings. Directive item 9 supplies the density premise and root/zero-power conventions; item 10 supplies the labelled no-loss forward premise.

| Knife | Single source mutation | Witness search |
|---|---|---|
| K1 | `kineticSymbol`: remove advective shift | local Fermat time polynomial |
| K2 | `metric2`: replace radial metric entry by one | local Fermat time polynomial |
| K3 | `dispersionSpeed`: replace local square by `c0^2` | local Fermat time polynomial |
| K4a | `rayDerivative`: freeze differentiated delta profile | pushed optical density |
| K4b | `rayDerivative`: freeze differentiated velocity profile | pushed optical density |
| K4c | `xiSlope`: differentiate constant instead of live xi ansatz | printed xi derivative |
| K5 | `returnDirection`: outgoing direction on return | odd one-form radial component |
| K6 | `radarSource`: one-way ER time instead of round trip | none found; full ray constraints unresolved |
| K7 | `quantifierMode`: fixed `b` instead of every `b` | none found; full ray constraints unresolved |
| K8 | `massDensity`: constant instead of live brane density | coordinate divergence |
| K9a | `localSoundSquared`: asymptotic density instead of local | exact fixed-ratio response |
| K9b | `localBulkDensity`: asymptotic density instead of local | exact power response |
| K10 | `branchGateLocal`: negate propagation gate | none found; full ray constraints unresolved |
| K11 | `vVector`: add azimuthal velocity component | compact local Fermat time polynomial |

Every construction is present. The transcript's `*_SITE` entries give literal replacements, unique occurrence counts and source lines. The validator confirms that each full copy differs only by that replacement (Text export omits the terminal newline).

K11 uses the unchanged mechanical prefix through `COMPACT_K11_BOUNDARY` for **both** source copies, with only the K11 mutation retained. It keeps untruncated local dispersion/Fermat and all twelve retained odd-one-form grades, no additional truncation, no samples/seeds or numerical approximation. Its symbolic sample-domain entry describes coverage only. **10 reached tags** are listed [here](S9b_repair_wolfram_round8_k11_evaluated_tags.txt); **364 not evaluated tags** are listed [here](S9b_repair_wolfram_round8_k11_not_evaluated_tags.txt). Each reached row prints the full baseline, the ordered compact **corrupted/unmutated** difference operands, and the ordered **compact-unmutated/full-baseline** method operands, each pair with its structural mark and coverage. All method operand pairs are structurally identical. The method residual supplies no witness. Unreached pairs print `NOT_EVALUATED` and `not compared`. The full corrupted K11 engine was not run.

## Commands and checks

All CAS invocations use `math -script`, one kernel at a time, `--pool s9b --memory-gib 13 --tasks-max 64`, zero swap, and no deadline. The run helper polls the launched guard until it ends; it does not schedule model polling. Guard stdout and stderr are retained under `_scratch/s9b_wl_repair8/`. The exact command arrays and working directories, guard outcomes and effective limits are in validation and the corresponding `_scratch/s11c/<run>/` directories.

Full harness, repository cwd:

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r8-full-final -- math -script research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Structural checks, repository cwd:

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r8-structural-checks5 -- math -script research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl
```

Repository engine:

```sh
python scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r8-engine-repository -- math -script research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

Standalone engine, cwd `/var/projects/toy_physics-s9b-r8-wl/_scratch/s9b_wl_repair8/standalone`:

```sh
python /var/projects/toy_physics-s9b-r8-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r8-wl/_scratch/s11c/s9b-r8-engine-standalone -- math -script S9b_light_bending_mathematica_audit.wl
```

Mechanical validation:

```sh
python research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_validate.py s9b-r8-full-final
```

The standalone and repository streams are nonempty and byte-identical, and also byte-identical to the harness baseline. The standalone directory contains only the engine copy. The capture was observed growing while its actual kernel PID was alive; observations are in validation. `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints nothing.

## Events and limits of the result

No guard refusal, kill, task-limit event or item 12 stop event occurred. Final CAS runs exit `0` with zero stderr bytes. Full-harness cgroup peak: **2507247616 bytes**; maximum sampled tasks: **32**; sampled swap: **0 bytes**.

Two early structural-check executions exposed author bugs in the checks' definition context and an unnecessary explicit `Global`K` shadow. Both were fixed before the full harness launch; these early outputs are excluded from results and remain in scratch. A scratch assembly script's source-slicing bug was also fixed before launch. None required a different physics construction. A preliminary full demonstration completed without messages; its witness bookkeeping omitted explicit ordered-case and joint auxiliary-solution records. Those author bugs were corrected without changing a physics construction, and a new full demonstration was run from the corrected final source. The preliminary source/captures remain pinned in scratch and its output is excluded from the reported results. The final pinned source is `_scratch/s9b_wl_repair8/harness-pinned-final.wl`.

Open: K6/K7/K10 finite searches found no certified witness. Their absence is reported, as amendment 11 permits. No mathematical equality follows from a structural mark or an unsuccessful search. Review and interpretation belong to the orchestrator. **Commit IDs: none.**

## Literal guard stdout from final-code runs

Full harness:

```json
{
  "exitCode": 0,
  "wallSeconds": 208.90914709097706,
  "unit": "s11c-guard-8e6581892259",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 208.83573423797498,
    "childPid": 758118
  }
}
```

Structural and bookkeeping checks:

```json
{
  "exitCode": 0,
  "wallSeconds": 1.4395652180537581,
  "unit": "s11c-guard-ee47bc6cbdff",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 1.3698718799860217,
    "childPid": 758018
  }
}
```

Repository engine:

```json
{
  "exitCode": 0,
  "wallSeconds": 15.840563025965821,
  "unit": "s11c-guard-952c6450ff52",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.774435869010631,
    "childPid": 756133
  }
}
```

Standalone engine:

```json
{
  "exitCode": 0,
  "wallSeconds": 15.790312759985682,
  "unit": "s11c-guard-ac98b055e07f",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.724726385029498,
    "childPid": 736455
  }
}
```
