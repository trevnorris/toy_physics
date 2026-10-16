# S9b Wolfram repair round 7 measurements

The final harness implements directive amendment 10 with the engine unchanged. No commits, reviews, production run or physics verdict. The user resumed after the initial admission refusal and replaced 16 GiB by **13 GiB**. Every result below comes from the final pinned harness/check source or the unchanged engine.

## Comparison method

The ordering is identity proof first, then symbolic existence or an exact witness, then UNDECIDED after a finite unsuccessful search. Inert integrals of zero normalize to zero before any evaluation. Polynomial/Boolean normalization and supported exact real quantifier elimination supply sufficient identity/change proofs; a failed search or agreeing evaluations supply no equality proof. Powers retain definedness charts. A chart-only failure is not global equality.

Unresolved parts are tested with four fixed smooth rational radial-function configurations. Each assigns one function, including fixed coefficients, to every live radial profile; derivatives come from those functions. These are witness inputs only, never engine substitutions. Candidate scalar points are exact, with exact algebraic FindInstance where supported. Nonfinite candidates and poles cannot certify a value change. The evaluator attempts complete rational integrals with their printed endpoints and exact real formulas; unsupported integrals and other unresolved expressions remain undecided. There is no numerical quadrature, cutoff, tolerance, seed or finite-precision result. Every certified witness prints its point, profiles, both evaluations, common-domain evaluation and zero error bound.

Value comparisons require both attached domains. Attached domains are compared separately; standalone relations and predicates remain values. Bound variables are renamed before substitution, irrelevant ambient assumptions are excluded, sets are compared by membership at an element, and statements with binders are compared only on free variables. Finite witness-search coverage is not a proof of completeness.

Final generic comparison demonstrations, including zero precedence, profile derivatives/integrals, predicates, quantified domains, renamed set variables, set membership, free-parameter statements, nested binders, undefined candidates and agreement without equality: [literal stdout](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_comparison.txt). Generating [script](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl), mechanically loading the final harness definitions. These test comparator behavior, not expected physics values.

## Final commands and literal streams

### comparison-final2

Working directory: `/var/projects/toy_physics-s9b-r7-wl`

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r7-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-comparison-final2 -- math -script /var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 2.2363789030350745,
  "unit": "s11c-guard-57371f76d446",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 2.166152346006129,
    "childPid": 516878
  }
}
```

[Complete literal stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-comparison-final2/stdout); [stderr](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-comparison-final2/stderr).

### full-final

Working directory: `/var/projects/toy_physics-s9b-r7-wl`

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r7-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-full-final -- math -script /var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 3883.652372639044,
  "unit": "s11c-guard-8311e9aaa724",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 3883.5863135259715,
    "childPid": 516932
  }
}
```

[Complete literal stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-full-final/stdout); [stderr](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-full-final/stderr).

### engine-repository

Working directory: `/var/projects/toy_physics-s9b-r7-wl`

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r7-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-engine-repository -- math -script /var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 15.6935867019929,
  "unit": "s11c-guard-3f9349ba979f",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.623570243013091,
    "childPid": 548535
  }
}
```

[Complete literal stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-engine-repository/stdout); [stderr](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-engine-repository/stderr).

### engine-standalone

Working directory: `/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/standalone-final4`

```sh
/usr/bin/python /var/projects/toy_physics-s9b-r7-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 13 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-engine-standalone -- math -script /var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/standalone-final4/S9b_light_bending_mathematica_audit.wl
```

Literal guard stdout:

```json
{
  "exitCode": 0,
  "wallSeconds": 15.793260288017336,
  "unit": "s11c-guard-b85cdc77cffa",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 15.724126130982768,
    "childPid": 549994
  }
}
```

[Complete literal stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-engine-standalone/stdout); [stderr](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final4-engine-standalone/stderr).

## Mechanical validation commands

```sh
python _scratch/s9b_wl_repair7/validate.py
python _scratch/s9b_wl_repair7/source_and_payload_checks.py
python _scratch/s9b_wl_repair7/resources.py
```

Generating scripts and literal stdout: [payload/tag validator](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/validate.py) → [stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/validation.stdout); [mutation/density-parts checker](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/source_and_payload_checks.py) → [stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/payload_checks.stdout); [resource summary](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/resources.py) → [stdout](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/resource_summary.stdout). These perform string, file and recorded-resource checks; they do not perform a physics/CAS calculation.

## Full payload and tag inventory

374 baseline tags = 330 shared + 44 local. No duplicate, missing or extra shared tag; local registry matches. Each of the thirteen full mutated copies, K11 compact comparison, and baseline-against-baseline control emits every baseline tag. Every full printed baseline/corrupted payload is checked against its corresponding capture. Every evaluated row carries separate value and domain differences; unreached K11 rows are explicitly NOT_EVALUATED. Counts: **164** rows containing CHANGED; **1875** rows containing UNDECIDED (these sets may overlap); **1908** rows whose printed pair is not `{0,0}` (including undecided rows). These are harness classifications, not physics conclusions.

[Validation](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_validation.json); [all tags](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_tags.txt); [changed rows](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_changed_rows.txt); [undecided rows](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_undecided_rows.txt); [every conditional NOT_ESTABLISHED tag](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_not_established_tags.txt) (**271**); [all harness rows carrying NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_harness_not_established_tags.txt). No new omitted item 10 object, shared name or engine profile family.

## Item 13 and supplied bindings

The engine is unchanged. General radial profiles and integral functionals remain live. Every-b conditions use the engine's Abel inversion and GM anchor, then signed velocity-squared reductions with their integration-constant domains. Mass balance uses Cartesian divergence on the supplied coordinate measure. No fold/LEDGER binding, new profile family or physics premise is introduced.

Supplied-input/spec-line bindings: dispersion/metric **45**; speed ratio **52**; polarization identification **59**; advection **64**; coordinate mass balance **69, 72–74**; embedding **87**; historical uniform anchor **95**; LAB_HELD **103–104**; reference speed **116**; delta/nonnegative root **124**; retained grades/counting **128, 135, 138–145**; reference objects **208–210**; pressure/sound/density **219–223**; responses **244–246**. Directive item 9 supplies the density premise and root/zero-power conventions. Item 10 supplies the labelled forward no-loss case. Its restrictions remain labelled substitutions, not profile-family selections. Witness configurations constrain only the tests and do not restrict any emitted engine object.

## Knives and K11

K1 acts at kinetic advection; K2 at the metric; K3 at dispersion speed; K4a at delta ray derivatives, K4b at velocity ray derivatives, K4c at xi slope; K5 at return direction; K6 at the radar slope's source; K7 at the every-b quantifier; K8 at mass density; K9a at the local sound response, K9b at local density in the power response; K10 at the branch gate; K11 at azimuthal velocity in advection. [Exact replacement and payload checks](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_payload_checks.json) verify one executable replacement per copy and the density-part object pairs. Wolfram text import/export drops the original terminal newline; the checker records that sole non-executable byte difference explicitly. The compact prefixes are byte-exact against their corresponding source prefixes.

K11 uses the unchanged exact symbolic mechanical prefix at `COMPACT_K11_BOUNDARY` in both copies. Its only source mutation is K11; the mutated full engine is not run. Configuration keeps untruncated local dispersion/Fermat and all twelve retained one-form grades, with symbolic constructor domains, no samples/seeds/precision approximation or additional truncation. **10** tags are reached: [complete list](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_k11_evaluated_tags.txt). Their compact baseline/corrupted payloads, differences, full baseline, method residual and coverage are printed. All method residuals are zero in value and domain. **364** unreached tags explicitly print NOT_EVALUATED: [complete list](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_k11_not_evaluated_tags.txt). Witness search on a compact difference is separate from the compact method's configuration and does not turn an unreached tag into an evaluated one.

## Executable and resource checks

The unchanged engine alone ran both in the repository and in an otherwise empty scratch directory. Both exited zero, streams were nonempty and byte-identical, and the scratch directory contains only the copied `.wl`. The harness baseline is also byte-identical to those streams. Capture growth was observed with a live kernel. `git status --porcelain research/pde_ledger_v3/mathematica/out/` printed nothing. [Executable checks](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_executable_checks.json); [resource evidence](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_resource_summary.json); [input hashes](/var/projects/toy_physics-s9b-r7-wl/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_final_source_hashes.json).

## Events and stopping

One earlier admission refusal in this round: `RuntimeError: unaccounted guarded service; refusing pooled overlap`. The user's resumed instruction states that the other run ended and the shared launcher was installed. Its literal output is retained in the adjacent measurements JSON.

Three incomplete resumed full diagnostic runs were author-terminated for comparator bugs, with sources/logs preserved and null guardReason. The first encountered `FindInstance::exvar` (nonfinite input) and `FindInstance::vlist` (empty variable list). The second exposed profile-name metadata in the solver inputs. Named profile heads are now excluded from scalar algebra, pure-function parameters retain their binding, and only formulas over the declared variables reach FindInstance. The third exposed a singular candidate during exact witness evaluation. Candidate arithmetic is now isolated: poles or divergent integrals reject only that candidate, and message-bearing evaluations are not cached. Other programming messages remain visible. The exact K6 row was checked after correction. None of these changes alters a physics construction. That incomplete diagnostic output is not a result here. The first source is [preserved draft](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/ablation_before_findinstance_fix.wl), with [literal messages](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-full-final/stderr) and [author termination record](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/author_bug_termination.json). The second [source](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/ablation_before_metadata_fix.wl), [messages](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final2-full-final/stderr) and [termination record](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/author_bug_termination2.json) are also retained. Exact diagnostic commands/outcomes are in [first diagnostic run record](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/before_findinstance_fix_final_runs.json) and [second diagnostic run record](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/before_metadata_fix_final_runs.json); these are excluded from final results. The third [run record](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/before_candidate_guard_final_runs.json), [source](/var/projects/toy_physics-s9b-r7-wl/_scratch/s9b_wl_repair7/ablation_before_candidate_guard.wl) and [messages](/var/projects/toy_physics-s9b-r7-wl/_scratch/s11c/s9b-r7-final3-full-final/stderr) are preserved as well.

No resumed guard admission refusal, resource kill, task-limit event or item 12 event. Final runs exited zero with zero stderr bytes, one task kernel at a time, 13 GiB hard caps and zero swap, without deadlines. All required build demonstrations are complete. Interpretation/review belongs to the orchestrator. Commit IDs: none.
