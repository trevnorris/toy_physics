# S15 reading build review, round 4: lookups

Mechanical lookups only (existence, verbatim retrieval, literal-match counts). Run from `/var/projects/toy_physics`.

## L1. The two leg reports

```text
$ sha256sum /var/projects/_review_scratch/s15-reading-r4-claude/REVIEW.md /var/projects/_review_scratch/s15-reading-r4-grok/REVIEW.md
b4cb51795008d8369f2137c544343bd5fee49d8b1f76ed1f950c4f082362f87b  /var/projects/_review_scratch/s15-reading-r4-claude/REVIEW.md
003b20bf8742a0c2bee4c63f16bd0a35be535345d19f7b3a644fec14a27954f4  /var/projects/_review_scratch/s15-reading-r4-grok/REVIEW.md
```

## L2. The reviewed baseline (hashes recorded before the first leg launched), and the review prompt both legs read

```text
$ sha256sum -c research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r4_baseline.sha256
/var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py: OK
/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl: OK
/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure_run.py: OK
/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_separation_departure.out: OK
/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out: OK
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md: OK
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S15_separation_departure_attempt_history.md: OK
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S15_separation_departure_attempt_history.json: OK
```

```text
$ sha256sum research/pde_ledger_v3/directives/S15_reading_build_review_prompt_r4.md
0b2dc58d821e675e7f4a23d83b0ceca710b9cad3d3c4a1f72f3ab73e34238692  research/pde_ledger_v3/directives/S15_reading_build_review_prompt_r4.md
```

## L3. Preserved copies of the reviewed transcripts

```text
$ sha256sum research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r4.out research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r4.out
10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6  research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r4.out
5779d92eb7820b11fa51068526145c3f61e04b69d2f10f86d6b1af923531abce  research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r4.out
```

## L4. The notes' text

```text
$ grep -n 'Every condition printed now that is true stays printed' research/pde_ledger_v3/directives/resume_s15_reading_repair_4.md
26:- no object of items 1, 2, 3, 5 or 6, and no item-4 object where `d` is absent, changes by value. Every condition printed now that is true stays printed.
```

```text
$ sed -n '14p;21p' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
- These six objects are not attempted again in this build. Each carries a not-obtained label, set by this decision, from the attempts already made on it.
- In each engine's final run, each of the six objects is printed with a not-obtained label that says it was set by this decision and that the object was not attempted in that run. This is a third kind of label. For these six objects it replaces the rules that a label is one of two kinds, the 6-hour stop or running out of memory (`resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_reading_repair_4h.md`). For every other object those rules stand, and no other object carries this label.
```

## L5. The runner: the decision label text, and the item-4 kernel run for each object

```text
$ sed -n '31,37p;206,211p' research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
DECISION_FAILURES = {
    f'region/component_{i}/{kind}_labels_{{drain1, drain2}}':
        'User decision 2026-11-06: complete item-4 conditions were not returned by the recorded prior attempts. '
        'This not-obtained label is set by that decision; the exponent-sensitive attempt was not made in this final run. '
        'The retained sufficient conditions are printed above. The run record lists the attempt history and the larger-memory follow-up.'
    for i in range(1,4) for kind in ('momentum','rate')
}
        if name in completed or name in args.exclude_object or (args.only_object and name not in args.only_object):continue
        failure = DECISION_FAILURES.get(name, '')
        if failure:
            state.setdefault('failed_objects', {})[name] = {'kind': 'user_decision_2026_11_06', 'label': failure}
            save()
        run_kernel('item4', name, failure)
```

## L6. The engine: conditions printed before the label check

```text
$ grep -n 'solo\[prefix<>"/condition_"\|solo\[prefix<>"/Q_under_condition_"\|attemptFailure=Environment\|emit\[prefix<>"/label_kind"\|nativeAnswer=S15NativeItem4' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-160
1364:   solo[prefix<>"/condition_"<>ToString[conditionIndex],If[valid,(branch && (rout[t]>rin[t]>0/.rules))/.inverse,Missing["NotObtained"]]];
1365:   solo[prefix<>"/Q_under_condition_"<>ToString[conditionIndex],If[valid,restricted/.inverse,Missing["NotObtained"]]],
1370: attemptFailure=Environment["S15_WL_ATTEMPT_FAILURE"];
1373:  emit[prefix<>"/label_kind",If[StringStartsQ[attemptFailure,"User decision 2026-11-06"],"user_decision_2026_11_06","resource_failure"]];
1377: nativeAnswer=S15NativeItem4[
1386:  solo[prefix<>"/condition_"<>ToString[conditionIndex],
1388:  solo[prefix<>"/Q_under_condition_"<>ToString[conditionIndex],
```

## L7. The six item-4 kernels in the final run, from the Wolfram run manifest

```text
$ python3 -c 'import json; j=json.load(open("research/pde_ledger_v3/directives/_measurements/S15_separation_departure_wolfram_run_manifest.json")); [print(k["object"], k["started"], k["ended"], k["returncode"], k["source_sha256"][:12], bool(k.get("attempt_failure_label"))) for k in j["execution"]["kernels"] if k.get("mode")=="item4" and k["object"].startswith("region/") and "{drain1, drain2}" in k["object"]]'
region/component_1/momentum_labels_{drain1, drain2} 2026-11-06T15:43:19.037823+00:00 2026-11-06T15:43:36.071677+00:00 0 366298d82e62 True
region/component_1/rate_labels_{drain1, drain2} 2026-11-06T15:43:36.123258+00:00 2026-11-06T15:44:47.160337+00:00 0 366298d82e62 True
region/component_2/momentum_labels_{drain1, drain2} 2026-11-06T15:44:47.161325+00:00 2026-11-06T15:45:02.164445+00:00 0 366298d82e62 True
region/component_2/rate_labels_{drain1, drain2} 2026-11-06T15:45:02.189738+00:00 2026-11-06T15:46:10.199189+00:00 0 366298d82e62 True
region/component_3/momentum_labels_{drain1, drain2} 2026-11-06T15:46:10.200403+00:00 2026-11-06T15:46:25.205326+00:00 0 366298d82e62 True
region/component_3/rate_labels_{drain1, drain2} 2026-11-06T15:46:25.206491+00:00 2026-11-06T15:47:33.214429+00:00 0 366298d82e62 True
```

## L8. The final Wolfram transcript: condition, label_kind and status lines of the six

```text
$ awk 'NR>=26890 && NR<=29360 && /\{drain1, drain2\}\/item4\/(condition_[0-9]+|label_kind|status):/ {s=$0; sub(/: .*/,"",s); sub(/^.*reading\/region\//,"",s); sub(/condition_[0-9]+$/,"condition_N",s); print s} NR>29360{exit}' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | uniq -c
      4 component_1/momentum_labels_{drain1, drain2}/item4/condition_N
      1 component_1/momentum_labels_{drain1, drain2}/item4/label_kind
      1 component_1/momentum_labels_{drain1, drain2}/item4/status
     16 component_1/rate_labels_{drain1, drain2}/item4/condition_N
      1 component_1/rate_labels_{drain1, drain2}/item4/label_kind
      1 component_1/rate_labels_{drain1, drain2}/item4/status
      4 component_2/momentum_labels_{drain1, drain2}/item4/condition_N
      1 component_2/momentum_labels_{drain1, drain2}/item4/label_kind
      1 component_2/momentum_labels_{drain1, drain2}/item4/status
     16 component_2/rate_labels_{drain1, drain2}/item4/condition_N
      1 component_2/rate_labels_{drain1, drain2}/item4/label_kind
      1 component_2/rate_labels_{drain1, drain2}/item4/status
      4 component_3/momentum_labels_{drain1, drain2}/item4/condition_N
      1 component_3/momentum_labels_{drain1, drain2}/item4/label_kind
      1 component_3/momentum_labels_{drain1, drain2}/item4/status
     16 component_3/rate_labels_{drain1, drain2}/item4/condition_N
      1 component_3/rate_labels_{drain1, drain2}/item4/label_kind
      1 component_3/rate_labels_{drain1, drain2}/item4/status
```

## L9. The run record's sentence on the six

```text
$ grep -n 'were not attempted in these final segments' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md | cut -c1-260
55:The six `user_decision_2026_11_06` objects were not attempted in these final segments; each retained its previously printed sufficient conditions and corresponding Q objects. All 42 nondeferred separation-dependent item-4 traversals returned their computed 
```
