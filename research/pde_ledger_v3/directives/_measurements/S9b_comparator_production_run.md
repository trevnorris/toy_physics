# Measurements — S9b comparator production output (generated 2026-10-15 23:36)

Generator: `_scratch/s9b_build/gen/s9b_comparator_production_record.sh` (sha256sum/cat/cmp/grep/head only; this file is written by it).
The comparator is the fourth repair, accepted in build review round 4 (`S9b_comparator_build_r4_review_disposition.md`).
The production output is the reviewed run `repair-04/full-02` itself, not a rerun: the domain-stratum inventory's order
differs between runs (disposition, observation c), so a rerun would not be the bytes both legs reviewed. The run was
made in the comparator author's export `/var/projects/toy_physics-s9b-comparator-r3` through the repository guard; its
inputs are the committed engine transcripts (below). Filing: `_scratch/s9b_build/gen/s9b_comparator_accept_file.sh` (cp only).

## The run
```
$ python3 -c "import json;d=json.load(open('_scratch/s9b_build/comparator_output/repair-04/full-02/guard/invocation.json'));print(d['command']);print(d['cwd'])"
['python', 'research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py', '--output', '_scratch/s11c/s9b-comparator-r3/full-02/comparison.jsonl', '--accounting', '_scratch/s11c/s9b-comparator-r3/full-02/accounting.jsonl']
/var/projects/toy_physics-s9b-comparator-r3
```

```
$ cat _scratch/s9b_build/comparator_output/repair-04/full-02/guard/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 11116.410525696992,
  "unit": "s11c-guard-c5425124acaf",
  "stderrBytes": 1732566,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 11116.238143314957,
    "childPid": 435791
  }
}
```

## Source run, and the committed comparator
```
$ sha256sum _scratch/s9b_build/comparator_output/repair-04/full-02/source/S9b_cross_engine_comparator.py _scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
326d204bcf76f81d40c141803106e42aa8bbb71da91414aef07a3cac710bb1f7  _scratch/s9b_build/comparator_output/repair-04/full-02/source/S9b_cross_engine_comparator.py
326d204bcf76f81d40c141803106e42aa8bbb71da91414aef07a3cac710bb1f7  _scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py
326d204bcf76f81d40c141803106e42aa8bbb71da91414aef07a3cac710bb1f7  research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
```

```
$ sha256sum _scratch/s9b_build/comparator_output/repair-04/full-02/source/test_S9b_cross_engine_comparator.py _scratch/s9b_build/comparator_r4/test_S9b_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py
efe7249bc18698fedb9d365aa1816a17014ab5967737af680744a5e31c989ad8  _scratch/s9b_build/comparator_output/repair-04/full-02/source/test_S9b_cross_engine_comparator.py
efe7249bc18698fedb9d365aa1816a17014ab5967737af680744a5e31c989ad8  _scratch/s9b_build/comparator_r4/test_S9b_cross_engine_comparator.py
efe7249bc18698fedb9d365aa1816a17014ab5967737af680744a5e31c989ad8  research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py
```

```
$ for f in S9b_comparator_build_report.md S9b_comparator_build_report_function_name_table.md S9b_comparator_build_report_coverage_reasons.md; do cmp _scratch/s9b_build/comparator_r4/$f research/pde_ledger_v3/directives/_measurements/$f && echo "$f IDENTICAL_TO_REVIEWED"; done
S9b_comparator_build_report.md IDENTICAL_TO_REVIEWED
S9b_comparator_build_report_function_name_table.md IDENTICAL_TO_REVIEWED
S9b_comparator_build_report_coverage_reasons.md IDENTICAL_TO_REVIEWED
```

## Inputs
```
$ head -1 _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl | grep -o '"input_sha256":{[^}]*}'
"input_sha256":{"py":"b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142","wl":"6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb"}
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142  research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb  research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
```

## Outputs filed, and their identity with the reviewed run
```
$ sha256sum _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out
ff62c9a3e55f19c1e18d380b92c3b8fd350158f63fe35d406f6a01c1295d9d38  _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl
ff62c9a3e55f19c1e18d380b92c3b8fd350158f63fe35d406f6a01c1295d9d38  research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out
```

```
$ sha256sum _scratch/s9b_build/comparator_output/repair-04/full-02/accounting.jsonl research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator_accounting.jsonl
4a0c4d098e24f1fb4be66d426f2927e92ff8464df4a114232fb0b652180e7f37  _scratch/s9b_build/comparator_output/repair-04/full-02/accounting.jsonl
4a0c4d098e24f1fb4be66d426f2927e92ff8464df4a114232fb0b652180e7f37  research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator_accounting.jsonl
```

```
$ cmp _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out && echo COMPARISON_IDENTICAL_TO_REVIEWED_RUN
COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```

```
$ cmp _scratch/s9b_build/comparator_output/repair-04/full-02/accounting.jsonl research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator_accounting.jsonl && echo ACCOUNTING_IDENTICAL_TO_REVIEWED_RUN
ACCOUNTING_IDENTICAL_TO_REVIEWED_RUN
```

```
$ grep -c '"S9b_cross_engine_comparator.py": "326d204b' _scratch/s9b_build/comparator_output/repair-04/full-02/verification.json
1
```

