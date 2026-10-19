# Measurements: scoped round-2 review of the S9b flow-only oracle check (generated 2026-10-19 16:10)

Generator: `_scratch/s9b_oracle/gen/check_r2_lookups.sh` (sha256sum/grep/sed/cut/wc/cat/git status only). The
legs' scratch directories are hard-linked at `_scratch/s9b_oracle/check_review/r2_{claude,grok}_evidence/`.

## Baseline, prompt and verdicts
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_baseline_r2.sha256)
research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check_values.wl: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/identification_map.json: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_visser_form.out: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_values.out: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_run.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_controls.md: OK
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r2.md | cut -c1-64
d25f4f4d640ef651c37e59a326af6e5235e08426e0d53af2b740c2c75b5acf07
```

```
$ grep -n '^- \*\*Z[1-4]\.\*\*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r2.md | cut -c1-120
27:- **Z1.** Whether each residual vanishes on its printed domain can be read from the output, without further algebra.
31:- **Z2.** The repair changes how residuals are printed, not what they are. Each residual's value under its printed
33:- **Z3.** The controls record reports FORM-control sensitivity as movement of values, not inequality of text.
34:- **Z4.** By the user's choice, the builder's Hamilton–Lisle FORM control is dropped. The controls record says that
```

```
$ grep -o 'Verdict: *[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude.md
Verdict: CLEAR
```

```
$ grep -o 'Verdict: [A-Z]*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_grok.txt | head -1
Verdict: CLEAR
```

```
$ grep -n -A2 '^\*\*Findings' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude.md
3:**Findings:** none that pass the filter.
4-
5-All work is in `/var/projects/_review_scratch/s9b-oracle-check-r2-claude/`. Every CAS run went through the pooled guard, one kernel at a time. Guard logs are under `/var/projects/toy_physics/_scratch/s11c/s9b-oracle-check-review-r2/claude_*`. Every run was admitted with exit 0, and none was refused or killed. The working tree was not modified. The baseline sha256 file checks OK on all 11 files.
```

```
$ grep -n -A2 '^## Findings' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_grok.txt
3:## Findings
4-
5-None.
```

## The repaired transcript: residuals and zero tests (literal counts)
```
$ grep -c '^RESIDUAL_ASSUMPTIONS' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
2806
```

```
$ grep -c '^RESIDUAL_ZERO_TEST' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
2806
```

```
$ grep -c '_RESIDUAL: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
253
```

```
$ grep '^RESIDUAL_ZERO_TEST' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -c '"value" -> True|>$'
361
```

```
$ grep '^RESIDUAL_ZERO_TEST' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -c '"value" -> False|>$'
40
```

## Each leg's own census and value comparison
```
$ cat /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/parse_triples.stdout
NEW_RESIDUALS 2806 ANOMALIES [] 0
OLD_RESIDUALS 2806 ANOMALIES [] 0
CLASS_COUNTS {'zero_test_not_established(conditional)': 2405, 'printed_0': 253, 'zero_test_False': 40, 'nonzero_form_zero_test_True': 108}
BY_TAG EXCHANGE_BALANCE_OPERATOR_RESIDUAL {'nonzero_form_zero_test_True': 8, 'printed_0': 10}
BY_TAG EXCHANGE_ORACLE_OPERATOR_SOURCE_RESIDUAL {'zero_test_not_established(conditional)': 8, 'printed_0': 10}
BY_TAG EXCHANGE_RESIDUAL {'zero_test_not_established(conditional)': 8, 'printed_0': 10}
BY_TAG FLOW_SIGNED_RESIDUAL {'zero_test_not_established(conditional)': 1172, 'zero_test_False': 40, 'printed_0': 40, 'nonzero_form_zero_test_True': 100}
BY_TAG FLOW_SQUARED_RESIDUAL {'zero_test_not_established(conditional)': 964, 'printed_0': 180}
BY_TAG OPTICAL_CHARACTERISTIC_RESIDUAL {'zero_test_not_established(conditional)': 253, 'printed_0': 3}
OLD_NEW_FIELD_DIFFERENCES {} first []
CONDITIONAL_SHAPES {'has_==': 2405}
DISTINCT_DOMAIN_IDS 22
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/2 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/4 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/5 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/6 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general 348
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/1 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/2 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/3 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/4 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/5 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/6 107
DOMAIN PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/general 156
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1 107
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/2 107
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3 107
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general 220
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/1 107
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/2 107
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/case/3 107
DOMAIN WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR/general 156
```

```
$ grep -E '^SUMMARY|^ANOMALY_COUNT' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/sample_full3.stdout
SUMMARY_KEY: "{class, has_witness, all_samples_zero, all_samples_nonzero} -> count"
SUMMARY: {"nonzero_form_test_True", True, True, False} -> 108
SUMMARY: {"printed_0", True, True, False} -> 253
SUMMARY: {"test_conditional", True, False, True} -> 2405
SUMMARY: {"test_False", True, False, True} -> 40
ANOMALY_COUNT: 0
```

```
$ wc -c /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/R1_vs_committed.diff /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/R2_vs_committed_visser.diff /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/run019_zero_tests.diff
0 /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/R1_vs_committed.diff
0 /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/R2_vs_committed_visser.diff
0 /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/run019_zero_tests.diff
0 total
```

```
$ cat /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/run019_zero_tests.stdout
rc=0 difflines=0 run019_tests=2806
```

```
$ sed -n '9,26p' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_grok_evidence/01_classify_and_diff.stdout
  "baseline": {
    "nlines": 17081,
    "nrecords": 2806,
    "classes": {
      "does_not_vanish": 2445,
      "printed_zero": 253,
      "nonzero_form_vanishes": 108
    },
    "kinds": {
      "conditional": 2405,
      "true": 361,
      "false": 40
    },
    "printed_zero_not_true_lines": [],
    "bad_count": 0,
    "message_count": 0,
    "lengths": {
      "does_not_vanish": {
```

## Zero tests respond to ablations (one line each)
```
$ grep -c '' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/A12_vs_baseline.diff /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/A3_vs_baseline.diff
/var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/A12_vs_baseline.diff:50
/var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/A3_vs_baseline.diff:2158
```

```
$ grep -n 'REVIEW_ABLATION_DOMAIN_REMOVED' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude_evidence/A12_vs_baseline.diff | head -2
16:> REVIEW_ABLATION_DOMAIN_REMOVED: {114, {(2*py`GM)/py`s9b$velocity$zero$radius == 0}}
41:> REVIEW_ABLATION_DOMAIN_REMOVED: {2058, {(2*GM)/r == 0}}
```

## Z3/Z4 and the two notes
```
$ grep -c '^FORM_VALUE_DELTA' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_values.out
18
```

```
$ grep -c '^FORM_VALUE_DELTA: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_values.out
10
```

```
$ grep -o 'EXCHANGE_[A-Z_]*RESIDUAL' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check_values.wl | sort | uniq -c
      2 EXCHANGE_BALANCE_OPERATOR_RESIDUAL
```

```
$ grep -n 'choices' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py
155: parser.add_argument('--control',choices=['baseline','visser_form'],default='baseline')
```

```
$ grep -c 'hl_spatial\|HL_SPATIAL' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl
0
```

```
$ sed -n '16p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_controls.md
**The Hamilton–Lisle FORM control is dropped.** With exact zero tests it did not complete at 50 GiB in `run_023` and `run_025`. These are historical commands for the pre-removal code, not commands supported by the final driver:
```

```
$ (cd /var/projects/toy_physics && git status --short research/pde_ledger_v3/scripts/out/)
 T research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
 D research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_hl_spatial_form.out
 T research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_visser_form.out
?? research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check_values.out
```

```
$ grep -n -A3 '^### Notes that change no claim' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r2_claude.md
124:### Notes that change no claim
125-- The value checker compares only `EXCHANGE_BALANCE_OPERATOR_RESIDUAL`. The Visser control also moves 8 `EXCHANGE_ORACLE_OPERATOR_SOURCE_RESIDUAL` tests (from `3ρ+2rρ'==0` to `False`); this is not misreported, only not counted.
126-- With the HL control dropped by the user's choice, the optical zero tests have no FORM control, and no record claims one.
127-- That the 2405 conditional tests are correct rests on this output as measured above, not on any structural guarantee from FullSimplify.
```

## Round 1's FORM change to Hamilton–Lisle (1) moved the optical residuals (Grok r1 evidence)
```
$ grep -n 'OPTICAL_CHARACTERISTIC_RESIDUAL' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_grok_evidence/ablate/form_cross/diff.stdout
49:OPTICAL_CHARACTERISTIC_RESIDUAL: base 256 other 256 differ 64 extra 0
```

```
$ grep -n 'FORM: ' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r1_grok.txt | cut -c1-200
39:- FORM: `(dr + β dt)^2` replaced by `(dr + β dt)(dr - β dt)`. Drift moves from `-(Sqrt[2]*Sqrt[gmH/r])` to `0`. The signed residual drops the `gmH` term. Exchange does not move (0 of 18). `ablat
```

