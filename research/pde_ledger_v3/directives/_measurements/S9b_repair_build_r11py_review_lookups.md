# Measurements: S9b build review, round 11, SymPy (generated 2026-10-19 01:43)

Generator: `_scratch/s9b_build/gen/s9b_build_r11py_lookups.sh` (sha256sum/grep/sed/head/wc/ls only). The legs'
evidence is hard-linked at `_scratch/s9b_build/s9b_repair_build_review_r11py_{claude,grok}_evidence/`.

## The reviewed files, the prompt, and the copy into this repository
```
$ (cd /var/projects/toy_physics-s9b-r11-py && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r11py.sha256 | grep -c ': OK$')
8
```

```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r11py.sha256 | grep -c ': OK$')
8
```

```
$ wc -l < /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r11py.sha256
8
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_prompt_r11py.md | cut -c1-64
6f45c0f9d6a58da6e9c7fb6e74091023751ff4f2d877af8d3b844de952fd95a3
```

```
$ git -C /var/projects/toy_physics diff --quiet 86097a44 -- /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_exports.py && echo engine-and-export-unchanged-since-86097a44
engine-and-export-unchanged-since-86097a44
```

## Verdicts
```
$ head -1 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude.md
**Verdict: NEEDS REVISION, for the records only.** The K10 repair, the other knives, the printed output and the inventory all hold up under my own computations and ablations. The single finding is an incomplete run record (item 8). Fixing it is a records edit and needs no rerun.
```

```
$ grep -o 'Verdict:\*\* [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_grok.txt
Verdict:** CLEAR
```

## K10: mutation unchanged, witness reproduced by both legs
```
$ sed -n '46p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K10': ('observable_gate = real_domain', 'observable_gate = sp.Not(real_domain)'),
```

```
$ sed -n '696p;770,771p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | cut -c1-200
    real_domain = everywhere(r,existence.subs(branch_map),traversed)
    observable_gate = real_domain  # K10: existence only
    branch_domain=sp.And(observable_gate,travel_domain)
```

```
$ grep -h 'emitted_objects_at_inputs' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/parse_k10_witness.stdout | head -2
emitted_objects_at_inputs :: (False, True)
```

```
$ grep -h 'FREE_SYMBOLS\|FUNCTIONS' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/k10_witness_step1.stdout | head -4
FREE_SYMBOLS base ['Z_E', 'Z_R', 'b', 'c_0']
FUNCTIONS base ['mu_perp', 'rho_br_profile', 's9b_V', 's9b_xi']
FREE_SYMBOLS corr ['Z_E', 'Z_R', 'b', 'c_0']
FUNCTIONS corr ['mu_perp', 'rho_br_profile', 's9b_V', 's9b_xi']
```

```
$ grep -h 'SET ConditionSet\|PART_VALUE' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/k10_witness_step2.stdout | head -6
  SET ConditionSet -> EmptySet | empty: True
SET ConditionSet -> {(12, 12)} | empty: False
  SET ConditionSet -> EmptySet | empty: True
SET ConditionSet -> {(12, 12)} | empty: False
  SET ConditionSet -> EmptySet | empty: True
SET ConditionSet -> {(12, 12)} | empty: False
```

```
$ grep -h 'REAL_ROOTS\|ACCEPTED\|INDEPENDENT_\|MATCHES_PRINTED' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_grok_evidence/k10_independent_eval.stdout
RAY0 REAL_ROOTS s9b_momentum_FLYBY {12}
RAY0 REAL_ROOTS s9b_turn_FLYBY {-12, 12}
RAY0 ACCEPTED [(12, 12)]
RAY1 REAL_ROOTS s9b_turn_ER {-s9b_momentum_ER, s9b_momentum_ER}
RAY1 REAL_ROOTS s9b_momentum_ER {12}
RAY1 ACCEPTED [(12, 12)]
RAY2 REAL_ROOTS s9b_turn_RE {-s9b_momentum_RE, s9b_momentum_RE}
RAY2 REAL_ROOTS s9b_momentum_RE {12}
RAY2 ACCEPTED [(12, 12)]
INDEPENDENT_BASELINE_DOMAIN False
INDEPENDENT_K10_DOMAIN True
MATCHES_PRINTED True
```

## Form ablations
```
$ grep -h '^[<>].*status\|^[<>].*emitted_objects_at_inputs' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/ablation_diff_mut_M2_single_root.txt | head -4
< Tuple(Str('status'), Str('witness')), 
> Tuple(Str('status'), Str('undecided')), 
< Tuple(Str('emitted_objects_at_inputs'), Tuple(false, true)), 
> Tuple(Str('emitted_objects_at_inputs'), Tuple(true, true)), 
```

```
$ grep -h '^[<>].*status' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/ablation_diff_mut_M3_integrals_unevaluated.txt | head -2
< Tuple(Str('status'), Str('witness')), 
> Tuple(Str('status'), Str('undecided')), 
```

```
$ grep -h 'tuple_solver_retains_all_roots\|ValueError' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_grok_evidence/form/form_ablation.diff | head -4
-PY_LOCAL_S9B_ABLATION_SOFTWARE_FIXTURE_tuple_solver_retains_all_roots: Tuple(Tuple(Str('evaluation'), FiniteSet(Tuple(Integer(-1), Integer(-1)), Tuple(Integer(1), Integer(1)))), Tuple(Str('trace'), Tuple(Tuple(Tuple(Str('binder'), Tuple(Symbol('self_check_x', real=True), Symbol('self_check_y', real=True))), Tuple(Str('input_set'), ConditionSet(Tuple(Symbol('self_check_x', real=True), Symbol('self
-PY_LOCAL_S9B_ABLATION_SOFTWARE_RESULT_tuple_solver_retains_all_roots: true
-PY_LOCAL_S9B_ABLATION_SOFTWARE_CHECKS: Tuple(Tuple(Str('syntax_whitespace'), true), Tuple(Str('argument_order'), true), Tuple(Str('no_algebra'), true), Tuple(Str('float_decimal_tail'), true), Tuple(Str('float_precision'), true), Tuple(Str('python_float_token'), true), Tuple(Str('binder_names_literal'), true), Tuple(Str('density_labels_are_not_conditions'), true), Tuple(Str('derivative_before_poin
+PY_LOCAL_S9B_ABLATION_SOFTWARE_FIXTURE_tuple_solver_retains_all_roots: Tuple(Tuple(Str('evaluation'), FiniteSet(Tuple(Integer(-1), Integer(-1)))), Tuple(Str('trace'), Tuple(Tuple(Tuple(Str('binder'), Tuple(Symbol('self_check_x', real=True), Symbol('self_check_y', real=True))), Tuple(Str('input_set'), ConditionSet(Tuple(Symbol('self_check_x', real=True), Symbol('self_check_y', real=True)), And(Equ
```

## Other knives against round 10
```
$ grep -h 'COMMON\|IDENTICAL\|DIFFERENT' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/compare_r10_r11.stdout | head -4
COMMON 5462 IDENTICAL 5433 DIFFERENT 29
DIFFERENT PY_LOCAL_S9B_ABLATION_SOFTWARE_FIXTURE_finite_attempt_can_certify_mathematical_value r10 (43, 7124) r11 (51, 6713)
DIFFERENT PY_LOCAL_S9B_ABLATION_SOFTWARE_FIXTURE_finite_attempt_can_certify_domain r10 (45, 7330) r11 (53, 6908)
DIFFERENT PY_LOCAL_S9B_ABLATION_SOFTWARE_FIXTURE_false_attached_condition_blocks_witness r10 (47, 6166) r11 (59, 5755)
```

```
$ grep -c 'equal_after_point_sort True' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/compare_modulo_point_order.stdout
16
```

```
$ grep -c 'equal_after_point_sort False' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/compare_modulo_point_order.stdout
0
```

## Claude finding 1: three admitted runs not in the record
```
$ grep -o '"guard_runs"' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_run_measurements.json | wc -l
1
```

```
$ grep -o 's9b-round11-[a-z0-9-]*' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_run_measurements.json | sort -u
s9b-round11-checks2-1
s9b-round11-full-1
s9b-round11-full-2
s9b-round11-full-3
s9b-round11-probe-1
s9b-round11-probe1
s9b-round11-probe-10
s9b-round11-probe-11
s9b-round11-probe-12
s9b-round11-probe-13
s9b-round11-probe-14
s9b-round11-probe-15
s9b-round11-probe-16
s9b-round11-probe-17
s9b-round11-probe-18
s9b-round11-probe-2
s9b-round11-probe-3
s9b-round11-probe-4
s9b-round11-probe-5
s9b-round11-probe-6
s9b-round11-probe-7
s9b-round11-probe-8
s9b-round11-probe-9
```

```
$ grep -c 's9b-round11-probe-19\|s9b-round11-checks1-1\|s9b-round11-gate1-1' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_run_measurements.json /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_run_measurements.json:0
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md:0
```

```
$ grep -n 'admission refusals' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md | cut -c1-300
19:**Events/open questions:** **21** memory-contention admission refusals, all recorded and retried at 18 GiB as authorized; no other refusal, guard kill or item 12 event. `Phi` versus `GM` and unsupplied material relations remain open. No production run or review performed.
```

```
$ head -14 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/run_records.stdout
ADMISSION_ROWS 25 KEYS ['attempt', 'capture', 'command', 'exit_code', 'finished_utc', 'log', 'memory_contention_retry', 'name']
RUN probe exit_code 0 log _scratch/s11c/s9b-round11-probe-19 outcome {'exit': 0, 'wall': 44.6}
  command ['python', '-B', '_scratch/s9b-repair-round11/probe.py']
  named_in_measurements False  named_in_report False
RUN checks1 exit_code 0 log _scratch/s11c/s9b-round11-checks1-1 outcome {'exit': 0, 'wall': 1.1}
  command ['python', '-B', 'research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py', '--scratch', '_scratch/s9b-repair-round11/checks1', '--checks-only']
  named_in_measurements False  named_in_report False
RUN gate1 exit_code 0 log _scratch/s11c/s9b-round11-gate1-1 outcome {'exit': 0, 'wall': 108.5}
  command ['python', '-B', '_scratch/s9b-repair-round11/probe_gate.py']
  named_in_measurements False  named_in_report False
RUN checks2 exit_code 0 log _scratch/s11c/s9b-round11-checks2-1 outcome {'exit': 0, 'wall': 1.1}
  command ['python', '-B', 'research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py', '--scratch', '_scratch/s9b-repair-round11/checks2', '--checks-only']
  named_in_measurements True  named_in_report True
RUN full exit_code 0 log _scratch/s11c/s9b-round11-full-3 outcome {'exit': 0, 'wall': 20247.4}
```

```
$ sed -n '423,433p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_repair_build_directive.md
17. **Report** (in the final message):
    - the tags emitted;
    - each deferred item (item 13) and how it was handled;
    - each binding and its spec line;
    - each family restriction and what it leaves out;
    - each knife, and where it acts;
    - every `NOT_ESTABLISHED`;
    - every guard refusal or kill;
    - every stop-and-report event (item 12);
    - for K11, the compact method, its coverage, and every tag printed as not evaluated;
    - every shared-vocabulary name (item 4) the engine does not emit, and why.
```

```
$ grep -n 'superseded' /var/projects/toy_physics/_scratch/s9b_build/repair_build/brief_py_r11.md
20:4. Rewrite the report as item 17 asks, so that nothing in it describes a superseded run, then stop. Commit nothing.
```

## The final full run's memory
```
$ cat /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11py_claude_evidence/out/measurements_memory.stdout | head -6
full {'memory.current': 19326566400, 'memory.peak': 19327352832, 'memory.swap.current': 0} {'low': 0, 'high': 0, 'max': 14250, 'oom': 0, 'oom_kill': 0, 'oom_group_kill': 0, 'sock_throttled': 0} 19327352832
checks2 {'memory.current': 136175616, 'memory.peak': 218890240, 'memory.swap.current': 0} {'low': 0, 'high': 0, 'max': 0, 'oom': 0, 'oom_kill': 0, 'oom_group_kill': 0, 'sock_throttled': 0} 19327352832
```

