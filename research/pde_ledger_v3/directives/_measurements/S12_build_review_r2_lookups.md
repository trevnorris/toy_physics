# S12 build review r2: mechanical lookups

Generator: `_scratch/s12/gen/build_r2_lookups.sh`. Verbatim retrieval, hashes, literal-match counts and existence checks only.

## D1. Baseline, prompt and working tree
```
$ sha256sum -c _scratch/s12/build_review/s12_build_baseline_r2.sha256
research/pde_ledger_v3/scripts/S12_order_source.py: OK
research/pde_ledger_v3/mathematica/S12_order_source.wl: OK
research/pde_ledger_v3/scripts/out/S12_order_source.out: OK
research/pde_ledger_v3/mathematica/out/S12_order_source.out: OK
research/pde_ledger_v3/directives/_measurements/S12_build_run_record.md: OK
$ sha256sum _scratch/s12/build_review/s12_build_review_prompt_r2.md
4a8b9ac830cfaf63f12d42625e39d8bb4cc8b3e937570cefb01549d4d6fe1d9f  _scratch/s12/build_review/s12_build_review_prompt_r2.md
$ git status --porcelain | diff _scratch/s12/build_review/status_before_r2.txt -  (empty = unchanged)
6a7,8
> ?? research/pde_ledger_v3/directives/_measurements/S12_build_review_r2_disposition.md
> ?? research/pde_ledger_v3/directives/_measurements/S12_build_review_r2_lookups.md
diff exit: 1
$ cat /var/projects/_review_scratch/s12-build-r2-grok/DONE
EXIT=0
```

## D2. Verdict lines
```
$ grep -n '^## Verdict' /var/projects/_review_scratch/s12-build-r2-claude/REVIEW.md
3:## Verdict: CLEAR
$ grep -n -o '\*\*Verdict: [A-Z]*\.\*\*' /var/projects/_review_scratch/s12-build-r2-grok/grok_stdout.txt
1:**Verdict: CLEAR.**
$ grep -n -A2 '^## Findings' /var/projects/_review_scratch/s12-build-r2-grok/grok_stdout.txt
3:## Findings
4-
5-None.
```

## D3. The round-1 face-concentrated flux: literal-match counts in the final engines and transcripts
```
grep -c 'supported_flux' research/pde_ledger_v3/scripts/S12_order_source.py: 0
grep -c 'supported_flux' research/pde_ledger_v3/mathematica/S12_order_source.wl: 0
grep -c 'supported_flux' research/pde_ledger_v3/scripts/out/S12_order_source.out: 0
grep -c 'supported_flux' research/pde_ledger_v3/mathematica/out/S12_order_source.out: 0
grep -c 'face_j_heads' research/pde_ledger_v3/scripts/S12_order_source.py: 0
grep -c 'face_j_heads' research/pde_ledger_v3/mathematica/S12_order_source.wl: 0
grep -c 'face_j_heads' research/pde_ledger_v3/scripts/out/S12_order_source.out: 0
grep -c 'face_j_heads' research/pde_ledger_v3/mathematica/out/S12_order_source.out: 0
grep -c 'faceFlux' research/pde_ledger_v3/scripts/S12_order_source.py: 0
grep -c 'faceFlux' research/pde_ledger_v3/mathematica/S12_order_source.wl: 0
grep -c 'faceFlux' research/pde_ledger_v3/scripts/out/S12_order_source.out: 0
grep -c 'faceFlux' research/pde_ledger_v3/mathematica/out/S12_order_source.out: 0
grep -c 'faceFluxHeads' research/pde_ledger_v3/scripts/S12_order_source.py: 0
grep -c 'faceFluxHeads' research/pde_ledger_v3/mathematica/S12_order_source.wl: 0
grep -c 'faceFluxHeads' research/pde_ledger_v3/scripts/out/S12_order_source.out: 0
grep -c 'faceFluxHeads' research/pde_ledger_v3/mathematica/out/S12_order_source.out: 0
grep -c 'face_supported' research/pde_ledger_v3/scripts/S12_order_source.py: 0
grep -c 'face_supported' research/pde_ledger_v3/mathematica/S12_order_source.wl: 0
grep -c 'face_supported' research/pde_ledger_v3/scripts/out/S12_order_source.out: 0
grep -c 'face_supported' research/pde_ledger_v3/mathematica/out/S12_order_source.out: 0
```

## D4. Run record against the guard records
```
$ grep -n -o -E '[0-9a-f]{64}|log-directory [^ ]+' research/pde_ledger_v3/directives/_measurements/S12_build_run_record.md
13:f116251b1f47fb31d2b67bac6867f653683fa6a76216b09feab72c4a73bcdda9
16:log-directory /var/projects/toy_physics/_scratch/s12/guard/wolfram_revision4_01
21:ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055
25:871495c15867462eb766da203bd338786d11fdb45474c4abc86b2df79e1cf926
28:log-directory /var/projects/toy_physics/_scratch/s12/guard/sympy_revision4_01
33:988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97
$ python3 -c ... _scratch/s12/guard/wolfram_revision4_01/outcome.json (exitCode, limitsVerified, guardReason)
0 True None
$ cmp _scratch/s12/guard/wolfram_revision4_01/stdout <transcript>
$ python3 -c ... _scratch/s12/guard/sympy_revision4_01/outcome.json (exitCode, limitsVerified, guardReason)
0 True None
$ cmp _scratch/s12/guard/sympy_revision4_01/stdout <transcript>
wolfram stdout == transcript
sympy stdout == transcript
```

## D5. Each leg's own derivation and controls exist
```
present, 8864 bytes: /var/projects/_review_scratch/s12-build-r2-claude/derive/derive_nGammaB.py
present, 14794 bytes: /var/projects/_review_scratch/s12-build-r2-claude/derive/derive_nGammaB.stdout
present, 1021 bytes: /var/projects/_review_scratch/s12-build-r2-claude/derive/weakform_control.py
present, 213 bytes: /var/projects/_review_scratch/s12-build-r2-claude/derive/weakform_control.stdout
present, 19701 bytes: /var/projects/_review_scratch/s12-build-r2-grok/derive_n_gamma_B.py
present, 37539 bytes: /var/projects/_review_scratch/s12-build-r2-grok/derive_n_gamma_B.stdout
present, 8110 bytes: /var/projects/_review_scratch/s12-build-r2-grok/compare_sympy_to_derivation.stdout
present, 697 bytes: /var/projects/_review_scratch/s12-build-r2-grok/compare_item2_and_wolfram.stdout
$ ls /var/projects/_review_scratch/s12-build-r2-claude/ablate/ | head -20
make_ablations.py
run_ablation.sh
S_A1.diff
S_A1_J_no_jump.guard.json
S_A1_J_no_jump.py
S_A1.out
S_A2.diff
S_A2_G_one_offset.guard.json
S_A2_G_one_offset.py
S_A2.out
S_C.diff
S_C_order_density_corrupt.guard.json
S_C_order_density_corrupt.py
S_C.out
term_diff.py
term_diff_S_A1.stdout
W_A1.diff
W_A1.out
W_A1_u_no_jump.guard.json
W_A1_u_no_jump.wl
$ ls /var/projects/_review_scratch/s12-build-r2-grok/diffs/
source_sympy_corrupt.diff
source_sympy_form.diff
source_wolfram_form.diff
sympy_corrupt__conditional_fast_ordering_region_2_source.diff
sympy_corrupt__general_order_region_2_source.diff
sympy_corrupt__material_face_1_trace_coefficient.diff
sympy_corrupt__material_face_2_trace_coefficient.diff
sympy_corrupt__supplied_material_balance.diff
sympy_form__conditional_fast_ordering_face_1_trace_coefficient.diff
sympy_form__conditional_fast_ordering_region_2_source.diff
sympy_form__general_order_face_1_trace_coefficient.diff
sympy_form__general_order_face_2_trace_coefficient.diff
sympy_form__general_order_region_2_source.diff
sympy_form__material_face_1_trace_coefficient.diff
sympy_form__material_face_2_trace_coefficient.diff
wolfram_form__conditional_fast_ordering_face_1_trace_coefficient.diff
wolfram_form__general_order_face_1_trace_coefficient.diff
wolfram_form__general_order_face_2_trace_coefficient.diff
wolfram_form__general_order_region_1_source.diff
wolfram_form__general_order_region_3_source.diff
```
