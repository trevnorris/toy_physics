# Measurements: S15 reading decision list, review pass 3 (scoped) lookups

Generator: `_scratch/s15/gen/s15_reading_p3_lookups.sh` (GNU grep, sed, cat, cmp, find, git, sha256sum, `LC_ALL=C`).
Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts and hashes. The
legs' scripts were not re-run here; their filed stdout is retrieved. `X/` is `/var/projects/_review_scratch/s15-readdl-p3-codex/`; `G/` is `/var/projects/_review_scratch/s15-readdl-p3-grok/`.

## M0. Baseline, prompt, legs finished, tree unchanged, no new files in the repository root

```
$ sha256sum -c _scratch/s15/reading_dl/reading_dl_baseline_p3.sha256
research/pde_ledger_v3/directives/S15_reading_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_decision_list_review_p2_disposition.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_decision_list_review_p2_lookups.md: OK
```

```
$ sha256sum _scratch/s15/reading_dl/reading_dl_review_prompt_p3.md
be6588865f13425a225e9b88142a2ab76065835603367390f9737db6424a1dac  _scratch/s15/reading_dl/reading_dl_review_prompt_p3.md
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p3-codex/DONE /var/projects/_review_scratch/s15-readdl-p3-grok/DONE; ls -la /var/projects/_review_scratch/s15-readdl-p3-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p3-grok/REVIEW.md | cut -c24-
exit 0
exit 0
 trevnorris 3595 Oct 28 12:28 /var/projects/_review_scratch/s15-readdl-p3-codex/REVIEW.md
 trevnorris 7357 Oct 28 12:43 /var/projects/_review_scratch/s15-readdl-p3-grok/REVIEW.md
```

```
$ head -12 /var/projects/_review_scratch/s15-readdl-p3-codex/codex.log | grep -E '^(model|reasoning effort|sandbox):'
model: gpt-6.1-sol
sandbox: danger-full-access
reasoning effort: xhigh
```

```
$ diff <(git status --porcelain | grep -v -F -e 'S15_reading_decision_list_review_p3_lookups.md' -e 'S15_reading_decision_list_review_p3_disposition.md'; git rev-parse HEAD) _scratch/s15/reading_dl/status_before_reading_dl_p3.txt && echo 'status unchanged apart from the p3 disposition and lookups files'
status unchanged apart from the p3 disposition and lookups files
```

```
$ find . -maxdepth 1 -newer _scratch/s15/reading_dl/reading_dl_review_prompt_p3.md ! -name .git | wc -l
0
```

## M1. The verdicts

```
$ grep -n -F -e '**FINDINGS**' -e '**CLEAR**' /var/projects/_review_scratch/s15-readdl-p3-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p3-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-140
s15-readdl-p3-codex/REVIEW.md:30:**CLEAR**: nothing outstanding within the revision 3 repair changes what is computed or what may be claimed
s15-readdl-p3-grok/REVIEW.md:60:**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
```

## M2. The repaired text (revision 3)

```
$ sed -n '52,55p;73,74p;86,88p;112,113p;124,125p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
- **Where they exist.** S15's domain and pressure conditions do not by themselves give every separation derivative.
  Each departure object carries a condition under which the separation derivative of its object exists, labelled
  where it is printed as a condition under which it exists. Such a condition is not a condition without which the
  derivative fails to exist, and not the complete condition for the object, unless that is shown.
2. `E_Q^geo` and `E_Q^resp`, for general profiles, exact on the domain under the condition each carries, with every
   live quantity kept.
Every object is exact on S15's domain under the conditions it carries: S15's pressure conditions, and for each
departure object a condition under which it exists, as above. No limit or expansion is taken in place of an exact
object. Where an engine does not obtain an object exactly, it prints a label saying so, beside what it did print.
  It reads each departure object where its printed condition holds, and reads nothing about it where that condition
  does not hold.
For each object S15 prints: its departure from the inverse square of the separation, exact, for general profiles,
with every quantity's change along the family live, where the condition it carries holds. On the power-law profile
```

## M3. Codex: probes (guarded run) and integrity check

```
$ grep -E '^(stationary_expansion_Eresp|psecond_at_base_left_right|whole_kink_separation_derivative|moving_remainder_derivative_right_left|C2_not_C3_EOS_derivative_order_3_left_right|C2_not_C3_expansion_derivative_right_left|whole_minus_expansion/separation_derivative/residual|smooth_sum_geo/residual|smooth_sum_resp/residual):' /var/projects/_review_scratch/s15-readdl-p3-codex/pressure_probes.stdout
whole_minus_expansion/separation_derivative/residual: 0
smooth_sum_geo/residual: 0
smooth_sum_resp/residual: 0
psecond_at_base_left_right: (-1, 1)
stationary_expansion_Eresp: 8*pi*Bp*R**3*d/3
whole_kink_separation_derivative: 8*pi*R**3*d/3
psecond_at_base_left_right: (-1, 1)
moving_remainder_derivative_right_left: (0, -8*pi*R**3/3)
C2_not_C3_EOS_derivative_order_3_left_right: (-6, 6)
C2_not_C3_expansion_derivative_right_left: (0, 0)
```

```
$ cmp /var/projects/_review_scratch/s15-readdl-p3-codex/pressure_probes.stdout /var/projects/_review_scratch/s15-readdl-p3-codex/_scratch/s11c/pressure_probes/stdout && echo 'top-level stdout byte-identical to the guarded stdout'; wc -c < /var/projects/_review_scratch/s15-readdl-p3-codex/_scratch/s11c/pressure_probes/stderr
top-level stdout byte-identical to the guarded stdout
0
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p3-codex/final_integrity.stdout
REPOSITORY_STATUS_UNCHANGED: True
BASELINE_MATCH: research/pde_ledger_v3/directives/S15_reading_decision_list.md True
BASELINE_MATCH: research/pde_ledger_v3/directives/_measurements/S15_reading_decision_list_review_p2_disposition.md True
BASELINE_MATCH: research/pde_ledger_v3/directives/_measurements/S15_reading_decision_list_review_p2_lookups.md True
SOURCE_OR_TRANSCRIPT_UNCHANGED: research/pde_ledger_v3/scripts/S15_two_drain_region.py True
SOURCE_OR_TRANSCRIPT_UNCHANGED: research/pde_ledger_v3/mathematica/S15_two_drain_region.wl True
SOURCE_OR_TRANSCRIPT_UNCHANGED: research/pde_ledger_v3/scripts/out/S15_two_drain_region.out True
SOURCE_OR_TRANSCRIPT_UNCHANGED: research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out True
GUARD_COPY_BYTE_IDENTICAL: True
GUARDED_PROBE_EXIT_CODE: 0
GUARDED_PROBE_LIMITS_VERIFIED: True
GUARDED_PROBE_STDERR_BYTES: 0
TOP_LEVEL_PROBE_STDOUT_BYTE_IDENTICAL: True
REPORT_EXISTS: True
OPEN_QUESTIONS: none within the scoped repair review
```

```
$ grep -o -E 'Different sufficient conditions remain permissible engine disagreement[^.]*\.[^.]*\.' /var/projects/_review_scratch/s15-readdl-p3-codex/REVIEW.md
Different sufficient conditions remain permissible engine disagreement. They must travel with the objects being compared; value agreement on shared coverage does not establish equal coverage or nonexistence elsewhere.
```

## M4. Grok: probes and phrase counts

```
$ grep -E '^(kink_(p_density_quotient_gap|composite_quotient_gap|ball_DQ)|expansion_(pprime_gap|pprime_derivative_gap|fixed_base_Eresp_over_volume)|whole_minus_(p_differentiable_gap_at_total|p_differentiable_gap_at_base|response_quotient_gap)|smooth_whole_minus_quotient_gap|path_(geo_quotient_gap|full_quotient_gap)|sphere_mean_(Q|DQ|Egeo)_residual):' /var/projects/_review_scratch/s15-readdl-p3-grok/existence_and_agreement.stdout
sphere_mean_Q_residual: 0
sphere_mean_DQ_residual: 0
sphere_mean_Egeo_residual: 0
kink_p_density_quotient_gap: 2
kink_composite_quotient_gap: 0
kink_ball_DQ: 0
expansion_pprime_gap: 0
expansion_pprime_derivative_gap: 2
expansion_fixed_base_Eresp_over_volume: Bprime
whole_minus_p_differentiable_gap_at_total: 0
whole_minus_p_differentiable_gap_at_base: 0
whole_minus_response_quotient_gap: -2
smooth_whole_minus_quotient_gap: 0
path_geo_quotient_gap: 2
path_full_quotient_gap: 2*Abs(bprime + 1)
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p3-grok/phrase_lookup.stdout
COUNT 'wherever it exists': 0 lines []
COUNT 'wherever they exist': 0 lines []
COUNT 'needs `p_br`': 0 lines []
COUNT 'p_br′': 0 lines []
COUNT 'without which': 1 lines [54]
COUNT 'complete condition': 1 lines [55]
COUNT 'unless that is shown': 1 lines [55]
COUNT 'condition each carries': 1 lines [73]
COUNT 'condition it carries': 1 lines [125]
COUNT 'printed condition': 2 lines [105, 112]
COUNT 'reads nothing': 1 lines [112]
COUNT 'do not by themselves': 1 lines [52]
COUNT 'a condition under which': 3 lines [53, 54, 87]
```

```
$ wc -c < /var/projects/_review_scratch/s15-readdl-p3-grok/existence_and_agreement.stderr; wc -c < /var/projects/_review_scratch/s15-readdl-p3-grok/phrase_lookup.stderr
0
0
```

## M5. The legs' reports, scripts and stdout (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-readdl-p3-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p3-codex/artifact_audit.py /var/projects/_review_scratch/s15-readdl-p3-codex/artifact_audit.stdout /var/projects/_review_scratch/s15-readdl-p3-codex/pressure_probes.py /var/projects/_review_scratch/s15-readdl-p3-codex/pressure_probes.stdout /var/projects/_review_scratch/s15-readdl-p3-codex/final_integrity.py /var/projects/_review_scratch/s15-readdl-p3-codex/final_integrity.stdout /var/projects/_review_scratch/s15-readdl-p3-grok/REVIEW.md /var/projects/_review_scratch/s15-readdl-p3-grok/existence_and_agreement.py /var/projects/_review_scratch/s15-readdl-p3-grok/existence_and_agreement.stdout /var/projects/_review_scratch/s15-readdl-p3-grok/phrase_lookup.py /var/projects/_review_scratch/s15-readdl-p3-grok/phrase_lookup.stdout | sed 's|/var/projects/_review_scratch/||'
5177c00b83d0f6f947ad18a2b9fa9a0bbb836747dc5f6b3d1c648ed57f04088e  s15-readdl-p3-codex/REVIEW.md
e3bfb1579ac43088ef1571236997ee0496b39cd045510f0941821edd5f705c6a  s15-readdl-p3-codex/artifact_audit.py
a1bf350a24ec53fbd575c265f1da96051a2149e410dbbc8c84ec6a3fa80fcd80  s15-readdl-p3-codex/artifact_audit.stdout
054f3fe3d1601e957c7fb0c3f726874382efffbce38ec4603599b21d1895f8cc  s15-readdl-p3-codex/pressure_probes.py
1f2f6eeefce5d349683d674e948e800c4cb6c543801e5765375e90f5ebe6fb40  s15-readdl-p3-codex/pressure_probes.stdout
3e070d09182424ba4e94b414ae20051340870cf16d3afd2811f582aaa9beabc6  s15-readdl-p3-codex/final_integrity.py
de7f12456cc7be6747258ef93f2435413b1119c721f09b1b26a0d5b9eb7481a1  s15-readdl-p3-codex/final_integrity.stdout
c0873eae95e696ffb8b2b50ff9df9e87f37e426dd29cb32fb5289fcd34a40ca6  s15-readdl-p3-grok/REVIEW.md
14e532790bad9044468cebc8e939ac6208d58535bc1cc3a1909e3de9878d240d  s15-readdl-p3-grok/existence_and_agreement.py
bf404ced1acaffd98b48675644c9a507c48125f7e7c7e824f9e5a9bc51355932  s15-readdl-p3-grok/existence_and_agreement.stdout
954014eef468af4069aec11aae4a3569c4d53dc24f9ab5f4b66ac8d9f995a645  s15-readdl-p3-grok/phrase_lookup.py
06dac267e3cef7ff7ebb36b9b16916532bc4e686fd2d58198ff89590887c29bc  s15-readdl-p3-grok/phrase_lookup.stdout
```
