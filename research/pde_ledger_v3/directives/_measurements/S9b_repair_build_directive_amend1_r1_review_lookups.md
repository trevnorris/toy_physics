# Measurements — S9b repair build directive amendment 1, review round 1 (generated 2026-10-08 19:02)

Generator: `_scratch/s9b_build/gen/s9b_amend1_r1_lookups.sh` (sha256sum/grep only).

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend1_review_baseline_r1.sha256
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r1.md
50f3c6aec1e2d98519322283b3c5d0a116145265d6b750fc973d62307c3228a0  _scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r1.md
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend1_review_r1_codex_final.txt _scratch/s9b_build/s9b_amend1_review_r1_grok.txt
_scratch/s9b_build/s9b_amend1_review_r1_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_amend1_review_r1_grok.txt:1:Verdict: CLEAR
```

```
$ grep -n 'Findings' _scratch/s9b_build/s9b_amend1_review_r1_codex_final.txt _scratch/s9b_build/s9b_amend1_review_r1_grok.txt
_scratch/s9b_build/s9b_amend1_review_r1_codex_final.txt:3:**Findings:** None within the requested scope.
_scratch/s9b_build/s9b_amend1_review_r1_grok.txt:17:**Findings.** None.
```

