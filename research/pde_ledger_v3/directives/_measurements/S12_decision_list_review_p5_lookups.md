# S12 decision list review pass 5: mechanical lookups

Generator: `_scratch/s12/gen/dl_p5_lookups.sh`. Verbatim retrieval and hashes only.

## R1. Baseline, prompt, working tree, exits
```
$ sha256sum -c _scratch/s12/dl_baseline_p5.sha256
research/pde_ledger_v3/directives/S12_decision_list.md: OK
$ sha256sum _scratch/s12/dl_review_prompt_p5.md
8bc4e5dd49b0499dfcf41a79b38f04a672ff1e7194c269d286f90e19dddc88e6  _scratch/s12/dl_review_prompt_p5.md
$ git status --porcelain | diff _scratch/s12/status_before_dl_p5.txt -  (empty = unchanged)
2a3,4
> ?? research/pde_ledger_v3/directives/_measurements/S12_decision_list_review_p5_disposition.md
> ?? research/pde_ledger_v3/directives/_measurements/S12_decision_list_review_p5_lookups.md
diff exit: 1
$ cat /var/projects/_review_scratch/s12-dl-p5-codex/DONE /var/projects/_review_scratch/s12-dl-p5-grok/DONE
EXIT=0
EXIT=0
```

## R2. Verdicts and property lines
```
$ tail -1 /var/projects/_review_scratch/s12-dl-p5-codex/REVIEW.md
**CLEAR**: nothing outstanding in the scoped revision-4 change changes what is computed or what may be claimed.
$ grep -n -o 'S5-1, S5-2, and S5-3 all hold[^.]*\.' /var/projects/_review_scratch/s12-dl-p5-codex/codex_review.txt | tail -1
7371:S5-1, S5-2, and S5-3 all hold.
$ grep -n '^## S5-[123] holds' /var/projects/_review_scratch/s12-dl-p5-grok/grok_stdout.txt
3:## S5-1 holds
28:## S5-2 holds
34:## S5-3 holds
$ grep -n -o '^\*\*CLEAR\*\*[^.]*\.' /var/projects/_review_scratch/s12-dl-p5-grok/grok_stdout.txt
134:**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
```

## R3. Each leg's own face derivation exists
```
present, 20652 bytes: /var/projects/_review_scratch/s12-dl-p5-grok/evidence/04_face_content.py
present, 3457 bytes: /var/projects/_review_scratch/s12-dl-p5-grok/evidence/04_face_content.stdout
present, 4667 bytes: /var/projects/_review_scratch/s12-dl-p5-codex/05_independent_face_calculus.py
present, 2440 bytes: /var/projects/_review_scratch/s12-dl-p5-codex/05_independent_face_calculus.stdout.txt
```
