# Measurements: S12 decision list, review pass 3, scoped (generated 2026-10-20 05:34)

Generator: `_scratch/s12/gen/dl_p3_lookups.sh`. Literal retrieval only.

## Q1. The pass: prompt, baseline, exits, verdict and property lines
```
$ sha256sum /var/projects/toy_physics/_scratch/s12/dl_review_prompt_p3.md | cut -c1-64
c426ea065fee1055da5bb2b62b17d96689863bccd475b11ae4f35631a64e910a
```

```
$ sha256sum -c /var/projects/toy_physics/_scratch/s12/dl_baseline_p3.sha256 | sed 's#/var/projects/toy_physics/##'
research/pde_ledger_v3/directives/S12_decision_list.md: OK
```

```
$ cat /var/projects/_review_scratch/s12-dl-p3-codex/DONE /var/projects/_review_scratch/s12-dl-p3-grok/DONE
EXIT=0
EXIT=0
```

```
$ grep -n 'CLEAR\|FINDINGS' /var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md /var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt | cut -c1-160
/var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md:16:**CLEAR**
/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt:3:**CLEAR**
/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt:92:**CLEAR**
```

```
$ grep -n 'S3-[1-4][^,]*holds' /var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md /var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt | cut -c1-140
/var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md:7:- **S3-1 (N1): holds.** U1 says, “Bulk material that settles into the slab re-or
/var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md:8:- **S3-2 (N2): holds.** The disposition requires arguments “so that no single ch
/var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md:9:- **S3-3 (N3–N6): holds in the requested physics scope.** The inventoried conten
/var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md:10:- **S3-4: holds.** The exact diff has four hunks containing eight changed blocks.
/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt:46:## S3-1 holds
/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt:61:## S3-2 holds
/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt:73:## S3-3 holds
/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt:88:## S3-4 holds
```

## Q2. The repository is unchanged by the legs
The first command drops this pass's own disposition and lookups, written after both legs exited.
```
$ git -C /var/projects/toy_physics status --porcelain | grep -v 'S12_decision_list_review_p3_' | diff /var/projects/toy_physics/_scratch/s12/status_before_dl_p3.txt - && echo IDENTICAL
IDENTICAL
```

