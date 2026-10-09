# Measurements — S9b repair build directive amendment 2, review round 0 (generated 2026-10-08 21:32)

Generator: `_scratch/s9b_build/gen/s9b_amend2_r0_lookups.sh` (sha256sum/grep/sed only).

## Reviewed version and verdicts
```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend2_reviewed_r0.md
a8d9cee57f132849c446aeb277f83377fd54a4aceafdc40ea97f409ba0036efa  _scratch/s9b_build/S9b_repair_build_directive_amend2_reviewed_r0.md
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend2_review_baseline_r0.sha256
a8d9cee57f132849c446aeb277f83377fd54a4aceafdc40ea97f409ba0036efa  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend2_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend2_review_r0_grok.txt
_scratch/s9b_build/s9b_amend2_review_r0_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_amend2_review_r0_grok.txt:1:Verdict: NEEDS REVISION
```

## The reviewed K11 paragraph (C1, C2, G1, G2)
```
$ grep -n 'K11 evaluation' -A 11 _scratch/s9b_build/S9b_repair_build_directive_amend2_reviewed_r0.md
164:    **K11 evaluation (amendment 2).** The guard killed K11's copy at 8 GiB and again at 16 GiB, both times after the
165-    same 31 tags. For K11 only, the harness may evaluate the corrupted copy by a compact method of the builder's
166-    choice that completes within item 11's limit.
167-    - The mutation is unchanged and stays in force in the evaluated copy.
168-    - The same method is applied to an unmutated copy of the engine. The difference is taken between those two
169-      evaluations.
170-    - For every tag, the harness prints the full baseline payload, both compact evaluations and their difference.
171-      It names the method, every truncation order, and every sampled value with its seed.
172-    - A tag the method does not reach is printed as not evaluated.
173-    - The engine and its baseline run are unchanged.
174-15. **Handoff.** Each builder works in its own fresh repository, exported from the commit that holds this
175-    directive, with no git history. The other engine's S9b files are absent from it. Builders are Codex
```

## G1: item 14's general rule and item 11's kill rule
```
$ grep -n 'exactly that one mutation\|for \*\*every\*\* tag\|Never answer a kill' _scratch/s9b_build/S9b_repair_build_directive_amend2_reviewed_r0.md
122:      limit needs the user. ⛔ Never answer a kill by narrowing or cheapening the requested object.
138:      below, it runs a copy with exactly that one mutation at the named construction site.
139:    - It prints the baseline payload, the corrupted payload and their difference for **every** tag, including tags
```

## G2: emission may not depend on a payload value (item 5); zero-diff under FORM (pipeline §4); reproducer provenance
```
$ grep -n 'Emission never depends' _scratch/s9b_build/S9b_repair_build_directive_amend2_reviewed_r0.md
58:   test returned. Emission never depends on a payload's value. Outside the branch-existence conditions the spec
```

```
$ grep -n 'Zero-diff under a claimed FORM knife\|mechanically extracted from' docs/development_pipeline.md
169:  mechanically extracted from — or calls — the named production construction; an independently reimplemented
181:  a `PASS` tag is the residual-asserted-zero defect in a new file. Zero-diff under a claimed FORM knife is a
```

## C2/G2: item 17's report list
```
$ grep -n '^17\.' -A 9 _scratch/s9b_build/S9b_repair_build_directive_amend2_reviewed_r0.md
179:17. **Report** (in the final message):
180-    - the tags emitted;
181-    - each deferred item (item 13) and how it was handled;
182-    - each binding and its spec line;
183-    - each family restriction and what it leaves out;
184-    - each knife, and where it acts;
185-    - every `NOT_ESTABLISHED`;
186-    - every guard refusal or kill;
187-    - every stop-and-report event (item 12).
188-
```

