# Measurements: S14a part 1 record review round 1 lookups (generated 2026-10-21 15:59)

Generator: `_scratch/s14a/gen/rec_r1_lookups.sh` (grep/sed/sha256sum/diff/git/cat/ls only; this file is
written by it). Mechanical lookups only; every computed result quoted here is a leg's saved stdout.

## L1. The round: repair, baseline, prompt, tree, reports

```
$ cat /var/projects/toy_physics/_scratch/s14a/rec/record_author_r1.DONE; grep -A1 'tokens used' /var/projects/toy_physics/_scratch/s14a/rec/record_author_r1.log | tail -1
exit 0
331,514
```

```
$ cd /var/projects/toy_physics && git diff --stat 25cb3a86 -- research/pde_ledger_v3 STATUS.md
 research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md   |   6 +-
 .../scripts/S14a_1_record_measurements.py          |  26 ++
 research/pde_ledger_v3/steps/S14a_1_projection.md  |  37 ++-
 .../_measurements/S14a_1_record_measurements.md    | 273 ++++++++++++++++-----
 4 files changed, 270 insertions(+), 72 deletions(-)
```

```
$ cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s14a/rec/s14a1_record_review_baseline_r1.sha256
research/pde_ledger_v3/steps/S14a_1_projection.md: OK
research/pde_ledger_v3/steps/_measurements/S14a_1_record_measurements.md: OK
research/pde_ledger_v3/scripts/S14a_1_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
research/pde_ledger_v3/steps/S12_order_source.md: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_build_runs.md: OK
STATUS.md: OK
```

```
$ cd /var/projects/toy_physics && git status --porcelain | diff - /var/projects/toy_physics/_scratch/s14a/rec/status_before_record_r1.txt && echo 'status unchanged'
6d5
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_record_r1_review_lookups.md
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s14a/rec/s14a1_record_review_prompt_r1.md
98897157501b959b646797861b764d5838eaaad15bcbf92cf7b7e2fdd467a30f  /var/projects/toy_physics/_scratch/s14a/rec/s14a1_record_review_prompt_r1.md
```

```
$ cat /var/projects/_review_scratch/s14a1-record-r1-grok/DONE; grep -nE 'Verdict: ' /var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md /var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md | cut -c1-200
exit 0
/var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md:3:**Verdict: CLEAR.** I found nothing that changes what may be claimed from S14a part 1, or what the record, the
/var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md:3:**Verdict: CLEAR**
```

```
$ grep -nE '^## Findings' -A2 /var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md; grep -nE '^\*\*Findings\.\*\*' /var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md | cut -c1-200
38:## Findings
39-
40-None.
7:**Findings.** None that change what may be claimed from S14a part 1, or what the record, the register, or the status lines hand on.
```

## L2. T1–T4 in the repaired files

```
$ grep -nE '18 quadratic|quadratic-form samples|quadratic samples' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-200
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:103:18 per form. Only the 18 quadratic-form samples can move when a derivative is taken at another
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:127:| Storage | `K=∂t∫shellΣ` | b171 exact residual unresolved; 54 stored samples `0`, 18 per form. Only its 18 quadrat
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:128:| Shell balance and its residual | `K+Φ+J=0`; residual is the printed lhs | b173 lhs and b174 unresolved exact residual
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:135:36 cubic. Only 18 of those quadratic samples belong to its unresolved lhs; the other 18 belong
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:206:| b50/b72/b128/b148 closed differential-operand integrals; b171 storage; b173 shell-balance lhs; b174 residual | **Open 
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:286:18 coupled quadratic-form samples can move when a derivative is taken at another point; the
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:287:18 affine and 18 cubic samples cannot make that check. b173's 36 quadratic samples across both
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:443:| `SUBSTRATE_REQUIREMENTS.md` | S12 acceptance in status and pass heading/provenance; J1 account names the S9b oracle st
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1528:and printed sample evidence: of each unresolved part's 54 stored samples, only the 18 quadratic-form
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1531:b173 has 36 quadratic samples across both sides, of which only 18 belong to its unresolved lhs.
```

```
$ grep -nE 'affine/cubic forms are unchanged|affine and cubic forms are unchanged' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-220
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:289:"The affine and cubic forms are unchanged" describes what its repair left alone, not sensitivity
```

```
$ grep -nF 'The bound `0≤f≤1` enters no computed object' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md | cut -c1-120
59:| U1 interior `0≤f≤1`, permitted jumps and no concentrated input fields; S12 encoding and user "Leave it out for 
```

```
$ grep -cF 'Commit nothing' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md; grep -cF 'material-loss' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md; grep -nF 'channels by which material leaves the slab' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md | cut -c1-160
0
0
10:specialization**, the record reports the two channels by which material leaves the slab, their sum `J`, the signed
```

```
$ grep -nE 'byte-identical|CMP_EXIT|matches the filed' /var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md /var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md | cut -c1-200
/var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md:85:- **Measurements (item 4).** The generator does retrieval only. My rerun into scratch is byte-identical to the filed
/var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md:67:**Measurements.** Commands are retrievals. The generator, pointed at scratch, matches the filed file:
/var/projects/_review_scratch/s14a1-record-r1-grok/REVIEW.md:71:CMP_EXIT:0
```

## L3. Notes below the filter (Claude)

```
$ sed -n '/^## Notes below the filter/,/^## Checked and found sound/p' /var/projects/_review_scratch/s14a1-record-r1-claude/REVIEW.md | cut -c1-240
## Notes below the filter (no change to claims or handoffs)

1. **R-S12-02 placement.** The by-target row now owes "printed regional/face and projected/shell compatibility". The
   entry's `requirement` field still names only S12's objects. The S14a objects are named in the new scope bullet:
   "The required objects also include …". The entry and its row say the same thing, unlike S12 r1 F3, where the entry
   contradicted its row. So the handoff is intact. Placing the objects in the `requirement` field would match how the
   S12 merge was written.
2. Record l.59 ends "(M17, T2)" after a sentence about concentrated source parts. T2 does not address that sentence;
   its source is the S12 decision list item 1 ("The rule does not apply to `n Γ_B`…"), which is in M01.
3. Record §6 quotes the earlier population directive ("Run it only on a prior result …"). M13 has that quote only
   through the accepted S12 record (l.405) and the S12 pass, not from the directive itself.
4. **My rerun left a file in the repository.** My rerun imported the generator through `importlib` so that its output
   went to scratch. That write side effect created, or replaced, the ignored bytecode cache
   `research/pde_ledger_v3/scripts/__pycache__/S14a_1_record_measurements.cpython-312.pyc` (mtime 15:43:14). It is not
   tracked (`git status` is unchanged). I did not delete it, because I cannot tell whether a copy existed before. All
   later scripts ran with `-I -B`.

## Checked and found sound
```

```
$ sed -n '59p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md | grep -oE '.{0,200}\(M17, T2\).{0,40}'
object. Permitted jumps and no concentrated input fields define the regional/traced function class of the same objects. The no-concentrated-input restriction does not remove concentrated source parts (M17, T2). |
```

```
$ git -C /var/projects/toy_physics check-ignore -v research/pde_ledger_v3/scripts/__pycache__/S14a_1_record_measurements.cpython-312.pyc; ls -la /var/projects/toy_physics/research/pde_ledger_v3/scripts/__pycache__/S14a_1_record_measurements.cpython-312.pyc
.gitignore:1:__pycache__	research/pde_ledger_v3/scripts/__pycache__/S14a_1_record_measurements.cpython-312.pyc
-rw-rw-r-- 1 trevnorris trevnorris 25575 Oct 21 15:43 /var/projects/toy_physics/research/pde_ledger_v3/scripts/__pycache__/S14a_1_record_measurements.cpython-312.pyc
```

## L4. Status text that describes this work as written for review

```
$ grep -nE 'written for source-first review|written for$' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md /var/projects/toy_physics/STATUS.md | cut -c1-200
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:3:**Status, 2026-10-21:** record and register pass written for source-first review. The constructions
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:8:The S14a part 1 record and register pass (2026-10-21) are written for source-first review;
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1483:### S14a part 1 pass — closed projection and two-channel shell compatibility · 2026-10-21; written for source-first r
/var/projects/toy_physics/STATUS.md:11:[projection record](research/pde_ledger_v3/steps/S14a_1_projection.md) and register pass are written for
```
