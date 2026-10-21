# Measurements: S14a part 1 record review round 0 lookups (generated 2026-10-21 15:37)

Generator: `_scratch/s14a/gen/rec_r0_lookups.sh` (grep/sed/sha256sum/diff/git/cat only; this file is written
by it). Mechanical lookups only; every computed result quoted here is a leg's saved stdout.

## K1. The round: author, baseline, prompt, tree, reports

```
$ cat /var/projects/toy_physics/_scratch/s14a/rec/record_author.DONE; grep -A1 'tokens used' /var/projects/toy_physics/_scratch/s14a/rec/record_author.log | tail -1
exit 0
239,180
```

```
$ cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s14a/rec/s14a1_record_review_baseline_r0.sha256
research/pde_ledger_v3/steps/S14a_1_projection.md: OK
research/pde_ledger_v3/steps/_measurements/S14a_1_record_measurements.md: OK
research/pde_ledger_v3/scripts/S14a_1_record_measurements.py: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
research/pde_ledger_v3/steps/S12_order_source.md: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_build_runs.md: OK
STATUS.md: OK
```

```
$ cd /var/projects/toy_physics && git status --porcelain | diff - /var/projects/toy_physics/_scratch/s14a/rec/status_before_record_r0.txt && echo 'status unchanged'
6d5
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_record_r0_review_lookups.md
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s14a/rec/s14a1_record_review_prompt_r0.md
eeaeebef5aab3a2a27d5dedcc3727324f54e38a5f0807432cf5b2efccbf3ea13  /var/projects/toy_physics/_scratch/s14a/rec/s14a1_record_review_prompt_r0.md
```

```
$ cat /var/projects/_review_scratch/s14a1-record-r0-grok/DONE; grep -nE 'Verdict: ' /var/projects/_review_scratch/s14a1-record-r0-claude/REVIEW.md /var/projects/_review_scratch/s14a1-record-r0-grok/REVIEW.md | cut -c1-200
exit 0
/var/projects/_review_scratch/s14a1-record-r0-claude/REVIEW.md:3:**Verdict: NEEDS REVISION** (one finding that changes what may be claimed from the sample check, three minor ones).
/var/projects/_review_scratch/s14a1-record-r0-grok/REVIEW.md:3:**Verdict: CLEAR**
```

```
$ grep -nE '^### [0-9]\. ' /var/projects/_review_scratch/s14a1-record-r0-claude/REVIEW.md | cut -c1-200
11:### 1. The carried "quadratic form only" limit does not travel with the seven unresolved exact rows
60:### 2. Minor: the U1 bound `0≤f≤1` is given a printed dependence, but no computed object uses it
82:### 3. Minor: a quotation attributed to the user is not in any retrieved source
100:### 4. Minor: the headline calls the channels "material-loss"
```

```
$ grep -nE 'byte-identical|byte for byte' /var/projects/_review_scratch/s14a1-record-r0-claude/REVIEW.md /var/projects/_review_scratch/s14a1-record-r0-grok/REVIEW.md | cut -c1-220
/var/projects/_review_scratch/s14a1-record-r0-claude/REVIEW.md:147:  - My rerun, with its output sent to scratch, is byte-identical (sha256 `c66be735…`, `regen_compare.stdout`).
/var/projects/_review_scratch/s14a1-record-r0-grok/REVIEW.md:26:- The generator, rerun with its output path pointed at scratch and the repository otherwise unread for writing, matches the filed measurements byte for byte
```

## K2. Claude 1: which samples can see a derivative at another point

```
$ sed -n '100,105p;122,125p;197,199p;276,280p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md | cut -c1-260
Closed integrals of the differential operands are also printed: b50/b72 and spherical b128/b148.
Their exact residuals are **not resolved to literal zero**; each has 54 stored sample residuals `0`.
They stay open exact cross-engine differences, **no owner named** (§5). The independently bound
integration-check operands b60/b61/b80/b81/b138/b139/b156/b157 do print exact residual `0`.
Those distinct comparisons do not close the unresolved integral rows. (Retrieve M07/M08; M02,
comparator round-1 "Carried".)
| Signed order lost | `Lχ=−∫shell Sχ`, including both face parts | b169 exact `0`. Positive denotes a net drain, negative net re-ordering; neither sign is enforced. |
| Outward outer material flux | `Φ(R,t)=4πR²F_r(R,t)` | b170 exact `0`; evaluated at a free outer radius, with no radius limit. |
| Storage | `K=∂t∫shellΣ` | b171 exact residual unresolved; 54 stored samples `0`. |
| Shell balance and its residual | `K+Φ+J=0`; residual is the printed lhs | b173 lhs and b174 unresolved exact residuals, 54 samples `0` each; b173 rhs exact `0`. |
| Difference | Classification and permitted reading | Owner / evidence |
|---|---|---|
| b50/b72/b128/b148 closed differential-operand integrals; b171 storage; b173 shell-balance lhs; b174 residual | **Open exact cross-engine differences.** Each exact expression stays printed; each has 54 stored sample residuals `0`. Finite sample agreement is t
the distinct projected divergence/source-law residual objects stay printed under their own names.
No supplied equation is applied by this record to erase the distinction. The comparator's repaired
nested binders retain variable identity, and its coupled quadratic probes detect derivative-point
changes where the unresolved exact residual cannot resolve them; affine/cubic forms are unchanged.
These are the accepted instrument's reported controls, not new runs here. The comparator has three
```

```
$ grep -nE 'quadratic form only|affine and cubic forms are unchanged' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_build_review_r1_disposition.md | cut -c1-220
44:  partial depends on every other argument. The affine and cubic forms are unchanged (E3).
96:    samples through the quadratic form only (E3).
```

```
$ grep -nE 'additively separable|changes no sample' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_build_review_r0_disposition.md | cut -c1-220
60:  - The test functions are additively separable: each form is a constant plus a sum of one-argument terms. The faces
64:    changes no sample.
```

```
$ sed -n '1,2p' /var/projects/_review_scratch/s14a1-record-r0-claude/testfn.stdout | cut -c1-200
functions per form {'affine': 396, 'quadratic': 396, 'cubic': 396}
functions containing a literal cross-argument product a_i*a_j (i!=j), per form {'quadratic': 396}
```

```
$ grep -nE '^b(50|72|128|148|171|173|174) ' /var/projects/_review_scratch/s14a1-record-r0-claude/samples_check.stdout
3:b50 {'affine:0': 18, 'quadratic:0': 18, 'cubic:0': 18}
4:b72 {'affine:0': 18, 'quadratic:0': 18, 'cubic:0': 18}
5:b128 {'affine:0': 18, 'quadratic:0': 18, 'cubic:0': 18}
6:b148 {'affine:0': 18, 'quadratic:0': 18, 'cubic:0': 18}
7:b171 {'affine:0': 18, 'quadratic:0': 18, 'cubic:0': 18}
8:b173 {'affine:0': 36, 'quadratic:0': 36, 'cubic:0': 36}
9:b174 {'affine:0': 18, 'quadratic:0': 18, 'cubic:0': 18}
```

## K3. Claude 2: the U1 bound

```
$ sed -n '58p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md | cut -c1-400
| U1 interior `0≤f≤1`, permitted jumps and no concentrated input fields; S12 encoding and user "Leave it out for now", 2026-10-20 | The list says U1 and its encoding are **"as the S12 decision list states them"**. | The regional/traced function class of the same objects. The no-concentrated-input restriction does not remove concentrated source parts. |
```

```
$ cut -c1-200 /var/projects/_review_scratch/s14a1-record-r0-claude/chi_bound.stdout
scripts/S14a_1_projection.py | lines mentioning chi: 7 | chi lines with an inequality/assumption token: [] | all Assumptions/positive declarations: [(46, 'r, a, R = sp.symbols("radius inner_radius out
mathematica/S14a_1_projection.wl | lines mentioning chi: 2 | chi lines with an inequality/assumption token: [] | all Assumptions/positive declarations: [(58, 'angularDensity = FullSimplify[Abs[spatial
scripts/out/S14a_1_projection.out | lines mentioning chi: 108 | chi lines with an inequality/assumption token: [] | all Assumptions/positive declarations: []
mathematica/out/S14a_1_projection.out | lines mentioning chi: 54 | chi lines with an inequality/assumption token: [] | all Assumptions/positive declarations: []
```

```
$ grep -nE 'enters no computed object' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_record_decision_list.md
93:   The record gives no standing the list does not give. One that enters no computed object is listed as such.
```

## K4. Claude 3: the quoted instruction

```
$ sed -n '352,355p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md
- The record contract calls the measurements generator **"a committed generator script"** and
  also excludes **"No commit and no push"**; the current user expressly says **"Commit nothing."**
  The durable generator is supplied in the working tree, uncommitted, with the final regenerated
  measurements. No acceptance/commit is invented. (Retrieve M01/M16.)
```

```
$ grep -nF 'Commit nothing' /var/projects/toy_physics/_scratch/s14a/rec/brief_s14a1_record_author.md | cut -c1-140
5:Other jobs may run on this machine and write under `_scratch/`; leave their files alone. Put your generator's working files under `_scratc
```

```
$ grep -cF 'Commit nothing' /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S14a_1_record_measurements.md
0
```

## K5. Claude 4: the channels' name

```
$ sed -n '10p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md
specialization**, the record reports the two material-loss channels, their sum `J`, the signed
```

```
$ grep -nE 'material-loss|material loss' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_decision_list.md
/var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md:10:specialization**, the record reports the two material-loss channels, their sum `J`, the signed
```

```
$ sed -n '54p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_decision_list.md
`J` counts **both channels** by which material leaves the slab around a throat:
```

```
$ sed -n '87p' /var/projects/toy_physics/research/pde_ledger_v3/../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md
re-ordering (return). Total material is conserved; brane-ORDER is not.
```
