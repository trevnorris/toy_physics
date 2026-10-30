# Measurements: S15 reading build review round 2 (round 3 of the build), lookups

Generator: `_scratch/s15/gen/s15_reading_r2_lookups.sh` (GNU grep, sed, awk field selection and literal-match counts,
git, sha256sum, `LC_ALL=C`). Each block is a command and its literal output: existence, verbatim retrieval, literal-match
counts, hashes. The generator first copies, by literal tag match, the item-4 lines of the objects whose item-4 status
is the finite-identity text into `_scratch/s15/reading_review/r2_lookup_extracts/`; M3 counts in those copies.
Computed results are quoted from the legs' saved stdout; their scripts were not re-run. The legs' closed forms and
conditions, and their lines that name a stratum where an identity holds, are not quoted here.

## M0. The round

```
$ sha256sum _scratch/s15/reading_review/s15_reading_build_review_prompt_r2.md
be51c0983bb968627b27db0ed0f33cd3342ccb318dcca3a337fa13cb24b63d00  _scratch/s15/reading_review/s15_reading_build_review_prompt_r2.md
```

```
$ sha256sum -c _scratch/s15/reading_review/s15_reading_build_baseline_r2.sha256
research/pde_ledger_v3/scripts/S15_separation_departure.py: OK
research/pde_ledger_v3/mathematica/S15_separation_departure.wl: OK
research/pde_ledger_v3/scripts/out/S15_separation_departure.out: OK
research/pde_ledger_v3/mathematica/out/S15_separation_departure.out: OK
research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md: OK
```

```
$ cat _scratch/s15/reading_review/r2_reports.sha256; sha256sum -c _scratch/s15/reading_review/r2_reports.sha256
87eb7516be3fec64a5911e2d4abf593b24e56957b66ab90b01c9bf9365954c61  /var/projects/_review_scratch/s15-reading-r2-grok/REVIEW.md
ed049e648266082031fc226223e8d9c04378ab58902f1371099bc3246c2368c5  /var/projects/_review_scratch/s15-reading-r2-claude/REVIEW.md
/var/projects/_review_scratch/s15-reading-r2-grok/REVIEW.md: OK
/var/projects/_review_scratch/s15-reading-r2-claude/REVIEW.md: OK
```

```
$ cat /var/projects/_review_scratch/s15-reading-r2-grok/DONE; wc -w /var/projects/_review_scratch/s15-reading-r2-claude/REVIEW.md /var/projects/_review_scratch/s15-reading-r2-grok/REVIEW.md
exit 0
 1721 /var/projects/_review_scratch/s15-reading-r2-claude/REVIEW.md
 1244 /var/projects/_review_scratch/s15-reading-r2-grok/REVIEW.md
 2965 total
```

```
$ git rev-parse HEAD; git status --short
6c0fc66bab0156714e48ca8ce7c52f3bec7fedc4
 M .claude/skills/prep-for-compact/SKILL.md
 M research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
 T research/pde_ledger_v3/scripts/out/S15_separation_departure.out
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair2_brief_lookups.md
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r2_disposition.md
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r2_lookups.md
```

```
$ diff <(git status --short | grep -v -F S15_reading_build_review_r2_; git rev-parse HEAD) _scratch/s15/reading_review/status_before_r2.txt && echo 'status unchanged apart from the r2 disposition files'
status unchanged apart from the r2 disposition files
```

```
$ find . -path ./.git -prune -o -path ./_scratch -prune -o -newer _scratch/s15/reading_review/s15_reading_build_review_prompt_r2.md -type f -print | grep -v -F S15_reading_build_review_r2_
```

```
$ ls -la --time-style=full-iso /var/projects/_review_scratch/s15-reading-r2-claude/stray_repo_root/x02_PYFORMA.stdout /var/projects/_review_scratch/s15-reading-r2-claude/x02_PYFORMA.stdout | awk '{print $5, $6, $7, $9}'; cat /var/projects/_review_scratch/s15-reading-r2-claude/stray_repo_root/NOTE.txt
0 2026-10-30 00:47:36.645963566 /var/projects/_review_scratch/s15-reading-r2-claude/stray_repo_root/x02_PYFORMA.stdout
400 2026-10-30 00:47:47.304865459 /var/projects/_review_scratch/s15-reading-r2-claude/x02_PYFORMA.stdout
Moved by the orchestrator on 2026-10-30 from /var/projects/toy_physics/x02_PYFORMA.stdout (untracked, 0 bytes, mtime 2026-10-30 00:47:36 -0600; the leg's own x02_PYFORMA.stdout in this directory was written 11 s later). mv preserved the mtime.
```

```
$ for j in $(ls _scratch/s11c/s15-reading-review-r2-guard); do if [ -f _scratch/s11c/s15-reading-review-r2-guard/$j/outcome.json ]; then echo "$j: $(tr -d ' \n' < _scratch/s11c/s15-reading-review-r2-guard/$j/outcome.json | cut -c1-110)"; else echo "$j: no outcome.json"; fi; done
claude-d01: {"exitCode":1,"wallSeconds":29084.17429874302,"unit":"s11c-guard-0213e015ad83","stderrBytes":0,"limitsVerified
claude-d01b: {"exitCode":1,"wallSeconds":3018.9931481759995,"unit":"s11c-guard-8cf5e6eae511","stderrBytes":0,"limitsVerifie
claude-d01c: {"exitCode":0,"wallSeconds":12.951186239952222,"unit":"s11c-guard-e35b0d4ce73b","stderrBytes":0,"limitsVerifie
claude-d01c2: {"exitCode":0,"wallSeconds":942.435704635107,"unit":"s11c-guard-b91d1647eef2","stderrBytes":0,"limitsVerified"
claude-d01d: {"exitCode":0,"wallSeconds":866.8338197730482,"unit":"s11c-guard-7f0369c167db","stderrBytes":0,"limitsVerified
claude-d02b: {"exitCode":0,"wallSeconds":896.6318780459696,"unit":"s11c-guard-04b5512fa1c4","stderrBytes":0,"limitsVerified
claude-e04: {"exitCode":0,"wallSeconds":231.44767329806928,"unit":"s11c-guard-6540c9242e48","stderrBytes":0,"limitsVerifie
claude-e04b: {"exitCode":0,"wallSeconds":231.2068539840402,"unit":"s11c-guard-d53d0d3ca423","stderrBytes":0,"limitsVerified
claude-e05b: {"exitCode":0,"wallSeconds":10408.938541523064,"unit":"s11c-guard-cecb79e19c70","stderrBytes":0,"limitsVerifie
claude-e05r: {"exitCode":0,"wallSeconds":1684.1772257740377,"unit":"s11c-guard-37e317abe420","stderrBytes":328,"limitsVerif
claude-e06: {"exitCode":0,"wallSeconds":1332.0548322820105,"unit":"s11c-guard-61878813f7e1","stderrBytes":0,"limitsVerifie
claude-e07: {"exitCode":0,"wallSeconds":5.2408667099662125,"unit":"s11c-guard-fa383f983a62","stderrBytes":0,"limitsVerifie
claude-e07b: {"exitCode":0,"wallSeconds":1166.5472621570807,"unit":"s11c-guard-cd76f9f13dc5","stderrBytes":0,"limitsVerifie
claude-e07r1: {"exitCode":0,"wallSeconds":5.084179715020582,"unit":"s11c-guard-9f068aa89012","stderrBytes":0,"limitsVerified
claude-e07r1b: {"exitCode":0,"wallSeconds":2240.9607567499625,"unit":"s11c-guard-33683c78594e","stderrBytes":0,"limitsVerifie
claude-e09: {"exitCode":0,"wallSeconds":40.951067727059126,"unit":"s11c-guard-60cb36e0082e","stderrBytes":0,"limitsVerifie
claude-e10: {"exitCode":1,"wallSeconds":4834.816985265003,"unit":"s11c-guard-58ceb89ad4b7","stderrBytes":0,"limitsVerified
claude-e10b: {"exitCode":1,"wallSeconds":548.1918297590455,"unit":"s11c-guard-774b8ccc91ff","stderrBytes":0,"limitsVerified
claude-e10c: {"exitCode":0,"wallSeconds":23.595504031982273,"unit":"s11c-guard-952c787c8cb6","stderrBytes":0,"limitsVerifie
claude-e10d: {"exitCode":0,"wallSeconds":59.30550533905625,"unit":"s11c-guard-1696ac3f680a","stderrBytes":0,"limitsVerified
claude-e10e: {"exitCode":0,"wallSeconds":24.54095304198563,"unit":"s11c-guard-cd38103db6d3","stderrBytes":0,"limitsVerified
claude-py_base: {"exitCode":0,"wallSeconds":456.8729247279698,"unit":"s11c-guard-5cac9f00bc34","stderrBytes":0,"limitsVerified
claude-py_depB: {"exitCode":1,"wallSeconds":2399.100346658961,"unit":"s11c-guard-6268de18dd14","stderrBytes":0,"limitsVerified
claude-py_depC: {"exitCode":0,"wallSeconds":569.644652256975,"unit":"s11c-guard-24ca8bdbaa3e","stderrBytes":0,"limitsVerified"
claude-py_famA: {"exitCode":0,"wallSeconds":1347.7307644890388,"unit":"s11c-guard-b8b79efcea0d","stderrBytes":0,"limitsVerifie
claude-py_formA: {"exitCode":0,"wallSeconds":1132.011970143998,"unit":"s11c-guard-1320789c6d70","stderrBytes":0,"limitsVerified
claude-rerun_py2: {"exitCode":0,"wallSeconds":12708.204215775011,"unit":"s11c-guard-b228c2bec280","stderrBytes":0,"limitsVerifie
claude-rerun_wl: {"exitCode":0,"wallSeconds":1513.448981732945,"unit":"s11c-guard-d5322bd3820e","stderrBytes":0,"limitsVerified
claude-t00: {"exitCode":0,"wallSeconds":0.08791353495325893,"unit":"s11c-guard-98bcf63b29cc","stderrBytes":0,"limitsVerifi
claude-t01: {"exitCode":1,"wallSeconds":0.686930222902447,"unit":"s11c-guard-d964bc56fef7","stderrBytes":1214,"limitsVerif
claude-t01b: {"exitCode":1,"wallSeconds":0.9854768960503861,"unit":"s11c-guard-fceecc2abc17","stderrBytes":1251,"limitsVeri
claude-t01c: {"exitCode":0,"wallSeconds":4.483181147952564,"unit":"s11c-guard-9fb85f9ad2d1","stderrBytes":0,"limitsVerified
claude-t02: {"exitCode":0,"wallSeconds":11.090511915041134,"unit":"s11c-guard-5dde88d9a192","stderrBytes":0,"limitsVerifie
claude-wl_balC: {"exitCode":0,"wallSeconds":1634.651872509974,"unit":"s11c-guard-2444438d2d2c","stderrBytes":0,"limitsVerified
claude-wl_depB: {"exitCode":0,"wallSeconds":1744.8813198480057,"unit":"s11c-guard-6f5bbbb38fee","stderrBytes":0,"limitsVerifie
claude-wl_formA: {"exitCode":0,"wallSeconds":4268.064500155044,"unit":"s11c-guard-27beaf21c1cc","stderrBytes":0,"limitsVerified
claude-x04g: {"exitCode":0,"wallSeconds":4.190993643016554,"unit":"s11c-guard-3a20ee196db7","stderrBytes":0,"limitsVerified
claude-x11: {"exitCode":1,"wallSeconds":26.948316542897373,"unit":"s11c-guard-d88bde586731","stderrBytes":519,"limitsVerif
claude-x11b: {"exitCode":0,"wallSeconds":378.28763447806705,"unit":"s11c-guard-3eff52d25ef3","stderrBytes":0,"limitsVerifie
grok-derive-01: {"exitCode":0,"wallSeconds":397.7550922851078,"unit":"s11c-guard-ae4dd8a0c51d","stderrBytes":0,"limitsVerified
grok-derive-02: {"exitCode":0,"wallSeconds":17.946922942996025,"unit":"s11c-guard-0d437f280dec","stderrBytes":0,"limitsVerifie
grok-derive-03: {"exitCode":1,"wallSeconds":1.3392503350041807,"unit":"s11c-guard-d2b8f474116f","stderrBytes":2257,"limitsVeri
grok-derive-04: no outcome.json
```

```
$ for j in $(ls _scratch/s11c/s15-reading-review-r2-guard | grep -E '^claude-(d01|e10|e10b|d01b|py_depB|e05b|e07r1b|rerun_py2)$'); do echo "$j $(grep -o -E '\"(memoryMax|startedUtc)\": [^,]*' _scratch/s11c/s15-reading-review-r2-guard/$j/invocation.json | tr '\n' ' ')"; done
claude-d01 "memoryMax": 21474836480 "startedUtc": "2026-10-29T17:50:04.651859+00:00" 
claude-d01b "memoryMax": 26843545600 "startedUtc": "2026-10-30T02:08:10.904492+00:00" 
claude-e05b "memoryMax": 26843545600 "startedUtc": "2026-10-29T19:40:27.927997+00:00" 
claude-e07r1b "memoryMax": 26843545600 "startedUtc": "2026-10-29T23:04:41.591666+00:00" 
claude-e10 "memoryMax": 26843545600 "startedUtc": "2026-10-30T00:34:13.976469+00:00" 
claude-e10b "memoryMax": 26843545600 "startedUtc": "2026-10-30T01:55:14.611571+00:00" 
claude-py_depB "memoryMax": 26843545600 "startedUtc": "2026-10-30T06:47:34.677914+00:00" 
claude-rerun_py2 "memoryMax": 26843545600 "startedUtc": "2026-10-30T08:00:13.917473+00:00" 
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; systemctl --user list-units 's11c-guard-*' --all --no-legend | wc -l
{}

0
```

```
$ sha256sum _scratch/s15/reading_review/r2_registry_before_reconcile_grok-derive-04.json; grep -c '"s11c-guard-' _scratch/s15/reading_review/r2_registry_before_reconcile_grok-derive-04.json; grep -o '"logDirectory": "[^"]*"' _scratch/s15/reading_review/r2_registry_before_reconcile_grok-derive-04.json
29f5bb2519dd46b97759b4c85b920ede6a5a1835a3c7748e4eb3874bef9a4686  _scratch/s15/reading_review/r2_registry_before_reconcile_grok-derive-04.json
2
"logDirectory": "/var/projects/toy_physics/_scratch/s11c/s15-reading-review-r2-guard/claude-d01"
"logDirectory": "/var/projects/toy_physics/_scratch/s11c/s15-reading-review-r2-guard/grok-derive-04"
```

```
$ sed -n '1,6p' _scratch/s15/reading_review/r2_claude_amendment_sendmessage.json; python3 -c "import json;print(json.load(open('_scratch/s15/reading_review/r2_claude_amendment_sendmessage.json'))['input']['message'])"
{
  "timestamp": "2026-10-30T01:49:31.086Z",
  "input": {
    "to": "a82b3a4672a508796",
    "summary": "Permit the Claude leg to stop its own stalled runs",
    "message": "Operational note from the orchestrator (no findings, no other leg's output): your guarded runs claude-d01 (unit s11c-guard-0213e015ad83, 20 GiB, no stdout since 12:08) and claude-e10 (unit s11c-guard-58ceb89ad4b7, 25 GiB, no stdout since 18:34) together hold enough of the pool that your rerun cannot be admitted. This amends the prompt's \"Do not stop guarded runs\" for your own runs only: you may stop any of your own guarded runs that you judge stalled, with `systemctl --user stop <unit>` while its launcher process is still alive. Then confirm `systemctl --user is-active <unit>` reports inactive and that its entry is gone from ~/.local/state/s11c-guard/resource-pool/reservations.json. Do not edit the registry, and do not touch runs that are not yours. Whether to stop either run is your call. Report in REVIEW.md what you stopped and what it had printed. Everything else in the prompt stands."
Operational note from the orchestrator (no findings, no other leg's output): your guarded runs claude-d01 (unit s11c-guard-0213e015ad83, 20 GiB, no stdout since 12:08) and claude-e10 (unit s11c-guard-58ceb89ad4b7, 25 GiB, no stdout since 18:34) together hold enough of the pool that your rerun cannot be admitted. This amends the prompt's "Do not stop guarded runs" for your own runs only: you may stop any of your own guarded runs that you judge stalled, with `systemctl --user stop <unit>` while its launcher process is still alive. Then confirm `systemctl --user is-active <unit>` reports inactive and that its entry is gone from ~/.local/state/s11c-guard/resource-pool/reservations.json. Do not edit the registry, and do not touch runs that are not yours. Whether to stop either run is your call. Report in REVIEW.md what you stopped and what it had printed. Everything else in the prompt stands.
```

## M1. The contract and the round 1 obligations

```
$ git show 7f52b4d2:research/pde_ledger_v3/directives/S15_reading_decision_list.md | sed -n '78,82p;86,88p'
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
   - `Q` on the family under each condition.
Every object is exact on S15's domain under the conditions it carries: S15's pressure conditions, and for each
departure object a condition under which it exists, as above. No limit or expansion is taken in place of an exact
object. Where an engine does not obtain an object exactly, it prints a label saying so, beside what it did print.
```

```
$ sed -n '84,91p' research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r1_disposition.md
- **What must be true.**
  - Part B item 4 holds for each of the 48 objects in each engine. The conditions are computed exactly from the
    printed closed form, with every special exponent value treated. They are printed as equations and inequalities
    in the named quantities, and `Q` on the family is printed under each.
  - A return that gives conditions is not filed as not obtained. A label appears only where an exact attempt that
    could have returned the conditions did not, and the label says what was attempted, beside what was printed.
  - Each engine reaches its conditions by its own construction.
  - No item-1, item-2, item-5 or item-6 object changes, and no item-3 object changes except as W2 requires.
```

```
$ sed -n '3,10p' _scratch/s15/reading_build/resume_s15_reading_repair_2.md
- W1, both engines: item 4 is not computed for the 48 closed forms that depend on the separation, and neither engine's attempt can return it.
  - SymPy abbreviates the kernels that hold `d` before the solver (`research/pde_ledger_v3/scripts/S15_separation_departure.py:776`). It prints each `solve_undetermined_coeffs` return and collects it in `attempted` (`:782–787`), and nothing else reads it. Then `complete=(not reduced.has(d,…))` (`:788`), and for a remainder that holds `d`, `solutions=[]` (`:820`). Whatever the solver returns, no condition is printed.
  - Wolfram poses `SolveAlways` on an `Implies` form (`research/pde_ledger_v3/mathematica/S15_separation_departure.wl:492`, `:499`). A returned list is kept as `answer` (`:501–507`), and `:510` files an unevaluated return and a returned list alike as not obtained.
  - Independent review reached the conditions, and `Q` under them, for representative objects of these classes by exact routes, and checked them against a second route. The gaps are not limits of the objects.
  - What must be true:
    - Part B item 4 holds for each of the 48 objects in each engine. The conditions are computed exactly from the printed closed form, with every special exponent value treated. They are printed as equations and inequalities in the named quantities, and `Q` on the family is printed under each.
    - A return that gives conditions is not filed as not obtained. A label appears only where an exact attempt that could have returned the conditions did not, and the label says what was attempted, beside what was printed.
    - Each engine reaches its conditions by its own construction.
```

```
$ sed -n '90,91p;113,120p' research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r0_disposition.md
**V2. Item 4's condition is the closed object set equal to zero, and "Q under the condition" is the family's `Q`.**
(Claude F2, Grok 2.) **CONFIRMED.**
- **What must be true.**
  - Part B item 4 holds for every closed form that item 3 prints, including those that depend on the separation. The
    conditions are reached by computation from the closed form and printed as equations and inequalities in the
    named quantities, and `Q` on the family is printed under each.
  - Where an engine does not obtain the conditions exactly after attempting them, it prints a label saying so beside
    what it did print (M1a).
  - No status line or run-record sentence states a property of a printed object that no computation in that engine
    produced.
```

## M2. The engines at the r2 baseline

```
$ sed -n '30,36p;293,303p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-160
  The scalar output solution space has dimension r. Its Wronskians differ by
  constant factors. A nonzero output of zero order k makes their order at
  least k-r+1. Thus k<=r*m+g*r*(r-1)/2+r-1, bounded by the same expression
  with the full system dimension N. This includes rank drops and the zero
  output without dividing by any parameter-dependent minor.
  Consequently the first B+1 exact coefficient equations at any regular
  point are necessary AND sufficient for the output to vanish identically.
    T=IdentityTransfer(transfer_name,index,size)
    R=S.MatrixSymbol(projection_name,1,size)
    Y=S.MatrixSymbol(initial_name,size,1)
    product=OrderedMatrixProduct(T,index,0,power-1)
    scalar=MatrixElement(S.MatMul(R,product,Y,evaluate=False),0,0)
    # Positive reference powers give every summand the same dimension.
    equation=S.Eq(S.Sum(S.Abs(reference**power*scalar)**2,(power,0,bound)),0,evaluate=False)
    output('matrix_aliases',{transfer_name:dict(shape=(size,size),entry_table=label+'/transfer_entries',default_entry=0,index=index),projection_name:dict(shape=
    output('ordered_product_definition','OrderedMatrixProduct(T,j,l,u)=T(u) ... T(l), with the identity matrix for u<l. The index is an integer; all sums and pr
    equation=S.And(regularity,equation)
    output('condition',equation)
```

```
$ sed -n '1093,1110p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
        if reduced.has(d):
            condition_parts=[]
            for powers,coefficient in reduced_coefficients.items():
                scalar=coefficient.xreplace(backward).xreplace(back_jets)
                scalar_tag=tag+'/finite_identity_'+''.join(map(str,powers))
                def identity_progress(step):
                    print('WORKING_ON: '+step,flush=True)
                    emit('WORKING_ON',step)
                condition,certificate=finite_identity_condition(scalar,d,identity_reference_point,scalar_tag,emit,identity_progress)
                condition_parts.append(condition)
            condition=S.And(rout>rin,rin>0,*condition_parts)
            solo(tag+'/condition',condition)
            reference_value=qclosed.subs(d,identity_reference_point)
            solo(tag+'/restriction_reference_operand',reference_value)
            restricted=identity_euler_template.subs({identity_distance:d,identity_reference:identity_reference_point,identity_value:reference_value},simultaneous=True)
            solo(tag+'/Q_under_condition',restricted)
            solo(tag+'/Q_on_family',qclosed)
            emit(tag+'/status','Finite parameter equations computed from the closed-form differential system and unit-sphere coefficient reduction. Q under the condition is computed from the Euler equ
```

```
$ sed -n '605p;610,611p;613,624p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-220
 productName=Unique["identityTransfer"];projectionName=Unique["identityProjection"];startName=Unique["identityInitial"];
 condition=Inactive[Sum][Abs[rv^ord matrixEntry[Inactive[Dot][projectionName,orderedMatrixProduct[productName,idx,0,ord-1],startName],1,1]]^2,{ord,0,bound}]==0;
 solo[prefix<>"/condition",condition];And[regularity,condition]];
(* Multiplicity certificate. For H Y'=M Y and f=P Y put
 P_0=P, P_(j+1)=H P_j'+P_j M-j H' P_j. Then H^j f^(j)=P_j Y.
 With m=deg P and g=max(0,deg H-1,deg M), deg P_j<=m+j g.
 If their rational rank is r, the first r rows span all subsequent rows.
 Every minor of their wedge has degree <=r m+g r(r-1)/2. Multiplication
 by an invertible analytic fundamental matrix preserves the minimum order
 of vanishing of that wedge. Its nonzero Wronskians differ by constants.
 A scalar output with zero order k therefore has k<=r m+g r(r-1)/2+r-1.
 The printed bound uses size>=r and remains valid when any parameter causes
 a rank drop. Thus the finite coefficient equations are equivalent to the
 identity, with no division by a parameter-dependent rank minor. *)
emit["item4_kernel/finite_identity_certificate","For the printed polynomial system H Y'=M Y and scalar numerator P Y, the finite multiplicity bound is N deg(P)+max(0,deg(H)-1,deg(M)) N(N-1)/2+N-1. The source gives the we
```

```
$ sed -n '655,667p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-220
 coefficients=(Last/@CoefficientRules[reduction[[2]],directionDifferentiationSymbols])/.Thread[slots->values];
 solo[prefix<>"/scalar_identity_operands",coefficients];
 Do[result=Block[{currentObject=prefix<>"/scalar_"<>ToString[k]},Catch[finiteIdentity[coefficients[[k]],dist,identityReferencePoint,prefix<>"/scalar_"<>ToString[k]],identityFailure]];
  If[Head[result]===List,solo[prefix<>"/attempt_return",result];
   emit[prefix<>"/status","NOT OBTAINED EXACTLY: differential closure and finite multiplicity elimination were attempted; the displayed unsupported primitive or algebraic step prevented completion."];
   AppendTo[exactGaps,{name,"item4",result}];Return[Null]];
  AppendTo[conditions,result],{k,Length[coefficients]}];
 solo[prefix<>"/condition",rout[t]>rin[t]>0 && And@@conditions];
 referenceValue=q/.dist->identityReferencePoint;
 solo[prefix<>"/restriction_reference_operand",referenceValue];
 restricted=identityRestrictionTemplate/.{identityD->dist,identityReference->identityReferencePoint,identityValue->referenceValue};
 solo[prefix<>"/Q_under_condition",restricted];solo[prefix<>"/Q_on_family",q];
 emit[prefix<>"/status","Finite parameter equations computed from the closed-form differential system and exact direction coefficient elimination. Q under that condition is computed from the Euler equation with its symbo
```

```
$ for c in 6c0fc66b WORKTREE; do if [ $c = WORKTREE ]; then py=$(cat research/pde_ledger_v3/scripts/S15_separation_departure.py); wl=$(cat research/pde_ledger_v3/mathematica/S15_separation_departure.wl); else py=$(git show $c:research/pde_ledger_v3/scripts/S15_separation_departure.py); wl=$(git show $c:research/pde_ledger_v3/mathematica/S15_separation_departure.wl); fi; echo "$c: py lines matching finite_identity|IdentityTransfer|OrderedMatrixProduct: $(echo "$py" | grep -c -E 'finite_identity|IdentityTransfer|OrderedMatrixProduct'); wl lines matching finiteIdentity|orderedMatrixProduct|identityTransfer: $(echo "$wl" | grep -c -E 'finiteIdentity|orderedMatrixProduct|identityTransfer')"; done
6c0fc66b: py lines matching finite_identity|IdentityTransfer|OrderedMatrixProduct: 0; wl lines matching finiteIdentity|orderedMatrixProduct|identityTransfer: 0
WORKTREE: py lines matching finite_identity|IdentityTransfer|OrderedMatrixProduct: 9; wl lines matching finiteIdentity|orderedMatrixProduct|identityTransfer: 5
```

## M3. Item 4 of the separation-dependent objects in the transcripts

```
$ grep -E '^reading/.*/item4/status: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | sed 's/^[^:]*: //' | cut -c1-60 | sort | uniq -c
     96 Computed factor zero-set conditions and Q after their substi
     48 Finite parameter equations computed from the closed-form dif
    192 Not applicable: the printed item-3 form retains a general fi
```

```
$ grep -E '^reading/.*/item4/status: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | sed 's/^[^:]*: //' | cut -c1-60 | sort | uniq -c
     48 "Finite parameter equations computed from the closed-form di
    192 "Not applicable to the general-field family integral."
     96 "Reduce conditions and Q simplified under each computed cond
```

```
$ wc -l _scratch/s15/reading_review/r2_lookup_extracts/py_prefixes.txt _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt _scratch/s15/reading_review/r2_lookup_extracts/py_Q_under_condition.txt _scratch/s15/reading_review/r2_lookup_extracts/py_condition_numbered_or_union.txt _scratch/s15/reading_review/r2_lookup_extracts/wl_prefixes.txt _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt _scratch/s15/reading_review/r2_lookup_extracts/wl_Q_under_condition.txt _scratch/s15/reading_review/r2_lookup_extracts/wl_condition_numbered_or_union.txt
     48 _scratch/s15/reading_review/r2_lookup_extracts/py_prefixes.txt
     48 _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt
     48 _scratch/s15/reading_review/r2_lookup_extracts/py_Q_under_condition.txt
      0 _scratch/s15/reading_review/r2_lookup_extracts/py_condition_numbered_or_union.txt
     48 _scratch/s15/reading_review/r2_lookup_extracts/wl_prefixes.txt
     48 _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt
     48 _scratch/s15/reading_review/r2_lookup_extracts/wl_Q_under_condition.txt
      0 _scratch/s15/reading_review/r2_lookup_extracts/wl_condition_numbered_or_union.txt
    288 total
```

```
$ sed 's#/item4$##' _scratch/s15/reading_review/r2_lookup_extracts/py_prefixes.txt | sed -E 's#^reading/##; s#/component_[123]/#/component_n/#' | sort | uniq -c
      3 inner/carried/component_n/labels_drain1
      3 inner/carried/component_n/labels_drain1+drain2
      3 inner/transport_centre/component_n/labels_drain1
      3 inner/transport_centre/component_n/labels_drain1+drain2
      3 inner/transport_radius/component_n/labels_drain1
      3 inner/transport_radius/component_n/labels_drain1+drain2
      3 outer/carried/component_n/labels_drain1
      3 outer/carried/component_n/labels_drain1+drain2
      3 outer/transport_centre/component_n/labels_drain1
      3 outer/transport_centre/component_n/labels_drain1+drain2
      3 outer/transport_radius/component_n/labels_drain1
      3 outer/transport_radius/component_n/labels_drain1+drain2
      3 region/component_n/momentum_labels_drain1
      3 region/component_n/momentum_labels_drain1+drain2
      3 region/component_n/rate_labels_drain1
      3 region/component_n/rate_labels_drain1+drain2
```

```
$ echo "py condition lines containing 'Eq(Sum(': $(grep -c -F 'Eq(Sum(' _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt)"; echo "py condition lines whose only Eq( are Eq(Sum(: $(awk '{a=gsub(/Eq\(/,"&"); b=gsub(/Eq\(Sum\(/,"&"); if (a==b) n++} END{print n+0}' _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt)"; for p in '(^|[^A-Za-z0-9_])[ab]_[12]([^A-Za-z0-9_]|$)' '[AB]_[12]\(t\)' 'Derivative\('; do echo "py condition lines matching $p: $(grep -c -E "$p" _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt)"; done
py condition lines containing 'Eq(Sum(': 48
py condition lines whose only Eq( are Eq(Sum(: 48
py condition lines matching (^|[^A-Za-z0-9_])[ab]_[12]([^A-Za-z0-9_]|$): 0
py condition lines matching [AB]_[12]\(t\): 0
py condition lines matching Derivative\(: 0
```

```
$ echo "wl condition lines containing 'Inactive[Sum]': $(grep -c -F 'Inactive[Sum]' _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt)"; echo "wl condition lines whose == count equals their Inactive[Sum] count: $(awk '{a=gsub(/==/,"&"); b=gsub(/Inactive\[Sum\]/,"&"); if (a==b) n++} END{print n+0}' _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt)"; for p in '(aa|bb)\[[12]\]' '(AA|BB)\[[12]\]' 'Derivative\['; do echo "wl condition lines matching $p: $(grep -c -E "$p" _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt)"; done
wl condition lines containing 'Inactive[Sum]': 48
wl condition lines whose == count equals their Inactive[Sum] count: 48
wl condition lines matching (aa|bb)\[[12]\]: 0
wl condition lines matching (AA|BB)\[[12]\]: 0
wl condition lines matching Derivative\[: 0
```

```
$ echo 'py Sum upper bounds (largest five):'; grep -o -E '\(identity_order, 0, [0-9]+\)' _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt | grep -o -E '[0-9]+\)$' | tr -d ')' | sort -n | tail -5 | tr '\n' ' '; echo; echo 'wl Sum upper bounds (largest five):'; grep -o -E '\{ord\$[0-9]+, 0, [0-9]+\}' _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt | grep -o -E '[0-9]+\}$' | tr -d '}' | sort -n | tail -5 | tr '\n' ' '; echo
py Sum upper bounds (largest five):
358488 358488 825513 825513 825513 
wl Sum upper bounds (largest five):
817596 817596 2249810 2249810 2249810 
```

```
$ grep -F 'reading/outer/carried/component_1/labels_drain1/item4/condition: ' _scratch/s15/reading_review/r2_lookup_extracts/py_condition.txt | cut -c1-600
reading/outer/carried/component_1/labels_drain1/item4/condition: (r_in(t) > 0) & (R(t) <= R(t)) & (R(t) > r_in(t)) & (-R(t) <= R(t)) & Eq(Sum(Abs((2*R(t))**identity_order*(reading_outer_carried_component_1_labels_drain1_item4_finite_identity_100_P*OrderedMatrixProduct(IdentityTransfer(reading_outer_carried_component_1_labels_drain1_item4_finite_identity_100_T, identity_index, 42), identity_index, 0, identity_order - 1)*reading_outer_carried_component_1_labels_drain1_item4_finite_identity_100_Y)[0, 0])**2, (identity_order, 0, 83)), 0) & Eq(Sum(Abs((2*R(t))**identity_order*(reading_outer_carried
```

```
$ grep -F 'reading/outer/carried/component_1/labels_{drain1}/item4/condition: ' _scratch/s15/reading_review/r2_lookup_extracts/wl_condition.txt | cut -c1-600
reading/outer/carried/component_1/labels_{drain1}/item4/condition: rout[t] > rin[t] > 0 && Inactive[Sum][2^(2*Re[ord$1311198])*Abs[matrixEntry[Inactive[Dot][identityProjection1704, orderedMatrixProduct[identityTransfer1703, idx$1311198, 0, -1 + ord$1311198], identityInitial1705], 1, 1]*rout[t]^ord$1311198]^2, {ord$1311198, 0, 1259}] == 0 && Inactive[Sum][2^(2*Re[ord$1312555])*Abs[matrixEntry[Inactive[Dot][identityProjection1707, orderedMatrixProduct[identityTransfer1706, idx$1312555, 0, -1 + ord$1312555], identityInitial1708], 1, 1]*rout[t]^ord$1312555]^2, {ord$1312555, 0, 0}] == 0 && Inactive
```

```
$ awk '{i=index($0,": "); print substr($0,1,i-1)}' _scratch/s15/reading_review/r2_lookup_extracts/py_Q_under_condition.txt | sort | uniq -c | awk '{print $1}' | sort | uniq -c; awk '{i=index($0,": "); print substr($0,1,i-1)}' _scratch/s15/reading_review/r2_lookup_extracts/wl_Q_under_condition.txt | sort | uniq -c | awk '{print $1}' | sort | uniq -c
     48 1
     48 1
```

```
$ grep -F 'reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition: ' _scratch/s15/reading_review/r2_lookup_extracts/py_Q_under_condition.txt | cut -c1-300; grep -F 'reading/outer/carried/component_1/labels_{drain1}/item4/Q_under_condition: ' _scratch/s15/reading_review/r2_lookup_extracts/wl_Q_under_condition.txt | cut -c1-300
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition: 4*(-4*N_1**2*(Piecewise((-pi*N_1*(3*R(t))**(2*a_1 + b_1 + 2)/(4*(2*a_1 + b_1 - 2)*R(t)**4 + 16*R(t)**4) + pi*N_1*R(t)**(2*a_1 + b_1 + 2)/(4*(2*a_1 + b_1 - 2)*R(t)**4 + 16*R(t)**4), Ne(2*a_1 + b_1, -2)), (pi*N_1*log(R(t))/(4*R(t
reading/outer/carried/component_1/labels_{drain1}/item4/Q_under_condition: (4*rout[t]^2*(-4*directionSeparation[1]^2*(Piecewise[{{(Pi*directionSeparation[1]*Log[rout[t]])/(4*rout[t]^4) - (Pi*directionSeparation[1]*Log[3*rout[t]])/(4*rout[t]^4), aa[1] + bb[1]/2 == -1}}, (Pi*directionSeparation[1]*rou
```

```
$ grep -m2 -E '^item4_kernel/(Euler_equation_return|finite_identity_certificate): ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-300; grep -m2 -E '^item4_kernel/(Euler_equation_return|finite_identity_certificate|restriction_template)[^:]*: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-300
item4_kernel/Euler_equation_return: Eq(identity_function(identity_distance), identity_reference**2*identity_value/identity_distance**2)
item4_kernel/finite_identity_certificate: For H dY/dd=M Y and f=P Y the source proves a zero-multiplicity bound N deg(P)+max(0,deg(H)-1,deg(M)) N(N-1)/2+N-1, including every rank degeneration. The finite transfer equations test exactly those coefficients at symbolic reference 2 R(t). Only finite pol
item4_kernel/finite_identity_certificate: "For the printed polynomial system H Y'=M Y and scalar numerator P Y, the finite multiplicity bound is N deg(P)+max(0,deg(H)-1,deg(M)) N(N-1)/2+N-1. The source gives the wedge/Wronskian proof, including rank drops. Transfer entries are the exact coefficient 
```

```
$ grep -F -m1 'reading/outer/transport_centre/component_1/labels_drain2/item4/condition_union: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -F -m3 'reading/outer/transport_centre/component_1/labels_{drain2}/item4/condition: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
reading/outer/transport_centre/component_1/labels_drain2/item4/condition_union: Eq(A_2(t), 0) | Eq(B_2(t), 0) | Eq(Derivative(X_2_1(t), t), 0)
reading/outer/transport_centre/component_1/labels_{drain2}/item4/condition: Derivative[1][xx[2, 1]][t] == 0 && rout[t] > 0 && 0 < rin[t] && rin[t] < rout[t]
reading/outer/transport_centre/component_1/labels_{drain2}/item4/condition: AA[2][t] == 0 && rout[t] > 0 && 0 < rin[t] && rin[t] < rout[t]
reading/outer/transport_centre/component_1/labels_{drain2}/item4/condition: BB[2][t] == 0 && rout[t] > 0 && 0 < rin[t] && rin[t] < rout[t]
```

```
$ for f in research/pde_ledger_v3/scripts/out/S15_separation_departure.out research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; do echo "$f: NOT OBTAINED $(grep -c 'NOT OBTAINED' $f); EXACT_OBJECT_GAP $(grep -c 'EXACT_OBJECT_GAP' $f)"; done
research/pde_ledger_v3/scripts/out/S15_separation_departure.out: NOT OBTAINED 0; EXACT_OBJECT_GAP 0
research/pde_ledger_v3/mathematica/out/S15_separation_departure.out: NOT OBTAINED 0; EXACT_OBJECT_GAP 0
```

## M4. The legs' evidence on item 4 (residuals and existence counts only)

```
$ cat /var/projects/_review_scratch/s15-reading-r2-claude/x06_item4_condition_lookup.stdout | cut -c1-260
PY: objects 48; with an object-level item4/condition line 48; containing a Sum jet equation 48; containing a literal exponent-only equality 0; condition line length min/max 525/2676; Q_under_condition lines per object [1]
PY: distinct item4/status texts: ['Finite parameter equations computed from the closed-form differential system and unit-sphere coefficient reduction. Q under the condition is computed from the Euler equation and its symbolic reference value. Exponent branches
WL: objects 48; with an object-level item4/condition line 48; containing a Sum jet equation 48; containing a literal exponent-only equality 0; condition line length min/max 522/4016; Q_under_condition lines per object [1]
WL: distinct item4/status texts: ['"Finite parameter equations computed from the closed-form differential system and exact direction coefficient elimination. Q under that condition is computed from the Euler equation with its symbolic reference value. Branches
```

```
$ grep -E '^## ' /var/projects/_review_scratch/s15-reading-r2-claude/e09_condition_py.stdout | cut -c1-120; grep -E '^OFF .*residual/scale=' /var/projects/_review_scratch/s15-reading-r2-claude/e09_condition_py.stdout | grep -o -E 'residual/scale=[^ ]+' | sort -t= -k2 -g | tail -3; echo "OFF residual lines: $(grep -c -E '^OFF .*residual/scale=' /var/projects/_review_scratch/s15-reading-r2-claude/e09_condition_py.stdout)"; grep -E '^OFF 2a1\+b1=37/97 scalar \(1, 0, 0\)(: S| n=0)' /var/projects/_review_scratch/s15-reading-r2-claude/e09_condition_py.stdout | head -4 | cut -c1-200; echo 'the scalar whose OFF residual/scale reaches 1.0:'; sed -n '/^OFF 2a1+b1=37\/97 scalar (1, 2, 0)/p' /var/projects/_review_scratch/s15-reading-r2-claude/e09_condition_py.stdout | cut -c1-200
## scalar (1, 0, 0): bound B=83, state size=42, reference 2*R(t)
## scalar (1, 2, 0): bound B=27, state size=14, reference 2*R(t)
## scalar (1, 0, 2): bound B=27, state size=14, reference 2*R(t)
residual/scale=1.0
residual/scale=1.0
residual/scale=1.0
OFF residual lines: 42
OFF 2a1+b1=37/97 scalar (1, 0, 0) n=0: c_n(transfer)=-1581.7921834966667612 t_n(direct Taylor of H*f)=-1581.7921834966667612 residual/scale=5.6272e-42
OFF 2a1+b1=37/97 scalar (1, 0, 0): S = sum_(n<=B) |d0^n c_n|^2 = 7.95821378079735e+21;  max_n |d0^n c_n| = 6.28673444e+10
the scalar whose OFF residual/scale reaches 1.0:
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=0: c_n(transfer)=0 t_n(direct Taylor of H*f)=4.7483553265588286407e-47 residual/scale=0.5098
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=1: c_n(transfer)=0 t_n(direct Taylor of H*f)=-3.3257799783261202507e-47 residual/scale=0.35707
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=2: c_n(transfer)=0 t_n(direct Taylor of H*f)=-9.3140816020961620654e-47 residual/scale=1.0
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=3: c_n(transfer)=0 t_n(direct Taylor of H*f)=-5.9823959180588496127e-47 residual/scale=0.6423
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=4: c_n(transfer)=0 t_n(direct Taylor of H*f)=-2.9919873529170464813e-47 residual/scale=0.32123
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=5: c_n(transfer)=0 t_n(direct Taylor of H*f)=-1.7383880789006704905e-47 residual/scale=0.18664
OFF 2a1+b1=37/97 scalar (1, 2, 0) n=6: c_n(transfer)=0 t_n(direct Taylor of H*f)=-5.1431076979257552148e-48 residual/scale=0.055219
OFF 2a1+b1=37/97 scalar (1, 2, 0): S = sum_(n<=B) |d0^n c_n|^2 = 0.0;  max_n |d0^n c_n| = 0.0
```

```
$ for f in e10d_condition_wl e10e_condition_wl; do echo "$f: lines $(grep -c . /var/projects/_review_scratch/s15-reading-r2-claude/$f.stdout); OFF residual lines $(grep -c -E '^OFF .*residual/scale=' /var/projects/_review_scratch/s15-reading-r2-claude/$f.stdout)"; grep -E '^OFF .*residual/scale=' /var/projects/_review_scratch/s15-reading-r2-claude/$f.stdout | head -2 | cut -c1-230; done; for j in claude-e09 claude-e10d claude-e10e; do echo "$j: $(tr -d ' \n' < _scratch/s11c/s15-reading-review-r2-guard/$j/outcome.json | cut -c1-60)"; done
e10d_condition_wl: lines 274; OFF residual lines 42
OFF 2a1+b1=37/97 n=0: c_n(transfer)=NumberForm[0, 20]  t_n(direct Taylor of H*f)=NumberForm[0``50.671949383027574, 20]  residual/scale=ScientificForm[0``-249.32805061697243]
OFF 2a1+b1=37/97 n=1: c_n(transfer)=NumberForm[0, 20]  t_n(direct Taylor of H*f)=NumberForm[0``50.23386771769388, 20]  residual/scale=ScientificForm[0``-249.76613228230616]
e10e_condition_wl: lines 92; OFF residual lines 14
OFF 2a1+b1=37/97 n=0: c_n(transfer)=NumberForm[443639.72387814900012220583836845990695237685`20., 20]  t_n(direct Taylor of H*f)=NumberForm[443639.72387814900012220583836845990695237685`20., 20]  residual/scale=ScientificForm[0``2
OFF 2a1+b1=37/97 n=1: c_n(transfer)=NumberForm[603884.83336592191641162987189201906424762038`20., 20]  t_n(direct Taylor of H*f)=NumberForm[603884.83336592191641162987189201906424762038`20., 20]  residual/scale=ScientificForm[0``2
claude-e09: {"exitCode":0,"wallSeconds":40.951067727059126,"unit":"s11c-
claude-e10d: {"exitCode":0,"wallSeconds":59.30550533905625,"unit":"s11c-g
claude-e10e: {"exitCode":0,"wallSeconds":24.54095304198563,"unit":"s11c-g
```

```
$ ls -la /var/projects/_review_scratch/s15-reading-r2-claude/d01_derive.py /var/projects/_review_scratch/s15-reading-r2-claude/d01c_momentum_rate.py /var/projects/_review_scratch/s15-reading-r2-claude/d02b_quadrature.py /var/projects/_review_scratch/s15-reading-r2-claude/d01_derive.stdout /var/projects/_review_scratch/s15-reading-r2-claude/d01c_momentum_rate.stdout /var/projects/_review_scratch/s15-reading-r2-claude/d02b_quadrature.stdout | awk '{print $5, $9}'; grep -E '^MAX ' /var/projects/_review_scratch/s15-reading-r2-claude/d02b_quadrature.stdout | cut -c1-140; sed -n '1p' /var/projects/_review_scratch/s15-reading-r2-claude/d01c_momentum_rate.stdout; grep -E '/not_obtained: ' /var/projects/_review_scratch/s15-reading-r2-claude/d01c_momentum_rate.stdout
21975 /var/projects/_review_scratch/s15-reading-r2-claude/d01_derive.py
18106 /var/projects/_review_scratch/s15-reading-r2-claude/d01_derive.stdout
8612 /var/projects/_review_scratch/s15-reading-r2-claude/d01c_momentum_rate.py
21904 /var/projects/_review_scratch/s15-reading-r2-claude/d01c_momentum_rate.stdout
4308 /var/projects/_review_scratch/s15-reading-r2-claude/d02b_quadrature.py
28231 /var/projects/_review_scratch/s15-reading-r2-claude/d02b_quadrature.stdout
MAX O1: scaled residual 2.859e-07; route-2 convergence estimate 8.343e-04
MAX O3: scaled residual 2.259e-08; route-2 convergence estimate 1.934e-04
MAX O4: scaled residual 4.907e-15; route-2 convergence estimate 8.059e-14
MAX O6: scaled residual 4.248e-09; route-2 convergence estimate 2.586e-05
MAX O7: scaled residual 4.343e-09; route-2 convergence estimate 1.952e-04
# integer strata -5..3; per-stratum time limit 900 s
O6O7/q=-5/not_obtained: closed form or identity solve not finished within 900 s
```

```
$ ls -la /var/projects/_review_scratch/s15-reading-r2-grok/derive_item34.py /var/projects/_review_scratch/s15-reading-r2-grok/derive_special.py /var/projects/_review_scratch/s15-reading-r2-grok/derive_conditions.py /var/projects/_review_scratch/s15-reading-r2-grok/derive_item34.stdout /var/projects/_review_scratch/s15-reading-r2-grok/derive_special.stdout /var/projects/_review_scratch/s15-reading-r2-grok/derive_conditions.partial.stdout | awk '{print $5, $9}'; grep -h -o -E 'abs=[0-9.e+-]+' /var/projects/_review_scratch/s15-reading-r2-grok/derive_special.stdout /var/projects/_review_scratch/s15-reading-r2-grok/derive_conditions.partial.stdout | sort | uniq -c; for j in grok-derive-01 grok-derive-02 grok-derive-03 grok-derive-04; do if [ -f _scratch/s11c/s15-reading-review-r2-guard/$j/outcome.json ]; then echo "$j: $(tr -d ' \n' < _scratch/s11c/s15-reading-review-r2-guard/$j/outcome.json | cut -c1-80)"; else echo "$j: no outcome.json; stdout lines $(wc -l < _scratch/s11c/s15-reading-review-r2-guard/$j/stdout)"; fi; done
14493 /var/projects/_review_scratch/s15-reading-r2-grok/derive_conditions.partial.stdout
9472 /var/projects/_review_scratch/s15-reading-r2-grok/derive_conditions.py
48025 /var/projects/_review_scratch/s15-reading-r2-grok/derive_item34.py
101809 /var/projects/_review_scratch/s15-reading-r2-grok/derive_item34.stdout
6376 /var/projects/_review_scratch/s15-reading-r2-grok/derive_special.py
7119 /var/projects/_review_scratch/s15-reading-r2-grok/derive_special.stdout
     16 abs=0.0
      1 abs=2.9404421517044706e-33
grok-derive-01: {"exitCode":0,"wallSeconds":397.7550922851078,"unit":"s11c-guard-ae4dd8a0c51d","
grok-derive-02: {"exitCode":0,"wallSeconds":17.946922942996025,"unit":"s11c-guard-0d437f280dec",
grok-derive-03: {"exitCode":1,"wallSeconds":1.3392503350041807,"unit":"s11c-guard-d2b8f474116f",
grok-derive-04: no outcome.json; stdout lines 188
```

```
$ wc -l /var/projects/_review_scratch/s15-reading-r2-grok/sympy_item4_conditions.txt /var/projects/_review_scratch/s15-reading-r2-grok/wolfram_quote_lines.txt; grep -c -F 'Eq(Sum(' /var/projects/_review_scratch/s15-reading-r2-grok/sympy_item4_conditions.txt; grep -c -F 'Inactive[Sum]' /var/projects/_review_scratch/s15-reading-r2-grok/wolfram_quote_lines.txt
   24 /var/projects/_review_scratch/s15-reading-r2-grok/sympy_item4_conditions.txt
    3 /var/projects/_review_scratch/s15-reading-r2-grok/wolfram_quote_lines.txt
   27 total
6
2
```

## M5. The claims

```
$ sed -n '15,18p;30,33p;39p' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md
- Inventory: 336 Q sections; 336 item-3 sections; 336 item-4 statuses; 48 separation-dependent objects with finite-identity conditions. Transcript end marker present: `True`.
- Finite-identity objects missing a condition: `[]`; missing restricted Q: `[]`. Remaining separation occurrences in the printed condition definitions: `[]`.

Requested-object gaps indexed by this final engine: none. General-field family integrals remain unevaluated where item 3 expressly permits them.
- Inventory: 336 Q sections; 336 item-3 sections; 336 item-4 statuses; 48 separation-dependent objects with finite-identity conditions. Transcript end marker present: `True`.
- Finite-identity objects missing a condition: `[]`; missing restricted Q: `[]`. Remaining separation occurrences in the printed condition definitions: `[]`.

Requested-object gaps indexed by this final engine: none. General-field family integrals remain unevaluated where item 3 expressly permits them.
Each engine solves its computed geometric scalar equation and inserts its closed Q at a regular symbolic reference, with the same instantaneous coefficients, to obtain Q under each condition. The family responses remain live.
```

```
$ grep -F -m1 'reading/outer/carried/component_1/labels_drain1/item4/status: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-400; grep -F -m1 'reading/outer/carried/component_1/labels_{drain1}/item4/status: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-400
reading/outer/carried/component_1/labels_drain1/item4/status: Finite parameter equations computed from the closed-form differential system and unit-sphere coefficient reduction. Q under the condition is computed from the Euler equation and its symbolic reference value. Exponent branches remain in the coefficient and initial-state definitions.
reading/outer/carried/component_1/labels_{drain1}/item4/status: "Finite parameter equations computed from the closed-form differential system and exact direction coefficient elimination. Q under that condition is computed from the Euler equation with its symbolic reference value. Branches are retained in the printed coefficient and initial-state definitions."
```

## M6. W2, and items 1, 2, 5 and 6 against round 1

```
$ grep -E '^SUMMARY' /var/projects/_review_scratch/s15-reading-r2-claude/e08_compare.stdout; echo "rows beyond 10x reference convergence: $(grep -c 'beyond 10x reference convergence' /var/projects/_review_scratch/s15-reading-r2-claude/e08_compare.stdout)"
SUMMARY PY: max scaled residual over evaluable values 1.361e-03; not evaluable/nonfinite 0
SUMMARY WL: max scaled residual over evaluable values 1.361e-03; not evaluable/nonfinite 0
SUMMARY WL_r1: max scaled residual over evaluable values 1.890e-07; not evaluable/nonfinite 30
rows beyond 10x reference convergence: 0
```

```
$ grep -n -o '^reading/region/component_[123]/rate_labels_{drain1}/item3/[A-Za-z_]*' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
57255:reading/region/component_1/rate_labels_{drain1}/item3/family_integral
57256:reading/region/component_1/rate_labels_{drain1}/item3/family_integral
57257:reading/region/component_1/rate_labels_{drain1}/item3/existence_condition
57259:reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand
57260:reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand
57262:reading/region/component_1/rate_labels_{drain1}/item3/object
57263:reading/region/component_1/rate_labels_{drain1}/item3/object
57264:reading/region/component_1/rate_labels_{drain1}/item3/exponent_branch_inventory
57265:reading/region/component_1/rate_labels_{drain1}/item3/retained_Appell_parameters
57266:reading/region/component_1/rate_labels_{drain1}/item3/exponent_divisor_supports
57267:reading/region/component_1/rate_labels_{drain1}/item3/exponent_divisor_supports
57268:reading/region/component_1/rate_labels_{drain1}/item3/exponent_divisor_loci
57269:reading/region/component_1/rate_labels_{drain1}/item3/exponent_divisor_loci
57274:reading/region/component_1/rate_labels_{drain1}/item3/status
62753:reading/region/component_2/rate_labels_{drain1}/item3/family_integral
62754:reading/region/component_2/rate_labels_{drain1}/item3/family_integral
62755:reading/region/component_2/rate_labels_{drain1}/item3/existence_condition
62757:reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand
62758:reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand
62760:reading/region/component_2/rate_labels_{drain1}/item3/object
62761:reading/region/component_2/rate_labels_{drain1}/item3/object
62762:reading/region/component_2/rate_labels_{drain1}/item3/exponent_branch_inventory
62763:reading/region/component_2/rate_labels_{drain1}/item3/retained_Appell_parameters
62764:reading/region/component_2/rate_labels_{drain1}/item3/exponent_divisor_supports
62765:reading/region/component_2/rate_labels_{drain1}/item3/exponent_divisor_supports
62766:reading/region/component_2/rate_labels_{drain1}/item3/exponent_divisor_loci
62767:reading/region/component_2/rate_labels_{drain1}/item3/exponent_divisor_loci
62772:reading/region/component_2/rate_labels_{drain1}/item3/status
67454:reading/region/component_3/rate_labels_{drain1}/item3/family_integral
67455:reading/region/component_3/rate_labels_{drain1}/item3/family_integral
67456:reading/region/component_3/rate_labels_{drain1}/item3/existence_condition
67458:reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand
67459:reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand
67461:reading/region/component_3/rate_labels_{drain1}/item3/object
67462:reading/region/component_3/rate_labels_{drain1}/item3/object
67463:reading/region/component_3/rate_labels_{drain1}/item3/exponent_branch_inventory
67464:reading/region/component_3/rate_labels_{drain1}/item3/retained_Appell_parameters
67465:reading/region/component_3/rate_labels_{drain1}/item3/exponent_divisor_supports
67466:reading/region/component_3/rate_labels_{drain1}/item3/exponent_divisor_supports
67467:reading/region/component_3/rate_labels_{drain1}/item3/exponent_divisor_loci
67468:reading/region/component_3/rate_labels_{drain1}/item3/exponent_divisor_loci
67473:reading/region/component_3/rate_labels_{drain1}/item3/status
```

```
$ for n in $(grep -n -E '^reading/region/component_[123]/rate_labels_\{drain1\}/item3/(object|Q_closed_operand): ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -d: -f1); do L=$(sed -n "${n}p" research/pde_ledger_v3/mathematica/out/S15_separation_departure.out); echo "line $n $(echo "$L" | cut -d: -f1)"; for f in '(2 + aa[1] + bb[1])' '(3 + aa[1] + bb[1])' '(4 + aa[1] + bb[1])' '(5 + aa[1] + bb[1])'; do echo "  literal $f: $(echo "$L" | grep -o -F "$f" | wc -l)"; done; done
line 57259 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand
  literal (2 + aa[1] + bb[1]): 0
  literal (3 + aa[1] + bb[1]): 0
  literal (4 + aa[1] + bb[1]): 0
  literal (5 + aa[1] + bb[1]): 0
line 57262 reading/region/component_1/rate_labels_{drain1}/item3/object
  literal (2 + aa[1] + bb[1]): 0
  literal (3 + aa[1] + bb[1]): 0
  literal (4 + aa[1] + bb[1]): 0
  literal (5 + aa[1] + bb[1]): 0
line 62757 reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand
  literal (2 + aa[1] + bb[1]): 0
  literal (3 + aa[1] + bb[1]): 0
  literal (4 + aa[1] + bb[1]): 0
  literal (5 + aa[1] + bb[1]): 0
line 62760 reading/region/component_2/rate_labels_{drain1}/item3/object
  literal (2 + aa[1] + bb[1]): 0
  literal (3 + aa[1] + bb[1]): 0
  literal (4 + aa[1] + bb[1]): 0
  literal (5 + aa[1] + bb[1]): 0
line 67458 reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand
  literal (2 + aa[1] + bb[1]): 0
  literal (3 + aa[1] + bb[1]): 0
  literal (4 + aa[1] + bb[1]): 0
  literal (5 + aa[1] + bb[1]): 0
line 67461 reading/region/component_3/rate_labels_{drain1}/item3/object
  literal (2 + aa[1] + bb[1]): 0
  literal (3 + aa[1] + bb[1]): 0
  literal (4 + aa[1] + bb[1]): 0
  literal (5 + aa[1] + bb[1]): 0
```

```
$ grep -E '^reading/region/component_[123]/rate_labels_\{drain1\}/item3/status: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-300
reading/region/component_1/rate_labels_{drain1}/item3/status: "Closed-form return from the angular, distance and finite-interval integration construction. The exponent-branch and retained-function inventories are printed separately."
reading/region/component_2/rate_labels_{drain1}/item3/status: "Closed-form return from the angular, distance and finite-interval integration construction. The exponent-branch and retained-function inventories are printed separately."
reading/region/component_3/rate_labels_{drain1}/item3/status: "Closed-form return from the angular, distance and finite-interval integration construction. The exponent-branch and retained-function inventories are printed separately."
```

```
$ cat /var/projects/_review_scratch/s15-reading-r2-claude/x02_PY.stdout | cut -c1-160; cat /var/projects/_review_scratch/s15-reading-r2-claude/x02_WL.stdout | cut -c1-160
PY ('construction/item1', 'only_r1') 1
PY ('construction/item1', 'same') 835
PY ('item1(reading Q)', 'same') 672
PY ('item2', 'differ') 12
PY ('item2', 'same') 3348
PY ('item3', 'differ') 150
PY ('item3', 'only_r2') 336
PY ('item3', 'same') 2154
PY ('item3_kernel', 'same') 67
PY ('item4', 'differ') 48
PY ('item4', 'only_r1') 432
PY ('item4', 'only_r2') 8712
PY ('item4', 'same') 5280
PY ('item4_kernel', 'only_r2') 28
PY ('item5', 'same') 4704
PY ('item6', 'same') 3246
PY ('construction/item1', 'only_r1') EXACT_OBJECT_GAP
PY ('item2', 'differ') reading/inner/pressure/component_1/whole_minus_expansion/E_geo
PY ('item2', 'differ') reading/inner/pressure/component_1/whole_minus_expansion/E_resp
PY ('item2', 'differ') reading/inner/pressure/component_2/whole_minus_expansion/E_geo
PY ('item2', 'differ') reading/inner/pressure/component_2/whole_minus_expansion/E_resp
PY ('item2', 'differ') reading/inner/pressure/component_3/whole_minus_expansion/E_geo
PY ('item2', 'differ') reading/inner/pressure/component_3/whole_minus_expansion/E_resp
PY ('item2', 'differ') reading/outer/pressure/component_1/whole_minus_expansion/E_geo
PY ('item2', 'differ') reading/outer/pressure/component_1/whole_minus_expansion/E_resp
PY ('item2', 'differ') reading/outer/pressure/component_2/whole_minus_expansion/E_geo
PY ('item2', 'differ') reading/outer/pressure/component_2/whole_minus_expansion/E_resp
PY ('item2', 'differ') reading/outer/pressure/component_3/whole_minus_expansion/E_geo
PY ('item2', 'differ') reading/outer/pressure/component_3/whole_minus_expansion/E_resp
WL ('construction/item1', 'only_r1') 1
WL ('construction/item1', 'same') 811
WL ('item1(reading Q)', 'same') 672
WL ('item2', 'same') 3360
WL ('item3', 'differ') 168
WL ('item3', 'only_r2') 1824
WL ('item3', 'same') 2136
WL ('item3_kernel', 'only_r2') 15
WL ('item3_kernel', 'same') 51
WL ('item4', 'differ') 330
WL ('item4', 'only_r1') 1680
WL ('item4', 'only_r2') 25578
WL ('item4', 'same') 4518
WL ('item4_kernel', 'only_r2') 31
WL ('item4_kernel', 'same') 14
WL ('item5', 'same') 5376
WL ('item6', 'same') 3242
WL ('construction/item1', 'only_r1') EXACT_OBJECT_GAP
```

```
$ tail -5 /var/projects/_review_scratch/s15-reading-r2-claude/x04_py_value_compare.guarded.child.stdout | cut -c1-200; cmp /var/projects/_review_scratch/s15-reading-r2-claude/x04_py_value_compare.unguarded.stdout /var/projects/_review_scratch/s15-reading-r2-claude/x04_py_value_compare.guarded.child.stdout && echo 'x04 unguarded and guarded stdout identical'
reading/outer/pressure/component_2/whole_minus_expansion/E_geo: parsed-equal=False; expand(r1-r2) = 0
reading/outer/pressure/component_2/whole_minus_expansion/E_resp: parsed-equal=False; expand(r1-r2) = 0
reading/outer/pressure/component_3/whole_minus_expansion/E_geo: parsed-equal=False; expand(r1-r2) = 0
reading/outer/pressure/component_3/whole_minus_expansion/E_resp: parsed-equal=False; expand(r1-r2) = 0
END x04
x04 unguarded and guarded stdout identical
```

```
$ head -30 /var/projects/_review_scratch/s15-reading-r2-grok/diff_sympy.txt | cut -c1-160; echo ...; head -20 /var/projects/_review_scratch/s15-reading-r2-grok/diff_wolfram.txt | cut -c1-160
only_r1 433
only_r2 9076
differ 210
ONLY_R2_KINDS Counter({'item4': 8712, 'item3': 336, 'other': 28})
DIFFER_KINDS Counter({'item3': 150, 'item4': 48, 'item2': 12})
--- differ sample ---
item3 reading/inner/carried/component_1/labels_drain1+drain2/item3/status
item3 reading/inner/carried/component_1/labels_drain1/item3/status
item3 reading/inner/carried/component_1/labels_drain2/item3/status
item3 reading/inner/carried/component_2/labels_drain1+drain2/item3/status
item3 reading/inner/carried/component_2/labels_drain1/item3/status
item3 reading/inner/carried/component_2/labels_drain2/item3/status
item3 reading/inner/carried/component_3/labels_drain1+drain2/item3/status
item3 reading/inner/carried/component_3/labels_drain1/item3/status
item3 reading/inner/carried/component_3/labels_drain2/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_drain1+drain2/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_drain1+remainder/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_drain1/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_drain2/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_remainder/item3/status
item2 reading/inner/pressure/component_1/whole_minus_expansion/E_geo
item2 reading/inner/pressure/component_1/whole_minus_expansion/E_resp
item3 reading/inner/pressure/component_1/whole_minus_expansion/item3/family_integral
item3 reading/inner/pressure/component_1/whole_minus_expansion/item3/object
item3 reading/inner/pressure/component_2/expansion_labels_drain1+drain2/item3/status
item3 reading/inner/pressure/component_2/expansion_labels_drain1+remainder/item3/status
item3 reading/inner/pressure/component_2/expansion_labels_drain1/item3/status
item3 reading/inner/pressure/component_2/expansion_labels_drain2/item3/status
item3 reading/inner/pressure/component_2/expansion_labels_remainder/item3/status
item2 reading/inner/pressure/component_2/whole_minus_expansion/E_geo
...
only_r1 1681
only_r2 27448
differ 498
ONLY_R2_KINDS Counter({'item4': 25578, 'item3': 1824, 'other': 46})
DIFFER_KINDS Counter({'item4': 330, 'item3': 168})
--- differ sample ---
item3 reading/inner/carried/component_1/labels_{drain1, drain2}/item3/status
item3 reading/inner/carried/component_1/labels_{drain1}/item3/status
item3 reading/inner/carried/component_1/labels_{drain2}/item3/status
item3 reading/inner/carried/component_2/labels_{drain1, drain2}/item3/status
item3 reading/inner/carried/component_2/labels_{drain1}/item3/status
item3 reading/inner/carried/component_2/labels_{drain2}/item3/status
item3 reading/inner/carried/component_3/labels_{drain1, drain2}/item3/status
item3 reading/inner/carried/component_3/labels_{drain1}/item3/status
item3 reading/inner/carried/component_3/labels_{drain2}/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_{drain1, drain2}/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_{drain1, remainder}/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_{drain1}/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_{drain2}/item3/status
item3 reading/inner/pressure/component_1/expansion_labels_{remainder}/item3/status
```

```
$ cat /var/projects/_review_scratch/s15-reading-r2-grok/loci_sympy.txt; head -5 /var/projects/_review_scratch/s15-reading-r2-grok/loci_wolfram.txt | cut -c1-200
loci_lines 0
flagged 0
loci_lines 288
flagged 192
FLAG ('reading/inner/carried/component_1/labels_{drain2}/item3/exponent_divisor_loci', '{}\n')
FLAG ('reading/inner/carried/component_1/labels_{drain2}/Q_operand/item3/exponent_divisor_loci', '{}\n')
FLAG ('reading/inner/carried/component_2/labels_{drain2}/item3/exponent_divisor_loci', '{}\n')
```

## M7. Controls

```
$ for a in wl_formA py_formA py_famA wl_depB py_depC wl_balC; do echo "== $a"; head -3 /var/projects/_review_scratch/s15-reading-r2-claude/abl/$a.diff | cut -c1-200; done
== wl_formA
--- /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
+++ /var/projects/_review_scratch/s15-reading-r2-claude/abl/wl_formA/S15_separation_departure.wl
@@ -436,7 +436,7 @@
== py_formA
--- /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
+++ /var/projects/_review_scratch/s15-reading-r2-claude/abl/py_formA/S15_separation_departure.py
@@ -880,7 +880,7 @@
== py_famA
--- /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
+++ /var/projects/_review_scratch/s15-reading-r2-claude/abl/py_famA/S15_separation_departure.py
@@ -951,7 +951,7 @@
== wl_depB
--- /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
+++ /var/projects/_review_scratch/s15-reading-r2-claude/abl/wl_depB/S15_separation_departure.wl
@@ -799,7 +799,7 @@
== py_depC
--- /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
+++ /var/projects/_review_scratch/s15-reading-r2-claude/abl/py_depC/S15_separation_departure.py
@@ -1188,7 +1188,10 @@
== wl_balC
--- /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
+++ /var/projects/_review_scratch/s15-reading-r2-claude/abl/wl_balC/S15_separation_departure.wl
@@ -51,7 +51,7 @@
```

```
$ for s in WLFORMA PYFORMA PYFAMA WLDEPB PYDEPC WLBALC; do echo "== $s"; cut -c1-200 /var/projects/_review_scratch/s15-reading-r2-claude/x02_$s.stdout | head -16; done
== WLFORMA
WLFORMA ('construction/item1', 'same') 811
WLFORMA ('item1(reading Q)', 'same') 672
WLFORMA ('item2', 'same') 3360
WLFORMA ('item3', 'differ') 51
WLFORMA ('item3', 'same') 4077
WLFORMA ('item3_kernel', 'differ') 1
WLFORMA ('item3_kernel', 'same') 65
WLFORMA ('item4', 'differ') 2457
WLFORMA ('item4', 'same') 27969
WLFORMA ('item4_kernel', 'same') 45
WLFORMA ('item5', 'same') 5376
WLFORMA ('item6', 'same') 3242
== PYFORMA
PYFORMA ('construction/item1', 'same') 835
PYFORMA ('item1(reading Q)', 'same') 6
PYFORMA ('item2', 'same') 30
PYFORMA ('item3', 'differ') 6
PYFORMA ('item3', 'same') 21
PYFORMA ('item3_kernel', 'only_r1') 4
PYFORMA ('item3_kernel', 'same') 63
PYFORMA ('item4', 'differ') 93
PYFORMA ('item4', 'same') 217
PYFORMA ('item4_kernel', 'same') 28
PYFORMA ('item5', 'same') 42
PYFORMA ('item6', 'same') 582
== PYFAMA
PYFAMA ('construction/item1', 'only_r2') 1
PYFAMA ('construction/item1', 'same') 835
PYFAMA ('item1(reading Q)', 'same') 6
PYFAMA ('item2', 'same') 30
PYFAMA ('item3', 'differ') 3
PYFAMA ('item3', 'same') 24
PYFAMA ('item3_kernel', 'same') 67
PYFAMA ('item4', 'differ') 73
PYFAMA ('item4', 'only_r1') 1
PYFAMA ('item4', 'same') 236
PYFAMA ('item4_kernel', 'same') 28
PYFAMA ('item5', 'same') 42
PYFAMA ('item6', 'same') 582
PYFAMA ('construction/item1', 'only_r2') EXACT_OBJECT_GAP
== WLDEPB
WLDEPB ('construction/item1', 'same') 811
WLDEPB ('item1(reading Q)', 'same') 672
WLDEPB ('item2', 'same') 3360
WLDEPB ('item3', 'differ') 54
WLDEPB ('item3', 'same') 4074
WLDEPB ('item3_kernel', 'same') 66
WLDEPB ('item4', 'differ') 2865
WLDEPB ('item4', 'only_r1') 216
WLDEPB ('item4', 'only_r2') 1002
WLDEPB ('item4', 'same') 27345
WLDEPB ('item4_kernel', 'same') 45
WLDEPB ('item5', 'same') 5376
WLDEPB ('item6', 'same') 3242
== PYDEPC
PYDEPC ('construction/item1', 'same') 835
PYDEPC ('item1(reading Q)', 'same') 6
PYDEPC ('item2', 'same') 30
PYDEPC ('item3', 'differ') 3
PYDEPC ('item3', 'same') 24
PYDEPC ('item3_kernel', 'same') 67
PYDEPC ('item4', 'differ') 66
PYDEPC ('item4', 'same') 244
PYDEPC ('item4_kernel', 'same') 28
PYDEPC ('item5', 'same') 42
PYDEPC ('item6', 'same') 582
== WLBALC
WLBALC ('construction/item1', 'differ') 78
WLBALC ('construction/item1', 'same') 733
WLBALC ('item1(reading Q)', 'differ') 66
WLBALC ('item1(reading Q)', 'same') 606
WLBALC ('item2', 'differ') 330
WLBALC ('item2', 'same') 3030
WLBALC ('item3', 'differ') 138
WLBALC ('item3', 'same') 3990
WLBALC ('item3_kernel', 'differ') 5
WLBALC ('item3_kernel', 'same') 61
WLBALC ('item4', 'differ') 2100
WLBALC ('item4', 'same') 28326
WLBALC ('item4_kernel', 'same') 45
WLBALC ('item5', 'differ') 264
WLBALC ('item5', 'same') 5112
WLBALC ('item6', 'same') 3242
```

```
$ for s in WLFORMA WLDEPB; do echo "== $s"; cut -c1-200 /var/projects/_review_scratch/s15-reading-r2-claude/x09_$s.stdout | head -8; done
== WLFORMA
item4|differ: tags 2457; identical after renumbering generated symbols 1566; still different 891
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/direction_coefficient_operands', 'content')
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/scalar_1/decomposition', 'content')
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/scalar_1/output_row', 'content')
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/scalar_1/polynomial_output', 'content')
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/scalar_1/projection_definition', 'content')
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/scalar_4/decomposition', 'content')
   ('reading/outer/carried/component_1/labels_{drain1, drain2}/item4/scalar_4/output_row', 'content')
== WLDEPB
item4|differ: tags 2865; identical after renumbering generated symbols 2091; still different 774
objects with real item4 differences: 54
   reading/inner/carried/component_1/labels_{drain1, drain2} 14
   reading/inner/carried/component_1/labels_{drain1} 14
   reading/inner/carried/component_2/labels_{drain1, drain2} 14
   reading/inner/carried/component_2/labels_{drain1} 14
   reading/inner/carried/component_3/labels_{drain1, drain2} 14
   reading/inner/carried/component_3/labels_{drain1} 14
```

## M8. Reproducibility

```
$ sha256sum /var/projects/_review_scratch/s15-reading-r2-claude/rerun/py/out/S15_separation_departure.out; ls /var/projects/_review_scratch/s15-reading-r2-claude/rerun; for j in claude-rerun_py2 claude-rerun_wl claude-x11b; do echo "$j: $(tr -d ' \n' < _scratch/s11c/s15-reading-review-r2-guard/$j/outcome.json | cut -c1-80)"; done
93d3eba10711abc52e36751a429cd0fc5b71738ff8b7c72095df319d1d2b531d  /var/projects/_review_scratch/s15-reading-r2-claude/rerun/py/out/S15_separation_departure.out
py
rerun_py.done
rerun_py.stdout
rerun_py.stdout.guard.json
rerun_py.stdout.guard_stderr
rerun_py.stdout.retries
rerun_py.stdout.stderr
rerun_wl.stdout
rerun_wl.stdout.guard.json
rerun_wl.stdout.guard_stderr
rerun_wl.stdout.retries
rerun_wl.stdout.stderr
wl
claude-rerun_py2: {"exitCode":0,"wallSeconds":12708.204215775011,"unit":"s11c-guard-b228c2bec280",
claude-rerun_wl: {"exitCode":0,"wallSeconds":1513.448981732945,"unit":"s11c-guard-d5322bd3820e","
claude-x11b: {"exitCode":0,"wallSeconds":378.28763447806705,"unit":"s11c-guard-3eff52d25ef3",
```

```
$ cat /var/projects/_review_scratch/s15-reading-r2-claude/x02_RERUNPY.stdout | cut -c1-200; cat /var/projects/_review_scratch/s15-reading-r2-claude/x09_RERUNWL.stdout | cut -c1-200
RERUNPY ('construction/item1', 'same') 835
RERUNPY ('item1(reading Q)', 'same') 672
RERUNPY ('item2', 'differ') 6
RERUNPY ('item2', 'same') 3354
RERUNPY ('item3', 'differ') 15
RERUNPY ('item3', 'same') 2625
RERUNPY ('item3_kernel', 'same') 67
RERUNPY ('item4', 'same') 14040
RERUNPY ('item4_kernel', 'same') 28
RERUNPY ('item5', 'same') 4704
RERUNPY ('item6', 'same') 3246
RERUNPY ('item2', 'differ') reading/outer/pressure/component_1/whole_minus_expansion/E_geo
RERUNPY ('item2', 'differ') reading/outer/pressure/component_1/whole_minus_expansion/E_resp
RERUNPY ('item2', 'differ') reading/outer/pressure/component_2/whole_minus_expansion/E_geo
RERUNPY ('item2', 'differ') reading/outer/pressure/component_2/whole_minus_expansion/E_resp
RERUNPY ('item2', 'differ') reading/outer/pressure/component_3/whole_minus_expansion/E_geo
RERUNPY ('item2', 'differ') reading/outer/pressure/component_3/whole_minus_expansion/E_resp
item4|differ: tags 1848; identical after renumbering generated symbols 1848; still different 0
```

```
$ echo "parsed-equal=True: $(grep -c 'parsed-equal=True' /var/projects/_review_scratch/s15-reading-r2-claude/x11_py_rerun_value_compare.child.stdout)"; echo "route1 zero after polylog(0,x)->x/(1-x): $(grep -c 'route1 expand(prod - rerun) after polylog(0,x)->x/(1-x) = 0$' /var/projects/_review_scratch/s15-reading-r2-claude/x11_py_rerun_value_compare.child.stdout)"; echo "route2 zero: $(grep -c 'route2 .*prod - rerun = 0$' /var/projects/_review_scratch/s15-reading-r2-claude/x11_py_rerun_value_compare.child.stdout)"
parsed-equal=True: 12
route1 zero after polylog(0,x)->x/(1-x): 9
route2 zero: 9
```
