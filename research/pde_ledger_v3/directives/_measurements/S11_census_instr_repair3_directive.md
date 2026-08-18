# Measurements — S11 census-instrument repair round 3 directive

Orchestrator verification commands behind the directive's five defect claims, run 2026-08-17
against the round-2 census outputs (`~/.s11_build/census_build2/`, sha256s in
`artifact_sha256.stdout` there) and instruments at `89ed80c9`. Leg reports:
`~/.s11_build/census_repair_grok_leg.txt` and the fresh-agent report (scripts under
`/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/review/`).

## Defect 1 (union poisoning) — grok finding 1, orchestrator-confirmed

```
$ sed -n '34p' ~/.s11_build/census_build2/containment_py.stdout | cut -c1-700
CONTAINMENT_COMPLETENESS ... tag=PY_S11_MAIN_D2_ROOT1_STACKED_DROP_JOINT_SOLUTION ...
missing_memberships=['"((Symbol(\'k1\', real=True), Mul(I, Symbol(\'k2\', real=True))),)"', ...
$ python3 - <<'EOF'
import sympy as sp
k1,k2=sp.symbols('k1 k2',real=True)
print(sp.simplify((k2 + sp.I*k1).subs(k1, sp.I*k2)))   # emitted {k2:-I k1} under candidate {k1:I k2}
print(sp.simplify((k2 - sp.I*k1).subs(k1, -sp.I*k2)))  # emitted {k2:+I k1} under candidate {k1:-I k2}
EOF
0
0
```

Both "missing" candidates satisfy an emitted branch identically ⇒ covered, verdicted OMITTED.

## Defect 2 (pole candidates) — grok finding 2, orchestrator-confirmed

```
$ for L in 1132 1192 1926 2001; do sed -n "${L}p" ~/.s11_build/census_build2/containment_wl.stdout \
    | grep -o 'tag=[A-Z0-9_]*\|verdict=[A-Z_]*'; done
tag=WL_S11_XFORM_EXTRA_D2_ROOT1_KW_ZERO_LOCUS_SOLUTION ... verdict=OMITTED_BRANCH   (and 3 more)
$ grep -m1 "WL_S11_XFORM_EXTRA_D2_ROOT1_KW_ZERO_LOCUS_EQUATIONS" \
    research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out | cut -c1-300
... (bComp*k1^2*rhoBr + ... - Sqrt[...])/(2*cs0^2*rhoBr^2) == 0}
```

The as-written equation is undefined at the candidate `rhoBr -> 0` (denominator `2 cs0^2 rhoBr^2`).

## Defect 3 (sampler collision) — agent finding 1, orchestrator-reproduced

```
$ python3  # simulate s11_census_math.py:827 assignment on (bComp,k1,k2,k3,muR), 4 serials
serial=0 bComp=1/2 muR=1/2 equal=True (bComp-muR)|sample=0
serial=1 bComp=1  muR=1  equal=True (bComp-muR)|sample=0
serial=2 bComp=2  muR=2  equal=True (bComp-muR)|sample=0
serial=3 bComp=3  muR=3  equal=True (bComp-muR)|sample=0
```

Positive pool length 4; `bComp`/`muR` at sorted-index distance 4 collide in every sample.

## Defect 4 (witness validated extinct) — agent finding 2, orchestrator-confirmed

```
$ grep -c "verdict=WITNESS_VALIDATED" ~/.s11_build/census_build2/containment_wl.stdout \
    ~/.s11_build/census_build2/containment_py.stdout
containment_wl.stdout:0
containment_py.stdout:0
$ grep -c "verdict=WITNESS_VALIDATED" ~/.s11_build/census_build/containment_wl.stdout \
    ~/.s11_build/census_build/containment_py.stdout
containment_wl.stdout:107
containment_py.stdout:99
```

Verdict split showing the failure direction is separate and intact
(failures fire on membership `SPURIOUS_BRANCH*` or premise FALSE, `s11_containment_census.py:302`):

```
$ grep "verdict=WITNESS_FAILURE" ~/.s11_build/census_build2/containment_wl.stdout \
    | grep -o 'premise_truth="[A-Z_]*"' | sort | uniq -c
     74 premise_truth="FALSE"
     98 premise_truth="UNDECIDED"
```

## Defect 5 (reducer sheet-line double count) — both legs, orchestrator-reproduced

```
$ python3  # count verdict-token lines by line kind across the four round-2 census stdouts
total limitation-token lines: 885
-- BRANCH_MEMBERSHIP_UNDECIDED by line kind --
  CONTAINMENT_BRANCH: 38
  CONTAINMENT_BRANCH_SHEET: 68
  CONTAINMENT_WITNESS_SHEET: 10
```

Object-level count 807; 78 sheet-progress lines double-count their parent objects.
`failures=889` (= 561 + 328 decided) and `findings=358` (= 74+29 / 13+70 / 172+0) reproduce
at object level — both legs and orchestrator agree.
