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

---

# Post-leg fold measurements (directive legs: Codex + Grok, 2026-08-17)

Leg reports: `~/.s11_build/repair3_directive_codex_leg.txt` (scripts `/tmp/s11_repair3_review/`),
`~/.s11_build/repair3_directive_grok_leg.txt` (scripts `/tmp/s11_r3_review/`).

## Defect 1 fold (piecewise union — Codex blocker, orchestrator-confirmed)

Candidate `{k1: I*Abs(k2)}` from `PY_S11_MAIN_D2_ROOT1_STACKED_DROP_JOINT` (Codex
`actual_union_piecewise_counterexample.stdout`: as-written equation statuses all ZERO):

```
$ python3  # branch constraints under the candidate, k2 real
branch A residual: k2 + Abs(k2)   (0 at k2=-3, 6 at k2=3)
branch B residual: k2 - Abs(k2)   (-6 at k2=-3, 0 at k2=3)
union covers all k2: True
neither branch alone covers: True True
```

Per-branch OR is refuted; the defined-branch union product (A·B = k2²−|k2|² ≡ 0) decides.
Grok's leg endorsed per-branch OR — refuted by this counterexample; Codex wins on computation.
Identical-cover population per Codex's stricter full-constraint-set predicate: 2 records /
5 candidates (grok's 7/70 counted single-constraint matches).

## Defect 2 fold (pole surface)

Leg counts disagree on the surface (grok 10/13, codex 12/13 records; codex
`reproduce_counts.stdout` shows per-record `('nan','UNDEFINED')` rows incl. D3/D4
ROOT_COINCIDENCE and STRATUM5–10; grok holds D2/D3/D4 coincidence records retain defined
missing candidates). Orchestrator-verified floor: the four `{rhoBr: 0}` KW_ZERO records and six
STRATUM `{bComp: muR*sRho}` single-candidate records (candidate_count=1 each, checked against
`containment_wl.stdout`). The directive therefore states the RULE, not a count; the builder's
per-record exclusion lines settle the number.

## Defect 3 fold (sampler population — both legs)

```
$ grep -c "UNDECIDED_ZERO_SAMPLES" ~/.s11_build/census_build2/probe_wl.stdout
81
```

Codex `sampler_collision_audit.stdout`: ≥51 records explicitly refuted by one generic
`bComp≠muR` point; 30 payloads truncated (unexamined). "28" (previous draft) was the round-1→2
transition subset, not the affected population.

## Defect 4/6 fold (witness semantics — both legs structurally, orchestrator on soundness)

Both legs: partition by PRE-substitution binding is wrong; substitute-then-classify is required.
Both legs took the instrument's `premise_truth=FALSE` on
`WL_S11_XKIN_ANISO_D2_STRATUM{5,6}_ROOT_COINCIDENCE_COEFF_REAL_WITNESS` at face value.
Orchestrator refutation — the atom is identically TRUE at the witness:

```
$ python3  # direct principal-value arithmetic at {bComp:1, muR:2, rhoBr:1, sRho:2}
A inner value: -81*k2**2        # strict-inequality atom -> contingent (k2**2 > 0)
B sqrt-arg after subst: 0       # Element[Sqrt[0], Reals] -> TRUE
$ python3  # instrument path: parse_payload(..., assumption_free=True) then subs
parsed pre-subst: Eq(((-re(k2)**2 + im(k2)**2)*(...)...   # re/im expansion, pre-substitution
post-subst: False                                          # branch-unsound
```

Same on STRATUM6 (atom-by-atom: TRUE,TRUE,TRUE,UNDECIDED,FALSE*,UNDECIDED,UNDECIDED,TRUE — the
FALSE* is the unsound expansion). So ≥2 of the 172 round-2 witness failures are instrument
artifacts; the round-2 "172 unchanged" constraint was dropped from the directive.
Witness binding scan (grok `06_witnesses.stdout`, `08_witness_payload_scan.stdout`): 307/480 WL
and 38/99 PY proved-nonempty witnesses bind some `k*` — the earlier draft's "binds only solve
variables" evidence sentence was false and is removed.
PY concrete-atom gap (grok `09_validate_plant_and_py_truth.stdout`):
`Q.real(k1) |-> Q.real(1) truth=UNDECIDED` — the directive now orders concrete-number atoms to
decide.

## Defect 5 (unchanged) — both legs reproduce 885 = 807 + 78; failures 889, findings 358 exact.
