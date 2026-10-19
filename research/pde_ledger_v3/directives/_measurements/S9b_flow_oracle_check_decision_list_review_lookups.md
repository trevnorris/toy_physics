# Measurements: review of the S9b flow-only oracle-check decision list, pass 1 (generated 2026-10-19 07:57)

Generator: `_scratch/s9b_oracle/gen/dl_review_lookups.sh` (sha256sum/grep/sed/cut/wc/git only). The legs'
evidence is hard-linked at `_scratch/s9b_oracle/dl_review_{codex,grok}_evidence/`.

## The reviewed list and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S9b_flow_oracle_check_decision_list_review_baseline.sha256)
research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md: OK
```

```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S9b_flow_oracle_check_decision_list_review_prompt.md | cut -c1-64
baaf2cf01a20728aff62bba1829c44ca8ff51342d92f7dd1a1744015957356c7
```

## Verdicts
```
$ grep -o 'Verdict:[^A-Z]*[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_codex_evidence/review.md | head -1
Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict:[^A-Z]*[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_grok.txt | head -1
Verdict:** NEEDS REVISION
```

## Both legs: the WL flow-only exchange field
```
$ sed -n '28,29p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md
  components and `IMPLIED_JN`; and WL's `GeneralSquaredVelocity`, `SquaredVelocity`, constant domains, cases and
  `Exchange`. The record `steps/S9b_light_bending.md` at `ff3bac63`, §3, transcribes them, but the transcripts are
```

```
$ grep '^WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"ImpliedExchange" -> ' | wc -l
3
```

```
$ grep '^WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"Exchange" -> ' | wc -l
0
```

```
$ grep '^WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"ImpliedExchange" -> ' | wc -l
3
```

```
$ grep '^WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_RADAR' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"Exchange" -> ' | wc -l
0
```

```
$ grep '^WL_LOCAL_S9B_MASS_BALANCE' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"ImpliedExchange" -> ' | wc -l
0
```

```
$ grep '^WL_LOCAL_S9B_MASS_BALANCE' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"Exchange" -> ' | wc -l
1
```

## Both legs: the line element and the optical operands
```
$ grep -E '^(WL|PY)_S9B_R_FLOW_ONLY_B_' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | grep -c -i 'metric\|lineelement'
0
```

```
$ grep -c '^PY_LOCAL_S9B_VECTOR_DISPERSION' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
1
```

```
$ grep -c '^WL_LOCAL_S9B_FERMAT_OBJECT' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
1
```

```
$ grep '^WL_LOCAL_S9B_SUPPLIED_INPUTS' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"Metric" -> ' | wc -l
1
```

```
$ git -C /var/projects/toy_physics show ff3bac63:research/pde_ledger_v3/steps/S9b_light_bending.md | sed -n '68,73p'
The supplied local leading-eikonal dispersion is

```
(ω − Vⁱkᵢ)² = cγ² gⁱʲkᵢkⱼ,
gᵢⱼ = δᵢⱼ + ∂ᵢξw ∂ⱼξw,       cγ² = μ⊥/ρbr,       ξw = ℓh.
```
```

```
$ grep -n 'insensitive to any overall multiplicative' -B2 -A1 /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_codex_evidence/visser.txt
1180-equation
1181-                                 hµν ∂µ ϕ ∂ν ϕ = 0.                             (86)
1182:This Eikonal equation is blatantly insensitive to any overall multiplicative pref-
1183-actor (conformal factor).
```

```
$ sed -n '49,67p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md
## What the check computes
For each oracle above and each engine, the check prints the following.
1. **The oracle's object.** Hamilton–Lisle's river velocity profile and its mass-function form. Visser's flowing-fluid
   background: its velocity, its continuity and density conditions, and its acoustic line element as the paper
   states them. Each is quoted verbatim, with section and equation number.
2. **The identification map.** The correspondence between S9b's symbols and the oracle's, for each condition §5
   names:
   - velocity, and the frame it is measured in;
   - free-fall versus lab time;
   - units, and `c₀` against the oracle's light or sound speed;
   - spatial dimension;
   - mass calibration: S9b's orbital `GM` against the oracle's mass parameter;
   - for Visser, also continuity, exchange and density, against S9b's implied `j_n` and live `ρ_br`.

   Each entry either cites the source that fixes it or prints `UNMATCHED`.
3. **The residual.** The difference between each engine's printed object and the oracle's object, under each
   identification the sources leave open. It is computed at S9b's retained order. Where the oracle's object is exact
   or covers more than S9b's restriction, the check prints the order used and what is left out. The comparison
   covers the general solution with its integration constant live, as well as the retained branch.
```

## Grok: Visser's two backgrounds, and the conditions section 5 names
```
$ grep -n -i 'continuity\|r^{-3/2}\|-3/2\|incompressible\|(53)\|(55)' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_grok_evidence/24_visser_s7_s8.txt | head -16
32:Given the Painlevé–Gullstrand line element, it might seem trivial to force the acoustic metric into this form: simply take ρ\\rho and cc to be constants, and set v\=2​G​M/rv=\\sqrt{2GM/r}? While this certainly forces the acoustic metric into the Painlevé–Gullstrand form the problem with this is that this assignment is incompatible with the continuity equation ∇→⋅(ρ​v→)≠0
34:The best we can actually do is this: Pick the speed of sound cc to be a position independent constant, which we normalize to unity (c\=1c=1). Now set v\=2​G​M/rv=\\sqrt{2GM/r}, and use the continuity equation ∇→⋅(ρ​v→)\=0\\vec{\\nabla}\\cdot(\\rho\\vec{v})=0 to deduce ρ​|v→|∝1/r2\\rho|\\vec{v}|\\propto 1/r^{2} so that ρ∝r−3/2\\rho\\propto r^{-3/2}. Since the speed of 
36:ds2∝r−3/2\[−dt2+(dr±2​G​Mrdt)2+r2(dθ2+sin2θdϕ2)\].ds^{2}\\propto r^{-3/2}\\left\[-dt^{2}+\\left(dr\\pm\\sqrt{2GM\\over r}dt\\right)^{2}+r^{2}\\left(d\\theta^{2}+\\sin^{2}\\theta\\;d\\phi^{2}\\right)\\right\].
38:(53)
46:We can turn this argument around and ask: Given a spherically symmetric flow of incompressible fluid, what is the acoustic metric? What is the corresponding line element in Schwarzschild coordinates? If we start by assuming incompressibility and spherical symmetry, then since ρ\\rho is position independent the continuity equation implies v∝1/r2v\\propto 1/r^{2}. But if ρ\\rho is position in
56:(55)
78:We should only use this canonical metric for the fluid region outside the bubble, and only in the approximation that the ambient fluid is incompressible (e.g. water). For the typically gaseous and necessarily compressible medium inside the bubble (e.g. air) we should use a separate acoustic metric. The two acoustic metrics need not be continuous across the bubble wall.
```

```
$ git -C /var/projects/toy_physics show 36f59ff1:research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md | grep 'M. Visser (1998)' | cut -c330-1000
qc/9712010). | Candidate for checking a separately derived flowing-fluid limit, its background continuity/force accounting and acoustic/GR comparison. The main comparison assumes barotropic, inviscid, irrotational fluid dynamics in three spatial dimensions; individual sink examples have additional density/EOS conditions. There is no shear-bearing sheet or off-sheet exchange in those examples. Thes
```

