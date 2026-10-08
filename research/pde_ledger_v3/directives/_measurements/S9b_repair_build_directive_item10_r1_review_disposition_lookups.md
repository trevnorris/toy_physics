# Measurements — S9b repair build directive item 10, review round 1 dispositions (generated 2026-10-08 15:07)

Generator: `_scratch/s9b_build/gen/s9b_item10_r1_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed version is the directive at `c41d217f`.

```
$ git show c41d217f:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sha256sum
7f2091cb2210c8a9f88e2dff3c675a4a4e8a160159716749c293c6c8267b02fa  -
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_item10_review_baseline.sha256
7f2091cb2210c8a9f88e2dff3c675a4a4e8a160159716749c293c6c8267b02fa  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep -n -o 'Verdict: NEEDS REVISION' _scratch/s9b_build/s9b_item10_review_codex_final.txt _scratch/s9b_build/s9b_item10_review_grok.txt
_scratch/s9b_build/s9b_item10_review_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_item10_review_grok.txt:1:Verdict: NEEDS REVISION
```

## The forward case as reviewed
```
$ git show c41d217f:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '66,80p'
   not built from another row, or from the amplitude before the response is applied.
10. **Forward case and single-mechanism restrictions (user requests, 2026-10-08).**
    - **Forward case: no far-zone loss.** This is the user's premise (2026-10-06/07/08): brane material leaves the
      brane only in throats, so `j_n ≡ 0` in the far zone, and the same brane mass crosses every sphere around the
      body.
      - Solve the supplied mass balance for `V` under this premise, with `ρ_br` live. Its integration constant is
        a new live symbol: the signed far-zone brane mass flux, with no value and no sign. Its relation to `GM` is
        open.
      - Substitute that `V` into Part A's observables, first with `δ ≡ 0`, `ξ_w ≡ 0`, then with `δ` and `ξ_w`
        live.
      - Print the deflection, the round-trip excess time and its `ln(1/b²)` coefficient, both effective `γ`s and
        both Part B residuals against the references, each as a function of `b`. Also print both one-way excess
        times and their nonreciprocal part, with `Z_E` and `Z_R` kept.
      - Print no matching condition for this case. Label every object with the premise.
    - For the deflection and for the radar `ln(1/b²)` coefficient, print each restriction below of:
```

```
$ git show c41d217f:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep -n 'far-zone flux'
41:   - Symbols that originate in this step: `n`, `K`, `m`, `GM`, `s`, `f`, `j_n`, `ρ₀`, and item 10's far-zone flux
```

## C1 — the forward case omits the difference between the two γs
```
$ grep -n 'The difference between the two' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
225:  - The difference between the two `γ`s, printed as an object.
```

## G1 — the flux constant has two readings
```
$ grep -n -A1 'FLUX_OVER_ODE_CONSTANT\|V_ODE_OVER_V_FLUX\|NONRECIP_ODE_OVER_FLUX\|ROUND_DENSITY_ODE_OVER_FLUX\|DPHI_DR_PLUS_4PI_R2_JN' _scratch/s9b_build/s9b_item10_review_grok_evidence/derive_item10.stdout
16:FLUX_OVER_ODE_CONSTANT
17-4*pi
--
19:DPHI_DR_PLUS_4PI_R2_JN
20-0
--
28:V_ODE_OVER_V_FLUX
29-4*pi
--
100:NONRECIP_ODE_OVER_FLUX
101-4*pi
--
103:ROUND_DENSITY_ODE_OVER_FLUX
104-16*pi**2
```

## G2 — the forward case drops the branch conditions
```
$ sed -n '177,181p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
them.
```

```
$ sed -n '211,212p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
- **Part A.** The branch-existence and path-traversal conditions. Then, for the LAB_HELD profile, as expressions
  or functionals of `δ`, `V` and `ξ_w`, with each retained monomial printed separately:
```

```
$ grep -n -A1 'LEG_CONDITIONS_IDENTICAL' _scratch/s9b_build/s9b_item10_review_grok_evidence/derive_item10.stdout
130:LEG_CONDITIONS_IDENTICAL
131-False
```

