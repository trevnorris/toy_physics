# Measurements — S9b repair build directive item 10, review round 2 dispositions (generated 2026-10-08 15:28)

Generator: `_scratch/s9b_build/gen/s9b_item10_r2_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed version is the directive at `cb3809f6`.

```
$ git show cb3809f6:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sha256sum
2637cd52ef36afb367487cf181a885b19cab28ed2277102302dadb3f3d8b22c7  -
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_item10_review_baseline_r1.sha256
2637cd52ef36afb367487cf181a885b19cab28ed2277102302dadb3f3d8b22c7  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_item10_review_r2_codex_final.txt _scratch/s9b_build/s9b_item10_review_r2_grok.txt
_scratch/s9b_build/s9b_item10_review_r2_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_item10_review_r2_grok.txt:1:Verdict: NEEDS REVISION
```

## G1 — the forward-case round-trip time is listed 'as a function of b'
```
$ git show cb3809f6:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '80,82p'
      - Print the deflection, the round-trip excess time and its `ln(1/b²)` coefficient, both effective `γ`s, their
        difference and both Part B residuals against the references, each as a function of `b`. Also print both
        one-way excess times and their nonreciprocal part, with `Z_E` and `Z_R` kept.
```

```
$ sed -n '169,172p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
2. **Round-trip (radar) excess time.** An emitter at distance `Z_E` on one side of the mass and a reflector
   at `Z_R` on the other, along the line. Subtract the flat round-trip time. Here
   `r_E = √(b² + Z_E²)` and `r_R = √(b² + Z_R²)`. Part A prints the full excess time; the Part B comparison
   uses only the coefficient of `ln(1/b²)` in the regime `Z_E, Z_R ≫ b`, for every profile, including tails
```

```
$ grep -n 'leibniz_d_delta_excess_d_ZE\|leibniz_d_excess_d_ZE\|CHECK V2_excess_depends_on_ZE' _scratch/s9b_build/s9b_item10_review_r2_grok_evidence/derive_item10.stdout
26:leibniz_d_excess_d_ZE = Phi**2*Z_E**2/(8*pi**2*c_0**3*rho_0**2*(Z_E**2 + b**2)**3)
47:leibniz_d_delta_excess_d_ZE = -2*delta(sqrt(Z_E**2 + b**2))/c_0
63:CHECK V2_excess_depends_on_ZE = True
```

## G2 — the bulk-density-alone rows omit the implied j_n
```
$ git show cb3809f6:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '93,94p'
    - For Part C's fixed-ratio and power responses, also print each condition restricted to `V ≡ 0`, `ξ_w ≡ 0`
      (bulk density alone), with its domain.
```

```
$ sed -n '232,234p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.
```

```
$ grep -n 'j_n_when_V_r_0' _scratch/s9b_build/s9b_item10_review_r2_grok_evidence/derive_item10.stdout
52:j_n_when_V_r_0 = 0
```

