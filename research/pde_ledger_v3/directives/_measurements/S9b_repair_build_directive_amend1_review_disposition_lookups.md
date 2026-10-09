# Measurements — S9b repair build directive amendment 1, review round 0 (generated 2026-10-08 18:44)

Generator: `_scratch/s9b_build/gen/s9b_amend1_r0_lookups.sh` (sha256sum/grep/sed/count only).

## Reviewed version
```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r0.md
18b00243467976c649c42be58c7f1e4dbb4e7649d1767d85b3a0cca0e7b01950  _scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r0.md
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend1_review_baseline_r0.sha256
18b00243467976c649c42be58c7f1e4dbb4e7649d1767d85b3a0cca0e7b01950  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## Verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend1_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend1_review_r0_grok.txt
_scratch/s9b_build/s9b_amend1_review_r0_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_amend1_review_r0_grok.txt:1:Verdict: CLEAR
```

## C1: the forward case names Part B's condition only; Part C's forward-case rows are absent
```
$ grep -n "Part B's condition\|Part C" _scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r0.md
10:**Amendment 1** (user, 2026-10-08). The forward case in item 10 now prints Part B's condition, where it printed
39:   the supplied Part C responses and the supplied references. The radial profiles are the ansatz.
59:   - Each Part B and Part C condition is solved relative to `GM` for every far-zone `b`. The solving is case by
70:9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
88:      - At both stages, also print Part B's condition for each observable under the same substitution, reduced as
97:      - Part B's condition, reduced as in item 6, with its implied `j_n`.
102:    - For Part C's fixed-ratio and power responses, also print each condition restricted to `V ≡ 0`, `ξ_w ≡ 0`
150:    - **K9a–b, Part C responses:** in the fixed-ratio response, the local `c_s(x)` is replaced by `c_s0`; in the
181:  - The outgoing bind-set is the reduced Part B and Part C conditions and the `j_n` each implies.
```

```
$ sed -n '102,104p' _scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r0.md
    - For Part C's fixed-ratio and power responses, also print each condition restricted to `V ≡ 0`, `ξ_w ≡ 0`
      (bulk density alone), with its domain and the `j_n` it implies through the supplied mass balance, with
      `ρ_br` symbolic.
```

## C1: the spec's Part C keeps V live in every row
```
$ sed -n '226,234p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
- **Part C.** The Part B conditions rewritten, to first order in `f`, under three supplied responses of `c_γ`
  to the local bulk density, with `V` and `ξ_w` left live in each:
  - `c_γ(x) ≡ c₀` (`δ(x) ≡ 0`);
  - pointwise fixed ratio: `c_γ(x)/c_s(x) = c₀/c_s0`;
  - symbolic power response: `c_γ(x)/c₀ = (ρ(x)/ρ₀)^s`, with `s` symbolic.

  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.
```

## C1: the leg's script prints forward Part C rows that carry Phi (count only; payloads withheld)
```
$ grep -c '^Part_C_forward_.*Phi' _scratch/s9b_build/s9b_amend1_review_r0_codex_evidence/derive.stdout
12
```

```
$ grep -o '^Part_C_forward_[a-z_]*' _scratch/s9b_build/s9b_amend1_review_r0_codex_evidence/derive.stdout
Part_C_forward_constant_speed_theta_residual
Part_C_forward_constant_speed_radar_residual
Part_C_forward_constant_speed_theta_condition_
Part_C_forward_constant_speed_radar_condition_
Part_C_forward_fixed_ratio_theta_residual
Part_C_forward_fixed_ratio_radar_residual
Part_C_forward_fixed_ratio_theta_condition_
Part_C_forward_fixed_ratio_radar_condition_
Part_C_forward_power_theta_residual
Part_C_forward_power_radar_residual
Part_C_forward_power_theta_condition_
Part_C_forward_power_radar_condition_
```

