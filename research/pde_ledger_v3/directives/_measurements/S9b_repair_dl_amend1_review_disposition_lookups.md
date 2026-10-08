# Measurements — S9b repair DL amendment-1 review dispositions (generated 2026-10-08 12:39)

Generator: `_scratch/s9b_build/gen/s9b_repair_dl_amend1_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed version is the frozen `_scratch/s9b_build/S9b_repair_decision_list_amend1_reviewed_v0.md`.

```
$ cut -c1-64 _scratch/s9b_build/s9b_repair_dl_amend1_review_baseline.sha256; sha256sum _scratch/s9b_build/S9b_repair_decision_list_amend1_reviewed_v0.md | cut -c1-64
7586d8cfdaec4ef730a7b9b9ab81e6959b10c94f3bdd01d224af6593cbc4b682
7586d8cfdaec4ef730a7b9b9ab81e6959b10c94f3bdd01d224af6593cbc4b682
```

```
$ grep -n -o 'Verdict:[* ]*[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_dl_amend1_review_codex_final.txt _scratch/s9b_build/s9b_repair_dl_amend1_review_grok.txt _scratch/s9b_build/s9b_repair_dl_amend1_review_codex_blocked1_final.txt
_scratch/s9b_build/s9b_repair_dl_amend1_review_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_dl_amend1_review_grok.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_dl_amend1_review_codex_blocked1_final.txt:1:Verdict: BLOCKED 
```

## Leg launch — the first Codex leg could not run the guard
```
$ grep -n 'Failed to connect to bus' _scratch/s9b_build/s9b_repair_dl_amend1_review_codex_blocked1_final.txt
13:   Failed to connect to bus: Operation not permitted (consider using --machine=<user>@.host --user to connect to bus of other user)
```

```
$ grep -n '^sandbox:' _scratch/s9b_build/s9b_repair_dl_amend1_review_codex_blocked1.txt _scratch/s9b_build/s9b_repair_dl_amend1_review_codex.txt
_scratch/s9b_build/s9b_repair_dl_amend1_review_codex_blocked1.txt:8:sandbox: workspace-write [workdir, /tmp, $TMPDIR]
_scratch/s9b_build/s9b_repair_dl_amend1_review_codex.txt:8:sandbox: danger-full-access
```

## A1 (Codex 1 = Grok 1) — P5's selection includes stating what the OPEN pieces must supply
```
$ grep -n 'Part D then states what those open' _scratch/s9b_build/s9b_repair_dl_amend1_review_prompt.md
16:  stress). Any other in-plane momentum current stays OPEN as now, and Part D then states what those open pieces must
```

```
$ grep -c 'what those open pieces must supply\|what the OPEN pieces' _scratch/s9b_build/S9b_repair_decision_list_amend1_reviewed_v0.md
0
```

```
$ grep -n '^- which of\|^- for each condition' _scratch/s9b_build/S9b_repair_decision_list_amend1_reviewed_v0.md
43:- which of `δ`, `V`, `ρ_br` and `j_n` it determines relative to `GM`;
44:- which of them stay free;
45:- for each condition, the `j_n` it implies.
```

```
$ grep -n 'arbitrary_OPEN_substitution_residual\|OPEN_current_shift_residual' _scratch/s9b_build/s9b_repair_dl_amend1_review_codex_evidence/physics.stdout
19:P5_P6_arbitrary_OPEN_substitution_residual = 0
20:P5_P6_OPEN_current_shift_residual = 0
```

```
$ grep -n 'sum_minus_Jstar\|identity_residual_r' _scratch/s9b_build/s9b_repair_dl_amend1_review_grok_evidence/s9b_dl_amend1_physics.stdout
27:identity_residual_r 0
30:sum_minus_Jstar 0
```

## A2 (Grok 2) — the pressure term had no tensor, so its force sign was free
```
$ grep -o 'isotropic pressure `p_br(ρ_br)` plus a linear viscous stress: [^|]*' _scratch/s9b_build/S9b_repair_decision_list_amend1_reviewed_v0.md | cut -c1-200
isotropic pressure `p_br(ρ_br)` plus a linear viscous stress: `η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k`. `p_br` is a general function of `ρ_br` only. Its compressio
```

```
$ grep -n 'omega2_over_k2\|^s -1\|^s 1' _scratch/s9b_build/s9b_repair_dl_amend1_review_grok_evidence/s9b_dl_amend1_physics.stdout
20:omega2_over_k2 [-c_comp2*s]
21:s -1 [c_comp2]
22:s 1 [-c_comp2]
```

```
$ sed -n '112,113p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
state. No sharp-sheet or finite-slab material reduction, stress measure, relaxation family, passive
sign or physical support is selected. Sources: spec §§1–6; contract §§1, 3–8 (M1).
```

```
$ sed -n '333,340p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**P6 — steady in-plane pressure (2026-10-08; D2).** In Part D, the supplied steady in-plane stress
is isotropic pressure, as the P1 material relaxes under steady load. In the Cartesian Cauchy-stress
convention where stress contracted with a unit normal gives traction, the premise is

```
𝒯_br,inplane^live,ij ≡ −p_br(ρ_br) δ^{ij} ,
c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
```
```

## A3 (Grok 3) — the dissipated power was assigned to an O2 operand defined with no formula
```
$ grep -o 'Its dissipated power belongs to[^.]*\.' _scratch/s9b_build/S9b_repair_decision_list_amend1_reviewed_v0.md
Its dissipated power belongs to `𝒫_ref/relax^live`, whose supplier and budget stay OPEN.
```

```
$ grep -n 'the heating stays in the OPEN energy budget' _scratch/s9b_build/s9b_repair_dl_amend1_review_prompt.md
20:  is lost; the heating stays in the OPEN energy budget. The linear form excludes power-law creep such as glacier ice's
```

```
$ grep -n '^| `𝒫_ref/relax^live`' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
279:| `𝒫_ref/relax^live` | Explicit power associated with `ℛ_ref/strain^live`, with no formula, sign or vanishing assumed. | O2 requirement; relaxation/reference ownership unassigned. |
```

```
$ grep -n 'phi_minus_2eta_edev2_minus_zeta_theta2' _scratch/s9b_build/s9b_repair_dl_amend1_review_grok_evidence/s9b_dl_amend1_physics.stdout
32:phi_minus_2eta_edev2_minus_zeta_theta2 0
```

## What both legs found sound (D4)
```
$ grep -n 'reading_A_residual\|counterexample_j_coord\|difference_at_s_0' _scratch/s9b_build/s9b_repair_dl_amend1_review_grok_evidence/s9b_dl_amend1_physics.stdout
38:reading_A_residual 0
39:difference_at_s_0 0
42:counterexample_j_coord 0
```

