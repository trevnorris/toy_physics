# Measurements: S15 decision list, review pass 1 lookups

Generator: `_scratch/s15/gen/s15_dl_p1_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval by line, literal-match counts, hashes. The legs'
scripts were not re-run here; their filed stdout is retrieved.

## P0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum research/pde_ledger_v3/directives/S15_decision_list.md research/pde_ledger_v3/directives/_measurements/S15_setup_lookups.md
6cc19b6658ea09ea0b125ac1ebfde7275bb648880668fbe58792b6c07eff5bc2  research/pde_ledger_v3/directives/S15_decision_list.md
0076cc5d423a49ec3576de7cc242bb5663228d046dcf561f1ff01f3189d969e1  research/pde_ledger_v3/directives/_measurements/S15_setup_lookups.md
```

```
$ cat _scratch/s15/dl/dl_baseline_p1.sha256
6cc19b6658ea09ea0b125ac1ebfde7275bb648880668fbe58792b6c07eff5bc2  research/pde_ledger_v3/directives/S15_decision_list.md
0076cc5d423a49ec3576de7cc242bb5663228d046dcf561f1ff01f3189d969e1  research/pde_ledger_v3/directives/_measurements/S15_setup_lookups.md
```

```
$ sha256sum _scratch/s15/dl/dl_review_prompt_p1.md
d348351dc3e5e2f0fa3b3c434ca30426c80d128d50ae82366e58356a065b0a46  _scratch/s15/dl/dl_review_prompt_p1.md
```

```
$ cat /var/projects/_review_scratch/s15-dl-p1-codex/DONE /var/projects/_review_scratch/s15-dl-p1-grok/DONE
exit 0
exit 0
```

```
$ git log -1 --format='%h %s'
4f907e9d Accept the S16 step record and register pass (round 2 clear)
```

```
$ git status --porcelain | /usr/bin/grep -v -F 'S15_decision_list_review_p1_lookups.md' | diff - <(sed -n '1,3p' _scratch/s15/dl/status_before_dl_p1.txt) && echo 'status unchanged apart from this file'
status unchanged apart from this file
```

```
$ sed -n '4p' _scratch/s15/dl/status_before_dl_p1.txt; git rev-parse HEAD
4f907e9d7788a125e11e955ffdb4f6007e85b917
4f907e9d7788a125e11e955ffdb4f6007e85b917
```

```
$ head -1 /var/projects/_review_scratch/s15-dl-p1-codex/REVIEW.md | cut -c1-200; /usr/bin/grep -n -F '6cc19b6658ea09ea' /var/projects/_review_scratch/s15-dl-p1-grok/REVIEW.md | cut -c1-160
Three findings affect what the builder sees, whether an object exists, and what the record may claim. Reviewed baseline: `6cc19b6658ea09ea0b125ac1ebfde7275bb648880668fbe58792b6c07eff5bc2`; lookups als
4:**sha256:** `6cc19b6658ea09ea0b125ac1ebfde7275bb648880668fbe58792b6c07eff5bc2` (matches `_scratch/s15/dl/dl_baseline_p1.sha256`).
```

## P1. Codex 1: the option texts in the list; the G13 row; M2

```
$ sed -n '31,42p' research/pde_ledger_v3/directives/S15_decision_list.md
## User decisions (2026-10-24)
The user selected these options. The wording is the orchestrator's option text. The build directive carries the
choices, not the options' reasons (M2).
- **Stress law, "Use P5 + P6":** "S15 adopts P5 and P6 as labelled premises. It computes the momentum the two drains'
  combined flow carries into a region around one drain, piece by piece (carried flow momentum, pressure, viscosity).
  The result is conditional on P5/P6 and on S14a's open far-flow requirement. It is not yet called the force on the
  throat."
- **Inflow, "Keep it general":** "Any radial inflow and density profile around each drain. S15 then says which
  profiles give an inverse-square pull, and that becomes gravity's flow requirement, compared with light's in S21.
  S14a left the far flow open, so supplying it would freeze it."
- **Motion, "Let them move":** "Keep each drain's velocity live. A moving drain that swallows material picks up a
  drag, the classic objection to 'gravity as inflow' (prior-art row G13). Holding the drains still would hide it."
```

```
$ sed -n '423p' docs/medium_requirements_and_prior_art.md | cut -c1-300
| `G13` objections to the sink ether | **KNOWN** | A moving drain feels a drag: Wang's formula (32) gives *"a universal damping force"*. Energy disappears without explanation, some drag must arise, and no process that creates matter is known (Wikipedia, secondary only). For pulsation, every pulsatio
```

```
$ sed -n '81,82p;85,86p' CLAUDE.md
the object; let the engine hand over what it built. The spec says what to **compute** — ⛔ never what
anything equals, is expected, or was measured. Withhold exactly one thing: an acceptance criterion
Naming an object does not waive first-principles derivation, form ablation, or evidence. Prior art's result
for our object is never a supplied premise (M3).
```

```
$ for p in drag G13 'inverse-square' 'freeze it' 'hide it'; do echo "$p: $(/usr/bin/grep -c -F -- "$p" research/pde_ledger_v3/directives/S15_decision_list.md)"; done
drag: 1
G13: 2
inverse-square: 2
freeze it: 1
hide it: 1
```

## P2. Codex 2: the remainders, the expansion, the domain; the EOS in the S9b spec; the leg stdout

```
$ sed -n '85,87p;102,105p;111p' research/pde_ledger_v3/directives/S15_decision_list.md
- **The total fields** are `V = V_1 + V_2 + W` and `ρ_br = ρ_1 + ρ_2 + ϱ`. The remainders `W` and `ϱ` are general
  live fields, and they hold whatever the two profiles do not, including any background flow. This is a labelled
  split of the total fields, not a superposition law. Nothing makes the remainders small or zero.
   Within each piece, the terms that contain fields of both drains are printed apart from the terms that contain
   one drain's fields only, and every term that contains a remainder is printed apart. The pressure of the total
   density does not split this way. It is printed whole, and beside it, as a labelled choice, its expansion through
   first order in `ρ_1` about `ρ_2 + ϱ`, with its order printed.
Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded.
```

```
$ sed -n '367,370p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
```

`p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
```

```
$ /usr/bin/grep -n -F 'EOS' research/pde_ledger_v3/directives/S15_decision_list.md; echo "count: $(/usr/bin/grep -c -F 'EOS' research/pde_ledger_v3/directives/S15_decision_list.md)"
count: 0
```

```
$ cat /var/projects/_review_scratch/s15-dl-p1-codex/stdout/pressure_domain.stdout
geometry: separation=5, R=1, r_in=1/2, probe=(0,1,0)
rho1_at_probe = 27
rho2_at_probe = 1
remainder_at_probe = -1
base_density = 0
total_density_at_probe = 27
p(total)_at_probe = log(27)
dp_drho(total)_at_probe = 1/27
scaled_expansion_path_at_probe = log(27*lam)
limit_pressure_lambda_to_0_plus = -oo
limit_first_derivative_lambda_to_0_plus = oo
smooth_probe_exact_pressure_minus_first_order = a**2
```

```
$ cat /var/projects/_review_scratch/s15-dl-p1-codex/stdout/pressure_guard.stdout
{
  "exitCode": 0,
  "wallSeconds": 0.48389899800531566,
  "unit": "s11c-guard-d8169f4593c6",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 0.416975181986345,
    "childPid": 230477
  }
}
```

## P3. Codex 3: the transport pieces and the record reading; the S16 record's momentum object; the leg stdout

```
$ sed -n '95,100p;134,141p' research/pde_ledger_v3/directives/S15_decision_list.md
1. **The momentum crossing each surface.** For the inner and the outer surface, two parts, kept apart, in the
   coordinates `x` of the supplied equations:
   - **The rest part:** `−Π^{ij}` contracted with the surface normal and integrated over the surface, in four pieces
     kept apart: the carried current, the pressure, the viscous stress and `Π_o`.
   - **The transport of `G` by the surface's motion,** in two pieces: from the motion of `X_2`, and from the change
     of the surface's radius.
- **What the record reads from the printed objects.** From the printed objects alone, it states:
  - which conditions on the profiles, the remainders, `Π_o` and `f` give the crossing momentum a stated dependence on
    the separation, including which give an inverse-square one, as the user's "Keep it general" asks;
  - the terms that depend on each drain's velocity;
  - the sign of a printed piece only on its stated conditions, and only for that piece (`:666–668`).
- **Requirements.** The record files each requirement as the register's rest-on test admits it, or on a user decision
  that it quotes. The user's "Keep it general" files gravity's flow requirement for S21's comparison. It merges into
  `R-S21-01` where the register's merge test says the two cover the same requirement.
```

```
$ sed -n '152,160p' research/pde_ledger_v3/steps/S16_moving_region.md
`momentumReceived` is a three-vector made from the integrated `f_i` and the two **rest**
momentum-flux surface pieces. `momentumRate` prints the full regional rate, and
`momentumSurfaceMotion` prints the motion pieces separately. The shared twelve surface-motion
components are bound one by one. Py's latter payload additionally contains two `Volume(0)`
bulk-reference entries per component, which WL's payload omits; those values/keys remain
unpaired, as do their six extra Py dimension tags. The corresponding separately emitted
rate pieces are compared elsewhere. No equality of received momentum and `dP/dt` is imposed,
and no moving-frame reinterpretation replaces the supplied `x`-coordinate object. Momentum
received is not a passive factor or the force a drain exerts on another body. (Retrieve
```

```
$ /usr/bin/grep -n -F 'crossing momentum' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "count: $(/usr/bin/grep -c -F 'crossing momentum' research/pde_ledger_v3/directives/S15_decision_list.md)"
135:  - which conditions on the profiles, the remainders, `Π_o` and `f` give the crossing momentum a stated dependence on
162:- **Whether the crossing momentum acts on the throat.** That needs the momentum exchanged at the faces (`R-S16-03`, no
count: 2
```

```
$ cat /var/projects/_review_scratch/s15-dl-p1-codex/stdout/control_surface_witness.stdout
total_velocity = Matrix([[beta*(t*v2 - x) + 1, -beta*y, -beta*z]])
independent_center_velocities = (v1, v2)
G = Matrix([[-beta*(-t*v2 + x) + 1, -beta*y, -beta*z]])
Pi = Matrix([[3*beta + (-beta*(-t*v2 + x) + 1)**2 + 1, -beta*y*(-beta*(-t*v2 + x) + 1), -beta*z*(-beta*(-t*v2 + x) + 1)], [-beta*y*(-beta*(-t*v2 + x) + 1), beta**2*y**2 + 3*beta + 1, beta**2*y*z], [-beta*z*(-beta*(-t*v2 + x) + 1), beta**2*y*z, beta**2*z**2 + 3*beta + 1]])
local_balance_lhs_minus_f = Matrix([[0, 0, 0]])
outward_material_flux_about_drain2_radius_b = -4*pi*b**3*beta
rest_crossing_outer = Matrix([[16*pi*R**3*beta/3, 0, 0]])
rest_crossing_inner = Matrix([[-16*pi*a**3*beta/3, 0, 0]])
integrated_f = Matrix([[R**3*(4*pi*beta*v2/3 - 16*pi*beta/3) - a**3*(4*pi*beta*v2/3 - 16*pi*beta/3), 0, 0]])
rest_plus_integrated_f = Matrix([[4*pi*beta*v2*(R**3 - a**3)/3, 0, 0]])
center_motion_both_surfaces = Matrix([[4*pi*beta*v2*(-R**3 + a**3)/3, 0, 0]])
outer_radius_transport = Matrix([[4*pi*R**2*Rd, 0, 0]])
inner_radius_transport = Matrix([[-4*pi*a**2*ad, 0, 0]])
region_momentum = Matrix([[4*pi*R**3/3 - 4*pi*a**3/3, 0, 0]])
stored_rate_minus_all_balance_terms = Matrix([[0, 0, 0]])
common_instant_geometry: d>1, R=1, r_in=1/2
same_fields_case_A_Rdot=1/d^2_crossing = Matrix([[4*pi/d**2, 0, 0]])
same_fields_case_B_Rdot=1_crossing = Matrix([[4*pi, 0, 0]])
case_A_d2_times_crossing_x = 4*pi
case_B_d2_times_crossing_x = 4*pi*d**2
rest_plus_integrated_f_same_in_both_cases = Matrix([[4*pi*beta*v2*(R**3 - a**3)/3, 0, 0]])
test_member_v2=4_integrated_f = Matrix([[0, 0, 0]])
test_member_v2=4_full_surface_crossing_A = Matrix([[4*pi/d**2, 0, 0]])
test_member_v2=4_full_surface_crossing_B = Matrix([[4*pi, 0, 0]])
```

```
$ cat /var/projects/_review_scratch/s15-dl-p1-codex/stdout/control_complete_guard.stdout
{
  "exitCode": 0,
  "wallSeconds": 1.5851471979985945,
  "unit": "s11c-guard-b771785412bc",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 1.518383560993243,
    "childPid": 230845
  }
}
```

## P4. Grok 1: the normal in the list; P6's traction; the leg stdout, sections 3 and 4

```
$ /usr/bin/grep -n -F 'normal' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160; echo "normal: $(/usr/bin/grep -c -F 'normal' research/pde_ledger_v3/directives/S15_decision_list.md)"; echo "outward: $(/usr/bin/grep -c -F 'outward' research/pde_ledger_v3/directives/S15_decision_list.md)"; echo "into: $(/usr/bin/grep -c -F 'into' research/pde_ledger_v3/directives/S15_decision_list.md)"
17:    with `{J, Θ_Q, N_∞,3, Q}` and the normalization in S22's debt table (`:1306–1307`);
23:  - No force law, coefficient or normalization is printed, and no printed object is called the force on a throat.
97:   - **The rest part:** `−Π^{ij}` contracted with the surface normal and integrated over the surface, in four pieces
113:The build computes no force law, coefficient, normalization or mass. It gives no relation of any printed object to
164:- **S22:** the source map `{J, Θ_Q, N_∞,3, Q}`, the normalization and `I_F,12` (`V3_STEP_PLAN.md:1306–1307`), and
176:step gives no force law, no normalization, no mass, no map from `J` to a drain strength, and no force on a throat.
normal: 6
outward: 0
into: 3
```

```
$ /usr/bin/grep -n -F 't_br,inplane' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | cut -c1-200
366:t_br,inplane^i ≡ T^{ij} n_j ,      F_br,inplane^i ≡ ∂_j T^{ij} ,
```

```
$ /usr/bin/grep -n -E 'volume_integral_of_div|surface_integral_region_outward|surface_integral_centre_radial_both_faces|outward_minus_volume|radial_both_minus_volume|convention_difference|_piece_factor' /var/projects/_review_scratch/s15-dl-p1-grok/s15_dl_checks.out
38:volume_integral_of_div = 4*pi*(R**2*h(R) - r_in**2*h(r_in))
39:surface_integral_region_outward = 4*pi*R**2*h(R) - 4*pi*r_in**2*h(r_in)
40:surface_integral_centre_radial_both_faces = 4*pi*R**2*h(R) + 4*pi*r_in**2*h(r_in)
41:outward_minus_volume = 0
42:radial_both_minus_volume = 8*pi*r_in**2*h(r_in)
43:convention_difference_radial_minus_outward = 8*pi*r_in**2*h(r_in)
51:outer_centre_piece_factor = Xdot1*er1 + Xdot2*er2 + Xdot3*er3
52:outer_radius_piece_factor = Rdot*(er1**2 + er2**2 + er3**2)
53:inner_centre_piece_factor = -Xdot1*er1 - Xdot2*er2 - Xdot3*er3
54:inner_radius_piece_factor = rdot*(-er1**2 - er2**2 - er3**2)
```

## P5. Grok 2: the grouping sentence; the leg stdout, section 7

```
$ sed -n '102,104p' research/pde_ledger_v3/directives/S15_decision_list.md
   Within each piece, the terms that contain fields of both drains are printed apart from the terms that contain
   one drain's fields only, and every term that contains a remainder is printed apart. The pressure of the total
   density does not split this way. It is printed whole, and beside it, as a labelled choice, its expansion through
```

```
$ /usr/bin/grep -n -E 'terms_in_both_drains_and_remainder|partition_|_each_monomial_in_one_bin|naive_overlapping' /var/projects/_review_scratch/s15-dl-p1-grok/s15_dl_checks.out
96:terms_in_both_drains_and_remainder = [2*V1*V2*rem_rho, 2*V1*W*rho2, 2*V2*W*rho1]
97:naive_overlapping_sum_minus_product = 2*V1*V2*rem_rho + 2*V1*W*rho2 + 2*V2*W*rho1
98:partition_P_sum_minus_product = 0
99:partition_P_both = V1**2*rho2 + 2*V1*V2*rho1 + 2*V1*V2*rho2 + V2**2*rho1
100:partition_Q_sum_minus_product = 0
101:partition_Q_both = V1**2*rho2 + 2*V1*V2*rem_rho + 2*V1*V2*rho1 + 2*V1*V2*rho2 + 2*V1*W*rho2 + V2**2*rho1 + 2*V2*W*rho1
102:partition_Q_both_minus_P_both = 2*V1*V2*rem_rho + 2*V1*W*rho2 + 2*V2*W*rho1
103:P_each_monomial_in_one_bin = True
104:Q_each_monomial_in_one_bin = True
```

## P6. The legs' scripts and stdout files exist (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-dl-p1-codex/scripts/source_audit.py /var/projects/_review_scratch/s15-dl-p1-codex/scripts/pressure_domain.py /var/projects/_review_scratch/s15-dl-p1-codex/scripts/control_surface_witness.py /var/projects/_review_scratch/s15-dl-p1-codex/stdout/source_audit.stdout /var/projects/_review_scratch/s15-dl-p1-codex/stdout/pressure_domain.stdout /var/projects/_review_scratch/s15-dl-p1-codex/stdout/control_surface_witness.stdout /var/projects/_review_scratch/s15-dl-p1-grok/s15_dl_checks.py /var/projects/_review_scratch/s15-dl-p1-grok/s15_dl_checks.out | sed 's|/var/projects/_review_scratch/||'
615ce4e514b35dfa50abf1e8b1663e5ec8a6da0bd0e3d3d8df5fa85cd5b9fd4d  s15-dl-p1-codex/scripts/source_audit.py
959cd8da32dbac23e6d97b99159748fc0f551865f0ad146f465b02466acafc85  s15-dl-p1-codex/scripts/pressure_domain.py
3169f9db4d64217bcee1487e1bc01ceac5d1623265ffdff9345d98604cdf3541  s15-dl-p1-codex/scripts/control_surface_witness.py
de5492e483f191385068832c9d773f131f3c0caf86a2ff14cbe540fcca3c4dee  s15-dl-p1-codex/stdout/source_audit.stdout
12e5799a98ff5cff2ea4e679ab99db0de58f0ef4ffb7e7b46c9ece5d48600509  s15-dl-p1-codex/stdout/pressure_domain.stdout
166b908973757b2791befaa9a3418977378aecbfeb87a18e5ef771712bed6ed6  s15-dl-p1-codex/stdout/control_surface_witness.stdout
9c834c8f62862cd765366d38462a337217d20f90c65c66511783d4a0906df79a  s15-dl-p1-grok/s15_dl_checks.py
9fd730e780f95975e45da04fc882243dee964dc53c307505a5d1021f87004536  s15-dl-p1-grok/s15_dl_checks.out
```
