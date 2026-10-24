# Measurements: S15 decision list, review pass 3 lookups

Generator: `_scratch/s15/gen/s15_dl_p3_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval by line, literal-match counts, hashes. The legs'
scripts were not re-run here; their filed stdout is retrieved.

## R0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/dl/dl_baseline_p3.sha256
research/pde_ledger_v3/directives/S15_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_setup_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/S15_decision_list_review_p2_disposition.md: OK
research/pde_ledger_v3/directives/_measurements/S15_decision_list_review_p2_lookups.md: OK
```

```
$ sha256sum _scratch/s15/dl/dl_review_prompt_p3.md
b144f964030af23a92ae82b13336f1acafa197bc54407980974feb982babd92e  _scratch/s15/dl/dl_review_prompt_p3.md
```

```
$ cat /var/projects/_review_scratch/s15-dl-p3-codex/DONE /var/projects/_review_scratch/s15-dl-p3-grok/DONE
exit 0
exit 0
```

```
$ git status --porcelain | /usr/bin/grep -v -F 'S15_decision_list_review_p3_lookups.md' | diff - <(sed -n '1p' _scratch/s15/dl/status_before_dl_p3.txt) && echo 'status unchanged apart from this file'
status unchanged apart from this file
```

```
$ sed -n '2p' _scratch/s15/dl/status_before_dl_p3.txt; git rev-parse HEAD
e393430ae8d17fec0b371e3388e5a59210414e9c
e393430ae8d17fec0b371e3388e5a59210414e9c
```

## R1. The verdicts

```
$ /usr/bin/grep -n -F 'CLEAR' /var/projects/_review_scratch/s15-dl-p3-codex/REVIEW.md /var/projects/_review_scratch/s15-dl-p3-grok/REVIEW.md | cut -c1-200
/var/projects/_review_scratch/s15-dl-p3-codex/REVIEW.md:10:**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
/var/projects/_review_scratch/s15-dl-p3-grok/REVIEW.md:52:**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
```

```
$ /usr/bin/grep -n -F 'FINDINGS' /var/projects/_review_scratch/s15-dl-p3-codex/REVIEW.md /var/projects/_review_scratch/s15-dl-p3-grok/REVIEW.md | cut -c1-200; echo "FINDINGS lines: $(cat /var/projects/_review_scratch/s15-dl-p3-codex/REVIEW.md /var/projects/_review_scratch/s15-dl-p3-grok/REVIEW.md | /usr/bin/grep -c -F 'FINDINGS')"
FINDINGS lines: 0
```

## R2. D6: the title is gone; the outcome words left outside the record section

```
$ sed -n '15p' research/pde_ledger_v3/directives/S15_decision_list.md
- **The plan's step.** S15 (`V3_STEP_PLAN.md:652–669`) (L1). The plan:
```

```
$ /usr/bin/grep -n -F 'attract' research/pde_ledger_v3/directives/S15_decision_list.md; echo "attract: $(/usr/bin/grep -c -F 'attract' research/pde_ledger_v3/directives/S15_decision_list.md)"
attract: 0
```

```
$ /usr/bin/grep -n -E 'expect|outcome' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-160
126:A name may name an object, never its value or a comparison outcome. The script prints objects and never states
127:conclusions. No acceptance criterion references an expected value.
```

```
$ sed -n '1,8p;28,34p;44,46p' /var/projects/_review_scratch/s15-dl-p3-grok/stdout/text_audit.stdout
record_section_lines = 143..173
n_lines = 196
WORD 'attract' total=0 outside=0
WORD 'repel' total=0 outside=0
WORD 'inverse' total=1 outside=0
  RECORD 147:     a named sum of printed pieces, a stated dependence on the separation, including which give an inverse-square one,
WORD 'sign' total=4 outside=3
  OUTSIDE 22:   - limits the sign claim to "the **leading matter-stress** contribution only" (`:666–668`).
WORD 'expect' total=1 outside=1
  OUTSIDE 127: conclusions. No acceptance criterion references an expected value.
WORD 'outcome' total=1 outside=1
  OUTSIDE 126: A name may name an object, never its value or a comparison outcome. The script prints objects and never states
WORD 'pull' total=1 outside=0
  RECORD 154:   on the named object, with all its conditions, those on `p_br`, `η` and `ζ` included. It is not called a pull or the
WORD 'leading' total=1 outside=1
  OUTSIDE 191: step gives no force law, no normalization, no mass, no map from `J` to a drain strength, and no force on a throat.
WORD 'F_12' total=0 outside=0
WORD 'square' total=1 outside=0
```

## R3. D7: the repaired record sentences; the live quantities in each item-1 piece (leg stdout)

```
$ sed -n '145,149p;152,156p' research/pde_ledger_v3/directives/S15_decision_list.md
- **What the record reads from the printed objects.** From the printed objects alone, it states:
  - which conditions on the profiles, the remainders, `p_br`, `η`, `ζ`, `Π_o` and `f` give a named printed piece, or
    a named sum of printed pieces, a stated dependence on the separation, including which give an inverse-square one,
    as the user's "Keep it general" asks. Each such dependence is stated with its conditions on every live quantity
    the piece contains, and on the surfaces: both radii, their rates, and the motion of `X_2`;
- **Requirements.** The record files each requirement as the register's rest-on test admits it, or on a user decision
  that it quotes. The user's "Keep it general" files gravity's flow requirement for S21's comparison, as a condition
  on the named object, with all its conditions, those on `p_br`, `η` and `ζ` included. It is not called a pull or the
  force on a throat (Retained OPEN). It merges into `R-S21-01` where the register's merge test says the two cover the
  same requirement.
```

```
$ cat /var/projects/_review_scratch/s15-dl-p3-grok/stdout/live_quantities.stdout
MEMBERSHIP
carried_current: u1,u2,rho1,rho2,varrho,W,X1,X2
pressure: rho1,rho2,varrho,p_br,X1,X2
pressure_expansion_first_order: rho1,rho2,varrho,p_br,X1,X2
viscous: u1,u2,W,eta,zeta,X1,X2
Pi_o: Pi_o
transport_motion_of_X2: u1,u2,rho1,rho2,varrho,W,X1,X2,dX2_dt
transport_radius_rate: u1,u2,rho1,rho2,varrho,W,X1,X2,R,dR_dt
viscous_flux_Derivative_eta_count 0
viscous_flux_Derivative_zeta_count 0
expansion_Derivative_p_br_count 1
div_viscous_Derivative_eta_count 3
div_viscous_Derivative_zeta_count 3
oint_p_nx_p_equals_rho_squared -8*pi*R**3*d/3
list_pressure_piece_x_p_equals_rho_squared 8*pi*R**3*d/3
oint_p_nx_p_equals_rho_series_through_R3 4*pi*R**3*(R**2 - 5*d**2)/(15*d**2)
pressure_laws_differ True
pA_value_at_b b**2
pB_value_at_b b**2
pA_prime_at_b 2*b
pB_prime_at_b 4*b
expansion_flux_ratio_pB_over_pA 2
div_T_visc_x (4*eta/3 + zeta)*exp(-sx)
div_T_visc_x_eta1_zeta0 4*exp(-sx)/3
div_T_visc_x_eta0_zeta1 exp(-sx)
div_T_visc_x_both0 0
ball_integral_exp_minus_x 4*pi*exp(-1)
ball_integral_positive True
surface_Txx_difference_at_x_eq_pm_R 0 0
div_difference_eta_var_minus_const 4*a*(-R**2 + sx**2 - 2*sx)*exp(-sx)/3
```

```
$ cat /var/projects/_review_scratch/s15-dl-p3-codex/_scratch/s11c/piece_dependencies/stdout
Generic pointwise sensitivity audit; no surface integrals or separation law computed.
carried_current : density_profiles_and_remainder, velocity_profiles_and_remainder, surface_normal_and_measure
whole_pressure : density_profiles_and_remainder, surface_normal_and_measure, p_br_function_including_any_derivatives
viscous_stress : velocity_spatial_derivatives, eta, zeta, surface_normal_and_measure
other_current : Pi_o, surface_normal_and_measure
centre_transport : density_profiles_and_remainder, velocity_profiles_and_remainder, surface_normal_and_measure, X2_motion
radius_transport : density_profiles_and_remainder, velocity_profiles_and_remainder, surface_normal_and_measure, radius_rate_and_direction
pressure_linear_choice : density_profiles_and_remainder, surface_normal_and_measure, p_br_function_including_any_derivatives
pressure_omission : density_profiles_and_remainder, surface_normal_and_measure, p_br_function_including_any_derivatives
d(whole_pressure_x)/d(rho1) = -nx*surface_measure*Subs(Derivative(p_br(_xi_1), _xi_1), _xi_1, remrho + rho1 + rho2)
d(viscous_stress_x)/d(eta) = 4*A0*nx*surface_measure/3 + A1*ny*surface_measure + A10*ny*surface_measure - 2*A11*nx*surface_measure/3 + A2*nz*surface_measure + A20*nz*surface_measure - 2*A22*nx*surface_measure/3
d(viscous_stress_x)/d(zeta) = A0*nx*surface_measure + A11*nx*surface_measure + A22*nx*surface_measure
pressure_linear_choice_contains_p_br_derivative = True
pressure_omission_contains_p_br_derivative = True
c_comp_squared_source_definition = dp_br/d(rho_br); no additional independent coefficient used
```

## R4. Nothing else: the supplied equations and the computed objects are unchanged from revision 2

```
$ diff <(git show cd853e16:research/pde_ledger_v3/directives/S15_decision_list.md | sed -n '/^## Supplied equations/,/^## Model point/p') <(git show e393430a:research/pde_ledger_v3/directives/S15_decision_list.md | sed -n '/^## Supplied equations/,/^## Model point/p') && echo 'Supplied equations through What the build computes: identical'
Supplied equations through What the build computes: identical
```

```
$ git diff --stat cd853e16 e393430a -- research/pde_ledger_v3/directives/S15_decision_list.md
 research/pde_ledger_v3/directives/S15_decision_list.md | 18 ++++++++++--------
 1 file changed, 10 insertions(+), 8 deletions(-)
```

## R5. The legs' scripts and stdout files exist (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-dl-p3-codex/scripts/source_audit.py /var/projects/_review_scratch/s15-dl-p3-codex/scripts/piece_dependencies.py /var/projects/_review_scratch/s15-dl-p3-codex/stdout/source_audit.stdout /var/projects/_review_scratch/s15-dl-p3-codex/_scratch/s11c/piece_dependencies/stdout /var/projects/_review_scratch/s15-dl-p3-grok/scripts/text_audit.py /var/projects/_review_scratch/s15-dl-p3-grok/scripts/live_quantities.py /var/projects/_review_scratch/s15-dl-p3-grok/stdout/text_audit.stdout /var/projects/_review_scratch/s15-dl-p3-grok/stdout/live_quantities.stdout | sed 's|/var/projects/_review_scratch/||'
c5c3cad36114bb60370f6e0c48b1c73347dbfec8145fcdeec15c1e44d428e31f  s15-dl-p3-codex/scripts/source_audit.py
0663a68785f0294431ecc0193c875f5228ea2a842fe7bd4d9c8ea4275c1d5733  s15-dl-p3-codex/scripts/piece_dependencies.py
cb9be42a5d9098ed3d6fdb3766893ccc7b10924dd556b61bbe6aa7f0dca0fb3f  s15-dl-p3-codex/stdout/source_audit.stdout
922bac345b2a751dd3ad916788519008249ecf234713afa7d8d84b120c9b249d  s15-dl-p3-codex/_scratch/s11c/piece_dependencies/stdout
b2a3de01e82711a1f679ab8fde42634cbf704569b2d75673976c39f694beb882  s15-dl-p3-grok/scripts/text_audit.py
aed1e3d44d641817a3eb4d7ccbd3e88271c5e727e5b1ea315fb9456317abb284  s15-dl-p3-grok/scripts/live_quantities.py
2070b295604a2f2028263bc6608ff80911c084224ad08b1ab7ca274642b423cd  s15-dl-p3-grok/stdout/text_audit.stdout
bae8499da22b1fe38d161c5ab1f5f2674ef4880d112c953f3c703730e69cdc14  s15-dl-p3-grok/stdout/live_quantities.stdout
```

```
$ wc -c /var/projects/_review_scratch/s15-dl-p3-grok/stdout/text_audit.err /var/projects/_review_scratch/s15-dl-p3-grok/stdout/live_quantities.err /var/projects/_review_scratch/s15-dl-p3-codex/stdout/piece_dependencies_guard.stderr | sed 's|/var/projects/_review_scratch/||'
0 s15-dl-p3-grok/stdout/text_audit.err
0 s15-dl-p3-grok/stdout/live_quantities.err
0 s15-dl-p3-codex/stdout/piece_dependencies_guard.stderr
0 total
```
