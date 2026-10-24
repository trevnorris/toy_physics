# Measurements: S15 decision list, review pass 2 lookups

Generator: `_scratch/s15/gen/s15_dl_p2_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval by line, literal-match counts, hashes. The legs'
scripts were not re-run here; their filed stdout is retrieved.

## Q0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/dl/dl_baseline_p2.sha256
research/pde_ledger_v3/directives/S15_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_setup_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/S15_decision_list_review_p1_disposition.md: OK
research/pde_ledger_v3/directives/_measurements/S15_decision_list_review_p1_lookups.md: OK
```

```
$ sha256sum _scratch/s15/dl/dl_review_prompt_p2.md
8b3955403b6a6037e26be828e6fb6fefb85bc992550fd9e36a5d6c7a8d07cee8  _scratch/s15/dl/dl_review_prompt_p2.md
```

```
$ cat /var/projects/_review_scratch/s15-dl-p2-codex/DONE /var/projects/_review_scratch/s15-dl-p2-grok/DONE
exit 0
exit 0
```

```
$ git status --porcelain | /usr/bin/grep -v -F 'S15_decision_list_review_p2_lookups.md' | diff - <(sed -n '1p' _scratch/s15/dl/status_before_dl_p2.txt) && echo 'status unchanged apart from this file'
status unchanged apart from this file
```

```
$ sed -n '2p' _scratch/s15/dl/status_before_dl_p2.txt; git rev-parse HEAD
cd853e16246072d9fddd3d83dc65fcc95b8d4a9b
cd853e16246072d9fddd3d83dc65fcc95b8d4a9b
```

```
$ /usr/bin/grep -n -F 'd83e69e5046ec71b' /var/projects/_review_scratch/s15-dl-p2-codex/REVIEW.md /var/projects/_review_scratch/s15-dl-p2-grok/REVIEW.md | cut -c1-200
/var/projects/_review_scratch/s15-dl-p2-codex/REVIEW.md:1:Reviewed revision 2 at `cd853e16`, SHA256 `d83e69e5046ec71be88acd240b4167971b466a2d07ed5f52936f4bf6d4239a01`. Sources were read before the lis
/var/projects/_review_scratch/s15-dl-p2-grok/REVIEW.md:3:Artifact: `research/pde_ledger_v3/directives/S15_decision_list.md` at `cd853e16246072d9fddd3d83dc65fcc95b8d4a9b`. sha256 `d83e69e5046ec71be88ac
```

## Q1. Codex 1: the plan's title on the builder's side; the record section's start; M2; the leg stdout

```
$ sed -n '14p' research/pde_ledger_v3/directives/S15_decision_list.md
- **The plan's step.** S15, "Two drains attract" (`V3_STEP_PLAN.md:652–669`) (L1). The plan:
```

```
$ /usr/bin/grep -n -E '^## ' research/pde_ledger_v3/directives/S15_decision_list.md
13:## Filing and framing
34:## User decisions (2026-10-24)
47:## Supplied equations (unfalsifiable in this build)
74:## The drains, their fields and the region
92:## What the build computes
128:## Model point
142:## For the record and its register pass (not in the build directive)
172:## Retained OPEN
186:## Intended claim
191:## Next
```

```
$ /usr/bin/grep -n -F 'attract' research/pde_ledger_v3/directives/S15_decision_list.md; echo "count: $(/usr/bin/grep -c -F 'attract' research/pde_ledger_v3/directives/S15_decision_list.md)"
14:- **The plan's step.** S15, "Two drains attract" (`V3_STEP_PLAN.md:652–669`) (L1). The plan:
count: 1
```

```
$ sed -n '652p' research/pde_ledger_v3/V3_STEP_PLAN.md | cut -c1-200
### S15 · Two drains attract ⚠ RUN S16 FIRST for regime context — ⛔ NOT as S15's premise
```

```
$ sed -n '81,82p' CLAUDE.md
the object; let the engine hand over what it built. The spec says what to **compute** — ⛔ never what
anything equals, is expected, or was measured. Withhold exactly one thing: an acceptance criterion
```

```
$ cat /var/projects/_review_scratch/s15-dl-p2-codex/stdout/answer_boundary.stdout
builder_side_attraction_hit_lines = [14]
LIST 14: - **The plan's step.** S15, "Two drains attract" (`V3_STEP_PLAN.md:652–669`) (L1). The plan:
record_only_section_line = 142
M2: The spec says what to **compute** — ⛔ never what anything equals, is expected, or was measured.
D1: The list states each choice as what is computed or kept live, with no option's reason and no prior-art result.
G9: two sources of like sign attract, and two of unlike sign repel, with forces varying inversely as the square of the distance
```

## Q2. Grok 1: the record's conditions; the live constitutive quantities in the list and in the S9b spec; the leg stdout

```
$ sed -n '144,149p' research/pde_ledger_v3/directives/S15_decision_list.md
- **What the record reads from the printed objects.** From the printed objects alone, it states:
  - which conditions on the profiles, the remainders, `Π_o` and `f` give a named printed piece, or a named sum of
    printed pieces, a stated dependence on the separation, including which give an inverse-square one, as the user's
    "Keep it general" asks. Each such dependence is stated with the conditions on the surfaces under which it holds:
    both radii, their rates, and the motion of `X_2`;
  - the terms that depend on each drain's velocity;
```

```
$ sed -n '57,59p;64p' research/pde_ledger_v3/directives/S15_decision_list.md
  G^i    = ρ_br V^i
  Π^{ij} = ρ_br V^i V^j − T^{ij} + Π_o^{ij}
  T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k
  - `p_br` is a general function of `ρ_br`. `η` and `ζ` are general live fields (P6).
```

```
$ sed -n '370,371p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
`p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
The shear viscosity `η` and bulk viscosity `ζ` are live general profiles, including their gradients.
```

```
$ for p in 'p_br' 'η' 'ζ'; do echo "$p in the record-section bullet (144-149): $(sed -n '144,149p' research/pde_ledger_v3/directives/S15_decision_list.md | /usr/bin/grep -c -F -- "$p")"; done
p_br in the record-section bullet (144-149): 0
η in the record-section bullet (144-149): 0
ζ in the record-section bullet (144-149): 0
```

```
$ sed -n '138,170p' /var/projects/_review_scratch/s15-dl-p2-grok/check_s15_dl_p2.out
========================================================================
13. Same profiles, two equations of state, two separation dependences
========================================================================
pressure_flux_x_p_equals_rho_through_s3 = 4*pi*d**3*s**5/15 - 4*pi*d**3*s**3/3
pressure_flux_x_p_equals_rho_squared_through_s3 = -8*pi*d**4*s**3/3
fluxes_differ = True
p_equals_rho_flux_at_fixed_R_truncated = -4*pi*R**3*(-R**2 + 5*d**2)/(15*d**2)
p_equals_rho_squared_flux_at_fixed_R_truncated = -8*pi*R**3*d/3
exact_flux_p_equals_rho_squared = -8*pi*R**3*d/3

========================================================================
14. Viscous piece depends on eta for an allowed radial profile
========================================================================
div_viscous_stress_x_at_unit_point_eta_1 = -4*exp(-1)
div_viscous_stress_x_at_unit_point_eta_0 = 0
eta_changes_the_divergence True

========================================================================
15. Same smooth remainder, shell flux of viscous stress
========================================================================
shell_div_x_eta_1_zeta_0 = 4*exp(-sx)/3
shell_div_x_eta_0_zeta_1 = exp(-sx)
shell_div_x_both_off = 0
prefactor_eta = 4/3
prefactor_zeta = 1
shell_integral_lower_bound_eta = -16*pi*(-R + a)*(R**2 + R*a + a**2)*exp(-R)/9
shell_integral_lower_bound_zeta = -4*pi*(-R + a)*(R**2 + R*a + a**2)*exp(-R)/3
eta_bound_positive_when_R_gt_a = True
zeta_bound_positive_when_R_gt_a = True
eta_bound_at_1_2 = 112*pi*exp(-2)/9
zeta_bound_at_1_2 = 28*pi*exp(-2)/3

========================================================================
```

```
$ cat /var/projects/_review_scratch/s15-dl-p2-grok/check_s15_dl_p2.err | wc -c
0
```

## Q3. The legs' scripts and stdout files exist (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-dl-p2-codex/scripts/answer_boundary.py /var/projects/_review_scratch/s15-dl-p2-codex/scripts/physics_witnesses.py /var/projects/_review_scratch/s15-dl-p2-codex/scripts/source_audit.py /var/projects/_review_scratch/s15-dl-p2-codex/stdout/answer_boundary.stdout /var/projects/_review_scratch/s15-dl-p2-codex/_scratch/s11c/physics_witnesses/stdout /var/projects/_review_scratch/s15-dl-p2-grok/check_s15_dl_p2.py /var/projects/_review_scratch/s15-dl-p2-grok/check_s15_dl_p2.out | sed 's|/var/projects/_review_scratch/||'
5b2721e69f9bb2de55a9ed1ed6fbc401c46df320cfcb9ebc04177ce13c87cb60  s15-dl-p2-codex/scripts/answer_boundary.py
2ad8e4dbb47a37b370f6f0aa36dfc605a322d0dca4614fe280dc39ea4fc42ba8  s15-dl-p2-codex/scripts/physics_witnesses.py
3e69e0ba2f693b1ab471b2c9cb08d2e414b80aa3a0db20f8cdab6d78a3768db1  s15-dl-p2-codex/scripts/source_audit.py
1132c9e936a2d8c680ef275bde5835b62994607d27b68156999ef931dd05cb4d  s15-dl-p2-codex/stdout/answer_boundary.stdout
4d591f6ea8f13da2cf57a492ebe29db9c1198acfa047cbbcbad6fc7063786a70  s15-dl-p2-codex/_scratch/s11c/physics_witnesses/stdout
b3867ffac7f9da619191eb7d43feefb9a8590a0710b8668ad8d993e85eace538  s15-dl-p2-grok/check_s15_dl_p2.py
2d129a2adb556ee23163fd87bc20e67c97d31aa8bc01298b297668caac0f7c64  s15-dl-p2-grok/check_s15_dl_p2.out
```
