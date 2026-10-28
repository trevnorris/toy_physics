# Measurements: S15 reading decision list, review pass 1 lookups

Generator: `_scratch/s15/gen/s15_reading_p1_lookups.sh` (GNU grep, sed, cat, git, sha256sum, `LC_ALL=C`). Each block
is a command and its literal output: existence, verbatim retrieval, literal-match counts and hashes. The legs'
scripts were not re-run here; their filed stdout is retrieved. `X/` is `/var/projects/_review_scratch/s15-readdl-p1-codex/`; `G/` is `/var/projects/_review_scratch/s15-readdl-p1-grok/`.

## M0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/reading_dl/reading_dl_baseline_p1.sha256
research/pde_ledger_v3/directives/S15_reading_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_setup_lookups.md: OK
_scratch/s15/gen/s15_reading_setup_lookups.sh: OK
```

```
$ sha256sum _scratch/s15/reading_dl/reading_dl_review_prompt_p1.md
a005df3ec8b5f30284695593609106bf835ad99a554b5451345a517c04c8335d  _scratch/s15/reading_dl/reading_dl_review_prompt_p1.md
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p1-codex/DONE /var/projects/_review_scratch/s15-readdl-p1-grok/DONE; ls -la /var/projects/_review_scratch/s15-readdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p1-grok/REVIEW.md | cut -c24-
exit 0
exit 0
 trevnorris  5012 Oct 26 06:54 /var/projects/_review_scratch/s15-readdl-p1-codex/REVIEW.md
 trevnorris 15239 Oct 26 07:16 /var/projects/_review_scratch/s15-readdl-p1-grok/REVIEW.md
```

```
$ head -12 /var/projects/_review_scratch/s15-readdl-p1-codex/codex.log | grep -E '^(model|reasoning effort|sandbox):'
model: gpt-6.1-sol
sandbox: danger-full-access
reasoning effort: xhigh
```

```
$ diff <(git status --porcelain | grep -v -F -e 'S15_reading_decision_list_review_p1_lookups.md' -e 'S15_reading_decision_list_review_p1_disposition.md'; git rev-parse HEAD) _scratch/s15/reading_dl/status_before_reading_dl_p1.txt && echo 'status unchanged apart from the p1 disposition and lookups files'
status unchanged apart from the p1 disposition and lookups files
```

## M1. The verdicts

```
$ grep -n -F -e '**FINDINGS**' -e 'Verdict: **FINDINGS**' /var/projects/_review_scratch/s15-readdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p1-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-120
s15-readdl-p1-codex/REVIEW.md:12:**FINDINGS**
s15-readdl-p1-grok/REVIEW.md:5:Verdict: **FINDINGS**. Three defects change what item 2 computes or what the record may c
```

```
$ grep -n -E '^[0-9]\. \*\*HIGH|^### [0-9]\. HIGH' /var/projects/_review_scratch/s15-readdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p1-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-160
s15-readdl-p1-codex/REVIEW.md:14:1. **HIGH — The promised inverse-square conditions solve only the geometry part.**
s15-readdl-p1-codex/REVIEW.md:36:2. **HIGH — The inherited pressure conditions do not define every new separation derivative.**
s15-readdl-p1-grok/REVIEW.md:28:### 1. HIGH — the separation derivative of the whole pressure, and of the whole minus the expansion, is not defined on the con
s15-readdl-p1-grok/REVIEW.md:103:### 2. HIGH — `E_Q^resp` of the first-order pressure needs a second derivative at `ρ_2 + ϱ`, which the carried condition do
s15-readdl-p1-grok/REVIEW.md:125:### 3. HIGH — `E_Q^geo` and `E_Q^resp` are not unique under "every term that contains a change `F′`"
s15-readdl-p1-grok/REVIEW.md:175:1. **HIGH.** The separation derivative of the whole pressure, and of the whole minus the expansion, does not exist under "defin
s15-readdl-p1-grok/REVIEW.md:176:2. **HIGH.** `E_Q^resp` of the first-order expansion does not exist under "differentiable at `ρ_2 + ϱ`" when `F′` moves tha
s15-readdl-p1-grok/REVIEW.md:177:3. **HIGH.** "Every term that contains a change `F′`" does not define a unique split of an unexpanded integral. The two readi
```

## M2. F1 (Codex 1): the family conditions solve only the geometric part

```
$ sed -n '40,43p;62,68p;92,93p;105,107p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
- **The departure from the inverse square:** `E_Q = 𝒟Q + 2Q`. `Q` falls off as the inverse square of the separation
  along the family exactly where `E_Q = 0` for every `D` on the domain.
- **Its two parts:** `E_Q = E_Q^geo + E_Q^resp`. `E_Q^resp` holds every term that contains a change `F′`;
  `E_Q^geo` holds the rest.
3. `E_Q^geo` on the profile family, with its dependence on `d` explicit and exact. Where the integrand holds no general
   field besides the family's amplitudes, the region's quantities and the velocities, the integrals are evaluated in
   closed form. Otherwise the object is printed with the family substituted and left unevaluated.
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
  - for each `Q`, the printed conditions under which its departure vanishes, and `E_Q^resp`;
  - gravity's flow requirement for S21, stated as those printed conditions;
For each object S15 prints: its departure from the inverse square of the separation, exact, for general profiles,
with every quantity's change along the family live. On the power-law profile family, the conditions under which each
closed-form departure vanishes. The claim is conditional on P5, P6, `R-S14a-01` and the labelled choices. It gives
```

```
$ sed -n '148,149p' research/pde_ledger_v3/directives/S15_decision_list.md
    as the user's "Keep it general" asks. Each such dependence is stated with its conditions on every live quantity
    the piece contains, and on the surfaces: both radii, their rates, and the motion of `X_2`;
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/outcome.json | tr -d ' \n' | cut -c1-120; echo
{"exitCode":0,"wallSeconds":1.426123300974723,"unit":"s11c-guard-3d8009b78371","stderrBytes":0,"limitsVerified":true,"ch

```

```
$ sed -n '1,14p' /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/stdout
F1 domain: d > R > r > 0; X1=d*e_z; X2=0; both cores excluded
F1 Q: shell momentum, labels_drain1+drain2, z component; b1=b2=0
rho1*V2 angular integral: 0
a1=-2 angular integral: -4*pi/d**2
a1=-2 Q: 4*pi*A1*B2*(-R**3 + r**3)/(3*d**2)
A1_prime=A1, B2_prime=B2: (Egeo/Q, Eresp/Q, E/Q) = (0, 2, 2)
same family A1=alpha*d, B2=beta*d; Q = 4*pi*alpha*beta*(-R**3 + r**3)/3
a1=1, A1_prime=A1, B2_prime=-4*B2: (Egeo/Q, Eresp/Q, E/Q) = (3, -3, 0)
same family A1=alpha*d, B2=beta/d^4; d^2*Q = 4*pi*alpha*beta*(-R**3 + r**3)/3
transverse azimuthal integrals = (0, 0)
a1=-2 all-direction Egeo = Matrix([[0, 0, 0]])
a1=-2 all-direction Eresp-2Q = Matrix([[0, 0, 0]])
a1=1 all-direction Egeo-3Q = Matrix([[0, 0, 0]])
a1=1 all-direction Egeo+Eresp = Matrix([[0, 0, 0]])
```

## M3. F2 (Grok 3): the split by term

```
$ sed -n '42,43p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
- **Its two parts:** `E_Q = E_Q^geo + E_Q^resp`. `E_Q^resp` holds every term that contains a change `F′`;
  `E_Q^geo` holds the rest.
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p1-grok/03_split_not_unique.stdout
Q: W/d**2
Dscr_expression: Integral(-2*W/d**2 + Wprime/d**2, (s, 0, 1))
Dscr_top_level_term_count: 1
Dscr_contains_Wprime: True
Dscr_value: (-2*W + Wprime)/d**2
E_Q: Wprime/d**2
reading_top/E_geo: 2*W/d**2
reading_top/E_resp: (-2*W + Wprime)/d**2
reading_top/sum: Wprime/d**2
reading_taylor/Dscr_at_Fprime_0: -2*W/d**2
reading_taylor/E_geo: 0
reading_taylor/E_resp: Wprime/d**2
reading_taylor/sum: Wprime/d**2
E_geo_readings_differ: True
E_resp_readings_differ: True
both_sums_equal_E: True
sphere_profile/Dscr: Integral(A*R**2*(2*R*d*cos(theta) - 2*d**2)*W(theta)*sin(theta)/(R**2 - 2*R*d*cos(theta) + d**2)**2 + A*R**2*Wprime(theta)*sin(theta)/(R**2 - 2*R*d*cos(theta) + d**2), (theta, 0, pi))
sphere_profile/top_level_terms: 1
sphere_profile/top_term_contains_Wprime: True
sphere_profile/integrand_geo_piece_count: 2
sphere_profile/integrand_resp_piece_count: 1
sphere_profile/geo_piece_has_Wprime: False
sphere_profile/resp_piece_has_Wprime: True
```

```
$ wc -c < /var/projects/_review_scratch/s15-readdl-p1-grok/03_split_not_unique.stderr
0
```

```
$ sed -n '24,26p' /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/stdout
CHECK integrating-factor residual = 0
CHECK instantaneous jets (delta_X1,delta_v1,delta_X2,delta_v2) = (d, j, 0, 0)
CHECK decomposition change (delta_Egeo,delta_Eresp,delta_E) = (-A*K*d, A*K*d, 0)
```

## M4. F3 (Grok 1): the whole pressure and the whole minus the expansion

```
$ sed -n '23,25p;73,74p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
- Everything S15's list supplies, unchanged: the sections "Supplied equations" and "The drains, their fields and the
  region", and the domain and pressure conditions (`:120–121`). This is revision 3, unchanged since `91248b67` (L2).
- The objects of that list's "What the build computes", items 1 and 2, as it names them. The accepted S15 build
Every object is exact on S15's domain, and each pressure object carries its condition. No limit or expansion is taken
in place of an exact object. Where an engine does not obtain an object exactly, it prints a label saying so, beside
```

```
$ sed -n '110,114p;120,121p' research/pde_ledger_v3/directives/S15_decision_list.md
   **The pressure.** The pressure of the total density is printed whole, where `p_br` is defined at `ρ_br`. Beside
   it, as a labelled choice, its expansion through first order in `ρ_1` about `ρ_2 + ϱ`, grouped as above (a
   function of `ρ_2 + ϱ` carries the labels drain 2 and remainder), with its order printed, and the whole pressure
   minus the expansion. The expansion exists only where `p_br` is defined and differentiable at `ρ_2 + ϱ`; that
   condition is labelled where it is printed.
Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, and each pressure
object where its condition above holds.
```

```
$ grep -n -F -e 'p_br` is a general function' -e 'c_comp(ρ_br)' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | cut -c1-200
367:c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
370:`p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
```

```
$ grep -E '^(whole|whole_smooth|subtract)/' /var/projects/_review_scratch/s15-readdl-p1-grok/02_pressure_domain.stdout; wc -c < /var/projects/_review_scratch/s15-readdl-p1-grok/02_pressure_domain.stderr
whole/p: Abs(rho - 2)
whole/defined_at_total: 0
whole/defined_at_base: 1
whole/pprime_at_base_plus: -1
whole/pprime_at_base_minus: -1
whole/base_derivative_exists: True
whole/pprime_at_total_plus: 1
whole/pprime_at_total_minus: -1
whole/total_derivative_exists: False
whole/integrand_quotient_plus: area
whole/integrand_quotient_minus: -area
whole/separation_derivative_exists: False
whole_smooth/derivative: 4*area
whole_smooth/chain_rule: 4*area
whole_smooth/match: True
subtract/whole_value: 0
subtract/expansion_value: 0
subtract/value_at_0: 0
subtract/quotient_plus: 2
subtract/quotient_minus: 0
subtract/separation_derivative_exists: False
0
```

```
$ sed -n '23p' /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/stdout
whole pressure family derivative = -10*pi/3
```

## M5. F4 (Codex 2, Grok 2): the response part of the first-order expansion

```
$ sed -n '15,22p' /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/stdout
F2 domain: d > 1 > r_in > 0; rho1=rho2=1; varrho(eps,x)=eps*z_cart
p(rho) = z + Abs(z - 1)**(3/2)
(p(1), p_prime(1), p(2)) = (1, 1, 3)
p_prime left/right limits at 1 = (1, 1)
outer first_order_choice Q(0) = 0
outer first_order_choice Q(+h) = -4*pi*(9*sqrt(h) + 5*h)/15
outer first_order_choice Q(-h) = 4*pi*(9*sqrt(h) + 5*h)/15
family derivative right/left limits = (-oo, -oo)
```

```
$ grep -E '^expansion/' /var/projects/_review_scratch/s15-readdl-p1-grok/02_pressure_domain.stdout
expansion/p: Piecewise((-rho**2/2 + rho - 1/2, rho < 1), (rho**2/2 - rho + 1/2, True))
expansion/p_simplified: Piecewise((-rho**2/2 + rho - 1/2, rho < 1), (rho**2/2 - rho + 1/2, True))
expansion/p_defined_at_base: 0
expansion/pprime_at_base: 0
expansion/coeff_quotient_plus: 1
expansion/coeff_quotient_minus: -1
expansion/coeff_derivative_exists: False
expansion/piece_at_0: 0
expansion/piece_quotient_plus: rho_1
expansion/piece_quotient_minus: -rho_1
expansion/response_derivative_exists: False
expansion/geo_derivative: 0
expansion/geo_defined: True
```

## M6. Checks the Grok leg found sound (filed stdout)

```
$ grep -E '^(radial_vector|magnitude_only|one_point|no_separation|no_X1)/' /var/projects/_review_scratch/s15-readdl-p1-grok/01_inverse_square_operator.stdout | cut -c1-160
radial_vector/Q: Matrix([[C*D_x/(D_x**2 + D_y**2 + D_z**2)**(3/2)], [C*D_y/(D_x**2 + D_y**2 + D_z**2)**(3/2)], [C*D_z/(D_x**2 + D_y**2 + D_z**2)**(3/2)]])
radial_vector/E: Matrix([[0], [0], [0]])
radial_vector/E_is_zero: True
magnitude_only/abs: d**(-2)
magnitude_only/E: Matrix([[-sin(d)/d], [cos(d)/d], [0]])
magnitude_only/E_is_zero: False
one_point/E: 2 - 2*d0/d
one_point/E_at_d0: 0
one_point/E_at_2d0: 1
no_separation/Dscr: 0
no_separation/E: 2*Q_no_separation
```

```
$ grep -E '^(special|limit_matches|r1_min)' /var/projects/_review_scratch/s15-readdl-p1-grok/04_profile_closed_form.stdout | cut -c1-160
limit_matches_special: True
special/Dscr: pi*R*(4*R*d + (R - d)*(R + d)*log((R - d)**2/(R + d)**2))/(d*(R - d)*(R + d))
special/E: pi*R*(4*R*d - (R - d)*(R + d)*log((R - d)**2/(R + d)**2))/(d*(R - d)*(R + d))
special/E_identically_zero: False
special/series_through_x3: 4*pi*R**2/d**2
special/E_of_that_truncation: 0
r1_min: -R + d
r1_min_minus_d_plus_R: 0
```

```
$ grep -E '^(held_velocity|held_rate|jet|jet_vs_path|limits_scaled|split)' /var/projects/_review_scratch/s15-readdl-p1-grok/05_family_coherence.stdout | cut -c1-160
held_velocity/Q: C*V2dot/d**2
held_velocity/E: 0
held_velocity/V2_still_present: True
held_velocity/E_zero_while_V2_live: True
held_rate_constant_in_d/Q: C*Rdot
held_rate_constant_in_d/E: 2*C*Rdot
held_rate_constant_in_d/E_is_zero: False
jet/rate: C*V1dot/d**2 + V2dot
jet/Dscr: C*(-2*V1dot + V1dot_prime)/d**2
jet/V2_coefficient_in_Dscr: 0
jet/Fp_coefficient: C
jet/V1_and_V2_independent: True
jet_vs_path/jet_Dscr: -2*C*V1dot/d**2
jet_vs_path/path_Dscr: C*(-V1dot - V2dot)/d**2
jet_vs_path/differ: True
limits_scaled/Dscr: (R**3 + 2*r_in**3)/(3*d**2)
limits_scaled_minus_held: R**3/d**2
split/value_at_0_A: -A0/d**2
split/value_at_0_B: -A0/d**2
split/values_match: True
split/Dscr_A: 2*A0/d**2
split/Dscr_B: 0
split/derivatives_differ: True
```

```
$ for f in /var/projects/_review_scratch/s15-readdl-p1-grok/0[1-6]_*.stderr; do echo "$(basename $f) $(wc -c < $f)"; done
01_inverse_square_operator.stderr 0
02_pressure_domain.stderr 0
03_split_not_unique.stderr 0
04_profile_closed_form.stderr 0
05_family_coherence.stderr 0
06_tag_inventory.stderr 0
```

## M7. The stray files a leg run left in the repository root

```
$ ls -la _scratch/s15/reading_dl/stray_p1/ | tail -n +4 | cut -c24-; cat _scratch/s15/reading_dl/stray_p1/02_pressure_domain.stderr
 trevnorris  112 Oct 26 07:12 02_pressure_domain.stderr
 trevnorris    0 Oct 26 07:12 02_pressure_domain.stdout
python3: can't open file '/var/projects/toy_physics/02_pressure_domain.py': [Errno 2] No such file or directory
```

## M8. The legs' reports, scripts and stdout (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-readdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p1-codex/physics_checks.py /var/projects/_review_scratch/s15-readdl-p1-codex/vector_and_split_checks.py /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/stdout /var/projects/_review_scratch/s15-readdl-p1-codex/_scratch/s11c/vector_checks/stdout /var/projects/_review_scratch/s15-readdl-p1-grok/REVIEW.md /var/projects/_review_scratch/s15-readdl-p1-grok/0[1-6]_*.py /var/projects/_review_scratch/s15-readdl-p1-grok/0[1-6]_*.stdout | sed 's|/var/projects/_review_scratch/||'
0ee4237aaa5f3f1502e570ff9d67cbcff742fccc569061ed383321b5df9e229a  s15-readdl-p1-codex/REVIEW.md
c78f7d341fe63e9f4607a20725e31af269053b5b90cc5c3d77f7bd29c13cb0d3  s15-readdl-p1-codex/physics_checks.py
47d39ee9580a9134f2aed2f38586e47337fd87801258f1624ac9aecdb87c2a6b  s15-readdl-p1-codex/vector_and_split_checks.py
b3bfb08c1b249ce53b2839d4466c64a582ec78b06a862a3ff9b66970ed299dce  s15-readdl-p1-codex/_scratch/s11c/physics_checks_final/stdout
d6f37c8fb92f3182f7568acf9cc4ad94856e0c6d976707bc53ccca2d4270961c  s15-readdl-p1-codex/_scratch/s11c/vector_checks/stdout
99d7adec5d132abfb261ad7f596643f10a58f05a3b0ce7d9b8ebcc26d8a251c0  s15-readdl-p1-grok/REVIEW.md
350709a67ffe59265685142868cee948f0f2b6ab94bd40a45f4edf978778f0d2  s15-readdl-p1-grok/01_inverse_square_operator.py
20018ce7e1283288ea4d027f1dc4e63341c42bb0a157015a841985529584c0f0  s15-readdl-p1-grok/02_pressure_domain.py
6918559ff1c4523a34fbffaf30abd42f6a996033ad0647eb855fba630b878219  s15-readdl-p1-grok/03_split_not_unique.py
f3335e944af568e4b1d0c88dd7a1bb3a4cef9284d01539fbf909743d35bc7b34  s15-readdl-p1-grok/04_profile_closed_form.py
839db7c51b0bf8203370ae1cacf34d110997c5ca102fd57070e5a375dbbcb859  s15-readdl-p1-grok/05_family_coherence.py
55cb657691bae15022f83d51497e0d56250ca3724f87cecc28db94e74fb1bee3  s15-readdl-p1-grok/06_tag_inventory.py
02c7f737c02337c76c3cf4d3b1865a983b2508812978f5441e0c0db2d6b7301d  s15-readdl-p1-grok/01_inverse_square_operator.stdout
14d5b626e018cf8b8a7cd6f2529b5a491f015c55044266a9f5a0e3035c8229f5  s15-readdl-p1-grok/02_pressure_domain.stdout
e3553cc817f330d2f12e8f045e25450200370f95cb0d26928cd2b3dbeadcc3d9  s15-readdl-p1-grok/03_split_not_unique.stdout
ae7585469076c51cbbad0ee2523080565aefb8d73eb7d81f64a84486f30cc7b4  s15-readdl-p1-grok/04_profile_closed_form.stdout
5c10750c3e8203121eeb54b3140673ff57ca20c50498642e3141f47e015c9151  s15-readdl-p1-grok/05_family_coherence.stdout
888fb82231d73059f604bb0bb0c54c573c81ffc4e5c800150c838a1a36bfb16b  s15-readdl-p1-grok/06_tag_inventory.stdout
```
