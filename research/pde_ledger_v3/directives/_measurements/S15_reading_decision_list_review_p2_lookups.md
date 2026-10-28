# Measurements: S15 reading decision list, review pass 2 lookups

Generator: `_scratch/s15/gen/s15_reading_p2_lookups.sh` (GNU grep, sed, cat, find, git, sha256sum, `LC_ALL=C`). Each
block is a command and its literal output: existence, verbatim retrieval, literal-match counts and hashes. The legs'
scripts were not re-run here; their filed stdout is retrieved. `X/` is `/var/projects/_review_scratch/s15-readdl-p2-codex/`; `G/` is `/var/projects/_review_scratch/s15-readdl-p2-grok/`.

## M0. Baseline, prompt, legs finished, tree unchanged, no new files in the repository root

```
$ sha256sum -c _scratch/s15/reading_dl/reading_dl_baseline_p2.sha256
research/pde_ledger_v3/directives/S15_reading_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_setup_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_decision_list_review_p1_disposition.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_decision_list_review_p1_lookups.md: OK
```

```
$ sha256sum _scratch/s15/reading_dl/reading_dl_review_prompt_p2.md
db29c2f7b5369167b5b89070294a1cc0c22a9107ec51685c0a67acfc5f4710b2  _scratch/s15/reading_dl/reading_dl_review_prompt_p2.md
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p2-codex/DONE /var/projects/_review_scratch/s15-readdl-p2-grok/DONE; ls -la /var/projects/_review_scratch/s15-readdl-p2-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md | cut -c24-
exit 0
exit 0
 trevnorris  3697 Oct 28 11:57 /var/projects/_review_scratch/s15-readdl-p2-codex/REVIEW.md
 trevnorris 18521 Oct 28 12:13 /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md
```

```
$ head -12 /var/projects/_review_scratch/s15-readdl-p2-codex/codex.log | grep -E '^(model|reasoning effort|sandbox):'
model: gpt-6.1-sol
sandbox: danger-full-access
reasoning effort: xhigh
```

```
$ diff <(git status --porcelain | grep -v -F -e 'S15_reading_decision_list_review_p2_lookups.md' -e 'S15_reading_decision_list_review_p2_disposition.md'; git rev-parse HEAD) _scratch/s15/reading_dl/status_before_reading_dl_p2.txt && echo 'status unchanged apart from the p2 disposition and lookups files'
status unchanged apart from the p2 disposition and lookups files
```

```
$ find . -maxdepth 1 -newer _scratch/s15/reading_dl/reading_dl_review_prompt_p2.md ! -name .git | wc -l
0
```

## M1. The verdicts

```
$ grep -n -F -e '**FINDINGS**' -e '**CLEAR**' /var/projects/_review_scratch/s15-readdl-p2-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-120
s15-readdl-p2-codex/REVIEW.md:5:**FINDINGS**
s15-readdl-p2-grok/REVIEW.md:125:**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
```

```
$ grep -n -E '^[0-9]\. \*\*(HIGH|MEDIUM|LOW)' /var/projects/_review_scratch/s15-readdl-p2-codex/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-140
7:1. **MEDIUM — The new EOS conditions exclude separation derivatives that exist.**
```

## M2. N1 (Codex 1): the named pressure conditions are stated as necessary

```
$ sed -n '50,53p;71,72p;84,86p;119,122p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
- **Where they exist.** Each departure object is stated where the separation derivative of its object exists, and it
  carries the condition under which it does, labelled where it is printed. For instance, the separation derivative
  of the whole pressure, and of the whole minus the expansion, needs `p_br` differentiable at the total density. The
  response part of the expansion needs `p_br′` differentiable at `ρ_2 + ϱ`.
2. `E_Q^geo` and `E_Q^resp`, for general profiles, exact wherever they exist on the domain, with every live quantity
   kept.
Every object is exact on S15's domain wherever it exists, and it carries the conditions under which it exists: S15's
pressure conditions, and for each departure object the conditions above. No limit or expansion is taken
in place of an exact object. Where an engine does not obtain an object exactly, it prints a label saying so, beside
## Intended claim
For each object S15 prints: its departure from the inverse square of the separation, exact, for general profiles,
with every quantity's change along the family live, wherever it exists. On the power-law profile family, the
conditions under which each closed-form geometric departure vanishes: the conditions under which the object falls off
```

```
$ sed -n '110,114p' research/pde_ledger_v3/directives/S15_decision_list.md
   **The pressure.** The pressure of the total density is printed whole, where `p_br` is defined at `ρ_br`. Beside
   it, as a labelled choice, its expansion through first order in `ρ_1` about `ρ_2 + ϱ`, grouped as above (a
   function of `ρ_2 + ϱ` carries the labels drain 2 and remainder), with its order printed, and the whole pressure
   minus the expansion. The expansion exists only where `p_br` is defined and differentiable at `ρ_2 + ϱ`; that
   condition is labelled where it is printed.
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p2-codex/_scratch/s11c/domain_checks_final/stdout
domain: d>R>rin>0; B>0; rho0>0; 0<z0<R; cores disjoint from shell
expansion_EOS: rho + (rho - 1)*Abs(rho - 1)/2
expansion_pprime_at_base_left_right: (1, 1)
expansion_pprime_derivative_at_base_left_right: (-1, 1)
expansion_Qeps: 8*pi*R**3*d*(B + Bprime*eps)*(eps + 1)/3
expansion_Q: 8*pi*B*R**3*d/3
expansion_Egeo: 8*pi*B*R**3*d
expansion_Eresp: 8*pi*Bprime*R**3*d/3
expansion_Etotal: 8*pi*R**3*d*(3*B + Bprime)/3
expansion_response_regular_at_base_without_psecond: True
whole_density_family_residual: 0
whole_EOS_at_total: (d*eps - z + z0)**2
whole_EOS_derivative_at_rho0_left_right: (-1, 1)
whole_Qeps: 8*pi*R**3*(d*eps + z0)/3
whole_Q0: 8*pi*R**3*z0/3
whole_separation_derivative: 8*pi*R**3*d/3
whole_Etotal: 8*pi*R**3*(d + 2*z0)/3
whole_composite_pressure_spatial_gradient: -2*d*eps + 2*z - 2*z0
END: exact integrated counterexamples to necessary EOS conditions
```

```
$ cat /var/projects/_review_scratch/s15-readdl-p2-codex/domain_checks_final.guard.stdout | tr -d ' \n' | cut -c1-80; echo
{"exitCode":0,"wallSeconds":0.5744224630179815,"unit":"s11c-guard-b75cd5559f72",

```

```
$ grep -n -E '^(p|r1sq|linear|rem|remprime|total|target|composite_pressure|wholeEOS) = ' /var/projects/_review_scratch/s15-readdl-p2-codex/domain_checks.py
23:p = rho+(rho-1)*s.Abs(rho-1)/2
35:r1sq = R**2+d**2*(1+eps)**2-2*R*d*(1+eps)*mu
36:linear = pbase+(B+eps*Bp)*r1sq*base_derivatives[0]
57:rem = rho0+(z-z0)**2-1-r1base
58:remprime = -2*d*(d-z0)
59:total = s.expand(r1changed+1+rem+eps*remprime)
60:target = rho0+(z-z0-eps*d)**2
62:composite_pressure = s.Abs(s.factor(total-rho0))
64:wholeEOS = s.Abs(rho-rho0)
```

## M3. The same named condition, for the whole minus the expansion with the base moving (Grok stdout, not filed)

```
$ sed -n '/response of whole-minus/,/whole_minus_response_derivative_exists/p' /var/projects/_review_scratch/s15-readdl-p2-grok/03_pressure_domain.stdout
=== response of whole-minus also needs p'' when the base moves ===
smooth whole-minus response 0 matches True
contains p'' term False
p'' piece left -1
p'' piece right 1
full whole-minus left 2
full whole-minus right 0
whole_minus_response_derivative_exists False
```

```
$ sed -n '114,115p;135,138p;147,148p' /var/projects/_review_scratch/s15-readdl-p2-grok/03_pressure_domain.py
def whole_minus(p, pprime, base, amp):
    return p(base + amp) - p(base) - pprime(base) * amp
def p_kink(r):
    return S.Rational(1, 2) * (r - 1) * S.Abs(r - 1)
def pprime_kink(r):
    return S.Abs(r - 1)
def Q_move(shift):
    return whole_minus(p_kink, pprime_kink, 1 + shift, 1)
```

```
$ grep -o -E 'The two pressure sentences are marked .{0,15}' /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md; grep -o -E 'A whole-minus response whose base moves .{0,140}' /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md
The two pressure sentences are marked "For instance".
A whole-minus response whose base moves has to carry differentiability of `p′` at `ρ_2 + ϱ` as well as differentiability of `p` at the total density; a carried profile has to c
```

## M4. Grok on the fold (its section 0) and its F3/F4 reproduction

```
$ sed -n '/^## 0. The fold/,/^## 1\./p' /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md | grep -n -E '^\*\*F[1-4]' | cut -c1-120
5:**F1.** Disposition, what must be true: "nowhere does the list present item 4's conditions as conditions for `E_Q = 0`
9:**F2.** What must be true: "`E_Q^geo` is `𝒟Q` with every change `F′` set to zero, plus `2Q`. `E_Q^resp` is `E_Q �
13:**F3 and F4.** What must be true: each departure is stated where its separation derivative exists and carries that co
```

```
$ sed -n '1,/^=== smooth pressure/p' /var/projects/_review_scratch/s15-readdl-p2-grok/03_pressure_domain.stdout
=== F3: whole pressure defined but not differentiable at the total density ===
p defined 0 left -1 right 1
separation_derivative_exists False
=== whole minus expansion at the same kink ===
left 0 right 2
whole_minus_derivative_exists False
=== smooth pressure: derivative exists and matches the chain rule ===
```

## M5. The legs' reports, scripts and stdout (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-readdl-p2-codex/REVIEW.md /var/projects/_review_scratch/s15-readdl-p2-codex/domain_checks.py /var/projects/_review_scratch/s15-readdl-p2-codex/physics_checks.py /var/projects/_review_scratch/s15-readdl-p2-codex/artifact_audit.py /var/projects/_review_scratch/s15-readdl-p2-codex/_scratch/s11c/domain_checks_final/stdout /var/projects/_review_scratch/s15-readdl-p2-codex/_scratch/s11c/physics_checks_corrected/stdout /var/projects/_review_scratch/s15-readdl-p2-codex/artifact_audit.stdout /var/projects/_review_scratch/s15-readdl-p2-grok/REVIEW.md /var/projects/_review_scratch/s15-readdl-p2-grok/0[1-8]_*.py /var/projects/_review_scratch/s15-readdl-p2-grok/0[1-8]_*.stdout | sed 's|/var/projects/_review_scratch/||'
ef2e26c39472ccb1c903e7189dbce952f2481e1011ddb6e3b200340360021a0a  s15-readdl-p2-codex/REVIEW.md
cc61e1d6fdb053dc35adbb213c24f382b55f5618bd3559ad175911638b1733f5  s15-readdl-p2-codex/domain_checks.py
90e3df2d3ddb7342881ce5b485a109752e4c8e29300d60c76f1c1feac9f93362  s15-readdl-p2-codex/physics_checks.py
c19181990bb8db3852a5366457b5b6ee2422556dfc980da2e958f5b7287baa50  s15-readdl-p2-codex/artifact_audit.py
b3e8364f6b2d604a547440e7c170b2b287aecc7a2a540f65cbebc4c1eaf243a2  s15-readdl-p2-codex/_scratch/s11c/domain_checks_final/stdout
5c97da315913ab4520087cb5445ed04b455ed2552792d2c685db15c9712b353f  s15-readdl-p2-codex/_scratch/s11c/physics_checks_corrected/stdout
fde1761867428120a8e699600d76aab7ccf03ed102b6afeca0ae8f8995d31f7d  s15-readdl-p2-codex/artifact_audit.stdout
49eb48d2453e55dd4dfec472be91144aa4585a1f2bd44e074ca65082e84c4c00  s15-readdl-p2-grok/REVIEW.md
4751631c09ff5f49bcb0e6c3b14a0ff75a70759ae3086d1e7591f4246a412a56  s15-readdl-p2-grok/01_euler_inverse_square.py
5b23e8a4acef2598d5b6f291096b88177b9442174f8b83c63a42c435d3b06b15  s15-readdl-p2-grok/02_parts_and_split.py
f302810bf96009717e9579e4eab57e4081b66fea6996c126f086150ab935faf4  s15-readdl-p2-grok/03_pressure_domain.py
a3fbf8d86706ae685cfc3d4ac53e7b9854beffaa1832bcd22a3a90afde2b9643  s15-readdl-p2-grok/04_family_holds.py
8fa4b675409e553e16a9b6d3b8d37ee5730a6017308ed547292214034af9bb5b  s15-readdl-p2-grok/05_profile_closed_form.py
bc5f42d3e2e71e0782cb6d1f3c5cafde5e09bebe76036e43c57b976c7816143c  s15-readdl-p2-grok/06_velocity_rate_objects.py
2e56de40eb0cc871c9293218e888771922d9a60dae52fa48dae8f07f03deafd6  s15-readdl-p2-grok/07_domain_and_sums.py
a75f2e4aa6033f582f0cf61a92a6b563d959ceed6d4bae1f62b4852ff414ea50  s15-readdl-p2-grok/08_tag_inventory.py
4163f9f55e355ea5e2171b41c6941b61474610717a144f54a24c88059ede9505  s15-readdl-p2-grok/01_euler_inverse_square.stdout
a2fc70f10fb58ed1c920dd92c7235026b6985f2da6eaeff4dcf4debe73bcb404  s15-readdl-p2-grok/02_parts_and_split.stdout
9a3b3a76abd47f85029c513d5694da57fd9305f3e862ccc1f2e1987c81a4a9b4  s15-readdl-p2-grok/03_pressure_domain.stdout
1cd2de2d1f487193a3205b622fb6c95c236b1389ebf4853659e94b694926a3c4  s15-readdl-p2-grok/04_family_holds.stdout
2ef16874d38766d75f96d0a2c61d74d243fbdbd5cc7f6fb6a059410c86f4fe8c  s15-readdl-p2-grok/05_profile_closed_form.stdout
cdd265517045c96f9549028effc215605563ffbb6147d642b6fb35d22954367d  s15-readdl-p2-grok/06_velocity_rate_objects.stdout
335c0c1c4eaf20c72afddbc1a74a50f940f6a8ecdd20e63b3eecac67e8f546ba  s15-readdl-p2-grok/07_domain_and_sums.stdout
686fa8b4a5181b84591ad3188ab29da0d7a12a5a3c0c52039158d9fd70977225  s15-readdl-p2-grok/08_tag_inventory.stdout
```

```
$ for f in /var/projects/_review_scratch/s15-readdl-p2-grok/0[1-8]_*.stderr /var/projects/_review_scratch/s15-readdl-p2-codex/_scratch/s11c/domain_checks_final/stderr; do echo "$(basename $f) $(wc -c < $f)"; done
01_euler_inverse_square.stderr 0
02_parts_and_split.stderr 0
03_pressure_domain.stderr 0
04_family_holds.stderr 0
05_profile_closed_form.stderr 0
06_velocity_rate_objects.stderr 0
07_domain_and_sums.stderr 0
08_tag_inventory.stderr 0
stderr 0
```
