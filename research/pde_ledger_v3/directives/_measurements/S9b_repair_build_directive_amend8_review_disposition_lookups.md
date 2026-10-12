# Measurements — S9b directive amendment 8, review round 0 (generated 2026-10-11 20:55)

Generator: `_scratch/s9b_build/gen/s9b_amend8_r0_lookups.sh` (sha256sum/grep/sed/cat/ls only).

## The round-0 reviewed version
```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend8_review_baseline_r0.sha256
cf258688a4f193d61a67ac7ff2ab72cfc13d2fa07168e753d20b2b55a3f7032d  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md
cf258688a4f193d61a67ac7ff2ab72cfc13d2fa07168e753d20b2b55a3f7032d  _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r0.md
4fba4aca1e0ab26677952195d7cad72d7a5651b2458a03e45476e2ff9d947132  _scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r0.md
```

```
$ grep -n 'that positivity is' _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md
165:   `ρ = ρ₀(1 + f)` is positive, and that positivity is a supplied input of each row whose response reads `ρ`, so it
```

## Both legs ran to completion, and both verdicts
```
$ cat _scratch/s9b_build/amend8_r0_codex.DONE _scratch/s9b_build/amend8_r0_grok.DONE
EXIT=0
EXIT=0
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
_scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_amend8_review_r0_grok.txt:1:Verdict: NEEDS REVISION
```

## Codex's finding, and what it says of item 14
```
$ sed -n 3p _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt
1. **Item 9, lines 164–166:** strict positivity adds a premise the spec does not supply. The independent script constructs a nonnegative density with a local zero and positive asymptotic density. At `n = 1`, the fixed-ratio response remains exactly defined and positive there, but the amendment excludes that stratum. **Fix:** p
```

```
$ grep -n -o 'No additional findings on item 14, coverage, or M2' _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt
5:No additional findings on item 14, coverage, or M2
```

## Grok's findings, and what it says of item 14
```
$ sed -n 9p _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
**1. Item 9 states a sign the spec does not write, and conjoining it changes the Part C domains.**
```

```
$ sed -n 43p _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
**2. As written, "whose response reads `ρ`" does not fix the K9 copies.**
```

```
$ grep -n -o 'This stands only while finding 1.s sign remains' _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
45:This stands only while finding 1's sign remains
```

```
$ grep -n -o '\*\*Item 14\.\*\* The new sentences match C1 and C2' _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
74:**Item 14.** The new sentences match C1 and C2
```

## Grok's script output for finding 1
```
$ grep -n 'DENSITY_GT_0\|RESPONSE_AT_NEG_DENSITY\|ZERO_DENSITY' _scratch/s9b_build/s9b_amend8_review_r0_grok_evidence/amend8_domain_predicates.stdout
6:ZERO_DENSITY_IDENTITY_RESIDUAL: 0.0
7:ZERO_DENSITY_PRODUCT: 0
8:ZERO_DENSITY_PRODUCT_GT_0: False
9:ZERO_DENSITY_RHO0_GT_0: True
10:ZERO_DENSITY_RHO_GT_0: False
21:POWER_RESPONSE_AT_NEG_DENSITY: 1
24:POWER_DENSITY_GT_0: False
27:FIXED_RESPONSE_AT_NEG_DENSITY: 1
31:FIXED_DENSITY_GT_0: False
34:FIXED_LOCAL_CS2_AT_ZERO_DENSITY: 1
35:FIXED_FAR_CS2_AT_ZERO_DENSITY: 1
36:FIXED_RESPONSE_AT_ZERO_DENSITY: 1
38:FIXED_ZERO_DENSITY_GT_0: False
39:FIXED_ZERO_DENSITY_GE_0: True
```

## The two leftover Wolfram units: reported by the legs, stopped and reconciled
```
$ grep -n -o 'supplementary Wolfram run hit the task ceiling' _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt
9:supplementary Wolfram run hit the task ceiling
```

```
$ grep -n -o 'It was not stopped' _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
5:It was not stopped
```

```
$ cat _scratch/s11c/s9b-amend8-review-r0/independent-wolfram-FqQ9xRVa/RECONCILIATION.md
# Reconciliation (orchestrator, 2026-10-11 20:54)

Unit `s11c-guard-462fd9341920` was a supplementary Wolfram run by the Codex leg of the amendment 8 review, round 0. The leg exited (its
`.DONE` marker exists) and reported the run as stuck at the 32-task ceiling, not used as evidence. The unit's owner,
pid 2943850, was dead. Nothing was queued. I stopped the unit with `systemctl --user stop`, confirmed it inactive and no
WolframKernel left, backed the registry up to `resource-pool/reservations.json.stale-amend8-r0-legs.bak`, and wrote
`{}`.
```

```
$ cat _scratch/s11c/s9b-amend8-review-r0/amend8-k9-forms/RECONCILIATION.md
# Reconciliation (orchestrator, 2026-10-11 20:54)

Unit `s11c-guard-4e9c9677fa88` was a supplementary Wolfram run by the Grok leg of the amendment 8 review, round 0. The leg exited (its
`.DONE` marker exists) and reported the run as stuck at the 32-task ceiling, not used as evidence. The unit's owner,
pid 2952936, was dead. Nothing was queued. I stopped the unit with `systemctl --user stop`, confirmed it inactive and no
WolframKernel left, backed the registry up to `resource-pool/reservations.json.stale-amend8-r0-legs.bak`, and wrote
`{}`.
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
0
```


# Round 1 (generated 2026-10-11 21:35; generator `_scratch/s9b_build/gen/s9b_amend8_r1_lookups.sh`)

## The round-1 reviewed version
```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend8_review_baseline_r1.sha256
8127b59a07cca95164594ac351e146416e9aeca9971e06063ebc73435615f651  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r1.md
8127b59a07cca95164594ac351e146416e9aeca9971e06063ebc73435615f651  _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r1.md
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r1.md
9c29891f821d55545909bf85ab5997d419a033aacbaee6847c740a6c6e06b92d  _scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r1.md
```

## Both legs ran to completion, and both verdicts
```
$ cat _scratch/s9b_build/amend8_r1_codex.DONE _scratch/s9b_build/amend8_r1_grok.DONE
EXIT=0
EXIT=0
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend8_review_r1_codex_final.txt _scratch/s9b_build/s9b_amend8_review_r1_grok.txt
_scratch/s9b_build/s9b_amend8_review_r1_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_amend8_review_r1_grok.txt:1:Verdict: NEEDS REVISION
```

```
$ grep -n -o 'Findings:\*\* None within the requested physics scope' _scratch/s9b_build/s9b_amend8_review_r1_codex_final.txt
3:Findings:** None within the requested physics scope
```

## Grok's finding
```
$ sed -n 7p _scratch/s9b_build/s9b_amend8_review_r1_grok.txt | cut -c1-200
**1. Item 9, the ρ ≤ 0 part of the fixed-ratio and power rows.** “Each predicate follows from the row’s response being real and positive” does not pick out one set once a point with ρ = ρ�
```

```
$ cat _scratch/s9b_build/s9b_amend8_review_r1_grok_evidence/density_parts_sympy.stdout
== cs2 derivative vs written formula
0

== dP/drho at n=1
K

== written cs2 at n=1
K/m

== written cs2 at n=1, rho=0
K/m

== pow at n=1 and rho=0 simultaneous
1

== pow substitute n=1 first then rho=0
1

== pow substitute rho=0 first then n=1
1

== S.Zero**S.Zero
1

== S.Zero**0 integer
1

== Pow(0,0,evaluate=False)
0**0

== Pow(0,0) doit
1

== fixed-ratio non-negative root vs power identity
rho=rho0
  nonneg 1
  power  1
  difference 0
rho=rho0/2
  nonneg 2**(1/2 - n/2)
  power  2**(1/2 - n/2)
  difference 0
rho=0
  nonneg sqrt(0**(n - 1)*rho0**(1 - n))
  power  0**(n - 1)
  difference -0**(n - 1) + sqrt(0**(n - 1)*rho0**(1 - n))
rho=-rho0
  nonneg sqrt(rho0**(1 - n)*(-rho0)**(n - 1))
  power  (-1)**(n/2 + 3/2)
  difference (-1)**(n/2)*I + sqrt(rho0**(1 - n)*(-rho0)**(n - 1))
rho=-2*rho0
  nonneg sqrt(rho0**(1 - n)*(-2*rho0)**(n - 1))
  power  (-2)**(n/2 - 1/2)
  difference -(-2)**(n/2 - 1/2) + sqrt(rho0**(1 - n)*(-2*rho0)**(n - 1))

== concrete substitutions, assumptions K>0,m>0,rho0>0,c0>0
n=1 rho=0
  nonneg 1 rewrite 1
  power  1 rewrite 1
  nonneg numeric-like 1
  power  numeric-like 1
n=1 rho=-rho0
  nonneg 1 rewrite 1
  power  1 rewrite 1
  nonneg numeric-like 1
  power  numeric-like 1
n=3 rho=-rho0
  nonneg 1 rewrite 1
  power  -1 rewrite -1
  nonneg numeric-like 1
  power  numeric-like -1
n=3 rho=0
  nonneg 0 rewrite 0
  power  0 rewrite 0
  nonneg numeric-like 0
  power  numeric-like 0
n=2 rho=-rho0
  nonneg I rewrite I
  power  I rewrite I
  nonneg numeric-like I
  power  numeric-like I
n=4 rho=-rho0
  nonneg I rewrite I
  power  -I rewrite -I
  nonneg numeric-like I
  power  numeric-like -I

== positive-symbol substitution trap
== forced powdenest of cs ratio
(f + 1)**(n - 1)

== sqrt of that at f=-2, n=3
1

== forced response
sqrt(rho0**(1 - n)*(rho0*(f + 1))**(n - 1))

f=0 n=3 1 simplified 1
f=-2 n=3 1 simplified 1
f=-1 n=3 0 simplified 0
f=-1 n=1 1 simplified 1
f=-2 n=1 1 simplified 1
f=0 n=1 1 simplified 1

== power response real/positive on rho<=0
s=2 f=-2 value 1 real True positive True doit 1
s=2 f=-1 value 0 real True positive False doit 0
s=2 f=-1/2 value 1/4 real True positive True doit 1/4
s=0 f=-2 value 1 real True positive True doit 1
s=0 f=-1 value 1 real True positive True doit 1
s=1 f=-2 value -1 real True positive False doit -1
s=-1 f=-2 value -1 real True positive False doit -1
s=1/2 f=-2 value I real False positive False doit I
s=3 f=-2 value -1 real True positive False doit -1
s=4 f=-1 value 0 real True positive False doit 0

== first-order coefficient at f=0
constant coef 0
fixed_power coef n/2 - 1/2
fixed_nonneg_sqrt coef n/2 - 1/2
power coef s

== K9 replacements
== K9a
1

== K9b
1

K9a depends on rho False
K9b simplified depends on f False

== domain predicates: response>0 versus with rho0>0 conjunct
equivalent without rho0 conjunct vs with, over reals? skip equivalent
point f=-2 s=2 rho0=1 in pred True
point f=-2 s=2 rho0=1 in pred_local False
point f=-2 s=2 rho0=-1 in pred True
point f=-2 s=2 rho0=-1 in pred_rho0 False
Xor simplified pred vs pred_rho0 (rho0 <= 0) & Contains((f + 1)**s, Reals) & ((f + 1)**s > 0)

== fixed-ratio domain membership at the round-0 points, non-negative c_s
n=1 rho=0 (1, 1, True, True)
n=1 rho=-1 (1, 1, True, True)
n=3 rho=-1 (3, 1, True, True)
n=3 rho=0 (0, 0, True, False)
n=2 rho=-1 (-2, I, False, False)
n=5 rho=-1 (5, 1, True, True)

== n=3 rho<0: sqrt(cs2 ratio) vs (rho/rho0)**1 after refine
sqrt formula Abs(rho)
at rho=-2 2 2
power at n=3 rho=-2 -2
sqrt>0 at rho=-2 True
power>0 at n=3 rho=-2 False
```

```
$ cat _scratch/s9b_build/s9b_amend8_review_r1_grok_evidence/density_parts_wolfram.stdout
== 0^0 and n=1 substitutions ==

                                        0
Power::indet: Indeterminate expression 0  encountered.
Indeterminate

                                        0
Power::indet: Indeterminate expression 0  encountered.
written cs2 n=1 rho=0: Indeterminate
n=1 first then rho=0: K/m

                                        0
Power::indet: Indeterminate expression 0  encountered.

General::stop: Further output of Power::indet will be suppressed during this calculation.
rho=0 first then n=1: Indeterminate
dP at n=1: K
limit n=1 rho->0 of written cs2: K/m
== fixed ratio at concrete points, K=m=rho0=1 ==

Syntax::sntue: Unexpected end of file (probably unfinished expression). (line 58 of "/tmp/s9b_amend8_r1_physics/density_parts_wolfram.wl")
```

```
$ cat _scratch/s9b_build/s9b_amend8_review_r1_grok_evidence/density_parts_wolfram2.stdout
== fixed ratio points ==
n=1 rho=0 nonneg InputForm[Indeterminate] power InputForm[Indeterminate] nonneg>0 False power>0 False
n=1 rho=-1 nonneg InputForm[1] power InputForm[1] nonneg>0 True power>0 True
n=3 rho=-1 nonneg InputForm[1] power InputForm[-1] nonneg>0 True power>0 False
n=3 rho=0 nonneg InputForm[0] power InputForm[0] nonneg>0 False power>0 False

Greater::nord: Invalid comparison with I attempted.

Greater::nord: Invalid comparison with I attempted.
n=2 rho=-1 nonneg InputForm[I] power InputForm[I] nonneg>0 False power>0 False

Greater::nord: Invalid comparison with I attempted.

General::stop: Further output of Greater::nord will be suppressed during this calculation.
n=4 rho=-1 nonneg InputForm[I] power InputForm[-I] nonneg>0 False power>0 False
n=5 rho=-1 nonneg InputForm[1] power InputForm[1] nonneg>0 True power>0 True
n=3 rho=-2 nonneg InputForm[2] power InputForm[-2] nonneg>0 True power>0 False
== power ==
s=2 f=-2 InputForm[1] real True positive True
s=2 f=-1 InputForm[0] real True positive False
s=0 f=-1 InputForm[Indeterminate] real False positive False
s=0 f=-2 InputForm[1] real True positive True
s=1/2 f=-2 InputForm[I] real False positive False
s=4 f=-1 InputForm[0] real True positive False
== K9 ==
```

## What the spec writes for the responses and the sound speed
```
$ sed -n 9,10p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Amendment 2 to Parts A–C** (orchestrator, 2026-10-09; a repair under the user's standing approval). `c_γ` is the
non-negative root of `c_γ²`, and branch existence is stated through `c_γ²` without presupposing its sign. Dispositions:
```

```
$ sed -n 219,220p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  - `P = Kρ^n`, with `n` symbolic, `ρ` the bulk number density and `m` the particle mass (v3 convention). So
    `c_s² = nKρ^(n−1)/m`.
```

```
$ sed -n 244,247p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  - `c_γ(x) ≡ c₀` (`δ(x) ≡ 0`);
  - pointwise fixed ratio: `c_γ(x)/c_s(x) = c₀/c_s0`;
  - symbolic power response: `c_γ(x)/c₀ = (ρ(x)/ρ₀)^s`, with `s` symbolic.

```

## The engines' response forms at 6b6e2ce4
```
$ sed -n 985,987p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    cs_ratio=sp.powsimp(sp.powdenest(cs_squared.subs(rho_bulk,rho0*(1+f_symbol))/cs0_squared,force=True),force=True)
    fixed_ratio_response = sp.sqrt(cs_ratio)  # K9a
    power_response = (1+f_symbol)**s  # K9b
```

```
$ sed -n 74,75p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
responseInputs = {c0, c0 Sqrt[localSoundSquared/(bulkSoundSquared /. rho -> rho0)],
  c0 (localBulkDensity/rho0)^s};
```

## Grok's leftover Wolfram unit: reported, stopped and reconciled
```
$ grep -n -o 'It was not stopped' _scratch/s9b_build/s9b_amend8_review_r1_grok.txt
52:It was not stopped
```

```
$ cat _scratch/s11c/s9b-amend8-review-r1/density-parts-wolfram-2/RECONCILIATION.md
# Reconciliation (orchestrator, 2026-10-11 21:34)

Unit `s11c-guard-41815c328fae` was a Wolfram run by the Grok leg of the amendment 8 review, round 1. The leg exited
(its `.DONE` marker exists) and reported the run as stuck at the 32-task ceiling, using only the lines it printed
before. The unit's owner, pid 2975456, was dead. Nothing was queued. I stopped the unit with `systemctl --user stop`,
confirmed it inactive and no WolframKernel left, backed the registry up to
`resource-pool/reservations.json.stale-amend8-r1-grok.bak`, and wrote `{}`.
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
0
```

