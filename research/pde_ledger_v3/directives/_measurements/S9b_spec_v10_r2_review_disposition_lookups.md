# Measurements — S9b spec v10 review round 2 dispositions (generated 2026-10-08 13:19)

Generator: `_scratch/s9b_build/gen/s9b_spec_v10_r2_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed files are the working copies checked against the baseline below.

```
$ sha256sum -c _scratch/s9b_build/s9b_spec_v10_review_baseline_r1.sha256
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_v10_spec_lookups.md: OK
```

```
$ grep -n -o 'Verdict:[* ]*[A-Z][A-Z ]*' _scratch/s9b_build/s9b_spec_v10_review_r2_claude.md _scratch/s9b_build/s9b_spec_v10_review_r2_grok.txt
_scratch/s9b_build/s9b_spec_v10_review_r2_claude.md:11:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_spec_v10_review_r2_grok.txt:1:Verdict: CLEAR
```

## E1 — the induced-measure rule denies a bound that the slope counting supplies for one reading
```
$ sed -n '68,73p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. Any induced-measure claim
  names its reading: the same densities re-expressed per induced volume, or a mass law imposed on
  the induced measure. It keeps `∂_r[(∂ξ_w)²]` live and attaches no relative order to the correction
  to `j_n` unless it states a condition bounding that correction relative to `j_n` itself (amended
  D4). No such bound is supplied here. Neither the supplied slope counting nor the historical
  gradient-scale condition `∂_r[(∂ξ_w)²] = O(ε/r)` supplies that bound. The latter is not imposed.
```

```
$ grep -n 'No source supplies that order\|unless it states a condition' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
61:  correction to the implied `j_n`". No source supplies that order. The gradient-scale condition
66:  unless it states a condition that bounds that correction relative to `j_n` itself.
```

```
$ grep -c 'slope counting' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
0
```

```
$ sed -n '478,480p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
The mass-law qualification travels explicitly: `∇·(ρ_br V)=−j_n` is a **coordinate-`d³x` measure**
input. For a claim reading `j_n` or `ρ_br` per induced measure, or comparing an induced-metric mass
law, the recorded qualification is a **relative `O(ε)` correction to `j_n`**, not a live O6 law;
```

```
$ grep -n '(∂ξ_w)² = O(ε)' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
132:  δ = O(ε),      (∂ξ_w)² = O(ε),      V/c₀ = O(ε^{1/2})   ⇒   (V/c₀)² = O(ε) .
```

The leg's script output:
```
$ grep -n 'reading 1\|reading 2 residual\|limit B->0' _scratch/s9b_build/s9b_spec_v10_review_r2_claude_evidence/induced_measure_readings.child_stdout.txt
2:reading 1: (j_1 - j_c)/j_c = -1 + 1/sqrt(Derivative(xi_w(r), r)**2 + 1)
3:reading 1 relative correction as series in s=(dxi)^2: -s/2 + 3*s**2/8 + O(s**3)
4:reading 1 relative correction free symbols/functions: ['Derivative(xi_w(r), r)', 'xi_w(r)']
6:reading 2 residual vs -(1/2) rho V d_r[(dxi)^2]/det g : 0
9:witness reading 1 relative correction = sqrt(r)/sqrt(r + sigma) - 1
11:witness reading 2 relative correction, limit B->0+ : oo
12:witness reading 1 relative correction, limit B->0+ : -1 + 1/sqrt(1 + sigma/r)
```

## E2 — `U_material^w` is undefined; read as the material velocity at the faces it would close the drain's normal outflow
```
$ grep -n 'U_material' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
387:ξ_w(r) ≡ 0 ,      U_material^w(r) ≡ 0 .
```

```
$ grep -o 'So `ξ_w` and the material `w`-velocity vanish there' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
So `ξ_w` and the material `w`-velocity vanish there
```

```
$ sed -n '71,77p' research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
follows from those names (C §§2, 5–7). On the steady supplied graph `w = ξ_w` (§3.1), a brane
material point with in-plane coordinate velocity `V^i` has the bulk-direction coordinate velocity
that the graph and `V` determine. Construct that component from them; it stays live through `V` and
`ξ_w`, and there is no separate bulk-direction velocity operand. Normal velocity content that the
graph does not determine belongs to `𝒩_br^live` and `𝒥_map`. One example is different local
velocities at the faces of a finite-thickness realization, whose centre and thickness content `ξ_w`
does not fix. Such content stays an unevaluated action with those operands named. This material
```

```
$ grep -o 'graph material velocity with `U^w = V_r ξ′`' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
graph material velocity with `U^w = V_r ξ′`
```

The leg's script output:
```
$ grep -n 'u^w at\|net outward normal\|centre-graph' _scratch/s9b_build/s9b_spec_v10_review_r2_claude_evidence/p6_wparity_refs.child_stdout.txt
11:u^w at centre w=0: 0
12:u^w at faces w=+L/2, -L/2: L*q(x)/2 , -L*q(x)/2
13:net outward normal mass flux through both faces (per area): L*q(x)*rho0(x)
14:centre-graph w-velocity V^i d_i xi_w with xi_w=0: 0
```

## E3 — author and orchestrator process text in the engine-facing spec
```
$ grep -n 'Authoring STOP\|Spec review\.\*\*\|No CAS, build' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
513:- **Spec review.** A fresh non-author Claude agent and Grok review v10 until clear, after this
539:**Authoring STOP.** Repair v10 in place and report changes for repair items 1–4, source conflicts
540:and missing sourced pieces. No CAS, build, review launch, commit, push or spawned agent is part of this task.
```

```
$ grep -n 'imports nothing' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
512:- **Engines.** SymPy, plus a blind Wolfram engine that imports nothing. No Lean (CLAUDE.md L5).
```

