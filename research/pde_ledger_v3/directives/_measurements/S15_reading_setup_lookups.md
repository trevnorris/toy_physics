# Measurements: S15 reading decision list, setup lookups

Generator: `_scratch/s15/gen/s15_reading_setup_lookups.sh` (GNU grep, sed, cut, git, sha256sum, `LC_ALL=C`), run from
`research/pde_ledger_v3/`. Each block is a command and its literal output.

## L1. The user's S15 decisions, verbatim (pass 1 disposition)

```
$ sed -n '105,116p' directives/_measurements/S15_decision_list_review_p1_disposition.md
## The user's decisions (2026-10-24), verbatim
The user selected these options in conversation. The wording is the orchestrator's option text, as revision 1 of
the list carried it (`:34–42`). The record quotes them from here. The build directive carries the choices only.
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

## L2. The S15 decision list: revision, the reading it leaves to the record, and its sections

```
$ git log --format='%h %s' -- directives/S15_decision_list.md | head -6
e393430a Repair the S15 decision list for pass 2 (revision 3, not reviewed)
cd853e16 Fold the S15 decision list (G2: one pass, folded once)
d6e4f74a Preserve the S15 decision list as reviewed in pass 1 (not accepted)
```

```
$ git diff --stat 91248b67 -- directives/S15_decision_list.md | tail -1; echo "diff lines: $(git diff 91248b67 -- directives/S15_decision_list.md | wc -l)"
diff lines: 0
```

```
$ sed -n '45,46p' directives/S15_decision_list.md
Which profiles give which dependence on the separation is the record's reading of the printed objects. The build
prints the objects.
```

```
$ grep -n -E '^## ' directives/S15_decision_list.md
14:## Filing and framing
35:## User decisions (2026-10-24)
48:## Supplied equations (unfalsifiable in this build)
75:## The drains, their fields and the region
93:## What the build computes
129:## Model point
143:## For the record and its register pass (not in the build directive)
174:## Retained OPEN
188:## Intended claim
193:## Next
```

```
$ sed -n '120,127p' directives/S15_decision_list.md
Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, and each pressure
object where its condition above holds.

The build computes no force law, coefficient, normalization or mass. It gives no relation of any printed object to
S14a's `J` or to the plan's `Q`, `Θ_Q` or `N_∞,3`. It identifies no printed object with the force on a throat.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value.
```

```
$ sed -n '145,156p' directives/S15_decision_list.md
- **What the record reads from the printed objects.** From the printed objects alone, it states:
  - which conditions on the profiles, the remainders, `p_br`, `η`, `ζ`, `Π_o` and `f` give a named printed piece, or
    a named sum of printed pieces, a stated dependence on the separation, including which give an inverse-square one,
    as the user's "Keep it general" asks. Each such dependence is stated with its conditions on every live quantity
    the piece contains, and on the surfaces: both radii, their rates, and the motion of `X_2`;
  - the terms that depend on each drain's velocity;
  - the sign of a printed piece only on its stated conditions, and only for that piece (`:666–668`).
- **Requirements.** The record files each requirement as the register's rest-on test admits it, or on a user decision
  that it quotes. The user's "Keep it general" files gravity's flow requirement for S21's comparison, as a condition
  on the named object, with all its conditions, those on `p_br`, `η` and `ζ` included. It is not called a pull or the
  force on a throat (Retained OPEN). It merges into `R-S21-01` where the register's merge test says the two cover the
  same requirement.
```

## L3. The accepted S15 build and comparator

```
$ git log -1 --format='%h %s' 3192123b; git log -1 --format='%h %s' 5a019394
3192123b Accept the S15 build (review round 1 clear)
5a019394 Accept the S15 comparator (round 1) and file its production output
```

```
$ sha256sum scripts/out/S15_two_drain_region.out mathematica/out/S15_two_drain_region.out scripts/S15_two_drain_region.py mathematica/S15_two_drain_region.wl
4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c  scripts/out/S15_two_drain_region.out
811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6  mathematica/out/S15_two_drain_region.out
55f1716333c04b64ce0ea2745463cf180ceef3d692f2b01ddf769570611986dc  scripts/S15_two_drain_region.py
37fa7ec54b9c91c4a59b364f60f629c36180ac34f5836e8a3b7e0b9996d030ec  mathematica/S15_two_drain_region.wl
```

```
$ grep -n -o -F 'COMPARISON_IDENTICAL_TO_REVIEWED_RUN' directives/_measurements/S15_comparator_production_run.md | tail -1
121:COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```

## L4. The printed object families (tag counts per engine)

```
$ for s in inner outer; do for p in carried viscous transport_centre transport_radius; do echo "$s/$p py $(grep -c -E "^$s/$p/component_[123]/labels_[a-z0-9+]+: " scripts/out/S15_two_drain_region.out) wl $(grep -c -E "^$s/$p/component_[123]/labels_\\{[^}]*\\}: " mathematica/out/S15_two_drain_region.out)"; done; done
inner/carried py 21 wl 21
inner/viscous py 21 wl 21
inner/transport_centre py 21 wl 21
inner/transport_radius py 21 wl 21
outer/carried py 21 wl 21
outer/viscous py 21 wl 21
outer/transport_centre py 21 wl 21
outer/transport_radius py 21 wl 21
```

```
$ for s in inner outer; do echo "$s pressure whole py $(grep -c -E "^$s/pressure/component_[123]: " scripts/out/S15_two_drain_region.out) wl $(grep -c -E "^$s/pressure/component_[123]: " mathematica/out/S15_two_drain_region.out); expansion groups py $(grep -c -E "^$s/pressure/component_[123]/expansion_labels_[a-z0-9+]+: " scripts/out/S15_two_drain_region.out) wl $(grep -c -E "^$s/pressure/component_[123]/expansion_labels_\\{[^}]*\\}: " mathematica/out/S15_two_drain_region.out); other py $(grep -c -E "^$s/other/component_[123]: " scripts/out/S15_two_drain_region.out) wl $(grep -c -E "^$s/other/component_[123]: " mathematica/out/S15_two_drain_region.out)"; done
inner pressure whole py 3 wl 3; expansion groups py 21 wl 21; other py 3 wl 3
outer pressure whole py 3 wl 3; expansion groups py 21 wl 21; other py 3 wl 3
```

```
$ grep -o -E '^inner/carried/component_1/labels_[a-z0-9+]+: ' scripts/out/S15_two_drain_region.out; grep -o -E '^inner/carried/component_1/labels_\{[^}]*\}: ' mathematica/out/S15_two_drain_region.out
inner/carried/component_1/labels_drain1: 
inner/carried/component_1/labels_drain2: 
inner/carried/component_1/labels_remainder: 
inner/carried/component_1/labels_drain1+drain2: 
inner/carried/component_1/labels_drain1+remainder: 
inner/carried/component_1/labels_drain2+remainder: 
inner/carried/component_1/labels_drain1+drain2+remainder: 
inner/carried/component_1/labels_{drain1}: 
inner/carried/component_1/labels_{drain2}: 
inner/carried/component_1/labels_{remainder}: 
inner/carried/component_1/labels_{drain1, drain2}: 
inner/carried/component_1/labels_{drain1, remainder}: 
inner/carried/component_1/labels_{drain2, remainder}: 
inner/carried/component_1/labels_{drain1, drain2, remainder}: 
```

```
$ grep -o -E '^region/component_1/[a-z0-9_+]+: ' scripts/out/S15_two_drain_region.out | sort -u | tr '\n' ' '; echo; grep -o -E '^region/component_1/[a-z_]+(\{[^}]*\})?: ' mathematica/out/S15_two_drain_region.out | sort -u | tr '\n' ' '
region/component_1/balance_dimension_operand_a:  region/component_1/balance_dimension_operand_b:  region/component_1/balance_dimension_residuals:  region/component_1/balance_equation:  region/component_1/balance_parts:  region/component_1/integrated_force:  region/component_1/momentum:  region/component_1/momentum_labels_drain1+drain2+remainder:  region/component_1/momentum_labels_drain1+drain2:  region/component_1/momentum_labels_drain1+remainder:  region/component_1/momentum_labels_drain1:  region/component_1/momentum_labels_drain2+remainder:  region/component_1/momentum_labels_drain2:  region/component_1/momentum_labels_remainder:  region/component_1/momentum_rate:  region/component_1/rate_labels_drain1+drain2+remainder:  region/component_1/rate_labels_drain1+drain2:  region/component_1/rate_labels_drain1+remainder:  region/component_1/rate_labels_drain1:  region/component_1/rate_labels_drain2+remainder:  region/component_1/rate_labels_drain2:  region/component_1/rate_labels_remainder:  
region/component_1/balance_dimension_operand_a:  region/component_1/balance_dimension_operand_b:  region/component_1/balance_dimension_residuals:  region/component_1/balance_equation:  region/component_1/balance_parts:  region/component_1/integrated_force:  region/component_1/momentum:  region/component_1/momentum_labels_{drain1, drain2, remainder}:  region/component_1/momentum_labels_{drain1, drain2}:  region/component_1/momentum_labels_{drain1, remainder}:  region/component_1/momentum_labels_{drain1}:  region/component_1/momentum_labels_{drain2, remainder}:  region/component_1/momentum_labels_{drain2}:  region/component_1/momentum_labels_{remainder}:  region/component_1/momentum_rate:  region/component_1/rate_labels_{drain1, drain2, remainder}:  region/component_1/rate_labels_{drain1, drain2}:  region/component_1/rate_labels_{drain1, remainder}:  region/component_1/rate_labels_{drain1}:  region/component_1/rate_labels_{drain2, remainder}:  region/component_1/rate_labels_{drain2}:  region/component_1/rate_labels_{remainder}:  ```

## L5. Where the drains' centres enter: the charts and the radial distances

```
$ grep -E '^(inner/chart|outer/chart|region/chart|distances): ' scripts/out/S15_two_drain_region.out | cut -c1-260
distances: [sqrt((x_1 - X_1_1(t))**2 + (x_2 - X_1_2(t))**2 + (x_3 - X_1_3(t))**2), sqrt((x_1 - X_2_1(t))**2 + (x_2 - X_2_2(t))**2 + (x_3 - X_2_3(t))**2)]
inner/chart: Matrix([[X_2_1(t) + r_in(t)*sin(theta)*cos(phi)], [X_2_2(t) + r_in(t)*sin(phi)*sin(theta)], [X_2_3(t) + r_in(t)*cos(theta)]])
outer/chart: Matrix([[R(t)*sin(theta)*cos(phi) + X_2_1(t)], [R(t)*sin(phi)*sin(theta) + X_2_2(t)], [R(t)*cos(theta) + X_2_3(t)]])
region/chart: Matrix([[s*sin(theta)*cos(phi) + X_2_1(t)], [s*sin(phi)*sin(theta) + X_2_2(t)], [s*cos(theta) + X_2_3(t)]])
```

```
$ grep -E '^(inner/chart|outer/chart|region/chart|radial_distances): ' mathematica/out/S15_two_drain_region.out | cut -c1-260
radial_distances: {Sqrt[(x - xx[1, 1][t])^2 + (y - xx[1, 2][t])^2 + (z - xx[1, 3][t])^2], Sqrt[(x - xx[2, 1][t])^2 + (y - xx[2, 2][t])^2 + (z - xx[2, 3][t])^2]}
inner/chart: {Cos[phi]*rin[t]*Sin[theta] + xx[2, 1][t], rin[t]*Sin[phi]*Sin[theta] + xx[2, 2][t], Cos[theta]*rin[t] + xx[2, 3][t]}
outer/chart: {Cos[phi]*rout[t]*Sin[theta] + xx[2, 1][t], rout[t]*Sin[phi]*Sin[theta] + xx[2, 2][t], Cos[theta]*rout[t] + xx[2, 3][t]}
region/chart: {s*Cos[phi]*Sin[theta] + xx[2, 1][t], s*Sin[phi]*Sin[theta] + xx[2, 2][t], s*Cos[theta] + xx[2, 3][t]}
```

## L6. The register rows the reading bears on

```
$ grep -n -E '^### R-S21-01|^### R-S14a-01|^\| \*\*S21\*\*' SUBSTRATE_REQUIREMENTS.md | cut -c1-160
117:| **S21** | `R-S21-01` | conditional deflection and radar-slope profile requirements, with Part C rewrites and implied exchange; comparison of force sectors
714:### R-S21-01 — conditional deflection and radar-slope profile requirements
974:### R-S14a-01 — the dynamically selected slow limit and far behaviour
```
