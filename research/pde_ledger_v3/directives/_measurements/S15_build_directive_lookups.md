# Measurements: S15 build brief, leak gate and contract check

Generator: `_scratch/s15/gen/s15_build_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by
the orchestrator.

## G0. The brief and the cleared list

```
$ sha256sum _scratch/s15/build/brief_s15_build.md; git show 91248b67:research/pde_ledger_v3/directives/S15_decision_list.md | sha256sum
8a744115249c46fc3c184718cac6da5559b25f2a6ec91ca8fe5912a99b1114c6  _scratch/s15/build/brief_s15_build.md
ccabfa672e91c771e9f6a07ecbf6649fd33e2139c519cc7b2fadc0412d0fd8aa  -
```

```
$ git log -1 --format='%h %s'
91248b67 Clear the S15 decision list (revision 3, scoped review clear)
```

## G1. The contract is the cleared list's builder-side sections, byte for byte

```
$ diff <(sed -n '/^-\{100\}$/,/^-\{100\}$/p' _scratch/s15/build/brief_s15_build.md | sed '1d;$d' | sed '$d') <(git show 91248b67:research/pde_ledger_v3/directives/S15_decision_list.md | sed -n '/^## Supplied equations/,/^## For the record/p' | sed '$d') && echo 'contract identical to 91248b67 sections Supplied equations .. Model point'
contract identical to 91248b67 sections Supplied equations .. Model point
```

```
$ /usr/bin/grep -n -E '^## ' _scratch/s15/build/brief_s15_build.md
6:## Supplied equations (unfalsifiable in this build)
33:## The drains, their fields and the region
51:## What the build computes
87:## Model point
```

## G2. Outcome words and prior-art names (fixed strings)

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'repel' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'repel' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/r' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/r' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'r²' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'r²' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'r^2' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'r^2' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' _scratch/s15/build/brief_s15_build.md)"
3:You are the author of the S15 build, working in /var/projects/toy_physics. Read AGENTS.md first. Your contract is the text between the two rules below. It is copied byt
118:When both scripts are written, each transcript is produced by a guarded run of the final code. Write a short run record, research/pde_ledger_v3/directives/_measuremen
count: 2
```

```
$ /usr/bin/grep -n -F -- 'drag' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'drag' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'damp' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'damp' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'sign' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'sign' _scratch/s15/build/brief_s15_build.md)"
110:The only place the physical symbols may be combined by hand is in writing the supplied balance with P5's momentum density and carried current and P6's stress as the c
count: 1
```

```
$ /usr/bin/grep -n -F -- 'positive' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'positive' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'negative' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'negative' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Bjerknes' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Bjerknes' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lamb' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lamb' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Newton' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Newton' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Wang' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Wang' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Hicks' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Hicks' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Riemann' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Riemann' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'G7' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'G7' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'G9' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'G9' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'G13' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'G13' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/build/brief_s15_build.md)"
81:The build computes no force law, coefficient, normalization or mass. It gives no relation of any printed object to
count: 1
```

```
$ /usr/bin/grep -n -F -- 'Q1' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q1' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q_1' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q_1' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q₁' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q₁' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q2' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q2' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '4π' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '4π' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'I_F' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'I_F' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'F_12' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'F_12' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Keep it general' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Keep it general' _scratch/s15/build/brief_s15_build.md)"
36:- **Each drain's fields** (the user's "Keep it general") are a radial inflow and a density profile about its live
count: 1
```

```
$ /usr/bin/grep -n -F -- 'Let them move' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Let them move' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Use P5' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Use P5' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

## G3. Symbols in proximity (both on one line), read by the orchestrator

```
$ /usr/bin/grep -n -F 'ρ_1' _scratch/s15/build/brief_s15_build.md | /usr/bin/grep -F 'ρ_2' | cut -c1-170
44:- **The total fields** are `V = V_1 + V_2 + W` and `ρ_br = ρ_1 + ρ_2 + ϱ`. The remainders `W` and `ϱ` are general
62:   **Grouping by drain.** Each field of the split carries one label: drain 1 (`ρ_1`, `V_1`), drain 2 (`ρ_2`, `V_2`)
69:   it, as a labelled choice, its expansion through first order in `ρ_1` about `ρ_2 + ϱ`, grouped as above (a
110:The only place the physical symbols may be combined by hand is in writing the supplied balance with P5's momentum density and carried current and P6's stress as the c
```

```
$ /usr/bin/grep -n -F 'u_1' _scratch/s15/build/brief_s15_build.md | /usr/bin/grep -F 'u_2' | cut -c1-170
110:The only place the physical symbols may be combined by hand is in writing the supplied balance with P5's momentum density and carried current and P6's stress as the c
```

```
$ /usr/bin/grep -n -E 'X_1' _scratch/s15/build/brief_s15_build.md | /usr/bin/grep -E 'X_2' | cut -c1-170
78:Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, and each pressure
110:The only place the physical symbols may be combined by hand is in writing the supplied balance with P5's momentum density and carried current and P6's stress as the c
```

## G4. The builder stays in its lane (no review, agent or downstream instruction)

```
$ /usr/bin/grep -n -i -F -- 'review' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' _scratch/s15/build/brief_s15_build.md)"
1:Status: the S15 decision list has passed its review (commit 91248b67), so you are cleared to PROCEED.
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/build/brief_s15_build.md)"
3:You are the author of the S15 build, working in /var/projects/toy_physics. Read AGENTS.md first. Your contract is the text between the two rules below. It is copied byt
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/build/brief_s15_build.md)"
120:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/build/brief_s15_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/build/brief_s15_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/build/brief_s15_build.md)"
118:When both scripts are written, each transcript is produced by a guarded run of the final code. Write a short run record, research/pde_ledger_v3/directives/_measuremen
count: 1
```
