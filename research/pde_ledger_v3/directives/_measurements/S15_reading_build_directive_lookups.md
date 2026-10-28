# Measurements: S15 reading build brief, leak gate and contract check

Generator: `_scratch/s15/gen/s15_reading_build_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is
a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by
the orchestrator. The brief is assembled by `_scratch/s15/reading_build/assemble.sh` from `head.md`, the two cleared
lists and `foot.md`.

## G0. The brief and the cleared lists

```
$ sha256sum _scratch/s15/reading_build/brief_s15_reading_build.md _scratch/s15/reading_build/head.md _scratch/s15/reading_build/foot.md _scratch/s15/reading_build/assemble.sh
ea806f1995e334c7c17d3635157b69648206e2c6bea52f01f2fb69b05b434642  _scratch/s15/reading_build/brief_s15_reading_build.md
560afecbdc5f61588a9fb1363324598ea4ec8197c53ce6d8e509594edab2451b  _scratch/s15/reading_build/head.md
2bd2ddb6867e696c3dc6a9e8fbc3dff8161664c908488ec07d576c080a3756ed  _scratch/s15/reading_build/foot.md
fe636d19cd8abb6780df2e49dc02ffd7979d019255076c73064851038a052fbe  _scratch/s15/reading_build/assemble.sh
```

```
$ git show 91248b67:research/pde_ledger_v3/directives/S15_decision_list.md | sha256sum; git show 7f52b4d2:research/pde_ledger_v3/directives/S15_reading_decision_list.md | sha256sum; sha256sum research/pde_ledger_v3/directives/S15_decision_list.md research/pde_ledger_v3/directives/S15_reading_decision_list.md
ccabfa672e91c771e9f6a07ecbf6649fd33e2139c519cc7b2fadc0412d0fd8aa  -
8da34e5e672465319c399e1a016382801fbfacdaea2192b13db76e8e27b02acb  -
ccabfa672e91c771e9f6a07ecbf6649fd33e2139c519cc7b2fadc0412d0fd8aa  research/pde_ledger_v3/directives/S15_decision_list.md
8da34e5e672465319c399e1a016382801fbfacdaea2192b13db76e8e27b02acb  research/pde_ledger_v3/directives/S15_reading_decision_list.md
```

```
$ git log -1 --format='%h %s'
7f52b4d2 Clear the S15 reading decision list (scoped pass 3: both CLEAR)
```

## G1. The contract is the two cleared lists' builder-side sections, byte for byte

```
$ diff <(sed -n '/^# Part A$/,/^# Part B$/p' _scratch/s15/reading_build/brief_s15_reading_build.md | sed '1,2d;$d') <(git show 91248b67:research/pde_ledger_v3/directives/S15_decision_list.md | sed -n '/^## Supplied equations/,/^## For the record/p' | sed '$d') && echo 'Part A identical to 91248b67 sections Supplied equations .. Model point'
Part A identical to 91248b67 sections Supplied equations .. Model point
```

```
$ diff <(sed -n '/^# Part B$/,/^-\{100\}$/p' _scratch/s15/reading_build/brief_s15_reading_build.md | sed '1,2d;$d') <(git show 7f52b4d2:research/pde_ledger_v3/directives/S15_reading_decision_list.md | sed -n '/^## Supplied (unfalsifiable/,/^## For the comparator/p' | sed '$d') && echo 'Part B identical to 7f52b4d2 sections Supplied .. Model point'
Part B identical to 7f52b4d2 sections Supplied .. Model point
```

```
$ /usr/bin/grep -n -E '^(# Part|## |-{100}$)' _scratch/s15/reading_build/brief_s15_reading_build.md
9:----------------------------------------------------------------------------------------------------
10:# Part A
12:## Supplied equations (unfalsifiable in this build)
39:## The drains, their fields and the region
57:## What the build computes
93:## Model point
107:# Part B
109:## Supplied (unfalsifiable in this build)
117:## The separation and the family (supplied definitions)
139:## The profile family (a labelled choice, items 3 and 4 only)
144:## What the build computes
176:## Model point
183:----------------------------------------------------------------------------------------------------
```

## G2. Outcome words and prior-art names (fixed strings)

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'repel' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'repel' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/brief_s15_reading_build.md)"
127:- **The departure from the inverse square:** `E_Q = 𝒟Q + 2Q`. `Q` falls off as the inverse square of the separation
count: 1
```

```
$ /usr/bin/grep -n -F -- 'falls off' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'falls off' _scratch/s15/reading_build/brief_s15_reading_build.md)"
127:- **The departure from the inverse square:** `E_Q = 𝒟Q + 2Q`. `Q` falls off as the inverse square of the separation
count: 1
```

```
$ /usr/bin/grep -n -F -- '1/r' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/r' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'r²' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'r²' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'r^2' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'r^2' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/brief_s15_reading_build.md)"
120:  - **Held** (labelled): the region and its motion (`X_2`, `Ẋ_2`, `r_in`, `R` and their rates), and the form of
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= -2' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= -2' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'vanish' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'vanish' _scratch/s15/reading_build/brief_s15_reading_build.md)"
161:   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
count: 1
```

```
$ /usr/bin/grep -n -F -- 'zero' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'zero' _scratch/s15/reading_build/brief_s15_reading_build.md)"
52:  split of the total fields, not a superposition law. Nothing makes the remainders small or zero.
193:Keep every field the contract calls live as a live function of its arguments: `u_1`, `u_2`, `ρ_1`, `ρ_2`, `W`, `ϱ`, `p_br`, `η`, `ζ`, `Π_o`, `f`, and each chang
count: 2
```

```
$ /usr/bin/grep -n -F -- 'drag' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'drag' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'damp' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'damp' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/brief_s15_reading_build.md)"
174:normalization or mass, and identifies no printed object with a pull or with the force on a throat.
count: 1
```

```
$ /usr/bin/grep -n -F -- 'sign' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'sign' _scratch/s15/reading_build/brief_s15_reading_build.md)"
191:The only place the physical symbols may be combined by hand is in writing: the supplied balance with P5's momentum density and carried current and P6's stress, as Par
count: 1
```

```
$ /usr/bin/grep -n -F -- 'positive' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'positive' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'negative' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'negative' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Bjerknes' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Bjerknes' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lamb' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lamb' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Newton' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Newton' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Wang' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Wang' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Hicks' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Hicks' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Riemann' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Riemann' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'G7' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'G7' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'G9' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'G9' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'G13' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'G13' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'requirement' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'requirement' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'S21' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'S21' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/brief_s15_reading_build.md)"
87:The build computes no force law, coefficient, normalization or mass. It gives no relation of any printed object to
174:normalization or mass, and identifies no printed object with a pull or with the force on a throat.
count: 2
```

```
$ /usr/bin/grep -n -F -- 'Q1' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q1' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q_1' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q_1' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Q₁' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Q₁' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '4π' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '4π' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'F_12' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'F_12' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Keep it general' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Keep it general' _scratch/s15/reading_build/brief_s15_reading_build.md)"
42:- **Each drain's fields** (the user's "Keep it general") are a radial inflow and a density profile about its live
count: 1
```

```
$ /usr/bin/grep -n -F -- 'Let them move' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Let them move' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Use P5' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Use P5' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Reading build' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Reading build' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

## G3. Symbols in proximity (both on one line), read by the orchestrator

```
$ /usr/bin/grep -n -F 'ρ_1' _scratch/s15/reading_build/brief_s15_reading_build.md | /usr/bin/grep -F 'ρ_2' | cut -c1-170
50:- **The total fields** are `V = V_1 + V_2 + W` and `ρ_br = ρ_1 + ρ_2 + ϱ`. The remainders `W` and `ϱ` are general
68:   **Grouping by drain.** Each field of the split carries one label: drain 1 (`ρ_1`, `V_1`), drain 2 (`ρ_2`, `V_2`)
75:   it, as a labelled choice, its expansion through first order in `ρ_1` about `ρ_2 + ϱ`, grouped as above (a
193:Keep every field the contract calls live as a live function of its arguments: `u_1`, `u_2`, `ρ_1`, `ρ_2`, `W`, `ϱ`, `p_br`, `η`, `ζ`, `Π_o`, `f`, and each chang
```

```
$ /usr/bin/grep -n -F 'u_1' _scratch/s15/reading_build/brief_s15_reading_build.md | /usr/bin/grep -F 'u_2' | cut -c1-170
193:Keep every field the contract calls live as a live function of its arguments: `u_1`, `u_2`, `ρ_1`, `ρ_2`, `W`, `ϱ`, `p_br`, `η`, `ζ`, `Π_o`, `f`, and each chang
```

```
$ /usr/bin/grep -n -E 'X_1' _scratch/s15/reading_build/brief_s15_reading_build.md | /usr/bin/grep -E 'X_2' | cut -c1-170
84:Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, and each pressure
118:- **The separation:** `D = X_1 − X_2`, `d = |D|`. The domain requires `d > R` (L2).
193:Keep every field the contract calls live as a live function of its arguments: `u_1`, `u_2`, `ρ_1`, `ρ_2`, `W`, `ϱ`, `p_br`, `η`, `ζ`, `Π_o`, `f`, and each chang
```

```
$ /usr/bin/grep -n -E 'a_k|b_k' _scratch/s15/reading_build/brief_s15_reading_build.md | /usr/bin/grep -E '[=<>]' | cut -c1-170
140:`u_k(r, t) = A_k(t) r^{a_k}` and `ρ_k(r, t) = B_k(t) r^{b_k}`, for `k = 1, 2`. The amplitudes `A_k`, `B_k` are general
```

```
$ /usr/bin/grep -n -F 'E_Q' _scratch/s15/reading_build/brief_s15_reading_build.md | /usr/bin/grep -E '= *0|vanish' | cut -c1-170
128:  along the family exactly where `E_Q = 0` for every `D` on the domain.
130:  - `E_Q^geo = (𝒟Q)|_{F′=0} + 2Q` is the departure when no quantity other than `X_1` changes along the family;
```

## G4. The builder stays in its lane (no review, agent or downstream instruction)

```
$ /usr/bin/grep -n -i -F -- 'review' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' _scratch/s15/reading_build/brief_s15_reading_build.md)"
1:Status: the S15 reading decision list has passed its review (commit 7f52b4d2), so you are cleared to PROCEED.
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/brief_s15_reading_build.md)"
3:You are the author of the S15 reading build, working in /var/projects/toy_physics. Read AGENTS.md first. Your contract is the text between the two rules below, in two p
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/brief_s15_reading_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/brief_s15_reading_build.md)"
203:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/brief_s15_reading_build.md)"
113:  printed them (`3192123b`), and the accepted comparator compared them (`5a019394`) (L3).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/brief_s15_reading_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/brief_s15_reading_build.md)"
201:When both scripts are written, each transcript is produced by a guarded run of the final code. Write a short run record, research/pde_ledger_v3/directives/_measuremen
count: 1
```
