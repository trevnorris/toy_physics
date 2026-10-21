# Measurements: S14a part 1 comparator decision list, setup lookups (generated 2026-10-21 10:06)

Generator: `_scratch/s14a/gen/cmp_setup_lookups.sh` (grep/sed/sort/comm/wc/git/sha256sum only; this file is
written by it). Mechanical lookups only: existence, literal-match counts, names.

## L1. The accepted transcripts and engine sources at 94fddf2f

```
$ cd /var/projects/toy_physics && git log -1 --format='%h %s' 94fddf2f
94fddf2f Accept the S14a part 1 build (cleared in scoped round 2)
```

```
$ cd /var/projects/toy_physics && sha256sum research/pde_ledger_v3/scripts/out/S14a_1_projection.out research/pde_ledger_v3/mathematica/out/S14a_1_projection.out research/pde_ledger_v3/scripts/S14a_1_projection.py research/pde_ledger_v3/mathematica/S14a_1_projection.wl
3bfb17aa5205e7f0ddec1eea44c920d4f0bb4e3dffe4ccf6fcc74abe88fbe0db  research/pde_ledger_v3/scripts/out/S14a_1_projection.out
7911ff19f667a388a2a54a3eb0d076ed06143d791ddb2b5bccfe397588652d5b  research/pde_ledger_v3/mathematica/out/S14a_1_projection.out
3d70f411dc894f22dc416b5cae7f7a09cbf8a3603b92a9a8e70ac7639da3e2f9  research/pde_ledger_v3/scripts/S14a_1_projection.py
365be5d41d3c1f34ddb727bef8a7b49b02bc39e2122a85644aa3dcd2a5870905  research/pde_ledger_v3/mathematica/S14a_1_projection.wl
```

```
$ cd /var/projects/toy_physics && git diff --stat 94fddf2f -- research/pde_ledger_v3/scripts/S14a_1_projection.py research/pde_ledger_v3/mathematica/S14a_1_projection.wl research/pde_ledger_v3/scripts/out/S14a_1_projection.out research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
0
```

## L2. Tag names: count per stream, repeats, names shared by both streams

A tag name is the text before the first ": " on a line; `WORKING_ON` lines are excluded.

```
$ for e in scripts mathematica; do printf '%s tags: ' $e; grep -v '^WORKING_ON' /var/projects/toy_physics/research/pde_ledger_v3/$e/out/S14a_1_projection.out | grep -E '^[A-Za-z][^:]*: ' | sed 's/: .*//' | wc -l; done
scripts tags: 600
mathematica tags: 383
```

```
$ for e in scripts mathematica; do printf '%s repeated names: ' $e; grep -v '^WORKING_ON' /var/projects/toy_physics/research/pde_ledger_v3/$e/out/S14a_1_projection.out | grep -E '^[A-Za-z][^:]*: ' | sed 's/: .*//' | sort | uniq -d | wc -l; done
scripts repeated names: 0
mathematica repeated names: 0
```

```
$ comm -12 <(grep -v ^WORKING_ON /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | grep -E '^[A-Za-z][^:]*: ' | sed 's/: .*//' | sort -u) <(grep -v ^WORKING_ON /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | grep -E '^[A-Za-z][^:]*: ' | sed 's/: .*//' | sort -u) | wc -l
0
```

## L3. Lines that do not start with a tag (payloads that span lines)

```
$ grep -cvE '^[A-Za-z_][^:]*: ' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out
6
```

```
$ grep -nvE '^[A-Za-z_][^:]*: ' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | cut -c1-60
164:[radius*sin(theta)*cos(phi)],
165:[radius*sin(phi)*sin(theta)],
166:[         radius*cos(theta)]])
169:[1,         0,                       0],
170:[0, radius**2,                       0],
171:[0,         0, radius**2*sin(theta)**2]])
```

```
$ grep -cvE '^[A-Za-z_][^:]*: ' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out
0
```

## L4. Notation heads, literal-match counts per stream

```
$ printf 'SymPy %s ' 'Integral('; grep -oF 'Integral(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Integral('; grep -oF 'Integral(' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Integral( 567
Wolfram Integral( 0
```

```
$ printf 'SymPy %s ' 'Inactive[Integrate]'; grep -oF 'Inactive[Integrate]' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Inactive[Integrate]'; grep -oF 'Inactive[Integrate]' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Inactive[Integrate] 0
Wolfram Inactive[Integrate] 790
```

```
$ printf 'SymPy %s ' 'Inactive[Equal]'; grep -oF 'Inactive[Equal]' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Inactive[Equal]'; grep -oF 'Inactive[Equal]' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Inactive[Equal] 0
Wolfram Inactive[Equal] 12
```

```
$ printf 'SymPy %s ' 'Derivative('; grep -oF 'Derivative(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Derivative('; grep -oF 'Derivative(' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Derivative( 1690
Wolfram Derivative( 0
```

```
$ printf 'SymPy %s ' 'Derivative['; grep -oF 'Derivative[' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Derivative['; grep -oF 'Derivative[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Derivative[ 0
Wolfram Derivative[ 2794
```

```
$ printf 'SymPy %s ' 'DiracDelta('; grep -oF 'DiracDelta(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'DiracDelta('; grep -oF 'DiracDelta(' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy DiracDelta( 72
Wolfram DiracDelta( 0
```

```
$ printf 'SymPy %s ' 'DiracDelta['; grep -oF 'DiracDelta[' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'DiracDelta['; grep -oF 'DiracDelta[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy DiracDelta[ 0
Wolfram DiracDelta[ 112
```

```
$ printf 'SymPy %s ' 'Heaviside('; grep -oF 'Heaviside(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Heaviside('; grep -oF 'Heaviside(' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Heaviside( 202
Wolfram Heaviside( 0
```

```
$ printf 'SymPy %s ' 'HeavisideTheta['; grep -oF 'HeavisideTheta[' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'HeavisideTheta['; grep -oF 'HeavisideTheta[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy HeavisideTheta[ 0
Wolfram HeavisideTheta[ 142
```

```
$ printf 'SymPy %s ' 'Ne('; grep -oF 'Ne(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Ne('; grep -oF 'Ne(' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Ne( 5
Wolfram Ne( 0
```

```
$ printf 'SymPy %s ' 'ConditionalExpression['; grep -oF 'ConditionalExpression[' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'ConditionalExpression['; grep -oF 'ConditionalExpression[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy ConditionalExpression[ 0
Wolfram ConditionalExpression[ 4
```

```
$ printf 'SymPy %s ' 'Matrix('; grep -oF 'Matrix(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | wc -l; printf 'Wolfram %s ' 'Matrix('; grep -oF 'Matrix(' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | wc -l
SymPy Matrix( 2
Wolfram Matrix( 0
```

## L5. How each stream prints a condition on the projected velocity (first 60 characters of each line)

```
$ grep -nF 'Ne(' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | cut -c1-60
360:SphericalSpecialization/projected velocity: {domain: Ne(
362:SphericalSpecialization/outer projected velocity: {domai
364:SphericalSpecialization/outer volume flux: {domain: Ne(I
614:SphericalSpecialization/dimensions/projected velocity/di
617:SphericalSpecialization/dimensions/volume flux/dimension
```

```
$ grep -nF 'ConditionalExpression[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | cut -c1-90
283:spherical specialization / shell/projected velocity: ConditionalExpression[Inactive[In
285:spherical specialization / shell/projected velocity volume flux: ConditionalExpression
394:spherical specialization / dimensions/projected velocity: ConditionalExpression[{0, 1,
396:spherical specialization / dimensions/projected velocity volume flux: ConditionalExpre
```

## L6. Earlier comparators, as the decision list cites them

```
$ cd /var/projects/toy_physics && git log -1 --format='%h %s' fe485b73 && git ls-tree --name-only fe485b73 research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py
fe485b73 Accept the S12 comparator (scoped round 2, with the engine-row note) and file its production output
research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py
research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py
```

```
$ cd /var/projects/toy_physics && git log -1 --format='%h %s' 88bb1f41 && git ls-tree --name-only 88bb1f41 research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py
88bb1f41 Accept the S9b comparator (cleared in round 4) and file its production output
research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py
```

```
$ cd /var/projects/toy_physics && ls research/pde_ledger_v3/directives/_measurements/ | grep -E '^S12_comparator_build_review_r[0-2]_disposition.md$'
S12_comparator_build_review_r0_disposition.md
S12_comparator_build_review_r1_disposition.md
S12_comparator_build_review_r2_disposition.md
```

```
$ cd /var/projects/toy_physics && grep -nE '^\*\*(F[1-4]|G[12]|J1)\. ' research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r{0,1,2}_disposition.md | cut -c1-200
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r0_disposition.md:23:**F1. The same derivative, printed in the two engines' notations, is read as two different objects.** (
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r0_disposition.md:47:**F2. Distinct live functions can share a test function, so evaluation cannot tell them apart.** (Claud
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r0_disposition.md:63:**F3. No control goes through the current rows' component pairing.** (Claude F3.) **CONFIRMED.**
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r0_disposition.md:76:**F4. Dimension-ledger terms are bound by being identical, then compared with each other and counted as
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r1_disposition.md:32:**G1. A label that is a list of strings is read as raw text, so a notation difference is printed as a t
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r1_disposition.md:50:**G2. The changed-label control does not go through the production label path.** (Claude 2.) **CONFIRME
research/pde_ledger_v3/directives/_measurements/S12_comparator_build_review_r2_disposition.md:41:**J1. The repaired `engine` row pairs components that the engines build as different objects.** (Claude
```

## L7. The comparator's object paths do not exist yet

```
$ cd /var/projects/toy_physics && for f in research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py; do test -e $f && echo "present $f" || echo "absent $f"; done
absent research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
absent research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
```
