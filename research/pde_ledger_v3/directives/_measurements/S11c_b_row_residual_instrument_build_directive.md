# S11c-b requested-truncation row-residual instrument — builder report

## Deliverable and run

The deliverable is `research/pde_ledger_v3/scripts/S11c_b_row_residual.py`.

Command:

```text
cd /var/projects/toy_physics && wc -l research/pde_ledger_v3/scripts/S11c_b_row_residual.py
```

Output:

```text
1064 research/pde_ledger_v3/scripts/S11c_b_row_residual.py
```

The transcript-only run command returned exit 0:

```text
cd /var/projects/toy_physics && python3 research/pde_ledger_v3/scripts/S11c_b_row_residual.py > /tmp/S11c_b_row_residual.out 2> /tmp/S11c_b_row_residual.err
```

Command:

```text
wc -c -l /tmp/S11c_b_row_residual.out /tmp/S11c_b_row_residual.err
```

Output:

```text
     547 16570943 /tmp/S11c_b_row_residual.out
       0        0 /tmp/S11c_b_row_residual.err
     547 16570943 total
```

Thus stdout is 16,570,943 bytes and stderr is empty.

## Requested truncation

The production function is `requested_truncation` at
`scripts/S11c_b_row_residual.py:171`; its scalar derivative-at-origin route is
at lines 156-168. The independent fixture implementation is
`independent_series_projection` at line 177 and uses sequential direct
`sp.series(..., 0, 2)` calls at lines 180-182.

Command:

```text
nl -ba scripts/S11c_b_row_residual.py | sed -n '155,182p'
```

The independent cross-check and one-sided corruption were emitted before the
guard. Command:

```text
rg '^(REQUESTED_TRUNCATION_CROSSCHECK|ONE_SIDED_CORRUPTION_(BEFORE|AFTER)) = ' /tmp/S11c_b_row_residual.out
```

Literal stdout:

```text
REQUESTED_TRUNCATION_CROSSCHECK = true
ONE_SIDED_CORRUPTION_BEFORE = Add(Mul(Pow(Symbol('eta_bg', real=True), Integer(2)), Symbol('fixture_eta_dropped')), Mul(Symbol('eta_bg', real=True), Symbol('fixture_kept'), Symbol('sigma_W', real=True)), Mul(Symbol('fixture_sigma_dropped'), Pow(Symbol('sigma_W', real=True), Integer(2))))
ONE_SIDED_CORRUPTION_AFTER = Mul(Symbol('eta_bg', real=True), Symbol('fixture_kept'), Symbol('sigma_W', real=True))
```

The before/after shows that the `eta_bg**2` and `sigma_W**2` terms are removed
while the `eta_bg*sigma_W` term survives.

## Equivalence dispatch

The named dispatch is `equivalence_dispatch` at
`scripts/S11c_b_row_residual.py:384`:

- strong exact branch: lines 393-399;
- weak modulo-exact-total-divergence branch, gated through
  `A._is_weak_scalar_density` and `A.classify_total_divergence`: lines 400-411;
- ordered componentwise branch: lines 412-413.

Command:

```text
nl -ba scripts/S11c_b_row_residual.py | sed -n '384,414p'
```

The run emitted 24 exact strong leaves, 20 weak coupling objects, and 52
ordered admissibility components. Command:

```text
rg '^EMITTED_ROW_ACCOUNTING = ' /tmp/S11c_b_row_residual.out
```

Literal stdout:

```text
EMITTED_ROW_ACCOUNTING = Association({'STRONG_EXACT': Integer(24), 'WEAK_MODULO_EXACT_TOTAL_DIVERGENCE': Integer(20), 'COMPONENTWISE_ORDERED_PAIRING': Integer(52)})
```

## Assembly accounting

The mass and thickness/kinetic accounting guards were emitted after all row
measurements. Command:

```text
rg '^ASSEMBLY_ACCOUNTING_GUARD_(MASS_ROW|THICKNESS_KINETIC_ROW) = ' /tmp/S11c_b_row_residual.out
```

Literal stdout:

```text
ASSEMBLY_ACCOUNTING_GUARD_MASS_ROW = Association({'EXPECTED_LAYER_OPERANDS': Integer(12), 'ACCOUNTED_LAYER_OPERANDS': Integer(12), 'COUNT_DIFFERENCE': Integer(0)})
ASSEMBLY_ACCOUNTING_GUARD_THICKNESS_KINETIC_ROW = Association({'EXPECTED_LAYER_OPERANDS': Integer(20), 'ACCOUNTED_LAYER_OPERANDS': Integer(20), 'COUNT_DIFFERENCE': Integer(0)})
```

The complete extraction-leaf multiset also closed. Command:

```text
rg '^ASSEMBLY_ACCOUNTING_GUARD = ' /tmp/S11c_b_row_residual.out
```

Literal stdout:

```text
ASSEMBLY_ACCOUNTING_GUARD = Association({'EXPECTED_LAYER_OPERANDS': Integer(248), 'ACCOUNTED_LAYER_OPERANDS': Integer(248), 'COUNT_DIFFERENCE': Integer(0), 'MULTISET_DIFFERENCE': Association({'MISSING': (), 'EXTRA': ()})})
```

The WL `FACE`/`FLUX` leaves are emitted as `ROW_FACE_ATTRIBUTED` and have the
accounting role `FACE_ATTRIBUTED_AND_EXCLUDED_FROM_RESIDUAL`; kinetic-origin
leaves are passed through `A._kinetic_pairs` only as
`PROVENANCE_ONLY_ALREADY_IN_COMPLETE_ROW`, never additively assembled.

## Both coupling blocks and relabelled adjoints

Command:

```text
rg '^COUPLING_BLOCK_ACCOUNTING = ' /tmp/S11c_b_row_residual.out
```

Literal stdout:

```text
COUPLING_BLOCK_ACCOUNTING = Association({'ADJOINTNESS_OPERAND_FORWARD:TRANSVERSE_TO_THICKNESS': Integer(4), 'ADJOINTNESS_OPERAND_REVERSE:THICKNESS_TO_TRANSVERSE': Integer(4), 'ADJOINTNESS_RESIDUAL:NO_SECTOR': Integer(4), 'BLOCK:THICKNESS_TO_TRANSVERSE': Integer(4), 'BLOCK:TRANSVERSE_TO_THICKNESS': Integer(4)})
```

All 20 weak objects have a certificate emission. Command:

```text
rg -c '^ROW_DIVERGENCE_CERTIFICATE = ' /tmp/S11c_b_row_residual.out
```

Output:

```text
20
```

## Residual emission and guard ordering

All 96 residuals use the single unconditional emission site at
`scripts/S11c_b_row_residual.py:443`; the arithmetic and assembly guards are
after the measurement loops at lines 1032-1052. There is no `assert` in the
file, so no assertion precedes any residual emission.

Command:

```text
rg -n 'assert|"ROW_RESIDUAL"|Arithmetic and assembly guards' scripts/S11c_b_row_residual.py
```

Output:

```text
443:        "ROW_RESIDUAL",
1032:        # Arithmetic and assembly guards intentionally follow every measured
```

Command:

```text
rg -c '^ROW_RESIDUAL = ' /tmp/S11c_b_row_residual.out
```

Output:

```text
96
```

Final static checks:

```text
python3 -m py_compile scripts/S11c_b_row_residual.py
git diff --check -- scripts/S11c_b_row_residual.py directives/_measurements/S11c_b_row_residual_instrument_build_directive.md
```

Both commands returned exit 0 with no output.
