# S11c-b row-residual coupling fix — implementation measurement

Artifact: `scripts/S11c_b_row_residual.py`
Committed comparison layer: `scripts/S11c_b_adjudicated_comparison.py` (unchanged)
Preserved pre-fix stdout: `/tmp/S11c_b_row_residual.out`
Fixed stdout: `/tmp/S11c_b_row_residual.fixed1.out`, `/tmp/S11c_b_row_residual.fixed2.out`

This is an instrument-build record only. It makes no residual-value, route,
in-scope, or which-engine judgment.

## Implementation evidence

### (a) Full-order classification uses the pre-bridge A-B residual

The coupling path forms `residual_pre_bridge_d` directly from the two retained
pre-`_bridge_d` operands, then calls the reviewed layer classifier with
`apply_bridge_d=True`. The result is emitted separately as
`FULL_PREBRIDGE_ROUTE` and `EULER_SIGNATURE`.

Command:

```bash
nl -ba scripts/S11c_b_row_residual.py | sed -n '849,860p;492,503p'
```

Relevant output:

```text
849  residual_pre_bridge_d = A._arithmetic_residual(
850      aligned.py_pre_bridge_d,
851      aligned.wl_pre_bridge_d,
852  )
856  full = A.classify_total_divergence(
857      residual_pre_bridge_d,
858      A.PRODUCTION_FIELD_REGISTRY,
859      apply_bridge_d=True,
860  )
493  "FULL_PREBRIDGE_ROUTE",
497  "EULER_SIGNATURE",
501  "EULER_SIGNATURE_SHA256",
```

### (b) `ROW_RESIDUAL` is the bridged, requested-truncated A-B density

The two emitted coupling operands and their A-B residual are all truncated
only after `_bridge_d`. The pre-bridge operand is not truncated.

Command:

```bash
nl -ba scripts/S11c_b_row_residual.py | sed -n '862,875p;504,515p'
```

Relevant output:

```text
862  py_trunc = requested_truncation(
863      A._bridge_d(aligned.py_pre_bridge_d)
865  wl_trunc = requested_truncation(
866      A._bridge_d(aligned.wl_pre_bridge_d)
868  residual = requested_truncation(
869      A._bridge_d(residual_pre_bridge_d)
```

The emitted labels at lines 504-515 are `ROW_OPERAND_PY_TRUNC`,
`ROW_OPERAND_WL`, and `ROW_RESIDUAL`, respectively.

### (c) The instrument computes `T(bridge(R_pre - div V))`, with fallback

The explicit homotopy call, formal divergence, normalization, bridge, and
requested truncation are inside the instrument. The `except Exception` path
marks the remainder unavailable; emission always retains `ROW_RESIDUAL` and
emits `NO_CLEAN_QUOTIENT`.

Command:

```bash
nl -ba scripts/S11c_b_row_residual.py | sed -n '877,900p;512,530p'
```

Relevant output:

```text
879  try:
880      vector = A._homotopy_vector(
881          residual_pre_bridge_d,
882          A.PRODUCTION_FIELD_REGISTRY,
884      remainder_pre_bridge_d = A._normalise_exact(
885          residual_pre_bridge_d
886          - A.formal_divergence(
891      in_scope_weak_remainder = requested_truncation(
892          A._bridge_d(remainder_pre_bridge_d)
898  except Exception as error:
899      in_scope_weak_remainder = None
900      no_clean_quotient_error = type(error).__name__
513  "ROW_RESIDUAL",
517  _emit("NO_CLEAN_QUOTIENT", ...)
520  "IN_SCOPE_WEAK_REMAINDER",
```

No Euler or homotopy operator is applied to a bridged residual. A static check
also finds no old `apply_bridge_d=False` call:

```bash
rg -n 'apply_bridge_d=False|classify_total_divergence|_homotopy_vector|formal_divergence' scripts/S11c_b_row_residual.py
```

```text
856:        full = A.classify_total_divergence(
880:            vector = A._homotopy_vector(
886:                - A.formal_divergence(
```

### (d) Witness canonicalization changes presentation only

The emission copy sorts association entries, applies pinned `signsimp`, and
uses lexicographically ordered expression text. SHA-256 is computed from that
presentation copy. The copy is made only after route, row residual, and
homotopy remainder computation; it is never assigned to the residual,
remainder, or `vector`.

Command:

```bash
nl -ba scripts/S11c_b_row_residual.py | sed -n '381,415p;877,904p'
```

Relevant output:

```text
381  def _deterministic_witness_presentation(value: object) -> object:
391      for label, item in sorted(
405      canonical_copy = sp.signsimp(value)
408      sp.sstr(canonical_copy, order="lex"),
413  def _witness_sha256(value: object) -> str:
414      payload = C.serialise(value).encode("utf-8")
415      return hashlib.sha256(payload).hexdigest()
877  in_scope_weak_remainder: sp.Expr | None
879  try:
902  euler_signature_presentation = (
903      _deterministic_witness_presentation(full.euler_signature)
```

## Full runs, coverage, and guards

Commands (independent processes; the second pins a hash seed):

```bash
/usr/bin/time -v python scripts/S11c_b_row_residual.py \
  > /tmp/S11c_b_row_residual.fixed1.out \
  2> /tmp/S11c_b_row_residual.fixed1.err
PYTHONHASHSEED=1 /usr/bin/time -v python scripts/S11c_b_row_residual.py \
  > /tmp/S11c_b_row_residual.fixed2.out \
  2> /tmp/S11c_b_row_residual.fixed2.err
rg -n 'Elapsed|Maximum resident|Exit status' \
  /tmp/S11c_b_row_residual.fixed1.err \
  /tmp/S11c_b_row_residual.fixed2.err
```

```text
fixed1 Elapsed: 2:02:28; Maximum resident: 838900 kbytes; Exit status: 0
fixed2 Elapsed: 2:03:14; Maximum resident: 834272 kbytes; Exit status: 0
```

Emission-count command:

```bash
for f in /tmp/S11c_b_row_residual.fixed{1,2}.out; do
  rg -c '^FULL_PREBRIDGE_ROUTE = ' "$f"
  rg -c '^EULER_SIGNATURE = ' "$f"
  rg -c '^EULER_SIGNATURE_SHA256 = ' "$f"
  rg -c '^IN_SCOPE_WEAK_REMAINDER = ' "$f"
  rg -c '^NO_CLEAN_QUOTIENT = false$' "$f"
  rg --include-zero -c '^NO_CLEAN_QUOTIENT = true$' "$f"
done
```

Both runs returned, in order: `20, 20, 20, 20, 20, 0`. Each run also emitted
24 `STRONG_EXACT`, 20 `WEAK_MODULO_EXACT_TOTAL_DIVERGENCE`, and 52
`COMPONENTWISE_ORDERED_PAIRING` cases.

Final guard command:

```bash
tail -n 5 /tmp/S11c_b_row_residual.fixed1.out
```

```text
COUPLING_BLOCK_ACCOUNTING = Association({'ADJOINTNESS_OPERAND_FORWARD:TRANSVERSE_TO_THICKNESS': Integer(4), 'ADJOINTNESS_OPERAND_REVERSE:THICKNESS_TO_TRANSVERSE': Integer(4), 'ADJOINTNESS_RESIDUAL:NO_SECTOR': Integer(4), 'BLOCK:THICKNESS_TO_TRANSVERSE': Integer(4), 'BLOCK:TRANSVERSE_TO_THICKNESS': Integer(4)})
ASSEMBLY_ACCOUNTING_GUARD_MASS_ROW = Association({'EXPECTED_LAYER_OPERANDS': Integer(12), 'ACCOUNTED_LAYER_OPERANDS': Integer(12), 'COUNT_DIFFERENCE': Integer(0)})
ASSEMBLY_ACCOUNTING_GUARD_THICKNESS_KINETIC_ROW = Association({'EXPECTED_LAYER_OPERANDS': Integer(20), 'ACCOUNTED_LAYER_OPERANDS': Integer(20), 'COUNT_DIFFERENCE': Integer(0)})
ASSEMBLY_ACCOUNTING_GUARD = Association({'EXPECTED_LAYER_OPERANDS': Integer(248), 'ACCOUNTED_LAYER_OPERANDS': Integer(248), 'COUNT_DIFFERENCE': Integer(0), 'MULTISET_DIFFERENCE': Association({'MISSING': (), 'EXTRA': ()})})
EMITTED_ROW_ACCOUNTING = Association({'STRONG_EXACT': Integer(24), 'WEAK_MODULO_EXACT_TOTAL_DIVERGENCE': Integer(20), 'COMPONENTWISE_ORDERED_PAIRING': Integer(52)})
```

The guards remain after every residual emission:

```bash
nl -ba scripts/S11c_b_row_residual.py | sed -n '1148,1185p'
```

The comment and guards begin at lines 1162-1169; none tests a route, residual,
remainder, in-scope verdict, or engine judgment.

## Two-invocation byte-identical witness demonstration

Command:

```bash
cmp \
  <(rg '^(EULER_SIGNATURE|EULER_SIGNATURE_SHA256) = ' \
      /tmp/S11c_b_row_residual.fixed1.out) \
  <(rg '^(EULER_SIGNATURE|EULER_SIGNATURE_SHA256) = ' \
      /tmp/S11c_b_row_residual.fixed2.out)
echo "all_20_witness_lines_cmp_exit=$?"
sha256sum \
  <(rg '^(EULER_SIGNATURE|EULER_SIGNATURE_SHA256) = ' \
      /tmp/S11c_b_row_residual.fixed1.out) \
  <(rg '^(EULER_SIGNATURE|EULER_SIGNATURE_SHA256) = ' \
      /tmp/S11c_b_row_residual.fixed2.out)
```

```text
all_20_witness_lines_cmp_exit=0
937a841b3fc72d56236952a3404823240d3f7c235b86eda1cfe7008cb6157315  /dev/fd/61
937a841b3fc72d56236952a3404823240d3f7c235b86eda1cfe7008cb6157315  /dev/fd/60
```

The compared streams contain 20 `EULER_SIGNATURE` lines and 20 matching
per-witness hash lines each.

## Literal pre-fix vs fixed non-coupling diff

The diff helper streams a selected complete `CASE` block through its terminal
`ROW_ASSEMBLY_ACCOUNTING` line. `cmp` was run on the resulting byte streams;
no symbolic reparsing or normalization was used.

Command summary:

```bash
extract_equivalence() {
  python -c 'import sys
wanted = sys.argv[2]
emit = False
for line in open(sys.argv[1], encoding="utf-8"):
    if line.startswith("CASE = "):
        emit = ("TextAtom(\047" + wanted + "\047)") in line
    if emit:
        sys.stdout.write(line)
    if emit and line.startswith("ROW_ASSEMBLY_ACCOUNTING = "):
        emit = False' "$1" "$2"
}
extract_non_coupling() {
  python -c 'import sys
skip = False
for line in open(sys.argv[1], encoding="utf-8"):
    if line.startswith("CASE = "):
        skip = "TextAtom(\047WEAK_MODULO_EXACT_TOTAL_DIVERGENCE\047)" in line
    if not skip:
        sys.stdout.write(line)
    if skip and line.startswith("ROW_ASSEMBLY_ACCOUNTING = "):
        skip = False' "$1"
}
for eq in STRONG_EXACT COMPONENTWISE_ORDERED_PAIRING; do
  cmp \
    <(extract_equivalence /tmp/S11c_b_row_residual.out "$eq") \
    <(extract_equivalence /tmp/S11c_b_row_residual.fixed1.out "$eq")
  echo "$eq cmp_exit=$?"
  extract_equivalence /tmp/S11c_b_row_residual.out "$eq" | sha256sum
  extract_equivalence /tmp/S11c_b_row_residual.fixed1.out "$eq" | sha256sum
done
cmp \
  <(extract_non_coupling /tmp/S11c_b_row_residual.out) \
  <(extract_non_coupling /tmp/S11c_b_row_residual.fixed1.out)
echo "ALL_NON_COUPLING cmp_exit=$?"
```

Literal diff summary:

```text
STRONG_EXACT cmp_exit=0 old_cases=24 new_cases=24 old_sha256=3884bd18a7fb3fc04606dea771aa7b60510f022fb18dc928f447b18be0162313 new_sha256=3884bd18a7fb3fc04606dea771aa7b60510f022fb18dc928f447b18be0162313
COMPONENTWISE_ORDERED_PAIRING cmp_exit=0 old_cases=52 new_cases=52 old_sha256=fc06f133ceddbe280b874b6a8040999059f709a0ebccaad5ed6ec801297106c5 new_sha256=fc06f133ceddbe280b874b6a8040999059f709a0ebccaad5ed6ec801297106c5
ALL_NON_COUPLING cmp_exit=0 old_sha256=fcd679d9c4e6794c93c0f91e5b3a8affde329b642ae232224a354e7912109130 new_sha256=fcd679d9c4e6794c93c0f91e5b3a8affde329b642ae232224a354e7912109130
```

The all-non-coupling stream removes only coupling `CASE` blocks and retains
the fixtures, strong rows (including mass), face-only rows, admissibility
rows, and final accounting emissions.

## Layer immutability and syntax

Commands:

```bash
python -m py_compile scripts/S11c_b_row_residual.py
git diff -- scripts/S11c_b_adjudicated_comparison.py
```

The compile command exited 0. The layer diff printed no output.
