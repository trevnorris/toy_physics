# Phase 0 Schema and Benchmark Validation

Date: 2026-08-24
Working directory: `/var/projects/toy_physics`

This read-only check supports only the structural statement that the benchmark's
named source paths are Git-tracked, its TeX labels and representative non-TeX
anchors resolve, and the two new Markdown files pass Git's whitespace check. It
does not validate any scientific conclusion.

## Source paths and whitespace

Command:

```bash
set -e
benchmark=memory/_meta/retrieval-benchmark.md
total=0
missing=0
while IFS= read -r path; do
  total=$((total + 1))
  if ! git ls-files --error-unmatch -- "$path" >/dev/null 2>&1; then
    printf 'MISSING_PATH=%s\n' "$path"
    missing=$((missing + 1))
  fi
done < <(rg -o '`(research|software)/[^`[:space:]]+' "$benchmark" | tr -d '`' | sort -u)
printf 'UNIQUE_SOURCE_PATHS=%s\nMISSING_SOURCE_PATHS=%s\n' "$total" "$missing"
git diff --check -- memory/_meta/schema.md memory/_meta/retrieval-benchmark.md
printf 'DIFF_CHECK=PASS\n'
```

Literal output:

```text
UNIQUE_SOURCE_PATHS=34
MISSING_SOURCE_PATHS=0
DIFF_CHECK=PASS
```

## Anchors

Command:

```bash
set -e
benchmark=memory/_meta/retrieval-benchmark.md
total=0
missing=0
while IFS= read -r anchor; do
  total=$((total + 1))
  if ! rg -q -F -- "$anchor" research software; then
    printf 'MISSING_TEX_ANCHOR=%s\n' "$anchor"
    missing=$((missing + 1))
  fi
done < <(rg -o '\\label\{[^}]+\}' "$benchmark" | sort -u)
printf 'UNIQUE_TEX_LABELS=%s\nMISSING_TEX_LABELS=%s\n' "$total" "$missing"
for target in \
  'research/pde_ledger_v3/steps/S11bA_interface_response.md|## ⭐⭐⭐ The headline — the leak is FREQUENCY-DEPENDENT' \
  'research/pde_ledger_v3/steps/S11bB_interface_assembly.md|## ⭐⭐⭐ THE HEADLINE — the velocity-coupled leak costs an ENERGY RESERVOIR' \
  'research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py|def compare_nullspace_basis(' \
  'research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out|CATEGORY: SUMMARY' \
  'research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl|buildPackage[' \
  'software/stage1_solver/STAGE1_VERDICT.md|## Interpretation — what it means (and does NOT mean)' \
  'software/em_charge_attribute/puncture_deflection_electric_sign_check.py|def section4_adjudicate(' \
  'software/em_charge_attribute/magnetism_moving_throat_check.py|def derive_route_b('; do
  path=${target%%|*}
  anchor=${target#*|}
  rg -q -F -- "$anchor" "$path"
done
printf 'REPRESENTATIVE_NON_TEX_ANCHORS=8\nMISSING_REPRESENTATIVE_NON_TEX_ANCHORS=0\n'
```

Literal output:

```text
UNIQUE_TEX_LABELS=25
MISSING_TEX_LABELS=0
REPRESENTATIVE_NON_TEX_ANCHORS=8
MISSING_REPRESENTATIVE_NON_TEX_ANCHORS=0
```
