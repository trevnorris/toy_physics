# Phase 0 Corpus Scan Measurement

Date: 2026-08-24
Working directory: `/var/projects/toy_physics`

This record supports the structural counts in `memory/_meta/CORPUS.md`. The
scan is read-only. It distinguishes the Git index from the much larger raw
filesystem.

## Tracked-content and project-size scan

Command:

```bash
git rev-parse HEAD
printf 'tracked research: '; git ls-files research | wc -l
printf 'tracked software: '; git ls-files software | wc -l
printf 'tracked research TeX: '; git ls-files research | rg '\.tex$' | wc -l
printf 'tracked research Markdown: '; git ls-files research | rg '\.md$' | wc -l
printf 'tracked research Python: '; git ls-files research | rg '\.py$' | wc -l
printf 'tracked research Wolfram: '; git ls-files research | rg '\.(wl|wls)$' | wc -l
printf 'tracked software Markdown: '; git ls-files software | rg '\.md$' | wc -l
printf 'tracked software Python: '; git ls-files software | rg '\.py$' | wc -l
printf 'tracked software Wolfram: '; git ls-files software | rg '\.(wl|wls)$' | wc -l
printf 'paper entry points: '; rg -l -F '\documentclass' research --glob '*.tex' | wc -l
for path in research/* software/*; do
  if [ -d "$path" ]; then
    tracked=$(git ls-files "$path" | wc -l)
    printf '%5d %s\n' "$tracked" "$path"
  fi
done | sort -k2,2
printf 'tracked scratch paths: '; git ls-files research software | rg '/_scratch/' | wc -l
printf 'tracked directive paths: '; git ls-files research software | rg '/directives?/' | wc -l
printf 'tracked red-team paths: '; git ls-files research software | rg '/redteam' | wc -l
printf 'tracked PDFs: '; git ls-files research software | rg '\.pdf$' | wc -l
printf 'filesystem files under research excluding cache: '; find research -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' -not -path '*/.pytest_cache/*' | wc -l
printf 'filesystem files under software excluding cache: '; find software -type f -not -path '*/.git/*' -not -path '*/__pycache__/*' -not -path '*/.pytest_cache/*' | wc -l
printf 'filesystem scratch files under software: '; find software -type f -path '*/_scratch/*' | wc -l
```

Literal output:

```text
30e96ee22245d4a7d0e873ad228cf5ce33de76f0
tracked research: 7235
tracked software: 604
tracked research TeX: 392
tracked research Markdown: 2774
tracked research Python: 518
tracked research Wolfram: 423
tracked software Markdown: 230
tracked software Python: 211
tracked software Wolfram: 69
paper entry points: 20
   11 research/1pn_hybrid
   16 research/1pn_optics
   13 research/1pn_orbital_dynamics
    5 research/1pn_spin_and_nbody
   34 research/4d
   10 research/4d_1pn_bridge
    5 research/4d_1pn_full
   10 research/4d_2_5pn
   74 research/4d_2pn
   34 research/4d_3pn
   42 research/4d_4pn
   11 research/4d_em_fields
    4 research/4d_plasma
   10 research/brane_bulk_ontology
    8 research/em_fields
    3 research/pde
  185 research/pde_audit
 5565 research/pde_ledger
  650 research/pde_ledger_v2
  544 research/pde_ledger_v3
  160 software/em_charge_attribute
   36 software/force_visualizer
  407 software/stage1_solver
tracked scratch paths: 20
tracked directive paths: 863
tracked red-team paths: 3752
tracked PDFs: 19
filesystem files under research excluding cache: 12535
filesystem files under software excluding cache: 122521
filesystem scratch files under software: 120979
```

## Paper-structure scan

Command:

```bash
rg -l -F '\documentclass' research --glob '*.tex' | sort |
while IFS= read -r f; do
  lines=$(wc -l < "$f")
  inputs=$(rg -c -F '\input{' "$f" 2>/dev/null || true)
  printf '%6d %3d %s\n' "$lines" "${inputs:-0}" "$f"
done
```

The columns are line count, direct `\input{...}` count, and path.

Literal output:

```text
  3112   0 research/1pn_hybrid/paper/1pn_hybrid.tex
  3895   0 research/1pn_optics/paper/1pn_optics.tex
  1316   0 research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex
  1921   0 research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex
  2039   0 research/4d/paper/4d.tex
  4074   0 research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
  6145   0 research/4d_1pn_full/paper/4d_1pn_full.tex
  6146   0 research/4d_2_5pn/paper/4d_2_5pn.tex
  7684   0 research/4d_2pn/paper/4d_2pn.tex
  8983   0 research/4d_3pn/paper/4d_3pn.tex
  7920   0 research/4d_4pn/paper/4d_4pn.tex
  1832   0 research/4d_em_fields/paper/4d_em_fields.tex
  5048   0 research/4d_plasma/paper/4d_plasma.tex
  2935   0 research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
  1972   0 research/em_fields/paper/em_fields.tex
  4024   0 research/pde/paper/pde.tex
    46  22 research/pde_ledger/paper/pde_ledger.tex
    32   7 research/pde_ledger/paper/pde_ledger_reader.tex
    44  20 research/pde_ledger_v2/paper/pde_ledger.tex
    46  10 research/pde_ledger_v3/paper/pde_ledger_v3.tex
```

## Ignore-policy inspection

Command:

```bash
find . -name .gitignore -type f -not -path './.git/*' | sort
sed -n '1,320p' .gitignore
for f in $(find research software -name .gitignore -type f | sort); do
  printf '\nFILE %s\n' "$f"
  sed -n '1,320p' "$f"
done
```

The authoritative literal policy remains in [`.gitignore`](../../.gitignore)
and `software/em_charge_attribute/.gitignore`; it is not duplicated here.
Notable routing consequences are recorded in `memory/_meta/CORPUS.md`.

Representative check:

```bash
printf 'tracked .out evidence paths: '; git ls-files research software | rg '\.out$' | wc -l
git check-ignore -v \
  software/stage1_solver/_scratch/probe.txt \
  software/stage1_solver/runs/probe.json \
  research/pde_audit/simulation/output/probe.json \
  research/pde_ledger_v3/paper/probe.pdf \
  software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production/probe.yaml
if git check-ignore -q research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out; then
  echo ignored
else
  echo not_ignored
fi
```

Literal output:

```text
tracked .out evidence paths: 86
.gitignore:102:**/_scratch/	software/stage1_solver/_scratch/probe.txt
.gitignore:93:software/**/runs/	software/stage1_solver/runs/probe.json
.gitignore:67:research/pde_audit/simulation/output/	research/pde_audit/simulation/output/probe.json
.gitignore:141:research/pde_ledger_v3/paper/*.pdf	research/pde_ledger_v3/paper/probe.pdf
software/em_charge_attribute/.gitignore:19:reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production/	software/em_charge_attribute/reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production/probe.yaml
not_ignored
```
