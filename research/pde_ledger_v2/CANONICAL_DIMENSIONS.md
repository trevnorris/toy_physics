<!--
OUTPUT ONLY: this table observes independent per-stage artifacts.
No script may ever import this file; that would make cross-checks tautological.
Never hand-edit this file; always regenerate it with
`python3 scripts/generate_canonical_dimension_table.py`.
-->

# Canonical dimension table

> **Coverage is incomplete:** this is the generated view of converted stages,
> not a complete corpus view.

- Dimension-bearing stage corpus: **30 stages**.
- Converted and represented: **5 of 30** — stage004, stage011, stage012, stage013, stage018.
- Not yet represented: stage002, stage003, stage005, stage006, stage007, stage008, stage009, stage010, stage016, stage021, stage023, stage027, stage030, stage031, stage032, stage034, stage035, stage036, stage037, stage038, stage039, stage040, stage041, stage042, stage044.
- Total quantity rows: **72** (one per exact emitted name and stage).
- Candidate-same-quantity groups: **2**.
- `NEEDS_ADJUDICATION` groups: **0**.

Values come only from committed `scripts/*.dimensions.txt` and
`mathematica/out/*.out` artifacts. Exact emitted names are primary; candidate
groups below are review flags, never automatic merges. Axis tuples are shown in
each stage's declared order, while canonical renderings are always labelled and
ordered `M`, `L`, `T`. Per-stage status uses the parser, labelled-axis comparison,
and artifact-name waiver registry in `scripts/compare_dimension_artifacts.py`.
A Python sidecar is rejected unless its source-digest assertions match
independent hashes of the current stage source and `ledger_dimensions.py`;
this is freshness, not source coverage.

## Stage coverage

| Stage | Axis order | Python quantities | Wolfram quantities | Quantity rows |
|---|---:|---:|---:|---:|
| stage004 | `(L,T,M)` | 20 | 20 | 20 |
| stage011 | `(L,M,T)` | 12 | 10 | 12 |
| stage012 | `(L,M,T)` | 19 | 18 | 19 |
| stage013 | `(L,M,T)` | 15 | 15 | 15 |
| stage018 | `(L,M,T)` | 6 | 6 | 6 |

## Quantities

| Scope | Quantity (scope-stripped emitted name) | Candidate key | Stage | Axis order | Exponents | Canonical labelled dimension | Per-stage engine status |
|---|---|---|---|---|---|---|---|
| `(none)` | `G3` | `G3` | stage004 | `(L,T,M)` | `{3, -2, -1}` | M⁻¹ L³ T⁻² | AGREE |
| `(none)` | `G4` | `G4` | stage004 | `(L,T,M)` | `{4, -2, -1}` | M⁻¹ L⁴ T⁻² | AGREE |
| `(none)` | `L` | `L` | stage004 | `(L,T,M)` | `{1, 0, 0}` | L | AGREE |
| `(none)` | `M` | `M` | stage004 | `(L,T,M)` | `{0, 0, 1}` | M | AGREE |
| `(none)` | `T` | `T` | stage004 | `(L,T,M)` | `{0, 1, 0}` | T | AGREE |
| `(none)` | `Tw` | `Tw` | stage004 | `(L,T,M)` | `{1, -2, 1}` | M L T⁻² | AGREE |
| `(none)` | `USigmaRR` | `USigmaRR` | stage004 | `(L,T,M)` | `{-1, -2, 1}` | M L⁻¹ T⁻² | AGREE |
| `(none)` | `action` | `action` | stage004 | `(L,T,M)` | `{2, -1, 1}` | M L² T⁻¹ | AGREE |
| `(none)` | `electricField` | `electricField` | stage004 | `(L,T,M)` | `{1, -2, 1}` | M L T⁻² | AGREE |
| `(none)` | `energy` | `energy` | stage004 | `(L,T,M)` | `{2, -2, 1}` | M L² T⁻² | AGREE |
| `(none)` | `force` | `force` | stage004 | `(L,T,M)` | `{1, -2, 1}` | M L T⁻² | AGREE |
| `(none)` | `lagrangianDensity` | `lagrangianDensity` | stage004 | `(L,T,M)` | `{-2, -2, 1}` | M L⁻² T⁻² | AGREE |
| `(none)` | `magneticField` | `magneticField` | stage004 | `(L,T,M)` | `{0, -1, 1}` | M T⁻¹ | AGREE |
| `(none)` | `maxwellCoeff` | `maxwellCoeff` | stage004 | `(L,T,M)` | `{-4, 2, -1}` | M⁻¹ L⁻⁴ T² | AGREE |
| `(none)` | `muWall` | `muWall` | stage004 | `(L,T,M)` | `{-1, 0, 1}` | M L⁻¹ | AGREE |
| `(none)` | `qA0` | `qA0` | stage004 | `(L,T,M)` | `{2, -2, 1}` | M L² T⁻² | AGREE |
| `(none)` | `qAi` | `qAi` | stage004 | `(L,T,M)` | `{1, -1, 1}` | M L T⁻¹ | AGREE |
| `(none)` | `rho3` | `rho3` | stage004 | `(L,T,M)` | `{-3, 0, 0}` | L⁻³ | AGREE |
| `(none)` | `velocity` | `velocity` | stage004 | `(L,T,M)` | `{1, -1, 0}` | L T⁻¹ | AGREE |
| `(none)` | `zero` | `zero` | stage004 | `(L,T,M)` | `{0, 0, 0}` | 1 | AGREE |
| `(none)` | `CorruptCsSquaredDim` | `CorruptCsSquaredDim` | stage011 | `(L,M,T)` | `{3, 0, -2}` | L³ T⁻² | AGREE |
| `(none)` | `CorruptKDim` | `CorruptKDim` | stage011 | `(L,M,T)` | `{19, 1, -2}` | M L¹⁹ T⁻² | AGREE |
| `(none)` | `CsSquaredDim` | `CsSquaredDim` | stage011 | `(L,M,T)` | `{2, 0, -2}` | L² T⁻² | AGREE |
| `(none)` | `EnergyDim` | `EnergyDim` | stage011 | `(L,M,T)` | `{2, 1, -2}` | M L² T⁻² | AGREE |
| `(none)` | `ExpectedCsSquaredDim` | `ExpectedCsSquaredDim` | stage011 | `(L,M,T)` | `{2, 0, -2}` | L² T⁻² | AGREE |
| `(none)` | `FourVolumeDim` | `FourVolumeDim` | stage011 | `(L,M,T)` | `{4, 0, 0}` | L⁴ | AGREE |
| `(none)` | `KDim` | `KDim` | stage011 | `(L,M,T)` | `{18, 1, -2}` | M L¹⁸ T⁻² | AGREE |
| `(none)` | `LengthDim` | `LengthDim` | stage011 | `(L,M,T)` | `{1, 0, 0}` | L | AGREE |
| `(none)` | `MassDim` | `MassDim` | stage011 | `(L,M,T)` | `{0, 1, 0}` | M | ONE_SIDED_PY (WAIVED) |
| `(none)` | `OmegaDim` | `OmegaDim` | stage011 | `(L,M,T)` | `{0, 0, -1}` | T⁻¹ | ONE_SIDED_PY (WAIVED) |
| `(none)` | `PressureDim` | `PressureDim` | stage011 | `(L,M,T)` | `{-2, 1, -2}` | M L⁻² T⁻² | AGREE |
| `(none)` | `RhoDim` | `RhoDim` | stage011 | `(L,M,T)` | `{-4, 0, 0}` | L⁻⁴ | AGREE |
| `(none)` | `K_dim` | `KDim` | stage012 | `(L,M,T)` | `{18, 1, -2}` | M L¹⁸ T⁻² | AGREE |
| `(none)` | `alpha_dim` | `alphaDim` | stage012 | `(L,M,T)` | `{-1, 0, 0}` | L⁻¹ | AGREE |
| `clean_walk` | `cs_dim` | `csDim` | stage012 | `(L,M,T)` | `{1, 0, -1}` | L T⁻¹ | AGREE |
| `clean_walk` | `cs_squared_dim` | `csSquaredDim` | stage012 | `(L,M,T)` | `{2, 0, -2}` | L² T⁻² | AGREE |
| `clean_walk` | `k_dim` | `kDim` | stage012 | `(L,M,T)` | `{-1, 0, 0}` | L⁻¹ | AGREE |
| `clean_walk` | `tan_argument_dim` | `tanArgumentDim` | stage012 | `(L,M,T)` | `{0, 0, 0}` | 1 | AGREE |
| `clean_walk` | `z00_dim` | `z00Dim` | stage012 | `(L,M,T)` | `{-1, 0, 0}` | L⁻¹ | AGREE |
| `clean_walk` | `z00_prefactor_dim` | `z00PrefactorDim` | stage012 | `(L,M,T)` | `{-1, 0, 0}` | L⁻¹ | AGREE |
| `(none)` | `corrupt_K_dim` | `corruptKDim` | stage012 | `(L,M,T)` | `{19, 1, -2}` | M L¹⁹ T⁻² | AGREE |
| `corrupt_walk` | `cs_dim` | `csDim` | stage012 | `(L,M,T)` | `{3/2, 0, -1}` | L³⁄² T⁻¹ | AGREE |
| `corrupt_walk` | `k_dim` | `kDim` | stage012 | `(L,M,T)` | `{-3/2, 0, 0}` | L⁻³⁄² | AGREE |
| `corrupt_walk` | `tan_argument_dim` | `tanArgumentDim` | stage012 | `(L,M,T)` | `{-1/2, 0, 0}` | L⁻¹⁄² | AGREE |
| `corrupt_walk` | `z00_prefactor_dim` | `z00PrefactorDim` | stage012 | `(L,M,T)` | `{-3/2, 0, 0}` | L⁻³⁄² | AGREE |
| `(none)` | `energy_dim` | `energyDim` | stage012 | `(L,M,T)` | `{2, 1, -2}` | M L² T⁻² | AGREE |
| `(none)` | `four_volume_dim` | `fourVolumeDim` | stage012 | `(L,M,T)` | `{4, 0, 0}` | L⁴ | AGREE |
| `(none)` | `mass_dim` | `massDim` | stage012 | `(L,M,T)` | `{0, 1, 0}` | M | ONE_SIDED_PY (WAIVED) |
| `(none)` | `omega_dim` | `omegaDim` | stage012 | `(L,M,T)` | `{0, 0, -1}` | T⁻¹ | AGREE |
| `(none)` | `pressure_dim` | `pressureDim` | stage012 | `(L,M,T)` | `{-2, 1, -2}` | M L⁻² T⁻² | AGREE |
| `(none)` | `rho_dim` | `rhoDim` | stage012 | `(L,M,T)` | `{-4, 0, 0}` | L⁻⁴ | AGREE |
| `(none)` | `K_eta` | `KEta` | stage013 | `(L,M,T)` | `{-1, 1, -2}` | M L⁻¹ T⁻² | AGREE |
| `k_dims` | `LL` | `LL` | stage013 | `(L,M,T)` | `{0, 1, -2}` | M T⁻² | AGREE |
| `k_dims` | `aL` | `aL` | stage013 | `(L,M,T)` | `{0, 1, -2}` | M T⁻² | AGREE |
| `k_dims` | `aa` | `aa` | stage013 | `(L,M,T)` | `{0, 1, -2}` | M T⁻² | AGREE |
| `(none)` | `k_shared` | `kShared` | stage013 | `(L,M,T)` | `{0, 1, -2}` | M T⁻² | AGREE |
| `m_dims` | `LL` | `LL` | stage013 | `(L,M,T)` | `{0, 1, 0}` | M | AGREE |
| `m_dims` | `aL` | `aL` | stage013 | `(L,M,T)` | `{0, 1, 0}` | M | AGREE |
| `m_dims` | `aa` | `aa` | stage013 | `(L,M,T)` | `{0, 1, 0}` | M | AGREE |
| `(none)` | `m_shared` | `mShared` | stage013 | `(L,M,T)` | `{0, 1, 0}` | M | AGREE |
| `(none)` | `ratio_dim` | `ratioDim` | stage013 | `(L,M,T)` | `{0, 0, -2}` | T⁻² | AGREE |
| `symbol_dims` | `L0` | `L0` | stage013 | `(L,M,T)` | `{1, 0, 0}` | L | AGREE |
| `symbol_dims` | `Tw` | `Tw` | stage013 | `(L,M,T)` | `{1, 1, -2}` | M L T⁻² | AGREE |
| `symbol_dims` | `beta` | `beta` | stage013 | `(L,M,T)` | `{-1, 0, 0}` | L⁻¹ | AGREE |
| `symbol_dims` | `muEta` | `muEta` | stage013 | `(L,M,T)` | `{-1, 1, 0}` | M L⁻¹ | AGREE |
| `symbol_dims` | `rAL` | `rAL` | stage013 | `(L,M,T)` | `{0, 0, 0}` | 1 | AGREE |
| `(none)` | `a` | `a` | stage018 | `(L,M,T)` | `{1, 0, 0}` | L | AGREE |
| `(none)` | `c_s0_dim` | `cS0Dim` | stage018 | `(L,M,T)` | `{1, 0, -1}` | L T⁻¹ | AGREE |
| `(none)` | `corrupted_u2_dim` | `corruptedU2Dim` | stage018 | `(L,M,T)` | `{-1, 0, 3}` | L⁻¹ T³ | AGREE |
| `(none)` | `u2` | `u2` | stage018 | `(L,M,T)` | `{0, 0, 2}` | T² | AGREE |
| `(none)` | `u4` | `u4` | stage018 | `(L,M,T)` | `{0, 0, 4}` | T⁴ | AGREE |
| `(none)` | `v5` | `v5` | stage018 | `(L,M,T)` | `{0, 0, 5}` | T⁵ | AGREE |

## Candidate-same-quantity groups

A dotted prefix is parsed as scope and is not part of the quantity name.
Candidate keys normalize separator boundaries in that scope-stripped name
to CamelCase while preserving initial and all other case. Case is meaningful:
`K_dim` is the EOS constant K while scope
`clean_walk`, quantity `k_dim` is the wavenumber k; likewise, documented
`c_E` (L T⁻¹) and `C_E` (M⁻¹ L⁻⁴ T²) are distinct quantities.
A one-sided group is not counted as an agreement even when its visible
dimensions are the same. A differing group is flagged `NEEDS_ADJUDICATION`;
the generator never chooses a winner.

| Case-sensitive candidate key | Members | Status |
|---|---|---|
| `KDim` | stage011 scope `(none)`, `KDim` [M L¹⁸ T⁻²; AGREE]<br>stage012 scope `(none)`, `K_dim` [M L¹⁸ T⁻²; AGREE] | AGREE |
| `Tw` | stage004 scope `(none)`, `Tw` [M L T⁻²; AGREE]<br>stage013 scope `symbol_dims`, `Tw` [M L T⁻²; AGREE] | AGREE |

## GROUPING LIMITATIONS

The general rule actually finds `KDim` (stage011) and `K_dim`
(stage012). It does not group the following known-same cross-stage
quantities because their initial case differs after separator-boundary
normalization:

- `CsSquaredDim` (stage011) and scope `clean_walk`, quantity
  `cs_squared_dim` (stage012);
- `CorruptKDim` (stage011) and `corrupt_K_dim` (stage012);
- `EnergyDim` (stage011) and `energy_dim` (stage012);
- `FourVolumeDim` (stage011) and `four_volume_dim` (stage012);
- `MassDim` (stage011) and `mass_dim` (stage012);
- `OmegaDim` (stage011) and `omega_dim` (stage012);
- `PressureDim` (stage011) and `pressure_dim` (stage012);
- `RhoDim` (stage011) and `rho_dim` (stage012).

These misses are reported only; they do not feed candidate grouping.
