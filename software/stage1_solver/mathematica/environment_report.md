# Mathematica 15 Environment Check

**Generated:** 2026-06-17 11:14:30 GMT-06:00

**Working invocation:** `timeout 600 wolframscript -script software/stage1_solver/mathematica/mt15_00_environment_check.wls`

**Kernel:** `$VersionNumber = 15.`; `15.0.0 for Linux x86 (64-bit) (May 19, 2026)`

**Machine/license:** `$MachineName = pop-os`; `$ProcessorCount = 8`; `$MaxLicenseProcesses = 2`

## Compact Status

| Item | Route area | Status |
| --- | --- | --- |
| 2 | FEM / NDSolve | PRESENT |
| 3 | Coordinate / curvilinear support | PRESENT |
| 4 | Eigensolvers | PRESENT |
| 5 | Sparse LinearSolve methods | PRESENT (caveat: small sparse smoke only; unknown-method sanity probe is excluded from this aggregate) |
| 6 | Matrix decompositions | PRESENT (caveat: some named symbols exist but did not pass the tiny-matrix decomposition smoke) |
| 7 | GPU | ABSENT |
| 8 | Export | PRESENT |
| 9 | ModelFit | PRESENT |

## Hoped-For Capabilities Requiring Re-Plan

- `PopovDecomposition on a tiny matrix`: ABSENT - defined = yes; tiny matrix smoke = no; Global shadow = no; result head = PopovDecomposition; result = Short[PopovDecomposition[{{4., 1.}, {2., 3.}}], 3]
- `PartialFractions on a tiny matrix`: ABSENT - defined = yes; tiny matrix smoke = no; Global shadow = no; result head = Symbol; result = Short[$Failed, 3]
- `PolynomialReduction on a tiny matrix`: ABSENT - defined = yes; tiny matrix smoke = no; Global shadow = no; result head = Symbol; result = Short[$Failed, 3]
- `LinearSolve[A, b, Method -> "HybridCPUGPU"]`: ABSENT - result = Short[$Failed, 3]

## Surprising / Useful Present Capabilities

- `LinearSolve[A, b, Method -> "Pardiso"]`: PRESENT - relative error = 0.
- `LinearSolve[A, b, Method -> "MUMPS"]`: PRESENT (caveat: small sparse smoke only; no production-scale or linked-library benchmark) - relative error = 1.1737147695332612*^-16
- `LDLDecomposition on a tiny matrix`: PRESENT - defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List
- `BunchKaufmanDecomposition on a tiny matrix`: PRESENT - defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List
- `OrderedSchurDecomposition on a tiny matrix`: PRESENT - defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List
- `PolarDecomposition on a tiny matrix`: PRESENT - defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List

## Detailed Probe Results

| ID | Route need | Probe | Status | Details |
| --- | --- | --- | --- | --- |
| 1.1 | gating prerequisite | `$VersionNumber >= 15.0` | PRESENT | $VersionNumber = 15. |
| 1.2 | gating prerequisite | `$Version / $MachineName / $ProcessorCount` | PRESENT | $Version = 15.0.0 for Linux x86 (64-bit) (May 19, 2026); $MachineName = pop-os; $ProcessorCount = 8 |
| 1.3 | gating prerequisite | `$MaxLicenseProcesses == 2` | PRESENT | $MaxLicenseProcesses = 2 |
| 2.1 | stationary open-throat + linearized BVPs | `Needs["NDSolve`FEM`"]` | PRESENT | load result = Short[True, 3] |
| 2.2 | stationary open-throat + linearized BVPs | `1-D FEM BVP u''[x] == -1 with Dirichlet endpoints` | PRESENT | Head = InterpolatingFunction; u(0.5) = Short[0.12500000000000053, 3] |
| 2.3 | stationary open-throat + linearized BVPs | `2-D unit-square Poisson FEM smoke` | PRESENT | Head = InterpolatingFunction; u(0.5,0.5) = Short[0.07367112522324072, 3] |
| 3.1 | keep throat coordinates natural | `CoordinateChartData["Spherical", ...]` | PRESENT | System symbol exists = True; sample = Short[#1[[1]] > 0 && 0 < #1[[2]] < Pi && Inequality[-Pi, Less, #1[[3]], LessEqual, Pi] & , 3] |
| 3.2 | keep throat coordinates natural | `Laplacian/Grad/Div with "Spherical" coordinate spec` | PRESENT | Laplacian ok = True; Grad ok = True; Div ok = True |
| 4.1 | BdG / mixed modal extraction | `Eigensystem on a small dense matrix` | PRESENT | eigenvalues = Short[{3.618033988749895, 1.381966011250105}, 3] |
| 4.2 | BdG / mixed modal extraction | `Generalized eigenproblem Eigenvalues[{A, B}]` | PRESENT | eigenvalues = Short[{2., 1.4999999999999998}, 3] |
| 4.3 | BdG / mixed modal extraction | `Sparse Arnoldi Eigenvalues[sparseA, 3, Method -> {"Arnoldi"}]` | PRESENT | eigenvalues = Short[{3.9776616524502586, 3.911145611572275, 3.801937735804832}, 3] |
| 5.Multifrontal | scalable linear algebra bottleneck | `LinearSolve[A, b, Method -> "Multifrontal"]` | PRESENT | relative error = 1.131782766116681*^-16 |
| 5.Pardiso | scalable linear algebra bottleneck | `LinearSolve[A, b, Method -> "Pardiso"]` | PRESENT | relative error = 0. |
| 5.Krylov | scalable linear algebra bottleneck | `LinearSolve[A, b, Method -> "Krylov"]` | PRESENT | relative error = 1.1737147695332612*^-16 |
| 5.Banded | scalable linear algebra bottleneck | `LinearSolve[A, b, Method -> "Banded"]` | PRESENT | relative error = 4.663869570635836*^-17 |
| 5.MUMPS | scalable linear algebra bottleneck | `LinearSolve[A, b, Method -> "MUMPS"]` | PRESENT (caveat: small sparse smoke only; no production-scale or linked-library benchmark) | relative error = 1.1737147695332612*^-16 |
| 5.DefinitelyNotAMethod | scalable linear algebra bottleneck | `LinearSolve[A, b, Method -> "DefinitelyNotAMethod"]` | ABSENT | result = Short[$Failed, 3] |
| 6.LUDecomposition | stability classification + pole extraction | `LUDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.CholeskyDecomposition | stability classification + pole extraction | `CholeskyDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.SchurDecomposition | stability classification + pole extraction | `SchurDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.JordanDecomposition | stability classification + pole extraction | `JordanDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.QRDecomposition | stability classification + pole extraction | `QRDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.SingularValueDecomposition | stability classification + pole extraction | `SingularValueDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.LDLDecomposition | stability classification + pole extraction | `LDLDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.BunchKaufmanDecomposition | stability classification + pole extraction | `BunchKaufmanDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.OrderedSchurDecomposition | stability classification + pole extraction | `OrderedSchurDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.PolarDecomposition | stability classification + pole extraction | `PolarDecomposition on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.PopovDecomposition | stability classification + pole extraction | `PopovDecomposition on a tiny matrix` | ABSENT | defined = yes; tiny matrix smoke = no; Global shadow = no; result head = PopovDecomposition; result = Short[PopovDecomposition[{{4., 1.}, {2., 3.}}], 3] |
| 6.JordanReduce | stability classification + pole extraction | `JordanReduce on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.FrobeniusReduce | stability classification + pole extraction | `FrobeniusReduce on a tiny matrix` | PRESENT | defined = yes; tiny matrix smoke = yes; Global shadow = no; result head = List |
| 6.PartialFractions | stability classification + pole extraction | `PartialFractions on a tiny matrix` | ABSENT | defined = yes; tiny matrix smoke = no; Global shadow = no; result head = Symbol; result = Short[$Failed, 3] |
| 6.PolynomialReduction | stability classification + pole extraction | `PolynomialReduction on a tiny matrix` | ABSENT | defined = yes; tiny matrix smoke = no; Global shadow = no; result head = Symbol; result = Short[$Failed, 3] |
| 7.1 | nice-to-have GPU acceleration | `Needs["CUDALink`"]` | PRESENT | load result = Short[True, 3] |
| 7.2 | nice-to-have GPU acceleration | `GPUArray symbol and small construction smoke` | ABSENT | System symbol exists = True; construction result = Short[$Failed, 3] |
| 7.3 | nice-to-have GPU acceleration | `LinearSolve[A, b, TargetDevice -> "GPU"]` | ABSENT | result = Short[$Failed, 3] |
| 7.4 | nice-to-have GPU acceleration | `LinearSolve[A, b, Method -> "HybridCPUGPU"]` | ABSENT | result = Short[$Failed, 3] |
| 8.1 | V2-22B packet handoff | `Export/Import round trip through .json` | PRESENT | Export result = Short["/var/projects/toy_physics/software/stage1_solver/mathematica/mt15_00_json_probe.json", 3]; Import result = Short[<\|"alpha" -> 1, "beta" -> {True, "x"}\|>, 3] |
| 8.2 | human/LLM-readable adjacent exports only | `YAML listed in $ExportFormats and Export[..., "YAML"]` | PRESENT | $ExportFormats contains YAML = True; Export result = Short["/var/projects/toy_physics/software/stage1_solver/mathematica/mt15_00_yaml_probe.yaml", 3] |
| 9.1 | post-freeze diagnostics only | `LinearModelFit smoke` | PRESENT | Head = FittedModel; fitted value at x=3 = Short[7., 3] |
| 9.2 | post-freeze diagnostics only | `NonlinearModelFit smoke` | PRESENT | Head = FittedModel; fitted value at x=1 = Short[2.7182818284590455, 3] |

## Notes

- This is an environment capability report, not a performance benchmark.
- Linear algebra `PRESENT` rows mean the method was accepted and returned a correct answer on the stated small smoke problem.
- GPU rows are nice-to-have only; failures there do not block the CPU route.
- The `DefinitelyNotAMethod` LinearSolve row is an internal sanity check that unknown methods are reported as absent.
