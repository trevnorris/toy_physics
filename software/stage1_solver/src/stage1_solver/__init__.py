"""Stage-1 branch-realization solver skeleton.

The package is intentionally small at step 1: reusable finite-volume radial
operators, cubic-GPE residuals, Newton/JVP machinery, independent references,
and a validation harness.
"""

from .config import (
    BackendConfig,
    CubicGPEConfig,
    BranchSmokeConfig,
    CoupledBranchMMSConfig,
    HarnessConfig,
    LinearEigenConfig,
    ManufacturedSolutionsConfig,
    MaxwellMMSConfig,
    NewtonConfig,
    CurrentMMSConfig,
    QuinticMatterMMSConfig,
    RadialGridSpec,
    TensorLaplacianMMSConfig,
    TensorGridSpec,
    WallMMSConfig,
    WallGridSpec,
)

__all__ = [
    "BackendConfig",
    "BranchSmokeConfig",
    "CoupledBranchMMSConfig",
    "CubicGPEConfig",
    "HarnessConfig",
    "LinearEigenConfig",
    "ManufacturedSolutionsConfig",
    "MaxwellMMSConfig",
    "NewtonConfig",
    "CurrentMMSConfig",
    "QuinticMatterMMSConfig",
    "RadialGridSpec",
    "TensorLaplacianMMSConfig",
    "TensorGridSpec",
    "WallMMSConfig",
    "WallGridSpec",
]
