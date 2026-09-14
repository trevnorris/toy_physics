# S11c thickness-coordinate regeneration status

Serial native producers and recorded physical/artifact checks. Each completed producer is committed locally; transcripts use DataLad/git-annex. No push.

## Completed b

```json
{
  "stage": "b",
  "exportRows": 2441,
  "changedValueSerializations": [
    "slab_operator",
    "slab_operator_term_origins"
  ],
  "addedKeys": [],
  "removedKeys": [],
  "nativeOutputBytes": 183361030,
  "nativeOutputSha256": "693f0e804ca5fc3e418406cd049aa3d4c56263e101a7f1ae3078172d674cf98f",
  "exportSha256": "3c6555f90a85382a5e96f3f961f0e36a4875d9259fe1be947fa8e683669b7084",
  "failures": [],
  "residualScalars": 2634,
  "nonzeroResidualScalars": 0
}
```

## Completed c1

```json
{
  "stage": "c1",
  "exportRows": 44,
  "changedValueSerializations": [],
  "addedKeys": [],
  "removedKeys": [],
  "nativeOutputBytes": 90722854,
  "nativeOutputSha256": "eb4bd18893804a76b15218628fd994834c6bce63b9655f2c42f0f14d855d6da4",
  "exportSha256": "a5353b86dc526fd1b56734fdbdf7dc7f42387a3195ab4de84e8994f751871cdc",
  "failures": []
}
```

## Completed c2

```json
{
  "stage": "c2",
  "exportRows": 70,
  "changedValueSerializations": [
    "s11cc2ClosedSlabOperator"
  ],
  "addedKeys": [],
  "removedKeys": [],
  "nativeOutputBytes": 530883300,
  "nativeOutputSha256": "9712191e3af5e7bbc2eb65824ca283f0c5cb3a35d5e21632d0821412ec018864",
  "exportSha256": "2ba5ed487b84ecb98290df6f1650fa6210fe2784cfa02a0f3cff96cd718411cd",
  "failures": [],
  "residualScalars": 2846,
  "nonzeroResidualScalars": 0
}
```

## Completed d

```json
{
  "nativeWallSeconds": 4166.187866047025,
  "nativePeakRssKiB": 2476192,
  "nativeArtifact": {
    "bytes": 83989687,
    "sha256": "19822837896b084fb8bcbd11b16aa5120e2c5b01c1c3084357aed279b7b70ce8"
  },
  "spectrum": {
    "packets": 24,
    "candidates": 432,
    "nullities": {
      "1": 336,
      "2": 96
    },
    "coverage": {
      "(11, 9, 11, 0, True, True)": 24
    },
    "maximumResidualByDimension": {
      "PULLBACK_RESIDUAL|(-5, -2, 1)": 0.0,
      "PULLBACK_RESIDUAL|(-3, -2, 1)": 0.0,
      "PULLBACK_RESIDUAL|(-5, -1, 1)": 0.0,
      "PULLBACK_RESIDUAL|(-3, -1, 1)": 0.0,
      "PULLBACK_RESIDUAL|(-1, -2, 1)": 0.0,
      "PULLBACK_RESIDUAL|(0, 0, 0)": 0.0,
      "RIGHT_RESIDUAL|(-2, -2, 1)": 4.7231752266122105e-15,
      "RIGHT_RESIDUAL|(-3, -1, 1)": 1.5118541179209328e-15,
      "RIGHT_RESIDUAL|(-1, -2, 1)": 2.6947576617039376e-15,
      "LEFT_RESIDUAL|(-1, 0, 0)": 3.773782390822243e-15,
      "LEFT_RESIDUAL|(0, 0, 0)": 3.2023728339893768e-15,
      "PROJECTOR_RESIDUAL|(0, 0, 0)": 4.449299500566018e-14,
      "PROJECTOR_RESIDUAL|(1, 0, 0)": 1.452356741593167e-14,
      "PROJECTOR_RESIDUAL|(-1, 0, 0)": 6.485071087474307e-15
    }
  },
  "jointSheet": {
    "packets": 24,
    "pathStatuses": {
      "TRANSPORTED": 480,
      "BRANCH_LOCUS_ON_PATH": 24
    }
  },
  "outstandingConstructions": [
    "FULL_END_SPECTRA_BEYOND_REFERENCE_MODE_JETS",
    "GENERIC_DOMAIN_SHEET_CONTINUATION",
    "CLOSED_NONLOCAL_BULK_CURRENT_AND_FLUX_NORMALIZATION",
    "COMPLETE_TWO_ENDED_SCATTERING",
    "POLES_RIESZ_OVERLAP",
    "SURVIVAL",
    "FLUX_BOOKKEEPING",
    "WEAK_COEFFICIENTS",
    "SECTION_5_CONTROLS",
    "OWN_ROWS_EXPORT"
  ],
  "failures": []
}
```

Remaining native producers: none in this regeneration queue.

Fresh endpoint/reference sources, two-frequency pairing and full current/adjoint normalization remain separate next steps. Native point/path/stratum records do not establish global coverage, scattering or section 3b profile-frequency bound poles.
