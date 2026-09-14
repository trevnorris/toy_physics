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

Remaining native producers: d.

The c2 artifact checkpoint is 537d78fd. Its post-save hash check stopped the
queue on truncated annex content before d started. The complete producer output
restored the existing annex key; full size/SHA256 and git-annex fsck now pass
for all five b/c1/c2 transcripts. See
[recovery evidence](S11c_thickness_coordinate_c2_annex_recovery_checkpoint.json).
The truncation cause remains unresolved. Resume at d with the existing checks.

Fresh endpoint/reference sources, two-frequency pairing and full current/adjoint normalization remain separate next steps. Native point/path/stratum records do not establish global coverage, scattering or section 3b profile-frequency bound poles.
