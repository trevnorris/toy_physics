# S11c thickness-coordinate endpoint checkpoints

Fresh calculations against the regenerated native producer. Each stage uses its own validator; numerical residuals retain their recorded domains.

## right_source

```json
{
  "tagCount": 3484,
  "checkedSourceObjects": 85,
  "checkedMetadataPaths": 368,
  "literalSourceResidualScalars": 33,
  "retainedNonzeroScalars": 0,
  "cancellationIdentityCount": 127,
  "nonzeroCancellationIdentities": 0,
  "wallSeconds": 453.55444050673395,
  "peakRssKiB": 201604
}
```

## left_source

```json
{
  "tagCount": 2314,
  "checkedSourceObjects": 85,
  "checkedMetadataPaths": 368,
  "literalSourceResidualScalars": 33,
  "retainedNonzeroScalars": 0,
  "cancellationIdentityCount": 82,
  "nonzeroCancellationIdentities": 0,
  "wallSeconds": 126.23547302465886,
  "peakRssKiB": 201724
}
```

## reference_source

```json
{
  "tagCount": 146,
  "residualScalarCount": 30,
  "nonzeroResidualCount": 0,
  "wallSeconds": 40.5872370749712,
  "peakRssKiB": 200784
}
```

## right_frequency

```json
{
  "tagCount": 1018,
  "objectCount": 508,
  "metadataPaths": 6573,
  "residualScalars": 2886
}
```

## left_frequency

```json
{
  "tagCount": 1018,
  "objectCount": 508,
  "metadataPaths": 6573,
  "residualScalars": 2886
}
```

## right_pairing

```json
{
  "tagCount": 374,
  "objectCount": 187,
  "metadataPaths": 4332,
  "retainedResidualScalars": 802,
  "retainedNonzeroScalars": 0,
  "wallSeconds": 962.0330720329657,
  "peakRssKiB": 202116
}
```

Remaining planned stages: left_pairing, reference_pairing.

Full endpoint current/adjoint maps and outward orientations precede variable-profile matching. No global exceptional coverage, complete scattering or profile-frequency bound-pole claim.
