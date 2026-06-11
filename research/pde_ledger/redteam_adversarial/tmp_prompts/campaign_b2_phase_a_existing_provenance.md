# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `033, 034, 035, 036, 037, 038, 039, 040, 041, 042, 043, 044, 045, 046, 047, 048, 049, 050, 051, 052, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_033.md
- redteam/pass2/reports/stage_034.md
- redteam/pass2/reports/stage_035.md
- redteam/pass2/reports/stage_036.md
- redteam/pass2/reports/stage_037.md
- redteam/pass2/reports/stage_038.md
- redteam/pass2/reports/stage_039.md
- redteam/pass2/reports/stage_040.md
- redteam/pass2/reports/stage_041.md
- redteam/pass2/reports/stage_042.md
- redteam/pass2/reports/stage_043.md
- redteam/pass2/reports/stage_044.md
- redteam/pass2/reports/stage_045.md
- redteam/pass2/reports/stage_046.md
- redteam/pass2/reports/stage_047.md
- redteam/pass2/reports/stage_048.md
- redteam/pass2/reports/stage_049.md
- redteam/pass2/reports/stage_050.md
- redteam/pass2/reports/stage_051.md
- redteam/pass2/reports/stage_052.md
- redteam/pass2/reports/stage_053.md
- redteam/pass2/reports/stage_054.md
- redteam/pass2/reports/stage_055.md
- redteam/pass2/reports/stage_056.md
- redteam/pass2/reports/stage_057.md
- redteam/pass2/reports/stage_058.md
- redteam/pass2/reports/stage_059.md
- redteam/pass2/reports/stage_060.md
- redteam/pass2/reports/stage_061.md
- redteam/pass2/reports/stage_062.md
- redteam/pass2/reports/stage_063.md
- redteam/pass2/reports/stage_064.md`
- Checkpoint provenance seed: `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`

Task:
Every pass-2 value-reconciliation entry and checkpoint-provenance value that maps to these stages must become either a fit-insertion-point candidate or an explicit pure-identity classification. Do not exempt a value because it was previously audited.

Emit only YAML:

```yaml
modality: existing_provenance
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    citation:
      path:
      line:
      excerpt:
    reason:
pure_identities:
  - stage:
    value:
    citation:
      path:
      line:
      excerpt:
    reason:
```
