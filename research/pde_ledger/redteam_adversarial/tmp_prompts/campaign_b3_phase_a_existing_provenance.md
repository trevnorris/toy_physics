# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 082, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 093, 094, 095, 096`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_065.md
- redteam/pass2/reports/stage_066.md
- redteam/pass2/reports/stage_067.md
- redteam/pass2/reports/stage_068.md
- redteam/pass2/reports/stage_069.md
- redteam/pass2/reports/stage_070.md
- redteam/pass2/reports/stage_071.md
- redteam/pass2/reports/stage_072.md
- redteam/pass2/reports/stage_073.md
- redteam/pass2/reports/stage_074.md
- redteam/pass2/reports/stage_075.md
- redteam/pass2/reports/stage_076.md
- redteam/pass2/reports/stage_077.md
- redteam/pass2/reports/stage_078.md
- redteam/pass2/reports/stage_079.md
- redteam/pass2/reports/stage_080.md
- redteam/pass2/reports/stage_081.md
- redteam/pass2/reports/stage_082.md
- redteam/pass2/reports/stage_083.md
- redteam/pass2/reports/stage_084.md
- redteam/pass2/reports/stage_085.md
- redteam/pass2/reports/stage_086.md
- redteam/pass2/reports/stage_087.md
- redteam/pass2/reports/stage_088.md
- redteam/pass2/reports/stage_089.md
- redteam/pass2/reports/stage_090.md
- redteam/pass2/reports/stage_091.md
- redteam/pass2/reports/stage_092.md
- redteam/pass2/reports/stage_093.md
- redteam/pass2/reports/stage_094.md
- redteam/pass2/reports/stage_095.md
- redteam/pass2/reports/stage_096.md`
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
