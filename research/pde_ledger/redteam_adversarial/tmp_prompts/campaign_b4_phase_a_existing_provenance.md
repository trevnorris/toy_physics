# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_097.md
- redteam/pass2/reports/stage_098.md
- redteam/pass2/reports/stage_099.md
- redteam/pass2/reports/stage_100.md
- redteam/pass2/reports/stage_101.md
- redteam/pass2/reports/stage_102.md
- redteam/pass2/reports/stage_103.md
- redteam/pass2/reports/stage_104.md
- redteam/pass2/reports/stage_105.md
- redteam/pass2/reports/stage_106.md
- redteam/pass2/reports/stage_107.md
- redteam/pass2/reports/stage_108.md
- redteam/pass2/reports/stage_109.md
- redteam/pass2/reports/stage_110.md
- redteam/pass2/reports/stage_111.md
- redteam/pass2/reports/stage_112.md
- redteam/pass2/reports/stage_113.md
- redteam/pass2/reports/stage_114.md
- redteam/pass2/reports/stage_115.md
- redteam/pass2/reports/stage_116.md
- redteam/pass2/reports/stage_117.md
- redteam/pass2/reports/stage_118.md
- redteam/pass2/reports/stage_119.md
- redteam/pass2/reports/stage_120.md
- redteam/pass2/reports/stage_121.md
- redteam/pass2/reports/stage_122.md
- redteam/pass2/reports/stage_123.md
- redteam/pass2/reports/stage_124.md
- redteam/pass2/reports/stage_125.md
- redteam/pass2/reports/stage_126.md
- redteam/pass2/reports/stage_127.md
- redteam/pass2/reports/stage_128.md`
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
