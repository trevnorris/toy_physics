# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_129.md
- redteam/pass2/reports/stage_130.md
- redteam/pass2/reports/stage_131.md
- redteam/pass2/reports/stage_132.md
- redteam/pass2/reports/stage_133.md
- redteam/pass2/reports/stage_134.md
- redteam/pass2/reports/stage_135.md
- redteam/pass2/reports/stage_136.md
- redteam/pass2/reports/stage_137.md
- redteam/pass2/reports/stage_138.md
- redteam/pass2/reports/stage_139.md
- redteam/pass2/reports/stage_140.md
- redteam/pass2/reports/stage_141.md
- redteam/pass2/reports/stage_142.md
- redteam/pass2/reports/stage_143.md
- redteam/pass2/reports/stage_144.md
- redteam/pass2/reports/stage_145.md
- redteam/pass2/reports/stage_146.md
- redteam/pass2/reports/stage_147.md
- redteam/pass2/reports/stage_148.md
- redteam/pass2/reports/stage_149.md
- redteam/pass2/reports/stage_150.md
- redteam/pass2/reports/stage_151.md
- redteam/pass2/reports/stage_152.md
- redteam/pass2/reports/stage_153.md
- redteam/pass2/reports/stage_154.md
- redteam/pass2/reports/stage_155.md
- redteam/pass2/reports/stage_156.md
- redteam/pass2/reports/stage_157.md
- redteam/pass2/reports/stage_158.md
- redteam/pass2/reports/stage_159.md
- redteam/pass2/reports/stage_160.md`
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
