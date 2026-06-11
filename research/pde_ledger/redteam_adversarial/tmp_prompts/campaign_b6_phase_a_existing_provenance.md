# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_161.md
- redteam/pass2/reports/stage_162.md
- redteam/pass2/reports/stage_163.md
- redteam/pass2/reports/stage_164.md
- redteam/pass2/reports/stage_165.md
- redteam/pass2/reports/stage_166.md
- redteam/pass2/reports/stage_167.md
- redteam/pass2/reports/stage_168.md
- redteam/pass2/reports/stage_169.md
- redteam/pass2/reports/stage_170.md
- redteam/pass2/reports/stage_171.md
- redteam/pass2/reports/stage_172.md
- redteam/pass2/reports/stage_173.md
- redteam/pass2/reports/stage_174.md
- redteam/pass2/reports/stage_175.md
- redteam/pass2/reports/stage_176.md
- redteam/pass2/reports/stage_177.md
- redteam/pass2/reports/stage_178.md
- redteam/pass2/reports/stage_179.md
- redteam/pass2/reports/stage_180.md
- redteam/pass2/reports/stage_181.md
- redteam/pass2/reports/stage_182.md
- redteam/pass2/reports/stage_183.md
- redteam/pass2/reports/stage_184.md
- redteam/pass2/reports/stage_185.md
- redteam/pass2/reports/stage_186.md
- redteam/pass2/reports/stage_187.md
- redteam/pass2/reports/stage_188.md
- redteam/pass2/reports/stage_189.md
- redteam/pass2/reports/stage_190.md
- redteam/pass2/reports/stage_191.md
- redteam/pass2/reports/stage_192.md`
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
