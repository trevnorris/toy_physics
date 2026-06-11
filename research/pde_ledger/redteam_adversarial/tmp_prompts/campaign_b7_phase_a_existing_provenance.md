# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_193.md
- redteam/pass2/reports/stage_194.md
- redteam/pass2/reports/stage_195.md
- redteam/pass2/reports/stage_196.md
- redteam/pass2/reports/stage_197.md
- redteam/pass2/reports/stage_198.md
- redteam/pass2/reports/stage_199.md
- redteam/pass2/reports/stage_200.md
- redteam/pass2/reports/stage_201.md
- redteam/pass2/reports/stage_202.md
- redteam/pass2/reports/stage_203.md
- redteam/pass2/reports/stage_204.md
- redteam/pass2/reports/stage_205.md
- redteam/pass2/reports/stage_206.md
- redteam/pass2/reports/stage_207.md
- redteam/pass2/reports/stage_208.md
- redteam/pass2/reports/stage_209.md
- redteam/pass2/reports/stage_210.md
- redteam/pass2/reports/stage_211.md
- redteam/pass2/reports/stage_212.md
- redteam/pass2/reports/stage_213.md
- redteam/pass2/reports/stage_214.md
- redteam/pass2/reports/stage_215.md
- redteam/pass2/reports/stage_216.md
- redteam/pass2/reports/stage_217.md
- redteam/pass2/reports/stage_218.md
- redteam/pass2/reports/stage_219.md
- redteam/pass2/reports/stage_220.md
- redteam/pass2/reports/stage_221.md
- redteam/pass2/reports/stage_222.md
- redteam/pass2/reports/stage_223.md
- redteam/pass2/reports/stage_224.md`
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
