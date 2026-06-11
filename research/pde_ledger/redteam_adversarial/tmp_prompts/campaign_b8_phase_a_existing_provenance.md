# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_225.md
- redteam/pass2/reports/stage_226.md
- redteam/pass2/reports/stage_227.md
- redteam/pass2/reports/stage_228.md
- redteam/pass2/reports/stage_229.md
- redteam/pass2/reports/stage_230.md
- redteam/pass2/reports/stage_231.md
- redteam/pass2/reports/stage_232.md
- redteam/pass2/reports/stage_233.md
- redteam/pass2/reports/stage_234.md
- redteam/pass2/reports/stage_235.md
- redteam/pass2/reports/stage_236.md
- redteam/pass2/reports/stage_237.md
- redteam/pass2/reports/stage_238.md
- redteam/pass2/reports/stage_239.md
- redteam/pass2/reports/stage_240.md
- redteam/pass2/reports/stage_241.md
- redteam/pass2/reports/stage_242.md
- redteam/pass2/reports/stage_243.md
- redteam/pass2/reports/stage_244.md
- redteam/pass2/reports/stage_245.md
- redteam/pass2/reports/stage_246.md
- redteam/pass2/reports/stage_247.md
- redteam/pass2/reports/stage_248.md
- redteam/pass2/reports/stage_249.md
- redteam/pass2/reports/stage_250.md
- redteam/pass2/reports/stage_251.md
- redteam/pass2/reports/stage_252.md
- redteam/pass2/reports/stage_253.md`
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
