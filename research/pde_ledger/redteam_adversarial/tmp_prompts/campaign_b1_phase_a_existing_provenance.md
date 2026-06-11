# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `001, 002, 003, 004, 005, 006, 007, 008, 009, 010, 011, 012, 013, 014, 015, 016, 017, 018, 019, 020, 021, 022, 023, 024, 025, 026, 027, 028, 029, 030, 031, 032`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_001.md
- redteam/pass2/reports/stage_002.md
- redteam/pass2/reports/stage_003.md
- redteam/pass2/reports/stage_004.md
- redteam/pass2/reports/stage_005.md
- redteam/pass2/reports/stage_006.md
- redteam/pass2/reports/stage_007.md
- redteam/pass2/reports/stage_008.md
- redteam/pass2/reports/stage_009.md
- redteam/pass2/reports/stage_010.md
- redteam/pass2/reports/stage_011.md
- redteam/pass2/reports/stage_012.md
- redteam/pass2/reports/stage_013.md
- redteam/pass2/reports/stage_014.md
- redteam/pass2/reports/stage_015.md
- redteam/pass2/reports/stage_016.md
- redteam/pass2/reports/stage_017.md
- redteam/pass2/reports/stage_018.md
- redteam/pass2/reports/stage_019.md
- redteam/pass2/reports/stage_020.md
- redteam/pass2/reports/stage_021.md
- redteam/pass2/reports/stage_022.md
- redteam/pass2/reports/stage_023.md
- redteam/pass2/reports/stage_024.md
- redteam/pass2/reports/stage_025.md
- redteam/pass2/reports/stage_026.md
- redteam/pass2/reports/stage_027.md
- redteam/pass2/reports/stage_028.md
- redteam/pass2/reports/stage_029.md
- redteam/pass2/reports/stage_030.md
- redteam/pass2/reports/stage_031.md
- redteam/pass2/reports/stage_032.md`
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
