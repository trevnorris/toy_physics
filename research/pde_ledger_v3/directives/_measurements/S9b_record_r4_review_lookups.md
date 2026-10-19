# Measurements: S9b record scoped review, round 4 (generated 2026-10-19 07:39)

Generator: `_scratch/s9b_build/gen/s9b_record_r4_lookups.sh` (sha256sum/grep/sed/cut/wc/git only). The legs'
evidence is hard-linked at `_scratch/s9b_build/record/s9b_record_review_r4_{claude,grok}_evidence/` (the Claude
leg's isolated clone is left in its scratch directory).

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_baseline_r4.sha256 | grep -c ': OK$')
4
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_prompt_r4.md | cut -c1-64
ecd7e2e579b14e15ab42afe8f9c6b88caaf01a05a7d7488546c8f236a63e026c
```

```
$ git -C /var/projects/toy_physics diff --stat d877e857 -- research/pde_ledger_v3/steps research/pde_ledger_v3/scripts/S9b_record_measurements.py research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
 research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md   |  30 +-
 .../scripts/S9b_record_measurements.py             |  47 +-
 research/pde_ledger_v3/steps/S9b_light_bending.md  |  38 +-
 .../steps/_measurements/S9b_record_measurements.md | 527 ++++++++++++++++++++-
 4 files changed, 595 insertions(+), 47 deletions(-)
```

## Verdicts
```
$ grep -o 'Verdict:[^A-Z]*[A-Z][A-Z]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_claude.md | head -1
Verdict: CLEAR
```

```
$ grep -o 'Verdict:[^A-Z]*[A-Z][A-Z]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_grok.txt | head -1
Verdict: CLEAR
```

## Both legs: the generator reproduces the filed measurements
```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md | cut -c1-64
158fd122647cae25f9e793ae2e17abb7f874a7fadcab45c7949bbe54bcde01ba
```

```
$ grep -h 'BYTE-IDENTICAL' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_claude_evidence/regen_compare.stdout.txt | head -1
BYTE-IDENTICAL
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_grok_evidence/regen/S9b_record_measurements.md | cut -c1-64
158fd122647cae25f9e793ae2e17abb7f874a7fadcab45c7949bbe54bcde01ba
```

## Item 1: the user decision's written source
```
$ git -C /var/projects/toy_physics show 4760e6e2:research/pde_ledger_v3/directives/_measurements/S9b_record_r2_review_disposition.md | sed -n '19,21p'
- **A user decision this round.** Asked whether a check of the flow-only condition against the river/acoustic-flow
  prior art should be queued, the user chose "Queue it after record". The check is new computation, done after the
  record is accepted and outside it: Codex writes the check and it is reviewed under G1.
```

```
$ sed -n '771,778p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
**User decision, 2026-10-19:** the written source is
`directives/_measurements/S9b_record_r2_review_disposition.md` at `4760e6e2`, lines 19–21,
retrieved verbatim in M21. The user chose **"Queue it after record"** for the check of the
flow-only condition against the river/acoustic-flow prior art. The check is new computation,
done **after the record is accepted and outside it: Codex writes the check and it is reviewed
under G1**. This queues the check without selecting an oracle's applicability, importing a flow
profile, or answering its second-route question. The identified river/inflow/acoustic candidates
in the table remain open with that handoff; no oracle is applied here.
```

```
$ sed -n '1076,1083p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
**User decision, 2026-10-19:** the written source is
`directives/_measurements/S9b_record_r2_review_disposition.md` at `4760e6e2`, lines 19–21,
retrieved verbatim in M21. The user chose **"Queue it after record"** for the check of the
flow-only condition against the river/acoustic-flow prior art. The check is new computation,
done **after the record is accepted and outside it: Codex writes the check and it is reviewed
under G1**. This queues the check without selecting an oracle's applicability, importing a flow
profile, or answering its second-route question. The identified river/inflow/acoustic candidates
in the table remain open with that handoff; no oracle is applied here.
```

```
$ grep -h 'IDENTICAL' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_claude_evidence/m21_check.stdout.txt | head -1
M21 == source lines 19-21: IDENTICAL
```

```
$ grep -h 'diff exit' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_claude_evidence/record_vs_register_oracle.stdout.txt | head -1
diff exit 0
```

## Item 2: the rows that carry the unresolved brane-zero value samples
```
$ grep -h 'mine == M22\|named - mine\|all my rows covered' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_claude_evidence/check_cells_vs_list.stdout.txt
mine == M22 row map: True | rows 41 | samples 123
   named - mine: []
   named - mine: []
   named - mine: []
all my rows covered by some cell: True | uncovered: []
```

```
$ head -12 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_grok_evidence/item2/cell_vs_retrieval.txt
A cell rows 31
subtraction rows 35
A equal True
A only cell []
A only retrieval []
B equal True
B cell ['F_FLOW_B_RESIDUAL_DEFLECTION', 'F_FLOW_B_RESIDUAL_RADAR', 'F_LIVE_B_RESIDUAL_DEFLECTION', 'F_LIVE_B_RESIDUAL_RADAR']
B retrieval ['F_FLOW_B_RESIDUAL_DEFLECTION', 'F_FLOW_B_RESIDUAL_RADAR', 'F_LIVE_B_RESIDUAL_DEFLECTION', 'F_LIVE_B_RESIDUAL_RADAR']
BR equal True
BR cell ['BRANCH_EXISTENCE', 'F_FLOW_BRANCH_EXISTENCE', 'F_FLOW_PATH_TRAVERSAL', 'F_LIVE_BRANCH_EXISTENCE', 'F_LIVE_PATH_TRAVERSAL', 'PATH_TRAVERSAL']
BR retrieval ['BRANCH_EXISTENCE', 'F_FLOW_BRANCH_EXISTENCE', 'F_FLOW_PATH_TRAVERSAL', 'F_LIVE_BRANCH_EXISTENCE', 'F_LIVE_PATH_TRAVERSAL', 'PATH_TRAVERSAL']
---M21---
```

## The Claude leg's observation: the A-class cell's wording
```
$ grep -n 'not assigned to every' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
561:| A_DEFLECTION, A_ROUND_TRIP, A_ONE_WAY_ER/RE, A_NONRECIPROCAL, A_RADAR_LOG; conditional F_FLOW/F_LIVE copies | Formed value samples and markers; all twelve applicable grades retained, four radar grades | No full along-ray domain comparison; no all-profile equality; printed numeric differences remain. The unresolved `brane_zero_at_point` value samples in this row class belong to `F_FLOW_A_DEFL
```

```
$ grep -h 'rows with brane-zero value samples' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r4_claude_evidence/general_restricted_A.stdout.txt
general A_* rows with brane-zero value samples: 44 | rows with any unresolved: 0 | outcome totals: {'residual formed': 129, 'boolean comparison formed': 3, 'unresolved coverage': 0}
restricted R_*_A_* rows with brane-zero value samples: 30 | rows with any unresolved: 0 | outcome totals: {'residual formed': 90, 'boolean comparison formed': 0, 'unresolved coverage': 0}
```

