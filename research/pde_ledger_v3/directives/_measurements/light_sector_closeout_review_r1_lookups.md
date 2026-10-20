# Measurements: light-closeout review round 1 (generated 2026-10-19 19:15)

Generator: `_scratch/light_closeout/gen/review_r1_lookups.sh`. Literal retrieval only. Packet lines are read
from the reviewed revision, commit `4c5b4c31` plus the r1 working-tree repair preserved in the next commit.

## E1. The round: prompt, baseline and verdict lines
```
$ sha256sum /var/projects/toy_physics/_scratch/light_closeout/review_prompt_r1.md | cut -c1-64
38d2211a702c23947f11546644c961bc2bbb0b9a65c730107f92c0c2c5aea76e
```

```
$ cat /var/projects/toy_physics/_scratch/light_closeout/packet_baseline_r1.sha256 | sed 's#/var/projects/toy_physics/##'
4df5da71e12f855cd4724faf57c0172332caaefebdc1959798f3ea196d577544  research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md
c008fa1785302357d3a7ed3f87feb40584a67fecac414e645aa2e838ba42846a  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_lookups.md
1ae4b63efffd9f52c43a00a653a311059178e860ef3ed53b656d1950e320eec0  research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_outcome_lookups.md
b42dc759d6973ddfb562bc5d79a7796b06e8ff7815846c0ce9d436741a5dd968  research/pde_ledger_v3/directives/_measurements/light_sector_reassessment_lookups.md
313bc3d3452d6ff931eb7e2fdcf915b2a8068350d817b7cc088d4d11e0d3bb18  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r0_disposition.md
0067a5e204ea57bbe4f4a572442ae53375e2c3f92b4e411857e1f4b55596065f  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r0_lookups.md
6b8b7ee0ee812953178682471e507d1f22deda7ddda2bb490bfb56b0dac5ad42  STATUS.md
a44d0b9c3334659288318bc0e240ae3de7956245e02c2acf93f02cc6fb9664ad  research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
0a25e5b39e20a379804590b2dd48e12f7bbf0600f98ab2d36f9568e2a74a6f9c  research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
```

```
$ grep -n '^FINDINGS\|^## FINDINGS\|^1\. ' /var/projects/_review_scratch/light-closeout-r1-codex/review_report.txt /var/projects/_review_scratch/light-closeout-r1-grok/grok_review.txt | cut -c1-200
/var/projects/_review_scratch/light-closeout-r1-codex/review_report.txt:396:FINDINGS
/var/projects/_review_scratch/light-closeout-r1-codex/review_report.txt:398:1. Medium — machine-capacity conclusion remains stronger than the measured resource failures.
/var/projects/_review_scratch/light-closeout-r1-grok/grok_review.txt:127:## FINDINGS
/var/projects/_review_scratch/light-closeout-r1-grok/grok_review.txt:129:1. Medium. `R-S22-01`'s source assigns E38 the radiating-mode condition that the survey states for the blackbody rows.
```

## E2. Grok's finding: the E38 condition
```
$ grep -n -A1 'acts as an additional' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
768:    under the condition that the extra longitudinal branch couples to matter and acts as an additional
769-    radiating mode. (C6, C7.)
```

```
$ sed -n '280p' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | grep -o 'E32–E34 on matter coupling.\{220\}'
E32–E34 on matter coupling and, for the blackbody comparison, ordinary thermal-radiation participation; E38 additionally on the extra branch's cosmological relativistic population, energy/perturbation mapping and decoupling history. E22's equal-s
```

```
$ sed -n '137p' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | grep -o 'An extra branch constrains this count.\{200\}'
An extra branch constrains this count only through its cosmological population, dispersion, gravitating energy/perturbations and temperature/decoupling history; a mode census alone predicts no \(\Delta N_{\rm eff}\). | [Primary Planck pa
```

## E3. Codex's finding: the opening of DEFERRED_HEAVY_RUNS
```
$ grep -n -A1 'OOM.d the earlier' /var/projects/toy_physics/research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
6:Some in-band engine self-checks OOM'd the earlier **30 GB** dev box, and several still do not fit the 64 GB box
7-in use since 2026-10-18 (see *Parked* below). The un-reduced variable-coefficient
```

```
$ grep -n 'memory needed to finish is unknown' /var/projects/toy_physics/research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
37:the OOM kills at 54 GiB, and the memory needed to finish is unknown. The guard's pool budget was 62 GiB, and
```

