# Measurements: S9b record review, round 2 (generated 2026-10-19 05:58)

Generator: `_scratch/s9b_build/gen/s9b_record_r2_lookups.sh` (sha256sum/grep/sed/head/wc only). The legs'
evidence is hard-linked at `_scratch/s9b_build/record/s9b_record_review_r2_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_baseline_r2.sha256 | grep -c ': OK$')
4
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_prompt_r2.md | cut -c1-64
ac4306755fb598aba4720728cee465608f0a2b16c74aaf0a9f74c84c7facbe9b
```

## Verdicts
```
$ head -1 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r2_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r2_grok.txt
Verdict: CLEAR
```

## Both legs: the generator reproduces the filed measurements
```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md | cut -c1-64
ff14384c1250b500f9b09435ad61ba7252bae9da4a9f3adb81a6f99e0e4936c6
```

```
$ grep -h 'WHOLE_FILE_LITERAL_EQUAL' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r2_claude_evidence/compare_regen.stdout | head -1
WHOLE_FILE_LITERAL_EQUAL True
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r2_grok_evidence/S9b_record_measurements.rerun.md | cut -c1-64
ff14384c1250b500f9b09435ad61ba7252bae9da4a9f3adb81a6f99e0e4936c6
```

## Claude 1: the unperformed oracle applicability checks
```
$ sed -n '152p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_record_directive.md
- No new CAS, derivation, oracle check, closure, normalization or profile choice.
```

```
$ sed -n '835p;1099p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
doubt"* — and says *"⇒ ⭐ **Run this route for every sector with identified prior art**, ⛔ not only the
⇒ ⭐ **Run this route for every sector with identified prior art**, ⛔ not only the light sector.
```

```
$ grep -n -i 'oracle\|applicab' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
32:or oracle was rerun. Harness evidence concerns its own engine only. (M1/M8/M12.)
91:traversal; a pointwise inequality is not the full applicability domain. PY declares
108:  nonnegative-root `speed_identity`/`root_domain` enters the observable applicability `gate`
141:applicability domains, the γ gates and the individually named branch tags must not be conflated.
287:Their attached domain keeps the solved ray/applicability gate, `2*GM/r≥0` throughout
510:| A_DEFLECTION, A_ROUND_TRIP, A_ONE_WAY_ER/RE, A_NONRECIPROCAL, A_RADAR_LOG; conditional F_FLOW/F_LIVE copies | Formed value samples and markers; all twelve applicable grades retained, four radar grades | No full along-ray domain comparison; no all-profile equality; unresolved brane-zero samples and printed numeric differences remain. Shared nonreciprocal objects retain §2/M14: K11 supports l
692:not fresh literature verification. **No oracle is applied in this record.**
696:| GR/PPN deflection and logarithmic radar references, A–C spec (M1) | Supplied comparison inputs in the accepted constructions; no new oracle check or GR solution reproduced here. |
701:| MacCullagh (1839) curl-only light / Whittaker Ch. V; Cauchy–Navier elasticity — survey attribution table; register R-S8-04/05 (M10) | None applied as an oracle. Existing light-carrier, longitudinal/compression, angular-momentum and objectivity obligations remain OPEN; the supplied optical ratio does not re-certify them. |
706:This route supplies no new entry: no additional sourced requirement follows from an oracle result
709:M13's bounded condition comparison do not discharge them. Applicability checks of these candidate
710:oracles remain unperformed, not negative results. The register's earlier pass-2 disagreement stays
```

```
$ sed -n '602,662p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md | grep -c -i 'oracle\|applicab'
0
```

```
$ sed -n '736,745p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
Questions requiring work beyond retrieval remain OPEN: whether differing within-WL condition fields
state the same condition (distinct from the item-5 count in M17); the engine-specific branch
bindings and every same-named branch-row difference in M19; complete cross-engine equivalence or
physical interpretation of numeric/domain/condition/γ/exchange differences; the free turning-radius,
constant and uncovered stratum questions; a substrate derivation and physical holder for the local
wave equation; optical/relaxation compatibility; whether any joint live profiles and a physical
source satisfy the conditions relative to orbital GM; Φ–GM and jn–vdr relations; O1–O7 closure and
O2's returned derivative interpretation. Ownership remains as in §7; no computation to answer these
questions is part of this record.

```

```
$ sed -n '736,745p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md | grep -c -i 'oracle\|applicab'
0
```

## Claude 2: the comparator's stratum inventory
```
$ grep -c '^{"kind":"stratum_inventory"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out
660
```

```
$ grep '^{"kind":"stratum_inventory"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -c 'unresolved coverage unless sampled in paired comparison'
658
```

```
$ grep -h 'INVENTORY_LINES\|EMPTY_INVENTORY_LINES' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r2_claude_evidence/enum/inventory.stdout
INVENTORY_LINES 660 total strata 19534
EMPTY_INVENTORY_LINES 2
```

```
$ grep -c 'stratum_inventory' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_record_measurements.py
/var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md:0
/var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md:0
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:0
/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_record_measurements.py:0
```

```
$ grep -c 'unresolved coverage unless sampled' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
/var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md:0
/var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md:0
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:0
```

```
$ sed -n '553,554p;585p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
**Every remaining coverage reason and every affected object is retained in M7**, verbatim from the
accepted coverage appendix, including the full no-value/local lists. The two remaining `not_formed`
- Domain-stratum inventory order can differ across runs; the accepted full-02 itself is filed.
```

```
$ grep -h 'stratum_inventory' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r2_grok_evidence/recount_coverage_stdout.txt | head -2
KINDS {'configuration': 1, 'operands': 330, 'comparison': 330, 'accounting': 423, 'stratum_inventory': 660, 'completion': 1}
```

```
$ grep -n -i 'inventory' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r*_review_disposition.md | cut -c1-200
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r4_review_disposition.md:35:| b | No test covers the Wolfram reader's operator precedence or chained comp
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r4_review_disposition.md:36:| c | The order of the domain-stratum inventory differs between runs; the mul
```

