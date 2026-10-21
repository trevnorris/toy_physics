# Measurements: S14a part 1 record decision list setup lookups (generated 2026-10-21 14:48)

Generator: `_scratch/s14a/gen/rec_setup_lookups.sh` (grep/sed/git/ls/test only; this file is written by it).
Mechanical lookups only. Paths are relative to `research/pde_ledger_v3/` unless they start with `../`.

## L1. The accepted commits the list cites

```
$ cd /var/projects/toy_physics && for c in e9eed614 94fddf2f de81b9fc ac909aaf e070c7f3 b08673a2; do git log -1 --format='%h %s' $c; done
e9eed614 Add the S14a part 1 decision list (G2: Codex + Grok, folded once; cleared)
94fddf2f Accept the S14a part 1 build (cleared in scoped round 2)
de81b9fc Accept the S14a part 1 comparator (round 1 clear) and file its production output
ac909aaf Accept the S12 record and its register pass (user: accept, note both)
e070c7f3 Accept the light-sector closeout (review clear)
b08673a2 Clear the S12 decision list revision 4 (scoped pass 5)
```

## L2. Sources: present or absent

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3 && for f in directives/S14a_1_decision_list.md scripts/S14a_1_projection.py mathematica/S14a_1_projection.wl scripts/out/S14a_1_projection.out mathematica/out/S14a_1_projection.out directives/_measurements/S14a_1_build_runs.md directives/_measurements/S14a_1_build_review_r0_disposition.md directives/_measurements/S14a_1_build_review_r1_disposition.md directives/_measurements/S14a_1_build_review_r2_disposition.md directives/S14a_1_comparator_decision_list.md scripts/S14a_1_cross_engine_comparator.py scripts/out/S14a_1_cross_engine_comparator.out directives/_measurements/S14a_1_comparator_production_run.md directives/_measurements/S14a_1_comparator_run_record.md directives/_measurements/S14a_1_comparator_builder_report.md directives/_measurements/S14a_1_comparator_build_review_r0_disposition.md directives/_measurements/S14a_1_comparator_build_review_r1_disposition.md steps/S12_order_source.md directives/_measurements/S12_record_r3_review_disposition.md directives/S12_record_decision_list.md SUBSTRATE_REQUIREMENTS.md ../../STATUS.md V3_STEP_PLAN.md DEFECT_REGISTER.md directives/O2_SHARED_PHYSICS.md steps/S9b_light_bending.md; do if [ -e "$f" ]; then echo "present $f"; else echo "ABSENT $f"; fi; done
present directives/S14a_1_decision_list.md
present scripts/S14a_1_projection.py
present mathematica/S14a_1_projection.wl
present scripts/out/S14a_1_projection.out
present mathematica/out/S14a_1_projection.out
present directives/_measurements/S14a_1_build_runs.md
present directives/_measurements/S14a_1_build_review_r0_disposition.md
present directives/_measurements/S14a_1_build_review_r1_disposition.md
present directives/_measurements/S14a_1_build_review_r2_disposition.md
present directives/S14a_1_comparator_decision_list.md
present scripts/S14a_1_cross_engine_comparator.py
present scripts/out/S14a_1_cross_engine_comparator.out
present directives/_measurements/S14a_1_comparator_production_run.md
present directives/_measurements/S14a_1_comparator_run_record.md
present directives/_measurements/S14a_1_comparator_builder_report.md
present directives/_measurements/S14a_1_comparator_build_review_r0_disposition.md
present directives/_measurements/S14a_1_comparator_build_review_r1_disposition.md
present steps/S12_order_source.md
present directives/_measurements/S12_record_r3_review_disposition.md
present directives/S12_record_decision_list.md
present SUBSTRATE_REQUIREMENTS.md
present ../../STATUS.md
present V3_STEP_PLAN.md
present DEFECT_REGISTER.md
present directives/O2_SHARED_PHYSICS.md
present steps/S9b_light_bending.md
```

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3 && for f in ../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md ../pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md ../pde_ledger_v2/notes/parameter_register.md ../pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md ../pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md ../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md; do if [ -e "$f" ]; then echo "present $f"; else echo "ABSENT $f"; fi; done
present ../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md
present ../pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md
present ../pde_ledger_v2/notes/parameter_register.md
present ../pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md
present ../pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md
present ../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md
```

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3 && git annex whereis scripts/out/S14a_1_projection.out mathematica/out/S14a_1_projection.out scripts/out/S14a_1_cross_engine_comparator.out 2>&1 | grep -E '^whereis|gin|here' 
whereis scripts/out/S14a_1_projection.out (2 copies) 
  	0449db84-a845-4397-ae2d-e50cac28a376 -- git@c44f82230a0c:/data/repos/trevnorris/toy_physics.git [gin]
   	654ba0a9-e5b3-4634-b74d-dcfd28d91b45 -- trevnorris@pop-os:/var/projects/toy_physics [here]
whereis mathematica/out/S14a_1_projection.out (2 copies) 
  	0449db84-a845-4397-ae2d-e50cac28a376 -- git@c44f82230a0c:/data/repos/trevnorris/toy_physics.git [gin]
   	654ba0a9-e5b3-4634-b74d-dcfd28d91b45 -- trevnorris@pop-os:/var/projects/toy_physics [here]
whereis scripts/out/S14a_1_cross_engine_comparator.out (2 copies) 
  	0449db84-a845-4397-ae2d-e50cac28a376 -- git@c44f82230a0c:/data/repos/trevnorris/toy_physics.git [gin]
   	654ba0a9-e5b3-4634-b74d-dcfd28d91b45 -- trevnorris@pop-os:/var/projects/toy_physics [here]
```

## L3. The deliverable paths do not exist yet

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3 && for f in steps/S14a_1_projection.md steps/_measurements/S14a_1_record_measurements.md scripts/S14a_1_record_measurements.py; do if [ -e "$f" ]; then echo "present $f"; else echo "absent $f"; fi; done
absent steps/S14a_1_projection.md
absent steps/_measurements/S14a_1_record_measurements.md
absent scripts/S14a_1_record_measurements.py
```

## L4. The step's decision list: its sections

```
$ grep -nE '^## ' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_decision_list.md
14:## Filing and framing
26:## Supplied equations (committed; unfalsifiable in this build)
53:## User decision on what `J` counts (2026-10-21)
61:## What the build computes
107:## Model point
114:## Inventories (for the record, under its own review; not computed, not used by the builder)
134:## Retained OPEN
149:## Intended claim
154:## Next
```

```
$ sed -n '/^## Intended claim/,/^## Next/p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_decision_list.md
## Intended claim
What the committed law gives for the slab-projected balances. For a throat with spherical symmetry: its drain rate,
the two channels that make up that rate, the order lost over the same shell, and its flux at any outer radius, all
with every field live; and, where the projected density is nonzero, the volume flux of its projected velocity. It is not a mass, a drain law, a throat solve or a flow.

## Next
```

## L5. The plan's S14a row

```
$ sed -n '619,621p;626,641p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
### S13 · The drain rate `J`, and what it is not
`J` is the invariant label. ⛔ **`m_defect` is a GAP** — `INFLOW_MASS_SOURCE_MISSING`; only a dimensional
bridge `α_J ħJ/c_γ²` exists, and pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED`.
### S14a · ⭐ RUN AFTER S13, BEFORE S14 · ⛔⛔ THE DRAIN BRIDGE — BLOCKING, and created by fixing S12

**Round 3 found that correcting S12 severed the chain to the imported gravity results.**

The committed dynamical `Γ_B` law **conserves total material and converts only order**
(`ledger_stage006_two_phase_chiB_ontology.md:81`). But the imported force derivation assumes a **Gauss
drain removing number flux** (`ledger_stage002_matter_stress_force_assembly.md:53,:111`), and S14/S20
use frozen DC sink/return completions. `Q = Θ_Q J / N_∞,3` *labels* the far-field flux; it does **not**
derive it from `Γ_B`.

⇒ **Derive, from the dynamical system:** the projected order-loss source · the definition of `J` · the
profile-dependent `J → Q` map · the controlled `ω → 0` return law.

⛔ **Until this succeeds, S14, S15 and S20 are CONDITIONAL, not carried results.** ⚠ Mark them so in
their records. ⭐ This is a real physics obligation produced by removing a frozen-wall shortcut — exactly
the kind of debt this ledger exists to surface, not a bookkeeping chore.
```

## L6. The build's acceptance: limits and N1

```
$ grep -nE 'not tested|shared by both operands' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_build_review_r1_disposition.md /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_build_review_r2_disposition.md | cut -c1-220
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_build_review_r2_disposition.md:36:    bulk-direction pieces and the face-concentrated parts are shared by both operands and not tested (K4)
```

```
$ sed -n '/^## Known item, carried open/,/^## Not findings/p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_build_review_r2_disposition.md
## Known item, carried open

**N1. The run record's construction-order sentence is stale.** (Orchestrator, from K8; neither leg raised it.)
- **The sentence.** Line 5 of `directives/_measurements/S14a_1_build_runs.md`: "The Wolfram source was written before
  the SymPy source." It is present at `0ee93e9b`, at `eaccd417` and now (K8).
- **What is true.** The first build's author wrote it, about the first Wolfram engine, and there it was true. The
  current Wolfram engine was rewritten by a separate author after repair 1 of the SymPy engine had finished (K8). Its
  log shows no read of the SymPy script (round 1's lookups, K1).
- **Standing.** It changes no computed object. It is a provenance claim, so the S14a part 1 record states the actual
  authorship and order from these dispositions and does not repeat the sentence. The record's decision list also
  orders the run record's line corrected.

## Not findings
```

```
$ grep -nE 'Three authors contributed' -A4 /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_build_review_r2_disposition.md
82:Three review rounds preceded acceptance: two full, one scoped. Three authors contributed:
83-- the first builder wrote both engines;
84-- the same author repaired the SymPy engine twice;
85-- a fresh author rewrote the Wolfram engine, then repaired it once.
86-
```

```
$ sed -n '5p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_build_runs.md
The Wolfram source was written before the SymPy source. Neither engine reads the other engine or its output. No cross-engine comparison or scientific interpretation was performed. Formal integrals retain the live fields; velocity objects carry their nonzero-density domain. All contract OPEN items remain OPEN.
```

## L7. The comparator's acceptance: what is carried to the record

```
$ sed -n '/^## Carried to the S14a part 1 record/,/^## Operations/p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_build_review_r1_disposition.md
## Carried to the S14a part 1 record (not findings)

- **Exact residuals not literally 0 or False** (E9).
  - Seven integral parts: b50, b72, b128, b148, b171, b173 and b174. For each, the exact residual as computed is
    printed, and the samples beside it are the resolved check. A derivative moved to another point moves these
    samples through the quadratic form only (E3).
  - b165, the shell domain: dispositioned in round 0.
  - b188 and b189, the projected-velocity and volume-flux dimensions: their "external printed domain contexts" part
    is not literally zero. SymPy's operand there is `True`. Wolfram's is `Ne(0, Integral(...))`: the projected
    density at the outer radius is nonzero (E9).

  The record interprets these.
- **Completeness.** That every object both engines print and the S14a part 1 decision list names is compared rests on
  the binding table's curation. The 343 curated reasons say so (E4).
- **The endpoint-delta convention** (item 4) is checked by the tests only: the streams contain no closed-endpoint
  integral (Grok, E9).

## Operations
```

## L8. The S12 acceptance: known items and stale status text

```
$ sed -n '/^\*\*Known items, carried open/,$p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S12_record_r3_review_disposition.md
**Known items, carried open.** Neither changes a claim:
- **J1.** The record's account of its own edits does not mention two status-only edits: the S9b oracle record's
  status line and the register header's light-closing sentence. Both are as the record decision list's item 9
  orders, and both are visible in `git diff 560c6c51`.
- **J2.** The "18 → 18" sentence for the round-3 repair cites M19/M20. Its count is in M12.

**Stale status text.** The accepted files still describe S12's record and register pass as written for, or under,
review: the record header, the register's status and S12 pass heading, and `STATUS.md`. The S9b record's acceptance
left the same kind of text, and S12's item 9 then repaired it for the light closeout. The next step's record list
updates these lines to this acceptance, and picks up J1 and J2.
```

```
$ sed -n '3,7p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md
**Status, 2026-10-20:** record and register pass written for review. The constructions are accepted at
`a3da43d9`; the comparator and its filed production output are accepted at `fe485b73`, with the engine-row
note below. This is the round-3 repair of reviewed baseline `fa5d1396`, awaiting the narrow check of H1/H2
(`directives/_measurements/S12_record_r2_review_disposition.md`, retrieved in M19). No drain law, throat solution or flow is delivered.
The light-sector closeout, including the flow-only oracle-comparison record, is accepted at `e070c7f3`.
```

```
$ grep -nE 'written for review|for review|awaiting' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-200
6:split, is accepted at `e070c7f3`. The S12 record and register pass (2026-10-20) are written for review;
229:  proof. ⚠ This is a **premise awaiting execution**, not a postulate: it is the kind of claim a CAS can
1152:additional separate γ test. This is not a count awaiting WL field-equivalence work (M17).
1298:### S12 pass — local order-source compatibility and separate inventories · 2026-10-20; written for review
1302:`fa5d1396`, awaiting the narrow check of H1/H2 (M19).
```

```
$ grep -nE 'written for review|for review|awaiting' /var/projects/toy_physics/STATUS.md | cut -c1-200
1:# Project status — light sector closed for now; S12 record written for review
94:**S12's record is written for review:** its accepted constructions and comparator report the source
```

```
$ grep -nE '18 → 18|M19/M20' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-200
/var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md:333:The round-2 repair of `4c37473c` was also 18 → 18, with no additional source gain. It restored the sourced
/var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md:340:The present repair of `fa5d1396` is also 18 → 18, with no additional source gain. It carries the
/var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md:343:S1's gate decides. It also supplies the accepted-pass edit account below (M19/M20).
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1310:The round-2 repair of `4c37473c` was also 18 → 18, with no additional source gain. It restored the sourced
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1315:The present repair of `fa5d1396` is also 18 → 18, with no additional source gain. It adds the
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1317:it changes no owner, computed object or comparison claim (M19/M20).
```

## L9. The register: entries and passes

```
$ grep -nE '^## |^### ' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-150
23:## Why this exists
38:## The failure this prevents
49:## Entry schema
72:## Entries
103:### R-S6-01 — `B_comp` must be retired or re-affirmed
118:### R-S7-01 — the slab-width flat direction, under a non-uniform width
139:### R-S1-01 — the brane's spatial dimension
177:### R-S1-02 — the substructure's shear response, as a function of the order parameter
216:### R-S1.5-01 — the GNLS's own mode content
246:### R-S6-02 — the reference state must be an equilibrium
279:### R-S8-01 — the form of the brane's quadratic stiffness functional
324:### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order
368:### R-S8-03 — the SIGN of the physical transverse stiffness
406:### R-S8-04 — what carries the brane's internal angular momentum
440:### R-S8-05 — the frame the brane's rotational stiffness is measured against
467:### R-S1-03 — the substructure's microscopic time-reversibility
484:### R-S8-06 — material displacement and the slab's quadratic inertia
518:### R-S12-01 — the reservoir and its power budget
550:### R-S12-02 — background drain/return and separate boundary data
612:### R-O2-01 — the net-power supplier and its stated budget
638:### R-S21-01 — conditional deflection and radar-slope profile requirements
829:### R-S22-01 — how the energy a throat releases splits among light, the extra wave and the bulk
894:## Population passes
915:### Pass 2 — S11, S11b-A/B and S11c PARTIAL · populated 2026-10-05; accepted `22be0137`
1021:### O2 pass — conditional live steady balance · populated 2026-10-07; accepted with the O2 record, `72866fcf`
1122:### S9b pass — conditional light bending and radar slope · populated 2026-10-19; accepted with the S9b record, `ff3bac63`
1250:### Light-closing pass — the oracle comparison and the extra wave · 2026-10-19; accepted `e070c7f3`
1298:### S12 pass — local order-source compatibility and separate inventories · 2026-10-20; written for review
1443:### Method for population passes
1480:## When this is done
```

## L10. The S12 record's owners line the step's decision list cites

```
$ sed -n '292p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md
| Interior-order law/history/strain; faces' law; source/flux forms | No separate owner named for the first two; source/flux forms remain S12's. Slab material completion remains separately S8/S22. None of the forms is delivered here (M01, M11, M18). |
```
