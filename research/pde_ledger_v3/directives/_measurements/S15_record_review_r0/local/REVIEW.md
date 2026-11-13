# S15 step record, independent source-first review, frozen round 0

**Verdict: changes required.** Six candidate files frozen in `../candidate.sha256` against `83237f4a`. I formed the source reading before opening the proposed record. The evidence supports a signed control-surface momentum balance, the geometric and live-response separation tests, and a conditional S21 filing. The sources supply neither the momentum delivered to the excluded core nor an inertial response. They do not determine two-drain acceleration.

1. **Momentum density and current, §1 line 11.** Assigning possible stress-carried momentum to `Π_o` misstates P5. A stress contribution to momentum **density** would change `G`, hence the stored momentum and its moving-boundary transport. The open other current is separate. S15 assigns a departure from the extended steady P6 **stress** to that other current. See `directives/S15_decision_list.md:62–73` and `directives/S9b_SHARED_PHYSICS.md:345–355`.
2. **Physical implication, opening line 5.** The phrase “has therefore been established as neither necessary nor sufficient” asserts that both implications were disproved. Neither necessity nor sufficiency was established. The supplied shell and free-mouth observables have no established conversion.
3. **Velocity reading, §2 line 27 and M3.** The candidate defines item 5 and lists cross-engine statuses but omits what the velocity differences actually give for the net crossing. At fixed fields, drain 1's explicit contribution is zero in all three components. Drain 2's contribution is the centre transport `∑∮Gⁱ(Ẋ₂·n)dS`. Drain 1 can still affect the flow and regional momentum rate. The instantaneous shell momentum and integrated force density themselves contain no explicit centre velocity when fields are held fixed. These distinctions are central to the requested reading and the comparison with drag. The original decision lists the velocity-dependent terms among the record’s deliverables (`directives/S15_decision_list.md:145–153`).
4. **Prior-art and reference qualifications, §4.** The original S15 decision asks to quote both sides of G7/G9/G13. Quote a short literal oracle and the actual corresponding S15 observable for each, especially the velocity subtraction opposite G13. S14a already prints `v_proj(R,t)=F_r(R,t)/Σ(R,t)` on the nonzero-density domain. Whether density permits area dilution to carry to velocity, and how that maps to S15, are open. Of its 54 samples for the unresolved b174 residual, only 18 can discriminate a changed derivative point.

**Requirements register and ancillary files.** Compared the actual by-target row, source, entry, pass, S16 record and project status to `83237f4a`. Only `R-S21-01` gains S15 provenance. The added test is conditional, contains the response, retains general fields and leaves the free acceleration open. The census stays 25 OPEN and the S16 acceptance update matches `4f907e9d`.

**Reproduction.** From the repository root, `(cd _scratch/s15/record_build/r0/candidate && sha256sum -c ../candidate.sha256)` returned six `OK` lines. The replay tree contains byte-identical copies of the frozen script, register and generated file; its `directives`, `scripts/out` and `mathematica` reference the original inputs. The frozen script independently reproduced the evidence file:

```text
$ python -B _scratch/s15/record_build/r0/local/replay/research/pde_ledger_v3/scripts/S15_record_measurements.py --check
S15_record_measurements.md 11025775 bytes
```

I independently matched every M3 operand locator, tracked final disposition, single-engine object label and effective refinement. For M4 I reset on each native kernel’s named-object start and used only the completed pass for that object. This excludes partial records from failed attempts. Every condition, corresponding restricted quantity, status and unresolved-list value agreed, in same-key occurrence order. The first audit parser retained failed-pass fragments and disagreed on a retried regional rate. The completed-pass parser below resolves that audit error.

```text
$ python -uB _scratch/s15/record_build/r0/local/audit.py
M4 144 eligible native objects {'condition': 5586, 'Q_under_condition': 5586, 'unresolved_parameter_branches': 42, 'status': 144}
M4 Wolfram full sha256 7d8044652316b879e7d0c61495bcf6fd241ed243009b251687cf5353b2b710fb
both_surfaces/sum_crossing/component_1 drain1 py_reading:177374202 0
both_surfaces/sum_crossing/component_1 drain2 py_reading:177466790 2750-byte expression
both_surfaces/sum_crossing/component_2 drain1 py_reading:184458547 0
both_surfaces/sum_crossing/component_2 drain2 py_reading:184551137 2750-byte expression
both_surfaces/sum_crossing/component_3 drain1 py_reading:191518452 0
both_surfaces/sum_crossing/component_3 drain2 py_reading:191610668 2714-byte expression
M3 336 objects 4176 role/field/operand rows {'precision_checks': 52, 'branch_checks': 7560} {'literal-zero-surface-drain1': 285, 'literal-zero-surface-drain2': 225}
register {'OPEN': 25} -> {'OPEN': 25}
```
