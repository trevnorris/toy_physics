# Measurements — O2 record directive review dispositions (generated 2026-10-07 21:56)

Generator: `_scratch/s9b_build/gen/o2_record_directive_lookups.sh` (sed/grep/sha256sum only; lines cut at 400 characters).

```
$ sha256sum _scratch/s9b_build/O2_record_directive_reviewed_v0.md
ddcb5ac0164366469db54ad5b7fa161c5e24e93707c8397438b56fbe336cb10c  _scratch/s9b_build/O2_record_directive_reviewed_v0.md
```

## D1 — a build directive in the record's source packet
```
$ grep -n 'build directive' _scratch/s9b_build/O2_record_directive_reviewed_v0.md
28:- **The accepted comparator** (`48c9bf33`) and its build directive `directives/O2_comparator_build_directive.md`.
```

```
$ grep -n -F '| **Step record /' CLAUDE.md
58:| **Step record / `.tex` card / physics-bearing prose** | Two; O or Cx; ⛔ never chosen by extension | Source-first fidelity; quote both sides; ⛔ no build directive in the packet; a measured physics claim still needs its command + stdout; for cards check suppressed macro fields | **Review-until-clear** about what may be claimed | Own artifact/version review required; script review alone is n
```

## D2 — the momentum-density premise in item 6, against the O2 sources
```
$ sed -n '/What S9b Part D may use/,/OPEN (/p' _scratch/s9b_build/O2_record_directive_reviewed_v0.md
6. **What S9b Part D may use.** The record identifies which O2 objects and conditional inputs a later S9b Part D
   calculation may use, and under which conditions. It also records the user's premise of 2026-10-07: the flowing
   brane's momentum density is `ρ_br V`. That premise is adopted for the S9b repair, not for O2, whose spec keeps
   `ℐ_br^live` OPEN (`O2_SHARED_PHYSICS.md:133`).
```

```
$ git show 4680e251:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | sed -n '133p'
| `ℐ_br^live`; **OPEN**, C §3 | Flowing/embedded inertial and kinetic response, entering momentum storage and transport. No identification with `ρ_br V`, a quadratic live kinetic energy, or a constant inertia is supplied. | S8, with S1.5 antecedents when available. |
```

```
$ grep -n -F 'momentum density' research/pde_ledger_v3/directives/O2_premise_decision_list.md research/pde_ledger_v3/directives/_measurements/O2_build_r3_review_disposition.md
research/pde_ledger_v3/directives/_measurements/O2_build_r3_review_disposition.md:23:- the bulk carry, the momentum density and flux, and the internal force;
```

## Fold applied
```
$ sha256sum research/pde_ledger_v3/directives/O2_record_directive.md
17fa670608576d413a2f1bfc5e92072bc1ac44a41cef146182cb6d70e5d9433c  research/pde_ledger_v3/directives/O2_record_directive.md
```

```
$ grep -c 'build directive' research/pde_ledger_v3/directives/O2_record_directive.md
0
```

```
$ grep -c 'ρ_br V' research/pde_ledger_v3/directives/O2_record_directive.md
0
```

