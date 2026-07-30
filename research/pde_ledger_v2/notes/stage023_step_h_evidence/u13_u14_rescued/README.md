# Rescued step-(h) artifacts — the sole evidence for U13 and U14

Copied verbatim (2026-07-30) out of `research/pde_ledger_v2/_scratch/stage023_h/`, which `.gitignore:96`
(`research/**/_scratch/`) excludes and which is being deleted. Nothing here was edited; the `.diff` is the
only generated file (`diff -u` between two of the copied inputs).

⚠ These back **prose-only** claims. Before this rescue, U13 and U14 were narrated in the stage note with no
committed artifact behind either — the artifacts existed only in scratch.

## U14 — joint perturbation; declaration freedom larger than one-at-a-time ablation shows

Cites: `notes/stages/ledger_stage023_nullspace_underdetermination.md` §1.7(4) **U14** (~`:701-710`).

| file | was | is evidence for |
|---|---|---|
| `u14_i7_joint.patch.py` | `_scratch/stage023_h/adversarial/patches/i7_joint.py` | the mutation itself — 16 declarations changed together (`a, c_s, omega, D0, K0c, K_eta, T_Omega, Z0_ret, Z1_ret, Omega_U, Omega_W, R_mix, g_U, g_W, R0, R1`) |
| `u14_i7_joint.stage.out` | `…/adversarial/runs/i7_joint/stage.out` | the stage staying green under it — `AUDIT_STATUS=PASS`, `TALLY sympy: 111 pass + 0 fail`, `OVERALL PASS`, `PHYSICS_VERDICT=FAIL_UNDERDETERMINED_NOT_PREDICTIVE` unchanged |
| `u14_i7_joint.cmp.out` | `…/adversarial/runs/i7_joint/cmp.out` | ⭐ the **scope limit** on U14: `RESULT\|stage=stage023\|status=FAIL\|mismatches=16`. The stage's own gate is blind; the cross-engine comparator is **not**, because the `.wl` carries an independent literal `baseDims` table. U14 is a claim about the Python's dimensional gate, not about the workstream's detection as a whole. |

## U13 — same-class mis-binding invisible to every instrument in this workstream

Cites: `notes/stages/ledger_stage023_nullspace_underdetermination.md` §1.7(4) **U13** (~`:692-700`), `:879`
(*"the same-class case is U13"*), `:904`; and `notes/stage023_py_rewrite_prediction.md:68,81,112` (P2's
falsification and P3's exclusivity half both rest on U13).

| file | was | is evidence for |
|---|---|---|
| `u13_expE_same_class_rebinding.diff` | generated: committed `scripts/ledger_stage023_nullspace_underdetermination_sympy_audit.py` vs `_scratch/stage023_h/fid/expE/…py` | the five same-class rebindings in `dimension_records()` — `M0←R0`, `K0c←Z1ret`, `g_U←g_W`, `eta_null←q_free`, `T0←P0_physical`. ⚠ The diff also carries a one-line `print_verdict_labels()` move: `expE` predates the R2 remediation (`../ABLATION_SUMMARY.md:36`), it is not part of the mutation. |
| `u13_expE.stage.out` | `…/fid/expE/out.txt` | the stage staying green — `AUDIT_STATUS=PASS`, `TALLY sympy: 111 pass + 0 fail`, `OVERALL PASS` |
| `u13_expE.dimensions.txt` | `…/fid/expE/…dimensions.txt` | the emitted payload under the mutation |
| `u13_baseline.dimensions.txt` | `…/fid/before.dimensions.txt` | the unmutated payload to compare it against |

**The byte-identical claim, re-verifiable from these two files:**

```
diff <(grep '^DIM|' u13_baseline.dimensions.txt) <(grep '^DIM|' u13_expE.dimensions.txt)   # empty; 29 records
```

Both reduce to md5 `61b45bda8972601fb19f8abe2820ea92`. Only the `DIMENSIONS|` header differs, and only in its
`source_sha256`/`ledger_dimensions_sha256` fields — which necessarily change when the source is mutated at
all. That is why the mis-binding is invisible: the 29 records a reader or comparator sees are the same bytes.

## Leg reports

| file | was | cited by |
|---|---|---|
| `REPORT_ADVERSARIAL.md` | `_scratch/stage023_h/REPORT_ADVERSARIAL.md` | `manifests/DIMENSION_REWRITE.md:469`, which **quotes** it (*"Positive control on the Mathematica side"*) as the measurement retiring the *"verification agents cannot confirm the reference side"* premise. Also the origin record for U14 (findings 10b / `$S/runs/i7_joint`). |
| `REPORT_REMEDIATE_H.md` | `_scratch/stage023_h/REPORT_REMEDIATE_H.md` | `../ABLATION_SUMMARY.md:36` cites its **R2** for the `print_verdict_labels()` defect. |

## Deliberately NOT rescued

The 51 per-record `.cmp`/`.stdout` captures under `_scratch/stage023_h/captures/` (832K). `../ABLATION_SUMMARY.md`
§Files already declares them not retained, and every per-row conclusion read from them is committed in
`../results.tsv` (51 rows) and `../include_list.tsv`. Conclusion already recorded ⇒ no rescue.
