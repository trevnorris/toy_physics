# S11c-b P2b build report — 2026-10-18

Built the §3a expression-valued coefficient bridge, invariant-object bijection certificate, argument-sensitive parse, and required controls. **The family factors and measured facts supplied by the decision list are premises; this build cannot falsify them.** Commits made: **none** (authorization cites `69bcc18d`). No agents, reviews, Wolfram kernels, or committed `.out` rewrites.

Paths below are relative to `/var/projects/toy_physics-p2b`. Script aliases: **H** = `research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py`; **A** = `research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py`; **C** = `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`.

Every run used this exact command template, replacing `LOG` and `SCRIPT` with its table entries:

```sh
python3 /var/projects/toy_physics-p2b/scripts/s11c_guarded_run.py --pool p2b --memory-gib 6 --log-directory /var/projects/toy_physics-p2b/_scratch/s11c/LOG -- python3 research/pde_ledger_v3/scripts/SCRIPT --p2b
```

| Run | SCRIPT | LOG; literal output |
|---|---|---|
| H | `S11c_b_handcoded_comparison.py` | `p2b-h-03`; `_scratch/s11c/p2b-h-03/stdout` |
| A | `S11c_b_adjudicated_comparison.py` | `p2b-a-01`; `_scratch/s11c/p2b-a-01/stdout` |
| C | `S11c_b_cross_engine_comparator.py` | `p2b-c-01`; `_scratch/s11c/p2b-c-01/stdout` |

All three exited **0**, with verified limits and empty stderr; their directories contain invocation, limit-validation, child-outcome, and outcome JSON. No deadlines were imposed. H includes the final complete-candidate emission and full-residual parse checks; A/C preceded those additions to H, which did not change their measured bridge images or routing.

| Requirement/control | Implementation lines | Run and literal stdout location |
|---|---|---|
| B1: expression image; coefficient-only substitution, closed factors | H:269–311,493–535 | H:1090–1869, `coefficient_image`, `family_factor` |
| B2: object matching, both degrees/counts | H:548–624 | H:1–1089, `P2B_PAIRING_*` |
| B3: unfolded invariants and emitted-action certificate on shared apply path | H:633–655; A:1318–1343 | H:1090–1869; A:970–1329, action/Hessian operands and residuals |
| B4: head classification, argument retention, jet decoding | H:198–206,258–311,658–719; C:315–352,391–403 | H:4153–4225; C:1–67 |
| B5: stale protection removed; actual classifier exercised | A:54,1345–1378 | A:1330–1385, `P2B_PROTECTION_*` |
| B6: CAS-only P2b emission, emit before guard; S11c-a unchanged | H:496–500,650–655,658–789; A:1309–1379 | All three stdout files; `_scratch/s11c/p2b-build/output-syntax.json` |
| Wrong pairing, locked factor | H:722–756 | H:4057–4152, `P2B_CONTROL_pairing_*` |
| Missing coefficient scale | H:735–757 | H:4057–4152, `P2B_CONTROL_coefficient_scale_*` |
| Folded family factor | H:735–757 | H:4057–4152, `P2B_CONTROL_family_factor_*` |
| String image instead of product | H:741–757 | H:4057–4152, `P2B_CONTROL_string_image_*` |
| Base/shifted width, modulus, eW; inert positive fixtures; thetaWave preservation | H:658–692 | H:4153–4213; C:1–55 |
| Both background-jet decode paths | H:693–719 | H:4214–4225; C:56–67 |
| Restored 07/10 protection | A:1345–1378 | A:1330–1385, current/restored route predicates |
| Independent PY/WL permutations; duplicate keys | H:761–787 | H:1870–4056, `P2B_PERMUTATION_*`, `P2B_DUPLICATE_*` |

Measured: **15/15 per source, 30/30 per branch, 60/60 overall**; each forward/reverse degree is **1**. All 60 certificate invariant/term residuals and adjudicated action/Hessian residuals are zero. Wrong-pair invariant residuals and missing-scale/string-image term residuals are nonzero. Folded-action term residuals are `Integer(0)` while factor locks are `false`. Both permutation equalities and duplicate rejection predicates are `true`. Each live-background residual remains nonzero; inert and thetaWave residuals are zero. All eight 07/10 planted residuals route with predicate `false` under current protection and `true` when protection is restored.

Changed durable files: **H, A, C, this report**. Generated work/log files remain under `_scratch/s11c/p2b-build/` and `p2b-{h-01,h-02,h-03,a-01,c-01}/`; their exact inventory is `_scratch/s11c/p2b-build/files.txt`. H-01 exited 1 on an internal-name/WL-pattern fixture mismatch; corrected engine spelling is covered by the successful runs. H-02 exited 0 before final emission additions. Committed input transcripts and S11c-a source match HEAD; `git diff --check` passed.

Scope/open questions: these runs validate the **P2b package**, not the legacy default CLI reports or P2a's final slab-row validation. No independent build-review clearance is claimed. Those later activities remain outside this authorization.
