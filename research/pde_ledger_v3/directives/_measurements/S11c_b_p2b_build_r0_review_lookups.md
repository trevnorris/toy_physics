# Measurements: S11c-b P2b build review, round 0 (generated 2026-10-18 17:13)

Generator: `_scratch/backlog_runs/gen_p2b_r0_lookups.sh` (sha256sum/grep/sed/awk line selection/wc only).
Paths: the build checkout `/var/projects/toy_physics-p2b`; the legs' evidence copied to `/var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics-p2b && sha256sum -c /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_baseline_r0.sha256 | grep -c ': OK$')
8
```

```
$ sha256sum /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_prompt_r0.md | cut -c1-64
ce922f61d5f69be8c8594a221f0437ad3cecdb31ab9037dc3fda845335850e2d
```

```
$ git -C /var/projects/toy_physics-p2b diff --stat | tail -1
 3 files changed, 435 insertions(+), 52 deletions(-)
```

## Verdicts
```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_grok.txt
Verdict: CLEAR
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude.md
Verdict: NEEDS REVISION
```

## Finding 1: a certificate disagreement or count mismatch stops the run before any family
```
$ sed -n '594,595p;600,601p;605,615p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
        if len(matches) != 1:
            raise C.InputError("invariant pairing degree")
        if degrees[id(py)] != 1:
            raise C.InputError("reverse invariant pairing degree")
            p2b_emit(f"{tag}_source_counts_{branch}_{source}", counts)
            if counts != sp.Tuple(15, 15):
                raise C.InputError("source cardinality")
        counts = sp.Tuple(*(sp.Integer(sum(row.branch == branch for row in records)) for records in (wl_records, py_records)))
        p2b_emit(f"{tag}_branch_counts_{branch}", counts)
        if counts != sp.Tuple(30, 30):
            raise C.InputError("branch cardinality")
    counts = sp.Tuple(sp.Integer(len(wl_records)), sp.Integer(len(py_records)))
    p2b_emit(f"{tag}_overall_counts", counts)
    if counts != sp.Tuple(60, 60):
        raise C.InputError("overall cardinality")
```

```
$ sed -n '655,656p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
            raise C.InputError("coefficient bridge certificate")

```

```
$ sed -n '1427,1429p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py
        H.install_gamma_bridge(py_tags, wl_tags)
        bridge_a_state = H.gamma_bridge()
        H.emit_gamma_certificate(bridge_a_state)
```

```
$ grep -n 'no assert on the residual' /var/projects/toy_physics-p2b/research/pde_ledger_v3/directives/S11c_b_p2b_gamma_bridge_directive.md
63:B6. Output states operands + residual; no verdict token, no assert on the residual (rule 5). S11c-a NOT regressed
```

```
$ grep -n 'def test_disagreement_is_measurement_exit_zero' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/test_S11c_b_cross_engine_comparator.py
209:    def test_disagreement_is_measurement_exit_zero(self) -> None:
```

```
$ grep -E '^(ind-wlinvariant|ind-wlaction)-a-default' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/EVIDENCE_INDEX.txt | cut -d'|' -f1,2
ind-wlaction-a-default | exit=2 
ind-wlinvariant-a-default | exit=2 
```

```
$ grep -c '' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/runs/ind-wlinvariant-a-default.stdout
482
```

```
$ grep -vc '^P2B_' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/runs/ind-wlinvariant-a-default.stdout
0
```

```
$ grep -vc '^P2B_' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/runs/ind-wlaction-a-default.stdout
0
```

## Finding 2: the controls do not pass through the production certificate
```
$ awk 'NR>=722 && NR<=757' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py | grep -c 'gamma_certificate_objects\|\.apply('
0
```

```
$ sed -n '739,740p;749p;751p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
                    action = pair.factor * pair.wl.action if package == "family_factor" else pair.wl.action
                    factor = sp.S.One if package == "family_factor" else pair.factor
                invariant_residual = sp.expand(pair.py.invariant - factor * invariant)
                factor_lock = sp.sympify(factor == P2B_FACTORS[source])
```

```
$ sed -n '645p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
        factor_lock = sp.sympify(image == pair.image and isinstance(image, sp.Mul))
```

```
$ grep -h 'P2B_CONTROL_family_factor_.*factor_lock' /var/projects/toy_physics-p2b/_scratch/s11c/p2b-h-03/stdout | head -4
P2B_CONTROL_family_factor_LAB_HELD_MU_R_BG_factor_lock = false
P2B_CONTROL_family_factor_LAB_HELD_W_BG_factor_lock = false
P2B_CONTROL_family_factor_MATERIAL_ADVECTED_MU_R_BG_factor_lock = false
P2B_CONTROL_family_factor_MATERIAL_ADVECTED_W_BG_factor_lock = false
```

```
$ grep -E '^mut-foldfactor' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/EVIDENCE_INDEX.txt | cut -d'|' -f1,2
mut-foldfactor-h-p2b | exit=1 
```

## Finding 3: the B4 fixture and the coverage check are built from the tables they test
```
$ sed -n '197,206p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
INERT_APPLIED = {"mu_theta_L", "mu_theta_M"}
LIVE_BACKGROUND_APPLIED = frozenset(C.LIVE_BACKGROUND_HEADS)
# Keep the complete decode domain separate from the removable applications.
BARE_APPLIED = frozenset(TRIAL_TEST_SUPPORT_APPLIED | INERT_APPLIED)
JET_APPLIED = BARE_APPLIED | LIVE_BACKGROUND_APPLIED
APPLIED_HEAD_CLASSES = {
    **{name: "INERT" for name in INERT_APPLIED},
    **{name: "LIVE_BACKGROUND" for name in LIVE_BACKGROUND_APPLIED},
    **{name: "TRIAL_TEST_SUPPORT" for name in TRIAL_TEST_SUPPORT_APPLIED},
}
```

```
$ sed -n '664p;689p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
    for head in sorted(LIVE_BACKGROUND_APPLIED | INERT_APPLIED | {"thetaWave"}):
    coverage = sp.sympify(set(APPLIED_HEAD_CLASSES) == set(JET_APPLIED))
```

```
$ sed -n '315p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py
LIVE_BACKGROUND_HEADS = frozenset({"widthBackground", "modulusBackground", "eWBackground"})
```

```
$ grep -E '^mut-dropwidth' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/EVIDENCE_INDEX.txt | cut -d'|' -f1,2
mut-dropwidth-c-p2b | exit=0 
mut-dropwidth-h-p2b | exit=0 
```

```
$ grep -h 'P2B_HEAD_classification_coverage' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_r0_claude_evidence/runs/mut-dropwidth-h-p2b.stdout
P2B_HEAD_classification_coverage = true
```

## Finding 4: the builder's runs predate its last script edits
```
$ stat -c '%y %n' /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py | cut -c1-19,36-
2026-10-18 16:05:23 /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py
2026-10-18 16:05:23 /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_adjudicated_comparison.py
2026-10-18 16:00:08 /var/projects/toy_physics-p2b/research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py
```

```
$ for d in a-01 c-01 h-01 h-02 h-03; do stat -c '%y %n' /var/projects/toy_physics-p2b/_scratch/s11c/p2b-$d/outcome.json | cut -c1-19,36-; done
2026-10-18 16:03:24 /var/projects/toy_physics-p2b/_scratch/s11c/p2b-a-01/outcome.json
2026-10-18 16:02:34 /var/projects/toy_physics-p2b/_scratch/s11c/p2b-c-01/outcome.json
2026-10-18 16:02:05 /var/projects/toy_physics-p2b/_scratch/s11c/p2b-h-01/outcome.json
2026-10-18 16:03:09 /var/projects/toy_physics-p2b/_scratch/s11c/p2b-h-02/outcome.json
2026-10-18 16:04:19 /var/projects/toy_physics-p2b/_scratch/s11c/p2b-h-03/outcome.json
```

## Guard path in the prompt
```
$ grep -n 'toy_physics-p2b/scripts/s11c_guarded_run.py' /var/projects/toy_physics/_scratch/backlog_runs/p2b_build_review_prompt_r0.md | cut -c1-120
77:  `python3 /var/projects/toy_physics-p2b/scripts/s11c_guarded_run.py --pool p2b-review --memory-gib 6 --log-directory
```

```
$ sha256sum /var/projects/toy_physics-p2b/scripts/s11c_guarded_run.py /var/projects/toy_physics/scripts/s11c_guarded_run.py | cut -c1-16
c747c266e78ed519
c747c266e78ed519
```

