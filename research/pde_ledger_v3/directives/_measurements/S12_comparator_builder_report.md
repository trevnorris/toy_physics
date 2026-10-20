# S12 comparator builder report — 2026-10-20

Built [scripts/S12_cross_engine_comparator.py](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py) and [scripts/test_S12_cross_engine_comparator.py](/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py); contract `3822b64e`, accepted inputs `a3da43d9`. No commit.
All evidence below is under `/var/projects/toy_physics/_scratch/s12/cmp_work/`; the final output is [comparison_final_02.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/comparison_final_02.jsonl).
Name table: [name_table.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/name_table.jsonl) lists all 51 injective rows, each object, Python line citation and Wolfram line citation; these are also printed in the output configuration.
Component correspondences: [component_correspondences.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/component_correspondences.jsonl) lists all 369 joined-row bindings, both engine-line citations and their semantic schemas.
Every joined row without a formed comparison: [uncompared_rows.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/uncompared_rows.jsonl) names each row and reason; the output retains both original tagged operands.
Uncompared rows comprise 46 expression-tree-dependent nested dimension ledgers; their cross-engine node correspondence is not supplied.
The two `*/material_elimination_operand` rows could not be formed: a SymPy scalar bulk multiple versus Wolfram time-derivative rules.
Per-row accounting and both parsed/actually-compared leaf counts: [row_accounting.tsv](/var/projects/toy_physics/_scratch/s12/cmp_work/row_accounting.tsv) and [row_accounting.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/row_accounting.jsonl); every unjoined occurrence has its engine and reason.
Counts (literal summary): `{"unjoined_py":84,"unjoined_wl":81,"joined":369,"joined_with_not_compared":46,"joined_with_not_formed":2,"joined_below_parsed_count":144}`.
Unpaired current time-density components and term-ledger components are all listed in [component_coverage.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/component_coverage.jsonl); no position-only or differently named pairing is made.
Distributions are extracted from each printed operand into three regular-region parts and normalized coefficients restricted to the two face supports. Unrestricted coefficient operands are labelled representation diagnostics.
Exact operands/residuals are emitted before domain restriction; the as-computed and supplied-domain residuals remain side by side, with rational evaluations in addition.
Test functions and points: [test_functions_and_points.jsonl](/var/projects/toy_physics/_scratch/s12/cmp_work/test_functions_and_points.jsonl) prints every live function for affine, quadratic and rational forms, three rational base points, all three regions and both faces.
Strata: generic; each face with density, velocity and flux jumps vanishing; each face with interior order vanishing; and each face with both conditions. Profiles stay live and argument derivatives are separated, including at the order-zero faces.
Printed evaluation rows: 139797; unresolved evaluation rows: 0. Every printed probe-domain value is retained in the output and manifest.
Machinery reused from S9b: the InputForm tokenizer/Pratt reader and independent parse-census design; native sstr parsing, multiargument derivatives/Subs, bindings, distribution extraction, domain restriction and probes are S12-specific.
Re-verified on the accepted streams: native framing/grammar, unique-name joins, semantic current schemas, term-keyed dimension ledgers, support extraction, original operands, complete tag accounting and parsed-versus-compared leaves.
Commands and literal guard outcomes: [research/pde_ledger_v3/directives/_measurements/S12_comparator_run_record.md](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S12_comparator_run_record.md); the manifest also records both full commands and output hashes.
Tests: 30 synthetic tests, 0 errors, 0 failures; guard exit 0, verified limits, no kill. Full comparison: guard exit 0, verified limits, no kill.
Tests: runtime 7.076304531015921 s, process peak RSS 72688 KiB, cgroup peak 71036928 bytes.
Full comparison: runtime 225.47244495898485 s, process peak RSS 145548 KiB, cgroup peak 195858432 bytes.
Leaf counts include syntax heads; comparing unit entries does not credit their matrix wrapper head.
Open coverage is the explicitly listed component/ledger correspondence gaps. Interpretation of residuals belongs to the S12 record; no new physics or physical rate/flow profile was chosen.
No residual target was given. No review rounds were launched. Both source engines and their accepted transcripts remain read-only.
