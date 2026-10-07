O2 comparator builder report — 2026-10-07

Built [comparator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py) and [synthetic tests](/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py). No commits were made; starting HEAD was `fdaa575b76df8746daef1372ea06308b81513fb9`.
The accepted engine sources were read, never run or edited. `datalad get` reported both transcripts already present.

The complete **159-row join table** and **43-row name table**, with both construction-line citations and spec objects, are in [catalog.json](/var/projects/toy_physics/_scratch/s11c/o2-comparator-build/catalog.json).
[accounting.jsonl](/var/projects/toy_physics/_scratch/s11c/o2-comparator-build/accounting.jsonl) lists every row's independently parsed and actually subtracted leaf counts, plus all unjoined occurrences and their reasons.
All 159 declared joins resolved. The census covers 42 SymPy and 12 Wolfram tagged objects, with no unaccounted leaves. There are 48 SymPy and 154 Wolfram unjoined occurrences; these include publication metadata, extra copies, and representation-specific content. These counts are accounting, not a physics conclusion.
Complete operands, mapped OPEN structures/orientations, structural differences, and three-valued comparisons are in [comparison.jsonl](/var/projects/toy_physics/_scratch/s11c/o2-comparator-build/comparison.jsonl). A compared count below its parsed count remains explicit.

Machinery adapted from S11c-b: next-tag multiline Wolfram framing, per-row accounting, and bounded-memory execution. No old comparator residual, parser name maps, time alarm, or verdict logic was imported.
Re-verified on these streams: complete framing/grammar, lossless function metadata, every selected path, injectivity, and full leaf partition. The new readers retain held heads, binders and applied arguments; scalar algebra uses exact subtraction. No numeric sampling or constitutive closure was introduced.

Final test command:
```sh
python scripts/s11c_guarded_run.py --pool o2-comparator --memory-gib 6 --log-directory /var/projects/toy_physics/_scratch/s11c/o2-comparator-build/tests-repoint -- python research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py
```
Full comparison command:
```sh
python scripts/s11c_guarded_run.py --pool o2-comparator --memory-gib 6 --log-directory /var/projects/toy_physics/_scratch/s11c/o2-comparator-build/full -- python /var/projects/toy_physics/_scratch/s11c/o2-comparator-build/pinned/O2_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/O2_live_balance_sympy_audit.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/O2_live_balance_mathematica_audit.out --output /var/projects/toy_physics/_scratch/s11c/o2-comparator-build/comparison.jsonl --accounting /var/projects/toy_physics/_scratch/s11c/o2-comparator-build/accounting.jsonl --catalog /var/projects/toy_physics/_scratch/s11c/o2-comparator-build/catalog.json
```
All 26 synthetic tests passed, including repoints of all 43 name rows between existing declared objects. Both final guarded commands exited 0; the full run had zero stderr bytes. Guard receipts verify 6 GiB memory containment, zero swap and no runtime/inactivity deadline. No kill or admission refusal occurred.
Final tests: guard wall `0.49434850085526705` s; cgroup peak `48275456` bytes.
Full run: guard wall `118.84675127104856` s; comparator runtime `118.37262136489153` s; process peak RSS `1328056` KiB; cgroup peak `1362120704` bytes.
Commands, literal child stdout/stderr, resource receipts and source hashes: [measurements](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/O2_comparator_builder_measurements.txt); compact [summary](/var/projects/toy_physics/_scratch/s11c/o2-comparator-build/summary.json).

**No residual target was given.** Unsupported, boolean and structural content remains a coverage finding, not a parse failure or an agreement claim.
Open questions: interpretation of the printed differences and coverage findings belongs to sub-step 7; independent build review belongs to the orchestrator. No review rounds or downstream work were launched. Stopping here.
