# Rescued stage044 verdict artifacts

`verdict_py.json` · `verdict_wl.json` — copied verbatim (2026-07-30) from
`research/pde_ledger_v2/_scratch/stage044/` before that gitignored tree was deleted.

**Cited by** `notes/stages/ledger_stage044_parent_action_reconciliation.md:394-396`, under **Governing**,
alongside two artifacts that were **already gone** by the time of this rescue and could not be recovered:
`stage044_synthesis_directive.md` and `OUT_stage044_ablation.txt`.

**Evidence for:** the stage's emitted reconciliation verdict — `P_retirement` (`stage006_drift 6->5`,
`stage007_drift 11->7`, `stage007_DOF 8->4`, net routeless `-5`), the operative/retired stage007 action
summands, and the `S_cons_G0` summand list.

⭐ The two files are **byte-identical**, which is the cross-engine (SymPy vs Mathematica) agreement itself,
not a duplicate.

⚠ **Weaker rescue than the U13/U14 set — this is insurance, not sole evidence.** `verdict_py.json` is a
generated output: `notes/rewrite_reference_table.md:112` records the stage044 `.py` as the corpus's *only*
filesystem writer (`:1343-1346`, `Path(__file__).resolve().parents[3]` → `_scratch/stage044/verdict_py.json`,
`mkdir(parents=True)`), and both engines are committed, so it is regenerable by re-running them.

⛔ **But the write target is the tree being deleted, and relocating the script moves it.** Stage044 is also
frozen/paused pending the 044-v2 un-freeze, so nothing is scheduled to regenerate these. Kept at 4.4K rather
than relying on a re-run that no one currently plans to perform.
