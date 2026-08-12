# S11 SymPy engine — build directive

## Authority and boundary

Rewrite `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` in full. Its products are
the flushed stdout tag stream and `research/pde_ledger_v3/scripts/S11_exports.py`. The current S11 engine,
its as-built copy and prior stdout are not build premises.
Those products are its only writes; every other file, including `S10_exports.py`, is read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics authority
and wins every conflict. Implement its §§1–10 directly. In particular, point to rather than duplicate §4;
§5, including its corollaries, no-verdict rule and locus protocol; §6 Q1–Q11; §7; §8; §9; and §10. Add no
expected value or acceptance criterion (`CLAUDE.md` rule 5).

Do not add parallel machinery for those contracts: no per-cell completion registry and no directive-local
exit policy. Section 5 corollary 4 applies as a property, without named exceptions.

## Accumulated export

Import `LEDGER` from `research/pde_ledger_v3/scripts/S10_exports.py` for §Q6r and as the incoming chain
record. Write `S11_exports.py` as the accumulated model-definition and knob record that the next step can
import. Its S11 candidate population is every tag emitted for the package §7 identifies as primary. Carry
forward that `LEDGER`; for every written S11 row, `class` is a property of the §8 named object the row
carries, classified from §§8–9's declarations. Do not author a per-quantity class map.

`research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` **F9** is the sole key-write rule.
Apply it; do not restate it or substitute an earlier export proposal.

Do not specify a tag-to-key transform. Emit as computed objects each candidate row's engine-derived key and
the imported `LEDGER` keys having those spellings. They are `CLAUDE.md` rule 2 operands, not a claim or
count; add no guard, threshold or assertion over them.

Choose **F6**'s first branch: publish `S11_exports.py` only if every declared cell for the package §7
identifies as primary completed. Use §7's observed run record; add no completeness field.

F9 is not an acceptance harness. Its file's **OWED TO THE BUILD REVIEW** checklist belongs to the review
legs run against the built script. Do not encode that checklist as guards, registries or invariants in the
engine, and do not make a check audit the input from which it was built.

The S10 exporter is mechanical precedent, not authority. Two inherited facts prevent verbatim reuse:

- `S10_brane_mode_spectrum_sympy_audit.py:1906-1910` makes `exact_reconstruction_match` a `type` plus
  `srepr` spelling comparison; it is not F9's equality proof.
- `S10_brane_mode_spectrum_sympy_audit.py:2087` hard-codes the predecessor step name.

Use §10 to report any conflict, ambiguity or unavailable construction. Do not fill such a gap with new
physics, an expected result, or a self-review mechanism.
