# Independent review — the ledger fold + under-export guard module

You are one of two independent review legs. The artifact is **Codex-written infrastructure code**. Derive the
contract yourself, then adversarially try to break it. This is a **guard**: its whole value is that it lets **no
under-export or mis-merge through** and that its tests **fail on a broken implementation**. A guard that passes
its own tests but has a hole is worse than none.

## Artifacts
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/ledger_fold.py`
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_ledger_fold.py`

## The contract it must meet (the authority)
`research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md` (`c04e071f`) §D2 (fold) + §D3 (guard),
and the build directive `research/pde_ledger_v3/directives/S11_ledger_fold_build_directive.md`. In one line:
`load_model(base, *deltas)` = chronological **last-wins on exact key**, ⛔ **never re-applies F9**; the guard
takes a consumer's `IMPORT_KEYS` (F9 write-keys) and must catch every under-export via **exact-key manifest
resolution**, **recursive symbol + `dimension_key` edge closure**, an **F9c-ambiguity raise-and-name**, and a
**bidirectional** access-recording smoke-test; plus minimum-mode and a promotion-delta helper. It is
**physics-free** — judge it as chain infrastructure.

## Required method — this is a SCRIPT, so ABLATE and try to hole the guard
Write your own probe scripts; save each script and its literal stdout to named `/tmp` paths and report them. ⛔
A prose "it looks correct" is discarded.

1. **Ablate each load-bearing function; confirm a test FAILS (decisive, not tautological).** In a `/tmp` COPY:
   make `load_model` **first-wins** instead of last-wins; make the closure **non-recursive** (one level only);
   make manifest resolution accept a **bare** key for a declared **prefixed** one; make the smoke-test proxy
   record nothing. Re-run `test_ledger_fold.py` after each and report the literal diff. ⛔ **Any ablation that
   leaves the tests still passing is a tautological test — report which test failed to move.**
2. **Try to construct an under-export the guard MISSES (the real probe).** Build a small synthetic base+delta
   where a consumer's kept row genuinely depends on a row that is **absent**, but the guard passes. Candidates
   to probe: a symbol that appears only **inside a nested container / Matrix / Tuple / relational** in the
   `value` srepr (does `free_symbols` surface it, or does the closure miss it?); a `dimension_key` chain **two
   deep** (does the recursion follow it?); a row referenced by a **string/relational** payload rather than a
   sympy `Basic`; an F9c-collided symbol whose two producers are present so the name is **not** ambiguous but
   the consumer still means the wrong one. If any such under-export passes the guard, that is a must-fix hole.
3. **Try to make `load_model` mis-merge.** Deltas out of order; a delta re-declaring a base key (F9b overwrite —
   is the audit correct?); a key present in two deltas; verify it never re-prefixes or re-compares (no F9
   re-apply). Report a case where the merged result is wrong.
4. **The smoke-test proxy completeness.** Can a lookup evade the access-recording proxy — `.get(k)`, `in`,
   iteration, `.keys()` then index, a cached reference? If a real dereference is not recorded, the bidirectional
   check is defeated. Probe it.
5. **F9c-ambiguity honesty.** Confirm it **raises** (does not silently pick) when a symbol name matches >1 fold
   key, and is a **no-op** when there is no collision. Find a collision shape it fails to detect.

## Physics/impact filter
Report a finding only if it catches a concrete way the module would (a) let an under-export reach a downstream
build, (b) mis-merge the fold, (c) pass a test on a broken implementation, or (d) silently mis-route an F9c
edge. Not style. If a probe finds nothing, say so and give the script+stdout that shows it.

## Ablation sandbox
⛔ Copy the artifacts to `/tmp` and ablate the COPY; never modify the working tree. Save every probe script and
its literal stdout to named absolute paths and report them.

## Output
Findings ranked most-severe first (must-fix = a real hole or a tautological test; nit otherwise), each with the
probe script path, its literal stdout, and the exact line in `ledger_fold.py`. End with a one-line verdict: is
the guard sound (no reachable under-export, all tests decisive), or are there must-fix holes.
