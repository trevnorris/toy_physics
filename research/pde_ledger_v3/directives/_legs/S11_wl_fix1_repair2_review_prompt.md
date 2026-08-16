# Scoped review: post-legs provenance repair to a Wolfram engine (script branch)

## Artifact
The delta ONLY: `/home/trevnorris/.s11_build/fix1_build/repair2/post_legs_only.diff` — the
post-review repair applied to the working-tree
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`.
The larger round-1 repair beneath it has already had two full review legs; do NOT re-review it.
Your scope is the delta and its interaction with the engine.

## What to check
The delta claims: (a) operand slots that name the primary QE attempt's return now carry the
primary's ACTUAL outcome (including its budget-expiry failure object), never the fallback's answer
under the primary's name; (b) whenever the route ran more than one attempt, the affected locus
records, component-count records, and constancy certificates carry the full attempt sequence;
(c) status tokens and test objects are computed exactly as before; (d) records whose primary never
expires are byte-identical to before. Governing text: the "Post-legs fold" section of
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round1_brief.md`.
Builder evidence to check against (not to trust): logs under
`/home/trevnorris/.s11_build/fix1_build/repair2/`.

## Required method — script branch, scoped
Ablate; don't just read. Each probe: copy the .wl to /tmp, mutate the COPY, run, report the
literal tag diff.

1. **Dead-payload probe on the NEW emission paths**: corrupt the attempt-sequence payload the
   delta emits (replace it with a sentinel) under a tiny primary budget — the output MUST move
   now; byte-identical output is a BLOCKING regression of the defect this delta exists to fix.
   Check all three record kinds the delta touches (locus operands, count records, certificates).
2. **Provenance honesty probe**: tiny primary budget, uncorrupted — verify the primary-named slot
   carries the expiry failure object and never the fallback's boolean; verify status tokens match
   the normal-budget run's tokens on records the fallback decides.
3. **Byte-identity probe**: normal budgets, one small cell (`MAIN 2`) — byte-identical to the
   committed baseline cell record.
4. **Conditionality**: the attempt sequence's PRESENCE may depend only on how many attempts ran —
   confirm no dependence on the outcome's value beyond that, and no new tag names.

## Operational constraints (identical for both legs)
- ⛔ Every kernel run wrapped in `timeout 600`; never more than ONE kernel at a time.
- ⛔ Never run `XKIN_ANISO 4` or `XKIN_ANISO 2` (they exceed the machine's memory — measured).
  `MAIN 2` (~15 s) suffices for every probe here.
- ⛔ Mutate only /tmp copies; never the working tree.
- ⭐ Save every probe script AND its literal stdout to named absolute paths and cite them.

## Physics filter
Report a finding only if it catches a way the emitted record could be wrong or dishonest.

## Report format
Verdict line; BLOCKING findings with literal evidence; non-blocking; probe table (probe → run →
literal outcome path).
