# Independent review — S11 census-instrument round-4 repair diff and its census

## Artifact

The repair diff `git diff fd9a5835..3f229bbf -- research/pde_ledger_v3/reduction/` (working tree
matches `3f229bbf`), the modules themselves under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/`, and the round-4 census outputs under
`/home/trevnorris/.s11_build/census_build4/` (build report `build_report.md`).

## Scope — SCOPED review of the round-4 diff

Full reviews already ran on the round-1 instruments and the round-2 and round-3 diffs; the round-4
acceptance surface is the folded directive
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_census_instr_repair4_directive.md`
(commit `fd38ebba` — read it and its `../_measurements/S11_census_instr_repair4_directive.md` twin; its
single defect class and its two numbered Required semantics were verified by two directive legs plus
orchestrator computation), under the brief (`S11_locus_census_instrument_brief.md`, `fa8c58b3`). Do not
re-litigate what earlier legs cleared outside this diff. Your questions:

1. **Does the repair implement the directive's two Required semantics in the executed census path?**
   Ablate each one, one at a time, in a /tmp copy against a real-record slice and report the literal
   verdict diff.
   - **Required 1 — zero-status of a free-symbol-less constant is decided by PROOF, the exact route
     mandatory.** ZERO only with an exact zero certificate (minimal-polynomial / exact equals-zero on
     algebraic numbers); NONZERO only with an exact nonzero certificate or a rigorous numeric enclosure
     that excludes zero; UNDECIDED only after the exact route was attempted and failed (resource-guarded,
     the guard's expiry recorded). NEVER NONZERO by failure-to-simplify; NEVER UNDECIDED merely because a
     fixed working precision could not separate a nonzero value from zero. This binds EVERY consumer of
     the status (coverage, membership, spurious refutation, probe refutation, definedness) — verify the
     binding reaches all of them, not just coverage. Ablate: swap the exact-certificate classifier back
     to the structural `simplified != 0` test and show the verdict diff.
   - **Required 2 — sampled coverage refutes only on a certified-uncovered point, under per-branch AND.**
     A DEFINED branch covers a sample point iff EVERY one of its substituted constraints is certified
     ZERO under Required 1; the point is certified-uncovered iff EVERY defined branch has at least one
     constraint certified NONZERO. Constraints belonging to ONE branch are NEVER multiplied together —
     that turns the branch's AND into an OR and manufactures coverage. `COVERAGE_UNDECIDED` keeps its own
     token and is never collapsed; a sample refutes coverage, it NEVER proves coverage (so `COVERED_SAMPLED`
     / `COMPLETE_FACTOR_COVER_SAMPLED` must not survive as production coverage verdicts). Ablate: multiply
     the within-branch constraints (reintroduce AND->OR) and show which candidate becomes falsely covered.
2. **Did any repair weaken a certified direction?** The directive's closing constraint lists them (the six
   round-3 classes, the seven round-2 classes, reducer object-level arithmetic, population reconciliation,
   certified-nonzero refutations, `PROVED_CONSISTENT`, the two-route residual comparisons). The genuine
   refutation direction MUST stay: a point where the union product is certified nonzero is a witness of
   non-coverage and must still be verdicted `NOT_COVERED_SAMPLED`; the round-3 spurious and probe
   refutations both legs upheld must not change. Run one-sided corruption on any two-route comparison the
   diff touches.
3. **Is the round-4 calibration real?** Four new plants, each demonstrated ABLE TO FAIL before the census,
   through production entrypoints and byte-shaped:
   - a free-symbol-less nested-radical exact zero routed through the production status classification —
     must NOT be NONZERO (specified to fail at `fd9a5835`);
   - a free-symbol-less algebraic constant nonzero but smaller than any fixed working precision (shape
     `Sqrt[10^200 + 1] - 10^100`) — must be certified NONZERO by the exact route (fails under a
     fixed-precision-only implementation);
   - a solution/equation pair whose candidates are inverse charts of the emitted branches with a union
     product too large for the symbolic simplify route — must NOT be omitted (fails at `fd9a5835`);
   - a genuinely uncovered candidate whose union product is certified nonzero at a sampled point — must
     still be `NOT_COVERED_SAMPLED` (guards the refutation direction).
   Verify at least the exact-zero radical plant and the tiny-nonzero plant actually fail as specified:
   run the round-4 calibration case against an `fd9a5835` copy of the instruments in /tmp (and, for the
   tiny-nonzero one, against a fixed-precision-only variant) and show the miss.
4. **Are the changed verdicts trustworthy?** The build report claims the round-4 census changed EXACTLY
   four verdicts, all completeness-level, and nothing else (per-key comparison against round 3):
     - terminal metrics R3 -> R4: WL omitted records 3 -> 2; reducer findings 348 -> 347 (omitted_branch
       73 -> 72; spurious 171 and witness 104 unchanged); limitations 815 -> 816; completeness_undecided
       6 -> 10; complete_factor_cover_sampled 3 -> 0; failures 917 unchanged; final verdict `ROUND_FAIL`
       unchanged. Every branch, spurious, witness, probe, `PROVED_CONSISTENT`, resource-expiry, parser,
       and residual-comparison verdict claimed per-key IDENTICAL.
     - the four changed completeness verdicts (re-decide EACH independently from the record's own operands):
       (a) WL `..._D2_ROOT_COINCIDENCE_COEFF_SOLUTION`: `OMITTED_BRANCH` -> `COMPLETENESS_UNDECIDED`;
       (b) WL `..._D3_ROOT2_RANK_DROP_K_SOLUTION`: `COMPLETE_FACTOR_COVER_SAMPLED` -> `COMPLETENESS_UNDECIDED`;
       (c) WL `..._D3_ROOT3_RANK_DROP_K_SOLUTION`: `COMPLETE_FACTOR_COVER_SAMPLED` -> `COMPLETENESS_UNDECIDED`;
       (d) PY `..._D4_ROOT1_RANK_DROP_K_SOLUTION`: `COMPLETE_FACTOR_COVER_SAMPLED` -> `COMPLETENESS_UNDECIDED`.
   Sample at least five changed/at-risk verdicts per census and re-decide each with your own independent
   script. A verdict changed WITHOUT the round-4 repair class behind it is a finding of the first order.
   Check specifically:
   - (a) the D2 flip: re-decide from the record's own emitted branches and candidates — are its candidates
     inverse charts of the emitted branches, and does the exact zero-status leave the union product
     genuinely undecided (not falsely omitted, and not falsely covered)?
   - (b)-(d): are these demotions justified — i.e. was each round-3 `COMPLETE_FACTOR_COVER_SAMPLED` a
     coverage claim resting only on samples (which never prove coverage), so `COMPLETENESS_UNDECIDED` is
     the honest verdict, with no genuine algebraic cover silently lost?
   - the D3/D4 `ROOT_COINCIDENCE` records must STAY `OMITTED_BRANCH`, and each carries TWO genuinely
     omitted families. Verify BOTH families of each are genuinely uncovered under per-branch AND — the
     round-3 agent leg's "cand[1] covered" claim was exactly the AND->OR artifact (it multiplied
     within-branch constraints so `bComp*0 = 0`); re-decide each family with EVERY branch constraint
     required to vanish simultaneously at a generic point of the family, and confirm no defined branch
     covers it.
5. **Reducer**: recompute failures/findings/limitations from the four round-4 census stdouts with your own
   script (object-level counting — sheet lines are evidence, not objects); any discrepancy is a finding.
   Verify the taxonomy stays closed over every token actually emitted, and that `findings = 171 spurious +
   72 omitted + 104 witness = 347` reconciles at object level.

## What you are handed

The round-4 directive + twin; the brief; the spec
(`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` §5, Q8a/Q8b); both
committed records (WL `a4cf6539`, PY `19591194`); round-3 outputs under
`/home/trevnorris/.s11_build/census_build3/` (delta baseline; `build_report.md` §5-§9); round-4 outputs
under `/home/trevnorris/.s11_build/census_build4/` (the canonical run is the unsuffixed `*.stdout`; the
`.parallel.*` and `.pre_exact_predicate.*` families are the builder's own preserved audit runs and are
NOT the reviewed outputs); `DEFECT_REGISTER.md` (7 entries, untouched — census rediscovery is
corroboration, not a repair target).

## Required method

These artifacts are SCRIPTS. Write your own re-decision scripts BEFORE reading the repaired code, from the
directive's Required sentences, the brief, the spec, and the committed records alone; save every script
and its literal stdout to named absolute paths under `/tmp/` and report those paths — a prose derivation
is worth nothing and will be discarded.

ABLATION IS MANDATORY. One class at a time, /tmp copy, real-record slice, literal verdict diff. A
corruption that moves nothing is a finding of the first order. Run one-sided corruption on any two-route
comparison the diff touches.

Probe for: a verdict verified using the object that produced it; an expected value living inside the
instrument or its calibration; an `assert` before the print it guards; a per-record line whose verdict
token does not follow from its printed operands; a calibration plant that cannot fail under the final code.

## Physics filter

Report a finding only if it catches a way a round-4 census verdict could be wrong, a real record class
could be silently unmeasured, or a certified-correct direction was weakened. Style, performance, and
hypothetical-input robustness are out.

## Sandbox

Copy anything you modify to /tmp and ablate the copy; never modify the working tree; never commit. No
Wolfram kernel — everything is Python + the record files. Cap any single script run at `timeout 600`; a
timeout is a failed probe to report, not a reason to retry harder. The full censuses take hours — run
slices, never the full population.

## Report format

Numbered findings, most severe first: claim, script path, literal stdout excerpt, instrument file:line
(at `3f229bbf`). If nothing survives the filter, list the ablations, samples, and literal outputs you ran
instead.
