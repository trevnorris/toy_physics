# S11b Wolfram engine — build directive

## Authority and boundary

Write `research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl` in full. Its
product is the flushed stdout tag stream; that is its only write. The two historical S11b `.wl` engines
(`mathematica/S11bA_*`, `mathematica/S11bB_*`) and their stdout are **not** build premises — this is one
unified step, not a merge of those two.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` is the sole physics authority
and the sole physics input: every equation, ansatz, premise, supplied law, derivation route and tag
obligation enters from it alone, and it wins every conflict with this directive. Implement its §§0–13
directly. Point to rather than duplicate §4; §5 and its basis construction; §6, its virtual-displacement
rule, causality diagnostics, convention cross-checks and energy accounting; §7; §8, including its
corollaries, the no-verdict rule and the locus protocol; §9; and §10. Add no expected value or acceptance
criterion (`CLAUDE.md` rule 5).

⛔⛔ **This is the blind engine.** The build's inputs besides the spec carry no S11b physics: this directive,
`CLAUDE.md`, and the mechanical precedent in the next section. **The engine this build writes reads no file
at run time — it imports nothing and re-derives from the spec's equations alone.**
`S9_export_chain_rebuild_directive.md:16-18`: **this independence is the design's one blindness control, and
nothing else may be built pretending to be one.** ⛔ This engine imports no `LEDGER` and writes none; the
§11 inherited-input lookup is the **importing engine's** and has no counterpart here — re-derive every
dimension from §1–§6 and emit **no** import tag (that absence is by construction, `S11b_unified_decisions.md`
G7, ⛔ not a denylist). ⛔⛔ **There is no do-not-read list, and none may be added** (rule 12, G9): blindness
is enforced by the engine importing nothing, ⛔ never by a sentence forbidding a read.

Do not add parallel machinery for the spec's contracts: no per-cell completion registry beyond §9's task
structure, and no directive-local exit policy. §8 corollary 4 applies as a property, without named
exceptions. ⛔ **No `VERDICT`, `PASS` or `FAIL`** (§8): a physics finding exits 0, and a boolean-valued test
is emitted **as the CAS object the test returned**, ⛔ never a host-language native boolean. This is a change
from the two historical S11b stages, which ended in a `VERDICT` line.

## Mechanical precedent

`research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` is **mechanical precedent,
not authority.** It shows the working shape under this design: an internal tag mapped to a standard
emitted-name string at the emit site, a duplicate emitted name a hard stop, every payload printed as
`NAME: value` in `InputForm` as computed, and flushing established on a redirected stream. Its physics is
another step's. §10 standardises emitted-name strings only (the `WL_S11B_<QUANTITY>` grammar); symbol
spellings inside the engine are its own. ⛔ Do not import or transcribe any object from it.

## Run discipline — binds the build's own demonstration runs

- The full derivation, and with it the writing of
  `mathematica/out/S11b_interface_coupling_law_mathematica_audit.out`, belongs to the **orchestrator after
  review**. The build's demonstration runs go under the build's own scratch paths, ⛔ never under
  `mathematica/out/`.
- The licence has two seats: run **at most one kernel** at any time.
- Kill criteria for a demonstration kernel, both by PID, both recorded with the task they interrupted and
  reported, ⛔ never answered by narrowing or cheapening the requested object: **600 seconds with no new
  output** (the runtime rule is observable progress, not elapsed time — a visibly streaming task may run
  long), or resident memory past **6 GB**.

## Acceptance — executable, no expected values

1. **Blindness / no-hidden-dependency and no-stray-write.** Copy the finished `.wl` alone into an empty
   scratch directory and run it there and in-repo: both runs exit 0, the streams are non-empty and
   byte-identical, and after each run its scratch directory contains exactly the one file that was copied
   in. Any run-time dependence on a repository file, and any write other than stdout, breaks this — and
   would break the blindness control.
2. **Flush under redirected stdout.** Run with stdout redirected to a file; on at least one task the capture
   file must be observed **growing while the kernel is alive**. `Print` alone does not establish flushing on
   a redirected stream; satisfying that under redirection is this build's to demonstrate. If every task
   completes too quickly to observe growth, record that in the report — ⛔ do not claim the flush obligation
   was demonstrated.
3. After all demonstrations: `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints
   nothing.

## Conflicts

Use the spec's §10 (and the §13 report) for any conflict, ambiguity, unavailable construction or object you
cannot emit under one-tag-per-object. ⛔ Do not fill such a gap with new physics, an expected result, or a
self-review mechanism. A refusal — `NOT_ESTABLISHED` with what is missing — is a valid output (§0).
