# Close-out review — are `C17` and `C18` actually closed?

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree).
⭐ The spec **before any of this repair** is `git show 767828f4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.

⚠ **Shared spec both engines read.** An error makes both engines agree on the same wrong thing, or
**disagree for a reason that is not physics**.

⚠⚠ **This repair took EIGHT rounds and bred sixteen defects along the way.** ⛔ Do not assume the survivors
are sound. ⭐ **The question is not "did the last round work" — it is "are `C17` and `C18` CLOSED".**

You are one of two independent legs; the other is not visible to you.

## ⛔ DO NOT READ
Any `S11_C17_C18_*directive*` or `*review_prompt*` file; anything under `_scratch/`.

## ⭐⭐ READING ORDER — ⛔ artifact last
1. `DEFECT_REGISTER.md`, `C17` and `C18` **in full**, including their measured divergences.
2. `directives/S11_C17_C18_spec_repair_decisions_v2.md` — `T1`–`T7`.
3. The committed outputs `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` and
   `scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⭐ **Open them.**
4. **Write down** how *you* would close `C17` and `C18`. Keep it.
5. **Only now** the working-tree file.

## ⭐⭐⭐ THE TWO QUESTIONS THAT DECIDE THIS

### 1 · Is `C17` closed?
⚠ `C17`: a stratum's `Q3`/`Q4` were evaluated **at one point** and emitted under the **component's** scope,
so two builders could pick different valid points on one component and emit different physics.
⭐ **Trace the register's own witness through the CURRENT text**, on real committed data. ⛔ Is there still
any place where a mode count moves and nothing says so?

### 2 · Is `C18` closed AS FAR AS IT CLAIMS?
⚠ `C18`: one engine emitted a stratum at `XFORM_EXTRA, D = 2` and the other emitted none, both faithful to
the words. ⚠⚠ **The spec claims only PARTIAL closure**: the divergence is now typed, and `§9` states that
**no computation resolves it**.
⭐ **Trace that exact case through the CURRENT text.** ⛔ What does each engine emit, and what can a
comparator conclude? ⭐ Is the divergence genuinely **visible and typed**, or can it still present as
*"a list against an empty list"*? ⛔ **And is `§9`'s disclosure honest and sufficient** — could a reader of
the record mistake this for something the build established?

## ⭐⭐ What to attack

1. ⛔⛔ **Leak.** ⭐ Clause by clause over the diff. ⚠ The builders were given measured values — root
   counts, a factor-of-two dispersion error, which spellings differ, that amplitudes proved unnecessary.
   ⛔ **None may be derivable from the file.** ⭐ Test: derivability, ⛔ not presence.
2. ⭐⭐ **Buildable by two engines that never speak?** ⛔ Does any engine receive anything at all from the
   other? ⭐ Is the Wolfram engine still importing nothing? ⚠ A cross-engine exchange was **removed** this
   round — ⛔ check nothing survives that still assumes one.
3. ⛔⛔ **THE REMOVAL.** ⭐ A large section was cut. ⛔ **Grep for every dangling reference and every
   sentence asserting a CONSEQUENCE of something no longer in the file.** ⚠ A removal falsifies downstream
   prose exactly as an edit does, and this is the round most likely to have left several.
4. ⭐ **The component comparison.** A component-scoped symbolic payload is compared as the difference
   reduced modulo **each engine's own** defining equations, both carried as operands. ⛔ Does anything pin
   which variable an engine eliminates? ⚠ A previous round pinned it and **deleted a branch**.
5. ⭐ **Counts vs symbolic.** ⛔ Verify on real committed data that the counts are invariant under which
   variable is eliminated, and that the reduction catches a wrong dispersion coefficient the counts miss.
6. ⚠ **Stale prose.** ⭐ Grep for sentences quoting a **consequence** of anything the repair changed.
   ⛔ **Seven of the eight rounds left exactly one.**
7. ⭐ **What of `T1`–`T7` is still undischarged?**

## Physics filter
Report a finding only if it catches a way the **physics could be wrong**, a **claim unsupported**, or
**wrong or incomparable artifacts**. ⛔ Not style.

## Method
- ⭐⭐ Quote both sides of every finding.
- ⭐ State what you opened or ran for every claim about a file, line or count.
- ⭐⭐ Questions 1 and 2 and attacks 3 and 5 are settled by **running it** on committed data. ⛔ Prose is
  discarded there.
- ⭐ Any CAS run: `timeout 600`, ⛔ one kernel at a time (two seats), ⛔ copy to `/tmp`.
- ⛔ Read-only.
- ⭐ End with: **is `C17` closed? is `C18` closed?**, `T1`–`T7` status, whether the file leaks, whether a
  stale sentence survives, and what you checked that could have failed and did not.
