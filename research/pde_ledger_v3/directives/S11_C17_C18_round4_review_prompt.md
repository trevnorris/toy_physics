# Independent review — the four-item revision of `S11_SHARED_PHYSICS.md`

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree).
⭐ Diff against its base:
```
git -C /var/projects/toy_physics diff 4d3fd4ae -- research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
```

⚠ **Shared spec both engines read.** An error makes both engines agree on the same wrong thing, or
**disagree for a reason that is not physics**.

⚠⚠ **Five revisions; the first four bred thirteen defects between them, every one inside a mechanism its
author INVENTED to make two engines describe a locus the same way.** ⭐ This revision was restricted to
four items, each settled by a measurement supplied to the builder, with an explicit instruction to invent
nothing and to **report** anything that could not be written without inventing.

⭐ **It reported two such things. Adjudicating those is the most important part of this review** (§B).

You are one of two independent legs; the other is not visible to you.

## ⛔ DO NOT READ
`directives/S11_C17_C18_round4_directive.md` (the directive this was built from), any `…_fix_round1_*`,
`…_round2_review_prompt`, `…_round3_review_prompt`; anything under `_scratch/`.

## ⭐⭐ READING ORDER — ⛔ artifact last
1. `DEFECT_REGISTER.md`, `C17` and `C18` in full.
2. `directives/S11_C17_C18_spec_repair_decisions_v2.md` — decisions `T1`–`T7`.
3. The committed outputs `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` and
   `scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⭐ **Open them.**
4. **Write down** how *you* would let two blind engines compare a stratum's spectrum when each may describe
   the component in different variables. Keep it.
5. **Only now** the artifact and the diff.

## §A · The four items — ⛔ judge each discharged or not, quoting the sentence

| | what had to become true |
|---|---|
| `M1` | the witness point-coverage token is decided over **everything that slot evaluates** — every emitted locus's equations **and** `M` — ⛔ not `M`'s symbols alone |
| `M2` | a component-scoped **symbolic** payload is a comparison row again, compared as **the difference reduced modulo the defining equations** — ⛔ and this must NOT become a rule pinning how either engine describes the component |
| `M3` | the orchestrator delivers the point in the **receiver's own coordinate names**, naming an existing map — ⛔ authoring no second map |
| `M4` | one stale sentence in `§4` that still counted **two** input fields |

## §B · ⭐⭐⭐ THE TWO REPORTED GAPS — ⛔ adjudicate these, they are the point of this review

The builder wrote the rule it was told to write and then reported that **neither is fully determined**.
⭐ **Settle each by opening the committed outputs, ⛔ not by argument.**

### B1 · *"reduced modulo `STRATUM<s>_DEFINING_EQUATIONS`" — ⛔ whose?*

⚠ Each engine emits **its own** defining equations, in **its own** chart. ⭐ The text does not say whether
the reduction is taken modulo one engine's set, the other's, or both.
⛔ **Does that ambiguity change the answer?** ⭐ Take the real committed component from the outputs, form
the two engines' symbolic payloads, and reduce the difference **each way**. ⭐ Report the literal results.
⚠ If all readings agree on the measured case, say whether that is luck or structural — ⭐ and construct a
case where they disagree, or state that you could not.

### B2 · `§Q6r`'s map may not be able to do the job `M3` gives it

The builder reports: `§Q6r`'s map is *this file's notation → `LEDGER` key names*, engine-local to the
**importing** engine; its right-hand side is one ASCII spelling; it carries **no** Wolfram spelling, and
Wolfram symbols cannot contain `_`. It also covers the coefficients only — committed `STRATUM<s>_POINT`
payloads additionally carry wavevector and amplitude components and at least one engine-local symbol.
⭐ **Verify all of that against the file and the outputs.**
⛔ Then: can the map deliver a point in a **Wolfram** receiver's own names? ⚠ The builder notes the `k` and
`a` spellings happen to coincide in the as-built pair, so the gap is **latent**.
⇒ ⛔ **Is `M3` actually discharged, or does the spec now instruct something unperformable?** ⭐ Say which,
and what the smallest honest fix is — ⛔ do not design a second map.

## §C · What else to check

1. ⛔⛔ **Does the file state or imply an answer?** ⭐ Clause by clause on added lines. ⚠ The builder was
   given **literal measured values** — a factor-of-two dispersion error, specific counts, a symbolic
   residual. ⛔ **None may be derivable from the file.** ⭐ The test is derivability, ⛔ not presence.
2. ⭐⭐ **Did `M2` become a normalisation?** ⛔ Does anything now pin which variable an engine eliminates,
   or oblige an engine to perform the reduction itself? ⚠ A previous round pinned an elimination and it
   **deleted a branch**. ⭐ Quote the guards and say whether they hold.
3. ⭐ **Did the revision touch what it was forbidden to?** ⛔ `§4`'s quoted block must be byte-identical to
   `directives/S10_SHARED_PHYSICS.md:111-113`; `_CANONICAL_LOCUS` must be untouched; no `§5` locus object
   deleted or weakened. ⭐ Check each by command.
4. ⚠ **Three additional edits the builder made on its own initiative** — an `N6_BASIS` carve-out, a `§9`
   table row, and a note that the name map is not a pairing condition. ⭐ Are they entailed by `M1`–`M3`,
   or are they scope creep? ⛔ Does the `N6_BASIS` carve-out weaken `M2`?
5. ⚠ **Stale prose. ⛔ ALL FOUR previous rounds left exactly one stale sentence.** ⭐ Grep the whole file
   for sentences referring to anything this revision changed. ⛔ Find the fifth, or state that there is
   none.
6. ⭐ What of `T1`–`T7` is still not discharged?

## Physics filter
Report a finding only if it catches a way the **physics could be wrong**, a **claim unsupported**, or
**wrong or incomparable artifacts**. ⛔ Not style.

## Method
- ⭐⭐ Quote both sides of every finding.
- ⭐ State what you opened or ran for every claim about a file, line or count.
- ⭐⭐ `B1`, `B2` and check 3 are settled by **running it** on real committed data. ⛔ Prose is discarded.
- ⭐ Any CAS run: `timeout 600`, ⛔ one kernel at a time (two seats), ⛔ copy to `/tmp` and run the copy.
- ⛔ Read-only.
- ⭐ End with: `M1`–`M4` discharged or not; a verdict on `B1` and `B2`; whether the file leaks; whether a
  fifth stale sentence exists; and what you checked that could have failed and did not.
