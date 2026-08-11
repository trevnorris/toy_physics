# Build directive — close `N2` by removing the translation, not by mapping it

⭐ **Target:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.
⛔ **Edit ONLY that file.** ⛔ Do not commit. ⛔ Do not run either physics engine.

⚠ You reported last round that `N2` could not be written: `COEFFICIENT_ORDERING` is *"sorted by the
engine's own stated rule"*, so positions do not correspond, and every repair available was a mechanism you
would have had to invent. ⭐ **That report was correct and is why this round exists.**

⭐⭐ **The measurement below removes the problem instead of solving it.** ⛔ There is no map, no ordering
rule, and no translation step in what follows.

---

## ⭐ The measurement this round rests on

⭐ **What a transported point MUST carry** — taken from the free symbols of every **locus** `_EQUATIONS`
(`RANK_DROP`, `STACKED_DROP`, `ROOT_COINCIDENCE`, `KW_ZERO_LOCUS`, `STRATUM<s>_DEFINING`) and of `M`, in
the committed outputs of **both** engines:

| | measured |
|---|---|
| minimal transport set | **12 symbols**, each engine |
| wavevector components | **5**, and ⭐ **both engines already spell them identically** |
| coefficients | **7** |
| coefficients whose spellings differ | **5** |
| amplitudes `a_i` | ⛔ **NOT needed by any locus or by `M`, in either engine** |
| dimension-bookkeeping symbols | ⛔ **NOT needed** |
| coefficients **absent** from `§Q6r`'s existing table | ⛔ **none — all 7 are already on it** |

⚠ **The committed `STRATUM<s>_POINT` payloads nevertheless carry amplitudes and engine-local
dimension-slot symbols.** ⇒ ⭐ those engines emit more than a point is.

---

## ⭐⭐ P1 · A point's payload KEYS are the standard names — ⇒ nothing needs translating

⭐ **What must become true:** wherever this file requires an exact point — `STRATUM<s>_POINT`,
`§5`'s `_REAL_WITNESS`, and `§Q8c`'s witness input — the payload's **keys** are the **standard coefficient
names this file already uses**, together with the wavevector components. ⛔ Not an engine's internal
spelling.

⚠ ⭐ **This is `S9_REWRITE_PLAN.md` `D12` applied to a payload instead of a tag**, and `D12` is quoted
there as: *"Internal variables stay whatever WL needs (`rhoBr`, `muR`); only the string in the emit call
carries the standard name."* ⭐ **Read `D12` yourself and say whether it reaches a payload's keys.**
⛔ **If it does not, STOP AND REPORT** — ⛔ do not extend it on my say-so.

⭐ The standard names are the ones `§Q6r` already lists. ⛔ Do not author a second table, ⛔ do not add a
Wolfram column, ⛔ and do not require either engine to change an internal variable.

## ⭐ P2 · A point IS an assignment to the coefficients and the wavevector — ⛔ and nothing else

⭐ **What must become true:** the definition of an exact point says it assigns the package's coefficients
and the wavevector components, ⛔ and carries nothing else.
⚠ Apply it wherever a point is defined or required, ⛔ not only in `§Q8c` — ⚠ **last round you correctly
declined `P2`'s predecessor because narrowing it in one place would have manufactured a stale sentence
elsewhere.** ⭐ Find every place.

## ⭐ P3 · Retire the translation prose that `P1` makes unnecessary

⚠ Two places still instruct the orchestrator to re-spell the point using `§Q6r`'s map — ⭐ you identified
them as lines `183` and `822` before this round's edits. ⛔ Both legs measured that the map cannot do it.
⭐ Under `P1` **no re-spelling happens at all** ⇒ the instruction goes, ⛔ and is not replaced by another.

---

## ⛔ Constraints

- ⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ NEVER what anything equals or was measured.** ⚠⚠ **Every number
  in the table above is FOR YOU.** ⛔ No count, symbol list, spelling pair or measured fact may appear in
  the file in any form a builder could **derive** it from. ⭐ Test: **derivability**, ⛔ not presence.
  ⚠ In particular ⛔ do not write how many coefficients there are, ⛔ which spellings differ, or ⛔ that
  amplitudes were found unnecessary — ⭐ state the requirement, ⛔ never the finding behind it.
- ⛔⛔ **`§4`'s quoted block is shared VERBATIM with the other steps' specs — ⛔ do not edit inside it.**
  ⭐ Verify `sed -n '156,158p'` of the target stays byte-identical to `sed -n '111,113p'` of
  `directives/S10_SHARED_PHYSICS.md`, and report it.
- ⛔ Do not touch `_CANONICAL_LOCUS`. ⛔ Do not delete or weaken any `§5` locus object.
- ⛔ Do not re-open `N1` or `N3` — ⭐ they landed last round.
- ⛔ No admissibility algorithm, no recursive stratification, ⛔ no rule pinning how an engine **describes**
  a component.
- ⚠⚠ **AFTER EDITING, grep the whole file for sentences quoting a CONSEQUENCE of anything you changed.**
  ⭐⭐ **All five previous rounds left exactly one stale sentence, and every one was a downstream sentence
  asserting a consequence the edit falsified.** ⛔ Do not make the sixth.

## Deliverables

1. The edited file.
2. ⭐ Per `P1`–`P3`: lines changed and the sentence that now makes it true.
3. ⭐ Your reading of `D12` — ⭐ **does it reach a payload's keys, or not?** ⛔ Quote it.
4. ⭐ The `§4` byte-identity result, and every place you found that defines or requires a point.
5. ⛔ Anything you could not write without inventing, or any fact above that proved false ⇒ **STOP AND
   REPORT.** ⭐ Reporting is success.
