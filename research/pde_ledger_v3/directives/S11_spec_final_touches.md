# Four one-line fixes to `S11_SHARED_PHYSICS.md` — ⛔ nothing else

**Do not commit.** One file: `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.

⭐⭐ **THE GOVERNING INSTRUCTION, unchanged from the last round: prefer REMOVING or LOOSENING to adding.**
⛔ Do not attempt to make this file complete, closed, or self-consistent in every clause. Two review rounds
have shown that each attempt to make it airtight breeds the next round's defects. ⛔ Add no new tag, rule,
guard or check. ⛔ Improve nothing you were not asked to touch. ⭐ All four items below are **one-line**
edits; if any of your diffs is longer than a line or two, you have over-reached.

---

## F1 · ⛔ BLOCKING — `§Q11`'s sign test contradicts the token set it binds

Around `:670`, `ROOT<r>_KW_SIGN` is specified as *"the **three-way** symbolic test on `k_w²` … token set as
in Q3"*. `§Q3` (`:294-296`) specifies **four** tokens: `POSITIVE`, `NEGATIVE`, `ZERO`, `UNDECIDED`.

⚠ **This is not decorative, and both review legs found it.** §3 supplies no relation between `c_s0` and the
stiffness coefficients, so the fourth token is **reachable for this object**. An engine implementing a
literal three-branch test must either force a determination or emit a token outside the pinned set — and
two engines resolving that differently is exactly the false disagreement the pinned token set exists to
prevent. ⇒ The file's own `§Q5` note records eleven such rows at the previous step.

⭐ Fix the wording so `§Q11` matches `§Q3`. ⛔ Do not change the token list, and ⛔ do not add a token.

## F2 · REMOVE — a deleted claim surviving in the header

`:8` still says the file supplies *"the name each one is emitted under"*. That is the closed-vocabulary
claim, which was deleted from `§8` last round, and it now contradicts `:827` (*"Where you emit an object
for which this file supplies no tag name, choose one and list it in the §10 report"*).

⭐ Delete the clause from the header. ⛔ Do not add a replacement.

## F3 · FIX — two stranded pointers to `§8`'s deleted vocabulary

`:186` and `:706` refer to *"§8's vocabulary"* / *"the FULL §8 tag vocabulary"*, which `§8` no longer
defines. The operative content survives in the adjacent sentences in both places, so ⭐ **reword each
pointer to refer to what §8 now says** — one tag per named object, at the scope §8 gives it. ⛔ Do not
restore any vocabulary claim, and ⛔ do not enumerate anything.

## F4 · FIX — `§Q8a`'s rank-drop tags are per-root but carry no root scope

Around `:527-534`, `RANK_DROP_MINORS` and `RANK_DROP_{K,COEFF,JOINT}_*` are computed *"for each `M_r`"* —
i.e. once per root — but carry no `ROOT<r>` scope, so per-root emission collides on a single tag name.
⚠ This is the **same collision** the previous round fixed for `KW_ZERO_LOCUS`.

⭐ Give them `ROOT<r>` scope, exactly as `ROOT<r>_KW_ZERO_LOCUS_*` now has. ⛔ Add no tag; rename only.

---

## ⛔ Do not touch

Everything else. In particular: the three clauses and five corollaries · the no-verdict rule · `§5`'s
five-object locus protocol · the pinned sign and marker token **sets** · the pinned monomial ordering and
the RREF requirement · `§Q7`'s five tags · `§Q11`'s bulk field, ansatz and its ban on the
back-substitution residual · the package list and the `D` sweep · `§Q6d`'s vacuity explanation · `§7`'s
runtime rule.

⚠ `§Q6d`'s explanation has twice been flagged as a leak and that finding was **considered and rejected**
both times. ⛔ Do not remove it.

## Acceptance — paste literal output

1. `git diff --stat` if the file is tracked, otherwise the net line change, and **each of your four diffs
   quoted in full**.
2. `grep -n "three-way\|four-way"` — the counts and the lines.
3. `grep -n "the name each one is emitted under"` — expect no hits.
4. `grep -n "RANK_DROP"` — every hit, showing the scope each now carries.
5. `grep -n "§8's vocabulary\|FULL §8 tag vocabulary"` — expect no hits.

## Constraints

- ⛔ No `git add`, no `git commit`, no other git write.
- ⛔ Touch no other file. ⛔ Create no new file. ⛔ Do not run either engine; neither exists yet.

## Report back — under 12 lines

1. `F1`–`F4`: done or not, with the net line change.
2. The acceptance output.
3. ⭐ Anything in `F1`–`F4` you judge wrong. ⛔ Do not list other defects you noticed — ⭐ say only how
   many you saw and leave them.
