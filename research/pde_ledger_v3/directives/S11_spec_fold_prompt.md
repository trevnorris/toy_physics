# Fold the review findings into the repaired S11 specification

## What you are editing

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **999 lines**,
uncommitted in the working tree. ⛔ The **only** file you edit.

⚠ You applied the five-item repair that produced this file. It went to two independent review legs. ⭐ The
five fixes below are what they found, plus one you flagged yourself. ⛔ Apply these five and **nothing
else** — everything not named stays byte-identical.

⛔⛔ **THIS FILE MAY NEVER STATE WHAT ANYTHING COMES OUT TO BE** — no value, count, sign, rank, dimension,
or expected effect of any package. ⛔ No justification from which one could be derived.

---

## 1 · `§2:55-56` misdirects the coefficient inventory

It reads: *"§6 says to read them off that package's own `W`."* ⛔ After the repair that is **false** —
`Q6:439-451` builds the inventory from the declared additive terms of `L_pkg`, together with every
inertial coefficient the kinetic density carries.

⭐ Rewrite that sentence so it points at what `§6` now does. ⚠ Keep the surrounding warning intact: a
package's `W` is a functional, and ⛔ nothing may assume which stiffness coefficients it contains or how
many. ⭐ The correction is to the **inventory domain**, ⛔ not to that warning.

## 2 · `STRATUM<s>_SOURCES` has no pinned token vocabulary

`:611-612` describes the payload in unbackticked prose. ⚠ Every other enumerated payload in this file
pins its tokens (`:322`, `:409-410`, `:491`, `:296`, `:529`, `:723`) — this is the sole exception, it is a
**cross-engine compared row**, and unpinned spellings are recorded at `:413-415` as having produced
eleven disagreements that were two names for one determination.

⭐ Pin the payload to a list drawn from exactly these three tokens, spelled and backticked this way:

```
RANK_DROP · STACKED_DROP · ROOT_COINCIDENCE
```

⭐ Pin the collection shape: the tokens for that component, **in the order listed above**, each appearing
at most once. ⚠ `Q8b:604-605` calls one of these *"stacked degeneracy"* while `Q8a:588` emits
`ROOT<r>_STACKED_DROP_MINORS`; ⭐ make the prose use the pinned token so one source has one name.

## 3 · `Q8b` does not say which locus families are candidates

`Q8a` now emits **six** families per root — `RANK_DROP_{K,COEFF,JOINT}` (`:582-584`) and
`STACKED_DROP_{K,COEFF,JOINT}` (`:593-595`) — and `Q3` emits two coincidence families (`:335-342`), eight
in all. `:604-606` names three *sources* but never says which families feed the candidate pool. ⛔ Two
engines taking different subsets produce different strata and different `Q3`/`Q4` reruns.

⭐ State that **all eight families are candidates**: all three `RANK_DROP` solve-variable sets, all three
`STACKED_DROP` sets, and both `Q3` root-coincidence families.
⛔ Do **not** add a deduplication rule, and ⛔ do not forbid one — `:607-608`'s existing silence stands.

## 4 · Two passages enumerate split-breaking by stiffness only

`Q4:371-373` says *"§7 contains a package whose **stiffness functional** is not built from `S_curl` and
`S_div` alone"*, and `§7:864-866` says `§Q4`'s warning *"is aimed at this package"* — naming
`W_XFORM_EXTRA` alone. ⚠ A package may now break that split through its **kinetic** density instead.

⭐ Generalise both so the concern is scoped by **what a package's action does**, ⛔ not by which member of
the pair carries it, and so `§7:866` no longer reads as scoping the concern to one package.
⛔ Do **not** state how any package's null directions lie relative to `k`, and ⛔ do not name which
packages are affected — `Q4:374-375`'s *"compute `N3` everywhere, unconditionally"* is the operative
instruction and stays exactly as it is.

## 5 · `Q6r`'s no-row rule collides with Corollary 4

You flagged this: `:521-522` emits no row for a coefficient the `LEDGER` does not carry, while
Corollary 4 (`:176-181`) forbids emission conditional on anything but package, dimension and quantity, so
that *checked-and-absent* stays distinguishable from *never checked*.

⭐ Resolve it by making the absence a **computed object** rather than a silence: `Q6r` emits, once per
`(package, D)`, the coefficients of `COEFFICIENT_ORDERING` that **resolved** to a dimension row and those
that did **not**, both read from the live `LEDGER` lookup. ⭐ Give the two objects tag names in this file's
style. ⭐ The per-coefficient comparison rows stay as they are — only resolved coefficients have vectors to
compare. ⛔ Do not manufacture a placeholder vector for an unresolved coefficient.

---

## Report back — ⛔ under 20 lines

1. The new line count, and one line per fix saying where you applied it.
2. ⭐ Anything above you could not apply without choosing something it did not specify — **name the
   choice.**
3. ⭐ Any place a fix made an existing passage wrong or inconsistent.
4. ⛔ Do not report what any computation will produce, and ⛔ do not evaluate whether a fix is a good idea.
