# Techniques that earned their keep — 2026-07-30/31

⭐ Each of these **caught something real this session**, with the catch named. ⛔ This is not a process
manual; it is a list of moves that paid. If one stops paying, drop it — the governing test is
`[[feedback-physics-not-ceremony]]`.

---

## 1. ⭐⭐ Triage before editing — "a mention is not damage"

Sort every hit into damaged/benign **before touching anything**, against explicit criteria.

**Paid:** a corpus-wide grep for `a` collapsed to **17 loci**. Found `STATUS.md` **already repaired** and
`stage043`'s engines **clean** — two files that would have been "fixed" pointlessly. Also found the
decision doc's own tier-1 table was **stale**.

⇒ Editing on a grep hit is how a symbol-level change becomes a corpus-wide mess.

## 2. ⭐⭐ Take the reviewer's CLASS, verify the FACT yourself

External reviewers are reliable about *kinds* of error and unreliable about *specific loci*.

**Paid three times:**
- Grok cited `parameter_register:295` for the second healing length. It is at **`:294`** — and opening it
  found something *sharper* than reported (the register already flags it as a convention, and it explains
  why the pin survived).
- I doubted a sub-agent's finding because my own grep window truncated the evidence, and had to reverse —
  the agent was right about both class *and* loci.
- Grok's `:185` for a token that is actually near `:189`.

⇒ Take the class. **Open the file for the fact.** Both directions: reviewers are wrong about loci, and
*you* are wrong when your grep is lazy.

## 3. ⭐⭐ The tripwire — leave a stale WRONG prediction where a cribber would find it

The committed plan said `dim_before 10, dim_after 6, Δ 4`. The correct answer after the actual edit was
`9/5/4`.

**Paid immediately:** Codex self-disclosed that a mis-rooted `--exclude` had exposed the forbidden plan
file. Without the tripwire that disclosure would have been **unfalsifiable** and the whole leg would have
had to be re-run. Because it reported `9/5/4` and not `10/6/4`, the exposure was *provably* harmless.

⇒ ⭐ Build the tripwire **before** you need it. It converts "I have to trust the builder" into "I can
check."

## 4. ⛔ "Do not read X" does not survive a grep — exclude structurally

**Paid:** GLM began grepping the reviews directory it had been told to avoid; Codex hit the plan file via
a mis-rooted exclusion. Neither was disobedient — `rg` has no notion of a forbidden path.

⇒ **Move the contaminating artifact out of the repository** before launching. Keep the prose prohibition
too — it costs nothing — but never rely on it. Second round, after moving the reviews to the scratchpad:
contamination check came back **zero**.

## 5. ⭐ Record your own derivation BEFORE reading anyone else's

Wrote my predicted payload to a file, timestamped, then read the reviews.

**Paid:** when GLM, Grok and Codex all returned `9/5/4`, I knew I hadn't anchored on them — and could say
so honestly rather than hoping.

## 6. ⭐⭐ Agreement is worth something only when the METHODS differ

Four parties derived the post-edit numbers: me (from the two switch semantics), GLM (SymPy Gröbner on a
reconstructed registry), Grok (the live `Registry` object in memory), Codex (the edited tree). **Four
routes, one answer.**

⇒ ⛔ Four parties running the same script is one derivation with three witnesses. Vary the *method*, or
the agreement is decorative.

## 7. ⭐ Run the acceptance commands before handing them to a builder

Ran all four gates on the pre-edit tree: green, with recorded output.

**Paid:** the directive's acceptance block was a command set I had *executed*, not one I hoped worked —
and the recorded baseline is what made "did the edit take effect?" answerable.

## 8. ⭐⭐ The Δ column is the tell

`baseline.dim_after` was **5 before and 5 after**. That column alone reads as *"nothing happened."*
**Δ went 5 → 4** — one fewer genuine constraint.

⇒ Predicted this in advance and used it as the acceptance criterion. ⛔ When a bookkeeping change and a
physics change would look identical in the headline number, find the column where they differ **before**
you look at the result.

## 9. ⭐ Ask "is this the same KIND of thing?", not "is this value right?"

The pin's value `ħ/(m c_s0)` **is** a healing length in a standard convention — correct number, correct
dimensions, plausible provenance. A value-check passes it.

⇒ The kind-check kills it in one line: `ξ_h` is **one number for the whole medium**; throat radii are
**one per particle**. ⛔ Run this on every identification of two same-dimension quantities.

## 10. ⭐ Scope rule: fix what computes, quarantine what only narrates

Cut step ① from ~21 files to **7**, with the omissions named explicitly in a do-not-touch list so the
builder wouldn't "helpfully" widen.

⇒ Repairing prose in a ledger you are superseding catches nothing about the physics. Repairing a registry
the next ledger will **compute with** catches plenty. ⚠ The honest cost — known-false statements left in
v2 — is acceptable **only** because they are recorded in the defect register and v3 derives its own count.

## 11. ⭐ Review the DIRECTIVE before it runs, not just the output

Two independent legs read the directive pre-execution. They found **five** defects, including two
**self-contradictions of mine** — a banner instruction that violated my own line-preservation rule, and a
"retire the probe" option that was a guaranteed dead end via a hardcoded count in an out-of-scope test.

⇒ Each would have cost a full build cycle. ⛔ The gate must **precede** launch — a builder snapshots the
prompt, so editing it afterwards does nothing.

## 12. ⭐ Ask reviewers to derive the answer, not just critique the question

The directive review prompt asked for an **independent derivation** of the post-edit numbers.

⇒ Turned two review legs into two extra independent references, for free. A reviewer who only critiques
gives you an opinion; one who derives gives you a *check*.

## 13. ⛔ Long jobs need `setsid` — the harness reaps background jobs

First review round: both legs **killed** mid-run (one at 606 bytes). Relaunched with
`setsid nohup … > file 2>&1 < /dev/null &` plus a `___DONE___` sentinel and an `until grep -q` watcher.
Second round completed cleanly.

⛔ Never wrap in shell `timeout` — SIGKILL wastes the whole run.

## 14. ⭐ A step whose result is "not determined" is a real step

`W_slab` selects no width. `m_defect` has no route. `g_phys` was never mapped.

⇒ Bank these as **completed steps with a negative result**, not as failures to retry. ⛔ Grinding on them
is how the interior swallowed a year.

## 15. ⭐ Verify the load-bearing claim yourself, even from a trusted source

Checked by hand: the three-pin determinant; the entailment `2m h0 ξ_h² − ħ² = 0` on the constraint
variety **and** its ideal combination; `INV2`'s `a`-cancellation; the `stage004` anchors at `:65`,
`:76`, `:102` after edit; the `test_registry`↔`CASES` coupling; `_validate_loci` being bounds-only.

**Paid:** `_validate_loci` turned out to check only that a range *fits inside the file* — never that the
cited lines still say what is cited. ⇒ **A line shift passes every gate silently.** That is the mechanism
behind the corpus-wide citation rot, and nobody had reported it.
