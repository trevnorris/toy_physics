# Independent physics review — the S10 DIRECTIVE PACKET, before any build runs

You are reviewing a **specification**, not a result. Nothing has been built yet. ⭐ **This is the one
artifact both engines will share, so an error here makes dual-engine agreement certify wrong physics.**
That is why this review happens first.

## Artifact — read all three, in this order

1. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
2. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_wl_directive.md`
3. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_py_directive.md`

## ⛔⛔ DO NOT READ — these contain the answers this build exists to produce

- `/var/projects/toy_physics/research/pde_ledger_v3/steps/` — **the entire directory**, every file. ⛔ Do
  not list it, do not grep it, do not open any file in it, and ⛔ do not read it via `git show` or
  `git log`.
- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`
  ⚠ These two are the **previous** engines, being replaced. Reading them anchors you to the framing under
  test.
- `/var/projects/toy_physics/research/pde_ledger_v3/paper/` — the TeX cards and the built PDF.

⭐ **You MAY read** `REBUILD_HANDOFF.md`, `CHARTER.md`, the `reduction/` registry, and the skills under
`.claude/skills/`.

## ⛔⛔ REQUIRED METHOD — A PROSE DERIVATION IS WORTH NOTHING

> *"trusting grok and codex and even yourself on how they 'rederive' is not trustworthy. Unless it's in CAS
> and we can see the output from the inputs, it's not to be trusted."*

⇒ ⭐ **Every derivation claim you make must be backed by a script you wrote and its literal stdout.** Save
both to named absolute paths **under `/tmp/`, ⛔ never inside the repository**, and give me the paths.
⛔ **A derivation claim with no script and no output will be discarded.**

⚠ Two hand derivations that agree are two unverifiable claims that happen to coincide — which is exactly
how two cancelling errors survive.

---

## What to check

### ⭐⭐⭐ A · Is `N3` really the transverse mode count? — THE LOAD-BEARING METHODOLOGICAL CLAIM

Shared physics §6 Q4 requires the transverse count to be computed as
**`nu_T = D − rank( M_r stacked on top of the single row kᵀ )`**, and it **forbids** obtaining that count
by classifying the vectors a `NullSpace` call returns.

⭐ **Settle this by computation, and it is the single most valuable thing you can do here:**

1. Prove or refute that `D − rank([M_r ; kᵀ])` equals `dim( null(M_r) ∩ k^⊥ )` for an arbitrary square
   `M_r` and nonzero `k`. ⚠ Is symmetry of `M_r` required? Is anything required?
2. ⭐⭐ **Construct an explicit counterexample matrix** on which the forbidden method (classify each
   `NullSpace` basis vector as parallel / perpendicular / neither) gives a **different or unusable** answer
   than `N3`, and **show the literal output of both methods** on it. If you cannot construct one, say so —
   ⛔ that is also a finding, and it would mean the directive is defending against nothing.
3. Is `N3` **sufficient**? Does the packet request everything needed to characterise each root's null
   space, or is there a case `N1`–`N7` cannot distinguish?

### B · Can every requested object actually be computed from what is supplied?

Work through §6 Q1–Q8 against §1–§3. ⛔ For each requested object, is it **computable from the supplied
action and ansatz alone**? Name anything that is not, and anything that is **ill-posed**, **ambiguous**, or
**tautological**.

⭐ **Specifically interrogate Q2's two-route claim.** Route A takes the Euler–Lagrange equations in
coordinate space and substitutes the plane wave; route B takes the second derivative of the plane-wave
quadratic form. ⛔ **Are these genuinely independent, or does one reduce to the other by construction?**
If the residual `M_A − M_B` is zero **for any input whatsoever**, it is operand theatre and the packet
should say so instead of dressing it as a check. ⭐ Answer this with a computation.

### C · Is the physics specification correct and complete?

- Is the Euler–Lagrange equation as written in Q1 correct, including signs?
- Is the stiffness density `S_curl` in §2 well defined in every `D`, and is the `1/2` consistent with
  summing `i` and `j` each over the full range?
- Is the ansatz in §3 consistent with taking a **real** quadratic form — does using a complex exponential
  while treating `a` as real introduce an error the packet does not address?
- ⚠ Is `omegaSquared` as a single symbol safe everywhere it is used?
- Does anything in the packet **assume** what it asks to be computed?

### D · ⚠ The static-or-instantaneous check — run it yourself

> **What timescale or rate does this specification send to `0` or `∞`, and what would it have governed?**

§3 declares `v₀ = 0`, no dissipation, and exactly-quadratic `L`. ⛔ **Are those the right premises, is the
list complete, and does any of them answer a question the step is trying to settle?** ⚠ A limit that
removes the very rate governing the quantity under study answers by assumption, and a closed-form result
looks equally healthy either way.

### E · Do the two engine directives specify the SAME physics?

⛔ Any asymmetry other than the deliberate one — engine 2 reads the registry, engine 1 does not — makes one
engine better-informed and quietly destroys the value of their agreement. Name every asymmetry you find.

### F · Does the packet leak an expected value?

⭐ The packet is supposed to supply **every equation** and withhold **only** an acceptance criterion
referencing an expected value. ⛔ Report anything that states, implies, or telegraphs a result: a count, a
root, a sign, a ratio, a dimension, or the *shape* of an answer — ⚠ **including inside a sentence that
forbids it**, which is a failure mode this project has measured repeatedly.

### G · ⚠ Runtime feasibility — a specific concern I want settled

§7 requires `XFORM_ANISO` (kinetic term `½ Σ_j ρ_j (∂_t u_j)²` with **distinct** `ρ_j`) at `D = 3` and
`D = 4`, and `MAIN` up to `D = 6`. ⛔ **Is solving `Det[M] = 0` for `omegaSquared` tractable there**, or
does it produce unsolvable radicals / `RootObject`s that make the downstream rank work meaningless?
⭐ **Test it** — actually run the determinant and solve for the anisotropic case at `D = 3` and `D = 4` and
report wall-clock and the form of the result. If it is intractable, ⭐ propose the **smallest** change to
the control that still breaks the isotropy the control exists to break.

⚠ No script in this project may exceed **10 minutes**. ⛔ Raising a timeout is never the fix.

---

## Physics filter

⛔ Report a finding only if it catches a way the **physics** could be wrong, or a way the **specification**
could make an engine compute the wrong thing. ⛔ Do not report style, and do not report "the script would
be wrong on a different input".

## Report format — ⛔ under 40 lines

1. **Findings, most severe first.** For each: what is wrong, where (file and line), why it matters
   physically, and the **script path + literal output** that establishes it.
2. Your explicit answer to **A.2** — the counterexample, or a statement that none exists.
3. Your explicit answer to **B**'s two-route question and to **G**'s tractability question.
4. ⭐ Anything wrong that this prompt did not ask about. ⭐ This is wanted.
