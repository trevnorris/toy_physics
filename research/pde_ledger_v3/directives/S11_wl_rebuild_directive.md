# S11 — build the Mathematica engine (engine 1, BLIND)

**Write:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
**Run it into:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out`

**Read first, in full, and treat as binding — it is the ONLY artifact you read for physics:**
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`

⛔ **Do not commit.** ⛔ **Do not modify any file** other than the two named above.

---

## ⛔⛔ WHAT THIS ENGINE MAY NOT READ

This is **engine 1** and it is **blind by construction**. Everything it computes comes from
`S11_SHARED_PHYSICS.md`. It must not read, import, grep, `git show`, or otherwise consult:

- `research/pde_ledger_v3/scripts/` — ⛔ **the sibling engine**, and every other SymPy engine in the ledger.
- `research/pde_ledger_v3/reduction/` — ⛔ **the registry, in any form.** This engine derives dimensions from
  the action alone. The registry comparison is engine 2's job and it is only a real comparison if this
  engine never saw it.
- `research/pde_ledger_v3/steps/` — ⛔ **all of it, every file, without exception.** ⚠ Some filenames in that
  directory state results; ⛔ **do not read the directory listing either.**
- `research/pde_ledger_v3/paper/` — the TeX cards and the built PDF.
- ⛔⛔ **`research/pde_ledger_v2/` — ALL OF IT, every file, without exception.** ⚠ That tree is this
  ledger's **predecessor**: it holds earlier engines *and their committed output files* for the same
  physics, under different names. ⛔ Do not read it, do not list it, do not `git show` it.
- ⛔⛔ **Every orientation, status and planning document in the repository**, at the repository root and at
  the ledger root — including `STATUS.md`, `REBUILD_HANDOFF.md`, `DEFECT_REGISTER.md`, `V3_STEP_PLAN.md`
  and any `SESSION_*.md`. ⚠ **These state results outright.** They are written for the orchestrator, and
  one of them exists specifically to tell the orchestrator things a builder must not be told.
- ⛔⛔ **`research/pde_ledger_v3/_asbuilt/` — all of it.** ⚠ It holds a frozen copy of an earlier engine for
  the same physics, kept only so the registry's line-range loci stay resolvable across a rebuild.
- ⛔ **Any other step's shared-physics file, engine, committed output, or checks config** — including the
  earlier-step engines and `out/` files sitting **in the very directory you are writing to**. ⚠ Several
  steps in this ledger share machinery with this one, and their committed outputs carry values for objects
  you are being asked to compute. ⭐ Compute yours; ⛔ do not look at theirs.
- ⛔ Any build transcript, log, or review report from any previous build.

⚠ **This list is not a spelling test.** ⛔ It is not satisfied by avoiding the exact paths named while
reaching the same content another way — `git show`, `git cat-file`, `git archive`, `git log -p`, a
`__pycache__` artifact, a symlink, or a file that merely quotes one of the above. ⭐ **The rule is the
intent: you compute every number in this build yourself.**

⭐⭐ **Two engines exist so that they can DISAGREE.** A disagreement is the single most valuable output this
build can produce, and it is destroyed the moment either engine is written against an answer it read
somewhere.

⚠ **You are overwriting a file that already exists.** ⛔ **Do not read the existing
`S11_stray_longitudinal_mathematica_audit.wl`, do not `git show` it, and do not adapt it.** It carries the
defects this rebuild exists to remove — it ends in a `VERDICT` tag, it renders symbolic booleans as the
typed words `TRUE`/`FALSE`, and its tag decomposition is not the one §8 now requires. ⭐ Write a new file
from the shared physics.

## The mandate

Implement **`Q1` through `Q11`** of the shared physics, for **every package** in its §7 table, at **every
`D`** in that package's sweep, obeying:

- §4, the structural rule;
- §5, the three clauses, the five corollaries, **no verdict**, and ⭐⭐ **the locus protocol**;
- §6's ordering requirement — ⭐ **`Q9` is computed first at each `D`**, because one package's action is
  built from its output;
- §7's `P_D` rule;
- §8's tag grammar, with `<ENGINE>` = **`WL`**.

## Mathematica-specific requirements

### Emission

- ⭐ Emit payloads with `ToString[expr, InputForm]` so every value is re-parseable. One line per tag,
  `TAG: <payload>`, per §8.
- ⭐⭐ **Emit each tag the moment it is computed, and flush.** ⛔ Do not accumulate tags in a list and print
  them at the end of a package, of a dimension, or of the run. ⚠ From outside, a run producing no output is
  indistinguishable from a solve that will never return — and that, ⛔ not elapsed time, is the failure
  mode this guards against.
- ⛔ **No `VERDICT`, `PASS`, `FAIL`, or `checks` list anywhere in the file.**
- ⛔ **No `Exit[1]` on a physics outcome.** Exit non-zero **only** if the kernel cannot complete the
  computation. ⚠ A builder iterating to exit 0 can make a genuine disagreement disappear.
- ⛔ **Every emitted tag name must be unique.**

### Assumptions

- ⭐ Carry §3's assumption set — **as one joint `And[...]`, exactly as written there, with the package's own
  domain premises joined in when that package is selected** — into every `FullSimplify` / `Simplify` /
  `Refine` / `Reduce` / sign test. ⚠ A rank or a sign computed without it is a different quantity, and
  engine 2 carries the identical joint set. ⭐ Emit the joint assumption expression itself as a tag.
- ⛔ **Add no premise relating the stiffness coefficients to one another** — not genericity, not
  distinctness, not an ordering. §Q8 asks for exactly the locus such a premise would delete.
- ⛔ Add no assumption about the sign of `ω²` anywhere in the construction of `M`.

### Solving, and the locus protocol

- ⚠⚠ **`_SOLUTION` is emitted with the domain LEFT UNRESTRICTED**, per §5's protocol table — the reality
  question is answered by `_REAL_ADMISSIBLE`, which is a separate object with its own operands. ⛔ Do not
  fold a domain restriction into the solve and call the protocol satisfied.
- ⭐ `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` are computed **from `_EQUATIONS`**, ⛔ never read off an
  empty solver token. Emit all five objects for **every** locus the spec requests, unconditionally.
- ⚠ `Solve[…, omegaSquared]` may return `ConditionalExpression`. ⛔ Do not drop the condition — emit it as
  its own tag with the `_LOCAL_` infix per §8. If a payload would print as a bare
  `ConditionalExpression[value, condition]`, emit **both** parts separately.

### The computations

- ⚠ `MatrixRank` returns the **generic** rank — that is the subject of `Q8`, so compute the rank-drop loci
  explicitly rather than assuming genericity, and ⭐ **at S11 the solve-variable sets include the stiffness
  coefficients**, not the wavevector alone.
- ⭐ For `Q4` **N3**, build the stacked matrix explicitly, e.g. `Join[Mr, {kvec}]`, and take its
  `MatrixRank`. ⛔ Do not attempt to infer `nu_T` from the `NullSpace` basis — §Q4 says why, and §7 says
  which package that warning is aimed at.
- ⭐ For `Q6`, walk the expression tree. The coefficient dimensions are **unknowns to be solved for**, and
  `COEFFICIENT_ORDERING` is built from the action's own terms per §Q6, ⛔ never from a fixed list you type.
- ⭐ For `Q7`, compute `c_i` from the Levi-Civita definition (`LeviCivitaTensor` or `Signature`) and compute
  `c·c` from **that result**. ⛔ Do not write out its expanded polynomial.
- ⭐ For `Q9`, emit `V1_BASIS`, `V2_BASIS` and `V6_BASIS` as the **reduced row echelon form** of the
  coefficient matrix in the pinned `MONOMIAL_ORDERING` coordinates — `RowReduce` on exact rationals.
  ⛔ Do not substitute a numerical or floating-point solve for the exact one.
- ⭐ `V6_OPERATOR` is the map `p ↦ p(R₀ G R₀ᵀ)` **as a linear operator on the `V1` space** — emit the
  operator, then its `(−1)`-eigenspace. ⛔ Do not obtain `V6` by subtracting `V2` from `V1`, and ⛔ do not
  obtain `V2` as the `(+1)`-eigenspace of `V6_OPERATOR`; §Q9 says what either shortcut destroys.
- ⭐ Set `$HistoryLength = 0`. If a simplifier becomes the bottleneck, prefer `Together` / `Cancel` /
  `FullSimplify` with the joint assumptions, and ⭐ **record which simplifier was used** as an engine-local
  `_LOCAL_` tag per §8.

### ⚠⚠ One trap that would hit BOTH engines identically

§Q3 asks for `ROOT_SOLUTION_SET` **"as returned, retaining multiplicity"** and `ROOT_DISTINCT` **"after
de-duplication"** as **two different objects**, with `ROOT_COUNT_ALL` and `ROOT_COUNT_DISTINCT` each counted
from its own list.

⚠ **The obvious root-finding call in either CAS silently de-duplicates**, which would make those two objects
identical **by construction** and their two counts equal for reasons that have nothing to do with the
physics. ⛔ That is not a result; it is the solver's default.

⇒ ⭐ **Produce `ROOT_SOLUTION_SET` by whatever route in your CAS actually retains multiplicity.** ⛔ This
directive names no API for it — ⭐ if you cannot obtain a multiplicity-retaining solution set, ⛔ do not
fake one: emit what you have and **report it under §10 item 3**.

## ⭐ If an object cannot be built the way the spec asks

⚠ Some objects in `Q9` and `Q11` may have only one obvious route in a given CAS, and §Q9 and §5 corollary 3
**forbid** particular routes because they make a residual identically zero. ⛔ **If the only route your CAS
offers for an object is one the spec forbids, do not take it.** ⭐ Emit what you can build honestly and
**report it under §10 item 3** — that is a finding about the specification and it is worth more than a
complete tag set obtained by taking the banned route.

⚠ The same applies to any `<QUANTITY>` name in §6/§7 you cannot emit under §8's one-tag-per-object rule.

## Runtime — ⛔ there is no wall-clock deadline

⚠⚠ **The rule is OBSERVABLE PROGRESS, ⛔ not total elapsed time.** A run that visibly completes
`(package, D)` cells for hours is **acceptable and expected**; §Q9's unknowns grow as `D²(D²+1)/2`.

- ⛔ **Do not narrow the declared sweep**, drop a `(package, D)` cell, or reformulate a requested object
  into a cheaper one to save time.
- ⭐ `RUN_PAIRS` and `SKIPPED_PAIRS` are **observed, not declared** — accumulate a pair only **after** that
  package has finished emitting, and emit both **after** the sweep, with `SKIPPED = declared ∖ completed`.
- ⭐ If a cell genuinely cannot complete, record it in `SKIPPED_PAIRS` **and report it in §10**. ⛔ Never
  drop one silently.

## Verify before you finish

Run the script over the **complete declared sweep**, with stdout redirected to the `out/` path above. It
must exit **0**. ⭐ **Paste nothing from its output into your report.**

## Operational

- ⛔ **Never `pkill` or `killall`** — shared machine. Kill by explicit PID only.
- ⚠ **At most two Mathematica kernels at any time**, two-seat licence.
- ⛔ No `git add`, no `git commit`, no other git write.

⇒ Then answer **§10 of the shared physics, under 25 lines.**
