# Launch prompt — paste this after /compact or at the start of a new session

⚠ **This file is a POINTER, not a summary.** The last session's launch prompt was a lossy hand-compression
that lost a whole review round and cost the first hour. ⛔ Do not let this one become that: it names
documents and states rules; it does **not** restate findings. If it and a named document disagree, ⛔ **the
document is right.**

---

Branch `ledger-v2-rebuild`. Run `git log --oneline -5` and `git status` first —
⛔ do not trust any hash written in a doc, including this one.

▶ **READ FIRST, AND IN FULL: `research/pde_ledger_v3/SESSION_2026-08-01.md`.**
It is the latest record and it GOVERNS. It carries the direction change, the model's foundational
postulate, the five standing decisions, what was banked, the simulation track, the open physics, and
§8 — how the previous session failed, by mechanism. ⛔ Do not start work before reading it.

▶ **THEN the orientation budget — and it is EIGHT documents, not the five `NEXT_SESSION.md` lists.**
The three it omits all govern:
- `docs/derivation_walkthrough_plan.md` — the method; **§5a** is why the registry is grown forward
- `research/pde_ledger_v2/walkthrough/DECISIONS.md` — ⛔ overturns the method doc where they disagree
- `docs/development_pipeline.md` — the operative process doc

▶ **`NEXT_SESSION.md` is still worth reading for its ELEVEN-ITEM list** of what a fresh session gets
wrong — each item is something an oriented orchestrator actually believed. ⛔ But its ordering and status
claims are stale; the session record says which.

---

## ⭐⭐ THE RULES — these are not style guidance

1. **Every step is walked SIDE BY SIDE.** Derive in the open, one move at a time, reasoning BEFORE the
   result. ⛔ Never pre-derive and present. ⭐ Flag every identification BEFORE making it.
2. ⭐ **Show the equations and the dimensions at every step** (user, 2026-08-01). Derive dimensions from
   the Lagrangian/action, ⛔ do not read them off the register and agree.
3. ⛔⛔ **Scripts are the volatile surface. Anything Codex writes gets a FRESH AGENT leg AND a GROK leg.**
   ⚠ The 2026-07-29/30 scale-back cut the *directive bookend*; it did NOT licence skipping review of
   builder-written code. A fresh agent of yours may review a file you made a minor fix to.
4. ⭐ **Every step adds its quantities and relations to `research/pde_ledger_v3/reduction/` and the gate
   must pass before the step banks.** The registry IS the requirements list — ⛔ do not defer it.
   ⚠ **Both engines:** Python gate *and* `research/pde_ledger_v3/mathematica/registry_dimensional_gate.wl`.
5. ⛔ **v3 is SELF-CONTAINED.** Copy files from v2; ⛔ never link to v2's tree.
6. ⛔ **The input COUNT is not the objective — physical motivation is.** Never fudge physics to shorten
   the list.
7. ⛔ **Verify a background job by its ARTIFACT** (`hook: Stop` / `tokens used` / mtime), ⛔ never by a
   process check — `pgrep -f "codex exec"` matches its own waiter, and it cost 1h42m once already.
8. ⛔ **Check a document's STATUS LINE before quoting it.** Two superseded docs were cited as live in one
   session, one carrying an explicit *"Do NOT build on it"* banner.

---

## ▶ WHERE WE ARE

**Banked:** `83668b97` S0.5 · `361e8114` v3's own `reduction/` + **S9** · `f2f5d9af` S9 reviews folded +
the Mathematica gate · `df57fc76` session record.

**S9 is done** — light's requirement on the medium, the first sector step under requirements-first.
Registry: 10 continuous quantities, 4 relations, residue = the four GNLS primitives + light's two
postulated constants. Both engines agree on that residue.

**▶ NEXT: S10 and S11** — the two transverse photons, then the stray-longitudinal departure. Then
**S11a**, the packet + simulation specification, which does not exist yet.

⚠ **Read the session record's §6 before S10** — it carries the `μ_br` ≠ `μ_R` hazard (Cauchy vs
MacCullagh, near-identical names, different objects) and the three-step story of light's departure.

## ⛔ OWED

1. ⛔⛔ **THE `.wl` GATE HAS THREE BLOCKING DEFECTS — ⛔ DO NOT TRUST ITS `WL_REGISTRY_PASS`.** Reviewed
   2026-08-01; see `SESSION_2026-08-01.md` §7 item 3 for the full findings and fixes. In short:
   **(1)** `--registry-dir` is dead and the orchestrator's "fix" **masked** it — the tracked script
   pointed at a known-**inhomogeneous** registry prints `HOMOGENEOUS=4 … PASS`, exit 0;
   **(2)** `qidsIn` misses a level-0 `Q` node, so a bare `[Q, x]` RHS — an alias/identification relation,
   a shape this project makes constantly — is admitted with `input_qids: []` and its output leaves the
   residue derived from nothing (⭐ Python does **not** have this hole; the engines genuinely disagree
   and only the `.wl` is wrong);
   **(3)** `kind` is never validated, so a one-character typo inflates the residue silently at exit 0.
   ⚠ **Fix these before building any `.wl` ablation harness**, and ⛔ assert on the `WL_REGISTRY_FAIL`
   **text**, not on exit ≠ 0 — a licence-contention run exits 40 with no gate output.
   ⭐ **Grok's leg CONFIRMED (1) independently** — `BLOCKING FINDINGS`, and it adds the sharpest form of
   the point: able-to-fail could only be demonstrated *"by hardcoding the registry path in a scratch
   copy of the script."* ⇒ ⛔ **No honest ablation harness can be built on that flag until it works.**
   ⚠ Grok did **not** find the level-0 `Q` hole; the fresh agent did. Different scopes, both needed.
   Full logs: `_scratch/grok_wl_review.txt`, and the fresh leg's findings in `SESSION_2026-08-01.md` §7.

2b. ⛔⛔ **NO `.tex` EXISTS FOR v3 — the step unit is INCOMPLETE.** v2 carries **44** stage `.tex` cards
   under `research/pde_ledger_v2/paper/stages/` plus a full paper (`parts/`, `appendices/`, `macros.tex`,
   a built PDF). ⭐ The defined unit is **six artifacts** — *note · TeX card · SymPy audit · independent
   Mathematica audit · source map · register entry* (`docs/model_map.md:269`). **S9 produced three:**
   the note, the registry entry, and the Mathematica gate. ⛔ **Missing: the TeX card and the source
   map** — and they are missing for S0.5 too. ⚠ Do not walk S10 without deciding the v3 paper's
   structure; writing 20 cards retroactively is the retrofit this method exists to avoid.
2. **The restructure to requirements-first is NOT applied to the docs.** Attempted twice, failed review
   twice, rolled back twice. The charter still says *"this is not a third method change"*; the step plan
   still runs the substrate first. ⭐ **Do the minimum to walk the next sector, ⛔ not the whole document.**
   Decision lists: `_scratch/RESTRUCTURE_V2.md`.
3. **A13 has four mutually exclusive homes** in the step plan, and *"settled at S12"* is unsupported.
4. **Two `able_to_fail.py` defects, filed unfixed** — `demonstrate_verdict_spoof` ⛔ cannot fail; the
   crash/spoof flags trip the case-count guard before their mechanism.
5. **Two registry nits, recorded not fixed** — `linear-unstrained-brane` is an inert tag;
   `denominator_guards` are declarative only.
6. ⚠ **A second S0.5-shaped question, undecided (USER CALL):** S0.5's own criterion now indicts all eight
   surviving medium rows, since the substrate steps run last. `ħ`/`m` come standard and stay; `K` and the
   EOS exponent are explicitly *not* standard.

## ⚠ WHY THE `.wl` GATE MISBEHAVED WHEN v2's 44 NEVER DO — measured

⛔ **Do not read the `.wl` blockers as "Mathematica is flaky."** v2's 44 stage scripts declare it in their
own headers — *"Print-only, standalone, **no arguments**, no exports"* — and across all 44:

| | count |
|---|---|
| use command-line arguments | **0** |
| import YAML | **0** |
| read **any** external data file | **0** |

⇒ They are **self-contained symbolic verifiers**: hardcode, compute, assert, print. Deterministic.

⭐⭐ **The v3 registry gate is the FIRST `.wl` here that takes an argument, parses YAML, and reads mutable
shared state.** All three are novel, which is why nobody had ever found that `$ScriptCommandLine` is empty
under `math -script` in this environment — no script had ever needed an argument.

⇒ ⭐ **All three blocking defects sit in that novel surface** (dead argument flag · parser hole on a data
shape · unvalidated schema field). ⛔ **None** is in the symbolic-algebra surface where the 44 have their
record. ⇒ **A reader of shared mutable state has a failure class a standalone verifier structurally
cannot have** — wrong file, stale file, silently-defaulted path. ⚠ Treat every future registry-reading
script, in either engine, as carrying that class; ⛔ the 44's track record does not transfer to it.

## ▶ NEXT SESSION — in this order

1. ⛔ **Fix the `.wl`'s three defects** (Codex builds; then a fresh agent leg AND a Grok leg). ⚠ One of
   the three is the orchestrator's, so the orchestrator ⛔ cannot judge the repair.
2. ⛔ **Decide v3's paper structure and backfill the TeX cards + source maps for S0.5 and S9** — two
   cards while it is still two. ⚠ **USER CALL:** does v3 get its own paper skeleton, or write into v2's?
   Given v3 is now self-contained, its own is the coherent answer.
3. **S10** — two transverse photons. ⚠ The count comes from the brane being **3-dimensional** (a 3-vector
   about `k` splits 1 longitudinal + 2 transverse), ⛔ **not** from codimension.
4. **S11** — the stray longitudinal, framed as the **third of three failed attempts** at Maxwell's
   no-longitudinal demand, ⛔ not a bonus mode.
5. **S11a** — the packet + simulation specification. ⭐ The one with real thinking in it.

## ⭐ THE SIMULATION TRACK — the user wants this

A superfluid sim with a photon in it, time stepped slowly, watching the medium move — and whether
something stays **quantized in a packet**. ⚠ That is a soliton and it needs **nonlinearity**, which the
corpus does not have (**C6**: no closed parent action; the same missing object as the geon).
⭐ **Cheaper first target: Test 5** (`notes/brane_bulk_handoff.md:880`) — *"check that setting `χ_B = 0`
removes shear/light propagation"* — designed, **never run**, needs only the linear equation plus the
gate, and tests a live requirement of light's.
