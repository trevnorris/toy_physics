# Script rebuild — state and what's next

Read `CLAUDE.md` first (how we work), then this (where we are).
Process detail: `docs/development_pipeline.md`. Defects: `DEFECT_REGISTER.md`.

Last updated 2026-08-05, end of session.

---

## Where S10 stands

| | state |
|---|---|
| **The mode count** — `ν_T = D − rank([M_r; kᵀ])`, `D−1` transverse | **Established.** Both engines compute the stacked rank independently; a leg confirmed all 26 generic `N3` tags present in both with no payload mismatch. Not re-verified this session. |
| **Q7** — the `D=3` form-and-normalisation comparison | **Broken.** SymPy hardcodes a fixed curl density and ignores the package (`…sympy_audit.py:1538`); Mathematica uses the package's own (`…wl:1312`). The spec instructs *both*, which is why. Repair directive ready. |
| **Cross-engine comparison** | **Does not run.** The harness reader rejects 220 of the Mathematica engine's 2983 payloads, including the Lagrangian and the Euler–Lagrange system, so the dimension and ablation layers never see the action. |
| **The step record** | **Reverted.** Must be rewritten with the disclosures below. |

**The critical path is the harness, not the spec.** Everything else is downstream of being able to compare
the two engines at all — and the Q7 directive now explicitly says provenance rests on cross-engine
comparison, so that claim is load-bearing.

---

## Scheduled after the rebuild — the substrate requirements register

`SUBSTRATE_REQUIREMENTS.md` (seeded 2026-08-05, two entries). The light sector is built on **S1–S8, none
of which have been run**. That is the requirements-first design, and it pays only if what each sector
needs from the substrate is written down. Measured: two obligations are recorded, both found by accident;
the rest is implicit prose across four step records.

**Pass 1 — S9 and S10** (user decision, 2026-08-05): they were rebuilt carefully, so they calibrate the
entry shape. **Pass 2 — S11, S11b, then each sector as it closes.** At Phase 5 the file becomes the
checklist the substrate is tested against.

---

## What's next, in order

1. **Harness rebuild** — `directives/S10_harness_rebuild_directive.md` (commit `bed199a8`).
   Needs two legs (orchestrator-written → Codex + Grok), then build, then two legs on the build.
   `D8` is the item that matters: the reader hardcodes a derivative shape from S9, which swept one
   dimension, while S10 sweeps four.
2. **Q7 spec repair** — `directives/S10_spec_Q7_repair_directive.md` (commit `fcda865a`).
   Three review rounds done; the finding class has moved from "this check cannot fail" to "this sentence
   claims too much", which is convergence. Build, then two legs.
3. **Rebuild SymPy's Q7** so it emits the density its action used.
4. **Write the step record.** Orchestrator's job. Must carry every disclosure below.
5. **TeX card** — Codex, with its own two legs.

---

## Deferred deliberately — kept, not lost

Scope decision (user, 2026-08-05): **keep S10 to S10.** These buy future comparability, not S10's
correctness.

- **The contract rewrite** — `directives/S10_spec_rewrite_directive_v2.md` plus
  `directives/S10_contract_DEFERRED_findings.md`, which holds two legs' findings on it including two
  attacks that got through and an unconstructible oracle node. Neither applied. Revisit before S11.
- **`DEFECT_REGISTER.md#f7`** — the comparison kernel equates a boolean with any nonzero number.
  Pre-existing. Measured exposure: **zero** cross-engine pairs on either S9 or S10 put a boolean on either
  side, so no agreement number is contaminated. The harness directive asks for a **tripwire** rather than
  the fix, so the deferral stays safe by construction rather than by a measurement that ages.
- **§8's quantity registry, the payload serialiser, the sign classifier, the Q6d unknown count.** All
  known-defective, all recorded. One narrow exception is *in* scope: `§8` has no stratum scope token, so
  the two engines emit different raw names for the Q8 stratum mode counts (`E7` in the repair directive).

---

## Disclosures the step record must carry

- **Q7's verdict was supplied.** `S10_SHARED_PHYSICS.md:384-385` tells both engines the residual is
  expected to be nonzero for non-curl packages, and `:380-382` states the measured failure shape. The
  residual is a computed polynomial; the judgement of what it should be was given.
- **The `homogeneous` counter is largely vacuous** — 1699 of 2139 counted-homogeneous tags had zero
  comparisons performed. S9's cited `1219/1219` is that same counter: disclosed, not verified.
- **Cross-engine coverage is a subset** — 677 configured pairs against ~2900 candidates. §8 pinned the tag
  prefix and never the quantity names, so both engines obeyed it and still diverged.
- **Provenance rests on cross-engine comparison**, not on Q7's residual. Nothing inside a single engine
  establishes that the emitted density is the one its action used — a leg proved the mutation test
  establishes dependence, not identity.
- **Prior art.** The light-sector algebra is MacCullagh's rotational aether (1839); curl-only stiffness
  giving `D−1` transverse modes plus a zero-restoring-force longitudinal is textbook continuum mechanics
  (Whittaker Ch. V). The record must **cite** it — it currently reads as though the mode structure were
  earned from our own postulates. Prior art is an **oracle, never a premise**.
- **S10 linearises about an unstrained brane** (`v₀ = 0`), which *is* MacCullagh's regime. So S10's
  agreement with him is not evidence our medium differs from a 19th-century aether — the divergence
  (compressibility, the background flow, the out-of-plane channel) is scheduled, not demonstrated.
- **What is genuinely ours** and should go in the Departure field: the `D`-sweep and the statement that the
  bulk never enters (he had no bulk); the in-plane/out-of-plane split; the longitudinal kept deliberately
  as the charge anchor. And for each of the three things that killed the elastic-aether programme —
  longitudinal modes, preferred frames, matter–aether coupling — what this model does instead, and whether
  it escapes or just relocates the problem.

---

## Pinned

- **`s10-as-built`** (`9309da70`) — the spec, both engines and both outputs as they actually ran. The step
  record cites these blobs, not any rewritten spec.
- **`wip-2026-08-05-unreviewed`** (`92461853`) — two builds that were committed without review and
  reverted. Reachable, not lost. Do not cherry-pick unreviewed.

---

## Known limits of the instrument

- The harness's ablation layer keyed on an S9 naming convention and matched **zero** S10 tags, reporting
  `compared=0 responsive=0 invariant=0` — which reads as health. That is what the rebuild fixes.
- Its coverage guard is scored against **what the engine emitted**, so deleting output improves it.
- A declared control absent from an entire dimension is silently dropped from every counter.
- `reduction/harness_ablation.py` reads counters with an integer-only regex, so it truncates a fractional
  coverage value to `0` — it cannot read the number the coverage guard rests on.
- Runtime margin is thin: modest form ablations pushed both engines past 600 s, always in the same package.
