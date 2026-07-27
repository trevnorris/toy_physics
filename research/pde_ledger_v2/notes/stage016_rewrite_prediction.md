# stage016 rewrite — sealed pre-registration, and its adjudication

**Sealed** before the build directive was launched, held **outside the repo** (`/tmp`, unreadable by
Codex and by every review agent) until the build *and* both review legs had landed. Copied in here at
step (i) as the record, per `manifests/DIMENSION_REWRITE.md` §4-e/§4-h2.

⚠ **Custody caveat, stated honestly:** committing it afterwards does **not** prove it predated the
build. Its value is as a *working* pre-registration — it disciplines the orchestrator's own thinking
and records falsified predictions — not as cryptographic evidence.

**Score: 6 confirmed · 1 falsified · 1 falsified-in-part · 1 confirmed-with-a-caveat.**
The two failures are both about *me*, and the second is the more useful one.

| # | prediction | outcome |
|---|---|---|
| **P1** | `evalDimensional` runs **6×** (1 clean + 4 corrupted + 1 arity re-run); emission must go at top level, not in-body | ✅ **CONFIRMED** — build independently measured 6; emission sits in `runDimensionalGate[]`, called once |
| **P2** | axis order `(L,M,T)`, established by a **slot→label binding in code**, not by prose; no slot identically zero | ✅ **CONFIRMED** |
| **P3** | ⚠ "Codex emits only the ~9 computed vectors and **omits the ~12 rule-table entries** as not worth comparing" | ❌ **FALSIFIED** — it emitted **all 21** and justified each in the enumeration. My reasoning was wrong in the good direction |
| **P4** | ≥9 emitted records are declared literals in both engines; the largest literal fraction yet | ✅ **CONFIRMED — 12 of 21**, independently counted by the build and the physics leg |
| **P5** | `K_eta` collides cross-stage as a **reduction level**, not drift; same shape for `T_w`, `μ_η` | ✅ **CONFIRMED** — 013 line / 016 volume / 023 reduced scalar, none renamed apart |
| **P6** | `a_dim` (016) **is** stage018's `a` — *recorded as unverified, since this is the `c_s0`/`c_S` shape that fooled two parties at 018* | ✅ **CONFIRMED** by two blind parties (build + physics leg), each citing its own evidence |
| **P7** | PASS multiset unchanged at **91** | ✅ **CONFIRMED** |
| **P8** | my read-only survey's **21** dimension-valued objects = what an independent enumeration finds | ⚠ **FALSIFIED IN PART — the important one.** The build's table has **35** rows (21 emitted + 14 excluded), and the review legs then found **my survey had missed two objects entirely**: `lambdaRef` (`.wl:227`), a live factor in the walk whose dimensionlessness comes from a `NumericQ` fall-through rather than any declaration, and the clean 6th `evalDimensional` re-run (`.wl:616`). **Both are symmetric omissions — exactly the failure mode §4-a exists to catch, committed by the survey that was supposed to prevent it.** |
| **P9** | `REACHABLE` ⇒ a print-only stage under D2, no new derivation and so no D3 independence argument owed | ✅ **CONFIRMED**, with a caveat: no *derivation* was added, but a single axis structure **was** introduced and the pre-existing renderer repointed at it. That is a structural edit, not a bare print |

## What P8's failure means beyond this stage

The `.wl` reachability survey method — the same one used to produce the **023, 027 and 021** verdicts
now recorded in `manifests/DIMENSION_REWRITE.md` §8 — **demonstrably misses objects**. It missed a
numeric factor that carries a dimension by fall-through, and it missed a repeat invocation whose result
is discarded. Neither is exotic; both are the kind of thing that is invisible precisely because nothing
consumes it.

⇒ Treat those three surveys' object counts as **lower bounds**, not inventories. The per-stage §4-a
enumeration, written by the builder against the file and then attacked by two review legs, is what
actually closes the set — the survey only decides the *shape* of the work. Recorded so the next stage
does not inherit a false sense of completeness.

## What the sealing was for

At stage018 an untracked prediction note was left in the working tree; a review agent read it and cited
my own wrong claim (P9 there) back as authority. Holding this file in `/tmp` closed that channel —
verified by grepping the build directory for leaked answers before launch, and by quarantining the two
review reports out of `_scratch/stage016/` before the builder ran.
