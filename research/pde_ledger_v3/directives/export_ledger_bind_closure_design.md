# Export architecture redesign — the bind-closure LEDGER (design / decision list)

**Orchestrator-authored, 2026-09-03; v2 folds the two decision legs
(`_legs/export_ledger_bind_closure_design_review.md`; Codex + Grok, both "not sound as-is", convergent — outputs
`_legs/export_ledger_bind_closure_design_review_{codex,grok}.out`).** This is a **decision list** (rule 7). It
settles **what the plain-git `*_exports.py` LEDGER must contain and how it is stored**, folding the two
deliberation legs (`_legs/export_ledger_minimum_{codex,grok}.out`) and reconciling the prior author decision
`chain_accumulate_or_generate_decision.md` (2026-08-11) the implementation drifted from. It does **not** decide
expression representation (§D4) and adds no physics.

**v2 changed (folded from the review legs):** the availability mechanism — F9 is **write-time** and the fold
never re-applies it; **no write-time pass-through** (the consumer's `IMPORT_KEYS` discovers everything from the
full fold); `IMPORT_KEYS` are **F9 write-keys**; the closure guard resolves **every edge's route**, not just
roots; the smoke-test is **bidirectional**; a **missing manifest is a hard error**; `.out` recovery needs an
**explicit promotion-delta**; `coupling_kernel` is **dropped** (c2 re-extracts, does not bind it); the **c1
spec §7 and N1** are amended, not just the directive.

## 0 · The measured problem (both deliberation legs, identical file)

`S11c_b_exports.py` = **58.75 MiB (61,608,243 bytes), 2,441 rows**. Measured over-export:

- **~56% (~33 MiB) is provably unbound** by any declared later step: `*_term_origins` provenance
  (**17.35 MiB, 29.5%**), dead accumulated S9/S10/S11/S11b history rows nothing references (~15 MiB), and
  unconsumed intermediates (`energy_basis_*`, `admissibility_*`, unreferenced `projection_*`).
- The **bound remainder** is two rows — `slab_operator` (18.1 MiB) and, if c2 declares it, `coupling_kernel`
  (7.6 MiB) — in **expanded** form (271 distinct symbols, 365,820 occurrences). That is a **representation**
  problem, out of scope (§D4). The minimum bind-closure is **~319 rows / ~18.2 MiB** without `coupling_kernel`,
  ~508 rows / ~25.5 MiB with it.

Two root causes: (1) membership was "every §4 MAIN primary" (D1) instead of downstream use; (2) each write
**accumulates** the whole inherited file forward, so dead rows ride along and the file grows monotonically.
`F10` (2026-09-03) narrowed only a step's **own** new rows and left the carried-forward inherited file whole —
so it does not shrink the chain. (Note: after F10's dissipation trim, c1's *current* file is ~64 MiB — already
under the 100 MiB cap; this redesign stops the **monotonic growth** and the over-export, it is not a live cap
fire.)

## D1 · MEMBERSHIP — the bind-closure test (both legs converged independently)

**A row is EXPORTED iff it is in the BIND-CLOSURE**, defined at the boundary after step *i* as the recursive
closure of

> the set of **F9 write-keys** that some **already-declared later step** binds **by key** as a **construction
> operand** — composes, substitutes, solves, evaluates, differentiates, or binds as a **control/regression
> operand**,

where the closure additionally pulls in, **as resolved F9 write-keys**: every **symbol / knob / profile**
referenced in a kept row's `value`; and every **`dimension_key`** target, recursively.

- ⛔ **Membership is decided by BINDING ALONE.** `class`, `step`, "primary", "model-level vs diagnostic" do
  **not** decide it. A diagnostic-labelled row a later step binds as a regression operand (c1 binds S11b
  `z_impermeable`, `z_by_regime`, `z_by_parity`, `added_mass`, `grazing_behaviour`, `face_response`(flat),
  `face_response_coeffs`, `permeable_dissipative_by_regime_and_parity`, and the five `degenerate_loci_*` for
  §5c/§5d) **is exported**; a model-level-looking row nothing binds (`energy_basis_variable`) is **emit-only**.
- ⭐ **Keys are F9 WRITE-KEYS** (bare on F9a/F9b, producer-prefixed on F9c). A closure edge is resolved to the
  *actual* row it binds, ⛔ never to a bare name that F9c has re-pointed to a predecessor (§D3).
- Everything not in the bind-closure is **EMIT-ONLY**: printed to stdout → the annexed `out/*.out` (git-annex +
  GIN, no cap), read by the comparator and the review legs, ⛔ never imported.

This **replaces `D1`'s membership rule** ("export every MAIN; errs toward more"). What survives is its one real
control: **the builder makes no per-object judgement** — the *spec* names the bind-set; the engine exports
exactly its closure. ⛔ **A missing manifest is a HARD ERROR (§D3), not a fallback to export-every-primary** —
that fallback is precisely what put the 17 MiB of `term_origins` in plain git.

## D2 · TOPOLOGY — generate over a frozen base; F9 at write time; fold = last-wins

Implement the **2026-08-11 decision** (`chain_accumulate_or_generate_decision.md`: *"each step exports only the
rows it defines or re-derives; the whole-model list is generated on demand"*), which superseded the
accumulation rules (D5, D7's marker, F4) and was never built. `F10`/`N1`'s "carry the whole inherited LEDGER"
re-asserted the superseded topology; it is reversed here.

**Deltas.** From S11c-c1 onward, each `*_exports.py` contains **only the rows that step defines or re-derives**
(its own delta) — a few MiB, never the accumulated model.

**F9 is a WRITE-TIME rule; the fold never re-applies it.** When a writer emits a row whose object name collides
with a key already in the **chronological prior fold**, it applies F9 *then* (F9b bare on proved-equal, F9c
producer-prefixed otherwise) and stores the **final F9-resolved key** in its delta. The `.out`/step-record
records the route.

**The fold.** A consumer resolves keys via `load_model(base, *deltas)`: the **single atomic frozen base**
`S11c_b_exports.py` (already F9-resolved internally) merged with each later delta in **validated chronological
order, last-wins on exact key**. ⛔ The fold does **not** re-run F9 (no re-prefixing, no re-comparing resolved
pairs); it verifies the declared routes and reports a route mismatch. The consumer then takes its own
`IMPORT_KEYS` closure (§D1) over the full fold. ⛔ **There is no write-time pass-through**: a row a *later*
consumer needs is discovered by *that consumer's* `IMPORT_KEYS` over the fold (c1 never names `slab_operator`;
c2 lists it and the frozen base supplies it).

**Migration = freeze the base, do NOT re-run the chain.** The existing accumulated files through
`S11c_b_exports.py` (< 100 MiB) are the frozen base: they hold the whole model through S11c-b, so no upstream
engine re-runs and **no `BUILD_INPUT_DIGESTS` sha cascade fires** (each file pins its parent's content sha; a
byte-identical base breaks nothing). New steps write own-rows deltas; the fold is `base ⊔ deltas`. This stops
the growth with zero re-execution. ⚠ It does **not** eliminate future whole-file-digest coupling (the 2026-08-11
decision left that provenance contract unchanged) and does **not** reclaim the base's own ~33 MiB of dead rows
(an optional later cleanup that *would* trigger the cascade — out of scope).

## D3 · UNDER-EXPORT GUARD — a manifest + edge-wise route/closure check (NOT D3-round-trip)

⚠ **The `D3` round-trip does not catch under-export** (both legs, verified): `_restore` builds `Symbol(...)`
directly from each row's srepr (`S11c_b_exports.py:17`) with no LEDGER lookup, and `publish_export` only
compares rows that **were** written — so a value round-trips even when a referenced symbol's **row** is absent.
`F10`'s "D3 catches the closure gap" is **false** and is struck. The guard is mechanical, at generate/fold time:

1. **Manifest.** Each consumer publishes an exact **`IMPORT_KEYS`** list of the **F9 write-keys** it binds (as
   construction or control/regression operands). ⛔ **No `PASS_THROUGH_KEYS`** — a later, non-adjacent consumer
   declares its own binds against the fold. **A missing manifest blocks publication of the delta (hard error).**
2. **Edge resolution by FULL SYMBOL IDENTITY, not name (measured refinement, 2026-09-03).** Compute the
   recursive `dimension_key` + referenced-atom closure; resolve **each edge** against the fold by the atom's
   **full identity** (its `srepr` — name **plus assumptions**), ⛔ not its bare name. Three outcomes:
   - **exactly one** fold row shares that identity ⇒ **resolved** — this is why `omega=Symbol('omega',real=True)`
     resolves to the bare `omega` row, not to `s11b_omega=Symbol('omega')` (distinct identities);
   - **zero** fold rows ⇒ the atom is a **free/structural symbol or function** (a coordinate, or a structural
     function like the window `O_window` — 0 rows, referenced 3,284× as `Function('O_window')`) ⇒ it is **not a
     required dependency**; ⛔ do **not** demand a row for it (requiring one was the over-reach that
     false-positived);
   - **two or more** fold rows share the identity under **different keys** ⇒ **genuinely ambiguous** (the
     unresolved routed-key contract, `S11_export_chain_decisions_v2.md:205`) ⇒ **raise-and-name**.
   ⚠ A referenced-name→producer index keys on **identity**, ⛔ never on "any row that USES the symbol": a row
   like `f_hold_e_W_0` (class `PREMISE`) that merely *mentions* `Symbol('W_0')` is a **user**, not a producer of
   `W_0`. **Measured (2026-09-03):** 16 F9c routed pairs exist — **11 resolve by identity, 5 are genuine
   same-srepr ambiguities, and 0 of the 5 are referenced by any critical-path root** (only `omega`, resolvable,
   is). The 5 off-path genuine ambiguities are the **deferred residual routed-key contract** — documented,
   resolved when a consumer actually binds one, ⛔ not a blocker for c1/c2. ⭐ The real under-export guarantee is
   the smoke-test (step 3): a build that looks up a genuinely-absent required row **fails there**.
3. **Bidirectional smoke-test.** The consumer builds against an **access-recording ledger proxy**; assert the
   set of **observed** construction/regression lookups **equals `IMPORT_KEYS`** — catching both an undeclared
   lookup (under-export risk) and a declared-but-unused key (stale manifest).
4. **Minimum-mode.** A delta's exported key set equals its own bind-closure contribution, apart from explicitly
   named infrastructure — an accidental re-accumulation fails loudly. (With the hard-error manifest, "export
   every primary" is never a selectable policy.)

**Recovery is an explicit promotion-delta, not a manifest edit.** The `.out` is the data-recovery net, but
adding a key to a later `IMPORT_KEYS` does **not** materialize the row — the fold correctly reports it missing.
To promote an emit-only object into the chain, its producer emits a **promotion-delta**: the row reconstructed
with its full schema/evidence (`display`, `_restore` srepr, `class`, `f9_operands`, route), or the producer is
rebuilt. ⛔ Never a PY re-derivation of an inherited object at consume time (`N1`; only WL re-derives, importing
nothing). "Emit-only is never computed then discarded" — but recovery is a real artifact, not a name.

## D4 · REPRESENTATION — out of scope, but recorded

The minimal **set** still keeps `slab_operator` (18.1 MiB) because c2 folds into it; its size is the
**expanded** representation plus the `display`+`srepr` double-store `D3` requires. ⭐ The comparator
canonicalizes the *difference* itself — `sp.cancel(sp.together(sp.expand(A−B)))`
(`S11c_b_cross_engine_comparator.py:1206`) — so agreement does **not** require the LEDGER to store the expanded
form. Compacting the stored representative (unexpanded / named coefficient sub-objects / dropping `display` from
the LEDGER — which would require explicitly superseding `D3`'s human-readable-rendering clause) is a **separate
lever** and its own workstream; ⛔ not solved here. This design changes only the row **set** and the storage
**topology**.

⚠ **`coupling_kernel` is NOT in the minimum.** c2 **re-extracts** the closed off-diagonal kernel from the
**closed full operator** (`S11c_c1_SHARED_PHYSICS.md:34`); it does **not** bind S11c-b's **open**
`coupling_kernel` (substitution into the already-extracted open kernel is "at most a residual after proving
per-channel commutation" — `_measurements/S11c_c_spec_review.md:59`; the extract-then-close non-commutation is
the reason). By the bind test it is emit-only unless c2, when authored, explicitly declares it a **regression
operand** — at which point c2's `IMPORT_KEYS` promotes it and the frozen base supplies it. ⛔ Do not pass it
through c1.

## Reconciliation with the existing record

- **`D1`** — membership rule replaced by §D1; the "no per-object builder judgement" clause retained; the
  "every MAIN / errs toward more" default is **removed** (a missing manifest is a hard error, §D3).
- **`F10`** (2026-09-03) — **superseded.** Its correct half (the EMIT-vs-EXPORT channel split, the
  `.out`-is-the-comparator-channel fact) survives in §D1/§D3; its wrong halves (category membership, "carry the
  whole inherited LEDGER", "D3 catches the closure gap", and the export-every-primary fallback) are struck.
- **`chain_accumulate_or_generate_decision.md`** (2026-08-11) — **un-superseded and implemented** by §D2, with
  the frozen-base migration it left open.
- **`N1` and the c1 spec §7 carry sentences** — **amended** by §D2: a consumer imports the **fold**
  (`base ⊔ deltas`), not "the prior one's exports"/"carry the whole LEDGER". This amendment is part of this
  decision (the spec wins over the build directive, so the directive alone cannot enact it).
- **`F9`** — untouched: it decides **write-keys** (which spelling a collision writes), not which objects exist;
  §D2/§D3 make it strictly write-time and require the manifest/fold to carry the resolved route.
- **`D3` round-trip / `F6` / `F3`** — retained (round-trip only what is written; publish-completeness; a
  re-derived row carries its evidence). ⛔ `D3` is **not** a closure guard (§D3 supplies that).

## Build plan (after this list passes two legs)

1. **The fold + guard module** (orchestrator-directed, Codex-built, its own 2 legs): `load_model(base, *deltas)`
   (chronological last-wins, route-verifying, ⛔ no F9 re-apply); the §D3 manifest / edge-wise route / closure /
   bidirectional-smoke-test / minimum-mode assertions; the promotion-delta helper. Because a consumer that
   imports this module gains a **6th executable input**, either **inline it into each self-pinned audit** or
   **add it to `BUILD_INPUT_DIGESTS`** — decide and state it here; the c1 five-input digest must be corrected
   accordingly.
2. **Amend the c1 spec §7 + N1** (per Reconciliation) so the authority no longer says "carry the whole LEDGER".
3. **Update the c1 SymPy build directive** to §D1–§D3: emit **own-rows delta only**; declare c1's `IMPORT_KEYS`
   **explicitly and exactly** — `mu_theta_operator`, the nine T-substrate rows, the constants/closure knobs
   incl. `v_bulk_normal_0`, `q_out`, `Lambda_{A,V,X}_0`, `tau_{A,V,X}`, the varying profiles/density reps, and
   the S11b §5c/§5d regression rows (⛔ not "closure will find them"); its **outgoing** rows are the **exact five
   roots** `dtn_flat_symbol`, `dtn_operator`, `dtn_kernel`, and the F9c write-keys `s11c_c1_face_response`,
   `s11c_c1_face_response_coeffs` (⛔ **no `dtn_*`/`face_response*` wildcard** — that re-admits the diagnostics
   F10 removed); import via the fold over the frozen base.
4. **Rebuild c1** under the new architecture (export = a small own-rows delta) → then resume c1's dual-engine
   path (WL engine, comparator, reconcile, step record). The blind **WL engine writes no ledger** (unchanged) —
   this design is PY-chain only.
