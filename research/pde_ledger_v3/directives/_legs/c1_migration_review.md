# Independent review — the S11c-c1 EXPORT-ARCHITECTURE MIGRATION (orchestrator-written; document + real-artifact)

## What changed and why
Mid-c1-build the LEDGER was found to over-export ~56%. A two-leg-gated design
(`research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md`, committed `c04e071f`/`fd8c89d0`)
replaced the export chain's **accumulate** topology with **generate-over-a-frozen-base**: each sub-step from
S11c-c1 onward imports the fold `load_model(frozen base ⊔ deltas)` and writes only its **own-rows delta**;
membership is the **bind-closure** (export a row iff a later step binds it); an under-export **guard**
(`research/pde_ledger_v3/scripts/ledger_fold.py`, verified `1be64084`) resolves manifest edges by full symbol
identity. This design is SETTLED and COMMITTED — ⛔ do **not** re-litigate the design itself.

This review is of the **migration** that enacts the design into three committed directives. The migration must
be a faithful, topology-only transcription of the design — ⛔ it must change **no physics**.

## Artifacts (review the working-tree versions against HEAD `2f61627a`)
Diff each against its committed baseline: `git diff 2f61627a -- <path>`.
1. `research/pde_ledger_v3/directives/S11c_decisions.md` — N1 amended (topology + digest sentences).
2. `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` — §1 header, §1a, §7 chain-output + tag-grammar.
3. `research/pde_ledger_v3/directives/S11c_c1_sympy_build_directive.md` — the "Delta export" section, the
   EMIT-vs-EXPORT + §D3 guard block, and `BUILD_INPUT_DIGESTS` (5→6 inputs).

## What you are handed (read all)
- The three edited files above (working tree) and their committed baselines (via `git diff`/`git show`).
- The DESIGN authority: `research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md` (§D1–§D4, §Build plan).
- The fold module + its public API: `research/pde_ledger_v3/scripts/ledger_fold.py` (`load_model`,
  `check_consumer`, `assert_lookups_equal_manifest`, `assert_delta_is_minimal`, `promotion_delta`).
- The grounding measurement: `research/pde_ledger_v3/directives/_measurements/c1_migration_import_keys_scan.py`.
- The frozen base ledger: `research/pde_ledger_v3/scripts/S11c_b_exports.py` (2441 rows, plain git, present).

## Required method — DOCUMENT branch, but it makes a REAL-ARTIFACT claim, so RUN the scan
This is prose (a decision list + a spec + a build directive), so read the DESIGN and the SPEC first, form your
own view of what the migration must say, and only then read the directive. BUT the directive's `IMPORT_KEYS`
manifest is a claim about the real frozen base — so **run** `check_consumer` yourself and verify it, and
**independently re-derive** the bind-set from the spec. A prose "looks right" is discarded; show the command and
its literal stdout for every artifact claim (`CLAUDE.md` rule 2). Save any script you write and its stdout to
named absolute paths and report them.

## The load-bearing questions (each is a way the physics/chain could be wrong)

1. **IMPORT_KEYS COMPLETENESS (under-export — the highest-value check).** Read the c1 spec §1c, §1d, §2, §3a,
   §3b, §5a–§5e, §6 and enumerate **every** upstream object the engine must **bind** (look up in the fold) as a
   construction or control/regression operand. Cross-check against the directive's declared `IMPORT_KEYS`
   (build directive, "Delta export" section). **Report any bound object the spec requires that is MISSING from
   IMPORT_KEYS** — that is an under-export that fails the build. (E.g.: are all nine T-substrate rows there? all
   three density representatives? every `Lambda_*`/`tau_*` knob the closure uses? the branch object `q_out`?)

2. **IMPORT_KEYS MINIMALITY (stale manifest).** The build's bidirectional smoke-test requires observed lookups
   to equal IMPORT_KEYS **exactly** — a declared-but-unused key raises. **Report any IMPORT_KEY the spec does
   NOT require c1 to bind.** In particular adjudicate the five `degenerate_loci_*` rows and `z_by_parity`: does
   the §5c uniform-limit / §5d zero-jet regression of the DtN and the face response actually look them up, or
   are some of them over-declared? The directive already flags this as an open subset question routed to §8 —
   say which keys you believe are genuinely bound vs. speculative.

3. **The face_response F9c handling.** c1 imports the **bare** S11b `face_response`/`face_response_coeffs` as
   §5c/§5d regression operands, AND writes its OWN curved response, which collides → F9c producer-prefix
   `s11c_c1_face_response`/`s11c_c1_face_response_coeffs`, leaving the bare S11b predecessors intact. Verify:
   (a) both bare keys resolve to the S11b flat rows in the fold (run it); (b) the design/spec support the F9c
   split; (c) the directive never overwrites the bare predecessors. A wrong bare/prefixed routing is the exact
   "face_response predecessor trap" the design warns of.

4. **FIDELITY to the design (topology only).** Does each edit faithfully transcribe §D1 (membership by binding,
   not label; missing manifest = hard error), §D2 (own-rows delta; F9 write-time; fold last-wins, no F9
   re-apply; ⛔ no write-time pass-through), and §D3 (the guard = check_consumer + smoke-test + minimal-mode;
   the round-trip does NOT catch under-export)? Report any drift, any place the migration weakens or contradicts
   the design, or any residual "accumulate"/"carry forward"/"export every primary" instruction left standing.

5. **The BUILD_INPUT_DIGESTS decision.** The directive pins **6** inputs: own source, `S11b_exports.py`,
   `S11c_a_exports.py`, `S11c_b_exports.py`, the spec, and `ledger_fold.py`. But `load_model` reads only
   `S11c_b_exports.py` for data. Is pinning `S11b_exports.py`/`S11c_a_exports.py` still correct (byte-provenance
   of the frozen base's ancestry), or should c1 pin only what it reads (own source, `S11c_b_exports.py`, spec,
   `ledger_fold.py`)? Give a reasoned verdict; this was deliberately left for you to adjudicate.

6. **NO PHYSICS CHANGED.** The migration must alter only export/import topology. Confirm no equation, order,
   control, object name, regime, branch, or dimension in the spec §§1–6 or §4 tag list was changed by these
   edits. Any physics change here is a defect. (Diff the spec: only §1 header, §1a, §7 should differ.)

7. **RUN check_consumer against the real base.** Independently run
   `cd research/pde_ledger_v3/scripts && python3 ../directives/_measurements/c1_migration_import_keys_scan.py`
   (or write your own). Confirm: single atomic base (2441 rows, 0 overwrites), all declared roots present, and
   the closure resolves with no AmbiguousSymbolError/ClosureError. Then **ablate to prove the check bites**:
   remove one genuinely-bound root from the candidate list and confirm the closure changes / a downstream edge
   is affected; add a bogus non-existent key and confirm `ManifestError`. Report literal stdout.

## Physics filter
Report a finding only if it catches a way the chain, the physics, or the c1 build could go wrong — an
under-declared bind, an over-declared bind, a physics change smuggled into a topology edit, a design-fidelity
drift, a wrong F9 routing, or a digest/provenance error. Do not report style. A missing IMPORT_KEY and a changed
equation are the two most valuable findings.

## Ablation sandbox
Copy anything you execute-ablate to `/tmp` and ablate the copy; ⛔ never modify the working tree. The scan and
`ledger_fold.py` are read-only; running them does not spawn CAS kernels, so no seat/timeout constraints apply.
