# Decision 17 — Trust-apparatus trim: keep the anti-fake core + frozen-data pinning, drop script content-pinning + the reseal ceremony

**Decided by the user, 2026-07-17**, at the post-B2 scope review the user pre-registered during the B2 production saga. Status: OPERATIVE for U1 Phase C and all subsequent builds. **Sealed B2 artifacts are untouched** — contract `8b8ee113…` remains valid as sealed (runner windows stay 1200s; reverting would invalidate it); this decision governs *new* contracts only.

## The decision
For Phase C onward, build contracts:

**KEEP (the anti-fake core — what actually caught the fakes):**
1. **Fresh-agent adversarial recompute** on every computed verdict (imports the engine, independently recomputes the decisive quantity; the RIG lesson).
2. **Per-tooth mutation campaigns** — every able-to-fail check gets a mutation that dies AT its own assert.
3. **Dual-engine (SymPy + independent Mathematica)** wherever Mathematica *can* independently verify.
4. **Trace/network containment** (network=0; sandboxed runs) for production compute.
5. **Hash-pinning of frozen physics data** — the derived/UNRESOLVED disposition sets and any frozen upstream inputs a build consumes are pinned by content hash (cheap, load-bearing: B2's disposition sha `d9434117…` staying byte-identical across reseals is what let approval compose).

**DROP:**
6. **Script content-pinning in an Obs-style manifest.** Script identity = the git commit anchor (`STARTUP_CONTRACT_COMMIT`-style: the orchestrator-supplied commit whose tree contains the scripts). A script fix = commit + re-run.
7. **The reseal ceremony** (full `--stage0` rebuild to re-bind changed script bytes). No stage-0 rebuild is triggered by a script edit; rebuilds happen only when *frozen physics data* actually changes.

## Evidence backing it (the B2 arc's ledger)
- **What caught fakery across the whole U1 arc (5 RIGs + the last latent bug):** fresh-agent recompute (RIG #2's YAML-token verdict; B1's scalar⊗δ ansatz; B2's 183/183 G9 witness verification), per-tooth mutations (15/15 killed-at-own-assert; the stale completion-gate tooth id), dual-engine independence (the on-disk engine inconsistency; forcing the thrice-faked χ–u Hessian to become real), trace containment (network=0; the containment trip). None depended on the manifest.
- **What the manifest + reseal ceremony delivered:** it never caught a fake. Six sealed-run attempts; five deaths caused or amplified by the seal itself — the `/proc/self/fd` path-bug fix *rejected by the manifest*; the broken re-seal path that clobbered two pinned files; reseals #2/#3/#4 (~50 min `--stage0` each) to bind one-line reviewed script fixes. File identity is exactly what git already provides.
- **User steer (recorded 2026-07-17, pre-decision):** the immutability/consistency half of the trust apparatus largely duplicates git at heavy friction cost; its real value was catching in-process fakery, not immutability.

## Scope + guardrails
- The trim removes *ceremony*, not *verification*. Every computed verdict still requires the fresh-agent recompute leg; honest-UNRESOLVED still requires the v48-style computed unavailability witness + derivability challenge; the four-leg (arbiter re-run / fidelity / adversarial recompute / external review) stands.
- Frozen-data pinning is one-way: a build may *consume* pinned upstream data; nothing in the trimmed contract may *rewrite* it. Any change to frozen physics data is a first-class amendment (reviewed), not a re-seal.
- If a future build shows in-process script tampering that git anchoring fails to catch (i.e., the dropped mechanism would have caught a real fake), that is grounds to revisit — bring to the user, do not silently re-add ceremony.

## Consequences / follow-ups
- The U1 Phase C directive is authored against this trimmed contract (leaner stage-0; no Obs script manifest; git-anchored scripts + hash-pinned frozen inputs).
- `docs/em_analog_next_phase_handoff.md` step (1) is discharged; the B2 anchors recorded there remain the authoritative pins for the sealed B2 artifacts.
