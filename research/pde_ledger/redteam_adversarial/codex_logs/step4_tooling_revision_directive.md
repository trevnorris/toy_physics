# Step 4 tooling REVISION directive — fix dedup/family via a target-normalization layer (Codex applies, Claude reviews)

Context: the Phase B tooling you built (session `019eb706`) is code-correct and passed review, BUT the two read-only proposals it generated are unsound, confirmed by a clean-agent data review AND your own read-only consult (record: `codex_logs/step4_dedup_family_consult_codex.txt`; you AGREED Q1–Q4, refined Q2/Q3/Q5). This directive implements the agreed fix. The agreed root cause: Phase A `parameter_names` are polluted (whole-sentence symbol captures from the fit-signal scan at `core.py:642`) and the literal-value extractor misses common forms (`sp.Float("4.651...")`, JSON colon values), so dedup over-merges (value-compatible fallback over a parameter-family union, transitively — e.g. stage073 collapses 7 distinct parameters across lines 6–94; stage247 merges 13 across lines 39–497) and family-build over-splits (param-set + each literal value as the key — chi_q→18 families; the Σ0 …867/…876 drift fragments across separate families).

This is a TOOLING revision only. Do NOT apply the real maps, do NOT run agents, do NOT adjudicate, do NOT edit `paper/`/`notes/`/`scripts/`/`graph/`. Do NOT commit; leave the tree dirty.

## Binding constraints (identical to the prior builds)

- Write set: `.claude/skills/adversarial-audit/**` and `research/pde_ledger/redteam_adversarial/**`; additive-only `adversarial:` config keys if genuinely needed. `docs/adversarial_audit_directive.md` FROZEN. Zero edits to `.claude/skills/redteam-audit/`.
- Every campaign-state-mutating write through the existing flock; read-only proposal generators take NO lock and mutate nothing.
- `timeout 600` everywhere; exit 124 = reformulate, never raise the cap. No OS temp files agents read (`tmp_prompts/`). YAML/markdown for agent-facing files.
- The Phase A scanners, `phase-a-ingest`, `phase-c-render`, `phase-b-ingest`, and `benchmark-ingest` validation rules built last session are CORRECT — do not regress them. Keep the manifest at 922 `scanned` candidates; this revision changes how candidates are GROUPED, not the candidate set.

## Item 1 (root fix) — deterministic target-normalization layer

New read-only command `target-resolve [--out PATH]` (default `redteam_adversarial/provenance/_target_layer.yaml`). Read-only on the manifest; writes only the sidecar. For each `scanned` candidate it resolves:

- `primary_target_parameter`: the ONE conceptual parameter the candidate is about. Resolution priority, deterministic and documented: (a) parse from `candidate_key` (the keys are human-meaningful: `stage073_eta_pinned_37`→`eta`, `stage_105_chi_q`→`chi_Q`, `stage168_sigma0_num_digit_variant`→`Sigma0_num`); (b) fall back to `parameter_names` ONLY when the key is uninformative; (c) apply the concept-variant alias table (Item 4) to map variants to one canonical concept id.
- `target_values`: ALL literal values for that parameter, extracted from the citation excerpt(s) AND the cited source line, with an improved extractor that catches at least: `sp.Float("4.651033550168876", 30)` and `sp.Rational`/`sp.Integer`; JSON `"key": <number>` colon forms; chained `= 0.05 = 1/20` (record both, recognize equality where exactly computable); closed forms like `Sqrt[4107 - 100*Pi^2]/(10*Pi)`; plain decimals and `a/b` fractions. Normalize for comparison where exactly decidable (0.05 == 1/20). If no value can be confidently resolved, leave empty.
- `target_confidence: high | low` — `low` when the key is uninformative, multiple plausible targets exist, or the value could not be extracted. `resolution_basis`: short string (which rule fired, which source the value came from).

The sidecar is the single source of clean parameter identity for everything downstream (dedup, family, and `phase-b-build`). It must be deterministic (same manifest → same sidecar modulo timestamp) and reviewable (one record per candidate with enough basis to audit the call).

## Item 2 — rebuild `dedup-propose` to consume the target layer, strict and value-anchored

- Load `_target_layer.yaml` (generate it first if absent, or require it — your choice, document it).
- Alias two candidates ONLY when ALL hold: same `anchor_stage`; same `primary_target_parameter` (post concept-alias); same normalized `target_value`, **non-empty and confidently resolved on BOTH sides**; citation is identical `(path,line)` OR same normalized path within ≤3 lines OR a pure path-prefix variant of the same line (e.g. `paper/...` vs `research/pde_ledger/paper/...`).
- DROP `records_value_compatible` and the no-conflicting-value fallback entirely. NO merging across a `parameter_family` union. Transitive grouping is allowed ONLY within an identical `(primary_target_parameter, target_value)` equivalence class — never across differing values. If the value is unresolved/low-confidence on either side, do NOT merge.
- Keep the `ambiguous` section (same stage+parameter, differing value/site). **Invariant: a candidate is in exactly one of {aliased, ambiguous-separate, standalone canonical} — never both.** Add an assertion that the alias∩ambiguous id-overlap is ZERO and fail loudly if not (the prior proposal had 153 overlaps).
- Expect the canonical count to RISE substantially (toward ~800+). That is the correct direction — over-merge hides fits, under-merge only costs redundant Phase C audits. Report actuals: alias groups, total aliases, canonical count, ambiguous groups, overlap (must be 0).

## Item 3 — rebuild `family-build` to key on conceptual-parameter identity

- Family key = `primary_target_parameter` (concept-normalized) ALONE — NOT param-set, NOT param+value.
- Each member carries its `target_values` as per-member attributes. When a family's members hold divergent values, emit a `value_divergence` finding on that family listing the distinct values + their stages/candidate ids (this is exactly how the Σ0 …867 vs …876 transposition must surface — both must be in ONE family `Sigma0_can`/`Sigma0_num` concept with a `value_divergence` finding, NOT in separate families).
- One primary family per candidate. Allow >1 family ONLY for genuinely multi-target block candidates (those whose `candidate_key`/reason legitimately covers multiple parameters); flag those (`multi_target: true`).
- Singletons in the long tail are acceptable.
- Every canonical maps to ≥1 family; `unmapped_canonical_count` must be 0.

## Item 4 — concept-variant alias table

A small, reviewable, deterministic alias table (`redteam_adversarial/provenance/_concept_aliases.yaml`, or a config section — your choice) mapping symbol variants to one canonical concept id. Seed it from the D3 chains and the observed variants, at minimum: `Sigma_0_can`/`Sigma0_can`/`Sigma0_num` → one concept; `chi_q`/`chi_Q` → `chi_Q`; the 54/5 quadrupole carriers; the 37/20 aspect-ratio carriers; `m0_hat`/`mhat0`; `T_can`/`S_can` transport. Used by `target-resolve` (Item 1). Keep it data, not hardcoded logic, so it is reviewable and extensible.

## Item 5 — seed the D3 headline chains as named families + coverage warnings

- `family-build` seeds explicit NAMED families for the D3 chains (membership populated mechanically via concept matching against the target layer):
  - `chain_quad_54_5` (54/5 + gamma_GR quadrupole target), `chain_aspect_37_20` (37/20 → F1 incl. the 4107 closed form), `chain_chi_Q_norm` (chi_Q/m0_hat/N_Q), `chain_sigma0_transport` (Sigma0_can/T_can incl. 867/876), `chain_barrier_222_224` (incl. V_known), `chain_calibration_245_253` (incl. lambda_L + CODATA mass ratio), `chain_wall_action` (wall-action coeffs + m̂/P0).
- For each seeded chain, the expected stage span is known (from the D3 list / fit_insertion_points). Emit a HARD `coverage_warning` listing any expected stage/candidate the mechanical rule failed to place in its seeded chain, so the orchestrator can route those to agent review. Do NOT silently absorb misses.

## Item 6 — `phase-b-build` consumes the primary target

`phase-b-build` currently iterates `entry["parameter_names"]` (the polluted union), so a single candidate would spawn ~28 provenance files. Change it to build the genealogy for the resolved `primary_target_parameter` (from the target layer), plus any additional targets only for flagged `multi_target` candidates. Preserve the dry-run behavior and the `synthesis_status` pending/complete state split from last session. This keeps Phase B genealogies keyed to the real target, not the symbol soup.

## Verification contract (iterate until everything exits 0)

1. Build Items 1–6.
2. Run the read-only generators on the REAL manifest, in order: `target-resolve`, then `dedup-propose`, then `family-build`. Leave `_target_layer.yaml`, `_concept_aliases.yaml`, `_dedup_proposal.yaml`, `_family_map.yaml` in the tree (overwriting the prior unsound proposals). Report:
   - target-resolve: high vs low confidence counts; how many candidates had a value successfully extracted (sanity-check that `Sigma0_num = sp.Float("4.651033550168876")` now extracts the value, and that `fit_stage_105_chi_q` resolves to ONE primary target).
   - dedup: alias groups, total aliases, canonical count, ambiguous groups, alias∩ambiguous overlap (MUST be 0). Spot-confirm stage073 is NO LONGER one merged group (its 7 distinct parameters become separate canonicals or only same-value aliases).
   - family: family count, singleton count, unmapped (MUST be 0); for EACH seeded D3 chain, the member count and whether its `value_divergence`/coverage findings populated — explicitly confirm the Σ0 …867 and …876 carriers are in ONE family with a `value_divergence` finding.
3. Fixture-verify `apply-alias-map` still applies a (small, strict) map correctly using dry-run fixtures; purge. Confirm `dry-run --stages 003 104 105` still satisfies every documented dry-run expectation; purge.
4. Do NOT apply the real maps; do NOT run agents; do NOT commit.

Report: per-item file:line changes; the design choices (key-parse rules, value-extractor coverage, concept-alias entries, seeded-chain membership rule); the real proposal counts from step 2 with the three explicit spot-confirmations (stage073 unmerged, chi_q one target, Σ0 867/876 co-familied with divergence finding); the fixture/dry-run/purge transcript; anything unsatisfied and why.
