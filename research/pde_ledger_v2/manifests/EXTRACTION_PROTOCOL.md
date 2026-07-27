# Stage-manifest extraction protocol (v2.1) — REUSABLE

Extract one machine-readable v2.1 `StageManifest` for one physics
derivation-ledger stage. Project root: `/var/projects/toy_physics`.
The manifest is a derived artifact, never a second source of truth. Every
machine-readable fact must be traceable to an audit engine or dated record.

## Read first

- `research/pde_ledger_v2/manifests/stage_manifest_schema_v2.json`
- `research/pde_ledger_v2/manifests/MANIFEST_README.md`
- `research/pde_ledger_v2/manifests/DIM_ORDER_DECISION.md`

Then read every source that exists for the stage:

- stage note: `research/pde_ledger_v2/notes/stages/ledger_stage<NNN>_*.md`;
- paper card: `research/pde_ledger_v2/paper/stages/stage_<NNN>.tex`;
- SymPy audit: `research/pde_ledger_v2/scripts/ledger_stage<NNN>_*sympy*.py`;
- Mathematica audit:
  `research/pde_ledger_v2/mathematica/ledger_stage<NNN>_*.wl`;
- ablation runner, if present; and
- `research/pde_ledger_v2/notes/parameter_register.md`.

Do not infer a fact that the sources do not supply.

## Runtime and tracked-evidence constraints

- Manifest extraction executes no Mathematica. It cites the SymPy audit and a
  digest-pinned Mathematica companion/output, so the ≤2-concurrent-Mathematica-
  seat cap does **not** bind extraction itself. The cap applies only when
  Mathematica audits are actually re-run.
- Stage-audit `.out` evidence is intentionally tracked: a `.gitignore`
  negation un-ignores `research/pde_ledger_v2/mathematica/out/*.out`. The
  repository-wide `*.out` ignore rule targets LaTeX artifacts, not these saved
  audit outputs.

## Extraction rules

1. **Named dimensions and source order.** For every symbol, record the audit
   script's actual tuple order in `dim_source_order`, the script path+locus in
   `dim_source`, and the raw positional exponent vector exactly as printed in
   `dim_source_tuple`. Transpose that tuple into a named `dim` object. Omit zero
   keys. The checker independently recovers the order from the source
   `Dim(...)` signature/dataclass fields and order docstring; never infer it
   from the manifest enum. The named map is authoritative. Exact rational
   exponents are allowed, so `sqrt` and other rational powers remain exact.

2. **Quantity identity and ownership.** Give every symbol a global stable
   `quantity_id`, its
   stage-local SymPy `parse_alias`, and the exact `definition_ref` that defines
   it (`here/claim_id` or `stageNNN/claim_id`). If the same printed name denotes
   different quantities, emit an explicit shared `alias` declaration listing
   all distinct quantity ids. A quantity has exactly one owner: the first stage
   whose claim actually binds it. Only a relation with the symbol on its lhs,
   a spectrum with matching `kernel_symbol`, or an operator identity with
   matching `acts_on` is a definition. Every later use points to that owner and
   declares the matching consume; never mint a local wrapper definition for an
   upstream quantity. Callable functions/operators need domain/codomain
   dimension signatures.

3. **Typed payloads.** Every claim must use one non-empty payload kind:
   `relation`, `operator_identity`, `spectrum`, `adjudication`,
   `record_range`, `set_cardinality`, or `token`. Do not use `note` as the
   machine-readable claim. Use the minimal `holds_within` conditions AST for
   conjunctions, alternatives, negation, XOR branches, and domains. Claims
   asserting a dimensional result carry `expected_dim`. Respect the schema's
   claim-kind/payload-kind matrix; do not disguise spectra, ranges, sets, or
   adjudications as relations merely to obtain a CAS mode.

4. **Exact expressions.** Emit canonical SymPy strings with implicit
   multiplication off and exact rationals only. Floats are errors. Use the
   README whitelist and declared signatures. LaTeX is a mirror, never the
   parsed source.

5. **Import completeness.** Inventory every quantity referenced by every claim
   payload, operator `acts_on`, spectrum operator/kernel fields, recursive
   condition/domain, consume `as_consumed`, structured substitution, and
   consume specialization domain. For each nonlocal `definition_ref`, emit a matching
   `consumes` entry resolving to that defining claim. A used nonlocal quantity
   without one is `UNDECLARED_IMPORT`; do not rely on a prose citation or
   `note`. If the producer claim id is not yet available, record a provisional
   ref in the manifest extraction report.

6. **Structured consumption.** Record the producer `ref`, typed
   `as_consumed` payload (or `as_consumed_dim` for `dim_equal`), exact `check`
   mode, and `substitutions`. Each substitution is `{lhs,rhs,
   backed_by}` and `backed_by` must resolve to an exported convention claim.
   Never encode substitution instructions in prose. `specialization` requires
   the explicit flag and a structured domain. Use the producer-kind mode:
   `spectrum_match`, `range_match`, `adjudication_match`, or `set_match` for
   those non-relation kinds. `cas_equivalence` on them is unsupported.
   `token_match` is only for real status/verdict string tokens, never a
   quantity. Cite opaque quantities with `opaque_quantity_match`, the producer
   `quantity_id`, and its export `source_digest`.

7. **Exports and C7.** **Export-complete (mandatory, user-approved 2026-07-24):**
   export EVERY operative claim AND every ownership (`declare_*`) claim of the
   stage; only retired, superseded, or departed claims stay unexported. Exports
   are NOT a curated highlight reel — downstream stages legitimately cite
   intermediate lemmas and quantity definitions, and a stage that under-exports
   forces its consumers into opaque citations or re-derivation. (The 030/031/032
   pilot proved the failure mode: 031 did not export `S_gg`, so 032 fell back to
   an opaque `C_V`.) Export-list size is free; the `NON_EXPORTED_CLAIM` guard
   still bites where it matters, because citing a RETIRED or departed claim
   remains a real error. Pin each export to its source evidence digest.
   `c7_binding` and `c7_expect` are
   optional until a real mutation exists. When present, the export binding
   carries its producing symbol/knob id, mutation env, executable mutation
   command, and exported semantic facet; the consume expectation carries the
   corresponding injection point, facet used, and expected first failing
   predicate. That predicate must already exist in the consumer's
   `verification.teeth`. A causal slice remains PARTIAL until C7 has run on all
   of its edges. Every `mutation_command` must be an honest mutator whose
   behavior changes in response to `C7_FACET` (or the declared
   `mutation_env`); a command that ignores the facet and merely prints the
   expected teeth silently defeats C7. The checker trusts command output and
   cannot prove mutator faithfulness. A future `C7_MUTATOR_INERT` canary may
   warn when a declared sentinel produces output identical to an unmutated
   baseline, but that warning is advisory because a legitimately robust export
   can be facet-insensitive.

8. **Lifecycle census.** Model each knob as an event with stable `knob_id`,
   register row, origin claim, effective stage, low/high count effect, category,
   pending state, and evidence. Pending debt remains counted. Retired lifecycle
   events use the register-grounded endpoint effects; departures are not knob
   events and always contribute `{low:0,high:0}`. `discharged` requires a claim
   ref to DERIVED-and-executed discharge evidence.

9. **Evidence.** Every symbol, claim, knob, departure, and verification tooth
   gets `{source_path,locus,source_digest,engine,method}`. Use
   `method: prose_only` for prose claims not verified by a script, such as the
   stage030 `A_eff` relation; such evidence must use `engine: prose`, never
   SymPy/Mathematica. Recompute every evidence digest from `source_path` before
   writing. Persist the extraction report and source digests in `extraction`.

10. **Genesis.** POSTULATED, CONV, and CALIBRATED claims require genesis.
    Set `origin: independent` only when a dated record proves it and include
    genesis evidence with path+commit, date, span, and the later claim it
    predates. `coordinated`/`target_matched` require refs. Otherwise use
    `unknown`.

11. **Verification teeth.** Copy actual predicate tokens, mutations, defended
    local claim ids, and per-tooth evidence from the scripts/ablation runner.
    Do not invent a tooth for an unexecuted assertion.

12. **Ranges and adjudications.** For record ranges, encode every convention
    axis as choices with exact low/high deltas. Component endpoint sums plus the
    Cartesian choice deltas deterministically produce the reported low/high;
    `spread=high-low`. For adjudications, emit `bucket_counts` and/or axis
    cardinalities when the source supplies them; counts must sum and axis sizes
    must multiply to `domain_cardinality`.

## Mechanical completion checks

- Parse the JSON and validate it with the draft-2020-12 v2.1 schema.
- Confirm every payload is typed and every object has no unknown properties.
- Parse all SymPy expressions with the local `parse_alias` table; reject floats.
- Check local claim ids, export ids, definition refs, and tooth claim ids.
- Run import-completeness and verify every substitution's exported convention.
- Recover each dimension order from the script, compare it to
  `dim_source_order`, compare the printed tuple to `dim_source_tuple`, and
  transpose it into the named dimension.
- Confirm per-item evidence and extraction source digests are present.
- Do not run `$RT exec-sympy`; if execution is useful, run the stage script
  directly.

## Return

Report:

`MANIFEST_WRITTEN: stage<NNN>.json` plus counts for symbols, claims, consumes,
exports, knobs, departures, and teeth.

Then list `UNGROUNDED`, `PROVISIONAL CONSUMES`, `PROSE/SCRIPT MISMATCH`, and
`JUDGMENT CALLS`. Keep the report compact and persist the same findings in the
manifest's `extraction.report`.
