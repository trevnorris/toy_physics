# Step 4C tooling directive — bulk Phase B build + systemic YAML perf (Codex applies, Claude reviews)

Context: Step 4C must create the *pending* Phase B mechanical bundles for ALL 915 non-duplicate canonical candidates so the genealogy agents can synthesize them. The existing per-candidate `phase-b-build` costs ~28.7s/call → ~7.3h for 915. Profiling (orchestrator, one call) shows the cost is entirely per-candidate-INVARIANT setup repeated each process:

```
load_manifest          11.9s   (MANIFEST.yaml is 2.7 MB; yaml.safe_load pure-Python)
target_params (1st)     2.5s   (loads the 599 KB target layer)
graph_index (1st)       2.1s   (loads the 500 KB atlas graph)
render_fit_file         2.0s
render_batches          4.7s
save_manifest          ~5s
source_evidence(1 cand) 0.06s  <-- the ONLY genuinely per-candidate cost
```

Two changes fix this. Both are pure tooling under `.claude/skills/adversarial-audit/**`. Do NOT touch `paper/`/`notes/`/`scripts/`/`graph/`; `docs/adversarial_audit_directive.md` is FROZEN; zero edits to `.claude/skills/redteam-audit/`. Do NOT commit; leave the tree dirty.

## Binding constraints (identical to prior builds)
- Write set: `.claude/skills/adversarial-audit/**` only (plus the verification scratch you purge). Additive-only `adversarial:` config keys if genuinely needed.
- Every campaign-state-mutating manifest write goes through the existing flock (`require_manifest_lock` / `_manifest_locked`). Read-only/dry paths take NO lock.
- `timeout 600` on any internal script; exit 124 = reformulate, never raise the cap. The bulk build is ONE Python process and must finish well under that.
- YAML/markdown for all artifacts. No OS temp files an agent would read.

## Change 1 (systemic) — use libyaml C loader/dumper when available

`yaml.CSafeLoader` and `yaml.CSafeDumper` are present in this environment (orchestrator-verified: C loader reads the 2.7 MB manifest in 1.45s vs 11.9s). Apply it everywhere `core.py` loads/dumps YAML, behind a capability check with a pure-Python fallback:

- `load_yaml` (core.py:229): use `yaml.CSafeLoader` when present, else `yaml.SafeLoader`. (Currently `yaml.safe_load`.)
- `write_yaml` (core.py:240): add a C-based no-alias dumper and use it when present, else the existing `NoAliasDumper`:
  ```python
  if hasattr(yaml, "CSafeDumper"):
      class NoAliasCDumper(yaml.CSafeDumper):
          def ignore_aliases(self, data): return True
      _DUMPER = NoAliasCDumper
  else:
      _DUMPER = NoAliasDumper
  ```
  Keep the exact same `yaml.dump(..., default_flow_style=False, sort_keys=False, width=120, allow_unicode=True)` call kwargs.
- Audit core.py for any other `yaml.safe_load` / `yaml.load` / `yaml.dump` / `yaml.safe_dump` sites (target layer, family map, dedup, benchmarks, synthesis ingest, graph index) and route them through `load_yaml`/`write_yaml` or apply the same C-loader/dumper selection so the speed-up is uniform.

**Content-neutrality is REQUIRED and must be proven.** The committed artifacts must not change in content. Acceptance test (run it, paste the transcript):
1. Copy `redteam_adversarial/MANIFEST.yaml` and `redteam_adversarial/provenance/_family_map.yaml` to scratch paths.
2. Load each with the NEW `load_yaml`, re-dump with the NEW `write_yaml` to a second scratch path, and ALSO re-dump with the OLD pure-Python `NoAliasDumper` to a third scratch path.
3. Assert (a) the parsed data structures are equal (`==`) across old/new loaders, and (b) the C-dumped text equals the Python-dumped text (a `diff` that is empty, or differs only in trailing whitespace you then normalize). If the C dumper emits any YAML anchor/alias (`&`/`*`) the test FAILS.
4. Purge all scratch files. The real MANIFEST.yaml / _family_map.yaml must be byte-unchanged by this test (`git status` clean for them after the test).

If the C dumper cannot be made textually identical to `NoAliasDumper`, STOP and report — do not silently change artifact formatting.

## Change 2 — bulk command `phase-b-build-all`

Add a new subcommand `phase-b-build-all [--ids-file PATH] [--limit N] [--force] [--dry]`:

- **Refactor, don't fork.** Extract the per-candidate body of `build_provenance` (core.py:3281) into a helper, e.g. `build_provenance_for_entry(env, manifest, candidate_id, *, graph_index, target_layer, save=False, render=False)` that builds the pending bundle(s) and mutates the in-memory `manifest`/`entry` but does NOT save/render. The existing single-candidate `phase-b-build` calls it with `save=True, render=True` (so its behavior is unchanged). The bulk command calls it in a loop with `save=False, render=False`, then renders fit-file + batches ONCE and saves the manifest ONCE under the flock at the end. There must be exactly one copy of the bundle-building logic.
- **Load once:** manifest, `graph_source_index`, and the target layer are loaded a single time before the loop and reused for every candidate.
- **Target set:** all candidates with `status == scanned`, no `duplicate_of`, that do not already have a complete provenance slice. Idempotent: skip a candidate that already has its pending bundle(s) written (paths exist) unless `--force`. `--ids-file` (one candidate id per line, `#` comments allowed) overrides the auto set with an explicit list. `--limit N` caps the number built (for piloting). `--dry` computes and reports the target set and per-candidate target-parameter count and writes NOTHING (no bundles, no manifest save).
- **Report:** target-set size, built, skipped-already-built, skipped-not-eligible, total bundles written, distribution of bundles-per-candidate (to surface multi_target), and any per-candidate errors (collect and continue; do not abort the whole run on one bad candidate — report it). On a real run also print final manifest counts.
- Real candidates stay `status: scanned` with `phase_b_status: synthesis_pending` (matching the existing single-candidate behavior); do NOT transition them. Do NOT run agents, do NOT adjudicate, do NOT touch benchmarks beyond what `build_provenance` already does (the dry-run placeholder maintenance).

## Verification contract (iterate until everything exits 0)
1. Build Change 1 + Change 2.
2. Run the Change-1 content-neutrality acceptance test above; paste the transcript; confirm MANIFEST.yaml / _family_map.yaml are git-clean afterward.
3. Exercise the bulk command WITHOUT mutating the real campaign state: run `phase-b-build-all --dry` and report the target-set size (expect ~915) and the bundles-per-candidate distribution. Then prove a real build works on a TINY scope into a THROWAWAY copy of the artifact tree (or run `--limit 2` then immediately `git checkout`/`rm` the 2 bundles + restore MANIFEST/BATCHES/fit_insertion_points so the real tree is left clean) — your choice, but the committed state must be left with ONLY the `.claude/skills/adversarial-audit/**` code changes dirty and NO pending bundles and NO manifest mutation. The ORCHESTRATOR will run the full real `phase-b-build-all` after reviewing your code.
4. Confirm `dry-run --stages 003 104 105` still satisfies every documented dry-run expectation; purge.
5. Time the bulk path: report the wall-clock of `--dry` and of the `--limit 2` real micro-run so we can extrapolate the full-run cost.

Report: per-change file:line edits; the refactor shape (the shared helper signature); the content-neutrality transcript; the `--dry` target-set + distribution; the micro-run timing + the clean-tree confirmation; anything unsatisfied and why. Leave the tree dirty with only the code changes.
