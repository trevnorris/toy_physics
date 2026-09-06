# S11c-c2 export-repair directive — decision gate (2 legs, one pass)

**Artifact:** `directives/S11c_c2_export_repair_directive.md` — a publication-only export-repair build directive
(orchestrator-written) → **two decision legs = Codex `gpt-5.6-sol` xhigh + Grok** (G1), one two-leg pass, fold
once, go (G2 — ⛔ not iterate-to-green). Identical prompt: `directives/_legs/S11c_c2_export_repair_decision_prompt.md`.
Reports: `_measurements/S11c_c2_export_repair_decision_{codex.md,grok.txt}` (Codex raw 1.25 MB in
`ext_logs/c2_export_decision_codex.txt`, hygiene).

**Both legs converged** (no leg-vs-leg disagreement): scope + the three edit sites, R1 membership (S11c-d binds
**both** operators — closed slab operator for its closed spectral equation/resolvent, coupling kernel for
Born/mixing; increment genuinely EMIT-only), R3 hygiene (no consumer binds `display`; ~6 MB / ~5 MB sstr dup on the
two operator rows), R4 guards, and the overstep fence are all **SOUND**. The one blocker is **R2** (representation
+ verification), which I under-specified. I verified the load-bearing findings myself against the code (rule 13):
the export row IS a nested cased Tuple; line 848–851 IS a compact-only srepr roundtrip (separate from a
compact-vs-emitted check); `difference`/`zip`@100 does not check tuple lengths.

**Fold applied (one pass) to the directive:**
- **R2 rewritten** (MUST-FIX, both legs): compact **only the per-case `VALUE` operator leaves**, keep the cased
  `{VALUE,MULTIGRADE,DIMENSION,bindings}` tree intact (a `cases()`/`named(…,'VALUE')` consumer requires it). Allowed
  transforms `factor`/`cancel`/`together`/`collect` staying one evaluable `sp.Basic` of the emitted leaf's type;
  ⛔ **CSE forbidden** (dummy temporaries aren't declarations → incomplete-skeleton `diff`; only shrinks if the
  replacement table is stored = opaque wrapper); ⛔ no `UnevaluatedExpr`/hold/pickle/string-only. **Preserve the
  singular locus** (pole set feeds d's resonances — ⛔ no denominator-canceling transform without a pole-set check).
  **Two separate guards**: KEEP the structural roundtrip (compact ↔ `_restore(srepr(compact))`), ADD a separate
  strict-recursive **emitted↔compact** semantic check (identical case-key sets / arities / mapping keys / `Str`
  labels / matrix shapes; exact metadata equality; leafwise `expand(decoded−emitted)==0`, **Integral-aware** —
  don't let `Integral(0,…)` slip through — with explicit length/shape checks since `zip` doesn't; hard-stop on
  mismatch; run against the generated module).
- **Nits folded:** allow a publication helper adjacent to `publish`; require a **short bounded `display`** (not
  omission) + note it supersedes S9 D3 / bind-closure D4's human-readable clause for the giant rows; `.out`
  identity = tag/emit-key names, not a 499 MB payload byte-diff; the run's `S11c_c2_sympy_guard_evidence.json`
  checkpoint write is an expected harness side-effect (not a builder artifact); the "no S11c-d `IMPORT_KEYS`
  manifest yet" caveat (both operators are the declared prospective binds).

**Disposition:** directive **GATED**, folded once, ready for the builder. Next: astra build (detached) →
2 re-review legs (fresh Claude + Grok) on the diff → adjudicate → commit.
