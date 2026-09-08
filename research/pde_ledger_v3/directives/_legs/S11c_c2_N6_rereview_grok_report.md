I'll review the revised S11c-c2 N6 ablation-harness directive against the engines and the exact PRE→POST diff, starting with the artifact, the four live engines, and the governing compact-PIT standard.The directive and compact-PIT standard are in hand. Next I’ll check the PRE→POST diff and the four engines for coverage, budget-patch sites, and any knife-bearing hunks.The diff looks confined to contract and budget text. I’ll now read the engine PIT/SHA emitters and the WL budget sites against the coverage claims.I have the engine-side evidence. One more mechanical check on tag coverage and the case-gated junk path, then I’ll write the verdicts.# S11c-c2 N6 ablation-harness directive — focused re-review

**Artifact:** `research/pde_ledger_v3/directives/S11c_c2_N6_ablation_harness_directive.md`  
**Diff:** `directives/_legs/S11c_c2_N6_rereview.diff` (7 hunks; 21 minus / 111 plus)  
**Method:** read directive + PRE→POST diff + four engines; `grep`/line counts only. No harness execution.

---

## Q1 — Frozen knives untouched — **SOUND**

The diff touches only: global clauses 1/4/5/7/8; the DEAD-print sentence; the deliverable-transcript sentence; Harness 1’s run envelope (budget patches); one diagnostic wording change (“variant triples” → “compact fingerprint/digest variant triples”); CHANGE-LOG bullets.

No hunk rewrites a knife **scope/site, old fragment, FORM replacement, NO-OP/IDENTITY, COEFFICIENT ×2, or DEAD set** for any of the ten knives. All listed old fragments are still present in the POST file. Knife-adjacent wording is output-contract only.

---

## Q2 — Coverage is real — **finding**

**SymPy arithmetic objects: covered by `n.pit`.**  
`pit` at `diagnostic_sympy.py:668` emits, per `(name, probe)`, `columns` plus `nonzero_modular_numerator` (`:769–770`). Covariance builds `objects = {('N6COV_'+name, None): table ...}` and calls `n.pit` (`covariance_sympy.py:250–252`) for `R_COV`, `SOURCE_ACTUAL`, `SOURCE_PREDICTED`, `R_COV_INCREMENT`, `R_COV_CONTROL_DELTA` (and extras). Reconcile `add(...)` then `n.pit` (`reconcile_sympy.py:246–276`) covers the certified/DEAD `N6RC_*` tables. Diagnostic `add(...)` then `pit(objects, ...)` (`:817–899`) covers `REP_INVARIANCE_*`, both-route `N6_SLOT_GUARD_RESIDUAL` / `N6_CLOSURE_GUARD_RESIDUAL`, both control groups, and the DEAD `N6_MU_RECONSTRUCTION_*` / `N6_MU_AMPLITUDE` / `N6_FACE_VELOCITY` objects.

**`n.sha` exists** at `diagnostic_sympy.py:91` (`hashlib.sha256(sp.srepr(c.cas(e)))`). Covariance already uses it on `Phi` inside `PROVENANCE.operand_sha256` (`:238–239`) and on census pieces (`:106,114`), but **does not** emit a digest of the `FROZEN_PHI` / `PHI_DOMAIN_CENSUS` payloads themselves (`:105–124`). Worker-side `n.sha(obj)` on those dicts is a real fallback (`cas` handles `Mapping`).

**WL arithmetic objects: PIT-covered.** `addNumeric` then `probeCase` (`audit.wl:875–900, 931, 998`) for every named primary/DEAD arithmetic tag. Per-object `PROBE_NUMERATORS` is at `:808`. `FROZEN_PHI` / `PHI_DOMAIN_CENSUS` are `putMeta` only (`:850, 975–976`), not `addNumeric` — matching the claimed non-PIT pair.

**WL digest source is not in the engine.** Directive clause 1 / CHANGE-LOG (`:21–26, :559–566`) require `emitted_object_sha256` from `RUN_PROVENANCE`. Engine counts: `emitted_object_sha256` = 0, `RUN_PROVENANCE` = 0. Existing SHA-256s are stage fingerprints in `RC:PROVENANCE` (`:958–959`) and a few other hashes — **not** per certified/DEAD object. `n.sha` is a SymPy function; it is not a WL hasher.

**Uncovered as specified:** every WL certified/DEAD object has no engine-emitted per-object digest. The two non-PIT objects (`WL_S11CC2_N6COV_FROZEN_PHI`, `WL_S11CC2_N6COV_PHI_DOMAIN_CENSUS`) have neither PIT nor digest.

---

## Q3 — Faithful witness — **finding** (WL; SymPy only if hash-before-projection keeps sample/symbolic content)

The printed PIT object is a **support bitmap / nonzero tally**, not GF(p) values:

- SymPy: `nonzero_modular_numerator` = per-column `any(num != 0)` (`diagnostic_sympy.py:770`); `numerator_denominator` is forbidden in the transcript.
- WL: fingerprint is `{valid, rejected, nonzero_count, circuit_leaves}` derived from `PROBE_NUMERATORS` (directive `:19–21`); the SparseArray values are not kept.

A FORM (or any COEFFICIENT ×2) that **does not change column support** leaves that fingerprint still. Measured WL fact: `R_N6` is already live (18 stored-nonzeros in the rebuild record). Doubling or contaminating a fully supported object need not move `nonzero_count`. COEFFICIENT ×2 over these primes never turns a nonzero column zero.

The backup witness is the digest of the **pre-projection payload**. That works for SymPy **if** the builder hashes the emitted pit dict (including `numerator_denominator`) or the DAG, using `n.sha` / SHA-256 **before** dropping the matrix. It fails if the builder hashes the compact fingerprint (digest becomes a restatement of the bitmap).

`n.sha` on in-memory Node DAGs is not well-defined (`cas` → `sp.sympify` on a dataclass). The spool JSON is the implementable target.

**Invisible bite:** a genuine FORM on a generically-nonzero WL object, with no engine digest, is invisible on the specified compact record. That is the load-bearing hole.

**Spurious movement:** NO-OP/IDENTITY with pinned SymPy seeds (`--seed 110603/110602`) and the same WL restricted `caseOrdinal` should not move a DAG digest. Hashing the **sample matrix** can still move DEAD digests if a FORM/COEFFICIENT changes joint singular rejections and therefore the accepted points (`audit.wl:774–791` jointly evaluates all `flatNodes`). Seed instability across variants is not introduced by the ordinal shift (see Q4d). COEFFICIENT ×2 **should** move a real object digest; that is the companion, not a false FORM.

---

## Q4 — WL budget fix — **finding** (K_junk rendered inert)

**(a) SOUND.** Unique driver `Do` containing `buildCase[case]` at `audit.wl:1011`; unique adaptive assignment `drawCount = If[...]` at `:760`. Second `Tuples[...]` is SETUP metadata at `:582–583`, correctly excluded. No `Environment` / `$CommandLine` / `$ScriptCommandLine`. Module initializer `drawCount = 8` at `:726` is overwritten by `:760`; replacing `:760` is the right site.

**(b) NON-KNIFE for the retained case’s construction — with one harness-level exception.** The iterator only changes the argument to unchanged `buildCase`. `drawCount` is assigned after `numericObjects` is compiled (`:729–760`). Prime list (`:733`), `branchCells`, and knife sites (`:17, :20, :892`) are untouched.

The exception is **K_junk’s case gate**, which the four-case driver used to observe:

```21:21:research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl
actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"};
```

```845:847:research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl
  activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0];
  mat = materialAmplitude[..., activeJunk]; muM = mat["MU"];
```

Junk amplitude is `el[pulled] + jk junkMu eW` at `:247`. Pinning the driver to `{"LAB_HELD","RHOBR_CONSTANT"}` (directive `:113`) forces `activeJunk = 0` even when FORM sets `actualJunkCoefficient = 1` or COEFFICIENT sets `2`. WL K_junk FORM, IDENTITY, and ×2 become identical on every printed object. That is an invisible knife bite caused by the budget envelope, not a change to the frozen knife text. (SymPy covariance `ACTUAL_JUNK` is **not** case-gated and still bites on the same pin.)

**(c) Comparative draws.** Stored WL bound `PER_VALID_DRAW = 22601/499969833 ≈ 4.5×10⁻⁵`. Four draws give a practical FN ≪ 10⁻¹¹ for a newly-nonzero column — adequate **to move a zero→nonzero bitmap**. It is **not** a substitute for digest on same-support FORM / ×2 (Q3). Adaptive `Max[8, 2⁻⁸⁰ …]` is a floor of 8 at this bound; 4 is a real relaxation of the engine’s identity-test budget, not of the comparative DAG digest.

**(d) SOUND for variant comparison.** Four-case `Tuples` order makes `LAB_HELD`/`RHOBR_CONSTANT` ordinal **2**; the single-case driver makes it ordinal **1**. Seeds `731921 + 10000 caseOrdinal + …` (`:763`) therefore differ from production four-case seeds. Every compared label uses the same restricted schedule, so baseline vs corrupted stays aligned. Printing actual seeds (directive `:137–138`) is required for any later cross-check against production.

---

## Q5 — No new leak / no new defect — **finding**

**Error-path / spool clauses are otherwise the right shape.** `/tmp` spool, no embed of symbolic objects/DAGs/sample matrices/raw stdout, bounded failure diagnostics (`:41–43, :62–63, :89–91`), `timeout --kill-after=5 600`, never launch the four-case engine. That would have blocked the ~285 MiB dump **if** every certified object were forced through compact PIT+digest with no exception.

**Residual leak / builder ambiguity:** clause 1 then says small scalar residuals `` `R_N6`, `R_cov`, `SPLIT_CHECK`, slot/closure guard residuals` `` “remain direct payloads” (`:28–30`). In all four engines those names are the **heavy** PIT tables (WL emission also carries `"ARITHMETIC" -> circuitExpression` at `audit.wl:804`). A builder who treats “direct payload” as the engine tag reintroduces the volume failure the revision exists to stop. Charitable reading (compact bitmap as the direct payload) collides with that parenthetical.

**Second builder trap:** “use `emitted_object_sha256` from `RUN_PROVENANCE`” names a nonexistent field (Q2). Combined with “do not reconstruct / do not patch an extractor” (`:34–35`), a builder can hash only the compact fingerprint and ship a transcript that cannot see same-support bites.

The compact-PIT + no-symbolic-fallback language, KB-scale deliverable sentence, and kill-after timeout are sound as intent.

---

## Overall — **NOT-SOUND**

Frozen knife **text** is intact. The output-contract intent (compact PIT + digest, no symbolic dump) and the two WL edit sites are the right kind of fix.

They are not complete:

1. **WL K_junk is inert** under the new single-case envelope (`actualJunkCase` vs `LAB_HELD`/`RHOBR_CONSTANT`).
2. **WL per-object digest is specified against fields the engine does not emit**, so the compact fingerprint (support tally only) can miss a real FORM/COEFFICIENT bite.
3. **“`R_N6` / `R_cov` remain direct payloads”** is a path back to a full-object transcript.

Until (1) the junk case is observable or K_junk is explicitly dropped from this WL harness, (2) the WL digest is defined as a worker hash of the pre-projection object/`PROBE_NUMERATORS`/ARITHMETIC (not a fictional `RUN_PROVENANCE` field), and (3) the direct-payload exception cannot mean the engine’s symbolic tag, the revised contract+budget should not be treated as a cleared build spec.
