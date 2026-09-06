# S11c-c2 export repair — re-review adjudication (gate record)

**Artifact:** the astra publication-only export repair — the diff of
`scripts/S11c_c2_selfenergy_fold_sympy_audit.py` from the reviewed baseline `8f3a017f` (publish rewrite +
`publication_compact` helper + a strict-recursive semantic guard; `EXPORT_ROOTS`/`export_key` membership) and the
regenerated `scripts/S11c_c2_exports.py` (**60,516,900 → 22,441,522 bytes**, 63% smaller).
**Role:** script / publication-only build → review-until-clear. **Authorship → legs (G1):** astra(Codex)-written ⇒
fresh Claude agent + Grok. Identical prompt: `directives/_legs/S11c_c2_export_repair_rereview_prompt.md`.
Blindness by absence: astra's builder report relocated during the review.

## Leg outcomes
- **Grok** (`_measurements/S11c_c2_export_repair_rereview_grok.txt`): all 6 checks — **CLEAR to commit**. Physics
  unchanged (baseline vs current TRIAGE byte-identical; construction functions unchanged, only the export map
  differs in `run()`), membership minimal (increment absent, both operators, closure minimal), semantic
  equivalence clean, transparency clean (no Dummy/UnevaluatedExpr), **guards BITE** (its `out+1` compaction
  corruption hard-stopped the semantic guard; membership re-add/drop both caught), no overstep.
- **Fresh Claude agent** (incomplete — wedged twice on reaped background waiters; stopped by orchestrator): checks
  2 (membership clean), 4 (transparency clean), 6 (no overstep) corroborate Grok. BUT it reported two APPARENT
  defects — both **disproven** (see below): a semantic recompute FAIL and a "B2 corruption returned clean". Its
  cross-process re-run method introduced fresh Integral bound-dummies, producing artifacts
  (`_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out`: the nonzero leaves are
  `2·_Dummy_34·X − 2·_Dummy_44·X` — identical term, relabeled dummy = the artifact signature).

## The leg-vs-leg disagreement — SETTLED by my own same-process computation (rule 13)
Command: `_measurements/S11c_c2_export_repair_adjudication_verify.py` → `…_verify.out`. Built one case in-process
(no serialization/cross-process dummy mismatch), ran astra's `publication_compact` on every closed-slab-operator
leaf, and applied astra's own algebraic-leaf comparison (`normalize_integrals` + Integral-protected expand):
```
closed slab operator: 5 algebraic leaves
FAITHFUL (compact==emitted):      5/5  fails=[]
GUARD BITES x2 corruption:        5/5  slipped=[]
GUARD BITES +marker corruption:   5/5  slipped=[]
INTEGRAND corruption on /THETA: bites=True
```
⇒ `publication_compact` is **faithful** (loses no information — the Agent's semantic FAIL is a cross-process
dummy artifact, not a real difference), and the semantic guard **bites** every real value change — doubling, a
distinct marker term, and an in-integrand change (the Agent's "B2 returned clean" was its own no-op corruption,
not a guard gap). This independently corroborates Grok and refutes both of the Agent's apparent defects.

## Verdict — the export repair is TECHNICALLY SOUND; ⚠ ONE usable leg + a PROCEDURAL gap (corrected 2026-09-06)
**Technical soundness (well-supported):** faithful compaction (my 5/5 + Grok + astra's in-run guard), guards bite
(my check + Grok), membership minimal + correct (increment EMIT-only; both closed operators, all 4 cases;
exported==required closure), transparency (no Dummy/UnevaluatedExpr; `display` shrunk to the root name), singular
locus preserved (astra's pole-revert + `reciprocal_powers_unchanged`), physics unchanged (no construction touched;
TRIAGE byte-identical), no overstep.
⚠ **PROCEDURAL GAP (Codex-sol compact-prep verify): the 2-leg discipline was NOT met.** Only **Grok** produced a
complete usable CLEAR report; the fresh-Claude leg wedged and never produced a usable report, and my rule-13
verification — however thorough — is **NOT a leg** (`CLAUDE.md`: both reports usable; orchestrator verification is
never a leg; no commit before both legs report). ⇒ calling this "2 re-review legs / re-reviewed clear" was an
**over-claim**. The commit `aa76105a` stands as the technically-supported reviewed baseline (Grok + my verify), but
**one complete usable fresh-Claude publication re-review leg is OWED** before procedural clearance — run it (with
foreground-only instructions to avoid the wedging) BEFORE advancing to the WL build.

⚠ **Process note (for leg-quality tracking):** the fresh-Claude leg repeatedly wedged on reaped background waiters
and its cross-process semantic method produced false-positive FAILs; a same-process comparison (as Grok and I did)
is the correct method for Integral-bearing operators. Its methodology, not a real defect, drove the divergence.

**Next:** commit the repaired baseline; then the light §5e/§3c spec-wording clarification; then the blind Wolfram
engine → its 2 legs → the c2 T7 comparator + reconcile → the c2 step record. ⛔ c1 STANDS.
