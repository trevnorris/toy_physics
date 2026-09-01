# S11c-b #89 WL §3a-repair directive — decision-review outcome (2 legs) + fold + scope split

⚠ ANSWER-BEARING (records the withheld completed count). For the orchestrator + REVIEW legs only; ⛔ never
handed to the WL builder. The build directive `S11c_b_89_wl_3a_repair_directive.md` stays leak-clean.

## Legs (orchestrator-written directive ⇒ Codex + Grok, rule 7 TRIGGER, before any builder)
- **Codex** (`~/.s11_build/S11c_b_89_wl_decision_codex.txt`; evidence
  `~/.s11_build/S11c_b_89_wl_decision/bilinear_quotient.{py,out}`): VERDICT 7 findings (4 blockers).
- **Grok** (`~/.s11_build/S11c_b_89_wl_decision_grok.txt`; evidence
  `~/.s11_build/S11c_b_89_wl_dirleg_grok/raw_vs_quotient.{py,stdout}`): VERDICT 7 findings (2 blockers).
Both did independent abstract-Kronecker computation (imports neither engine). I verified each finding (rule
13): findings 2/3/5/7 against the WL code directly, finding 3 against Codex's booleans, the operator freeze
against Grok's frozen-vs-live EL ranks and the code path.

## The converged measurement (both legs' scripts agree)
- Complete **O(3)-Kronecker** first-jet field-bilinear family = **15/source**; raw rank = live-quotient rank
  = 15 (nullity 0). Frozen-coefficient quotient = 8/source. Two-source w/ uniform 10 ⇒ RAW/LIVE **40**,
  FROZEN **26**. (Grok stdout `RAW_TWO_SOURCE 40`, `FROZEN_QUOTIENT_TWO_SOURCE 26`; Codex `raw_equals_live_
  quotient=True`, `frozen_quotient_has_positive_nullity=True`.) ⇒ the directive's basis-count PREMISE HOLDS:
  WL's plain `MatrixRank` over first-jet `rankVariables` reaches the correct count once the family is
  complete; ⛔ no divergence-quotient recipe and no leaked integer needed. §3a "independence as field
  bilinears (B1 not applied)" forbids the *constraint* quotient, not the divergence quotient; the two
  readings AGREE on the count for the complete family (they diverge only under the *frozen* coefficient).
- **Levi-Civita adds 4 parity-ODD pseudoscalars** (`RAW_GAIN_FROM_LEVI 4`); the S11b group has parity +
  unbroken `w→−w`, so these are FORBIDDEN. The enumeration must be O(3)-Kronecker only.
- **Operator freeze is REAL** (the big miss): the slab operator `evaluatedModel` (L905–911) applies
  `applyProfile` (frozen 2nd/3rd jets, L297–310 `higherRules`) to the energy records BEFORE the EL variation
  `constrainedRows`→`variationalSource` — the opposite of §3b (276–278 "do not freeze a coefficient before
  differentiation"). Frozen-EL rank **8** vs live-EL rank **15**; 6 of the current 8 invariants carry live
  Hessian in EL (`WL8_LIVE_HESSIAN_ATOM_COUNTS [0,0,6,6,6,6,6,3]`). The Hessian is the non-absorbable
  operator content (matches #88). The live path exists (`applyBackgroundProfileWithGeneratedJets`, L360) but
  is used ONLY for admissibility (L1337–1341), not the operator. ⇒ Completing the four lists repairs the
  basis COUNT but leaves SLAB_OPERATOR/COUPLING_KERNEL at frozen-coefficient EL.
- My rule-17 control (frozen-vs-live `applyProfile` on the energy candidates) is TAUTOLOGICAL: candidates are
  first-jet only, so the energy-poly rank difference is 0 regardless (`ENERGY_POLY_RANK_DIFFERENCE 0`). The
  lethal freeze is in the EL, not the monomials.

## Findings adjudicated (all VALID; verification per rule 13)
1. **Operator freeze denied + mis-cited (Codex#5 SHOULD-FIX, Grok#1 BLOCKER — I take BLOCKER).** Directive
   §2's "no separate operator freeze (L568)" is false; L568/L576 are the admissibility constructor, not the
   slab operator (`evaluatedModel`). VERIFIED in code + Grok's EL ranks.
2. **Rule-17 control on the wrong object (Codex#3, Grok#2 — BLOCKER).** Replace the energy-poly frozen-vs-
   live control (keep it as a Hessian-in-energy guard, emits 0) with a load-bearing **operator-row
   frozen-EL vs live-EL rank** emit pointed at `evaluatedModel`.
3. **Enumeration not uniquely fixed (Codex#4 BLOCKER, Grok#3 — the real completeness control).** "admissible
   / symmetry-reduce / energy dimension / divU and the full shear" are recipe words a subset can satisfy.
   Name it as the object: every **O(3) Kronecker** scalar of one first-jet spurion × a quadratic in
   {u,∇u,θ,∇θ,e_W,∇e_W}, ⛔ NO Levi-Civita (parity), `w→−w` imposed, exact thickness map imposed, then raw
   rank. ∇u has THREE independent Kronecker contractions with (s,u), not one. The free coefficient supplies
   energy dimension (the bare invariant does not).
4. **Thickness map not imposed on the new invariants (Codex#2 — BLOCKER).** Uniform terms use mapped
   `localEw` (L463–467); the new invariants pass RAW `ewWave`/`ewVariation` (L505–506, L572–574), violating
   §3a (246–247, map before rank). Codex `exact_map_changes_retained_eta_sigma_grade=True`. VERIFIED in code.
5. **Dimension check self-validating (Codex#6).** `newInvariantDimensions` is a position-keyed hand list
   (L1750–1764); coeff dim = energy − inv dim (L1766–1768); homogeneity adds them back (L1828–1831) → an
   incorrect invariant dimension passes. Derive invariant dims independently, compare, then coeff dims.
6. **Anchoring agreement tautological (Codex#7).** `basisRepresentativeIndices` computed once from LAB_HELD
   (L621–623), reused for both branches (L1208–1211); counts are the length of the same slice → agree by
   construction. Compute rank/span independently per anchoring, emit the span-equivalence residual.
7. **`timeout 600` on the full regeneration (Grok#5).** 600 s is the ABLATION budget; the full ~97 MB run is
   silence/RSS-bounded, not elapsed. Keep 600 s for ablations only; full run uses the historical WL budget.
8. **Blindness prohibition is a denylist (Codex#1 BLOCKER vs Grok#6 NIT — RESOLVED to Grok).** rule 12: "do
   not read PY" is a request, not a control; the engine imports nothing anyway. The REAL control is the
   UNIQUE enumeration (finding 3) — a deterministic O(3)-Kronecker family removes any "tune to 15" freedom,
   so the leaked answer in #86/#87/this note is harmless. Remove the prohibition; do not cite answer-bearing
   measurements in the directive; keep the count withheld. (Filtered-worktree build is optional defense-in-
   depth; the program's precedent is in-repo blind builds + from-scratch review legs + form ablation.)

## SCOPE DECISION — SPLIT the operator unfreeze out (both legs bless deferral if honest)
#89 WL is scoped to the **BASIS** (enumeration completion + thickness map + the basis-level tautology fixes
5/6). The **operator unfreeze** (EL-before-freeze then live-path reduction) is a DIFFERENT mechanism on a
DIFFERENT surface (`evaluatedModel` + all operator consumers) — bundling two large unrelated changes is a
rule-15 defect risk. It becomes **#89b WL** (its own 2 decision legs + build), the WL analog of PY's #88
consumer fix, and a prerequisite of the integration re-adjudication. The #89 WL directive therefore: (a)
RETRACTS the false "no operator freeze" claim; (b) states the operator freeze is a KNOWN separate defect NOT
closed by #89, deferred to #89b; (c) emits the operator-row **frozen-EL vs live-EL rank** as a computed
DIAGNOSTIC (documents the freeze in the `.out`, asserts nothing). PY parity note: PY #89 fixed its operator;
WL reaches operator parity only after #89b, so the integration operator comparison waits on #89b.

## Disposition
Directive REWRITTEN folding 1–8 once (rule 7: one decision pass, fold, go — ⛔ not re-reviewed to green; the
build's 2 legs check the implementation). NEXT = leak-gate the rewrite → Codex build (WL, danger-full-access,
xhigh) → 2 build legs (fresh Claude + Grok; SERIALIZE Mathematica). Then #89b WL operator, then integration.
[[feedback_basis_independence_must_not_freeze_spurion]] [[feedback_decision_list_length_is_the_defect_rate]]
