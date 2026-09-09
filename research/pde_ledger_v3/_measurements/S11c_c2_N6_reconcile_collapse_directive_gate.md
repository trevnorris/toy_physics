# S11c-c2 N6 reconcile-collapse instrument — build directive decision gate (G2)

**Directive:** `directives/S11c_c2_N6_reconcile_collapse_build_directive.md` (orchestrator-written; astra will
build the instrument). **Legs (2, identical prompt `_legs/S11c_c2_N6_reconcile_collapse_directive_prompt.md`,
launched on sight, detached):** Codex-sol `gpt-5.6-sol` xhigh + Grok `grok-4.6` high (orchestrator-written → Codex
+ Grok). Reports `_legs/S11c_c2_N6_reconcile_collapse_directive_{codex,grok}.md`. Both EXIT=0.

## Round 1 verdict: BOTH **DIRECTIVE NOT-SOUND**, CONVERGENT (6 computation-changing MUST findings).
The gate paid off — the v1 "frozen dictionary" was descriptions of construction stages, three of which would
have made the collapse VACUOUS or reintroduced η, plus real gaps. Verified each finding (G4) against the cited
engine sites, then folded ALL into v2 (⛔ physics-bearing directive ⇒ review-until-clear, not one-pass).

### Convergent MUST findings (both legs) — verified + folded
1. **Map 6 (velocity/source factors) = whole-object equality of independently-computed subobjects + WRONG
   citations.** Both legs: would suppress a `SOURCE_*` component by assumption; RESOLVED.md keeps face-velocity
   correctness OPEN (`:41-43`). Grok: cited WL `:359,364` are `massSubstrate` constraint/Jacobian + a comment;
   real velocity `.wl:274,304,862-863`, source-solve `:366-381`; PY source factors `diagnostic:388-389`.
   **Folded:** map 6 REMOVED from the dictionary; velocity/source-solve factors stay in the residual (listed in
   the Excluded block).
2. **Map 2 "wave-projection degree 2" = a construction stage, vacuous on the surfaced leaves.** Baked into μ_M
   (PY `brane:1970-1981`; WL `:246`); emitted sources are wave-LINEAR (`diagnostic:368-374`; WL `waveScale[…,1]`
   `:381`). **Folded:** kept ONLY the Jacobian spelling `1+tr(∇u)` as a leftover-factor identification (map 7);
   the degree-2 projector is EXCLUDED.
3. **Map 4 `WBg→W0(1+η·w1)` reintroduces η into already-graded leaves.** Both engines expand profiles before
   grading (WL `:136`→`:141`; PY `diagnostic:247-248`). **Folded:** map 3 (v2) = leftover profile-jet
   NAME/scaling only (`w1ProfileJet{ij}`↔`w1_profile_d{ij}`, `σ_W/L_W`, `μ_R/W_0`); ⛔ no η re-expansion;
   enforced by the grade-combining tripwire.
4. **Φ is ungraded metadata — coefficientwise contract not implemented.** `extract_meta` → `MAP_VARIABLE`/
   `FIELD_PATH`, no grade axes (`comparator:688-765`); "already graded by comparator key" false for Φ; Φ was
   also dropped from the unchanged-dictionary sentence. **Folded:** new "Φ as a graded collapse operand" section
   (extract_meta → profile-expand → explicit 4-coefficient extraction → dictionary + witness); Φ added to the
   unchanged-dictionary set.
5. **Missing primitive: source-wave point-preserving map.** WL bare `sourceMap` jets (`.wl:707-711`) vs SymPy
   applied jets at Y (`diagnostic:392-396`); comparator deliberately does NOT bridge them (`:287-289,990-993`).
   **Folded:** map 2 (v2) = exact `jet[f,I] ↔ ∂_I s11cc2Field_f(Y,t)` point-preserving correspondence.
6. **Controls tautological + census crosswalk missing + operand-pair API + PIT bounds.** One-sided corruption
   needed a CONTROL-DELTA (a leaf already nonzero passes vacuously); maps with no surviving occurrence can't
   "move a leaf"; blanket control tautological if operand:=operand; `compare_family` releases operands
   (`:794-825,862-898`) so no pair API; census needs `kappa_a↔ADVECTION`/`MAX_RANK` crosswalk; PIT needs
   disjoint primes/bounds. **Folded:** controls redesigned (per-map defining-relation control moving a
   control-delta + production-reachability census + blanket ablation with the Control-1 leak dropped + grade
   tripwire); operand-pair extraction respecified (import comparator PRIMITIVES + group-by-key, no comparator
   edit, STOP+report if a helper is needed); census semantic crosswalk added; PIT bounds (disjoint primes/seeds,
   denominator rejection, degree bound, max attempts/time/RSS, undecided-on-exhaustion) added.

### Also folded (leg nits, verified)
Map 3 (live density) → leftover-NAME identification, ⛔ not constant↔field equality (c1 mandated re-adjudication,
not equality — `S11c_c1_comparator_reconcile.md:137-163`); map 8 (c1) → ε/ω conditional-if-atom-present only,
Fourier/on-shell dropped; corrected citations throughout (`diagnostic:382`, `.wl:839-840`, `brane:767-821`,
`:1727-1872`). The Control-1 "shipped run must NOT be all-collapsed" leak REMOVED (value-free).

### Parts both legs found SOUND (kept)
Independent retained rectangle (4 grades); the exclusion list (whole-object / cross-anchoring / σ_W-binding /
default on-shell); BASELINE as nominal-control duplicate; witness three-valuedness + residual-before-witness +
no-verdict; value-free intent; builder fence.

## Architecture note (⭐ the key design change in v2)
The frozen dictionary is now specified as **allowed map TYPES + exact engine SITES + a per-map defining-relation
control**; astra CONSTRUCTS the explicit rewrite tables from the sites and EMITS each as an object. This keeps
the orchestrator from hand-typing error-prone symbol tables (the route-2 failure mode / rule-15 precedent) while
keeping the directive orchestrator-owned; the re-review legs verify each constructed map is a mechanical
name/structure fact that bites. Every admitted map is a NAME/STRUCTURE correspondence — ⛔ no construction-stage
operation is an operand rewrite.

## Disposition
Folded ONCE (all findings convergent + verified). ⛔ physics-bearing directive ⇒ re-review (2 decision legs on
v2) until clear before the build. NEXT = relaunch Codex + Grok on v2.
