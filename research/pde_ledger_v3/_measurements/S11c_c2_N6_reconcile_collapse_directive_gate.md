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

## Round 1 disposition
Folded ONCE (all findings convergent + verified). ⛔ physics-bearing directive ⇒ re-review (2 decision legs on
v2). Relaunched Codex + Grok on v2.

## Round 2 verdict: BOTH **DIRECTIVE NOT-SOUND** again — but CONVERGENT + NARROWER (v1 hazards resolved).
Reports `_legs/S11c_c2_N6_reconcile_collapse_directive_{codex,grok}_r2.md`. Both legs explicitly confirm the
round-1 hazards are FIXED (no whole-object equality remains; map-2 degree-2 + map-4 η-re-expansion correctly
excluded) and that retained order, the three-valued witness, BASELINE-as-duplicate, value-free discipline, the
fence, and the import-primitives extraction reuse are all SOUND. Remaining convergent MUST findings — all about
**explicitly FREEZING the physics-bearing pieces** rather than letting astra construct them:
- **M-a (both) — map 4 (energy-basis coefficient pairing) is not frozen.** The pairing carries physics
  normalization/scale, NOT mechanical spelling: WL `bRho/2·WBg·θ²`, `cCoupling·WBg·θ·eW`, else
  `energyCoefficient<i>` (`.wl:211,217,222-226`) vs SymPy `B_rho_3·W_bg/(2·W_0)`, `C·W_bg`, `kappa_theta/2`
  (`brane:1584-1594`), first-jet `gamma_s11cb_*` (`brane:357,1820`). So the table needs leftover scales
  (`B_rho_3 = bRho·W_0`, `energyCoefficient↔kappa_theta/2`). "astra constructs from sites" delegates a
  residual-changing choice ⇒ the directive must SUPPLY the complete frozen contraction-ID table (orientation +
  scale per pair), or a deterministic reviewed basis-change algorithm + its expected table.
- **M-b (both) — map 6 (leftover density) has no surviving emitted atom name.** Both engines rebind live density
  INSIDE source construction (PY `diagnostic:378,382` → `inputs.density[(ρ,)][1]`; WL `sourceBind` `rhoFace→density3`,
  `density3=density4·WBg` `.wl:839,380`), so `density3`/`inputs.density[…]` are expressions/locals, not emitted
  names. Map 6 is therefore either inert (its control can't bite) or a new constant↔live-expression rewrite
  (reopens the c1 rule-17 hazard). ⇒ REMOVE from the production rewrite dictionary; keep live-density only as a
  separately-emitted premise/census comparison.
- **M-c (Codex) — typed map scopes.** The source-wave bare↔applied-`Y` map is justified ONLY at the source
  extraction boundary (`diagnostic:392`; WL `sourceMap` bare jets `.wl:707`); applying it to Φ (an abstract jet
  map with no `Y` evaluation, `cov:63`; WL `:236`) changes the claim from abstract-jet-map equality to
  after-point-evaluation equality. ⇒ freeze each map with a typed family/stage scope (source-wave → `SOURCE_*`
  only; Φ gets domain/multi-index spelling + pre-grade profile substitution, ⛔ not source-point evaluation).
- **M-d (both) — controls incomplete.** Corruption control needs LOCALITY (nonzero delta on ≥1 declared-use
  leaf AND zero delta OUTSIDE the use-set — else an over-broad map passes); add a real DROPPED-map control
  (delete each active primitive, require a residual on a predeclared defining-relation operand); run blanket
  through the same extraction/rewrite/witness path.
- **M-e (Codex) — census crosswalk must be FROZEN exactly, not "e.g."** Enumerate `DOMAIN`/`COVERAGE`/
  `UNCOVERED`/`MAX_RANK` ↔ SymPy counterparts + `kappa_a↔ADVECTION`, `kappa_j↔JUNK`; ⛔ `MATERIAL_NORMAL` is
  ONE-SIDED (WL-only, like `THICKNESS`) — do not synthesize a SymPy zero for it.

## Round 2 disposition — ⭐ RULE 15: CHANGE THE AUTHOR (2nd heavy round at the build-directive gate).
Two decision-gate rounds NOT-SOUND on the same material (the frozen dictionary). The remaining fixes are exactly
the error-prone explicit-table specification (map-4 scales, census crosswalk) that hand-authoring has repeatedly
gotten wrong (the N6 route-2 precedent: hand-written maps failed 3× → Codex authored → CLEAR). The rule-15-armed
condition ("a 2nd heavy round at the build-directive gate ⇒ hand re-author to Codex") is met. ⇒ **DELEGATE v3
authoring to Codex-sol**: fold the round-2 MUST findings into v3 with the explicit frozen map-4 table + census
crosswalk + typed scopes + map-6 removal + corruption-locality/dropped-map controls; keep everything both legs
found sound. Authorship changes O→Codex ⇒ v3's review legs are **fresh Claude + Grok** (Codex-written → not
Codex-reviewed). NEXT = Codex authors v3 → fresh-Claude + Grok review-until-clear → astra build.
