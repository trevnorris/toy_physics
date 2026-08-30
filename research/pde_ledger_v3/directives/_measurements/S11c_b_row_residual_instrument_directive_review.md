# S11c-b row-residual instrument — decision-list review (2 legs, orchestrator-written → Codex + Grok)

Directive reviewed: `directives/S11c_b_row_residual_instrument_build_directive.md` (v1).
Legs (raw): `~/.s11_build/S11c_b_row_residual_dirreview_codex.txt` (Codex, xhigh, 9 findings),
`~/.s11_build/S11c_b_row_residual_dirreview_grok.txt` (Grok, high, 8 findings). Both blind to my step-2 note.

## CONFIRMED (both legs, independent derivation): the requested-truncation reading
Both legs independently derived from §2a (`S11c_b_SHARED_PHYSICS.md:202-207`) + §3a (`:248-254`):
**retain a monomial `ε^c η^a σ_W^b` iff `c≤1 ∧ a≤1 ∧ b≤1`; grade by bookkeeper power, not spatial-derivative
order; Taylor-LINEARIZE coefficients** (`W_bg²→W_0²(1+2η w1)+O(η²)`, `1/(1+η w1)→1−η w1+O(η²)`).
§3b's "retain full spatial dependence" = do not FREEZE before differentiation, and does NOT override the
amplitude truncation. **Both: NOT a spec ambiguity.** ⇒ the pivot is settled (3 independent reads incl. mine):
kinetic (2,0)=η² and coupling a≥2/b≥2 are OUT of the requested truncation; coupling agrees within scope at
(0,1),(1,1) IF the in-scope expressions match.

## FINDINGS (folded once; rule-13 verdicts). C#=Codex, G#=Grok; merged where they overlap.
- **[LEAK] C1=G7 — VALID.** §0 L7 "twenty aligned cases carry an independently-nonzero per-bucket
  difference" + L9-11 "engines partition the same strong-form operator" pre-disclose nonzero residuals and
  pre-classify as representational; L130-135 "residual must move (advective/flux)" leaks an expected effect.
  Contradicts §382-383 "No residual value is supplied" / §479-482. FOLD: value-free §0 (only "per-bucket ≠
  operator agreement"); ablate on synthetic fixtures or emit effects without expected value.
- **[WRONG OBJECT — assembly] C2=G1 — VALID.** The step-1 twenty-case set is a SELECTED-DIFFERENCE registry
  (`S11c_b_background_multigrade.py:54-101`), NOT a partition registry; it omits BULK_ENERGY/FACE/FLUX/
  ACCUMULATION. FOLD: assemble each engine BY ITS OWN origin schema (PY L1573-1583: BULK_ENERGY/KINETIC/
  FACE_FLUX/ADVECTIVE; WL L851-862: KINETIC/BULK_ENERGY/FACE/FLUX/ADVECTIVE/ACCUMULATION); align the
  resulting SEMANTIC rows; do not use TARGET_IDENTITIES as the assembly universe.
- **[SIGN + non-additive] C4 — VALID, critical.** PY EXPANDED SUBTRACTS kinetic (`- rhobr·u_tt`,
  `- mu_W·W_bg²·e_tt`, L1470-1490) while KINETIC origin stores them POSITIVE (L1576-1580) ⇒ a literal
  bucket-sum flips the kinetic sign. LOCAL/DIVERGENCE_FLUX/EXPANDED are 3 REPRESENTATIONS of one
  contribution, not 3 additive origins. WL ACCUMULATION (L860) omitted from my sum. FOLD: use each engine's
  OWN complete emitted row with its OWN signs — PY EXPANDED (mechanical) / ADVECTIVE (mass); WL
  operator[row] — never a hand-specified bucket sum.
- **[μ_θ vs θ vs mass] C3=G5 — VALID.** PY `theta_expanded = ε·mu_theta_amplitude` (L1478-1484) is the
  CONSTITUTIVE μ_θ operand, mapped by the extractor to `MU_THETA` (`S11c_b_cross_engine_comparator.py:766-778`),
  NOT a θ EOM; `ADVECTIVE_MASS_OPERAND → MASS_EVOLUTION_ROW`. WL has NO θ row (rows: U_MOMENTUM/THICKNESS/
  MASS_EVOLUTION/CENTER_FACE, L839-847). FOLD: drop the "θ-balance" strong row; compare MU_THETA separately;
  compare the mass-evolution row (this is the ADVECTIVE family).
- **[PY faces not additive] C5=G3 — VALID.** PY FACE_FLUX = the whole S11c-a substrate bundle (traction/
  virtual-work/…, L1434-1455), attached as FACE_FLUX_BOUNDARY_OPERANDS (L1543-1550), NOT added to EXPANDED;
  WL folds faces into generalized U_FACE/EW_FACE rows (L807-814) then into operator[row]. The layer does NOT
  do this variational conversion. FOLD: faces are NOT among the 20 differences — scope the strong-row
  comparison to the commensurate bulk+kinetic(+mass) content; compare faces separately only if a reviewed
  PY-generalized-face adapter is built (not needed for the 20).
- **[COUPLING weak + both blocks] C6=G4 — VALID.** §3c is a WEAK restriction, IBP boundary term fixed to 0
  by compact support (§292-308) ⇒ residual defined modulo EXACT total in-plane divergence (full `∇·(cF)`,
  product rule; §1d: `c∇·F` alone is NOT discardable, the `−(∇c)·F` is physics). Layer already gates this:
  `_is_weak_scalar_density`→`classify_total_divergence` for COUPLING_KERNEL only (L779-782, 965-1018).
  Registry is one-sided (only FORWARD/TRANSVERSE_TO_THICKNESS) — §3c requires BOTH transverse→thickness AND
  thickness→transverse blocks; do NOT manufacture reverse as "± forward". FOLD: coupling → weak modulo exact
  total in-plane divergence via the layer's classifier, at requested truncation, BOTH blocks + relabelled
  adjoints.
- **[strong-row zero test] G8 — VALID.** Do NOT apply a total-in-plane-divergence quotient to the STRONG
  slab rows — §1d: variable-coeff IBP generates first-jet terms that are PHYSICS. Strong rows = EXACT
  difference (integral linearity only for ∫-splitting canon, NO divergence quotient); coupling = weak route.
- **[admissibility componentwise] C7 — VALID.** §3d admissibility = ε⁰ background first variation, ordered
  pairing (bulk-DOF body force + per-face traction); components live in different slots, MUST NOT be summed.
  Both engines emit U/THETA/E_W body + both tractions (PY L2230-2256, WL L1336-1349); the reused registry
  had only BODY_FORCE/DOF=THETA. FOLD: compare the full ordered admissibility association COMPONENTWISE
  (3 momentum, θ, e_W, each traction), all 4 anchoring/density cases, attach c=0, NOT via slab-row assembly.
- **[ε + no full residual] C8 — VALID.** WL applies `truncateBackground` BEFORE emission (L154-160, 994-1003)
  ⇒ WL operands ALREADY truncated to (≤1,≤1); an "untruncated ROW_RESIDUAL_FULL" compares a fuller PY vs a
  pre-truncated WL — asymmetric/meaningless. The extractor strips PY ε via `coeff_epsilon`
  (`…comparator.py:455-466`); WL carries ε as metadata. FOLD: DROP ROW_RESIDUAL_FULL; attach ε order c from
  family METADATA (c=1 wave objects, c=0 admissibility), never by symbol-counting; truncate PY to (≤1,≤1)
  (WL already is).
- **[tautological guards] C9=G6 — VALID.** WINDOW_CLEAN with `remainder := operand − polynomial` is 0 by
  construction (same-route Taylor bookkeeping, `S11c_b_background_multigrade.py:183-200`) — an arithmetic
  bookkeeping guard, NOT an independent truncation validation. WL BUCKET_PARTITION_CHECK is vacuous (rows and
  origins from the same local expressions in one constructor). FOLD: no hard-stop on the reconstruction
  identity as a physics check; validate the truncation against an INDEPENDENTLY-implemented CAS-series
  projection + one-sided synthetic corruption; keep a partition check only where an independent complete-row
  route exists (PY per-term vs full-density EXPANDED), described as assembly-accounting, not physics.
- **[DIVERGENCE_FLUX mis-bucket] C7'/G2 — VALID minor.** DIVERGENCE_FLUX is a FORM (LOCAL/DIVERGENCE_FLUX/
  EXPANDED) not an origin. FOLD: strike from the origin list.

## FOLD DECISION (rule 15 / [[feedback-delegate-build-when-prose-directive-repeats]])
The per-family assembly is BESPOKE and not pre-enumerable in prose (both legs proved a uniform recipe wrong).
⇒ revise to a SURGICAL directive: LOCK the verified physics (the (≤1,≤1) truncation reading; per-family
object definitions; strong-exact vs weak-modulo-divergence split; ε-from-metadata; componentwise
admissibility; both coupling blocks); REUSE the layer's existing per-family adapters (`_kinetic_pairs`,
`_admissibility_py_parts`, the weak coupling route, the extractor row mapping) with VERIFIED refs + line
numbers; DELEGATE per-engine assembly (own schema, own signs) to the builder; DEFINITION-OF-DONE with
accounting (every emitted origin accounted; both coupling blocks; admissibility componentwise; c from
metadata); value-free; genuine guards. Two legs done + folded ONCE → Codex build → 2 build legs gate the
artifact (rule 14 ablation). ⛔ do NOT re-leg the directive (iterating toward green).
