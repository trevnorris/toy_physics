# v3 directive review — fresh Opus Agent leg (Codex-authored directive → fresh Claude + Grok)

**Artifact:** `directives/S11c_c2_N6_reconcile_collapse_build_directive.md` (v3). Independent derivation-based review.

## VERDICT: DIRECTIVE SOUND (clear to build)
⚠ Adjudicated by the orchestrator as a FALSE-NEGATIVE on the map-4 graded-application stage — see the gate
record; the Grok leg caught a real MUST this leg missed (L-R13). This leg's positive derivations (the 19-row
table is correct UNGRADED) stand and were used to confirm Grok's finding is a STAGE error, not a table error.

## Map 4 — the 19-row energy-basis table: re-derived every row from the two constructors, CONFIRMED CORRECT (ungraded)
- R_W / uniform rows from `uniform_coefficient` (`brane:1586-1595`): E01 `θ²` catalogue `B_rho_3·W_bg/(2·W0)` ⇒
  `bRho → B_rho_3/W_0`, both ½ retained; E02 `θ·eW` with `E=W0·e_W/W_bg` (`:1631,1658`) ⇒ `cCoupling → R_W·C`;
  E04 `r·q` with `br=D(E)=R_W(∇e_W − e_W∇W_bg/W_bg)` (`:1626-1642`) ⇒ `energyCoefficient4 → R_W·kappa_theta_W`
  + a leftover feeding E14; E07 `q·q → kappa_theta/2`; E09 `θ·trG → G_theta_u`.
- E14 triangular mixing exactly right: the E04 leftover (`eW·(gW·q)`, one background jet, retained) adds to the
  direct `be·(gW·q)` term ⇒ `R_W·(gamma_s11cb_w_bg_14 − kappa_theta_W/W_bg)`; the E03/E13 two-background-jet
  chain pieces correctly dropped at retained order.
- First-jet rows (E03,E05,E06,E08,E10-E13,E15-E19): independently enumerated `enumerate_new_candidates`
  (`:1510-1527`) + `delta_contractions` (`:1334-1353`) + `perfect_matchings` (`:1313-1323`); reproduced the
  index→structure map exactly (04=θ·(g·u), 06/07=g_a·q_i·G_ai swaps, 08=(g·q)·trG, 12=θ·(g·q), 13=θ·(g·r),
  14=be·(g·q)); gamma/contraction pairing index-matched in `construct_energy:1819-1822`; MU_R_BG vs W_BG source
  assignment correct every row.

## Other maps (1,2,3,5,6,7,8): all grounded, primitive, no whole-object equality
Map 1 grad-θ dual (atom map); Map 2 source-wave at Y scoped S_A/S_P/S_B, Φ excluded; Map 3 profile scales
(σ_W·profile/L_W^(r-1)), numeric-leaf spelling only; Map 5 Φ spelling only; Map 6 correctly REMOVED (both engines
substitute density inside source construction — census-only `LIVE_DENSITY_PREMISE`); Map 7 Jacobian occurrence-
gated; Map 8 ε/ω occurrence-conditioned. Excluded list correct.

## Retained order / witness / controls / census / value-free / fence: all SOUND
Coefficientwise 4-grade; witness three-valued, residual-first, disjoint PIT primes, bounded, no verdict;
controls bite (corruption-locality bidirectional, dropped-primitive on defining-relation operand, blanket
unreachable from production, grade tripwire); census crosswalk exact (6 pairs), THICKNESS/MATERIAL_NORMAL/JUNK_*
WL-one-sided; BASELINE nominal-control duplicate; value-free; fence + DoD + comparator-primitive reuse all sound
(all named primitives exist).

## Non-blocking observations (FOLD both)
- **OBS-1:** could not hand-reproduce the WL `energyCoefficient<i>`↔contraction assignment for E03-E19 (derives
  from Mathematica `DeleteDuplicates` + `RowReduce` pivot order); verified the SymPy side exactly + each row
  names the same physical contraction. A WL-index error fails SAFE (mismatched rename → nonzero residual →
  `residual_remains`, never false `collapsed`). ⇒ recommend the builder `test_` ASSERT the WL
  `energyCoefficient↔contraction` map against the actual `constructEnergy` output, not only the SymPy side +
  rank/invertibility.
- **OBS-2:** rename-only rows (no R_W: E06-E12, E15-E19) collapse essentially by construction; the load-bearing
  discriminators are the R_W normalizations (E01-E05, E13, E14) + the contraction-structure match. ⇒ the
  orchestrator should weight the R_W rows + the surviving residual when adjudicating, NOT the count of collapsed
  rename-rows (a note for the disposition, not the directive).
