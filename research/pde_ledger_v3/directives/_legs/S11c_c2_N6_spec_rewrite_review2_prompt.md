# Spec rewrite RE-REVIEW (round 2) — corrected c2 §5c (N6). STATIC only (⛔ run nothing, load no `.out`).

You are ONE of two independent legs. **Round 2**: round 1 (Codex-sol + Grok) confirmed the N6 *axis* is now
correct (Eulerian vs material-coordinate of the SAME increment at fixed `α,ρ`; cross-anchoring reclassified; `Δρ`
forbidden as an anchoring bridge; M2 clean), but **blocked commit** on the items below. The orchestrator folded the
corrections. Verify each is now correctly and completely fixed and introduced no new defect. ⛔ Do not rubber-stamp;
a real problem is a finding.

## Artifact
Current `directives/S11c_c2_SHARED_PHYSICS.md` §5c + diff vs HEAD:
`git -C /var/projects/toy_physics --no-pager diff -- research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md`.
Consistency sources: c2 §3c (increment def), §3a (the fold = §3a δp_s-substitution), §7 (imports); parents
`S11c_a_SHARED_PHYSICS.md` §2c/§5a; `S11c_b_SHARED_PHYSICS.md` §5a (the route-2 material-slab template); sibling
`S11c_c1_SHARED_PHYSICS.md` §5a (the bulk-DtN global-scaling prohibition); `S11c_decisions.md` N4/N6.

## The round-1 blockers that were folded — verify each
1. **Route-2 object naming (the commit blocker).** §5c must now name the objects, not a recipe: `SLAB_M` = **S11c-b
   §5a's material-coordinate construction of `slab_operator` at the same `α`** (constructed in this control, not
   imported); `close(SLAB_M)` = **§3a substitution of the imported same-`α` c1 `s11c_c1_face_response`** into
   `SLAB_M` (⛔ not a second DtN, ⛔ not c1's Hanzawa operand); `T_{M→E}` = `Δρ` map. And it must ⛔ forbid
   reconstructing the bulk DtN by S11c-a's global scaling (c1 §5a: secular at infinity). **Check the current
   wording makes route 2 compute the SAME object as route 1** (a build cannot substitute a bulk-Hanzawa or a
   second DtN). Cross-check `SLAB_M` really is what S11c-b §5a defines, and the prohibition really matches c1 §5a.
2. **Route-1 gloss.** Round 1 flagged "direct Eulerian graph/level-set" as S11c-a substrate language, not c2's fold.
   Is route 1 now described as c2's §3c increment built in Eulerian variables (§3a substitution into the imported
   Eulerian `slab_operator` + §3c extraction)?
3. **One-sided corruption.** Is the advective (`u·∇Σ_E⁰`) probe now on the **material→Eulerian map-back** (⛔ not the
   `δp_s`-free slab base, which cancels per §3c), the tilt probe on the Eulerian route, the "each route move alone"
   wording fixed, and a printed emit package `S11CC2_CONTROL_INDEPENDENCE_{BASE,CORRUPTED,RESIDUAL}[α,ρ,probe]`
   added (print-not-assert, M2)? Any remaining "require it to move" that reads as a builder assertion?
4. **Cross-anchoring mandatory + precise name.** Is `S11CC2_ANCHORING_L_MINUS_M[ρ]` now **mandatory** (not "if
   emitted") and defined with the precise `S11CC2_SELF_ENERGY_INCREMENT[α,ρ]` objects (not generic `SelfEnergy`),
   still with ⛔ no zero target and ⛔ never labeled a rep-invariance residual?
5. **No new defect / consistency.** Any contradiction introduced with §3a/§3c/§4/§5d/§5e or the parents/sibling? Is
   §5e's uniform-limit smoke test still coherent beside the corrected §5c? Does the fold accidentally leak a target
   or manufacture a derivation path?

## Output
For each of 1–5: finding + quoted lines (CONFIRMED vs unsettled). End: is the corrected §5c **now correct to
commit**, or the exact remaining wording that must change.
