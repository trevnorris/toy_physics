# S11c-d wider-box one-/two-momentum refinement

Momentum-domain actions are published at 995d5909 and annex-verified at
ce220b19. The matching profile-order baseline changes the action by 6.94e-18;
cutoff 2→3 and 3→4 changes reach 2.25e-7 and 1.14e-7. The two-momentum rows
dominate those changes. Resolve their wider-box quadrature before interpreting
the finite differences as tails. This is the next planned numerical check,
with no change to the engine, physical input or retained solver/export contract.

1. Consume the accepted cutoff-3/4 actions and all bound native rows, sources,
   profiles, ordered limits and saved field assignments. Verify the publication,
   all current/frozen sources, packets, records, worker layouts and partials.
   Reconstruct each accepted complete action from its saved native terms.
2. Use four workers, one per field/cutoff pair, each with one native thread and
   a 2 GiB address-space ceiling. Each worker retains its accepted full baseline
   explicitly. Hold source/profile orders 256/512, bounds 48/14, regulator 0.2
   and its own momentum box fixed throughout.
3. Refine the 40 single-momentum rows at outer orders 324 then 432. With those
   values retained, refine all 30 paired rows at outer orders 324 then 432,
   keeping inner order 32; then raise only inner order to 48 and 64. The ten
   accepted three-momentum rows remain explicitly held at 144/24/24. No new
   three-momentum integration is claimed. Emit held-layout and native-term
   masks, both operands and literal differences at every stage.
4. Preserve the unchanged native streaming quadrature, current/source/profile
   evaluation, complete-action contraction and tested four-process dispatcher.
   Before production compare actual saved production prefixes at all four
   field/box combinations, isolate every settings change, and check the four
   coarse worker sweeps against serial native groups and independent cell
   contractions. Coarse rules are instrument evidence only.
5. Save every completed group and record and every 64-batch partial before
   later guards. Verify native row/source/profile/limit census, actual measure
   controls, finite masses, exact read-only rule caches, unchanged held values,
   full metadata/emission replay and pre/post packet/source identities. Use
   fresh write keys and the original source/profile cache dimensions. Publish
   only accepted production through DataLad/git-annex and commit checkpoints.
6. Inspect measured refinements before choosing independent outer quadrature,
   wider-box three-momentum refinement or further domain/tail work. This sweep
   cannot establish full wider-box convergence while three-momentum rows are
   held. Uniform/independent-grade and exceptional-locus coverage, infinite
   tails/interchange, Abel limits, scattering and poles remain separate.

Recover instrument or emission failures from saved operands with explicit
source/helper joins. Do not repeat the nine-hour complete-domain computation,
upstream factorization or physics for plumbing. Keep durable repository scratch,
one supervised job and silent local completion/error wake-ups; no polling or
recurring checks. The separately owned Q9 repair does not change these pinned
inputs; any future export rebase requires the recorded consumed-root joins.
