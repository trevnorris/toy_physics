# S14a part 1 decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_decision_list_review_p1_lookups.md`, beside this file
(generator `_scratch/s14a/gen/dl_p1_lookups.sh`). "K-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S14a_1_decision_list.md`, orchestrator-written, sha256 `ec48eb74…`, with
  `S14a_1_setup_lookups.md`. Both are preserved unchanged in the commit that adds this file. Both matched the
  baseline after the pass, and the working tree was unchanged (K1).
- **Role.** A pre-builder decision list (G2: one Codex + Grok pass, fold once). It also names the object the build
  computes, so that content is reviewed until clear (spec row).
- **Prompt and legs.** `_scratch/s14a/dl_review_prompt_p1.md`, sha256 `ad95e580…` (K1). Both legs received it
  identically, in parallel. Codex (xhigh, workspace-write in its scratch) wrote
  `/var/projects/_review_scratch/s14a1-dl-p1-codex/REVIEW.md`, and Grok (`grok-4.7`, effort high) wrote
  `/var/projects/_review_scratch/s14a1-dl-p1-grok/REVIEW.md`. Both exited 0 (K1).
- **Verdicts.** Codex: FINDINGS, 3. Grok: FINDINGS, 4 (K2). Each finding below carries its leg's script and
  filed stdout.
- **What Grok checked and found sound** (its "Checked, and not a separate finding"):
  - the supplied pair matches stage 006;
  - `J` is the sum of the two channels in the 2026-10-21 decision, with the order object printed separately;
  - no value, steady state, `ω → 0`, `J_χ = 0` or exchange sign is given to the build;
  - v2's labels are record-only;
  - "counted in full" fixes the full weight of a face's concentrated part.

## Findings

**P1. The face terms and the face channel are not fixed as material crossing moving, tilted faces.** (Grok 1.)
**CONFIRMED.**
- **The list.** Item 1 asks for "every term the integration leaves at each face", and item 2 for the face channel,
  "the net rate at which material leaves the shear-holding region through its faces". The model point keeps the
  faces moving (K3).
- **The source.** S12's material face factor keeps the face's motion and tilt:
  `a_{r,s}=u_r^w−∂t w_s−u_r^i∂i w_s` (S12 record l.79, K3).
- **The evidence.** Grok's three witnesses (moving, tilted, and both) compare two readings. Reading A is the
  endpoint coordinate flux. Reading B is material crossing each face relative to its motion. Only B closes storage,
  `J` and the outer flux: the closure prints `0` for B and `1`, `1`, `2` for A (K8).
- **Why it changes what is computed.** Item 1 does not say the projected balance is for the time rate of the
  projected density. So a builder may leave the faces' motion inside the integral, and the face channel and `J`
  then differ.
- **What must be true.** Item 1's projected balances are balances for the time rate of each projected density and
  the in-plane divergence of each projected flux, with the faces' motion and tilt kept. The face terms are whatever
  that leaves at each face. The face channel is the net material crossing the faces relative to their motion. S12's
  face factor may be cited as the supplied, verified object it is.

**P2. "Every profile a live function of radius and time" freezes the bulk coordinate.** (Codex 2, Grok 2.)
**CONFIRMED.** Both legs found it.
- **The list.** Item 1 says "Every field stays live, as S12 kept them". Item 2 says "with every profile a live
  function of radius and time" (K4).
- **The source.** S12's fields live on `(t,x¹,x²,x³,w)`, with order `f(t,x,w)` (S12 record l.62–63). The decision
  list puts `n` "on the 4D domain (in-plane `x^i`, bulk direction `w`)" (K4).
- **The evidence.**
  - Codex's spherical fixture keeps `w` dependence and satisfies both balances. Flattening its native fields in `w`
    moves the volume flux (`flattened_minus_preserved = 2/5`, K8).
  - Grok's `w`-dependent column has face leaving `exp(-t)`. The same column with `u^w` independent of `w` has face
    leaving `0` (K8).
- **What must be true.** The specialization imposes only spherical symmetry in the three in-plane coordinates.
  Native fields keep their dependence on radius, `w` and time. Faces depend on radius and time, and so do projected
  objects. No further restriction is imposed. If one were, it would be labelled as supplied and the claim limited
  to it.

**P3. "The order lost" has no sign.** (Grok 3.) **CONFIRMED.**
- **The list.** It prints "the order lost" and "the difference between the face channel and the order lost", with
  no sign stated (K5).
- **The source.** In stage 006, `Gamma_B < 0` is the drain and `Gamma_B > 0` is re-ordering (l.86–87). The S12
  decision list signs `Γ_B` the same way (K5).
- **The evidence.** Grok's sign probe prints `-1` if "order lost" is the source and `1` if it is the decrease
  (K8). The difference flips.
- **What must be true.** The order object over the shell is defined with an explicit sign relative to the committed
  source `n Γ_B` and its face parts, and it states which sign is the drain in stage 006's convention. The printed
  difference uses that definition. No exchange sign is imposed: U2 leaves the sign free.

**P4. The projected velocity divides by a density that can vanish.** (Codex 1, Grok 4.) **CONFIRMED.**
- **The list.** "The projected velocity is the projected material flux divided by the projected material density"
  (K6).
- **The source.** S12's record: "no division by a possibly jumping or vanishing density is performed" (l.161,
  K6).
- **The evidence.**
  - Codex's fixture has a nonnegative density (a square), moving faces and nonzero order flux, and satisfies both
    balances. Its projected density vanishes at a point, where the quotient is `0/0`. The volume flux's limits along
    two approaches differ by `-4/13` (K8).
  - Grok's empty slab (`n = 0`) gives `nan` (K8). Its other witness has a density that changes sign
    (`quotient_n_changes_sign True`, K8). A number density does not change sign, so that witness is outside the
    field class, and the finding does not rest on it.
- **What must be true.** The projected velocity and its volume flux are defined only where the projected density
  is nonzero. That condition is carried wherever they are printed. No value is filled in by continuity. The
  densities, fluxes and balances stay on the full domain.

**P5. The OPEN list assigns the faces' law to S12.** (Codex 3.) **CONFIRMED.**
- **The list.** "the conversion law's forms, the faces' law and the form of `J_χ` (S12, under "Keep S12")" (l.116,
  K7).
- **The source.**
  - The accepted S12 record lists the interior-order law (with its history and strain), the faces' law, and the
    source and flux forms. It says: "No separate owner named for the first two; source/flux forms remain S12's.
    Slab material completion remains separately S8/S22" (l.292, K7).
  - "Keep S12" keeps "every owner the sources and you gave" as written (l.15–18, K7).
  - S12 retains the drain/return functions, their boundary data and the throat's order-conversion source
    (l.18–20, K7).
  - The S12 decision list puts the slab's material law with S8/S22 (l.140, K7).
- **What must be true.** The OPEN list carries owners exactly as the accepted S12 record gives them:
  - source and flux forms: S12's, under "Keep S12";
  - the interior-order law and the faces' law: no separate owner named;
  - slab material completion: S8/S22.

  The list assigns no owner of its own.

## The leg disagreement

On P5, Grok found that ownership held, reading the faces' law as S12's. K7 bears out Codex: the accepted record
names no owner for the faces' law. Codex did not raise P1 or P3. Both legs raised P2 and P4, each with its own
witness.

## Adjudication

P1 to P4 change what is computed, and P5 changes what may be claimed. All five are folded once into the list (G2).
The fold changes the named object, so the folded content is reviewed until clear: a scoped Codex + Grok pass on the
folded content runs before any builder.
