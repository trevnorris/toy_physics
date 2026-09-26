# Centre motion (S11c-d): my recommendation

## Bottom line

I recommend **A-minimal done as a reflection-covariance lemma, plus a short scope note that confines the claims to the reflection-even sector.** The note is D-shaped but is not a downgrade. The centre gap is real for independent two-face drives, centre-sector channels and initial-value problems. It is not a real gap for the retained thickness-driven observables. Those observables don't need a centre equation, a centre inertia value, a centre=0 constraint or an odd-sector uniqueness theorem. They need one symmetry premise, which should be demonstrated from the source rather than assumed.

## Q1: Is the gap real for the retained observable?

Only partly. Take the mirror map M: w → −w, which swaps the two faces. Under M, ζ_c ≡ (ζ₊+ζ₋)/2 is odd, while δW, θ and u are even (b-shared-physics.md §1a lines 62–67).

The supplied inputs look M-covariant:
- **Face maps:** `R_s = (x, sW_bg/2 + ζ_s)` (a-shared-physics.md §3a).
- **Closure:** one set of `Λ_I, τ_I` shared by both faces (a §1c lines 115–119).
- **Normals:** the outward normal's in-plane tilt is the same on both faces (c1 §2a lines 191–193). That is exactly what the mirror predicts.
- **Bulk:** the two half-spaces are disconnected and carry the same radiation condition (S11b §1 line 37).
- **Energy basis:** reflection stays unbroken (b §1c lines 126–127).

If the assembled linear operator commutes with M, it splits exactly into even and odd blocks, with zero coupling in both directions. Then:

- **Existence:** For an even drive (thickness or transverse incident, as in C_{T→H}), the five-field solution with ζ_c ≡ 0 solves the full system. This holds for any M-respecting centre mechanics, whether massless, inertial or flexural. The completed diagnostic already supplies the needed half: the even-to-odd load `Fcentre = (−4/W₀)·D = 0` in all four cases (diagnostic-report.md).
- **Selection:** Even observables don't depend on any homogeneous odd solution when the odd-to-even block is zero. The retained claims therefore don't need odd-sector invertibility or a spectrum.

What the evidence actually establishes:
- **The centre-load zero** establishes even→odd decoupling, but only for the implemented chain and only by inheriting c2's `D=0`. It is not a second acoustic derivation. The one-face sign control shows the route could register an asymmetry.
- **The bare-DtN parity zeros** (saved-navigation.json, `PY_S11CC1_DTN_BY_PARITY` `OFF_DIAGONAL_BLOCKS` = `Integer(0)`, both anchorings) establish both-direction decoupling for the bare bulk operator at first shape order only.
- **Nothing implemented** yet shows odd→even decoupling through the slab, closure and virtual-work rows. b's `face_generalized_force_rows` only builds `CENTER_FACE_GENERALIZED_ROW` from the `(branch,"DELTA_W","ZETA_C")` work. c2's `build_face` feeds only the thickness-keyed face velocity (line 88).

That missing direction is the actual premise. It comes from a symmetry of the source definitions, not from new mechanics.

## Q2: Do the supplied laws already give a centre balance?

Literally, they do. a supplies a centre virtual direction (a §1c lines 152–159; `virtual_work_cases`), and neither the supplied T nor U contains ζ_c (b §1c lines 104–108; `kinetic_balance_from_energy` has velocities `(*u_t, e_t)` only). The balance-law method then gives a massless centre row: net centre face load = 0.

The provenance makes that reading doubtful as an intended choice. S11b was centre-fixed: `ζ_± = ±δW/2` and `δ_vx_± = n̂_± δ_v(δW)/2` (S11b-inherited-physics.md §6 lines 265–270). ζ_c was later added to a as geometry, without any mechanics. The basis enumerates `{u,∇u,θ,∇θ,e_W,∇e_W}` only (b §3a), so a symmetry-allowed centre-gradient term was never considered. Neither was the translational inertia that the same slab mass carries in-plane.

So the packet does not settle whether the massless centre is a deliberate model or an artefact of history. If a general centre sector is ever wanted, the smallest missing input is a **declaration**: centre inertia, and whether a centre-gradient stiffness joins the basis. I won't type a coefficient. The retained claims don't need that declaration if the lemma holds.

## Ranking

| Rank | Option | Verdict for retained claims | Cost |
|---|---|---|---|
| 1 | **A-minimal (symmetry lemma) + even-sector scope note** | Closes centre selection for C_{T→H}, even bulk-depth flux and weak FORM; invariant under any M-respecting centre mechanics | Small: a source argument; optionally one read of existing rows |
| 2 | D (conditional handoff) | Correct fallback if the lemma fails on a named premise | Export and review work remains |
| 3 | B (prescribed centre) | For even drives it gives the same answer with zero holding reaction, so it adds nothing but a scope revision; wrongly implies centre=0 matters | Small, but changes the question |
| 4 | C (centre mechanics) | Needed only for odd-sector or general-drive claims; unjustified by a missing code row | Large, uncertain |

Codex's route is right in spirit but frames the question too broadly. "Do the supplied laws determine the centre balance?" is the wrong first question for the retained observable. The answer to it is ambiguous, and the retained claims don't depend on it.

## Next work and stop condition (Q3, Q4)

**Necessary work:**

1. **Write the M-covariance lemma from the supplied definitions**, citing each premise:
   - face maps and heights (a §3a);
   - outward velocity, flux and traction definitions (a §3b);
   - face-independent closure constants (a §1c);
   - density representatives that depend only on in-plane position (a §3c line 390);
   - identical radiation conditions (S11b §1; c1 §2a);
   - the M-invariant basis rule (b §1c).

   It must cover both anchorings and both density representatives. It must note that the deliberately one-sided controls (b §5a upper-face jet reversal) break M by design and are excluded.
2. **Check the lemma against the saved results without recomputing anything:**
   - `D=0` in the four cases;
   - the bare-DtN off-diagonal zeros;
   - the face-velocity literals. Under the `ZETA_C` key, `ζ_c,t` enters with opposite outward signs on the two faces (odd), as it should. The MATERIAL_ADVECTED `−½σ_W u_t·∇w₁` term under the same key has the same outward sign on both faces (even). Classify by actual M-parity, not by the key label.
3. **Rewrite the retained A9/A11/A12 statements** as M-even-sector claims, with the face-drive map V₊ = V₋ (outward) and ζ_c ≡ 0 derived as the odd part of an even state. This satisfies accepted-scope.md lines 183–188, which ask for the map "including any elimination of the independent centre displacement".

**Optional robustness:** close the existing open odd→even rows (a's `(…,"ZETA_C",…)` virtual-work text through c1) and confirm the literal zero. That is an implemented receipt for a spec-level lemma, not a prerequisite. Global uniqueness, the odd-sector spectrum and centre dynamics are all unnecessary.

**Stop condition:** stop as soon as any premise fails.

The one ambiguity I found: a §3c (lines 394–401) supplies the window `𝒪(G₊, G₋)` without saying it is symmetric under swapping its two arguments. If it isn't, and its projection terms reach the retained rows, the lemma fails at that point. Report that premise once and move to D. Don't escalate to C.

**What would change my view:** if the S11c-e consumer or A9 truly requires odd-sector channels, "every open channel" in d §S-matrix (d-retained-scattering-context.md line 224), or independent two-face drives. Then B or C becomes a real model decision.

## Strongest objection

accepted-scope.md lines 184–185 say: "A parity assumption … is not a substitute."

My answer is that a covariance lemma derived from the stated definitions, with the one-face control showing the route can detect asymmetry, is a derivation, not an assumption. Still, the scope owner must accept that reading explicitly.

A second objection: the only implemented even→odd evidence inherits c2's `D=0` and its cross-engine debts. The lemma doesn't depend on that implementation. But if the implementation disagreed with the lemma, that would be a demonstrated defect in c2, not a centre-physics gap.

## Impact on existing claims (Q5)

- **Four-case results:** unchanged. No new defect is demonstrated. The c1/c2 fidelity and cross-engine debts carry forward as inherited debt.
- **Weak conversion FORM:** bindable for the even sector, with the scope note.
- **End-channel vs bulk-depth flux:** both are even-sector. The bulk-depth construction gets its per-face drive from the symmetric map. Cross terms between the even and odd sectors don't enter the bilinear flux.
- **Radiation coverage:** unaffected. The decaying-bulk vs radiating obstacle (options.md item 6) stays separate and unresolved, and the lemma does not touch it.
- **S11c-e handoff:** receives even-sector conversion that holds for any reflection-respecting centre mechanics. It must accept that centre-channel conversion, independent-face drives and initial-value problems are excluded.

## Current spec vs scope revision

- **Within the current spec:** the lemma, the even-sector claims and the ζ_c ≡ 0 derived as the odd part.
- **Scope revision:** only B, C, or dropping the full odd-sector channel coverage if the contract really demands it.