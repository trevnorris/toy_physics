I'll read the packet index and the options memo first, then check the source and saved evidence those documents point to.The native functions show a separate centre coordinate and a thickness-only kinetic caller. I’ll read the governing text next to see what that does to the retained scattering observable.**A-minimal, and then stop.** The retained driven observable does not need a centre equation, a centre solver, or a uniqueness theorem. It needs one written claim boundary that uses the force cancellation and the bare parity record at their actual strength.

The preferred route is the current specification, not a model change. Rank: **A-minimal**, then a **D-style limit** only on the two sentences the evidence does not support, then **B** only if the experiment is later redefined, and **C** last.

| Route | What it is for | What it would add | Place |
|---|---|---|---|
| **A-minimal** | Say what the closed five-field solution already is | A source disposition and a claim sentence. No new algebra | Do this, then stop |
| **D-limit** | Keep two sentences unclaimed | No construction | Attach to A. Do not replace A9 with a reduced handoff |
| **B** | Prescribe `ζ_c`, including a fixed centre | A different experiment, plus holding or reaction accounting if power is claimed | Only by an explicit scope revision |
| **C** | Add centre inertia, stiffness, or a new constraint | New physics and a new response channel | Not supported by an absent code row |

## What is actually established

Thickness and centre are independent geometric coordinates. In `a-shared-physics.md` §1a, `δW ≡ ζ₊ − ζ₋` and `ζ_c ≡ (ζ₊ + ζ₋)/2`, and no S11c-a computation may set `ζ_c = 0`. The same prohibition is repeated in `b-shared-physics.md` §1a and `c1-shared-physics.md` §1a, with one named exception: a centre-fixed uniform regression. Native `DOFS = ("DELTA_W", "ZETA_C")` and `dof_fields` implement that split. For `DELTA_W`, the centre slot is set to zero inside that one source direction (`build_material_face_source`). That is a direction, not an elimination.

The balance that is actually solved is narrower. `S11b-inherited-physics.md` §6 names the state variables `u`, `δW`, and `θ`, and its displayed virtual work tests `δ_v(δW)` only. `b-shared-physics.md` §1a distinguishes those internal fields from the face pair `{ζ₊, ζ₋}`. Section §3a builds the stored-energy basis from `{u, ∇u, θ, ∇θ, e_W, ∇e_W}` — one thickness coordinate — plus background jets. Section §3b asks for equations of `{u, θ, e_W}` only. The supplied kinetic energy, repeated from S11b §5 through a §1c and b §1c, is `½ ρ_br⁰ |∂_t u|² + ½ μ_W (∂_t δW)²`. Native `kinetic_balance_from_energy` uses `u_t` and `e_t` only. That is a supplied field list, not a missing coefficient and not a discovered massless law. `constraint_fold_from_source` likewise eliminates `θ` onto `δ_v u` and `δ_v e_W`.

c2 closes that operator. `c2-shared-physics.md` §1c describes `slab_operator` over `{u, θ, e_W}`, with the θ-row velocity contribution carried by `W_0 e_{W,t}`. Native `REPRESENTATION = 'DELTA_W'`, and `expanded_rows` / `build_case` close `U`, `THETA`, and `E_W`. c1’s response is a map from prescribed `V_s` and `μ_θ` (`c1-shared-physics.md` §0; native `response_operator_case`). Every saved coefficient input in `saved-navigation.json` is `Tuple(Str('V_S'), Str('MU_THETA'))`. Nothing there selects `ζ_c_t`.

Two checks constrain how far that specialization can be pushed.

The completed centre-load diagnostic (`diagnostic-report.md`) finds, in all four anchoring/density cases,

```text
Fcentre = (rplus+rminus) S + (rplus-rminus) D = (−4/W_0) D = 0
```

with `rplus = −2/W_0`, `rminus = 2/W_0`, and the zero inherited from the saved c2 half-difference `D = 0`. `face_generalized_force_rows` is where that centre row is extracted, beside `U` and `E_W`, from thickness virtual work tested on `ZETA_C`. This says the recorded thickness-driven closure has no centre-force obstruction. It does not select `ζ_c`, and the report says so.

The saved first-shape bare DtN parity record has `OFF_DIAGONAL_BLOCKS` equal to `Integer(0)` for both `LAB_HELD` and `MATERIAL_ADVECTED`. Native assembly puts the face kernels in `FACE_BASIS` and the transformed blocks in `THICKNESS_CENTRE_BASIS`. Vanishing off-diagonals mean those two face kernels agree, so that bare operator is diagonal in the thickness/centre basis. `c2-shared-physics.md` §1b records cross-engine agreement of “the parity matrix (= the kernel)” and, separately, leaves the whole-form `dtn_operator` undecided. The navigation limits are right: these zeros are not a slab, constraint, or current-sector theorem. The centre diagonal block is not in the supplied views, so its invertibility is unknown. The one-face sign control has `is_zero` undecided. The open-work byte equality between `(DELTA_W, ZETA_C)` and `(ZETA_C, ZETA_C)` is indexed text equality, not a closed centre equation.

S11b B0b does ask for `Z` separately on the `δW` and `ζ_c` combinations. That is a bulk-response split. It is not a slab unknown.

## The five questions

**1. The gap is real only for an unqualified two-face claim.** For the retained driven observable it has been overstated. That observable is the profile-conditioned transverse↔thickness response of the reduced closed operator (`d-retained-scattering-context.md` §2–§3a): end-channel `J_H` from the thickness-channel current of that pencil, plus a distinct bulk-depth flux driven by the solved state’s face data (`accepted-scope.md` §2). The solved state has no `ζ_c` unknown. `∂_t ζ_c` is absent because c2 consumes the `DELTA_W` face velocity, not because a centre solution was proved to vanish.

The live ambiguity is the status of `δ_v ζ_c`. Curved `a-shared-physics.md` §3b defines virtual work from the full face map, which includes `ζ_c`, while S11b §6 and b §3b vary only the internal fields. Read as a constraint, the stored centre row must vanish. Read as a probe, it is diagnostic. The computed `Fcentre = 0` satisfies both readings and selects `ζ_c` under neither. c1 §1b’s “no incoming waves” fixes the outgoing DtN branch. It is not centre initial data. The formal noninvertibility condition, `S11CC1ZeroInSpectrum` of the per-face resolvent, is not a centre uniqueness theorem.

**2. The supplied laws do not contain a centre balance to finish.** There is no inherited centre inertia or stiffness to type. The smallest missing input, if a later claim needs it, is a selection premise: the retained incident data include no independent centre drive, and the reported particular solution is the existing `DELTA_W` closure. A homogeneous centre addition would be another sector. Its exclusion is that premise, not a spectral proof. I would not invent the centre diagonal coefficient.

**3. Next work is a disposition, not a calculation.** Record the field list, the `DELTA_W` closure, `Fcentre = 0` with its `D = 0` provenance and the surviving c1/c2 debts, and the bare off-diagonal zeros. State that “`ζ_c` eliminated” and “both face drives determined for arbitrary centre motion” are not claimed. Then stop. No repeat of the centre diagnostic, no new mode in d, no radiating Green kernel.

I would change this if a later claim requires arbitrary independent `ζ_c`, or if an inspection of the closed rows — not done here — found a live `ζ_c_t` slot. c2 §1c’s description of the θ-row velocity term as `e_{W,t}` is against that, but this packet does not re-expand every closed term. The saved `MATERIAL_ADVECTED` `ZETA_C` face velocity does contain `u_t` profile terms; that is the kinematics of the centre source direction, and the four-case force check already returned zero on the closed thickness state.

**4. Necessary evidence is already saved.** Optional, and only if someone insists on the unqualified bulk-drive sentence: one reuse of the saved closed assembly with the existing `ZETA_C` face-velocity direction as a source perturbation, asking whether the five closed rows move. That is the missing back-reaction direction. It is not a spectrum, a radiation method, or a second derivation of `Fcentre`. End-channel `J_H` does not need it.

**5. Effects.** The four real-frequency cases are not newly defective. They remain the `DELTA_W` closure on the decaying-bulk domain, with the inherited c1 debts in `accepted-scope.md` §2 and the c2 sign conventions in `c2-shared-physics.md` §1a. The weak conversion FORM and end-channel `J_H` can be built from that reduced operator. Bulk-depth flux may use the thickness-direction drive only with the label that centre motion was not eliminated; `Fcentre = 0` is not that exterior map. A9 is not discharged, and it is not blocked as a whole: the centre issue qualifies the two-face sentence, while A11/A12 and nonempty radiating support stay blocked by `radiation-boundary-finding.md` item 2. Completing centre motion would not supply continuous-spectrum coverage. S11c-e may consume the supported thickness-sector coefficients inside their domain. It must not read them as a centre-resolved leakage total. Pole deferral is untouched.

## Objection, and the line between spec and revision

The strongest objection is that one-way compatibility plus a bare acoustic diagonal does not prove the closed slab rows are invariant under a centre velocity. Profile, constraint, and `Λ_X` traction could still mix the sectors. If they do, the five-field results are the `ζ_c`-absent specialization, and an invariant-sector sentence would be too strong. That objection limits the wording. It does not create a centre unknown inside the operator d already consumes, and it does not justify C.

A is the current-spec route: keep `ζ_c` geometrically independent, do not set it to zero, do not add a field. B or C would revise the experiment or the energy field list. The fixed-centre choice is allowed only as that revision, never as a derived outcome of this diagnostic.