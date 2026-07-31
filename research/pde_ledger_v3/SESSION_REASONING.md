# How a symbol repair became a foundations question — the reasoning chain

**Session 2026-07-30/31.** Written so the *argument* survives, not just the conclusions. Every link has a
locus; ⛔ if you doubt one, open the file rather than re-deriving from here.

---

## 1. The job was "repair the `a`-pin damage." The damage was not what we expected.

Triage rule (from D-01a): **a mention is not damage.** Damage is only (a) the pin presented as
determining a physical radius, (b) pin and throat radius conflated, (c) their dimensional coincidence
recorded as evidence.

Sorting every hit first — before editing anything — collapsed a corpus-wide grep into **17 loci**, and
found that `STATUS.md` was **already repaired** and `stage043`'s engines were **clean**. Two of the ~8
"tier 1" files needed nothing.

⭐ **But the load-bearing form of the damage was not narrative.** In 13 places `parameter_register.md`
uses `a` as a genuine physical radius inside a physics formula — `K̄₀a⁵/(27c_s⁵)`, `z = aω/c_s`,
`g_base ∝ a^(−7/2)`, `(c_s/a)²(6+(ma)²)` — and tags it **`CONV`**, so it is **not counted**. Several sit
directly under headline claims: *"ZERO new counted knobs"*, *"Part-II counted CALIB set UNCHANGED"*.

⇒ **A calibrated physical input was being excused as a units convention. The count was understated.**
The prose was a symptom; the accounting was the disease.

## 2. The check that should have caught it could not see the question

`stage043`'s medium block (`:606-615`) carries 10 symbols and 5 constraints — one of which is
`scale_a*mass*cs0 − hbar`, **the pin relation used as content**. That block is identical to the v2
registry's, which is exactly `acceptance_check.py`'s expected payload.

⇒ **Both sides of the acceptance comparison admitted the pin.** They agreed because they *shared the
assumption under test*. The `MATCH` carried no information about the thing being asked.

⚠ Worse: `stage043` tags the same pin `CAT_CONVENTION` in its register and consumes its equation as
content 270 lines later. **Classified as convention, used as content, in one artifact.**

## 3. Why the pin is not physics — and the sharper reason

`{c_s0, ħ, m_GNLS}` over `{L,T,M}` is `det = 1`: **rank 3, nullity 0**. Already a complete unit system.
So a *fourth* pin must produce exactly one dimensionless monomial, and since `a` is a pure length that
monomial is forced to be `a·m·c_s0/ħ`. **The form was never in question.**

⚠ **Codex's correction, which sharpens this:** rank forces the monomial's *form*, not its *value*.
Setting it to 1 is an extra claim unless `a` is *defined* as the unit length. ⇒ The pin was **not
contentless bookkeeping** — it was a **false claim**: that the throat radius equals the medium's natural
length.

⭐⭐ **And the deepest version, from investigating A7 at the very end:** `ħ/(m c_s0)` **is** a healing
length in a standard convention — that is *why* identifying `a` with "the GNLS healing core" felt
natural and survived. **The pin's number was never wrong. The category was.** A medium-wide constant was
used as a per-defect quantity.

⇒ **The test that catches this class:** ⛔ not *"is this value right?"* but ***"is this the same KIND of
thing — one number for the whole medium, or one per particle?"*** The first question passes the pin. The
second kills it instantly, because `ξ_h` is one number and leptons are a family.

## 4. Which raised the question the session actually turned on

If throat radii are species-indexed, **how does the radius scale with mass?**

**Route A — external.** `k e²/a = m c²` gives `a ∼ r_e ∝ 1/m`. ⛔ Doesn't apply here: it assumes the rest
mass is *entirely* electrostatic self-energy, and in this model charge is the **Z₂ ±w orientation**,
identical across leptons, while mass is drain/wave energy. The spec already de-rates it
(*"never an interior-PDE or force-sign input"*).

**Route B — my error, recorded so it is not revived.** I derived `a ∝ m^{1/3}` from
`F(a) = A/a + B/a² + Ca³` by dropping `B` and keeping `m = (18/11)A/a`. **Wrong**: the `18/11` comes from
the partition `E_w:E_f:E_PV = 11:2:5`, i.e. `E_f/E_w = 2/11`, which **requires `B ≠ 0`**. The partition
also *fixes* `C = 80A⁵/(11⁵B⁴)`, which I had treated as free. Three reviewers, two SymPy-verified.

**Route C — the corpus's own answer.** With `B` kept: `a_* = 11B/(2A)`, `F_* = 36A²/(121B)`, hence
`a ∝ m^{-1/2}` — heavier ⇒ **smaller**.

## 5. ⭐⭐ But route C is dead, and that is the session's real result

The corpus's family route gives `a_j ∝ 1/(2j+1)` **and** `F_j ∝ (2j+1)²`. The second is
**`1 : 9 : 25`** — and `lepton_mass_notes.md:424-473` says of it: *"decisively ruled out"* against
`206.77` and `3477.37`.

⇒ **The slope and the falsified masses are the same derivation.** You cannot keep one and discard the
other — *for the support-only route*.

⛔⛔ **CORRECTED 2026-07-31, and the correction matters more than the claim.** I then wrote *"there is no
surviving mass–radius slope in this model."* **That is false**, and an external review caught it.

`notes/lepton_mass_notes.md` is **3239 lines**. I had read about 135 of them. Section 6 says outright
(`:863`): *"the old support-only `1:9:25` falsifier was **too simple**, because once throughput and
geometry are allowed to respond dynamically, the family ladder changes."* Later sections develop
alternative routes — `φ_j = R_j^{3/2}/√ν_j` (`:740`), the low-harmonic benchmark `Wν = R^{9/5}`
(`:2199`).

⇒ **The honest statement:** *no unique, target-blind, empirically successful mass–radius law survives;
several mutually incompatible conditional slopes remain.* ⚠ And they are weak — `:857` concedes the
turbine route bridges the gap only *"in a regime very close to the critical point `s = 3/4`"* and *"does
not naturally produce"* the observed ratios.

⛔⛔ **The lesson is the real finding.** This is **§9's own failure mode — "asserting absence from a
partial search" — recurring in the same document that names it**, on the session's headline conclusion.
I read 4% of a file and wrote "there is no…". ⇒ ⛔ **A universal negative requires having read the whole
artifact, or it is not a finding — it is an impression.** `wc -l` before concluding absence.

## 6. And killing it removed the model's only family label

The falsified picture — leptons as support-mode excitations `j` of **one** throat — made **charge
universality automatic**: same throat, same puncture, same charge. That is presumably why it was
attractive.

Its death leaves the model needing a family label that **changes mass by 207× while leaving charge
exactly identical**, and it no longer has one.

⇒ ⭐⭐ **The question: what makes a muon a muon?** Everything else — `L/a`, `R*(E)`, the geon profile —
is downstream of it.

**And the corpus already records the tension.** The conceptual story says charge is pure direction with
universal magnitude (*"Charge ⊥ mass … which is why the charge is the same regardless of mass"*), but
the **computed** Coulomb magnitude is geometric — `3Q_E²ℓ tanh²(b/ℓ)/(8πRb²)` with `Q_E` `[CALIBRATED]`
— and the manifest carries `q_h/Q_E universality is required but not earned from the current b/ell
ledger` → **`FALSIFIABLE_TENSION`**.

⭐ **This is leverage, not just a problem.** Charge universality is a **held-out, target-blind fact**:
nothing is fitted to it. So it *constrains* how sleeve geometry scales with mass. A constant `L/a` makes
every sleeve length scale with particle size; universality demands the charge-carrying combination
**not** scale. That is a constraint equation on `L(a)`, not a free ratio — and it is a fourth independent
reason `L/a = 37/20` is unsupported (the others: admitted ansatz, adjudicated `free_choice`, contradicted
by the solved `0.9575`).

## 7. Why none of this was visible: the year built linear response

`two_throat_simulation_handoff_spec.md:324` lists what is committed — the bulk Madelung system and
*"the **quadratic** brane Lagrangian and its transverse/longitudinal sectors"*.

**Quadratic is the whole story.** It is small-oscillation physics about an assumed equilibrium — exactly
right for waves and far-field forces, and **structurally incapable** of describing a defect, which is a
large-amplitude nonlinear configuration where the order parameter goes to zero.

⇒ The year didn't fail to solve the brane–bulk. It solved the **linearised** one. Every sector's
unresolved item points at the same missing interior, which is why *one* throat solve is the shared R1
for gravity, electric sign and magnetism at once.

⚠ **Corollary that matters for v3:** gravity is clean **because it is insulated from the hard part**.
`4d_2_5pn.tex:613` — *"controlled by the worldline/worldtube multipoles rather than by arbitrary internal
details of the defect."* A tidy v3-gravity ledger is **not** evidence the model is sound.

## 8. What is actually missing, and two corrections to earlier belief

1. ⛔ **The geon — the mass carrier — has no equation.** *"its profile is a **declared OPEN input**
   [POSTULATE: mass mechanism; profile OPEN]"*. Not deferred. **Undeclared.**
2. ⛔ **No closed parent action.** The BVP is a *"template, not yet a closed BVP"*; well-posedness
   `[OPEN]`; status `UNRESOLVED(closure)`.
3. **Three unlinked mass-ish quantities**: `J` (the drain rate, gravity's source), the **geon** (what
   constitutes mass), and `m_defect` (a `GAP`, `INFLOW_MASS_SOURCE_MISSING`, only a dimensional bridge).

⚠ **Correction 1 — `V_χ` is NOT an unfilled slot.** I claimed it was and proposed "write the wall
potential" as the restart. **Wrong.** It is written with coefficients
(`g0_closure_card_v0.md:115-124`) and the kink **is** solved (`stage006`: `δ`, `σ_wall` derived). The
real problems are subtler: those rest on two *postulated* constants `{a_B, κ_B}`; `W_slab` is
`FREE-UNREDUCED` (*"double-well selects NO width"*); and the polar-smectic program that would have
grounded the medium returned **`FAIL_COUPLE_STRESS_NOGO`**, marked *"SUPERSEDED at the brane-existence
level"*.

⚠ **Correction 2 — `THROAT_DRAIN_DESTABILIZED` does not say the particle is dead.** I cited
`pathA_26_derrick.md:58-72` and missed `:101-123`, headed *"NOT a computed result; read before citing the
top-line"*, which states outright: *"must NOT be read as 'the gravitating particle is killed.'"*
Conservative existence is **generic** (interior positive-definite minimum in **75/75** scans); the
instability sits above `gcrit ≈ 0.006` while the tested box demanded `g ≈ 89`, a synthetic worst case.
⭐ **The real gap is that `g_phys` was never mapped** — and, deeper, that with `J_w = 0` on the exact
trapped mode, *whether "gravity = drain" is even nonzero* is unbuilt.

## 9. The pattern under all of it

**Ten pin-shaped identifications** found (register §A): two quantities sharing a dimension, silently
equated. `a`-pin · `ℓ = δ` · `L_W := L` · implicit `β·a = 1` · `c_w ≡ c_s` · the length chain ·
`r_e` body-vs-mouth · `W_slab` vs `L/a` · and one of mine (sonic horizon ≡ mouth radius, retracted
same session).

⇒ ⭐ **This is not a bug that recurs. It is *the* systematic failure mode of the project**, and it is why
`a` had a plausible value: every one of these is two *real* lengths, correctly computed, wrongly
identified. ⛔ Assume one is present in any sector not yet checked with the **kind-test** (§3).
