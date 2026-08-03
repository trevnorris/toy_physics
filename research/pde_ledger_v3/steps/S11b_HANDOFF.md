# S11b — where the work stands, 2026-08-03

⚠ **Read `steps/S11bA_interface_response.md` first.** It is the closed result and it supersedes anything
older about the interface.

---

## Status

| | state |
|---|---|
| **S11b-A** — the bulk's response to moving faces | ✅ **CLOSED** — two blind engines, five independent derivations, zero disagreements |
| **S11b-B** — the homogeneous assembly | **directive at rev 3, in review.** ⛔ No build has run |
| **S11b-C** — the non-uniform transverse coupling | **NEW, not started.** Split out of B, see below |

## ⭐⭐ What S11b-A established, in one place

```
q² = ω²/c_s0² − k²          Z = ω ρ_m / q_out
Z_perm = (ρ_m r + Λ_V⁰)/(y r − Λ_p⁰) ,   r = 1 − iωτ ,   y = q_out/ω
m_add  = ρ_m/α per face (evanescent, outward)
```

**The headline, and it was not expected:** evanescent dissipation is created **entirely by the interface
closure** — propagating loss is radiation resistance and needs no permeability. Within it, the **velocity**
coupling `Λ_V` dissipates **iff `ωτ ≠ 0`**: a velocity-coupled flux does no net work when it responds
instantly because it stays in quadrature; give the conversion a finite rate and it lags, and the lag does
the work. The channel **freezes out** as `ωτ → ∞`. ⇒ **slow disturbances see the leak; fast ones do not.**

⚠ **That result exists only because the user objected to a memoryless closure.** Under an instantaneous law
the term is identically zero and the pre-registration's `P-B` would have scored as **confirmed**, banking a
false mechanism into the step where the dissipation question gets decided.

## Why B was split, twice

**Whole-interface rev 1 and rev 2 were both rejected before any build.** Rev 1 specified a system with **no
linear coupling at all** — both engines would have returned a block-diagonal result and concluded *"the
longitudinal does not couple,"* an artifact invisible to their agreement. ⇒ split into A (bulk response)
and B (assembly). **A closed cleanly**, which is the evidence the cut was right.

**B rev 1 was then rejected** for mandating a non-uniform background while fixing every perturbation to a
plane wave — position-dependent coefficients mix wavevectors, so no global dispersion exists. ⇒ split again:
**B is the homogeneous assembly** (does the longitudinal radiate or stay bound — a question needing no
gradients), **C is the non-uniform coupling** (is light's confinement unconditional — which genuinely needs
them).

⭐ **The seam is real, not a retreat.** Each split turned structural failures into local ones.

## ⛔ Three traps S11b-A measured, carried into B and C

1. The per-face inertial loading is `+ρ_m/α` **against the outward acceleration on both faces**. The signed
   pair `(ρ_m/α, −ρ_m/α)` is a **convention artifact**, ⛔ not negative inertia.
2. ⛔ **Propagating `Re Z` is radiation resistance**, present with impermeable faces. It is ⛔ **not**
   evidence of transfer through the interface. Only **evanescent** `Re Z` is created by the closure.
3. The bulk's rest-frame linearisation discards a relative correction of order **`v₀|q_n|/ω`**, ⛔ not
   `v₀/c_s0`, and it **exceeds first order whenever `k c_s0 ≫ ω`** — the regime carrying the added mass and
   the evanescent dissipation.

## ⭐⭐ What S11b-C must handle that B does not

Both reviewers named these; ⛔ none is optional:
- **Tilted faces.** With `W₀(x)` varying, an in-plane displacement across `∇W₀` produces **normal** face
  motion of the same gradient order as the coupling being computed.
- **Eulerian vs material density.** The two differ by an advective `u·∇ρ_4D` term that ⭐ *directly changes*
  the transverse–thickness coupling.
- **Plane waves are not eigenmodes.** Either treat the inhomogeneity as a perturbation about a uniform state
  with mode conversion as an off-diagonal coupling, or say what restriction is being granted. ⛔ Do not
  demand a global dispersion relation.
- ⚠ **A uniform-limit control is empty unless the coupling can have a gradient-independent part.** If the
  coupling is *defined* from gradients, setting them to zero passes by construction.

## Where B's directive stands

`directives/S11bB_{SHARED_PHYSICS,wl_header,py_header,wl_directive,py_directive}.md`, rev 3,
`sha256 2e80bb49…`. Two rejections folded. Its distinctive features:

- ⭐ **No single in-plane compression modulus is supplied.** Compression is carried by `θ` and `e_W`, and
  where a frozen-thickness modulus sits is an **output**. ⇒ the double-count that killed an earlier revision
  is structurally impossible.
- ⭐ **A derivation prescription is mandated** (Lagrangian + multiplier + Rayleigh dissipation), because the
  constraint is **non-holonomic** — it has a source and memory — and substituting it into the energy gives
  a *different answer* from local stress plus continuity.
- ⭐ **A symmetry-allowed cross term `C` is carried**, because a diagonal energy silently imposes
  separability, which is part of what the step must determine.
- ⭐ **The branch is defined constructively** (`Im q_out ≥ 0`), because `Im ω` is the deliverable and two
  sheets give opposite signs.
- ⭐ **The registry is barred from BOTH engines**; insertion is a separate later pass.

**Pre-registration is committed and quarantined:** `a68245b4`. ⚠ It marks which predictions come from a
solid derivation and which from a sketch, and states that the transverse coupling's exponent is
**deliberately withheld** — with the note that the author privately expects `|g|²` and ⛔ must not read that
outcome as confirmation of a derivation never performed.

## ⚠ Known limits, ⛔ recorded not fixed — red-team items

- **19 of 23 gated tags in S11b-A's `.wl` cannot fail** under upstream corruption. Its `PASS` certifies less
  than it appears to. ⛔ **Not evidence against the values** — five independent derivations back them.
  ⚠ A repair pass attempted this and **made it worse**: it added 183 lines whose checks were `x === x`,
  because the repair directive mandated *"a check that reads the same expression that is emitted"* while
  **banning independent re-derivation**. ⭐ The gates that survive real corruption are exactly those
  comparing against a **physics-independent constant**.
- **The convective-order tag text still reads `O(v₀/c_s0)`** in both S11b-A scripts. Corrected in the step
  record; the script text is a red-team item.
- One engine renamed a tag against the spec (`RELATIVE_FLUX_CHANNELS` → `..._COMBINATIONS`); values agree.

## ⭐⭐ Process lessons banked today

- **`feedback_static_or_instantaneous_check`** — a standing check now in `/step-run`: *what timescale or
  rate did I just send to 0 or ∞, and what would it have governed?* ⛔ The limit is illegitimate when the
  removed rate is close to what the step is trying to determine. ⚠ It never looks like removing time; it
  looks like a simple constitutive law.
- **Reviewer set follows authorship** — orchestrator-written artifacts get **Codex + Grok**; Codex-written
  ones get **a fresh Claude agent + Grok**.
- **Blindness by absence, not instruction** — raw transcripts now write outside the repo; answer-bearing
  documents are moved out of the tree during a build.
- ⛔ **Serialize review legs that ablate Mathematica** — a 2-seat licence; one leg died mid-run at exit 144.
- ⛔ **A finding about a CHECK is not a finding about the PHYSICS**, and there is a **red-team phase** for
  hardening. ⭐ When tempted to harden one engine's internal checks, **build the other engine instead.**

## ▶ Next actions, in order

1. **Fold rev 3's review** (running at handoff time; outputs in the session scratchpad).
2. **Build B's blind `.wl`**, then quarantine it and build the `.py`. ⛔ Reviews launch **on sight**.
3. **Cross-engine comparison**, pre-registration scoring, step record.
4. **Registry insertion as a separate pass**, once B's physics is banked.
5. **Then S11b-C**, using the four requirements listed above.
