# S11b — where the work stands, 2026-08-03

⚠ **Read `steps/S11bA_interface_response.md` first.** It is the closed result and it supersedes anything
older about the interface.

---

## Status

| | state |
|---|---|
| **S11b-A** — the bulk's response to moving faces | ✅ **CLOSED** — two blind engines, five independent derivations, zero disagreements |
| **S11b-B** — the homogeneous assembly | ⛔ **directive at rev 3, REJECTED by both legs.** Rev 4 owed. ⛔ No build has run |
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
- ⚠ **A derivation prescription is mandated** (Lagrangian + multiplier + Rayleigh dissipation) because the
  constraint is non-holonomic. ⛔ **REJECTED — the Rayleigh term double-counts `Z_perm`'s dissipation.** See
  the rev-3 fold below for the suggested single non-redundant rule.
- ⭐ **A symmetry-allowed cross term `C` is carried**, because a diagonal energy silently imposes
  separability, which is part of what the step must determine.
- ⛔ **The branch rule `Im q_out ≥ 0` is WRONG** — it is not the retarded continuation and would flip the
  sign of the deliverable. ⚠ Rev 4 must not simply pick another rule; make the engines report and justify.
- ⭐ **The registry is barred from BOTH engines**; insertion is a separate later pass.

**Pre-registration:** committed at `a68245b4` and **restored to the tree for the compact**.
⛔⛔ **RE-QUARANTINE IT (move it out of the tree) BEFORE B's FIRST BUILD** — `steps/S11bB_PREREGISTERED_PREDICTION.md`. ⚠ It marks which predictions come from a
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

## ⛔⛔ B's REV-3 REVIEW — BOTH LEGS, FOLD THESE INTO REV 4

⚠⚠ **READ THE DIAGNOSIS BEFORE THE LIST.** Rev 1 and rev 2 failed structurally. **Rev 3's new findings are
all in the parts that were CHANGED to fix rev 2** — each revision was a substantial rewrite, and each
rewrite introduced fresh defects of the same class. ⇒ ⭐⭐ **REV 4 MUST BE MINIMAL AND SURGICAL: fix the
named findings and change nothing else.** ⛔ Do not rewrite the specification again.

⭐ **The sharpest instance of why:** rev 2's branch rule was *vague but safe* ("state which sheet");
rev 3 replaced it with a *precise and wrong* one (`Im q_out ≥ 0`), which is **not** the retarded
continuation and would flip the sign of the deliverable. ⇒ ⭐⭐ **When fixing an ambiguity you are not sure
of, a wrong precise answer is worse than the ambiguity — make the engines report and justify their choice
instead of picking badly.**

### ⭐⭐⭐ The finding that outlives the directive

**The scope sentence "propagate freely, decay, or fail to exist" OMITS INSTABILITY.** The quadratic energy
carries no boundedness conditions on its moduli or the cross term `C`, so an unbounded scalar sector could
produce a **growing mode** — and both builders were steered to misclassify it rather than report it as
⭐ **the charter-level falsifier it would be.** ⛔ A falsification channel was closed by accident, in a
ledger whose entire purpose is being able to fail. **Fix this first.**

### Convergent (both legs)

- ⛔ **`Z_perm` and the Rayleigh function double-count the closure's dissipation.** `Z_perm` already folds
  `Λ_p`, `Λ_V`, `τ` into the pressure–velocity law; a Rayleigh term built from the same `J_±` adds a second
  dissipative force from the same object — and its form is left free, so engines can also *disagree*.
  ⭐ **Suggested single non-redundant rule:** treat `J_±` as **determined** by the closure plus the bulk
  solution rather than as an independent field. The constraint then stops being non-holonomic in effect,
  substitution becomes legitimate, `Z_perm` carries all dissipation as one external face force, and ⛔ **no
  Rayleigh function is needed.** ⚠ Verify that claim rather than assuming it.
- ⛔ **The stored-energy list is not a closed symmetry basis.** Terms involving `∇·u` and mixed
  scalar-gradient stiffnesses are allowed at the same order. ⇒ ask for the **closed basis**, ⛔ do not give
  a list and say "report omissions" — `B8-B` can otherwise certify removal of gradient stiffness while an
  omitted allowed term still changes the dispersion.

### Grok only

- ⛔⛔ **`U` has NO `u`↔compression coupling; the only path is `B1`'s constraint, and that is left VERBAL.**
  ⇒ both engines could write the balance **without a `∇·u` term**, obtain a free longitudinal mode, and
  agree. ⭐ **This is whole-interface rev 1's false null rebuilt inside `B1` — the third time this coupling
  has been left implicit after being derived explicitly in the walk.** ⇒ **WRITE THE CONSTRAINT OUT.** It is
  linearised mass conservation: kinematics, not a result.
- §3b's variational steps never say where the integrated-out bulk enters, while `B2` still requires it.

### Codex only

- ⛔ **`Im q_out ≥ 0` is not the retarded continuation** — see the diagnosis above.
- ⛔ **"Assume passivity" is not computable** from `J = Λ_pδp + Λ_V V`: it needs the brane–bulk
  chemical-potential affinity and its reciprocal traction law, neither supplied. ⚠ A stationary face with an
  internal density perturbation has `δp = V = J = 0` here, so **pressure-driven conversion is excluded by
  construction.**
- ⛔ **Calling `Λ⁰ → 0` "reversible" contradicts** the supplied fact that propagating `Re Z` survives as
  radiation resistance with impermeable faces.
- ⛔ **Control A does not isolate thickness:** holding `δW = 0` also zeroes the face velocity and so removes
  pressure loading and permeability simultaneously.
- ⛔⛔ **The note that B6's controls "could only re-encode a predetermined result" constrains the answer and
  defeats blindness** — it instructs engines to discard ablation dependence rather than discover it. ⭐ Same
  defect class the rule was written to prevent. **Delete it.**
- ⛔ **Normalizations undefined** for the internal-DOF count, `B3`'s response, `B4`'s stress, `B6`'s
  coefficient — engines can emit incomparable dimensions under identical tags.
- ⛔ **`q_n` is never defined** (only `q_out` is), and `v₀|q_n|/ω` is not a real smallness measure for
  complex `ω`.

### Cleared by both legs — ⛔ do not re-open

Scope boundary is explicit and not leaked · header physics symmetric and `reduction/` barred from both ·
no pre-registration leak · `B8` controls B/C/D are genuine form cuts · the closed count is square **given a
correct `B1`** · tractability and orphan check clean.

## ▶ Next actions, in order

1. **Fold rev 3's review** (running at handoff time; outputs in the session scratchpad).
2. **Build B's blind `.wl`**, then quarantine it and build the `.py`. ⛔ Reviews launch **on sight**.
3. **Cross-engine comparison**, pre-registration scoring, step record.
4. **Registry insertion as a separate pass**, once B's physics is banked.
5. **Then S11b-C**, using the four requirements listed above.
