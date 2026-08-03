# S11b — where the work stands, 2026-08-03

⚠ **Read `steps/S11bA_interface_response.md` first.** It is the closed result and it supersedes anything
older about the interface.

---

## Status

| | state |
|---|---|
| **S11b-A** — the bulk's response to moving faces | ✅ closed on its own terms, ⚠ but see **THE CLOSURE FINDING** below — its `Z_perm` is **conditional on a closure now believed thermodynamically incomplete** |
| **S11b-B** — the homogeneous assembly | ⛔ **directive at rev 4, REJECTED by both legs.** Rev 5 owed, ⛔ **BLOCKED on a user decision** (see below). ⛔ No build has run |
| **S11b-C** — the non-uniform transverse coupling | **NEW, not started.** Split out of B, see below |

## ⛔⛔ THE CLOSURE FINDING — the most important thing on this page

**Three independent parties, none of whom saw the others' work, reached the same conclusion:** the
interface closure `J_± = Λ_p δp + Λ_V V_±` that S11b-A derived and handed to B is **not
Onsager-admissible**. The conjugate force for interfacial mass transfer is the **chemical-potential jump
across the face**, ⛔ not the bulk pressure alone. ⇒ the closure should read `J_± = Λ_p(δp − μ_θ) + Λ_V V_±`
or equivalent, with `Λ_p`, `Λ_V` tied to **one positive-definite Onsager matrix**.

| who | how they got there |
|---|---|
| Codex, rev-3 directive review | *"passivity is not computable from `J = Λ_pδp + Λ_V V`: it needs the brane–bulk chemical-potential affinity and its reciprocal traction law"* — plus: a face at rest with an internal density perturbation has `δp = V = J = 0`, so **pressure-driven conversion is excluded by construction** |
| route agent A | route (c) still admits growth; `max Im ω = +0.00357` **verified at 30 digits**, 151/400 random draws |
| route agent B | same, `max Im ω = +2.84` over its draws; identifies the sink `μ_θ𝒥` as **sign-indefinite** |

⭐⭐ **The consequence that matters: the instability channel cannot be read until this is settled.** §0 of the
directive now (correctly) makes a growing root a first-class falsifier — but with a non-Onsager closure a
growing root is **unreadable**, because it could be the model failing *or* the closure being inadmissible.
⛔ **Opening the falsification channel and leaving the closure unfixed makes the channel worthless.**

⚠⚠ **NOT established, ⛔ do not assume it:** neither agent verified that the corrected closure actually
**restores** passivity. Both asserted it would. ⇒ ⭐ **that must be COMPUTED by B, never assumed** — and
note that "fix the closure until the instability goes away" would be exactly the rescue the charter forbids.
⭐ The legitimate move is the reverse: impose Onsager because the **second law** requires it, then report
whatever stability follows.

⚠ **This does not make S11b-A *wrong*.** A's `Z_perm` is correct **given** its closure, and A's record
already says `Λ_p⁰`, `Λ_V⁰`, `τ` are coefficients of a law B would assemble. But a closure whose flux
depends on `μ_θ` couples the face response to a **brane** variable, so `Z_perm` would no longer be a pure
bulk-response function. ⇒ **A's deliverable changes shape if the closure changes.**

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

## B's REV-3 REVIEW — ✅ **ALL FOLDED INTO REV 4.** Kept for the diagnosis, ⛔ not as open work

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

## ⛔⛔ B's REV-4 REVIEW — REJECTED, and ⭐⭐ THE PATTERN IS THE POINT

Rev 4 (`40086811`, shared block `sha256 54e53ade`) made surgical fixes to every rev-3 finding. **Both legs
rejected it, and every single finding was NEW-IN-REV-4 — in material the fixes had changed.** ⛔ Codex:
*"No PRE-EXISTING finding survived the physics filter."*

⇒ ⭐⭐ **THREE CONSECUTIVE REVISIONS HAVE BRED NEW DEFECTS IN EXACTLY THE PARTS CHANGED TO FIX THE LAST
ROUND.** ⚠ But the defects are **not** scattered — across rev 2, 3 and 4 the blockers cluster in **two
subsystems only**: the **derivation route** and the **complex-frequency continuation**. The brane physics —
`B1`, the energy, the tasks, the controls — has been converging and is now largely clean.
⇒ ⭐ **That is a signal about WHERE the spec is weak, ⛔ not a reason to rewrite it again.** Both subsystems
were therefore extracted and settled as **standalone physics questions** before rev 5 is written.

### What rev 4 got wrong (fold into rev 5)

**Codex — 2 BLOCKERs:**
- ⛔ **Removing the branch rule left §1 UNDER-DETERMINED.** Requirements 1–3 fix the retarded germ in the
  upper half-plane and its real-axis boundary, ⛔ **not** the continuation to B5's lower-half-plane poles;
  the branch points sit **on** the real axis, so descent paths with different winding reach different
  sheets. ⇒ **SETTLED — see below.**
- ⛔ **The "`J_±` is determined" declaration broke the variational route.** One engine varies `J_±[δW,δp]`
  inside the multiplier term and gets `λ δJ` interface forces; another treats it as an external source and
  omits them. Both obey the directive; the dispersions differ. ⇒ **SETTLED — see below.**

**Grok — 4, of which two are mine to own:**
- ⛔ **The passivity block ASSERTS THE ANSWER.** I wrote *"passivity is not computable"* — but `Re Z_perm ≥ 0`
  **is** an explicit inequality on the supplied symbols. ⚠ I collapsed two different questions (face-port
  dissipativity vs thermodynamic admissibility of the closure) and answered both. ⭐ Same defect class as
  the `B6` note I deleted one section earlier. ⇒ ask for **both**, assert **neither**.
- ⛔⛔ **My "hypothesis to test" CANNOT FAIL.** I mandated the prescription, then had engines enumerate loss
  channels *from the equations that prescription produces* and check power balance — an identity for any
  system built that way. ⇒ **the `x === x` defect from A's repair pass, reproduced by me in a directive.**
- ⛔ The cross-check says *"decay rate"*, which fights the growth channel opened two sections earlier.
- ⛔ *"redundant modulo B1"* is ill-defined (also Codex): `B1` is sourced, memory-carrying, and changes rank
  at `ω = 0`. ⇒ **judge independence as FIELD BILINEARS with no constraint applied**, carry every
  independent invariant symbolically, and report separately which become redundant once the constraint is
  applied — and whether that differs between the impermeable and flux-on cases.

**Codex — also:** the §3 symmetry statement is incomplete (no translation invariance, no in-plane parity, no
equivalence modulo total divergences) ⇒ a literal engine adds a **pinning term** `½K|u|²` and **gaps the
modes**. ⇒ state the symmetry group completely.

## ⭐⭐ RESOLVED SEPARATELY — two subsystems, four independent agents, SymPy/numerics

### 1 · The branch — ⛔ Codex was right, Grok's clear was wrong

**Both branch agents independently returned CLAIM B: under-determined**, with explicit numbers. Two paths
from the same upper-half-plane start to the same lower-half-plane point give `q_out` differing by **exactly
`−1`** — one a decaying **normal mode**, the other a growing **leaky resonance**, same `ω`. Fed to a toy
radiative-loading secular equation the two sheets give `Im ω = −0.0745` vs `+0.0745`: ⭐ **the sign of the
deliverable is set by the branch, not by the requirements.** Disagreement is **200 % of the deliverable at
every damping strength**, and ⛔ **does not shrink as poles approach the real axis** — a straight-line root
search can cross the cut and swap sheets **silently, with no numerical warning**.

⭐ One agent reproduced Codex's **rev-3** finding without being told it existed: `Im q_out ≥ 0` for all
complex `ω` is non-analytic and disagrees exactly where it matters.

⭐⭐ **The prescription for rev 5** (both agents converged; second agent's clause 3 is load-bearing):
1. Requirements 1–2 fix `q_out` **on the real axis only**. ⚠ Requirement 1 must be read as an **energy-flux**
   condition, valid for **both signs of `ω`** — ⛔ not a phase-velocity condition, which breaks at `ω < 0`.
2. At complex `ω`, `q_out` is **defined** as the continuation of the real-axis branch reached from
   `ω + i0⁺` **downward along the ray of fixed `Re ω`** — equivalently, deform the inverse-Fourier-in-time
   contour downward from above the real axis with the branch points `ω = ±c_s0|k|` fixed and their cuts
   dragged **vertically down**.
3. ⛔⛔ **Requirements 1–2 must NOT be re-imposed at complex `ω`.** Whatever `|w| → ∞` behaviour the
   continuation produces is a **RESULT to report**, ⛔ never a criterion for re-selecting the root.
   ⚠⚠ **An engine that re-applies "must decay" at a complex pole lands on the wrong sheet and converts a
   damped resonance into an APPARENT INSTABILITY** — which, with §0 now declaring growth a first-class
   finding, would manufacture a fake falsifier.
4. Report the degenerate point `ω = ±c_s0|k|` (`q_out = 0`, the two bulk solutions coalesce, the second
   going linear in `w`); ⛔ neither requirement selects anything there — continuity does.
5. If a pole's trajectory crosses `Re ω = ±c_s0|k|` under parameter variation, **report that it has left
   the sheet**, ⛔ do not re-select it onto one.

⭐ **Why supplying a precise rule is right this time, having been wrong in rev 3:** rev 3 **asserted** a
precise rule; rev 5's was **derived and numerically verified twice, independently.** ⇒ the lesson is
**"be precise only where you have verified"**, ⛔ not "never be precise".

### 2 · The derivation route — ⛔ the action principle was the wrong tool

**Both route agents converged**, with SymPy determinants and root-locus numerics:

- ⭐ **Route (c) — write the BALANCE LAWS directly**, using `U` only as the constitutive potential. It is
  **identical** to the Lagrangian route restricted to **Lagrange–d'Alembert virtual work** — `J_±` and `δp`
  held **fixed** as prescribed external sources, bulk load as external virtual work — ⛔ **not** stationary
  action. One agent verified `det[(a)-fixed]/det[(c)]` is a bare monomial, i.e. the same dispersion.
- ⛔⛔ **Varying `J_±` inside the multiplier is INCORRECT — not merely ambiguous.** Variation transposes the
  retarded kernel into an **ADVANCED** one (`Λ⁰/(1−iωτ) → Λ⁰/(1+iωτ)`, pole `−i/τ → +i/τ`): an anti-causal,
  **energy-generating** kernel. ⇒ a spurious growing root at `≈ +i/τ` (`+3.089i` at `τ=1/3`; `+1.517i` at
  `τ=0.7` against `i/τ = 1.4286i`), ⚠ **whose growth rate tends to `1/τ` as `Λ⁰ → 0`** — a finite-rate
  instability from infinitesimal coupling, which is what falsifies it.
- ⛔ **Route (b), substituting the constraint into `U`, is wrong the same way** and **invents a mode**
  (determinant degree **6** vs **5**). One agent finds `(b) ≡ (a)+δJ` identically; the other finds them
  differing at `O(a)` vs `O(ā)`. ⚠ **Minor unresolved discrepancy between the agents — not load-bearing**,
  since both agree `(c)`/`(a)`-fixed is correct and both others are wrong.
- ⛔⛔ **The routes AGREE in the impermeable limit** ⇒ ⚠ **rev 4's `Λ⁰ → 0` consistency check would have
  passed for all three wrong routes.** One agent adds that the split survives at `τ = 0`, so it is caused by
  **dissipation**, ⛔ not by memory.
- ⭐ **The trap inside route (c):** the normal-force balance must include the slab's own internal pressure
  pushing the faces apart (`+μ_θ`, `μ_θ ≡ ∂U/∂θ`). One agent got this wrong first pass and the variational
  routes caught it. ⛔ Omit it and (c) disagrees with everything.

⭐⭐ **THE ACCEPTANCE CHECK TO ADOPT — mechanical, gradeable, and it does NOT foreclose instability:**
**if `Y(−ω)`, `Λ(−ω)` or `Λ*(ω)` appears anywhere in the final equations, the derivation is wrong.**
⛔⛔ **DO NOT adopt the other proposed gate — "every root must satisfy `Im ω ≤ 0`".** It would re-close the
falsification channel §0 just opened. ⚠ Use the **specific** diagnostics instead: a root at `+i/τ` means
someone varied `J`; an unattributable `a·μ_θ` term means someone substituted into `U`.

## ▶ Next actions, in order

0. ⛔⛔ **USER DECISION REQUIRED on the closure** (top of this page) before §2 of rev 5 can be written.
1. **Write rev 5** — fold the six rev-4 findings plus the two resolutions above. ⭐ Replace §3b's
   variational mandate with balance laws; ⛔ delete the un-failable accounting check.
2. **Review rev 5** — I author it ⇒ **Codex + Grok**. ⛔ Neither may be told the resolutions came from
   agents; they must judge the spec as written.
3. **Build B's blind `.wl`**, then quarantine it and build the `.py`. ⛔ Reviews launch **on sight**.
   ⛔⛔ **BUILD OUTSIDE THE REPOSITORY** — see the blindness-sweep result below.
4. **Cross-engine comparison**, pre-registration scoring, step record.
5. **Registry insertion as a separate pass**, once B's physics is banked.
6. **Then S11b-C**, using the four requirements listed above.

## ⛔⛔ BLINDNESS SWEEP — the bar list is INSUFFICIENT

An agent swept the tree for anything that could hand a builder S11b-B's answers. ⚠ **The two strongest
leaks are OUTSIDE the barred tree entirely:** `docs/two_throat_simulation_handoff_spec.md` §2.2 (the
constrained-field elimination, the renormalized inertia, the resulting longitudinal pole — B4 + B5's shape)
and `software/stage1_solver/reports/pathA_36_c5_phase_potential.md` (its source, plus a control table
spanning B1/B4/B5/B6/B8). Also unbarred: `_scratch/` (verbatim mirrors of **every** barred file **plus
prior reviewers' own worked derivations**), the built **PDF** of the paper, and `DEFECT_REGISTER.md`.

⇒ ⛔ **Do not answer this by extending the denylist.** The directive is self-contained by design, so
**run the build in a directory OUTSIDE the repository** and move the deliverable in afterward.
⚠ **Keep the empty output directory in place** — quarantining a whole `scripts/` dir once left a builder
nowhere to write and it exited 0 with the file parked in `/tmp`.
⛔ **Absence does not survive `git show`.** Bar history reads explicitly and audit the build log for them.
