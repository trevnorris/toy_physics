# Decision list — the unified S11b step ("the interface coupling law")

**Orchestrator-written, 2026-08-19.** Settles the structure, naming, provenance and scope questions the
S11b **shared spec** and build directive must honour. ⛔ **Not itself a build directive**, and ⛔ not the
shared spec. It gets **two legs (Codex + Grok)** before any builder launches (rule 7); one pass, then
fold and go.

⭐ IDs are the **G-series** ("G" = the unified S11b step). The established export-chain rules **`F1`–`F9`**
(`S11_export_chain_decisions_v2.md`) and the three script clauses (`.claude/skills/build/SKILL.md`) are
**inherited, not restated** — a G-item points at them.

⚠ Two facts frame everything: S11b was two execution stages (A = the bulk's response to moving faces,
B = the homogeneous interface assembly), built on the **old pre-chain pattern**; and the A→B transition
is a **generalisation**, established below by reading the specs, ⛔ not a rename.

---

## A · Structure

**G1 · S11b is ONE step on the export chain, subsuming the historical A and B stages.**
⭐ One shared spec, two blind engines, one `S11b_exports.py`, one comparator, one step record, one card.
The A/B split was **execution history, not ledger structure** (`steps/S11b_HANDOFF.md:9-11`); A earned
**no ledger rows by design** (`steps/S11bA_interface_response.md:198`). ⛔ S11b-C is a **separate later
step**, not in this build.

**G2 · The unified step is B-LED.** The canonical constitutive law is **B's affinity closure**; A's
regime is recovered as an **explicit reduction slice** (see `G4`,`G5`), ⛔ never discarded. B's law is a
generalisation that **must reproduce A's `Z_perm`** on that slice — a consistency the specs already
require (`directives/S11bB_SHARED_PHYSICS.md:203`), so unification is faithful, not a merge of rivals.

---

## B · Naming & provenance — the two resolved snags

**G3 · The two `v₀` are DIFFERENT physical quantities ⇒ DISTINCT standard names.**
S11's `v₀` is the background **brane's in-plane/tangential** velocity, pinned to zero
(`directives/S11_SHARED_PHYSICS.md`, supplied premise). S11b's `v₀` is the **bulk's normal drain across
a face** — it couples to the normal wavenumber `q_n` (`steps/S11bA_interface_response.md:99`). ⇒ S11b's
drain gets its **own** standard name (e.g. `v_bulk_normal_0`); ⛔ **never reuse the key `v_0`**, and ⛔
this is **not a premise override** (S11b does not set the brane moving).
⚠ **Why the name is load-bearing, not cosmetic:** the imported bare `v_0` key holds `Symbol('v_0')`
(`scripts/S11_exports.py:11730`), so an `F9` **object** comparison against a reused `v_0` would prove
them EQUAL and silently merge two quantities — the exact same-name-different-quantity failure `D5`
exists to prevent. The drain's ledger role is a **free parameter tagging a recorded neglected convective
correction** (a scope limit), ⛔ not an active term in the derived operator (both A and B leave it
uncarried).

**G4 · `Λ_A⁰` ≠ `Λ_p⁰` — a redefinition, kept distinct with a documented reduction.**
B's `Λ_A⁰` is the coefficient of the **full affinity** `𝒜 = μ_s − δp/ρ_m`
(`directives/S11bB_SHARED_PHYSICS.md:163,171`); A's `Λ_p⁰` is the **raw-pressure** coefficient
(`directives/S11bA_SHARED_PHYSICS.md:180`). ⭐ Canonical = `Λ_A⁰`. `Λ_p⁰` is a **distinct, derived
pressure-slice coefficient** carrying the reduction B documents as an acceptance check at `μ_s=0`
(`directives/S11bB_SHARED_PHYSICS.md:193,199`) — point at that line; ⛔ do **not** alias the two, and ⛔
do not reprint the reduction's value here (it is the acceptance value, withheld per `G13`).

**G5 · Three independent relaxation times `τ_A, τ_V, τ_X`; A's single `τ` is a restricted slice.**
B declares them free and independent (`directives/S11bB_SHARED_PHYSICS.md:167,187`); A imposes **one shared
`τ`** as its own scope limit (`directives/S11bA_SHARED_PHYSICS.md:194`). ⇒ the step derives the general
**independent-time** response and emits the equal-time A-slice (`τ_A=τ_V=τ`, `μ_s=0`) as a reduction
check. ⛔ `τ` is not an alias for any one B time.

**G6 · `Λ_X⁰`, `τ_X` are SUPPLIED constitutive inputs, ⛔ not inherited from A and ⛔ not derived into
existence.** B's reciprocal traction is a **prescribed** mechanical response
(`directives/S11bB_SHARED_PHYSICS.md:447`); A had no such channel. ⇒ introduce them with **S11b
`SUPPLIED`/`PREMISE` provenance**; the step DERIVES their consequences (power form, passivity region,
whatever relation between `Λ_X` and `Λ_V` is forced, if any — B says *assert no outcome in advance* (`directives/S11bB_SHARED_PHYSICS.md:296-303`), so ⛔ this list does not state it), ⛔ never the channel's existence. If a deeper derivation
of the channel is later wanted, that is **new physics beyond the current specs**, not a bookkeeping
repair.

---

## C · Chain wiring — inherit the pattern, with three named deviations

**G7 · The chain mechanism is S11's, inherited verbatim.** Point at: `F1`–`F9`
(`S11_export_chain_decisions_v2.md`); the blind-Wolfram control
(`S9_export_chain_rebuild_directive.md:16-18`); the exports spine, flat keys, `BUILD_INPUT_DIGESTS`, the
three-valued `F9` collision, `D1` iterate-the-emitted-collection, `D3` in-run round-trip. The SymPy
engine imports the **S11 LEDGER**; the consumed physics imports include `c_s0` **and `μ_R`** (B lists `μ_R` as an established input, `directives/S11bB_SHARED_PHYSICS.md:76-85`; `scripts/S11_exports.py:2182` carries it), bound to the imported
object, ⛔ never re-declared. ⚠ Whether B's `ρ_br⁰` equals S11's `ρ_br` is a mapping the spec must settle; `ρ_m` **originates in S11b** (absent upstream). `BUILD_INPUT_DIGESTS` pins exactly `{the S11b audit, S11_exports.py,
S11b_SHARED_PHYSICS.md}`.

**G8 · Three deviations from "copy S11 verbatim," each because S11 is itself incomplete or S11b differs:**
- **(a)** Author the **cross-engine comparator to the frozen `T7` contract** (`S11_C17_C18_spec_repair
  _decisions_v2.md:53-60`): join by name, residual the paired payloads, **reject a native boolean as a
  residual operand**, three-valued undecided, plus a repoint ablation. ⚠ **Net-new for the light
  sector** — S11 shipped without a comparator.
- **(b)** Restore the **`D3` in-run exec-and-residual round-trip** of the finished export module (S10 had
  it; S11 dropped it for serialization symmetry).
- **(c)** Include the `_RELATIONALS` reviver block (S11b emits interface-law relations/inequalities;
  S9 omitted it).

**G9 · The denylist is CUT.** `S11bA_SHARED_PHYSICS.md`'s "WHAT NEITHER ENGINE MAY READ" section is a
prohibition-based blindness control — a denylist, which rule 12 forbids. The **only** blindness control
is the blind Wolfram engine (`G7`). ⛔ Do not carry any do-not-read list into the S11b spec.

**G10 · The engines PRINT computed objects and residuals, then guard — ⛔ no terminal `VERDICT`.**
Point at the three script clauses and the "no VERDICT — a physics finding exits 0" rule. ⇒ the old
`.wl`'s 19/23 unfalsifiable gates and the `x===x` tautology
(`steps/S11bA_interface_response.md:164-172`) are **cut**; a check is the inherited **independent-routes + one-sided-corruption** obligation (`.claude/skills/build/SKILL.md` corollary 3 / rule 14) plus the blind second engine, ⛔ never a longer gate checklist and ⛔ never a definition-vs-its-own-substitution residual (which passes vacuously). The `S11bA_wl_repair_directive.md` (which mandated the
tautology) is **superseded, not folded**.

---

## D · Physics scope

**G11 · Carry forward (recompute; physics unchanged, only re-packaged).** The uniform-background interface
law: `q² = ω²/c_s0² − k²`; the three-regime bulk response `Z` and the grazing singularity; the per-face
added mass (its sign is a **convention artifact** — both faces load positively,
`steps/S11bA_interface_response.md:191-193`); the permeable response `Z_perm`; the projection identity;
the parity rule (**even window AND symmetric interval**); the flux channels; both validity inequalities;
the breathing-mode stability quantity and its separating condition (⚠ the closed root is a `k=0`/impermeable **slice**, `steps/S11bB_interface_assembly.md:140`, ⛔ not the general boundary — report the sign, ⛔ do not assume it); the transverse coupling on a uniform background — compute it (B6) and **if it vanishes, say so** (`directives/S11bB_SHARED_PHYSICS.md:750`), ⛔ do not assert zero; the closed quadratic energy-invariant basis — construct and check it (`directives/S11bB_SHARED_PHYSICS.md:345`), ⛔ the count is a result, not given. ⭐ Name the objects; the spec supplies the equations.

**G12 · Required fixes, folded in (each a measured defect, not a preference):**
- **(a)** The discarded convective order must be **DERIVED** from the convective dispersion relation `(∂_t+v₀∂_w)²φ = c_s0²∇²φ`, ⛔ not presupposed (an obligation-to-compute like `G13`; the known value is diffed on our side, ⛔ not stated here). ⚠ The old tag text was
  wrong in both engines (`steps/S11bA_interface_response.md:93-114,174`). ⚠ It was a **spec** defect
  (the spec presupposed a power of `v₀/c_s0`); every S11b spec leg must ask *"does the spec
  presuppose the FORM of any answer?"*
- **(b)** The velocity-coupled dissipation headline is **CONDITIONAL**: admissible **only with a named
  reservoir and a stated power budget** (`steps/S11b_HANDOFF.md:37-39`). ⛔ Never written as
  "forbidden/DEAD" — the engines emitted a **region tagged not-used-to-remove-the-root**, not a
  prohibition.
- **(c)** The two owed card items (`Λ_A⁰` used undefined; the dropped `Λ_p⁰=0` qualifier on the
  finite-memory departure) are fixed when the card is re-pointed.
- **(d)** B's **uncarried** background-flow (convective) correction is recorded as a standing scope
  limit, ⛔ not silently dropped again.

**G13 · Withhold the acceptance value (rule 5).** The `μ_s=0` reduction relating `Λ_p⁰` to `Λ_A⁰`, and
the A-slice reproduction of `Z_perm`, are **acceptance criteria referencing an expected value**: the
spec states them as an **obligation to compute-and-emit**, ⛔ never prints the target, and the **diff
happens on our side**. A supplied object stated in the spec is flagged as unfalsifiable within the build.

**G14 · The uniform interface questions S11 deferred to S11b are THIS step's; ⛔ only genuinely non-uniform or nonlinear work defers.**
⭐ **In scope here** (uniform background; these are B5's, and S11 handed them to "the S11b interface and full-slab spectrum", `steps/S11_stray_longitudinal.md:148-150`): the assembled **longitudinal mode's fate** — propagate / decay / grow / fail-to-exist as a normal mode — its dissipation origin, and behaviour at the sound-cone / `KW_ZERO_LOCUS` (`c_L=c_s0`, `q_n=0`) **grazing threshold** (`directives/S11bB_SHARED_PHYSICS.md:15-19,720-731`). ⚠ `Z`'s grazing singularity is already carried in `G11`, so the coupled-mode grazing question ⛔ cannot leave with it.
⛔ **Deferred to C** (non-uniform only): the full **variable-coefficient** slab spectrum, actual leakage rates in the non-uniform problem, and whether light's confinement is **unconditional**.
⛔ **Deferred to the nonlinear light program** (⛔ NOT C): the DC/harmonic/sideband radiation audit and nonlinear intensity coupling.
⇒ the step record names each deferral's correct owner, ⛔ never implies a uniform result closed a non-uniform question.

---

## E · What the legs must decide

⭐ For each G-item: is the **object named correctly** and is the **provenance** (imported / derived /
supplied / reduction-slice) right? ⛔ Two specific traps to hunt: **(1)** any place a name is reused
across a semantic boundary such that `F9`'s object comparison would prove a **false equal** (the `v_0`
hazard, `G3`); **(2)** any place the list **restates** an inherited obligation in weaker words instead of
pointing at it. ⚠ And the standing question of `G12(a)`: does any decision here **presuppose the form of
an answer** the engines are supposed to compute?
