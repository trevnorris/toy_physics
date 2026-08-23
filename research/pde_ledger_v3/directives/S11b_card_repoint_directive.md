# S11b card re-point directive (Codex writes the `.tex`; its own 2 legs after)

## Authority and boundary
Re-point `research/pde_ledger_v3/paper/steps/S11b_interface_coupling_law.tex` IN PLACE (committed baseline
`8ddccb74`-era; the card itself is the Aug-4 A/B version). `CLAUDE.md` binds. ⭐ **The source of truth is the
committed step record `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md`** — the card must state
nothing the step record does not, and must reflect its corrected physics. Read that step record and
`directives/S11b_SHARED_PHYSICS.md` first. This is a documentation re-point of already-computed, cross-checked
results (⛔ no engine runs, ⛔ no rule-5 withholding — the results are settled and may be stated).

⚠ **Preserve what already works:** the `\claimstatus`, the "region not prohibition" framing of the
admissibility section, the two refuted preregistered predictions (two surviving scalar dof; added inertia
`ρ_m/(2α)` per face, not `ρ_m/α`), the `\paragraph{Verification scope.}` handled as a **plain paragraph**
(⛔ NOT `\stagefield{Verification}{...}` — that field is suppressed in the default build, `paper/macros.tex:19`,
and the scope/limits must stay visible), and the `\paragraph{Known limits --- recorded, not fixed.}` list
structure. ⛔ Do not delete a load-bearing limit.

## The changes (each cites the step record / decision list)

### 1. Re-point the source records
The "Inputs and regime" paragraph currently names `S11bA_interface_response.md` + `S11bB_interface_assembly.md`
as "the only source records." ⭐ Make the **unified step record**
`research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` the source of record (A/B are archival
pre-unification history — name them as that, ⛔ not as the current source). S11b is now ONE unified
export-chain step (A + B subsumed; **C — the non-uniform transverse coupling — is a separate later step**).

### 2. G12(c) — define `Λ_A⁰`, and restore the dropped `Λ_p⁰=0` qualifier
- **Define `Λ_A⁰`** (used undefined in the admissibility region). It is the **DC (`ω→0`) affinity-channel
  flux coefficient**: the face flux is driven by the affinity `𝒜 = μ_s − δp/ρ_m` through `Λ_A(ω) 𝒜` with
  `Λ_A(ω) = Λ_A⁰/(1 − iωτ_A)` (parallel to the card's existing `Λ_p`, `Λ_V` kernels). Relate the affinity and
  raw-pressure channels by the **static slice map** `Λ_p⁰ = −Λ_A⁰/ρ_m` (the `μ_s=0` reduction; both engines
  now emit it — `ZPERM_SLICE_MAP`, step record §"cross-engine comparison"). Define `Λ_V⁰`, `Λ_X⁰` similarly as
  the DC velocity- and reciprocal-traction-channel coefficients so the admissibility region reads with all
  symbols defined.
- **Restore the `Λ_p⁰ = 0` qualifier** on the finite-memory departure. The "Departure --- finite-memory
  velocity channel" branch is the **pure velocity branch** (`Λ_p⁰ = 0`, i.e. the raw-pressure DC channel is
  off; the coupling is velocity-driven). State that qualifier where the branch is introduced (it was dropped
  in the current card).
- **Narrow the `τ_A` qualifier.** The card says "while `τ_A` is unrestricted." Correct it to the step record:
  `τ_A` is unrestricted **by the passivity conditions** — no relation beyond the supplied `τ_A ≥ 0`; ⛔ not
  that negative relaxation times are admitted.

### 3. G12(d) — record the background-flow correction as a STANDING scope limit
Add to "Known limits --- recorded, not fixed." the **uncarried background-flow (convective) correction**:
`O(v₀ |q_n| / ω)`. The engines emit the relative term `−2 q v₀ / ω` (`RELATIVE_TERM`), so first order fails
where `|q v₀ / ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` —
large `k c_s0/|ω|` is **necessary, not sufficient**. ⛔ Do not silently drop it again; it was recorded by A,
lost in B's first draft, and restored only after a review caught it. (`v₀` is the same steady drain named as
the departure's candidate reservoir.)

### 4. Transverse null — the stability is CONDITIONAL (⚠ write it in the card's retained representative)
The card states `ρ_br⁰ ω² = μ_R k²`, `Im ω = 0` **unconditionally**. Correct it to the step record's result,
but ⛔ **keep the card internally consistent with its `(∇·u)²` energy representative** (change 6): the card's
energy list carries `(∇·u)²`, no `st²`/`μ_S`, so the transverse working stiffness is `μ_R`, ⛔ **NOT
`μ_R + μ_S/2`** — displaying `μ_R + μ_S/2` beside a `μ_S`-free energy would make the EOM not follow from the
displayed energy and would use `μ_S` undefined (the very defect G12(c) fixes for `Λ_A⁰`; carrying both an
independent `st²` and an independent `(∇·u)²` would restore X-1's eleventh invariant).
- The transverse-to-thickness coupling is **identically zero** and the mode is **non-dissipative**
  (`TRANSVERSE_DISSIPATION ≡ 0`) — both unconditional.
- Write the oscillator `ρ_br⁰ ω² = μ_⊥ k²` with the transverse stiffness `μ_⊥ = μ_R` **in the card's retained
  `(∇·u)²` representative**. `Im ω = 0` (a stable real-frequency mode) holds **only where `μ_⊥ = μ_R ≥ 0`**;
  ⛔ §5 assumes no positivity (§0 — decay and growth are both admissible), so `μ_⊥ < 0` is admissible and gives
  a growing transverse root `ω = ±i k √(|μ_⊥|/ρ_br⁰)`.
- ⭐ **Cross-engine note (parenthetical, ⛔ not the working formula):** the same physical stiffness is `μ_R` in
  the blind-WL / card `(∇·u)²` representative and `μ_R + μ_S/2` in the SymPy `st²` representative (the
  energy-basis representative split); ⛔ do not introduce `μ_S` into the card's working equations.
- Keep the "dimension of the vanishing coupling is undetermined" note and the "does NOT settle unconditional
  confinement (that is C's non-uniform question)" caveat.

### 5. Verification scope — re-point to the frozen T7 cross-engine comparison
The current "Verification scope." paragraph describes single-engine `VERDICT: PASS` strings. Replace with the
current architecture: the SymPy engine (imports the S11 LEDGER) and the **blind** Wolfram engine (imports
nothing, re-derives — the only cross-engine control) were compared by the **frozen T7 comparator**
(`scripts/S11b_cross_engine_comparator.py`, run `scripts/out/S11b_cross_engine_comparison.out`). ⭐ **No
compared object is a physics contradiction**; the disagreements are emission-format (SymPy `Tuple` vs Wolfram
`Association`), coefficient-basis (the representative split), naming, a `W₀` normalization convention, a
denominator-clearing form/domain artifact (the `DEGENERATE_LOCI` loci coincide on the regular domain, differing
only where a cleared denominator vanishes — the grazing locus `q=0`), and a few **coverage gaps** (one engine
solves what the other leaves open — ⛔ not a conflict). ⚠ Keep the
load-bearing limit that **the independent physics verification comes from reviewer derivations, ⛔ not from
the engines agreeing** — a supplied object lands identically in both.

### 6. Energy basis — note the total-divergence quotient (X-1)
The card correctly reports **ten** independent quadratic invariants. Add the brief provenance: the count is
**ten under §5's equivalence modulo total in-plane divergences**; the basis representative is **not unique**
(the two engines retained different but equivalent representatives — one the traceless-symmetric `st²`, one
`(∇·u)²` — spanning the same 10-dimensional quotient). ⚠ One engine initially over-counted to eleven by
omitting the divergence quotient and was corrected (recorded in the step record; keep the card to one clause).
⛔ Do not restate the five new terms wrongly — keep the card's existing list, which uses the `(∇·u)²`
representative.

## Fidelity and scope
⭐ Every physics statement must match the step record; ⛔ add no new result. ⚠ Keep the card self-contained and
readable; the ledger card records the step's result, ⛔ it does not re-derive. Where the card must name a
computed object, name it as the step record does. ⛔ Do not touch other cards or `paper/macros.tex`.

## Report (§13) — under 15 lines
The edits (paragraph → what changed), a confirmation that every physics statement traces to a step-record
line, that the suppressed-`Verification` handling is preserved (plain paragraph), and that the "region not
prohibition" framing and the two refuted predictions are intact.
