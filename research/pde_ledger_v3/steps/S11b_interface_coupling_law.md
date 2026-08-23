# S11b — the linear brane–bulk interface coupling law (step record)

Slug `S11b_interface_coupling_law`. This is the unified step record for S11b, subsuming the historical A
(bulk face response) + B (homogeneous assembly) execution stages into **one export-chain step**; **C — the
non-uniform transverse coupling — is a separate later step** (deferred, scope ratified in the decision list
`directives/S11b_unified_decisions.md` G14, `ddd0ae4c`). The physics authority is the unified spec
`directives/S11b_SHARED_PHYSICS.md` (`1a2395a3`).

⭐ This record is the **interpretation** layer; every result below names the computed object and the commit
that produced it. The engines PRINT objects; they state no conclusions (the conclusions are here).

## What the step computes

The linear response of the brane–bulk interface: given a slab of finite thickness `W` in the `w` direction
with two faces meeting the bulk, derive how a normal face motion couples to the bulk, what drives material
across the interface, and the closed spectrum of the assembled system on a **uniform** background.

## The two engines and their agreement

- **SymPy engine** `scripts/S11b_interface_coupling_law_sympy_audit.py` — the packaging engine: imports the
  S11 LEDGER (1663 rows), binds `c_s0`/`μ_R`/`ρ_br⁰=rho_br` to the imported objects, writes
  `S11b_exports.py` (1958 rows). Committed `864d6f41`; X-1 repair `53fcd98d`.
- **Wolfram engine** `mathematica/S11b_interface_coupling_law_mathematica_audit.wl` — **blind**: imports
  nothing, re-derives every object from §1–§6 of the spec. Committed `ec89f9df`; repair `bd598ae7`. Its
  blindness is the only cross-engine control (it cannot transcribe the `.py` because it never reads it).
- **T7 comparator** `scripts/S11b_cross_engine_comparator.py` (reused from S10) — joins the two transcripts
  by emitted-object name, residuals paired payloads, rejects a native boolean as a residual operand, is
  three-valued, and was frozen before it saw either output. Re-run against both repaired engines: `fba6a34c`
  (`scripts/out/S11b_cross_engine_comparison.out`).

⭐ **No compared object is a physics contradiction.** After the WL repair (F-WL-1) and the SymPy X-1
repair, the comparator's headline is AGREE 21 / DISAGREE 108 / UNCOMPARED 25 / UNPAIRED 71; the 108
disagreements are **format, coefficient-basis, naming, convention, and a few coverage gaps — not conflicting
values** (adjudicated below). ⚠
**Finding — the engines are not emission-parallel:** a SymPy engine emits Python `Tuple`s and a Wolfram
engine emits `Association`s, so the comparator marks 102 objects `STRUCTURE_DISAGREE(tuple vs Association)`
even where the physics agrees (the parsed values coincide, or differ only by the invertible coefficient-basis
remap of the representative split — see the adjudication below). This is a property of the two CAS output
shapes and the basis choice, not a physics divergence; it is why `FINAL_OPERATIONAL_STATUS` reads `FAIL`
while the physics agrees.

## The result

⭐⭐ **The interface response is a REGION, ⛔ not a prohibition** (this was an early overclaim, corrected
2026-08-05: the engines emit a classification tagged *"not used to remove the root"*, under a spec stating
*"admissibility is a classification, not an acceptance gate"* — the verdict must not be manufactured in the
writeup).

**Passivity region** (the power form ≥ 0 — the interface is not an energy source):
```
Λ_A⁰ ≥ 0   and   [ Λ_V⁰ = Λ_X⁰ = 0   or   (Λ_X⁰ = −Λ_V⁰ , τ_V = τ_X = 0) ]
```
A separate, **conditional** Onsager–Casimir test gives `Λ_X(ω) = −Λ_V(ω)` by two independent routes; `τ_A`
is unrestricted by the passivity conditions (no relation beyond the supplied `τ_A ≥ 0`, ⛔ not that negative
relaxation times are admitted).

⚠⚠ **The passive region does NOT bound this model**, for two reasons named in the walk:
1. Onsager–Casimir requires microscopic time-reversibility, which the model **postulates nowhere** — the
   medium is a substructure; thermodynamic laws emergent from particles made *of* it are not inherited by it.
2. Passivity/reciprocity describe fluctuations about an **equilibrium**; the reference state carries `v₀ ≠ 0`
   (a driven steady state). ⇒ **Standing rule:** a non-passive coupling is admissible only with a NAMED
   reservoir and a STATED power budget; the real gate is **observational, not thermodynamic**. The model
   supplies a candidate reservoir — the background drain `v₀`.

**The chemical-potential drive.** Material crosses the interface driven not by bulk pressure alone but by an
affinity `𝒜 = μ_s − δp/ρ_m`. Under a pressure-only closure a compressed brane next to a calm bulk had
nothing driving conversion, which is wrong — a compressed brane *wants* to shed material. The `μ_s = 0`
reduction turns the flux pressure-driven but ⛔ does not remove the velocity channel `Λ_V V`.

**The velocity channel and its thermodynamic fate.** The velocity-driven conversion is an energy **source**
unless it is instantaneous (`τ → 0`) — inside the passive region the second law forbids it, so the model
refuses an unphysical process on its own terms. Outside the region it survives only against a named reservoir.

**Transverse mode — decoupled on a uniform background.** The transverse-to-thickness coupling
`∂²U/∂u_T ∂e_W` is **identically zero** and the transverse mode is **non-dissipative** (`TRANSVERSE_DISSIPATION`
≡ 0, unconditionally) — a decoupled oscillator `ρ_br⁰ ω² = μ_⊥ k²` (`S11B_TRANSVERSE_DISPERSION`,
`S11B_TRANSVERSE_COUPLING`). ⚠ It is a **stable real-frequency** mode (`Im ω = 0`) **only where the transverse
stiffness `μ_⊥ ≥ 0`**: §5's moduli are free and ⛔ no positivity is assumed (§0), so `μ_⊥ = μ_R + μ_S/2 < 0`
is admissible and gives a growing transverse root `ω = ±i k √(|μ_⊥|/ρ_br⁰)` — decay and growth are both
admissible outcomes here, ⛔ not excluded. ⛔ This is the **uniform** limit only and does NOT settle
unconditional confinement — that is the non-uniform question deferred to C. ⚠ The two
engines express the transverse stiffness `μ_⊥` in **different coefficient bases** (a consequence of the
energy-basis representative split, below): the blind WL engine emits `μ_⊥ = μ_R`, the SymPy engine
`μ_⊥ = μ_R + μ_S/2`, with identical roots under the invertible coefficient identification WL-`μ_R` ≡
SymPy-`(μ_R + μ_S/2)`. ⛔ They are NOT equal under the comparator's naive `muR↔mu_R` name-transliteration,
which is why `TRANSVERSE_DISPERSION` shows a cross-engine difference (an extra `k²μ_S/2`, masked under its
tuple-vs-Association STRUCTURE flag) — the same physical stiffness, not a physics disagreement.

**Breathing-mode stability.** On the `k = 0` slice with impermeable faces (`Λ_A⁰ = Λ_V⁰ = 0`) **and no
reciprocal traction (`Λ_X⁰ = 0`)** — all three cuts, since "impermeable" alone does not set `Λ_X⁰` — the
bulk load is kept and `K₀ = B_ρ⁽³⁾ − 2C W₀ + k_W W₀² > 0`; growth iff `C > (B_ρ⁽³⁾ + k_W W₀²)/(2W₀)`. ⚠ The growing root is **not** an energy-conservation violation — the
stored energy has no minimum in that direction and the accounting closes exactly. ⚠ This is a slice result,
⛔ not the general breathing boundary.

## The energy basis and the X-1 correction

⭐ The spec (§5) carried **five** stored-energy terms but instructed each engine to CONSTRUCT the symmetry-
allowed basis and check its closedness (`S11B_ENERGY_BASIS_COUNT`). The engines found more allowed
invariants than the spec named — including `(∇·u)²`, the coupling that fell out of the specification four
times while it lived in prose.

⛔ **X-1 — a one-engine over-count, caught by the sibling.** The SymPy engine first emitted **11** basis
invariants; the blind WL engine emitted **10**. §5's symmetry group states *"equivalence modulo total
divergences — two densities differing by a total in-plane divergence are the same term; do not count both."*
The SymPy independence test judged only pointwise polynomial independence over the field components and
omitted that quotient; two of its `∇u`-sector invariants differ by a total in-plane divergence
(`st² ≡ ½curl² + ⅔(∇·u)²`, since `(∇·u)² − tr((∇u)²)` is a total divergence). The spec-correct count is
**10**; the SymPy engine over-counted. Verified independently three ways (the two build legs' own
Euler–Lagrange enumerations, and the orchestrator's own `x1_independent_basis_count.py`).

The X-1 repair (`53fcd98d`) made the SymPy independence judgment honor §5's quotient (Euler–Lagrange-signature
rank) and eliminated the redundant invariant by **REWRITING** it into the retained curl²/strain invariants —
folding its coefficient in — ⛔ not by deleting the term (which would have deleted the stiffness). The two
routes are checked by an emitted `ENERGY_REEXPRESSION_RESIDUAL` whose Euler–Lagrange derivative is printed
(identically zero) and demonstrably goes nonzero under a wrong fold.

⚠ **Representative split (a computed non-uniqueness, not a disagreement).** The 10-dimensional quotient has
no canonical basis: the SymPy engine kept `st²` and dropped `(∇·u)²`; the blind WL engine kept `(∇·u)²` and
dropped `st²`. Both span the identical complement, and `ENERGY_BASIS_COUNT` AGREES (10 = 10). This split is
why the comparator flags `ENERGY_BASIS_OMISSIONS` and `ENERGY_BASIS_INDEPENDENT_TERMS` as differing — a
basis-choice artifact, ⛔ not a physics difference.

## Cross-engine comparison — the adjudication (rule 13)

⭐ **No compared object is a physics CONTRADICTION** — the two engines never emit conflicting values for the
same quantity. The frozen re-run `fba6a34c` records AGREE 21 / DISAGREE 108; the 108 fall into four benign
classes plus a small set of coverage gaps:
- **102 STRUCTURE**, of two kinds. (i) *Output format*: SymPy `Tuple` vs Wolfram `Association`/`Rule` shapes,
  values coinciding once the shape is stripped (e.g. `ZPERM_SLICE_MAP` = `−Λ_A⁰/ρ_m` in both). (ii)
  *Coefficient-basis remap from the energy-basis representative split*: because SymPy retains `st²(μ_S)` where
  WL retains `(∇·u)²`, any object carrying a ∇u-sector energy coefficient is in a different but
  invertibly-related basis (WL `μ_R ≡` SymPy `μ_R + μ_S/2`; WL's `(∇·u)²`-coefficient `≡` SymPy `B_div + 2μ_S/3`),
  so the comparator's `muR↔mu_R` transliteration leaves a residual (e.g. `TRANSVERSE_DISPERSION`'s extra
  `k²μ_S/2`) that the coefficient map removes.
- **3 KEY** — Association-key-name mismatches (`CONTROLS_ON_TRANSVERSE`, `DISSIPATION_ORIGIN`,
  `LIMITS_AND_PATH`), ⛔ not a value split (the readable one, `CONTROLS_ON_TRANSVERSE`, is coupling 0 on both).
- **1 DIM** — `DIM_THICKNESS_RESPONSE` differs by one length (`W₀`) in the **drive normalization** (WL's route
  is `thicknessDisplacement − generalizedPressure`, i.e. δW per pressure; PY uses an independent
  force-per-volume route) — a convention, ⛔ not a physics disagreement; all other dimensions AGREE.
- **2 CONTENT** — both representation: `ENERGY_BASIS_OMISSIONS` (the representative split above), and
  `DEGENERATE_LOCI_EQUATIONS` — the loci are the **same equation on the regular domain** (PY = `(−i)·`WL's
  polynomial after clearing PY's denominator `q_out ρ_m (ω τ_A + i)`, with `q_out ≡ q`); ⚠ they differ only
  where that denominator vanishes — the grazing locus `q = 0` and the kernel pole `ω τ_A + i = 0`, which PY's
  cleared form excludes and WL's polynomial retains. A denominator-clearing form/domain artifact.
- ⚠ **Coverage gaps (flagged STRUCTURE, but not mere format).** A few objects are computed by ONE engine and
  left open by the other — ⛔ a difference in what was *solved*, not a conflict: `DEGENERATE_LOCI_SOLUTION`
  (PY solves `Λ_A⁰ = (i ω q ρ_m τ_A − q ρ_m)/ω`, WL emits an empty `SOLUTION`) and `ONSAGER_DETERMINABLE`
  (PY tags `UNDECIDED`, WL solves both routes with residual 0). Neither overturns a physics conclusion — the
  reciprocity `Λ_X = −Λ_V` is independently carried by WL's solve *and* PY's generic branch — but the record
  ⛔ does not claim every disagreement is pure format.
- **F-WL-1 cleared** — the withheld B2c/G13 acceptance map now agrees cross-engine: both emit the static
  `Λ_p⁰ = −Λ_A⁰/ρ_m` (the WL repair fixed a wrong-sign, wrong-level extraction; `ZPERM_SLICE_MAP`).

## Static-or-instantaneous freezes, standing limits, and owed items

- **Frozen wall width.** `ρ_br⁰ ≡ ρ_4D⁰ W₀` **is** the imported `rho_br` (S11's wall width is frozen). This
  is a recorded **freeze**, ⛔ not a fix — named here per the static-or-instantaneous check.
- ⛔⛔ **The background-flow correction `O(v₀ |q_n| / ω)` is uncarried and unbounded.** The engines emit the
  relative term `−2 q v₀ / ω` (`RELATIVE_TERM`), so first order fails where `|q v₀ / ω| ≳ 1`; in the
  `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is
  **necessary, not sufficient** (it also needs `|v₀|/c_s0` non-negligible). `v₀` is the user's dark-energy
  leak; it was recorded
  by A, lost in B's first draft, and restored only after a review caught it. It is the top known limit and is
  carried forward to C (and to the nonlinear program for the DC/harmonic/sideband radiation audit).
- **Owed card items** (deferred at the original close; ⛔ do not silently drop): G12(c) `Λ_A⁰` used undefined
  and a dropped `Λ_p⁰ = 0` qualifier — fix at the card re-point; G12(d) B's uncarried background-flow
  convective correction `O(v₀|q_n|/ω)` recorded as a standing scope limit.
- **C is deferred** — the non-uniform variable-coefficient transverse spectrum (is light's confinement
  unconditional?). The uniform coupling is identically zero, so a uniform-limit control for C is
  known-vacuous; C's coefficient is bounded by **bench-top optics** (a slit edge converting an O(1) fraction
  of a photon into the thickness channel would break diffraction gratings) — falsifiable with no cosmology
  and no gravity sector.

## Operational note

⚠ Long Mathematica jobs were **spuriously killed** three times this session (not OOM, not the watchdog, not
cron). The blind WL engine's repair review was carried by the fresh-Agent ablation leg + a Grok source audit
+ orchestrator diff-verification; **Grok's full WL ablation leg was blocked by those kills** (the Agent path
completes under Mathematica; the grok CLI did not). Prefer the Agent path for any future Mathematica ablation
leg, and recover a spurious kill mechanically (the edits complete before the kill). The SymPy X-1 build and
the comparator re-run are pure Python and immune.

## Provenance

Decision list `ddd0ae4c` (G1–G14, 2 legs) · spec `1a2395a3` (2 legs) · per-engine build directives `9bd2f184`
(2 legs each) · SymPy engine `864d6f41` (2 legs, 2 emission repairs) · blind WL engine `ec89f9df` (2 legs) +
repair `bd598ae7` · T7 comparator build/freeze `17fe32c8` · SymPy X-1 directive `b4c02381` (2 legs) + build
`53fcd98d` (2 legs) · comparator re-run `fba6a34c`. Run checklist: `steps/S11b_RUN_CHECKLIST.md`.
