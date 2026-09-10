# RE-AUTHOR the S11c-d SHARED PHYSICS spec — v4 (rule-15 author change)

You (`gpt-5.6-sol`) are re-authoring an **orchestrator-owned physics spec** under `CLAUDE.md` rule 15: the previous
author (the orchestrator) folded two revisions (v2, v3) and each bred over-corrections, so authorship passes to you.
Your v4 will then be reviewed by a fresh Claude agent + Grok (⛔ not by you). Write the best physics spec you can; ⛔ do
not defer to the prior author's phrasings where a finding shows them wrong.

## The task
**Overwrite** `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md` with **v4** — a complete, self-contained
physics-authority spec (the document two blind CAS engines, one SymPy + one Wolfram-that-imports-nothing, read to
independently construct the same objects, plus their comparator). Keep the house format (§0–§8), tag prefix `S11CD_`,
and the **v3 structural frame** (both review legs confirmed it sound); fold **every** round-3 finding below.

## Read first
1. The current base: `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md` (v3 — structurally sound, folds
   v1/v2). Preserve what its legs confirmed; change what round-3 flagged.
2. The two round-3 leg reports you must fold in full (must-fix **and** nits):
   `research/pde_ledger_v3/directives/_legs/S11c_d_shared_physics_review_r3_codex.md` (10 must-fix) and
   `research/pde_ledger_v3/directives/_legs/S11c_d_shared_physics_review_r3_grok.md` (3 must-fix + 4 nits).
3. Ground truth to check the folds against: `directives/S11c_decisions.md` (N5/N6/N7/N10–N15),
   `directives/S11c_a_SHARED_PHYSICS.md` §2a (`:171-210`, the `w₁`/`m₁` ansatz + the branchwise density maps),
   `steps/S11c_c2_self_energy_fold.md` + `_measurements/S11c_c2_N6_reconcile_disposition.md` (the c2 status/DEBT), and
   the real `scripts/S11c_c2_exports.py` keys.

## The round-3 fold-set (fold ALL; the leg reports carry the full derivations)
1. **§2 — do not assert uniform-Born/DWBA coincidence** (it silently needs the off-diagonal baseline `K₀=0`, which is
   the *withdrawn F*). Decompose and require the engine to **compute** `K₀`, `K₋`, `K₊`; define asymptotic channels
   from the **full** asymptotic block operators, not the bare diagonal `L±`; permit the uniform-mode simplification
   only *conditionally* on computed vanishing baselines. [Codex F1]
2. **§3c — the induced-amplitude decomposition.** `ΔA = A − A₀` still contains the zero-jet **contrast** term
   `(η¹σ_W⁰)`; separate the zero-jet (`η`) and first-jet (`σ_W`) amplitude components. The physical conversion flux has
   **baseline `O(ε²)` and interference `O(ε²λ)`** terms; `B(ΔA,ΔA)/J_in` is an **induced-field quadratic form**, ⛔ not
   the physical total conversion fraction. Apply the `O(εη)`/`O(ε²η²)`/`O(η²)` `N12` labels **conditionally** on the
   computed baseline/interference disposition. [Codex F2]
3. **§1c/§3 — Fourier premises.** Require short-range (`L¹`) decay for **every profile derivative the retained kernel
   consumes** (Riemann–Lebesgue needs it); the full profiles do **not** decay (unequal asymptotes) so treat step
   profiles in coordinate space or after an **explicit asymptotic subtraction**; **define the exact reduced 1-D
   transform** and the **mapping from every c2 3-D Fourier carrier** (the `(2π)⁻³` convention, the tangential
   `δ²(Q_∥)`, the per-unit-edge-area amplitude) so two engines cannot differ by `2π`/dimensional factors. [Codex F3]
4. **§1d/§3d — kinematics.** `Q_nL_W → 0` is the **zero-transfer** limit only (⛔ not "sudden", ⛔ not "maximal
   conversion" — overshoot profiles need not peak there, and the flux weights/other channels move the maximum); large
   `|Q_nL_W|` suppresses *this derivative form factor* (⛔ not "WKB"). **WKB/adiabaticity is a separate local
   wavelength/gap condition** (`kL_W`), forbidden as an additional `σ_W→0` expansion, ⛔ not identified with large
   `Q_nL_W`. Keep `η`, `σ_W`, **and** kinematics (splitting the form-factor argument `Q_nL_W` from the wave-sharpness
   `kL_W`) as distinct live quantities. ⛔ Do not take `L_W→0` at fixed `η` (`σ_W→∞`). Replace "form-factor node
   `Q_nL_W=nπ`" with "a zero of the computed form factor, if one exists" (the representative `tanh` has no real
   nodes). [Codex F4 = Grok F2]
5. **§1c admissibility — the branchwise density maps, verbatim.** `RHO4_CONSTANT`: `ρ_4D,bg⁰` constant,
   `ρ_br,bg⁰ = ρ_4D,bg⁰·W_bg` (varies). `RHOBR_CONSTANT`: `ρ_br,bg⁰` constant, `ρ_4D,bg⁰ = rho_br/W_bg` (**varies** —
   the N4 advection probe needs this gradient). Carry both density asymptotes/gradients live; ⛔ do not imply
   `ρ_br,bg⁰` varies in the RHOBR branch. [Codex F5]
6. **§5a — honestly DOWNGRADE the d-level N6.** A coordinate rewrite of the *already-constructed imported closed
   kernel* is **not** an independent shape/coordinate derivation of its N3/N4 content, and **no exported provenance
   isolates tilt from advection** (flipping the `w₁′` factor moves several channels). ⇒ relabel §5a as a
   **scattering-coordinate covariance regression** that ⛔ does **not** discharge the kernel-level N6/N3/N4
   independence — that was **c2's** control and stays **carried as unclosed DEBT**; the one-sided atom mutations are
   **shape-sensitivity** probes, ⛔ not clean channel isolation. Name the actual chart (in-plane `x = X + u` at fixed
   anchoring). Keep: `∇W_bg→0`/`η→0` rejected as a corruption, `RHO4_CONSTANT` computed absence (⛔ no `A−A`), no
   anchoring corruption, no `Φ` on amplitudes/modes/flux, print-not-target. [Codex F6 + Grok nit3]
7. **§5b — the jet-zero regression.** `w₁′ = 0` everywhere ⇒ `Δw₁ = 0`, so it **cannot** retain unequal asymptotes.
   Run **two separate** uniform regressions — the left constant background everywhere, and the right constant
   background everywhere (equal ends in each) — plus the `η = σ_W = 0` reference; the coupling is a **computed**
   object, ⛔ not asserted to vanish (withdrawn F). [Codex F7]
8. **§3a/§3c/§3b — canonical S-matrix / flux / confinement.** Require the **complete left/right channel matrix** (⛔
   not "one end OR full matrix"). Derive an **explicit modal flux bilinear current for every asymptotic channel** (⛔
   not just group velocity — the multi-component, `ω`-dependent, possibly non-Hermitian operator needs the mode's
   bilinear current with left/right modes and `∂_ωL`). The **converted** flux uses the **converted (thickness)
   channel's own current** (⛔ not "outgoing transverse-channel flux"); define `C_{T→H} = J_{H,out}/J_{T,in}`. Emit the
   **transverse survival functional** (reflected + transmitted transverse blocks) for confinement (`N13` = survival of
   the transverse channel), ⛔ not a vague condition. [Codex F8]
9. **§3b/§2 — bound pole.** **State** the truncation limitation next to the §2 exemption: Evans/Jost zeros are the
   spectrum of the **first-shape-order** imported operator; omitted `O(η²,σ_W²)` terms are the **same order** as a weak
   bound state (`E ∼ λ²`), so existence/location are **truncated-model data**, ⛔ not a controlled parent-theory
   channel; give the threshold/separation domain in which a pole may be promoted to a physical claim. Make the
   canonical object the **pole set + normalized Riesz residue/projector** (`⟨l, ∂_ωL r⟩ = 1`) — an Evans/Jost
   determinant is non-canonical (defined up to a nonvanishing analytic factor; the nonlocal operator needs
   trace-class/Fredholm premises for a determinant). Spectral-overlap-only unless a **concrete capture protocol** is
   supplied; keep the non-Hermitian physical-sheet/normalizable/zero-width/all-channels-closed criteria; keep "not a
   Bloch band." [Codex F9 = Grok F3]
10. **§7 — comparator debt.** A residual on the final **projected** mixing amplitude cannot *surface* the c2 carrier
    (40)/source (76)/Φ (18) operand debt (projection + channel-sum may cancel them; the disposition says those are a
    **schema non-join**). ⇒ say the d amplitude **remains conditional on and PROPAGATES** the debt; claim direct
    surfacing only if d emits separately-defined common-basis projections of those operand families with a sound join.
    [Codex F10]
11. **§0/§1c/§5c — class membership.** The class is a **thickness** interface: `Δw₁ ≠ 0` ⇒ `W₋ ≠ W₊` (the N5 object).
    `m₁` stays **independent** (constant, bump, or interface — all in-class); ⛔ `Δm₁ ≠ 0` is **not** a class gate — it
    is a **per-profile discriminant of the modulus subchannel**. The edge-vs-bump discriminant (§5c) is `Δw₁` for the
    thickness (and analogously `Δm₁` for the modulus control), ⛔ not "either jump vanished." [Grok F1]

**Nits:** define `J` from the S11b energy current / the imported closed operator, ⛔ not c2's EMIT-only traction–slab
pairing tag [Grok nit2]; §5d — "emit the computed projection", ⛔ not a typed `∝ k·a` shape [Grok nit1]; cite the
`w₁`/`m₁` **definition** at `S11c_a_SHARED_PHYSICS.md:171` (not `:190`, which is the `η`/`σ_W` independence) [Codex
nit1]; fix the "both residuals" antecedent [Codex nit2].

## Governing constraints — ⛔ do NOT violate while folding
- **M2 (name the object; withhold only the acceptance criterion).** The spec says what to **compute**; ⛔ it states no
  component value, sign, order, parity, or grade **beyond** the supplied `(ε,η,σ_W)`/`λ` power-counting contract.
  Withhold exactly one thing: an acceptance criterion referencing an **expected value** (the falsification numeric
  bound / the `O(1)` reductio — kept orchestrator-side). ⛔ **No leaked target** anywhere: no "nonzero", no "must
  vanish", no "must move" — controls emit **baseline + operand + residual**, disposition adjudicated our side.
- **The c2-honesty block (§1b) is load-bearing — ⛔ do NOT weaken or "close" any of it.** Keep: the two closed
  operator/kernel VALUES per-engine SOUND; the cross-engine content = the N6 covariance thread only (Reading B, the
  matched zeros are `(0)−(0)`, ⛔ not operand agreement); `R_N6 = 18/288` is the **per-engine SymPy raw /
  schema-unmatched** result; the operand DEBT (carrier 40 / source 76 / Φ 18) UNADJUDICATED and **material to this
  consumer**; the leftover SHAPE ⛔ not "just thickness"; **F and G WITHDRAWN** (paused, not discharged, ⛔ not
  established); the two S11c-b sign conventions, six §3d re-adjudications, three N6 premise caveats — all named.
  Preserve **both** `R_N6=18/288` and `R_cov` no-nonzero.
- **No new upstream physics (`N15` at the right layer).** d consumes c2's operator/kernel **verbatim** and emits
  **profile moments / form factors derived from the imported kernel** — ⛔ never a new local constitutive constant; a
  missing invariant is recorded as an **upstream N15 debt**.
- **Consume only the REAL c2 export rows:** `s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`,
  `s11cc2Fieldtheta`, `s11cc2FieldeW`, `s11cc2Fieldu{1,2,3}`, `s11cc2Coefficientw1Profile`,
  `s11cc2Coefficientm1Profile`, `s11cc2FourierW1ProfileHatTransfer`, `s11cc2FourierW1ProfileJetHat*` (+ `*Dimension`),
  and the reachable constants. ⛔ There are **no** term-origin, parity, self-energy-increment, or §3d export rows — the
  increment is EMIT-only parent provenance, ⛔ not an import. The exact `IMPORT_KEYS` root set is fixed later at the
  build directive, ⛔ not enumerated-then-frozen in the spec.
- **Keep the settled frame — ⛔ do not re-litigate it, only fold the findings:** the localized-interface class + Born
  in contrast `η` with `σ_W`/kinematics live; the strong-edge as a *named downstream obligation* (a nonzero Born
  coefficient gives no lower bound — keep the `sin²(ηG)` counterexample; `C_strong(1)`, ⛔ not "`F(1)`"); the
  distorted-wave organization; the two-channel `N13`; `N11a` inert; the chain/`T7` comparator/blind-WL/`N14`
  reservations; supplied-vs-computed split. The house `[[wiki-link]]` cross-references may be kept where v3 used them.
- **This is orchestrator-context material** carrying parent physics values (`R_N6=18/288`, `R_cov=0`) — those are
  **parent facts already in the spec**, keep them; they are ⛔ not answers to withhold (the withheld thing is only the
  future falsification bound). The spec is a **physics-authority document**, not a blind build.

## Output
Overwrite `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md` with the complete v4. At the end of your run,
print a short **fold report**: per round-3 finding (1–11 + nits), one line naming the section you changed and how.
⛔ Do not print the whole spec back. ⛔ Make no other file changes.
