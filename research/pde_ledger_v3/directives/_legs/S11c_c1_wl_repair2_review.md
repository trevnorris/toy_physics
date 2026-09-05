# Independent physics re-review — S11c-c1 Wolfram engine, REPAIR-2 (Codex-written)

## Artifact
`research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` (working tree; repaired against
baseline commit `13f0bd2c`). It is a **Mathematica script**; its only product is a flushed stdout tag stream (51
tags). You are one of two serialized legs (the other ablates the same `.wl`; run ONE kernel at a time — the
licence has two seats).

## ⛔ METHOD — form your own target view BEFORE opening the artifact (a re-review gap being closed)
Prior re-review legs opened the artifact first, which weakened them. So, in order:
1. **First**, from the spec below, write your OWN short derivation (a `.wl` or `.py`) of the two objects under
   test — (a) what the Fredholm/noninvertibility operator's `Z` must be, and (b) what the per-parity dissipation
   Hermitian form must be — and save it + its literal stdout to named absolute paths. ⛔ A prose derivation is
   worth nothing; without a script + literal stdout your derivation claims are discarded.
2. **Then** open the artifact and ablate it (below). Quote the spec text and the artifact line for every finding.

## What was repaired (the two emit-only objects to re-review) — SUPPLIED
Repair-2 fixed exactly two objects; everything else must be **byte-identical** to `13f0bd2c`.
- **R1 — `NONINVERTIBILITY_CONDITION` / `fredholmFunctionSpaceOperator`** (`.wl` ~`:574-600`, consumer ~`:1520-1555`):
  the operator whose invertibility is in question must be `[I + (Λ_A/ρ_m²) Z]` with `Z` the **two-momentum DtN
  operator carrying BOTH legs** `q_out(momentumOutput)` and `q_out(momentumInput)` (the same construction as
  `DTN_OPERATOR`/`operatorCompositionRecordFromDerivation`) — the prior code froze both `N_0` factors to
  `momentumOutput` (a single-leg WKB freeze, spec §3a forbids). The only legitimate diagonal symbol is the
  already-emitted **flat** one (`DEGENERATE_LOCI_*`); a `k=k′` slice of the full operator via a fresh dummy is NOT
  a legitimate second object.
- **R2 — `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`** (emission ~`:1714+`): must emit the **parity-combination** memory
  Hermitian form per parity (`PARITY_DELTA_W`, `PARITY_ZETA_C`), formed by the SAME change of basis
  `PERMEABLE_PORT_HERMITIAN`/`portParityCombination` uses — the prior code emitted the per-face form under both
  parity keys (a dead parity axis). `UNIFORM_LIMIT_*` was permitted to stay byte-identical (a lighter NIT deferred
  to the step record).

## Sources of truth
`directives/S11c_c1_SHARED_PHYSICS.md`: §1 `:76-83` (`δW≡ζ₊−ζ₋`, `ζ_c≡(ζ₊+ζ₋)/2`, `ζ_s=ζ_c+sδW/2`),
`S11c_a_SHARED_PHYSICS.md:54-59` (`V_s=s∂_tζ_s`, the outward orientation), §3a `:247-278` (the two-momentum
operator + face-parity structure), §3b `:299-320` (the Fredholm condition + the three-object dissipation audit "per
regime pair and parity combination"). Siblings as needed. `S9…:16-18` is the blindness control. `CLAUDE.md` rules
2, 3, 5, 6, 11, 17. ⛔ Do NOT read the SymPy engine or the step records.

## The ablations (MANDATORY — a FORM change, not a coefficient rescale; report literal stdout)
Copy the `.wl` to `/tmp` and ablate the COPY (⛔ never the working tree).
1. **R1 both-legs.** Confirm the emitted `NONINVERTIBILITY_CONDITION` operator's `Z`-part contains BOTH
   `qOut[omega, momentumOutput]` and `qOut[omega, momentumInput]` (the probes `FREDHOLM_Z_HAS_Q_INPUT`/`_OUTPUT`).
   Then **constructor-level re-freeze** the input leg (set the right-leg momentum to `momentumOutput`, as
   `DTN_OPERATOR` re-freezes) and confirm `FREDHOLM_Z_HAS_Q_INPUT` **MOVES** (both-present → input-absent). ⛔ Check
   this is a genuine constructor change, not a syntactic `Count` a payload-wide sprinkle of `qOut[…momentumInput]`
   could fake. Confirm the flat diagonal symbol / `DEGENERATE_LOCI_*` are byte-identical to baseline.
2. **R2 parity combination.** Confirm the base `PARITY_DELTA_W` and `PARITY_ZETA_C` dissipation payloads are
   **distinct computed objects** in the curved case (⛔ not byte-identical, as they were pre-repair). Then
   **one-sided sign-flip the `+`-face memory Hermitian form at its construction** and confirm the two
   parity-combination payloads **MOVE DIFFERENTLY** (not identically) — a genuine per-face→parity dependence, not a
   dead axis.
3. **Protected core byte-identical.** Diff the emitted values (or the construction) against `13f0bd2c` for
   `PERMEABLE_PORT_HERMITIAN`, `DTN_KERNEL`, `DTN_FLAT_SYMBOL`, `DTN_OPERATOR`, `FACE_RESPONSE`/`_COEFFS`, the
   `ENERGY_*` audit, T-a..T-i, and confirm they are unchanged.
4. **⭐ Independently adjudicate `PERMEABLE_PORT_HERMITIAN` (the object repair-2 did NOT touch, on the ground that
   it is correct).** Derive, in your own script, whether its diagonal blocks are the CONGRUENCE
   `Aᵀ·diag(P₊,P₋)·A` (`A=faceToParityMatrix`) of the power-conjugate form given `V_s=s∂_tζ_s` — i.e. both
   diagonal blocks EVEN combinations, the odd combination the coupling — or whether a `DELTA_W`↔`(ζ₊−ζ₋)` reading
   is right. If your derivation disagrees with "leave it byte-identical," that is a MUST-FIX finding: say so with
   the computation.

## Rule-5 / independence checks
For R1's re-freeze and R2's corruption: confirm each control genuinely MOVES and is NOT a structural zero killed by
a trivial substitution; confirm the artifact states NO expected residual value / no "residual is zero for the
correct labelling" (the specific rule-5 leak this remediation removed). Report any `assert` preceding the value it
guards, any conclusion emitted as an unconditional literal, any A−A control.

## Mathematica run discipline (in-prompt so both legs share it)
```
⛔ Wrap EVERY kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on.
⛔ NEVER raise the timeout; ⛔ never run more than one kernel at a time (the licence has TWO seats).
⛔ Copy the artifact to /tmp and ablate the COPY. ⛔ Never modify the working tree.
⭐ Save every ablation script AND its literal stdout to named absolute paths, and report those paths.
```

## Physics filter
Report a finding only if it catches a way the physics could be wrong; do not report "the script would be wrong on
a different input." Severity: MUST-FIX (a repaired control does not bite, the object is wrong physics, a leak, or
the protected core moved) vs NIT. If the two repairs are sound and the core is byte-identical, say so and name the
two or three things you checked most closely (with artifact + spec lines and your ablation stdout paths).
