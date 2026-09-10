# FOCUSED refix — S11c-d SHARED PHYSICS §1c Fourier carrier + §8, to v6 (rule-15, the material keeps breeding)

You (`gpt-5.6-sol`) are making a **narrow, surgical** correction to an otherwise-SOUND physics spec. The §1c
Fourier-carrier passage has failed a review gate **twice** (v4 supplied a wrong `(2π)⁻³` map; v5's "compute from the
kernel's `DiracDelta³` piece" is non-executable and still leaks). ⛔ **Edit ONLY** the §1c Fourier-carrier passage, the
§8 supplied/computed lines about that convention, and the version label. ⛔ **Do NOT touch anything else** — every
other part of the spec is SOUND across both legs, and any collateral edit will force a needless re-review.

## The file
`research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`

## What is wrong now (round-5 Grok F1, orchestrator-VERIFIED against the real sources)
The §1c block beginning **"⚠⚠ The imported c2 carrier's normalization is NOT supplied — it must be COMPUTED…"** through
the **"Reconstruction round-trips against the DERIVED map…"** paragraph is defective:
1. It tells engines to reduce the insertion "against the imported kernel's own flat identity (its
   `DiracDelta³(k_out−k_in)` piece)" — but **`DiracDelta` does not occur in `scripts/S11c_c2_exports.py` (count 0)**.
   c2 **strips** the deltas before folding (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:373-374`:
   `deltas = diagonal.atoms(sp.DiracDelta); z0out = diagonal.xreplace({d: 1 …})`). The three-delta identity lives on
   **c1's** `dtn_kernel` (`scripts/S11c_c1_exports.py`), NOT on the consumed `s11cc2ClosedCouplingKernel`. So the
   prescribed computation **cannot execute** on the consume-set.
2. The sentence "The c2 engine's transform convention is an **unnormalised forward transform with a normalised
   inverse** … the `(2π)⁻³` sits on the inverse …" — together with the supplied `f̂_red` and the geometry
   `∂_{yᵢ}f = n̂ᵢ f′/L_W` — **uniquely fixes the 3-D→1-D factor without looking at any kernel**. That is the very
   numeric-convention **leak** the fold claimed to remove (`M2`).
3. It points the **blind Wolfram engine** at `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (the OTHER engine's
   construction script). There is **no WL self-energy engine**; WL must re-derive from the sibling **specs** only
   (`S9_export_chain_rebuild_directive.md:16-18`). Neither `S11c_c1_SHARED_PHYSICS.md` nor `S11c_c2_SHARED_PHYSICS.md`
   states any `2π`/`DiracDelta` convention (count 0/0) — it is **not shared physics**. A spec-level pointer into a
   construction script is the opposite of the blindness control.
4. **§8 contradicts §1c**: §8 lists "exact reduced Fourier/3-D carrier convention" under **SUPPLIED**
   (unfalsifiable), while §1c says the c2 carrier normalization is **not** supplied.

## The fix (Grok's exact minimal fix — implement this, ⛔ add no cleverness)
Rewrite the §1c Fourier-carrier passage so that:
- **`f̂_red` stays SUPPLIED.** Keep the coordinate/`Q`/`s` definitions and the reduced-transform definition
  `f̂_red(s) ≡ ∫dξ e^{−isξ}f(ξ)` and its inverse — that is S11c-d's own **1-D convention**, a named object, legitimately
  supplied. Keep the "distributional for a full step / coordinate-space or asymptotic-subtraction" caveat.
- **DELETE** the "unnormalised forward transform / normalised inverse / `(2π)⁻³`" sentence **and** the pointer to
  `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`. ⛔ Supply **no** numeric `(2π)`/`δ²(Q_∥)`/`[L_W/(2π)]` map for the
  imported carrier.
- ⛔ **Do NOT** tell engines to read a `DiracDelta³` piece out of `s11cc2ClosedCouplingKernel` (it is not there).
- **Require each engine to EMIT the 3-D→1-D reduction of ITS OWN closed-kernel Fourier symbols as a computed object
  with BOTH operands:** (i) the 3-D carrier **as it actually appears in that engine's own closed coupling kernel**, and
  (ii) the reduced 1-D kernel obtained by applying **that same engine's own Fourier convention** (the one realized in
  its own construction of the two-momentum identity / profile insertion — ⛔ not a typed `[L_W/(2π)]` or `(2π)²L_W`) to
  the §1c interface geometry `f = f(n̂·y/L_W)`. **The comparator joins those two engines' reduced kernels** (⛔ not a
  pre-factored coefficient). The `2π`, the `δ²(Q_∥)`, and the dimensional content **fall out of each engine's own
  computation** — the cross-engine agreement (or residual) on the reduced kernel is the measurement.
- Keep the anti-tautology point (⛔ `A_3D ≡ [L_W/(2π)]δ²A_edge` is `A−A`) and the flux point (strip the tangential
  `δ²(Q_∥)` before squaring) — but phrase them so they do **not** reintroduce a supplied `[L_W/(2π)]` map; each engine's
  reduced kernel and its reconstruction are its own computed objects.
- **Optional witness:** if a build later uses the **AGREE'd c1 `dtn_kernel`** flat identity as a convention witness,
  name **`dtn_kernel`** (⛔ not `s11cc2ClosedCouplingKernel`), and still require the 1-D interface-geometry step to be
  **computed**, not typed. (You may mention this as the legitimate route; the exact use is the build directive's.)
- **§8:** move the Fourier item so that **`f̂_red` (the 1-D reduced-transform convention) is SUPPLIED**, while **the
  3-D→1-D reduction factor and the 3-D↔1-D reconstruction are COMPUTED** (per engine). Fix the current wording that
  lists "exact reduced Fourier/3-D carrier convention" as SUPPLIED.

## Also
- **Bump the version label to v6.** The header (line 1) and the two "v4"/"Spec v4" lines (~16, ~21) currently say v4;
  set them to v6 and add one clause noting the lineage: Codex-authored v4 → orchestrator §1c/nit fold (v5) →
  this Codex §1c-Fourier refix (v6, rule 15, per Grok's round-5 F1). Keep the house §0–§8 structure.
- ⛔ Leave the §7 consume-set as is (the `s11cc2MiddleMomentum*`/`s11cc2OutgoingNormalMomentum` keys are correctly
  deferred to the build directive's `IMPORT_KEYS` root set — not a defect).

## Governing (unchanged, do not violate)
`M2`: name the object, supply no expected value; the only SUPPLIED convention here is the 1-D `f̂_red` definition.
Blindness: ⛔ no spec-level pointer into any engine's construction script; the blind WL re-derives from specs. Keep the
c2-honesty premises, the settled frame, and every other section **byte-for-byte unless the fix above requires the
edit**.

## Output
Edit only the passages named above. Print a short report: the exact §1c sentences deleted, the replacement's key
requirement (per-engine self-reduction + comparator-join), the §8 change, and the version-label change. ⛔ No other
file edits; ⛔ do not reprint the whole spec.
