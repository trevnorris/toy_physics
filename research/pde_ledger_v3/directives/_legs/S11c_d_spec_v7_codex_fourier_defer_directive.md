# FOCUSED refix — S11c-d SHARED PHYSICS §1a/§1c/§7 Fourier, to v7 (defer the mechanics to the build directive)

You (`gpt-5.6-sol`) are making a **narrow, surgical** correction to an otherwise-SOUND physics spec. The
Fourier-reduction description has now failed **three** review gates (v4 supplied a wrong `(2π)⁻³` map; v5 was
non-executable + leaked; v6 misdescribes the real kernel's Fourier content and does not wire the comparator). The
decision (user-approved) is to stop enumerating the kernel's Fourier micro-structure in the spec and instead **defer
the Fourier mechanics to the build directive** — exactly as §7 already defers the `IMPORT_KEYS` root set — while
keeping the OBJECT, the supplied `f̂_red` convention, and the leak-safe CONTROL in the spec.

⛔ **Edit ONLY** the §1a Fourier-content sentence(s), the §1c Fourier-carrier passage, the §7 comparator residual list,
and (if needed) the §8 line — plus the version label. ⛔ **Do NOT touch anything else** — every other part of the spec
is SOUND across many legs.

## The file
`research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`

## Verified facts about the REAL consume-set (ground the edit in these)
- `s11cc2ClosedCouplingKernel` carries **0** `DiracDelta` (c2 already applied the momentum-conservation identity and
  **strips** the deltas — `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:373-374`). The 3-D momentum deltas live on
  **c1's** `dtn_kernel` (`FLAT_DIAGONAL`); ⛔ they are **not** the 3-D→1-D operand and ⛔ `dtn_kernel` must **not** be
  bound to supply the reduction convention (binding it — delta without `(2π)³` — over a kernel whose Integrals already
  carry `π⁻³ = 1/(2π)³` double-counts `(2π)³`).
- The kernel's Fourier content is **not** only the transfer `Q = k_out − k_in`: it also carries **middle-leg** symbols
  (`ŵ(k_out − k_mid)`, `ŵ(k_mid − k_in)`; thousands of `s11cc2MiddleMomentum` occurrences) — the **required three-leg
  second-scattering** at the mixed retained grade (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:398-400`: "Second
  scattering is required at the mixed retained grade. A three-leg triangular representation … its middle leg
  integrated below. There is no single-momentum division.") — plus in-plane **position integrals** (over
  `s11cc2Y{1,2,3}`) already carrying `π⁻³`. Dropping the middle-leg symbols drops the retained `ησ_W` grade
  (S11c-a `:195-198` requires first order in **each** background bookkeeper, including the mixed grade).

## The fix (keep the control + `f̂_red`; make §1a accurate; DEFER the enumeration/`(2π)` to the build directive; wire §7)

**§1a — make the Fourier-content description ACCURATE, ⛔ not exhaustive.** Replace the current
"its Fourier content is on the momentum transfer `Q = k_out − k_in`, ⛔ not a single unspecified `k`" with an accurate
statement: the closed kernel's Fourier content appears as the profile-transfer carrier **and** the **middle-leg**
symbols of the required three-leg second-scattering, over the kernel's **in-plane position integrals**, and it carries
**no** 3-D momentum deltas (c2 stripped them; the `dtn_kernel` deltas are **not** the 3-D→1-D operand). ⛔ The exact
enumeration of which Fourier symbols appear (transfer vs middle-leg), the `(2π)` bookkeeping, and the operand
construction are **the build directive's** (against the real kernel, its own two decision legs), ⛔ not enumerated
here. (You may keep `dtn_kernel` in the reachable-import list **only** if you add, right there, that it is ⛔ not bound
as the reduction-convention operand.)

**§1c — keep the leak-safe control + `f̂_red`; DEFER the mechanics.** Keep: the supplied 1-D convention `f̂_red` (the
only supplied Fourier convention); the requirement that **each engine COMPUTE and EMIT the 3-D→1-D reduction of its
own closed-kernel Fourier content — every such symbol (transfer and middle-leg) plus the already-present in-plane
integrals — as an object with BOTH operands, using its own realized Fourier convention** (⛔ not a typed `[L_W/(2π)]`
or `(2π)²L_W` map, ⛔ not by binding `dtn_kernel`); the comparator joins the two engines' **reduced kernels**; the
anti-tautology point (⛔ `A_3D ≡ [L_W/(2π)]δ²A_edge` is `A−A`); and the flux-`δ²(Q_∥)`-stripping point. **Add** one
explicit deferral clause: the exact Fourier-symbol set, the `(2π)`/dimensional bookkeeping, and the operand
construction are fixed at the **build directive** against the real kernel (its own two decision legs verify them),
⛔ not typed in this spec — the spec supplies only `f̂_red`, the reduce-your-own-kernel requirement, and the
comparator-join control. ⛔ Remove any residual sentence that names "the two-momentum identity" in a way that invites
binding `dtn_kernel`, or that implies the transfer `Q = k_out − k_in` is the only Fourier content.

**§7 — wire the control into the comparator.** Add the **both-operand reduced-kernel join** to the load-bearing
comparator residual list (alongside the S-matrix / currents / conversion / survival / Riesz joins), so the §1c control
is not lost at comparator authoring. State the join key includes the profile-transfer and middle-leg reduced symbols.

**§8 — keep the split** already present: `f̂_red` SUPPLIED; each engine's 3-D→1-D reduction factor (over transfer +
middle-leg + integrals) and reconstruction COMPUTED.

## Also
- **Bump the version label to v7** (header line 1 and the two "v6"/"Spec v6" lines), extending the lineage note (…→ v6
  Codex Fourier refix → v7 Codex Fourier-defer, per Grok's round-6 F1).

## Governing (unchanged, do not violate)
`M2`: name the object; the only SUPPLIED Fourier convention is `f̂_red`; the leak-safety is the CONTROL (each engine
its own convention + comparator joins the reduced kernels), ⛔ never a supplied `(2π)`/`δ²` map. Blindness: ⛔ no
spec-level pointer into any engine's construction script; the blind WL re-derives from specs; ⛔ do not bind
`dtn_kernel` as the convention operand. Keep the c2-honesty premises, the settled frame, and every other section
**byte-for-byte unless the fix above requires the edit**.

## Output
Edit only the passages named above. Print a short report: the §1a change (accurate + deferral), the §1c change (kept
control + added deferral + removed the two-momentum-identity/only-transfer wording), the §7 comparator-join addition,
and the version-label change. ⛔ No other file edits; ⛔ do not reprint the whole spec.
