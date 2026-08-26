# Physics consultation — S11c-a T-f: a background bulk current, symbolic vs hardcoded-zero

## What this is
S11c-a builds the interface shape-derivatives of a toy superfluid brane in two INDEPENDENT computer-algebra
engines (a SymPy engine and a blind Wolfram engine), each re-deriving the same physics from one shared spec.
On the T-f object (the dynamic-window projection identity and its shape derivative), the two engines treat the
**rest-frame background bulk current** `j⁰` differently. We are NOT asking you to confirm a decision we have
already made — we have not made one. We want your independent physics read, backed by a real CAS computation.
Do not assume either engine is correct, do not assume the shared spec is correct, and do not infer a preferred
answer from how this is phrased. If you think both engines are wrong, or the spec is wrong, say so.

## Read the source of truth first (open the files; quote what you rely on verbatim)
`research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`:
- §1 — the supplied rest-frame bulk fields. Note line ~48: *"`v_bulk_normal_0` is the steady bulk-normal drain
  and remains only the inherited rest-frame scope limit; the convective bulk problem is not reopened."*
- §1b — the current law and continuity: `j = ρ_4D v_bulk` , `∂_tρ_4D + ∇₄·j = 0`.
- §2d — the supplied background state `𝔅⁰` and its zeros (`V_s⁰ = 0 , J_s⁰ = 0 , 𝒜_s⁰ = 0`).
- §3c — the shifted-trace law and the two-argument dynamic window `𝒪`, and this sentence VERBATIM:
  *"In this scope the traced bulk velocity, the perturbation pressure, and the bulk current have zero
  background — `V_s⁰=J_s⁰=0` (§2d), the drain `v_bulk_normal_0` is the inert rest-frame scope limit of §1,
  `δp` has background value zero (§3b), and the rest-frame background current `ρ_4D⁰v_bulk⁰` vanishes — and
  the supplied density background depends on the in-plane anchor, not on `w`."*
  and earlier in §3c: *"Every background face value or normal derivative appearing in this law is obtained by
  differentiating a member of the supplied background state `𝔅⁰` (§2d); none may be introduced as a free
  premise."*

## The two engine treatments (UNTRUNCATED — but re-read the raw files yourself before relying on them)
A truncated/mis-copied engine form has manufactured a false conclusion in this area before. Open both files.

**Engine A — SymPy** (`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`), lines ~659-664:
```python
    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    density_background = rho4_bg_exact.subs(parameter, 0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )
```
So A constructs `j⁰ = ρ_4D⁰ · v_bulk⁰` and, having set `v_bulk⁰ = 0`, `j⁰` is identically zero from the start.
Its projection (lines ~1152-1168) uses only the perturbation current `j_bulk` (= `delta_j_bulk_i`) and its
divergence.

**Engine B — Wolfram** (`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`),
lines 434-446 and 798-822:
```wolfram
pressureZero[coordinates_List, normal_] := 0;
...
currentWZero[coordinates_List, normal_] := currentWBackground @@ Append[coordinates, {normal, time}];
currentXZero[index_][coordinates_List, normal_] := Symbol["currentXBackground" <> ToString[index]] @@ ...;
...
  currentX = Table[currentXZero[index][..] + waveOrder currentXWave[index][..], {index, spatialDirections}];
  currentW = currentWZero[..] + waveOrder currentWWave[..];
  leftSpace = Total[Table[D[Inactive[Integrate][window currentX[[index]], {normalCoordinate,-Inf,Inf}],
      spatialCoordinates[[index]]], {index, spatialDirections}]];
  windowSpace = Inactive[Integrate][Total[Table[D[window, spatialCoordinates[[index]]] currentX[[index]],
      {index, spatialDirections}]], {normalCoordinate,-Inf,Inf}];
  windowNormal = Inactive[Integrate][D[window, normalCoordinate] currentW, {normalCoordinate,-Inf,Inf}];
```
So B carries the background current as **independent nonzero symbols** `currentWBackground` /
`currentXBackground{i}`, decoupled from `v_bulk⁰`. It DOES set `pressureZero := 0`. The shape derivative used by
T-f is `shapeDerivative[e] := Together[Expand[D[e, waveOrder] /. waveOrder -> 0]]` (wl:45); for the dynamic
route `window` depends on `waveOrder` through the perturbed height, so the shape derivative of `window·j⁰` is
`(∂_waveOrder window)·j⁰`. (`projectionShapeDerivative`, wl:215-241, differentiates the integrand under the
integral — verify this yourself.)

## A prior measurement, for you to verify independently — do NOT take it as a premise
A prior run found that B's background-current symbols survive in the emitted `PROJECTION_SHAPE_DERIV` /
`PROJECTION_DYNAMIC_OPERAND` / `PROJECTION_RESIDUAL` (they are absent in `PROJECTION_STATIC_OPERAND`, which uses
a flat window with no `waveOrder`), i.e. they do not cancel among themselves in B's assembly. Re-derive this;
do not trust the claim.

## The questions — answer each; Q1 REQUIRES a CAS computation with literal stdout
The projection objects are integrals `∫ (…) dw` over the whole normal line `w ∈ (−∞,∞)` with the window `𝒪`
and its derivatives decaying at `±∞`, so integration by parts in `w` (dropping vanishing boundary terms) is a
legitimate identity on them.

**Q1 (COMPUTE, show script + literal stdout).** Reconstruct B's background-current contribution to the T-f
projection shape derivative (the `leftSpace`/`windowSpace`/`windowNormal` terms proportional to
`currentWBackground` / `currentXBackground{i}`, at `waveOrder → 0`). Then test whether that contribution
**vanishes as an integral identity** once you impose the background continuity relation implied by §1b at the
static background — `∇₄·j⁰ = ∂_{x1}j⁰_{x1}+∂_{x2}j⁰_{x2}+∂_{x3}j⁰_{x3}+∂_w j⁰_w = 0` (since `∂_tρ_4D⁰ = 0`) —
allowing integration by parts in `w`. Report the literal residual: does it reduce to zero, to a nonzero
integral, or to a pure boundary term? Save the script and its stdout to named absolute paths and report them.

**Q2.** Which is the more informative computation of whether the background current contributes to the projected
evolution law: A's (set `v_bulk⁰=0` at the outset, so `j⁰≡0` before the projection is computed), or B's (carry
`j⁰` symbolic and let the projection reveal whether it cancels)? State the tradeoff plainly.

**Q3.** §1 calls `v_bulk_normal_0` a *steady bulk-normal drain*; §3c calls it *"inert … not reopened"* and
declares `ρ_4D⁰v_bulk⁰` vanishes. Is `j⁰=0` a **fundamental** zero, or a **scope choice** that sets aside a real
(drain-sourced) background current? If the latter and the drain is nonzero, is B's surviving `∫ j⁰·δΩ` term a
genuine physical effect (a background drift coupling to the interface shape change), or an artifact?

**Q4.** Given Q1–Q3, what is the correct treatment of the background current in these two engines — set it to
zero (accept the scope choice), carry it symbolic in BOTH engines (and record whatever term survives), or
something else? What does each choice cost in information? If your Q1 residual is zero, does that change your
answer?

## Method
- Q1 must be a runnable CAS script (SymPy or Wolfram) with its literal stdout saved to named absolute paths you
  report. A prose-only derivation is discarded. Copy any repo file you execute to `/tmp`; never modify the tree.
- Physics filter: report only what bears on whether the physics is right. No style/naming points.
- ⛔ Do not tune anything to a recorded value or to the sibling engine. We want the computed truth, whatever it is.
