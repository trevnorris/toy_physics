# Independent physics adjudication — S11c-a T-f: the background current in the dynamic-window projection

## Your task
S11c-a builds the interface shape-derivatives of a toy superfluid brane in two INDEPENDENT computer-algebra
engines (a SymPy engine and a blind Wolfram engine), each re-deriving the same physics from one shared spec.
On the T-f object (the dynamic-window projection identity and its shape derivative), the two engines treat
the **background bulk current** differently. Your job is to determine, **from the spec alone**, which
treatment is faithful — or whether both are, or neither. Derive it yourself; do not assume either engine is
correct, and do not try to guess a "expected" answer from how the question is phrased.

## Read the source of truth FIRST, before the engine forms
Open and read these, form your own view of what the spec requires, and only then read the engine treatments
below:
- `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`
  - §1 (the supplied rest-frame bulk fields and the current law `j = ρ_4D v_bulk`, `∂_tρ_4D + ∇₄·j = 0`),
  - §2d (the supplied background state `𝔅⁰`; what has zero background),
  - §3c (the shifted-trace law AND the two-argument dynamic window; the sentence beginning "In this scope
    the traced bulk velocity, the perturbation pressure, and the bulk current have zero background …"),
  - §4 T-f (`S11CA_PROJECTION_*` — the projection identity, its dynamic-window and static-flat routes, and
    its shape derivative).
- Quote the exact spec sentences you rely on. Do not paraphrase a premise you are going to lean on.

## The exact engine forms (UNTRUNCATED — but re-read the raw files yourself before comparing)
Do NOT trust these transcriptions as premises; open both engine files and read the raw definitions at the
cited lines. A truncated or mis-copied engine form has manufactured a false finding here before.

**Engine A — SymPy** (`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`), the traced
bulk fields, lines ~660-683:
```python
    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    density_background = rho4_bg_exact.subs(parameter, 0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )
    ...
    current_perturbation = tuple(
        affine_bulk_perturbation(j_bulk[i], dw_delta_j_bulk[face][i], face)
        for i in range(4)
    )
```
`j_bulk[i]` are the symbols `delta_j_bulk_{i}` (perturbation current). The projection (lines ~1150-1168)
integrates the window against `j_bulk` and `current_divergence = Σ_i grad_j_bulk[i][i]` only.

**Engine B — Wolfram** (`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`),
the current fields, lines 434-446:
```wolfram
pressureZero[coordinates_List, normal_] := 0;
...
currentWZero[coordinates_List, normal_] :=
  currentWBackground @@ Append[coordinates, {normal, time}];
currentWWave[coordinates_List, normal_] :=
  currentWPerturbation @@ Append[coordinates, {normal, time}];
currentXZero[index_][coordinates_List, normal_] :=
  Symbol["currentXBackground" <> ToString[index]] @@ Append[coordinates, {normal, time}];
currentXWave[index_][coordinates_List, normal_] :=
  Symbol["currentXPerturbation" <> ToString[index]] @@ Append[coordinates, {normal, time}];
```
and the projection assembly, lines ~798-822:
```wolfram
  currentX = Table[currentXZero[index][spatialCoordinates, normalCoordinate] +
      waveOrder currentXWave[index][spatialCoordinates, normalCoordinate], {index, spatialDirections}];
  currentW = currentWZero[spatialCoordinates, normalCoordinate] +
    waveOrder currentWWave[spatialCoordinates, normalCoordinate];
  leftSpace = Total[Table[D[Inactive[Integrate][window currentX[[index]],
    {normalCoordinate, -Infinity, Infinity}], spatialCoordinates[[index]]], {index, spatialDirections}]];
  windowSpace = Inactive[Integrate][Total[Table[
    D[window, spatialCoordinates[[index]]] currentX[[index]], {index, spatialDirections}]], {normalCoordinate, ...}];
  windowNormal = Inactive[Integrate][D[window, normalCoordinate] currentW, {normalCoordinate, ...}];
```
The shape derivative used by T-f is `shapeDerivative[e] := Together[Expand[D[e, waveOrder] /. waveOrder -> 0]]`
(the coefficient of the first power of `waveOrder`). `window` for the dynamic route depends on `waveOrder`
through the perturbed height; for the static route it is `staticFlatWindow` (no `waveOrder`).

Note the two engines name the current differently: A's `delta_j_bulk_i` is B's `currentXPerturbation{i}` /
`currentWPerturbation` (the perturbation current). The question below is about the OTHER piece — the
background part (`currentXZero`/`currentWZero` = `currentXBackground`/`currentWBackground` in B; the
`current_background` tuple in A).

## The questions to answer FROM THE SPEC
1. In the T-c…T-i scope of S11c-a (which includes T-f), does the spec supply the rest-frame background
   current `ρ_4D⁰v_bulk⁰` as **zero**, or as a nonzero independent field? Cite the exact sentence.
2. Given your answer to (1): in the shape derivative of the dynamic-window projection, should any term
   proportional to the **background** current (independent of the perturbation current) survive? Show the
   computation. In particular: if the background current is zero, does `∂_waveOrder [ window(waveOrder) ·
   j⁰ ] |_{waveOrder=0}` vanish? If it is nonzero, what term does it leave?
3. Therefore, which engine treatment is faithful to the spec on this point — A (background current
   constructed as `ρ_4D⁰·v_bulk⁰` with `v_bulk⁰=0`, hence zero), B (background current carried as an
   independent nonzero symbol `currentWBackground`/`currentXBackground`), both, or neither? State it plainly.

## Method requirements
- If you make any derivation claim, back it with a runnable CAS script AND its literal stdout, saved to named
  absolute paths that you report. A prose derivation with no script is discarded (this project's standing
  rule). Where the claim is purely a spec-reading, quote the verbatim spec text instead.
- Copy any file you execute to `/tmp` and work on the copy; never modify the working tree.
- Physics filter: report a conclusion only if it bears on whether the physics is right — i.e. whether a
  background-current term that the spec does (or does not) zero would appear in the emitted projection
  object. Do not report style, naming, or "it would be wrong on a different input."

## What to return
- Your spec-grounded answer to (1), (2), (3), with the verbatim spec citations and your script paths+stdout.
- If you conclude one engine diverges from the spec, name which, name the exact term that should not be
  there (or that is missing), and cite the spec line that decides it.
