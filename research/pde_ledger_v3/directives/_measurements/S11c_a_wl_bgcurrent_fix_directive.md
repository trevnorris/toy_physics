# Measurements — S11c_a_wl_bgcurrent_fix_directive

Every factual claim in the directive carries the command that produced it (rule 2). Literal outputs below.
(The physics justification for *why* zeroing is safe — the continuity-cancellation computation — deliberately
lives in `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §16c, NOT here: it is a downstream result the builder must
not see or engineer toward. §3c: "Which trace terms then survive is computed from these premises, not stated
here.")

## Claim: WL introduces the background current as free-premise symbols (defs at 438-446)
```
$ grep -n "currentWBackground\|currentXBackground\|currentWZero\|currentXZero" \
    research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl
438:currentWZero[coordinates_List, normal_] :=
439:  currentWBackground @@ Append[coordinates, {normal, time}];
442:currentXZero[index_][coordinates_List, normal_] :=
443:  Symbol["currentXBackground" <> ToString[index]] @@
615:  "NORMAL_CURRENT" -> {currentWZero, currentWWave, dimensionFlux},
616:  "CURRENT_X1" -> {currentXZero[1], currentXWave[1], dimensionFlux},
617:  "CURRENT_X2" -> {currentXZero[2], currentXWave[2], dimensionFlux},
618:  "CURRENT_X3" -> {currentXZero[3], currentXWave[3], dimensionFlux},
801:    currentXZero[index][spatialCoordinates, normalCoordinate] +
804:  currentW = currentWZero[spatialCoordinates, normalCoordinate] +
```
⇒ only defs (438-446), inventory register (615-618), projection assembly (801-804); NO rule sets them to 0.

## Claim: WL already sets its background bulk velocity to zero by construction (115-119)
```
$ sed -n '115,119p' research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl
bulkVelocityZero[coordinates_List, normalCoordinate_] := Through[
  {zeroVelocityOne, zeroVelocityTwo, zeroVelocityThree, zeroVelocityW}[
      Sequence @@ Append[coordinates, {normalCoordinate, time}]]] /.
    {zeroVelocityOne[___] -> 0, zeroVelocityTwo[___] -> 0,
      zeroVelocityThree[___] -> 0, zeroVelocityW[___] -> 0};
```
⇒ every component (including normal/drain `zeroVelocityW`) → 0. So `j⁰ = ρ⁰·v⁰` built from this is zero.

## Claim: the SymPy engine constructs j⁰ = ρ⁰·v⁰ (=0) — the property WL must mirror (660-664)
```
$ sed -n '660,664p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    density_background = rho4_bg_exact.subs(parameter, 0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )
```

## Claim: WL's free-symbol background current SURVIVES in the emitted T-f physics tags (motivates the fix)
```
$ python3 ~/.s11_build/S11c_a_bgcurrent_check.py
TAG                          cases bg_cases bg_hits pert_hits
PROJECTION_SHAPE_DERIV           8        8    1660       344
PROJECTION_STATIC_OPERAND        8        0       0        40
PROJECTION_DYNAMIC_OPERAND       8        8    1660       344
PROJECTION_RESIDUAL              8        8    1660       384
```
⇒ WL carries free background-current symbols into the emitted current-consuming objects (present in the
dynamic route, absent in the flat-static route). The fix removes the free premise at its source.
