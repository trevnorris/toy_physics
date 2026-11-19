(* Independent current-first derivation. Integrate the constant spherical conormal
   load inward from the reservoir, then impose each incremental mouth law. *)
ClearAll["Global`*"];
check[name_, got_, want_] := Module[{value = Factor[got]},
  Print[name, " = ", ToString[value, InputForm]];
  If[!TrueQ[Simplify[value == want]], Exit[1]]
];
drop = Integrate[1/(4 Pi kap r^2), {r, ri, ro},
                 Assumptions -> ro > ri > 0 && kap > 0];
check["outer compliance", drop, (1/ri - 1/ro)/(4 Pi kap)];
trace = bath + force drop;
(* Solve for the load from one condition. Substitute it into the field to
   obtain the two readouts at the inner radius. *)
reading[constraint_] := Module[{f},
  f = force /. First[Solve[constraint == 0, force]];
  Factor /@ {trace /. force -> f, f}
];
pinned = reading[trace];
fluxHeld = reading[force];
(* Fixed bare source plus retained local energy km h^2/2. *)
sourceHeld = reading[force + km trace];
robin = reading[force + kr trace];
normalize[v_] := Factor /@ (v /. kap -> (1/ri - 1/ro) stiff/(4 Pi));
check["fixed bend pair per drive", normalize[pinned/bath], {0, -stiff}];
check["fixed conormal pair per drive", normalize[fluxHeld/bath], {1, 0}];
check["fixed source with local curvature per drive", normalize[sourceHeld/bath],
      {stiff/(stiff+km), -stiff km/(stiff+km)}];
check["fixed ideal source coincides with fixed conormal",
      (sourceHeld /. km -> 0) - fluxHeld, {0,0}];
check["mixed slope coincides with same local curvature",
      (robin /. kr -> km) - sourceHeld, {0,0}];
all = {pinned, fluxHeld, sourceHeld, robin};
check["opposite drive residual", Flatten[(all /. bath -> -bath) + all], ConstantArray[0,8]];
check["reserved negative half residual", Flatten[(all /. bath -> -bath/2) + all/2],
      ConstantArray[0,8]];
shifted = reading[force + kr (trace-bath)];
check["co-moving reference", shifted, {bath,0}];
Exit[0];
