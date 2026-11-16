(* Separate Cartesian construction of density samples across a spherical shell
   of aperture points. Its radius r lies inside the three-dimensional mouth. *)
ClearAll["Global`*"];
cartDensity[x_, w_, s_] := base + g w + h s x/r + e (s x/r)^2 + k (s x/r)^3;
height[x_, s_] := eta s x/r;
read[x_, s_] := cartDensity[x, height[x, s], s];
diff[x_, s_] := read[x, s] - read[x, 0];
(* For surface area on the 2-sphere, an axial band dx has area 2 pi r dx. *)
band = 2 Pi r;
assumptions = r > 0 && s > 0;
project = FullSimplify[
    Integrate[band (x/r) diff[x, s], {x, -r, r}]/(s Integrate[band (x/r)^2, {x, -r, r}]),
    assumptions];
poles = Simplify[(diff[r, s] - diff[-r, s])/(2 s)];
Print["shell projection and opposite poles = ", InputForm[{project, poles}]];
Print["full and half speed projection = ", InputForm[{project, project /. s -> s/2}]];
limits = Limit[{project, poles}, s -> 0];
Print["zero-drive limits = ", InputForm[limits]];
fixed = D[cartDensity[x, 0, s], s] /. s -> 0;
moving = D[read[x, s], s] /. s -> 0;
Print["fixed and registered linear fields = ", InputForm[{fixed, moving}]];

(* Recover the same scalar from an independently supplied number current and
   a nonzero relative normal speed, varying both before taking the derivative. *)
inferred = D[(j0 + s j1)/(v0 + s v1), s] /. s -> 0;
Print["quotient response = ", InputForm[inferred]];
current = read[x, s] (v0 + s v1);
reconstructed = Simplify[inferred /. {j0 -> (current /. s -> 0),
                                      j1 -> (D[current, s] /. s -> 0)}];
Print["current reconstruction and difference = ", InputForm[{reconstructed, Simplify[reconstructed-moving]}]];
If[!TrueQ[Simplify[limits == {moving, moving} /. x -> r]] || Simplify[reconstructed-moving] =!= 0,
   Exit[1]];
