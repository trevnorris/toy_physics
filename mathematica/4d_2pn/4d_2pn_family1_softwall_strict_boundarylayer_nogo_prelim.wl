ClearAll[section, pass, fail, recordCheck, Kexpr, objective, vars];

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[s_String] := Print["PASS: " <> s];
fail[s_String] := Print["FAIL: " <> s];
recordCheck[label_String, expr_] := Module[{v = TrueQ[expr]}, If[v, pass[label], fail[label]]];

section["1) Strict two-moment Family-1 boundary-layer ansatz"];

Print["Surface profile:"];
Print["  sigma(mu) = 1 - q mu^2 + r mu^4"];
Print["Strict surface action:"];
Print["  E = 1/2 ∫ dΩ [ A J(mu) Ψ^2 + B (Fθ(mu) (∂θΨ)^2 + Fφ(mu) (∂φΨ)^2/sin^2θ ) ]"];
Print["with J, Fθ, Fφ obtained by pulling the steep isotropic wall layer back to the reference mouth sphere."];

ClearAll[A, B, q, r];
J = Expand[-7 mu^8 r^2 + 6 mu^6 q r + 8 mu^6 r^2 - mu^4 q^2 - 8 mu^4 q r + 2 mu^4 r + 2 mu^2 q^2 - 2 mu^2 q + 1];
Ftheta = Expand[8 mu^8 r^2 - 8 mu^6 q r - 8 mu^6 r^2 + 2 mu^4 q^2 + 8 mu^4 q r - 2 mu^2 q^2 + 1];
Fphi = Expand[-8 mu^8 r^2 + 8 mu^6 q r + 8 mu^6 r^2 - 2 mu^4 q^2 - 8 mu^4 q r + 2 mu^2 q^2 + 1];

Print["Quadratic flare truncation:"];
Print["  J(mu)    = ", J];
Print["  Fθ(mu)   = ", Ftheta];
Print["  Fφ(mu)   = ", Fphi];

section["2) Channel stiffness formulas"];

Kexpr = <|
  {0, 0} -> Expand[7 A q^2/15 - 26 A q r/35 - 2 A q/3 + 23 A r^2/63 + 2 A r/5 + A],
  {1, 1} -> Expand[11 A q^2/35 - 2 A q r/5 - 2 A q/5 + 13 A r^2/77 + 6 A r/35 + A + 8 B q^2/35 - 32 B q r/105 + 32 B r^2/231 + 2 B],
  {1, 0} -> Expand[27 A q^2/35 - 10 A q r/7 - 6 A q/5 + 25 A r^2/33 + 6 A r/7 + A - 16 B q^2/35 + 64 B q r/105 - 64 B r^2/231 + 2 B],
  {2, 0} -> Expand[13 A q^2/21 - 94 A q r/77 - 22 A q/21 + 6185 A r^2/9009 + 6 A r/7 + A - 16 B q^2/7 + 320 B q r/77 - 320 B r^2/143 + 6 B],
  {2, 1} -> Expand[13 A q^2/21 - 230 A q r/231 - 6 A q/7 + 205 A r^2/429 + 10 A r/21 + A + 8 B q^2/21 - 96 B q r/77 + 800 B r^2/1001 + 6 B],
  {2, 2} -> Expand[5 A q^2/21 - 58 A q r/231 - 2 A q/7 + 25 A r^2/273 + 2 A r/21 + A + 16 B q^2/21 - 64 B q r/77 + 320 B r^2/1001 + 6 B]
|>;

targets = <|
  {0, 0} -> 4/45,
  {1, 1} -> 2/7,
  {1, 0} -> 1/4,
  {2, 0} -> 4/9,
  {2, 1} -> 2/3,
  {2, 2} -> 8/3
|>;

Do[
  Print["K_{", key[[1]], ",", key[[2]], "} = ", Kexpr[key]],
  {key, Keys[Kexpr]}
];

section["3) Numerical best-fit against passed support data"];

objective = Total[(Kexpr[#] - targets[#])^2 & /@ Keys[Kexpr]] // Expand;
vars = {A, B, q, r};

nmin = NMinimize[
  {
    objective,
    -10 <= A <= 10 && -10 <= B <= 10 && -10 <= q <= 10 && -10 <= r <= 10
  },
  vars,
  WorkingPrecision -> 40,
  Method -> "DifferentialEvolution"
];

bestObj = First[nmin];
bestRules = Last[nmin];

Print["Best-fit rules:"];
Print[bestRules];
Print["\nBest objective = ", N[bestObj, 20]];

residuals = Association@Table[
  key -> N[(Kexpr[key] - targets[key]) /. bestRules, 20],
  {key, Keys[Kexpr]}
];

maxAbs = Max[Abs[Values[residuals]]];
Print["Max absolute channel residual = ", N[maxAbs, 20]];
Print["\nChannel residuals:"];
Do[
  Print["  ", key, " -> ", residuals[key]],
  {key, Keys[residuals]}
];

recordCheck["Strict two-moment model does not exactly fit the passed support data (numeric obstruction)", bestObj > 10^-6];

section["4) Conclusion"];

Print["Main result:"];
Print["- The strict two-moment isotropic boundary-layer pullback does NOT exactly reproduce the passed support-channel data."];
Print["- The best global fit still leaves a sizable positive residual."];
Print[""];
Print["Interpretation:"];
Print["- Pure isotropic penetration + pure geometry pullback is too small a model."];
Print["- At least one extra local wall-stress / traction / profile degree of freedom is required."];
Print["- This makes the reduced variational wall-Hessian a genuine new structure, not just a reparameterization of the strict isotropic layer."];
