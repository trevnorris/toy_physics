(* Fixed, round electromagnetic and Berry connections on an isolated sphere.
   All angular averages below use the unit sphere. Energy is in units
   B = 2 Pi R^2 alpha^2/u, fields in units Sqrt[alpha/u]. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];

(* For signed difference degree k the monopole spectrum is
   j(j+1) - (k/2)^2, j = |k|/2 + ell, ell a nonnegative integer.
   Use kAbs >= 0 to display the dependence on the excitation number. *)
angularMomentum = kAbs/2 + ell;
spectrum = Expand[angularMomentum (angularMomentum + 1) - (kAbs/2)^2];
bottom = spectrum /. ell -> 0;
emit["sphere spectrum", spectrum];
emit["excitation above bottom", Factor[spectrum - bottom]];
native[h_, s_] := -2 h s;
mismatch[n_, h_, s_] := n - native[h, s];
bestDegree[n_, s_] := Min[Abs[mismatch[n, -1, s]], Abs[mismatch[n, 1, s]]];
degreeIdentity = Resolve[ForAll[n, bestDegree[n, 1] == Abs[Abs[n] - 2]], Reals];
mirrorDifference = Expand[mismatch[n, h, -s] - mismatch[n, -h, s]];
conjugateDifference = Expand[mismatch[-n, h, s] + mismatch[n, -h, s]];
emit["minimum difference degree equals Abs[Abs[n]-2]", degreeIdentity];
emit["mirror exchange and sign conjugate-exchange residuals",
     {mirrorDifference, conjugateDifference}];

(* m is the total mean density of both helicities. The spectral estimate
   bounds the sum of both kinetic energies, and Jensen bounds the quartic
   average by m^2 before minimizing with m >= 0. *)
lowerMaterial = 4 x bottom m - 2 m + m^2;
vertex = m /. First[Solve[D[lowerMaterial, m] == 0, m]];
emit["Jensen polynomial and unconstrained density vertex", {lowerMaterial, vertex}];
lower[n_, degree_, x_, y_] := y n^2 - Max[0, 1 - x degree]^2;
incoherent[n_, y_] := y n^2;
emit["incoherent energy", incoherent[n, y]];
emit["matched n=2 constant field energy", lower[2, 0, x, y]];
emit["neutral flat-native constant field energy", lower[0, 0, x, y]];
(* For C_h=0 one has degree=|n|, so any charged state lies above -1.
   With native degrees +/-2, every |n| >= 2 lies above or at 4 y - 1.
   The matched degree-zero connection at |n|=2 realizes this lower bound. *)
allChargedMargin = Resolve[ForAll[{n, y}, y > 0 && Abs[n] >= 1,
                                  y n^2 - 1 - (-1) > 0], Reals];
tailMargin = Resolve[ForAll[{n, y}, y > 0 && Abs[n] >= 2,
                            y n^2 - 1 >= 4 y - 1], Reals];
emit["flat-native charged and native |n|>=2 universal comparisons",
     {allChargedMargin, tailMargin}];

(* For positive unit degree use the regular northern gauge A_phi=1-z.
   The southern gauge differs by d phi and the section has a simple zero
   at z=0. Its complex conjugate gives the negative-degree trial. The
   axisymmetric angular average is Integrate[.,{z,0,1}]. *)
profile = Sqrt[z] (1 + t (1 - z));
sinSquared = 4 z (1 - z);
radialMetric = sinSquared/4;
azimuthPotential = 1 - z;
gradDensity = Together[radialMetric D[profile, z]^2 +
                       azimuthPotential^2 profile^2/sinSquared];
average[expr_] := Expand[Integrate[Expand[expr], {z, 0, 1}]];
norm = average[profile^2];
gradient = average[gradDensity];
quartic = average[profile^4];
emit["unit trial {norm,angular kinetic,quartic}", {norm, gradient, quartic}];
emit["kinetic excess over lowest harmonic", Factor[gradient - norm/2]];
(* a is the squared overall amplitude. A polynomial trial can explore
   higher harmonics; minimize over nonnegative a at each fixed shape. *)
trialPolynomial = y + a (4 x gradient - 2 norm) + a^2 quartic;
trialVertex = Together[a /. First[Solve[D[trialPolynomial, a] == 0, a]]];
energyAtVertex = Factor[trialPolynomial /. a -> trialVertex];
emit["unit-trial constrained optimum amplitude squared",
     Max[0, trialVertex]];

xAt[r_] := 1/(3 r^2);
yAt[r_] := 2/(9 r^3);
shape = 1/2;
trialAt[r_, testShape_] := Module[{rules = {x -> xAt[r], y -> yAt[r], t -> testShape}},
  trialPolynomial /. a -> Max[0, trialVertex /. rules] /. rules
];
referenceRadii = {1, 101/100};
comparison[r_] := Module[{xVal = xAt[r], yVal = yAt[r], unit, zero, tail},
  unit = trialAt[r, shape];
  zero = lower[0, bestDegree[0, 1], xVal, yVal];
  tail = lower[2, bestDegree[2, 1], xVal, yVal];
  {r, xVal, yVal, unit, zero, tail, Min[zero, tail] - unit}
];
emit["shape t", shape];
emit["lowest-harmonic trial and positive-shape trial at R=1",
     {trialAt[1, 0], trialAt[1, shape]}];
rows = comparison /@ referenceRadii;
Scan[emit["{R,x,y,unit trial,neutral lower,|n|>=2 lower,certified gap}", #] &, rows];

(* Verify strict separation on the entire preparation interval via rational
   inequalities, retaining the same shape even as the optimum amplitude
   adjusts. Check both active branches before using their polynomials. *)
intervalRules = {x -> xAt[r], y -> yAt[r], t -> shape};
unitOnBranch = energyAtVertex /. intervalRules;
intervalTests = {
  (trialVertex /. intervalRules) > 0, 1 - 2 xAt[r] > 0,
  unitOnBranch < -(1 - 2 xAt[r])^2, unitOnBranch < 4 yAt[r] - 1
};
fullInterval = Resolve[ForAll[r, 1 <= r <= 101/100, And @@ intervalTests], Reals];
emit["closed preparation interval strict comparisons", fullInterval];

(* Fixed-material size controls. At small radius every charged sector has
   the common lower bound y-1; at large radius the matched pair is exact.
   For the remaining |n|>=3 use the universal lower bound 9y-1. *)
small = 1/2;
emit["R=1/2 {neutral lower,neutral incoherent,all charged lower}",
  {lower[0, bestDegree[0, 1], xAt[small], yAt[small]], incoherent[0, yAt[small]],
   yAt[small] - 1}];
large = 2;
largeBounds = {
  lower[2, 0, xAt[large], yAt[large]],
  lower[0, bestDegree[0, 1], xAt[large], yAt[large]],
  lower[1, bestDegree[1, 1], xAt[large], yAt[large]],
  9 yAt[large] - 1
};
emit["R=2 {matched |n|=2,neutral lower,unit lower,|n|>=3 lower}", largeBounds];

(* A separately postulated local splitting, with strongHelicity defined
   relative to the chosen local normal. Both signs and both channels remain
   admissible. For any allocation of the total mean density, weighting by the
   largest quadratic gain gives a lower bound, including mixed sections. *)
helicities = {-1, 1};
alphaRatio[h_, strongHelicity_] := If[h === strongHelicity, 1, 4/5];
splitLower[n_, s_, strongHelicity_, x_, y_] :=
  y n^2 - Max[0, Table[alphaRatio[h, strongHelicity] -
                       x Abs[mismatch[n, h, s]], {h, helicities}]]^2;
splitComparison[r_, s_, strongHelicity_] := Module[
  {selected = Sign[native[strongHelicity, s]], unit, opposite, zero, tail},
  unit = trialAt[r, shape];
  opposite = splitLower[-selected, s, strongHelicity, xAt[r], yAt[r]];
  zero = splitLower[0, s, strongHelicity, xAt[r], yAt[r]];
  tail = 4 yAt[r] - 1;
  {r, s, strongHelicity, selected, unit, opposite, zero, tail,
   Min[opposite, zero, tail] - unit}
];
splitRows = Flatten[
  Table[splitComparison[radius, sgn, strong], {radius, referenceRadii},
        {sgn, {-1, 1}}, {strong, {-1, 1}}], 2];
Scan[emit["split {R,s,strong h,selected n,trial,opposite lower,neutral lower,tail lower,gap}",
          #] &, splitRows];
(* The two pointwise symmetries act on signed differences as follows.
   Fixed weights: (n,s,h)->(-n,-s,h) conjugates each section.
   Swapped weights: (n,s,h)->(n,-s,-h) relabels the two sections. *)
fixedWeightMirror = Expand[mismatch[-n, h, -s] + mismatch[n, h, s]];
swappedWeightMirror = Expand[mismatch[n, -h, -s] - mismatch[n, h, s]];
swappedRatios = Table[alphaRatio[-h, -strong] - alphaRatio[h, strong],
                      {h, helicities}, {strong, helicities}];
emit["split mirror {fixed weights conjugation,swapped weights relabel,ratio residuals}",
     {fixedWeightMirror, swappedWeightMirror, swappedRatios}];
emit["split flat-native n=0 energy, charged lower excess",
     {-Max[alphaRatio[1, 1], alphaRatio[-1, 1]]^2, y n^2}];
checks = {degreeIdentity, mirrorDifference == 0, conjugateDifference == 0,
          allChargedMargin, tailMargin, fullInterval,
          And @@ (Positive[Last[#]] & /@ rows),
          And @@ (Positive[Last[#]] & /@ splitRows),
          fixedWeightMirror == 0, swappedWeightMirror == 0,
          swappedRatios == ConstantArray[0, {2, 2}],
          First[largeBounds] < Min[Rest[largeBounds]]};
If[!TrueQ[And @@ checks], Exit[1]];
