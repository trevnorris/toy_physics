(* Ordered round collar, r in [R,lambda R]. The radial and angular fields
   are unrestricted in the energy space, with free material edges. There
   is no coherence penalty in the exterior. Only the trial is constant
   in radius. n is an integer varied across sectors. Both helicities are
   available at each mouth. *)
ClearAll["Global`*"];
emit[label_, data_] := Print[StringJoin[label, " = ", ToString[data, InputForm]]];

(* tRad = r/R, f = R Sqrt[u/(K g)] Psi, and g is the stronger coefficient.
   B3 = 2 Pi K^2 g^2 intInv/(u R). The dimensionless material integrand
   with angular averages implicit and a common prefactor of 2/intInv is
     tRad^2/g Sum[|D_rad f_h|^2] + (1/g) Sum[|D_S f_h|^2]
       - Sum[(g_h/g)|f_h|^2] + tRad^2/2 (Sum[|f_h|^2])^2.
   D_S uses the unit sphere. *)
interval[expr_] := Integrate[expr, {tRad, 1, lambda},
                             Assumptions -> lambda > 1];
width = interval[1];
volume = interval[tRad^2];
intInv = interval[1/tRad^2];
radialRatio = Factor[width^2/(volume intInv)];
emit["radial {width,volume,inverse square integral,trial factor}",
     {width, volume, intInv, radialRatio}];

(* The monopole angular energy j(j+1)-q^2 has q=|n-C_h|/2,
   j=q+ell. Jensen applies to the total density of both components. *)
level = Expand[(degree/2 + ell) (degree/2 + ell + 1) - degree^2/4];
gap = level /. ell -> 0;
localMinimum = Factor[(tRad^2 m^2/2 - gain m) /.
                      m -> (m /. First[Solve[tRad^2 m == gain, m]])];
emit["sphere {angular bottom,level increment}", {gap, Factor[level - gap]}];
emit["radially unrestricted local bound for gain>=0",
     {localMinimum, Factor[2 interval[localMinimum]/intInv]}];
native[h_, s_] := -2 h s;
degreeAt[n_, h_, s_] := Abs[n - native[h, s]];
g = 3/2;
ratio[h_, strong_, weakRatio_] := If[h == strong, 1, weakRatio];
sectorLower[n_, s_, strong_, weakRatio_, y_] := Module[{gainMax},
  gainMax = Max[0, Table[ratio[h, strong, weakRatio] - degreeAt[n, h, s]/(2 g),
                         {h, {-1, 1}}]];
  y n^2 - gainMax^2
];

(* A north-regular degree-k section. In the southern patch its phase
   changes by k phi. Near the south pole the amplitude vanishes as
   (distance)^k for the two profiles used here. *)
avg[expr_] := Integrate[Expand[expr], {z, 0, 1}];
moments[profile_, k_] := Module[{kinetic},
  kinetic = Together[z (1-z) D[profile, z]^2 +
                     (k (1-z))^2 profile^2/(4 z (1-z))];
  {avg[profile^2], avg[kinetic], avg[profile^4]}
];
unitMoments = moments[Sqrt[z] (1 + (1-z)/2), 1];
neutralMoments = moments[z, 2];
flatMoments = moments[1, 0];
emit["{norm,kinetic,quartic} of {unit,neutral degree two,flat} profiles",
     {unitMoments, neutralMoments, flatMoments}];
(* a is squared amplitude. Trial profiles occupy the strong channel so
   the quadratic gain coefficient is one there. The optimum is feasible
   whenever norm>kinetic/g. *)
trialAction = Expand[2/intInv (width a (grad/g - norm) + volume a^2 quartic/2)];
amplitudeSquared = Factor[a /. First[Solve[D[trialAction, a] == 0, a]]];
trialMaterial = Factor[trialAction /. a -> amplitudeSquared];
trial[n_, moments_, lambdaVal_, y_] :=
  (y n^2 + trialMaterial /. Thread[{norm, grad, quartic} -> moments]) /.
      lambda -> lambdaVal;
amplitudes = Factor[amplitudeSquared /. Thread[{norm, grad, quartic} -> #]] & /@
  {unitMoments, neutralMoments, flatMoments};
emit["stationary amplitudes squared for {unit,neutral,flat}", amplitudes];
emit["trial material energy divided by radial factor", Factor[trialMaterial/radialRatio]];
amplitudesPositive = Resolve[ForAll[lambda, lambda > 1,
                                    And @@ Thread[amplitudes > 0]], Reals];

(* Fix k_e once in physical units. Dilation of the shell changes B3 but
   cancels out of y when its aspect ratio is held fixed. *)
referenceAspect = 5/4;
targetY = 2/9;
bindingScale = 2 Pi stiffness^2 g^2 intInv/(quarticU radius);
kElectric = Factor[(targetY radius bindingScale) /. lambda -> referenceAspect];
yGeneral = Factor[(kElectric/radius)/bindingScale];
ySelfSimilar = yGeneral /. lambda -> referenceAspect;
emit["{single k_e,general y,self-similar y,H}",
     {kElectric, yGeneral, ySelfSimilar, radialRatio /. lambda -> referenceAspect}];
unitUpper = trial[1, unitMoments, referenceAspect, ySelfSimilar];
tailLower = 4 ySelfSimilar - 1;
tailCheck = Resolve[ForAll[n, Abs[n] >= 2, ySelfSimilar n^2 - 1 >= tailLower],
                    Reals];
emit["{incoherent n=0,unit upper,incoherent unit,|n|>=2 lower}",
     {0, unitUpper, ySelfSimilar, tailLower}];
emit["unit trial strictly below incoherent at every R>0", unitUpper < 0];
emit["universal |n|>=2 lower inequality", tailCheck];
comparison[s_, strong_, weak_] := Module[{aligned, zero, opposite},
  aligned = Sign[native[strong, s]];
  zero = sectorLower[0, s, strong, weak, ySelfSimilar];
  opposite = sectorLower[-aligned, s, strong, weak, ySelfSimilar];
  {weak, s, strong, aligned, unitUpper, zero, opposite, tailLower,
   Min[zero, tailLower, If[weak == 1, Infinity, opposite]] - unitUpper}
];
rows = Flatten[Table[comparison[s, strong, weak], {weak, {1, 4/5}},
                     {s, {-1, 1}}, {strong, {-1, 1}}], 2];
Scan[emit["{weak ratio,s,strong h,aligned n,trial,neutral lower,opposite lower,tail lower,gap}",
          #] &, rows];
(* At equal weights both conjugate signs have the same energy. At split
   weights fixed h across a mirror conjugates the charge sector. Swapping
   the material preference with the mirror preserves the charge sector. *)
symmetry = {
  Expand[(-n - native[h, s]) + (n - native[-h, s])],
  Expand[(-n - native[h, -s]) + (n - native[h, s])],
  Expand[(n - native[-h, -s]) - (n - native[h, s])]
};
emit["{equal-weight conjugate exchange,fixed mirror,swapped mirror} residuals", symmetry];
flatUpper = trial[0, flatMoments, referenceAspect, ySelfSimilar];
chargedLower = ySelfSimilar - 1;
emit["flat-native {neutral upper,all charged lower,gap}",
     {flatUpper, chargedLower, chargedLower - flatUpper}];

(* At fixed absolute thickness R0/4 let p=R/R0. The trial for n=0
   uses a degree-two section in the strong channel at either mirror. *)
aspectFixedThickness = 1 + 1/(4 p);
yThin = Factor[yGeneral /. lambda -> aspectFixedThickness];
thinNeutral = Factor[trial[0, neutralMoments, aspectFixedThickness, yThin]];
thinCharged = Factor[yThin - 1];
thinRange = Reduce[p > 0 && thinCharged >= 0, p, Reals];
thinCheck = Resolve[ForAll[p, p > 0, thinNeutral < 0], Reals];
emit["fixed thickness {y,neutral upper,charged lower,charged lower nonnegative}",
     {yThin, thinNeutral, thinCharged, thinRange}];
emit["fixed thickness neutral trial strictly below incoherent", thinCheck];

(* With all g_h switched off the physical material energy is kinetic
   plus a nonnegative quartic. With +mu0 Sum[|Psi|^2] the least local
   mass is mu0-K g/r^2, increasing across the collar. In either case,
   a nonnegative quadratic cap is attained by the zero section for ANY
   n; the fixed Coulomb term still costs k_e n^2/R. *)
positiveMassRange = Reduce[radius > 0 && mu0 > 0 && stiffness > 0 &&
                           mu0 - stiffness g/radius^2 >= 0, radius, Reals];
emit["positive mass sufficient radius criterion",
     {radius^2 >= stiffness g/mu0, positiveMassRange}];
emit["physical {zero-binding,positive-mass-above-threshold} sector minima",
     {kElectric n^2/radius, kElectric n^2/radius}];
checks = {amplitudesPositive, unitUpper < 0, tailCheck,
          And @@ (Positive[Last[#]] & /@ rows),
          symmetry == {0, 0, 0}, flatUpper < chargedLower, thinCheck};
If[!TrueQ[And @@ checks], Exit[1]];
