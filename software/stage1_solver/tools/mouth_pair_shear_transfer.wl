(* Flat resting slab, linear calibrated drives and readouts, fixed nonzero
   three-dimensional wave vector. Reflection in the normal w exchanges the
   two faces without reversing the in-brane laboratory curl. Work in one
   transverse helicity with curl eigenvalue q = +k or -k. The projector
   argument permits polarization mixing within each reflection sector. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label <> " = " <> ToString[value, InputForm]];
checks = {};
test[condition_] := AppendTo[checks, TrueQ[condition]];

(* Generic retarded inverse response. Entries can be complex functions of
   frequency and wave vector; no reciprocal or instantaneous law assumed.
   Face coordinates are {u_t,u_b,A_t,A_b}. The parity coordinates are
   {u_even,A_even,u_odd,A_odd}. *)
reflect = IdentityMatrix[4][[{2, 1, 4, 3}]];
parity = {{1, 1, 0, 0}, {0, 0, 1, 1}, {1, -1, 0, 0}, {0, 0, 1, -1}}/Sqrt[2];
raw = Array[entry, {4, 4}];
invariant = (raw + reflect.raw.reflect)/2;
reduced = Expand[parity.invariant.Transpose[parity]];
offBlock = {reduced[[1 ;; 2, 3 ;; 4]], reduced[[3 ;; 4, 1 ;; 2]]};
emit["reflection eigenvalues, generic parity-projected off-diagonal blocks",
  {Diagonal[parity.reflect.Transpose[parity]], offBlock}];
test[offBlock == ConstantArray[0, {2, 2, 2}]];
blocks = {{me, fe, 0, 0}, {ge, ae, 0, 0}, {0, 0, mo, fo}, {0, 0, go, ao}};
de = Det[blocks[[1 ;; 2, 1 ;; 2]]];
do = Det[blocks[[3 ;; 4, 3 ;; 4]]];
response = Inverse[blocks];
flux = {0, 0, 0, q};
commonDrive = {1, 0, 0, 0}; differentialDrive = {0, 0, 1, 0};
faceDrives = (parity.#) & /@ IdentityMatrix[4][[1 ;; 2]];
channels = Simplify[flux.response.#] & /@ Join[{commonDrive, differentialDrive}, faceDrives];
emit["retarded determinant, even and odd determinants and factorization residual",
  {Factor[Det[blocks]], {de, do}, Expand[Det[blocks] - de do]}];
emit["odd-flux transfer from {common,differential,top,bottom} mechanical force",
  channels];
emit["{common mechanical,odd-flux} auto-responses",
  {commonDrive.response.commonDrive, flux.response.flux}];
test[channels[[1]] == 0 && Simplify[channels[[3]] - channels[[2]]/Sqrt[2]] == 0 &&
  Simplify[channels[[4]] + channels[[2]]/Sqrt[2]] == 0];

(* A two-coordinate degeneracy can produce identical frequencies in the
   separate auto-spectra. Their joint weight at that frequency has rank 2.
   z=(omega+i0)^2; a pole weight below multiplies 1/(lambda-z). *)
pick = DiagonalMatrix[{1, q}];
degenerate = pick.Inverse[(c k^2 - z) IdentityMatrix[2]].Transpose[pick];
degWeight = Simplify[(c k^2 - z) degenerate];
emit["equal-speed independent even shear and odd signed-field pole weight, determinant",
  {degWeight, Det[degWeight]}];
test[degWeight[[1, 2]] == 0 && Det[degWeight] == q^2];

(* Break exchange symmetry through a calibrated coefficient eta multiplying
   the two spatial gradients. Unit inertias, positive elastic matrix:
   a,b>0 and eta^2<a b. Both unconstrained coordinates remain mobile. *)
elastic = {{a, eta}, {eta, b}};
mixed = k^2 elastic - z IdentityMatrix[2];
sourceToReading = Simplify[pick.Inverse[mixed].Transpose[pick]];
disc = Sqrt[(a - b)^2 + 4 eta^2];
vel = {(a + b - disc)/2, (a + b + disc)/2};
projectors = {(vel[[2]] IdentityMatrix[2] - elastic)/disc,
              (elastic - vel[[1]] IdentityMatrix[2])/disc};
weights = Simplify[pick.#.Transpose[pick]] & /@ projectors;
reconstruction = Simplify[Total[MapThread[#1/(k^2 #2 - z) &, {weights, vel}]] -
                             sourceToReading];
emit["unequal-face elastic matrix, characteristic polynomial in squared speed, squared speeds",
  {elastic, Factor[Det[elastic - s IdentityMatrix[2]]], vel}];
emit["mixed odd-flux/common-shear transfer; pole weights from slow, fast poles",
  {sourceToReading[[2, 1]], weights[[All, 2, 1]]}];
emit["determinants of the two pole-weight matrices, full response reconstruction residual",
  {Simplify[Det /@ weights], reconstruction}];
equalBare = FullSimplify[{vel, weights[[All, 2, 1]]} /. b -> a, eta > 0];
emit["equal bare speeds with eta>0, squared speeds and cross weights", equalBare];
test[Simplify[Det /@ weights] == {0, 0} && reconstruction == ConstantArray[0, {2, 2}]];

(* Constitutive alternative: only one gapless odd material coordinate W in
   the observed odd sector, and both an actual differential mechanical
   coordinate and the measured odd flux are postulated to be proportional
   to W (the latter through a curl). Stiffness B and inertia rho positive.
   Calibrate g through differential face displacement under a known drive.
   Ordinary common drive has no coupling to W in the symmetric slab. *)
readout = {0, g, q};
oneWeight = Outer[Times, readout, readout]/rho;
oneResponse = Outer[Times, readout, readout]/(stiff k^2 - rho z);
emit["postulated odd mechanical and flux law: squared speed, pole weight for {common,differential,flux}, rank",
  {stiff/rho, oneWeight, MatrixRank[oneWeight]}];
emit["odd-field response to differential drive, common drive",
  oneResponse[[3, {2, 1}]]];
test[MatrixRank[oneWeight] == 1 && oneResponse[[3, 1]] == 0];
emit["checks", checks];
If[!TrueQ[And @@ checks], Exit[1]];
