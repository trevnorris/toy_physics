(* Unit cubic lattice. There is one directed edge registration per positive
   direction per cell, shared with all incident squares. Each square owns
   four corner angles. K = stiffness, jn = corner inertia, ie = edge
   inertia, s = angular frequency squared. The face exchange reverses all
   internal angles. Phasors are exp(i k.n-i omega t). *)
ClearAll["Global`*"];
emit[name_, data_] := Print[name, " = ", ToString[data, InputForm]];
checks = {};
zeroCheck[name_, expr_] := Module[{res},
  res = Map[Factor, Flatten[{expr}]];
  emit[name, DeleteDuplicates[res]];
  AppendTo[checks, AllTrue[res, # === 0 &]]
];

id4 = IdentityMatrix[4];
cycle = RotateLeft[id4] - id4;
slip = Array[b, 4];
minAngles = PseudoInverse[cycle].slip;
minStrain = Expand[cycle.minAngles - slip];
strainProjector = id4 - cycle.PseudoInverse[cycle];
uMin = Factor[stiff minStrain.minStrain/2];
emit["square minimum {mean-zero corners, strains, energy}", {minAngles, minStrain, uMin}];
zeroCheck["corner static equilibrium and free mean", Join[Transpose[cycle].minStrain, {Total[minAngles]}]];

(* The rows introduce the same shared edge into consecutive sides. The
   ordering is (0,mu,mu+nu,nu) with a positive oriented boundary. *)
basis = IdentityMatrix[3];
pairs = {{1, 2}, {2, 3}, {3, 1}};
squareMap[{mu_, nu_}, phase_] := {
  basis[[mu]], phase[[mu]] basis[[nu]], -phase[[nu]] basis[[mu]], -basis[[nu]]};
maps[phase_] := squareMap[#,phase]& /@ pairs;
zeroPhases = ConstantArray[1,3];
wavePhases = {zx, zy, zz};
boundaries = Total /@ maps[wavePhases];
grad = wavePhases - zeroPhases;
emit["loop sums on three Bloch edges", boundaries];
zeroCheck["exact gradient circulation", boundaries.grad];
emit["edge rank at zone center and at phases {i,1,1}",
  MatrixRank[boundaries /. Thread[wavePhases -> #]]& /@ {zeroPhases,{I,1,1}}];
(* Give the full physical strain map, retaining every corner mean. *)
strainMap[phase_] := Join[-(Join @@ maps[phase]), KroneckerProduct[IdentityMatrix[3],cycle],2];
fullRanks = MatrixRank[strainMap[#]]& /@ {zeroPhases,{I,1,1}};
emit["full {strain rank, mechanical nullity} at center and {i,1,1}", {#,15-#}& /@ fullRanks];
(* Physical angular momenta: changing a by a vertex gradient while changing
   each adjacent mark by that vertex angle preserves every spring strain.
   The corresponding momentum is conserved for unconstrained initial data. *)
vertexLift[{mu_,nu_},phase_] := {1,phase[[mu]],phase[[mu]] phase[[nu]],phase[[nu]]};
liftedGradient = Join[grad, Flatten[vertexLift[#,wavePhases]& /@ pairs]];
zeroCheck["vertex transformation on full strains", strainMap[wavePhases].liftedGradient];

(* Eliminate the corner equation at s unequal to its bare poles. The mean
   is decoupled exactly. Physical means at s=0 are counted above. *)
cornerResponse = Map[Factor,
   Inverse[stiff Transpose[cycle].cycle-jn s id4].(stiff Transpose[cycle]),{2}];
response = Map[Factor,id4-cycle.cornerResponse,{2}];
modeVectors = {{1,1,1,1},{1,0,-1,0},{0,1,0,-1},{1,-1,1,-1}};
modeValues = Factor[#.response.#/(#.#)]& /@ modeVectors;
emit["signed spring-response factors for cycle harmonics {0,1a,1b,2}",modeValues];
zeroCheck["cycle harmonic eigenvector residuals",
  MapThread[response.#1-#2 #1&,{modeVectors,modeValues}]];
effective[phase_] := Map[Factor, stiff Sum[
  Transpose[squareMap[pairs[[n]],1/phase]].response.squareMap[pairs[[n]],phase],
  {n,3}]-ie s IdentityMatrix[3], {2}];
axisKernel = effective[{z,1,1}];
diag = Diagonal[axisKernel];
zeroCheck["axis off diagonals and transverse splitting", {axisKernel-DiagonalMatrix[diag],diag[[2]]-diag[[3]]}];
(* Read the axis symbols using h=(2-z-1/z)/4. The product of the
   rational denominators at a pole is not used to classify modes there. *)
asH[expr_] := Factor[(1-h) (expr /. z->1)+h (expr /. z->-1)];
axisL = asH[diag[[1]]];
axisT = asH[diag[[2]]];
zeroCheck["axis factors after h substitution",
  ({axisL,axisT} /. h->(2-z-1/z)/4) - diag[[1;;2]]];
emit["axis longitudinal/transverse inverse response, h=sin(k/2)^2",{axisL,axisT}];
axisSlope = Factor[-D[axisT,h]/D[axisT,s] /. {s->0,h->0}];
emit["transverse acoustic {d(omega^2)/dh, speed squared}",{axisSlope,axisSlope/4}];
emit["longitudinal inverse response at zero frequency",axisL /. s->0];

(* The leading derivative of the full static form in any direction is
   obtained from the three face circulations. Static response already fixes
   the principal form and only the k=0 kinetic coefficient enters at order k^2. *)
qVec = {qx,qy,qz};
curlLeading = Coefficient[boundaries /. Thread[wavePhases->zeroPhases+t qVec],t,1];
inertiaZero = Map[Factor,-D[effective[zeroPhases],s] /. s->0,{2}];
emit["three-dimensional leading {inertia,curl symbol without i}",{inertiaZero,curlLeading}];
zeroCheck["principal spatial form residual",
  Transpose[curlLeading].curlLeading-(qVec.qVec IdentityMatrix[3]-Outer[Times,qVec,qVec])];
(* Compute all zone-center modes from the physical 15-angle Hessian, even
   at the poles of the preceding Schur complement. *)
atZero = strainMap[zeroPhases];
mass = DiagonalMatrix[Join[ConstantArray[ie,3],ConstantArray[jn,12]]];
emit["full zone-center characteristic polynomial det(M^-1 H - s)",
  Factor[Det[Inverse[mass].(stiff Transpose[atZero].atZero)-s IdentityMatrix[15]]]];

(* Energy conservation on the physical spring/rotor graph. One square
   receives power from each of its four edges. Each transfers the opposite
   power into its own kinetic energy. Edge speeds below are signed. *)
cornerV = Array[v,4]; slipV = Array[w,4]; ext = Array[x,4];
torque = -stiff Transpose[cycle].ext;
squarePower = Expand[cornerV.torque+stiff ext.(cycle.cornerV-slipV)];
emit["corner torques and square power",{torque,squarePower}];
zeroCheck["square/edge power cancellation and corner net torque",
  {squarePower+stiff ext.slipV,Total[torque]}];
(* On an x wave with y-edge amplitude A, use the exact spring torque on
   the right y-edge of its xy square. The energy per cell uses full corner
   kinetics. For an eigenmode the average kinetic and spring energies agree. *)
ampY = {0,1,0};
xySlip = squareMap[{1,2},{z,1,1}].ampY;
rightStrain = Expand[-(response.xySlip)[[2]]];
rightVelocityStar = I omega/z;
powerPhasor = stiff rightStrain rightVelocityStar amp^2/2;
(* All parameters and amp are real. The velocity conjugate supplied the
   only factor of i, so conjugating sends i -> -i and z -> 1/z. *)
powerZ = Factor[(powerPhasor-(powerPhasor /. z->1/z))/2];
powerSin = Factor[powerZ/( (z-1/z)/(2 I) )];
emit["right edge time-average energy current divided by sin(k)",powerSin];
entrained = Sum[
  (cornerResponse.squareMap[pairs[[n]],{1/z,1,1}].ampY).
    (cornerResponse.squareMap[pairs[[n]],{z,1,1}].ampY),{n,3}];
modalMass = asH[Factor[ie+jn entrained]];
zeroCheck["edge power equals group velocity times full energy",
  {modalMass+D[axisT,s],powerSin-omega amp^2 D[axisT,h]/8}];
emit["long-wave right-edge current coefficient divided by omega sin(k) A^2",
  Factor[powerSin/(omega amp^2) /. s->0]];

(* Controls. Uniform pinning costs anchor*phi_j^2/2 at every corner. A
   single mean constraint per square, by comparison, leaves all strains
   attainable. Cutting the fourth spring gives a three-link path. *)
pinned = id4-cycle.Inverse[stiff Transpose[cycle].cycle+anchor id4].(stiff Transpose[cycle]);
uniformStiffness = Factor[stiff Sum[
  With[{m = squareMap[pairs[[n]],zeroPhases]},(m.ampY).pinned.(m.ampY)],{n,3}]];
path = cycle[[1;;3]];
pathResidual = IdentityMatrix[3]-path.PseudoInverse[path];
emit["all-corner anchoring {uniform static stiffness, rigid-anchor limit}",
 {uniformStiffness,Limit[uniformStiffness,anchor->Infinity]}];
emit["cut square {path rank,remaining spring residual projector}",{MatrixRank[path],pathResidual}];
zeroCheck["pinning only the first corner changes no static strains",
  cycle.(minAngles-ConstantArray[minAngles[[1]],4])-cycle.minAngles];
emit["all identity checks",And @@ checks];
Exit[If[And @@ checks,0,1]];
