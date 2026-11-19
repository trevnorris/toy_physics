(* One parcel with linear right-moving dispersion vp traverses a static
   interval of length v, in units hbar = transit time = gap = hopping = 1.
   The middle Hamiltonian maintains the lead gap and adds the rim couplings.
   Positive torque denotes the force exerted by the parcel on F. *)
ClearAll["Global`*"];
emit[label_, values_] := Print[label, " = ", ToString[N[values, 13], InputForm]];
residuals = {};
record[res_] := AppendTo[residuals, Max[Abs[Flatten[{res}]]]];
id2 = IdentityMatrix[2];
sigma = Table[PauliMatrix[j], {j, 3}];
lead = KroneckerProduct[IdentityMatrix[3], (id2 - sigma[[3]])/2];
up = (IdentityMatrix[6] - lead)/3;
mixed = IdentityMatrix[6]/6;
right = IdentityMatrix[3];
left = Reverse[right];
achiral = right[[{1, 3, 1}]];

interaction[pins_, flux_] := Module[{hop, local},
  hop = {{0, -1, -Exp[I flux]}, {-1, 0, -1}, {-Exp[-I flux], -1, 0}};
  local = Sum[KroneckerProduct[DiagonalMatrix[UnitVector[3,j]],
      -2 pins[[j]].sigma], {j, 3}];
  local + KroneckerProduct[hop, id2]];

(* Independently integrate the edge force by advancing the density at 256
   equal time intervals. This uses the physical current operator with the
   sign fixed by the Peierls factor on 20 and its adjoint. *)
edgeForce[f_] := KroneckerProduct[
    {{0, 0, I Exp[I f]}, {0, 0, 0}, {-I Exp[-I f], 0, 0}}, id2];
mean[rho_, operator_] := Re[Tr[rho.operator]];
integrateForce[h_, force_, initial_] := Module[
  {n = 256, state = initial, step, adjoint, values},
  step = MatrixExp[-I h/n];
  adjoint = ConjugateTranspose[step];
  values = Prepend[Table[state = step.state.adjoint; mean[state,force], {n}],
     mean[initial,force]];
  (First[values] + Last[values] + 4 Total[values[[2 ;; n ;; 2]]] +
    2 Total[values[[3 ;; n-1 ;; 2]]])/(3 n)];

(* The off-diagonal block of the propagator is its exact first derivative
   in F. Integrating U(t)^dagger(-d_F h)U(t) gives -i U^dagger d_F U.
   At each interface -delta<h> is the change in drift energy vp. *)
read[pins_, f_, initial_, active_:True] := Module[
  {v, h, dq, force, zero, augmented, u, du, out, de, tau, dE, tauInt,
   entry, leave, driftIn, driftOut, driftDuring},
  v = If[active, interaction[pins,theta], ConstantArray[0,{6,6}]];
  dq = N[D[v,theta] /. theta -> f, 40];
  v = N[v /. theta -> f, 40];
  h = N[lead + v, 40];
  zero = ConstantArray[0,{6,6}];
  augmented = MatrixExp[-I ArrayFlatten[{{h,dq},{zero,h}}]];
  u = augmented[[1;;6,1;;6]];
  du = augmented[[1;;6,7;;12]];
  out = u.initial.ConjugateTranspose[u];
  de = mean[out - initial,lead];
  tau = mean[initial,-I ConjugateTranspose[u].du];
  dE = 2 Re[Tr[lead.du.initial.ConjugateTranspose[u]]];
  force = If[active, N[edgeForce[f],40], zero];
  tauInt = integrateForce[h,force,initial];
  entry = mean[initial,v]; leave = -mean[out,v];
  driftIn = -entry; driftOut = -leave;
  driftDuring = mean[initial-out,h];
  record[{de-entry-leave, de+driftIn+driftOut, driftDuring, tauInt-tau,
          Norm[ConjugateTranspose[u].u-IdentityMatrix[6]],
          Tr[out]-1, Norm[force+dq]}];
  {de, tau, dE, entry, leave, driftIn, driftDuring, driftOut, tauInt-tau}];

emit["columns {DeltaEint,tau_F,dEint/dF,Eentry,Eexit,Kentry,Kzone,Kexit,integrated-minus-derivative torque}", {}];
Do[
  plus[f] = read[right,f,up];
  minus[f] = read[left,f,up];
  emit["xyz F="<>ToString[f], plus[f]];
  emit["zyx F="<>ToString[f], minus[f]],
  {f,{0,2/25,-2/25,1/25,-1/25}}];
Do[
  record[plus[f][[{1,2,3}]] - {1,-1,-1} minus[-f][[{1,2,3}]]];
  emit["xyz half-differences {DeltaEint,tau_F} and DeltaEint secant slope at |F|=" <>
       ToString[f], Join[(plus[f]-plus[-f])[[{1,2}]]/2,
       {(plus[f][[1]]-plus[-f][[1]])/(2 f)}]],
  {f,{2/25,1/25}}];
controlP = read[achiral,2/25,up];
controlM = read[achiral,-2/25,up];
unpolarized = read[right,2/25,mixed];
off = read[right,2/25,up,False];
emit["xzx F=2/25 and -2/25 {DeltaEint,tau_F}", {controlP[[1;;2]],controlM[[1;;2]]}];
emit["mixed xyz and zero interaction F=2/25 {DeltaEint,tau_F}", {unpolarized[[1;;2]],off[[1;;2]]}];
record[{controlP[[1]]-controlM[[1]], controlP[[2]]+controlM[[2]],
   unpolarized[[1;;3]], off[[1;;3]]}];
emit["maximum absolute residual across all rows and controls", Max[residuals]];
Exit[If[TrueQ[Max[residuals] < 10^-7],0,1]];
