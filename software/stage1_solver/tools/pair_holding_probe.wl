(* Axisymmetric harmonics in unit-radius sphere frames. Evaluate each incident
   multipole by translation to regular harmonics at the other center. All
   boundary matching and angular force integration use exact arithmetic. *)
ClearAll["Global`*"];
SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
translation[s_, l_, dx_] :=
  (-1)^s Sign[dx]^(s+l) Binomial[s+l,l]/Abs[dx]^(s+l+1);

measure[centers_, degree_] := Module[
  {n = Length[centers], labels, shift, op, rhs, coeff, regular, inertia,
   flux, incident, angular, indexes, local, incoming, trace, vx, speed2,
   p, mu, integrate, mobility, load, hold, free, probes, observed, predicted,
   e, checks, results = {}, probe = 3/200, i, j, k, minus, plus},
  labels = Flatten[Table[{i,l},{i,n},{l,0,degree}],1];
  shift = Table[If[row[[1]] == col[[1]], 0,
      translation[row[[2]],col[[2]],centers[[row[[1]]]]-centers[[col[[1]]]]]],
    {row,labels},{col,labels}];
  op = DiagonalMatrix[-(Last /@ labels)-1] + DiagonalMatrix[Last /@ labels].shift;
  rhs = Table[If[col==0,-Boole[row[[2]]==0],Boole[row=={col,1}]],
    {row,labels},{col,0,n}];
  coeff = LinearSolve[op,rhs];
  regular = shift.coeff;
  angular = Table[LegendreP[l,mu],{l,0,degree}];
  integrate[expr_] := With[{c = CoefficientList[Expand[expr],mu]},
    Sum[2 c[[j]]/j,{j,1,Length[c],2}]];
  flux = incident = ConstantArray[0,n];
  inertia = IdentityMatrix[n];
  Do[
    indexes = Range[(i-1)(degree+1)+1,i(degree+1)];
    local = coeff[[indexes,All]];
    incoming = regular[[indexes,All]];
    trace = angular.(local[[All,1]]+incoming[[All,1]]);
    (* Exact boundary normal speed -1 and tangential velocity from the trace. *)
    vx = -mu+(1-mu^2) D[trace,mu];
    speed2 = 1+(1-mu^2) D[trace,mu]^2;
    flux[[i]] = (3/2) integrate[vx+mu speed2/2];
    incident[[i]] = incoming[[2,1]];
    (* Projection onto the sole acceleration harmonic carrying net force. *)
    inertia[[i]] -= (local+incoming)[[2,2;;]];
    ,{i,n}];
  mobility = Inverse[inertia];
  Do[
    load = flux-3 k incident;
    hold = -load;
    free = LinearSolve[inertia,load];
    probes = Table[e = UnitVector[n,j];
      minus = LinearSolve[inertia,load-probe e];
      plus = LinearSolve[inertia,load+probe e];
      {minus,plus}, {j,n}];
    observed = Transpose[(#[[2]]-#[[1]])/(2 probe)& /@ probes];
    predicted = observed.(-hold);
    checks = {observed===mobility, predicted===free, SymmetricMatrixQ[inertia],
      If[n==1,flux==={0} && inertia==={{3/2}}, Total[hold]===0 && Total[free]===0]};
    Print["k=",k," {H_target,a_minus,a_plus,a_free}=",
      ToString[N[{hold[[1]],probes[[1,1,1]],probes[[1,2,1]],free[[1]]},16],InputForm],
      " controls=",ToString[checks,InputForm]];
    If[!And@@checks,Exit[1]];
    AppendTo[results,{hold[[1]],probes[[1,1,1]],probes[[1,2,1]],free[[1]]}],
    {k,{0,1}}];
  Print["mobility_target_row=",ToString[N[mobility[[1]],16],InputForm]];
  results
];

Do[
  Print["centers=",ToString[centers,InputForm]," degree=5"]; low=measure[centers,5];
  Print["centers=",ToString[centers,InputForm]," degree=9"]; high=measure[centers,9];
  change=Max[Abs[Flatten[high-low]]];
  Print["max_refinement_change=",ToString[N[change,6],InputForm]];
  If[change>2*10^-10,Exit[1]],
  {centers,{{0},{0,10},{0,20}}}];
