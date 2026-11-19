(* What spin-blind hopping reads after the collar preparation. A single
   incoming spin-up carrier has a phase-registered path on either arm.
   The energy is a held-state response to a bias against the preparation.
   h is the positive survival amplitude of the excited component. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
checks = {};
check[value_] := AppendTo[checks, TrueQ[value]];
up = {1, 0};
id = IdentityMatrix[2];
sigma = Table[PauliMatrix[j], {j, 3}];
project[axis_] := (id + axis.sigma)/2;
outer[v_] := Outer[Times, v, Conjugate[v]];
inner[v_, w_] := Conjugate[v].w;

(* Each row is a spin amplitude for an orthogonal receiver record. For a
   fixed axis, after the first emission the ground state cannot emit again.
   For repeated encounters the vacuum map composes exactly. Excitations at
   different encounter numbers are orthogonal and have summed probability
   1-h^2. The effective excited receiver triplet below has that total norm.
   The common phase of the emitted amplitudes is unobservable here. *)
records[axis_, h_] := Module[{g = project[axis], e = id-project[axis]},
  Prepend[Table[Sqrt[(1-h^2)/2] g.sigma[[j]].e.up, {j, 3}], (g + h e).up]
];
siteState[amplitudes_] := Total[outer /@ amplitudes];
zeros = ConstantArray[0, {3, 2}];
leftRecords[amplitudes_] := Flatten[Join[amplitudes, zeros]];
rightRecords[amplitudes_, overlap_] := Flatten[Join[{First[amplitudes]},
  overlap Rest[amplitudes], Sqrt[1-overlap^2] Rest[amplitudes]]];
spinMarginal[joint_] := siteState[Partition[joint, 2]];

(* At the opposed pole both arms have the same local reduced density. Only
   the final receiver record distinguishes common and tagged preparations. *)
held = records[{0, 0, -1}, h];
left = leftRecords[held];
realH[value_] := FullSimplify[value, Assumptions -> 0 < h <= 1 && 0 <= u <= 1];
heldMarginals = realH[spinMarginal /@ {left, rightRecords[held,u]}];
heldOverlap = realH[inner[left, rightRecords[held,u]]];
crossingMarginals = realH[heldMarginals /. h -> 1/Sqrt[2]];
crossingFringes = (heldOverlap /. h -> 1/Sqrt[2]) /. {{u -> 0}, {u -> 1}};
emit["held marginals, record overlap u, spin-blind cross-arm amplitude",
  {heldMarginals, u, heldOverlap}];
emit["at zero polarization: marginals, fringe amplitudes {tagged,common}",
  {crossingMarginals, crossingFringes}];
check[heldMarginals[[1]] === heldMarginals[[2]] &&
  crossingMarginals === {id/2, id/2} && crossingFringes === {1/2,1}];

(* A numerical readout needs a path-phase reference. A spin-independent
   hopping energy at fixed state gives E/t= -Re[e^(-i a) <left|right>].
   Zero and quarter-turn biases determine both quadratures. Rephasing an
   arm rephases the corresponding link. *)
axes = IdentityMatrix[3];
finite = records[#, 1/2] & /@ axes;
vacuum = First /@ finite;
link[vectors_] := MapThread[inner, {vectors, RotateLeft[vectors]}];
taggedLinks = MapThread[inner, {leftRecords /@ finite,
  (rightRecords[#,0] & /@ RotateLeft[finite])}];
localPhases = {1, I, -1};
changed = MapThread[Times, {localPhases, vacuum}];
phaseFactors = link[localPhases /. z_?NumericQ :> {z}];
(* An explicit two-site hopping matrix acts on paths times records and spin. *)
hopping[bias_] := -KroneckerProduct[{{0, Exp[-I bias]}, {Exp[I bias], 0}},
    IdentityMatrix[Length[left]]];
joint = Join[leftRecords[finite[[1]]],rightRecords[finite[[2]],0]]/Sqrt[2];
energy = FullSimplify[inner[joint, hopping[a].joint], Assumptions->Element[a,Reals]];
quadratures = FullSimplify[-energy /. a -> #] & /@ {0,Pi/2};
emit["tagged triangle: arm overlaps, vacuum prediction, rephased loop/control",
  {taggedLinks, link[vacuum], Times@@link[changed], Times@@taggedLinks}];
emit["held energy / hopping and two bias quadratures", {energy,quadratures}];
check[taggedLinks === link[vacuum] && link[changed] === phaseFactors taggedLinks &&
  Times@@link[changed] === Times@@taggedLinks];
check[quadratures === ReIm[First[taggedLinks]]];

(* The vacuum amplitudes are never normalized branch probabilities in the
   energy. Normalize only to ask about their preferred connection. With
   x=cos(theta), azimuth f, differentiate its actual phase connection. *)
normal = {Sqrt[1-x^2] Cos[f], Sqrt[1-x^2] Sin[f], x};
realSphere[value_] := FullSimplify[value, Assumptions ->
    -1<x<1 && 0<h<=1 && Element[f,Reals]];
vac = realSphere[First[records[normal,h]]];
norm2 = realSphere[inner[vac,vac]];
normalized = vac/Sqrt[norm2];
aPhi = realSphere[-I inner[normalized,D[normalized,f]]];
poles = FullSimplify[{aPhi /. x -> 1, aPhi /. x -> -1}, Assumptions->0<h<=1];
emit["vacuum spinor, squared norm, azimuth connection, {north,south}, cut flux/(2 Pi)",
  {vac,norm2,aPhi,poles,Last[poles]-First[poles]}];
check[realSphere[norm2>=h^2 && Re[First[vac]]>=h] && poles === {0,0}];
(* A circulating spin whose transmission amplitude at each contact is the
   local density matrix provides a different, closed-loop probe. Each
   filter is a contraction for 0<=s<=1. An unpolarized beam and an
   unfiltered reference give Tr[F1 F2 F3]/2 for this traversal order. *)
filter[n_] := (id+s n.sigma)/2;
octants = (If[Det[#]>0,#,Reverse[#]] &) /@ (DiagonalMatrix /@ Tuples[{-1,1},3]);
loop[vertices_] := Expand[Tr[Dot@@(filter /@ vertices)]/2];
loopSet = Union[Factor /@ (loop /@ octants)];
degrees = Table[FullSimplify[Total[Arg /@ ((loop /@ octants) /. s -> p)]/(2 Pi)],
  {p,{0,1/2,1}}];
emit["filtered orthogonal octants: distinct fringes, summed preferred flux/(2 Pi) for s=0,1/2,1",
  {loopSet,degrees,N[degrees,9]}];
check[First[degrees]===0 && Last[degrees]===1 && 0<N[degrees[[2]]]<1];
emit["checks",checks];
If[!TrueQ[And@@checks],Exit[1]];
