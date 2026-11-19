(* Hard-core pairs on a periodic cubic lattice; energies and total face
   twists in lattice units. Number is held fixed during differentiation.
   Directed Peierls phases below multiply final-head, initial-tail entries. *)
ClearAll["Global`*"];
emit[label_, expr_] := Print[label, " = ", ToString[expr, InputForm]];
checks = {};
test[expr_] := AppendTo[checks, TrueQ[expr]];

(* Two sites, two distinguishable hard-core bosons per site at most
   (one per face). Only a separated pair is reached by one elementary hop
   from the on-site-pair sector, at energy cost gap. *)
eye = IdentityMatrix[2];
destroy = {{0, 1}, {0, 0}};
localOp[op_, slot_] := KroneckerProduct @@
  ReplacePart[ConstantArray[eye, 4], slot -> op];
directed[tail_, head_] := Transpose[localOp[destroy, head]].localOp[destroy, tail];
faceHop[tail_, head_, h_, phase_] := -h (phase directed[tail, head] +
  Transpose[directed[tail, head]]/phase);
bare = faceHop[1, 3, ht, ut] + faceHop[2, 4, hb, ub];
pairStates = 1 + (FromDigits[#, 2] & /@
    {{0, 0, 0, 0}, {1, 1, 0, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}});
other = Complement[Range[16], pairStates];
bond2 = Map[Factor, -bare[[pairStates, other]].bare[[other, pairStates]]/gap, {2}];
transfer = -bond2[[3, 2]]/(ut ub) // Factor;
vNatural = bond2[[4, 4]] + bond2[[1, 1]] - bond2[[2, 2]] - bond2[[3, 3]] // Factor;
extra = Factor[-2 transfer - vNatural];
emit["bare order-two bond in {empty,left,right,full}", bond2];
emit["pair hop T, induced V, induced V minus 2T, extra V for ferro point",
  {transfer, vNatural, Factor[vNatural - 2 transfer], extra}];
test[Factor[transfer - 2 ht hb/gap] == 0 && Factor[extra + 2 (ht + hb)^2/gap] == 0];

(* Adding T(n_i+n_j) per undirected bond merely shifts each number sector
   by a linear energy. At V=-2T the following local projector results. *)
exchange = IdentityMatrix[4][[{1, 3, 2, 4}]];
ferroBond = {{0, 0, 0, 0}, {0, hop, -hop, 0}, {0, -hop, hop, 0}, {0, 0, 0, 0}};
emit["shifted ferro bond, divided by T; T(1-swap) discrepancy",
  {ferroBond/hop, Simplify[ferroBond - hop (IdentityMatrix[4] - exchange)]}];
test[ferroBond == hop (IdentityMatrix[4] - exchange)];

(* Fixed-number ground state is constant: each exchange transposition fixes
   it, and nearest-neighbor transpositions generate all site permutations
   on any connected lattice. An oriented bond moves the pair in exactly
   Binomial[n-2,m-1] configurations. *)
forwardProb = FullSimplify[FunctionExpand[Binomial[n - 2, m - 1]/Binomial[n, m]],
    Assumptions -> n > m > 0 && Element[{n, m}, Integers]] // Factor;
deltaEnergy = hop n forwardProb 2 (1 - Cos[(thetaT + thetaB)/ell]);
(* deltaEnergy is the trial-state energy change. Its Hessian is the true
   ground-energy Hessian: every component of H'|Dicke> is proportional
   to Sum[n_(x+1)-n_x] on complete periodic rows and hence zero. This
   equality concerns the second derivative only. *)
angles = {thetaT, thetaB};
helicity = Map[Factor, (ell^2/n) Table[
    D[deltaEnergy, a, b] /. {thetaT -> 0, thetaB -> 0}, {a, angles}, {b, angles}], {2}];
bulk = Map[Factor[Limit[# /. m -> filling n, n -> Infinity]] &, helicity, {2}];
costs[mat_] := { {1, 1}.mat.{1, 1}, {1, -1}.mat.{1, -1} } // Factor;
emit["directed bond coherence, finite helicity matrix", {forwardProb, helicity}];
emit["bulk helicity, same/opposite directional curvature, normalized eigenvalues",
  {bulk, costs[bulk], Eigenvalues[bulk] // Factor}];
test[Simplify[bulk - 2 hop filling (1 - filling) ConstantArray[1, {2, 2}]]
    == ConstantArray[0, {2, 2}]];
scaling = Table[With[{ns = side^3, ms = Floor[side^3/2]},
   {side, ms, ns, First[First[helicity]] /. {n -> ns, m -> ms, hop -> 1}}],
   {side, {3, 4, 6, 8}}];
emit["L,M,N,scaled TT curvature for floor-half cubic sequence at T=1", scaling];
emit["empty and full matrices", {helicity /. m -> 0, helicity /. m -> n} // Simplify];
separate = DiagonalMatrix[{helicity[[1, 1]] /. hop -> topHop,
                            helicity[[1, 1]] /. hop -> botHop}];
emit["two independent exchange copies with M per face, matrix and directional curvatures",
  {separate, costs[separate]}];
test[Subtract @@ costs[separate] == 0];

(* Independent finite cubic matrix construction. In seam gauge the first
   derivative moves the ground-state vector and a nonzero current correction
   subtracts from the diamagnetic curvature. Include both contributions. *)
side = 3; ns = side^3; ms = 2;
sites = Tuples[Range[0, side - 1], 3];
idx[pos_] := 1 + FromDigits[pos, side];
bonds = Flatten[Table[
    {idx[pos], idx[Mod[pos + UnitVector[3, axis], side]], axis,
     Boole[axis == 1 && pos[[1]] == side - 1]},
    {pos, sites}, {axis, 3}], 1];
states = Subsets[Range[ns], {ms}];
codes = Total[2^(# - 1)] & /@ states;
lookup = AssociationThread[codes -> Range[Length[states]]];
(* One row per possible swap, including the source row occupation sign.
   A positive sign means the row is the head-occupied final state. *)
moves = Reap[Do[
    Do[
      sgn = BitGet[codes[[i]], b[[2]] - 1] - BitGet[codes[[i]], b[[1]] - 1];
      If[sgn != 0,
        next = BitXor[codes[[i]], 2^(b[[1]] - 1) + 2^(b[[2]] - 1)];
        Sow[{i, lookup[next], sgn, Boole[b[[3]] == 1]/side, b[[4]]}]],
      {b, bonds}], {i, Length[states]}]][[2, 1]];
dim = Length[states]; ones = ConstantArray[1, dim];
matrix[values_] := SparseArray[Thread[(moves[[All, {1, 2}]]) -> values], {dim, dim}];
adjacency = matrix[ConstantArray[1, Length[moves]]];
h0 = SparseArray[Band[{1, 1}] -> adjacency.ones, {dim, dim}] - adjacency;
score = Total[sites[[#, 1]]]/side & /@ states;
finiteResults = Table[Module[{weights, deriv, second, firstState, para, dia, residual},
    weights = moves[[All, col]];
    deriv = matrix[-I moves[[All, 3]] weights];
    second = matrix[weights^2];
    firstState = If[col == 4, 0 ones, -I (score - Mean[score])];
    residual = h0.firstState + deriv.ones;
    para = 2 Re[Conjugate[firstState].(deriv.ones)]/dim;
    dia = ones.(second.ones)/dim;
    {Max[Abs[deriv.ones]], Max[Abs[residual]], dia, para, dia + para}],
    {col, {4, 5}}];
emit["cube {L,M,dimension}, max norm H0 on uniform vector",
  {{side, ms, dim}, Max[Abs[h0.ones]]}];
emit["uniform/seam {max current action, eigenvector derivative residual, diamagnetic, current correction, total twist curvature} at T=1",
  finiteResults];
finiteExpected = ns/side^2 helicity[[1, 1]] /. {n -> ns, m -> ms, hop -> 1};
emit["count prediction for total twist curvature, discrepancies",
  {finiteExpected, finiteResults[[All, 5]] - finiteExpected}];
test[Max[Abs[h0.ones]] == 0 && finiteResults[[All, 2]] == {0, 0} &&
  finiteResults[[All, 5]] == {finiteExpected, finiteExpected} && finiteResults[[2, 4]] < 0];

(* Independent pinned control at half filling. Tile an even lattice by
   isolated nearest-neighbor dimers, keep finite intra-dimer transfer and
   the same -2T attraction there, and impose opposite site energies +/-pin.
   Endpoints have E0=0, E2=-2T. A phase on an isolated dimer is pure gauge. *)
dimer = {{pin, -hop/u}, {-hop u, -pin}};
poly = Factor[Det[en eye - dimer]];
dimerGround = -Sqrt[pin^2 + hop^2];
pinGap = Simplify[-2 hop - 2 dimerGround];
emit["pinned dimer one-pair polynomial, E0,E1,E2, particle-hole gap, twist curvature",
  {poly, {0, dimerGround, -2 hop}, pinGap, D[dimerGround, thetaT, thetaB]}];
test[Simplify[poly /. en -> dimerGround] == 0 && FreeQ[poly, u] &&
     FullSimplify[pinGap > 0, hop > 0 && pin > 0]];

(* At V=-2T+delta, only the delta*n_i*n_j term changes the shifted
   Hamiltonian. Endpoints have shifted energy densities 0 and 3 delta;
   a compact filled domain gives 3 delta*filling in the infinite cube.
   For delta<0 this attains the lower bound from Sum n_i n_j <= 3 M.
   The uniform state is a trial off the point, not an exact eigenvector. *)
joint = Factor[forwardProb (m - 1)/(n - m)];
uniformSlope = Factor[Limit[3 joint /. m -> filling n, n -> Infinity]];
coexistenceSlope = 3 filling;
emit["Dicke nearest-neighbor joint occupancy, homogeneous and separated trial d(e/site)/d(delta), their difference",
  {joint, uniformSlope, coexistenceSlope, Factor[uniformSlope - coexistenceSlope]}];
emit["unshifted tuned E0(M), inverse number susceptibility numerator",
  {-6 hop m, Expand[-6 hop (m + 1) - 6 hop (m - 1) + 12 hop m]}];
(* A plane wave in the one-pair sector; global spin lowering commutes
   with exchange and carries nonzero-momentum magnons into interior sectors. *)
directions = Join[IdentityMatrix[3], -IdentityMatrix[3]];
disp = hop Total[(1 - Exp[I {kx, ky, kz}.#]) & /@ directions];
magnon = ComplexExpand[disp] // Simplify;
emit["exchange one-pair dispersion, quadratic term along one axis",
  {magnon, SeriesCoefficient[magnon /. {ky -> 0, kz -> 0}, {kx, 0, 2}]}];
emit["checks", checks];
If[!TrueQ[And @@ checks], Exit[1]];
