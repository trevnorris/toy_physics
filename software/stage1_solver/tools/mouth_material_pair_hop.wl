(* Left/right patches, one mobile component on each face; the faces are
   distinguishable. Hard-core bottom occupancy in the one-vacancy control.
   Paired sites have energy zero and separated sites energy gap > 0.
   ht,hb are real positive bare hopping energies, ut,ub unit link phases.
   Rows are final states. Matrix entry (R,L) is the directed L->R hop. *)
ClearAll["Global`*"];
emit[name_, value_] := Print[name <> " = " <> ToString[value, InputForm]];
hop[t_, u_] := {{0, -t/u}, {-t u, 0}};
eye = IdentityMatrix[2];
swap = {{0, 1}, {1, 0}};
bottom["particle"] = hop[hb, ub];
(* A left vacancy means bottom material at right, and vice versa.
   Use the original material hopping in the complementary basis. *)
bottom["hole"] = swap . hop[hb, ub] . swap;
paired = {1, 4}; separated = {2, 3}; (* LL,LR,RL,RR *)
bare[kind_] := KroneckerProduct[hop[ht, ut], eye] +
  KroneckerProduct[eye, bottom[kind]] + DiagonalMatrix[{0, gap, gap, 0}];
schur[h_, energy_] := h[[paired, paired]] + h[[paired, separated]] .
  Inverse[energy eye - h[[separated, separated]]] . h[[separated, paired]];
(* At energy zero the Schur complement gives the order-hop^2 / gap
   Hamiltonian. Its energy-dependent version below is exact for this dimer. *)
transfer[h_] := Factor[schur[h, 0][[2, 1]]];
link[kind_, top_, bot_] := Factor[transfer[bare[kind]] /
  (transfer[bare[kind]] /. {ut -> 1, ub -> 1}) /. {ut -> top, ub -> bot}];
charTarget = (e (e - gap) - (ht + hb)^2) (e (e - gap) - (ht - hb)^2);
Do[
  h = bare[kind];
  paths = Table[Factor[-h[[4, mid]] h[[mid, 1]] / gap], {mid, separated}];
  exactResidual = Factor[Det[e IdentityMatrix[4] - h] - charTarget];
  emit[kind <> " bottom L->R bare hop", bottom[kind][[2, 1]]];
  emit[kind <> " forward via LR,RL", paths];
  emit[kind <> " pair diagonal, forward, reverse to second order",
       Factor /@ {schur[h, 0][[1, 1]], transfer[h], schur[h, 0][[1, 2]]}];
  emit[kind <> " exact energy-dependent forward", Factor[schur[h, e][[2, 1]]]];
  emit[kind <> " common,opposite phase response on one link",
       {link[kind, u, u], link[kind, u, 1/u]}];
  (* Four oriented edges around a plaquette. Each constituent's first
     three phases are arbitrary; the last fixes its total holonomy. *)
  topLoop = {t1, t2, t3, wt/(t1 t2 t3)};
  botLoop = {b1, b2, b3, wb/(b1 b2 b3)};
  emit[kind <> " normalized pair loop for face holonomies wt,wb",
       Factor[Times @@ MapThread[link[kind, #1, #2] &, {topLoop, botLoop}]]];
  emit[kind <> " exact open-dimer secular polynomial discrepancy", exactResidual];
  If[exactResidual =!= 0 || Factor[Total[paths] - transfer[h]] =!= 0, Exit[1]],
  {kind, {"particle", "hole"}}];
emit["open-dimer secular polynomial", charTarget];
(* Independent check on the magnitude using the exact two low levels. *)
low[amp_] := (gap - Sqrt[gap^2 + 4 lam^2 amp^2])/2;
splitting = FullSimplify[SeriesCoefficient[
  low[ht - hb] - low[ht + hb], {lam, 0, 2}], gap > 0];
emit["exact low-doublet splitting at order hopping squared", splitting];
If[!TrueQ[Simplify[splitting == 4 ht hb/gap, gap > 0]], Exit[1]];
