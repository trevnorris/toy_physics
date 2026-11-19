(* Independent wave boundary matrix, linear contact transfer and link counting. *)
ClearAll["Global`*"];
emit[name_, value_] := Print[name, " = ", ToString[value, InputForm]];
verify[test_] := If[!TrueQ[test], Print["Failed check: ", HoldForm[test]]; Exit[1]];

(* The lab depth runs from 0 to s L. Propagate displacement and lab slope. *)
prop[z_] := {{Cos[k z], Sin[k z]/k}, {-k Sin[k z], Cos[k z]}};
Do[
  nodeToEnd = prop[s len];
  endCondition = nodeToEnd[[2, 2]] /. k -> angle/len;
  firstK = (angle /. First[Solve[{endCondition == 0, 0 < angle < Pi}, angle, Reals]])/len;
  (* Use the end amplitude to normalize the two real tangent displacements. *)
  endMotion = amp {Cos[phase], p Sin[phase]};
  rimSlope = Simplify[endMotion / (nodeToEnd[[1, 2]] /. k -> firstK), len > 0];
  squaredTraction = Simplify[ten^2 rimSlope.rimSlope];
  disp = prop[z][[1, 2]] rimSlope /. k -> firstK;
  power = -ten speed firstK D[disp, phase].D[disp, z];
  emit[{"side", s, "polarization", p, "frequency, squared rim traction, mean power flux"},
       {speed firstK, squaredTraction, Simplify[Integrate[power, {phase, 0, 2 Pi}]/(2 Pi)]}],
  {s, {1, -1}}, {p, {1, -1}}];

(* General real symmetric elastic load in the circular basis. *)
basis = {{1, 1}, {I, -I}}/Sqrt[2];
emit["elastic end in circular basis",
  Simplify[ConjugateTranspose[basis].{{aa, bb}, {bb, dd}}.basis]];

(* An extra real, skew velocity coupling is a finite-mode sensitivity control.
   The positive-frequency phasor {1,-I p} uses Exp[I omega t]. *)
Do[
  impedance = mass (omega0^2 - omega^2) IdentityMatrix[2] +
              I omega s gamma {{0, -1}, {1, 0}};
  polarization = {1, -I p};
  eigen = Expand[(impedance.polarization)[[1]]];
  verify[Simplify[impedance.polarization == eigen polarization]];
  shifted = FullSimplify[omega /. Solve[eigen == 0, omega],
                         mass > 0 && omega0 > 0 && gamma > 0];
  positive = Select[shifted, FullSimplify[# > 0, mass > 0 && omega0 > 0 && gamma > 0] &];
  verify[Length[positive] == 1];
  emit[{"added gyroscopic term, side", s, "polarization", p}, First[positive]],
  {s, {1, -1}}, {p, {1, -1}}];

(* Mirror the kinetic operator and the domain of the contact.
   Plane-wave outgoing impedance gives transmission from either lab side. *)
op[f_, v_] := -cc D[f, {w, 2}] + v f;
emit["mirror operator residual", Simplify[op[f[-w], v[-w]] - (op[f[w], v[w]] /. w -> -w)]];
jump = {{1, 0}, {h/cc, 1}};
Do[
  waves = {{1, 1}, {s I k, -s I k}};
  transfer = Simplify[Inverse[waves].MatrixPower[jump, s].waves];
  sol = First[Solve[Thread[{t, 0} == transfer.{1, r}], {r, t}]];
  emit[{"side", s, "current transmission"},
       FullSimplify[Abs[t /. sol]^2, cc > 0 && k > 0 && Element[h, Reals]]],
  {s, {1, -1}}];
(* Match outgoing tails at imaginary wavenumber. *)
decay = kap /. First[Solve[Det[Transpose[{jump.{1, kap}, {1, -kap}}]] == 0, kap]];
emit["bound mode inverse length and energy, decaying sector", {decay, -cc decay^2}];
emit["bound mode current", ComplexExpand[Im[Conjugate[b Exp[-kap w]] D[b Exp[-kap w], w]]]];

(* Count sectors via the generating polynomial of oriented link fluxes.
   A shared link contributes opposite powers to its two incident vertices. *)
gen = Expand[(x/y + y/x) (x + 1/x)^3 (y + 1/y)^3 x^4 y^4];
mult[q_] := Coefficient[Coefficient[gen, x, 4 + 2 q], y, 4 - 2 q];
sectors = Select[Range[-2, 2], mult[#] > 0 &];
emit["neutral pair multiplicities", Table[{q, mult[q]}, {q, sectors}]];
pairEnergy[s_, q_, bias_] := q^2 - 2 bias s q;
Do[
  emit[{"lambda/U", bias, "left side", s, "ground charge"},
       MinimalBy[sectors, pairEnergy[s, #, bias] &]],
  {bias, {0, 1, 2}}, {s, {1, -1}}];
Do[
  rivals = DeleteCases[sectors, s];
  emit[{"strict unit-selection interval", s},
       Reduce[And @@ (pairEnergy[s, s, ratio] < pairEnergy[s, #, ratio] & /@ rivals), ratio, Reals]],
  {s, {1, -1}}];
