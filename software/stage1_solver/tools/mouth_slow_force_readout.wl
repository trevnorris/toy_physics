(* Use a velocity/displacement work decomposition at the turning and crossing
   points, then extrapolate the independent interior-trace carrier readout. *)
force[t_] := bias + amp ((k - mass w^2) Cos[t] - drag w Sin[t]);
atCos = Simplify[(force[0] - force[Pi])/(2 amp)];
atSin = Simplify[(force[Pi/2] - force[3 Pi/2])/(2 amp)];
signedMass = Simplify[(k - atCos)/w^2];
cycleWork = Integrate[force[t] D[amp Cos[t], t], {t, 0, 2 Pi}];
Print["turning force per amplitude, crossing force per amplitude, recovered mass, cycle work = ",
  ToString[{atCos, atSin, signedMass, cycleWork}, InputForm]];
If[signedMass =!= mass || Simplify[-atSin/w - drag] =!= 0 ||
   Simplify[cycleWork - Pi amp^2 w drag] =!= 0, Exit[1]];

manifest = Import["docs/_measurements/lateral_mouth_acceleration.json", "RawJSON"];
folder = "_scratch/s11c/lateral_accel_e24_c6_01";
matches = Select[manifest["runs"], MemberQ[StringSplit[#["command"]], folder] &];
If[Length[matches] != 1, Exit[1]];
record = First[matches];
If[FileHash[folder <> "/stdout", "SHA256", "HexString"] =!= record["stdout_sha256"], Exit[1]];
lines = Select[record["stdout"], StringStartsQ[#, "omega="] &];
(* These two stored columns contain nonnegative, ordinary decimal literals. *)
number[s_] := Module[{parts},
  If[!StringMatchQ[s, RegularExpression["[0-9]+(\\.[0-9]+)?"]], Exit[1]];
  parts = StringSplit[s, "."];
  FromDigits[StringJoin[parts]]/If[Length[parts] == 2, 10^StringLength[Last[parts]], 1]
];
rows = Table[With[{a = Association[Rule @@@ (StringSplit[#, "="] & /@ StringSplit[line])]},
      {number[a["omega"]], number[a["direct_compensating_per_A"]], number[a["compensating_per_A"]]}],
    {line, lines}];
If[rows[[All, 1]] =!= {0, 1/200, 1/100}, Exit[1]];
Do[
  fit = InterpolatingPolynomial[{#[[1]]^2, #[[column]]} & /@ Rest[rows], z];
  observed = rows[[1, column]];
  Print[{"interior-trace flux", "boundary-law flux"}[[column - 1]],
    " slow extrapolation, archived slow tangent, signed difference = ",
    ToString[N[{fit /. z -> 0, observed, (fit /. z -> 0) - observed}, 15], InputForm]],
  {column, 2, 3}];
