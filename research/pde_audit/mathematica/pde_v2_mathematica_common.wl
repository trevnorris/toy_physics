ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 78]];
  Print[title];
  Print[StringRepeat["=", 78]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 78]];
  Print[title];
  Print[StringRepeat["-", 78]];
);

scriptDir[] := DirectoryName[$InputFileName];
auditDir[] := DirectoryName[DirectoryName[$InputFileName]];
fixtureDir[] := FileNameJoin[{auditDir[], "scripts", "fixtures"}];
fixturePath[name_String] := FileNameJoin[{fixtureDir[], name}];
auditPath[name_String] := FileNameJoin[{auditDir[], name}];
artifactDir[] := FileNameJoin[{auditDir[], "scripts", "output", "artifacts"}];
artifactPath[name_String] := FileNameJoin[{artifactDir[], name}];

fileSHA256[path_String] := ToLowerCase[IntegerString[FileHash[path, "SHA256"], 16, 64]];

jsonImport[path_String] := Module[{},
  If[!FileExistsQ[path],
    Print["FAIL: missing JSON file: ", path];
    Exit[1];
  ];
  Import[path, "RawJSON"]
];

pass[name_String] := Print["PASS: ", name];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail: ", fmt[detail]]];
  Exit[1];
);

scalarSimplify[expr_] := Quiet[Cancel[Together[expr]]];

deepSimplify[expr_] := If[
  ListQ[expr],
  Map[scalarSimplify, expr, {ArrayDepth[expr]}],
  scalarSimplify[expr]
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# == 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr == 0]
];

checkTrue[name_String, condition_, detail_: Missing["NotAvailable"]] := If[
  TrueQ[condition],
  pass[name],
  fail[name, detail]
];

checkZero[name_String, expr_] := Module[{res},
  res = deepSimplify[expr];
  If[allZeroQ[res], pass[name], fail[name, res]];
];

checkEqual[name_String, lhs_, rhs_] := checkZero[name, lhs - rhs];

nearZeroQ[x_, tol_: 10^-9] := TrueQ[Abs[N[x]] <= tol];

weightedBar[x20_, x21_, x22_] := (x20 + 2 x21 + 2 x22)/5;
groupA[x20_, x21_, x22_] := (2 x20 - x21 - x22)/10;
groupB[x20_, x21_, x22_] := (x21 - x22)/2;

responseU2[D0_, D2_] := -D2/D0;
responseU4[D0_, D2_, D4_] := (D2^2 - D0 D4)/D0^2;
prefactorP0[D0_, N0_] := N0/D0;
prefactorP2[D0_, D2_, N0_, N2_] := (D0 N2 - 2 D2 N0)/D0^2;
prefactorP4[D0_, D2_, D4_, N0_, N2_, N4_] :=
  (D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0)/D0^3;
