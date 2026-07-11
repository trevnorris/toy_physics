(* MIDWAY knob-audit codimension dry-run: independent Mathematica CAD engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
assertTrue[condition_, msg_] := If[! TrueQ[condition], fail[msg]];

hbar = Unique["hbar"];
mass = Unique["mass"];
bigK = Unique["K"];
rho0 = Unique["rho0"];
cs0 = Unique["cs0"];
xiH = Unique["xiH"];
h0 = Unique["h0"];
scaleA = Unique["a"];
cGamma = Unique["cGamma"];
lambdaGamma = Unique["lambdaGamma"];

mediumVars = {hbar, mass, bigK, rho0, cs0, xiH, h0, scaleA, cGamma, lambdaGamma};
mediumBaseline = {
  mass*cs0^2 - 5*bigK*rho0^4,
  scaleA*mass*cs0 - hbar,
  xiH^2*mass^2*cs0^2 - 2*hbar^2,
  4*h0 - mass*cs0^2,
  lambdaGamma*cs0 - cGamma
};

muEta = Unique["muEta"];
tW = Unique["tW"];
beta = Unique["beta"];
kEta = Unique["kEta"];
vp0OverLc = Unique["vp0OverLc"];
tOmega = Unique["tOmega"];
aB = Unique["aB"];
kappaB = Unique["kappaB"];
delta = Unique["delta"];
sigmaWall = Unique["sigmaWall"];

wallVars = {muEta, tW, beta, kEta, vp0OverLc, tOmega, aB, kappaB, delta, sigmaWall};
wallBaseline = {
  kEta - tW*beta^2,
  2*aB*delta^2 - kappaB,
  36*sigmaWall^2 - 2*aB*kappaB
};

blocks = <|
  "M" -> <|
    "vars" -> mediumVars,
    "cases" -> <|
      "baseline" -> mediumBaseline,
      "C-M1" -> Join[Most[mediumBaseline], {lambdaGamma - lambdaGamma}],
      "C-M2" -> Join[mediumBaseline, {xiH^2 - 2*scaleA^2}],
      "C-M3" -> Join[mediumBaseline, {bigK - rho0}]
    |>
  |>,
  "Wchi" -> <|
    "vars" -> wallVars,
    "cases" -> <|
      "baseline" -> wallBaseline,
      "C-W1" -> Join[{kEta - kEta}, Rest[wallBaseline]],
      "C-W2" -> Join[wallBaseline, {muEta - tW}]
    |>
  |>
|>;

expectedPayload = <|
  "M" -> <|
    "baseline" -> <|"dim_before" -> 10, "dim_after" -> 5, "Delta" -> 5|>,
    "C-M1" -> <|"dim_before" -> 10, "dim_after" -> 6, "Delta" -> 4|>,
    "C-M2" -> <|"dim_before" -> 10, "dim_after" -> 5, "Delta" -> 5|>,
    "C-M3" -> <|"dim_before" -> 10, "dim_after" -> 4, "Delta" -> 6|>
  |>,
  "Wchi" -> <|
    "baseline" -> <|"dim_before" -> 10, "dim_after" -> 7, "Delta" -> 3|>,
    "C-W1" -> <|"dim_before" -> 10, "dim_after" -> 8, "Delta" -> 2|>,
    "C-W2" -> <|"dim_before" -> 10, "dim_after" -> 6, "Delta" -> 4|>
  |>
|>;

regionDimension[vars_, equations_] := Module[{formula, region, result},
  formula = (And @@ Thread[equations == 0]) && (And @@ Thread[vars > 0]);
  region = ImplicitRegion[formula, Evaluate[vars]];
  result = TimeConstrained[RegionDimension[region], 60, "DIMENSION_UNCERTIFIED"];
  If[! IntegerQ[result],
    Print["DIMENSION_UNCERTIFIED: equations=", equations];
    fail["CAD RegionDimension did not return an integer"]
  ];
  result
];

computeCase[vars_, equations_] := Module[{before, after},
  before = regionDimension[vars, {}];
  after = regionDimension[vars, equations];
  <|"dim_before" -> before, "dim_after" -> after, "Delta" -> before - after|>
];

payload = AssociationMap[
  Function[blockName,
    With[{block = blocks[blockName]},
      AssociationMap[
        Function[caseName, computeCase[block["vars"], block["cases"][caseName]]],
        Keys[block["cases"]]
      ]
    ]
  ],
  Keys[blocks]
];

(* The SymPy engine emits this identical canonical association as compact JSON.
   All values here were independently derived by real CAD RegionDimension. *)
assertTrue[SameQ[payload, expectedPayload], "dimension payload mismatch: " <> ToString[payload, InputForm]];

Print["BLOCK  CASE      dim_before  dim_after  Delta"];
KeyValueMap[
  Function[{blockName, cases},
    KeyValueMap[
      Function[{caseName, record},
        Print[StringRiffle[{
          blockName,
          caseName,
          ToString[record["dim_before"]],
          ToString[record["dim_after"]],
          ToString[record["Delta"]]
        }, "  "]]
      ],
      cases
    ]
  ],
  payload
];
Print["NOTE: beta_2(w) is a PROFILE and the R35 radial-scalar reductions " <>
  "{Mtilde,Ktilde,Ttilde_Omega}=Integral[density*beta_2^2] are integral " <>
  "functionals, not polynomializable; they are outside this computed check " <>
  "and remain reasoned-tally only."];
Print["COMPARISON_PAYLOAD: " <> ExportString[payload, "RawJSON", "Compact" -> True]];
Print["MIDWAY_CODIM_CHECK: CONSISTENT"];

Exit[0];
