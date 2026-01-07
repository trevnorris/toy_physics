(* ============================================================
   cylinder_throat.wl  — Paper 7 Step 1
   Cylinder throat DtN / impedance map with clean CLI output.

   CLI:
     math -script mathematica/inner_throat/cylinder_throat.wl -- -cli
     math -script mathematica/inner_throat/cylinder_throat.wl -- -cli -scan
   ============================================================ *)

ClearAll["Global`*"];

(* ----------------------------
   Robust CLI detection
   ---------------------------- *)
scriptArgs = Quiet@Check[$ScriptCommandLine, {}];
allArgs = DeleteDuplicates@Flatten@{scriptArgs, $CommandLine};

cliMode =
  ($FrontEnd === Null) || MemberQ[allArgs, "-cli"];

runScan =
  MemberQ[allArgs, "-scan"];

plotMode = !cliMode;

If[cliMode, Print["[CLI mode enabled]"];];

(* ----------------------------
   User parameters (edit these)
   ---------------------------- *)
c      = 1.0;         (* wave speed *)
a      = 1.0;         (* throat radius *)
Lovera = 1.85;        (* aspect ratio L/a *)
L      = Lovera*a;    (* throat length *)

wallBC   = "Neumann"; (* "Neumann" or "Dirichlet" at r=a *)
bottomBC = "Neumann"; (* "Neumann" or "Dirichlet" at z=-L *)

gamma  = 0.01;        (* small damping: ω -> ω + i γ *)
mMax   = 0;           (* azimuthal modes m=0..mMax *)
nMax   = 4;           (* radial mode count; Dirichlet: 1..nMax, Neumann: 0..nMax *)

ωmin = 0.0;
ωmax = 30.0;

(* picked mode *)
mPick = 0;
nPick = If[wallBC === "Dirichlet", 1, 0];

(* ----------------------------
   CLI-safe numeric formatting
   IMPORTANT: use OutputForm so it prints as plain text, not NumberForm[...]
   ---------------------------- *)
NumStr[x_, prec_:12] := Module[{nf},
  nf = NumberForm[
    N[x, prec],
    {Infinity, prec},
    NumberPadding -> {"", ""},
    ExponentFunction -> (Null &)
  ];
  ToString[nf, OutputForm]
];

PrintTable[header_List, rows_List] := Module[{},
  Print[StringRiffle[header, "\t"]];
  Scan[
    (Print[StringRiffle[(ToString[#, OutputForm] & /@ #), "\t"]]) &,
    rows
  ];
];

(* ----------------------------
   Bessel derivative (stable recurrence)
   J_m'(x) = (J_{m-1}(x) - J_{m+1}(x))/2
   ---------------------------- *)
BesselJPrime[m_Integer, x_?NumericQ] := (BesselJ[m - 1, x] - BesselJ[m + 1, x])/2;

(* ----------------------------
   Radial roots for wall boundary
   ---------------------------- *)
RadialRoot::bc  = "Unknown boundary condition `1`.";
RadialRoot::dir = "Dirichlet roots use n>=1; got n=`1`.";

RadialRoot[bc_, m_Integer, n_Integer] := Module[{x0, x1, guess, sol},
  Which[
    bc === "Dirichlet",
      If[n < 1, Message[RadialRoot::dir, n]; Return[$Failed]];
      BesselJZero[m, n],

    bc === "Neumann",
      If[m == 0 && n == 0, 0,
        (* interlacing midpoint guess *)
        x0 = BesselJZero[m, Max[n, 1]];
        x1 = BesselJZero[m, Max[n, 1] + 1];
        guess = (x0 + x1)/2;

        sol = FindRoot[BesselJPrime[m, x] == 0, {x, guess},
                       WorkingPrecision -> 30, AccuracyGoal -> 16, PrecisionGoal -> 16];
        x /. sol
      ],

    True,
      Message[RadialRoot::bc, bc];
      $Failed
  ]
];

Lambda[bc_, m_Integer, n_Integer] := RadialRoot[bc, m, n]/a;

(* ----------------------------
   κ_{mn}(ω) and DtN eigenvalue Z_{mn}(ω)
   ---------------------------- *)
Kappa[ω_?NumericQ, m_Integer, n_Integer, γ_:0.0] := Module[{w, k, λ},
  w = ω + I*γ;
  k = w/c;
  λ = Lambda[wallBC, m, n];
  Sqrt[k^2 - λ^2]
];

Zmn::bc = "Unknown bottom boundary condition `1`.";

Zmn[ω_?NumericQ, m_Integer, n_Integer, γ_:0.0] := Module[{κ},
  κ = Kappa[ω, m, n, γ];
  Which[
    bottomBC === "Dirichlet", κ*Cot[κ*L],
    bottomBC === "Neumann",  -κ*Tan[κ*L],
    True, (Message[Zmn::bc, bottomBC]; $Failed)
  ]
];

(* ----------------------------
   Resonance predictor (pole locations)
   ---------------------------- *)
ResonanceOmega[m_Integer, n_Integer, q_Integer] := Module[{λ, κres},
  λ = Lambda[wallBC, m, n];
  κres = Which[
    bottomBC === "Dirichlet", (q*Pi)/L,
    bottomBC === "Neumann",  ((q + 1/2)*Pi)/L,
    True, $Failed
  ];
  c*Sqrt[λ^2 + κres^2]
];

(* ----------------------------
   CLI summary
   ---------------------------- *)
If[cliMode,
  Print["=== Cylinder DtN Summary ==="];
  Print["wallBC=", wallBC, "  bottomBC=", bottomBC,
        "  a=", NumStr[a, 6], "  L=", NumStr[L, 6],
        "  (L/a)=", NumStr[L/a, 6]];
  Print["picked mode: m=", mPick, " n=", nPick];
];

(* ----------------------------
   Mode table
   ---------------------------- *)
modeRows = Flatten[
  Table[
    Table[
      Module[{lam = Lambda[wallBC, m, n]},
        {ToString[m], ToString[n], NumStr[lam, 12], NumStr[c*lam, 12]}
      ],
      {n, If[wallBC === "Dirichlet", 1, 0], nMax}
    ],
    {m, 0, mMax}
  ],
  1
];

If[cliMode,
  Print["--- Mode table (m, n, λ, ω_cutoff=c λ) ---"];
  PrintTable[{"m","n","lambda","omega_cutoff"}, modeRows];
];

(* ----------------------------
   Resonances for picked mode
   ---------------------------- *)
qMax = 6;
resList = Table[{q, ResonanceOmega[mPick, nPick, q]}, {q, 0, qMax}];

If[cliMode,
  Print["--- Approx resonance frequencies ω_res(q) for picked mode ---"];
  PrintTable[{"q","omega_res"},
    ({ToString[#[[1]]], NumStr[#[[2]], 15]} & /@ resList)
  ];

  If[wallBC === "Neumann" && mPick == 0 && nPick == 0,
    Print["Expected ω0 = (π/2) c / L = ", NumStr[(Pi/2)*c/L, 15]];
    Print["Expected spacing Δω = π c / L = ", NumStr[Pi*c/L, 15]];
    Print["Measured spacing (q=0→1) = ", NumStr[resList[[2,2]] - resList[[1,2]], 15]];
  ];
];

(* ----------------------------
   Power proxy diagnostics:
   print both ω Im(Z) and -ω Im(Z) (normal/sign convention)
   ---------------------------- *)
PproxyPlus[ω_?NumericQ] := ω*Im[Zmn[ω, mPick, nPick, gamma]];
PproxyMinus[ω_?NumericQ] := -ω*Im[Zmn[ω, mPick, nPick, gamma]];

If[cliMode,
  Print["--- Power-proxy sign check (sample points, γ=", NumStr[gamma, 6], ") ---"];
  samples = {0.5, 1.0, 2.0, 3.0, 5.0};
  ppRows = Table[
    {NumStr[ω, 6], NumStr[PproxyPlus[ω], 12], NumStr[PproxyMinus[ω], 12]},
    {ω, samples}
  ];
  PrintTable[{"ω","ω Im(Z)","-ω Im(Z)"}, ppRows];
];

(* ----------------------------
   Optional: quick L/a scan proxy (smoke test)
   ---------------------------- *)
If[cliMode && runScan,
  Print["--- Running quick L/a scan proxy (this can take a bit) ---"];

  scanLovera = Range[1.2, 2.6, 0.05];

  selectionProxy[LoveraVal_?NumericQ] := Module[{Ltmp, λ, ωgrid, vals, Ztmp},
    Ltmp = LoveraVal*a;
    λ = Lambda[wallBC, mPick, nPick];

    Ztmp[ω_?NumericQ] := Module[{wloc, kloc, κloc},
      wloc = ω + I*gamma;
      kloc = wloc/c;
      κloc = Sqrt[kloc^2 - λ^2];
      Which[
        bottomBC === "Dirichlet", κloc*Cot[κloc*Ltmp],
        bottomBC === "Neumann",  -κloc*Tan[κloc*Ltmp]
      ]
    ];

    ωgrid = Subdivide[ωmin + 10^-3, ωmax, 2500];
    vals = Re[Ztmp /@ ωgrid];
    Max[vals]
  ];

  scanData = Table[{x, selectionProxy[x]}, {x, scanLovera}];

  Print["L/a\tmax_ω Re(Z) (crude)"];
  Scan[
    (Print[NumStr[#[[1]], 4], "\t", NumStr[#[[2]], 10]]) &,
    scanData
  ];
];

If[cliMode && !runScan,
  Print["(scan skipped; add -scan to run L/a proxy scan)"];
];

(* ----------------------------
   Plotting (Notebook mode only)
   ---------------------------- *)
If[plotMode,
  zReImPlot =
    Plot[
      Evaluate[{Re[Zmn[ω, mPick, nPick, gamma]], Im[Zmn[ω, mPick, nPick, gamma]]}],
      {ω, ωmin, ωmax},
      PlotRange -> All,
      AxesLabel -> {"ω", "Z(ω)"},
      PlotLegends -> {"Re Z", "Im Z"},
      PlotLabel -> Row[{
        "value Z_{", mPick, ",", nPick, "}(ω),  wallBC=", wallBC,
        ", bottomBC=", bottomBC, ",  L/a=", Lovera
      }]
    ];

  zMagPhasePlot =
    Plot[
      Evaluate[{Abs[Zmn[ω, mPick, nPick, gamma]], Arg[Zmn[ω, mPick, nPick, gamma]]}],
      {ω, ωmin, ωmax},
      PlotRange -> All,
      AxesLabel -> {"ω", ""},
      PlotLegends -> {"|Z|", "Arg(Z)"},
      PlotLabel -> "Magnitude and phase"
    ];

  passivityPlotPlus =
    Plot[
      Evaluate[PproxyPlus[ω]],
      {ω, ωmin, ωmax},
      PlotRange -> All,
      AxesLabel -> {"ω", "ω Im(Z)"},
      PlotLabel -> "Power-proxy (potential-style):  ω Im(Z)"
    ];

  passivityPlotMinus =
    Plot[
      Evaluate[PproxyMinus[ω]],
      {ω, ωmin, ωmax},
      PlotRange -> All,
      AxesLabel -> {"ω", "-ω Im(Z)"},
      PlotLabel -> "Power-proxy (flipped normal):  -ω Im(Z)"
    ];

  Print[zReImPlot];
  Print[zMagPhasePlot];
  Print[passivityPlotPlus];
  Print[passivityPlotMinus];
];

(*"
Output:

[CLI mode enabled]
=== Cylinder DtN Summary ===
wallBC=Neumann  bottomBC=Neumann  a=1.  L=1.85  (L/a)=1.85
picked mode: m=0 n=0
--- Mode table (m, n, λ, ω_cutoff=c λ) ---
m	n	lambda	omega_cutoff
0	0	0.	0.
0	1	3.831705970208	3.831705970208
0	2	7.015586669816	7.015586669816
0	3	10.173468135063	10.173468135063
0	4	13.323691936314	13.323691936314
--- Approx resonance frequencies ω_res(q) for picked mode ---
q	omega_res
0	0.849079095564809
1	2.547237286694427
2	4.245395477824044
3	5.943553668953662
4	7.641711860083279
5	9.3398700512129
6	11.03802824234251
Expected ω0 = (π/2) c / L = 0.849079095564809
Expected spacing Δω = π c / L = 1.698158191129618
Measured spacing (q=0→1) = 1.698158191129618
--- Power-proxy sign check (sample points, γ=0.01) ---
ω	ω Im(Z)	-ω Im(Z)
0.5	-0.019388493462	0.019388493462
1.	-0.207819743789	0.207819743789
2.	-0.115344811521	0.115344811521
3.	-0.2744445218	0.2744445218
5.	-0.46803556568	0.46803556568
"*)
