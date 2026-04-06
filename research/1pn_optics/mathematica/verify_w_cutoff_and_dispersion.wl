(* ============================================================
   verify_w_cutoff_and_dispersion.wl

   Purpose:
     For transverse mode index m and bulk scale Lw, compute:
       kw (Neumann or Periodic),
       cutoff frequency fc = (c_s kw)/(2π),
       and the dispersion factor
         F(ω) = sqrt(1 + (kw/k)^2) = 1/sqrt(1 - (ωc/ω)^2),
       where ωc = c_s kw and ω = 2π f.

     This avoids the misleading "band by wavelength" issue for m>0.

   ============================================================ *)

ClearAll["Global`*"];

(* --- knobs --- *)
cS = 3.*10^8;                 (* use c for far-zone estimate; adjust if needed *)
mList  = {0, 1};              (* zero + first excited *)
LwList = {1.*10^-6, 1.*10^-3, 1., 10.};   (* meters *)

freqs = {
  {"1 GHz (radio)",        1.*10^9},
  {"10 GHz (microwave)",   1.*10^10},
  {"1 THz (far IR)",       1.*10^12},
  {"100 THz (IR)",         1.*10^14},
  {"500 THz (optical)",    5.*10^14},
  {"1 PHz (near UV)",      1.*10^15}
};

kwNeumann[m_, Lw_] := m Pi/Lw;
kwPeriodic[m_, Lw_] := 2. Pi m/Lw;

omega[f_] := 2. Pi f;

(* dispersion factor for group-index multiplier relative to the k_w=0 branch *)
Ffactor[om_, omc_] := If[om <= omc, Indeterminate, 1./Sqrt[1. - (omc/om)^2]];
deltaFromF[F_] := If[F === Indeterminate, Indeterminate, F - 1.];

fmt[x_] := If[x === Indeterminate, "NO-PROP",
  ToString[FortranForm[N[x, 6]]]
];

printBlock[modeName_, kwFun_] := Module[{},
  Print["\n======================================================="];
  Print[modeName];
  Print["======================================================="];
  Print["Columns: Lw, m, kw(1/m), fc(Hz), then Δ(f)=F-1 for each test frequency"];
  Print[""];

  Do[
    Do[
      Module[{kw = kwFun[m, Lw], omc, fc, deltas},
        omc = cS*kw;
        fc  = omc/(2. Pi);
        deltas = Table[
          deltaFromF[Ffactor[omega[freqs[[i,2]]], omc]],
          {i, Length[freqs]}
        ];

        Print[
          "Lw="<>fmt[Lw]<>
          "  m="<>ToString[m]<>
          "  kw="<>fmt[kw]<>
          "  fc="<>fmt[fc]<>
          "  Δ=["<>StringRiffle[fmt /@ deltas, ", "]<>"]"
        ];
      ],
      {m, mList}
    ],
    {Lw, LwList}
  ];

  Print["\nLegend (Δ list order): ", StringRiffle[freqs[[All,1]], " | "]];
];

printBlock["Neumann modes", kwNeumann];
printBlock["Periodic modes", kwPeriodic];

Print["\nInterpretation:"];
Print["  • m=0 => kw=0 => fc=0 => Δ=0 for all frequencies (clean photon branch)."];
Print["  • m>=1 => fc ~ c/(2Lw) (Neumann) or fc ~ c/Lw (Periodic): a cutoff exists."];
Print["  • For f <= fc: NO-PROP (evanescent along the brane)."];
Print["  • Just above fc: Δ becomes large (strong dispersion near cutoff)."];
Print["\nDone."];

(*"
Output:

=======================================================
Neumann modes
=======================================================
Columns: Lw, m, kw(1/m), fc(Hz), then Δ(f)=F-1 for each test frequency

Lw=1.e-6  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=1.e-6  m=1  kw=3.141592653589793e6  fc=1.5e14  Δ=[NO-PROP, NO-PROP, NO-PROP, NO-PROP, 0.04828483672191819, 0.01144347484834718]
Lw=0.001  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=0.001  m=1  kw=3141.592653589793  fc=1.5e11  Δ=[NO-PROP, NO-PROP, 0.01144347484834718, 1.1250018983055554e-6, 4.500000305718288e-8, 1.1250000264695359e-8]
Lw=1.  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=1.  m=1  kw=3.141592653589793  fc=1.5e8  Δ=[0.01144347484834718, 0.00011251898793540605, 1.1250000264695359e-8, 1.1251000131551336e-12, 4.5075054799781356e-14, 1.1324274851176597e-14]
Lw=10.  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=10.  m=1  kw=0.3141592653589793  fc=1.5e7  Δ=[0.00011251898793540605, 1.1250018983055554e-6, 1.1250000930829174e-10, 1.1324274851176597e-14, 4.440892098500626e-16, 2.220446049250313e-16]

Legend (Δ list order): 1 GHz (radio) | 10 GHz (microwave) | 1 THz (far IR) | 100 THz (IR) | 500 THz (optical) | 1 PHz (near UV)

=======================================================
Periodic modes
=======================================================
Columns: Lw, m, kw(1/m), fc(Hz), then Δ(f)=F-1 for each test frequency

Lw=1.e-6  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=1.e-6  m=1  kw=6.283185307179586e6  fc=3.e14  Δ=[NO-PROP, NO-PROP, NO-PROP, NO-PROP, 0.25, 0.04828483672191819]
Lw=0.001  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=0.001  m=1  kw=6283.185307179586  fc=3.e11  Δ=[NO-PROP, NO-PROP, 0.04828483672191819, 4.5000303752207316e-6, 1.8000004864404673e-7, 4.500000305718288e-8]
Lw=1.  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=1.  m=1  kw=6.283185307179586  fc=3.e8  Δ=[0.04828483672191819, 0.0004503039779921725, 4.500000305718288e-8, 4.5001780080156095e-12, 1.800781745942004e-13, 4.5075054799781356e-14]
Lw=10.  m=0  kw=0.  fc=0.  Δ=[0., 0., 0., 0., 0., 0.]
Lw=10.  m=1  kw=0.6283185307179586  fc=3.e7  Δ=[0.0004503039779921725, 4.5000303752207316e-6, 4.5000003723316695e-10, 4.5075054799781356e-14, 1.7763568394002505e-15, 4.440892098500626e-16]

Legend (Δ list order): 1 GHz (radio) | 10 GHz (microwave) | 1 THz (far IR) | 100 THz (IR) | 500 THz (optical) | 1 PHz (near UV)
"*)
