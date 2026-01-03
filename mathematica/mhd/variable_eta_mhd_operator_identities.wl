(* ::Package:: *)

(*
  variable_eta_mhd_operator_identities.wl

  Goal:
    Verify identities for variable resistivity/diffusivity η(x,t) in MHD-like induction:
      opGood = -Curl( η Curl(B) )
      opNaive =  η Laplacian(B)

    Show:
      opGood - opNaive == -(Grad η) × (Curl B) - η Grad(Div B)
      Div(opGood) == 0   (constraint preservation)
      Local magnetic energy identity:
        B · opGood == -Div( (η Curl B) × B ) - η |Curl B|^2
*)

Print["============================"];
Print["STEP 3 (FIXED): Variable resistivity operator identities"];
Print["============================\n"];

$RecursionLimit = 4096;

ClearAll[x,y,z,t, η, Bx,By,Bz, Grad3,Div3,Curl3,Lap3,Dot3,Cross3,Norm2];

vars = {x,y,z};

(* --- Basic vector calculus (componentwise) --- *)
Grad3[f_] := {D[f,x], D[f,y], D[f,z]};
Div3[v_List] := D[v[[1]],x] + D[v[[2]],y] + D[v[[3]],z];

Curl3[v_List] := {
  D[v[[3]],y] - D[v[[2]],z],
  D[v[[1]],z] - D[v[[3]],x],
  D[v[[2]],x] - D[v[[1]],y]
};

Lap3[v_List] := {D[v[[1]],x,x] + D[v[[1]],y,y] + D[v[[1]],z,z],
                 D[v[[2]],x,x] + D[v[[2]],y,y] + D[v[[2]],z,z],
                 D[v[[3]],x,x] + D[v[[3]],y,y] + D[v[[3]],z,z]};

Dot3[a_List,b_List] := a[[1]] b[[1]] + a[[2]] b[[2]] + a[[3]] b[[3]];
Cross3[a_List,b_List] := {
  a[[2]] b[[3]] - a[[3]] b[[2]],
  a[[3]] b[[1]] - a[[1]] b[[3]],
  a[[1]] b[[2]] - a[[2]] b[[1]]
};
Norm2[v_List] := Dot3[v,v];

(* --- Symbols --- *)
ηf = η[x,y,z,t];
B  = {Bx[x,y,z,t], By[x,y,z,t], Bz[x,y,z,t]};

curlB = Curl3[B];
divB  = Div3[B];

Print["Operator (recommended):   opGood = -Curl( η Curl(B) )"];
Print["Operator (naive):         opNaive =  η Laplacian(B)\n"];

opGood  = -Curl3[ ηf * curlB ];
opNaive =  ηf * Lap3[B];

Print["============================"];
Print["STEP 3A: Operator identity"];
Print["============================\n"];

expected = -Cross3[ Grad3[ηf], curlB ] - ηf * Grad3[ divB ];

diff = FullSimplify[ opGood - opNaive - expected ];

Print["opGood - opNaive - expected (should be {0,0,0}):"];
Print[diff, "\n"];

Print["If Div(B)=0, expected reduces to  -(Grad η)×(Curl B)."];
diffDivFree = FullSimplify[ (opGood - opNaive + Cross3[Grad3[ηf], curlB]) /. divB -> 0 ];
Print["Div-free check (should be {0,0,0}):"];
Print[diffDivFree, "\n"];

Print["============================"];
Print["STEP 3B: Divergence preservation"];
Print["============================\n"];

divGood  = FullSimplify[ Div3[opGood] ];
divNaive = FullSimplify[ Div3[opNaive] ];

Print["Div(opGood)  (should be 0):"];
Print[divGood, "\n"];

Print["Div(opNaive) (generally nonzero when η varies):"];
Print[divNaive, "\n"];

Print["============================"];
Print["STEP 3C: Magnetic energy dissipation identity"];
Print["============================\n"];

(* Identity used:
   B · (-Curl(η Curl B)) == -Div( (η Curl B) × B ) - η |Curl B|^2
*)
fluxTerm = -Div3[ Cross3[ ηf * curlB, B ] ];
dissTerm = -ηf * Norm2[curlB];

energyCheck = FullSimplify[ Dot3[B, opGood] - (fluxTerm + dissTerm) ];

Print["Check: B·opGood - (flux + diss)  (should be 0):"];
Print[energyCheck, "\n"];

Print["Result (local form):  ∂t(|B|^2/2) = fluxTerm + dissTerm"];
Print["Integrated (periodic/decaying boundary):  d/dt ∫|B|^2 dV = -2 ∫ η |Curl B|^2 dV <= 0\n"];

Print["============================"];
Print["STEP 3D: Vorticity analog"];
Print["============================\n"];
Print["Replace B -> ω (vorticity). Then -Curl(η Curl ω) is the diffusive/reconnective operator that breaks frozen-in topology."];
Print["Done."];

(*"
Output:

============================
STEP 3 (FIXED): Variable resistivity operator identities
============================

Operator (recommended):   opGood = -Curl( η Curl(B) )
Operator (naive):         opNaive =  η Laplacian(B)

============================
STEP 3A: Operator identity
============================

opGood - opNaive - expected (should be {0,0,0}):
{0, 0, 0}

If Div(B)=0, expected reduces to  -(Grad η)×(Curl B).
Div-free check (should be {0,0,0}):
{-(η[x, y, z, t]*(Derivative[1, 0, 1, 0][Bz][x, y, z, t] + Derivative[1, 1, 0, 0][By][x, y, z, t] + Derivative[2, 0, 0, 0][Bx][x, y, z, t])), -(η[x, y, z, t]*(Derivative[0, 1, 1, 0][Bz][x, y, z, t] + Derivative[0, 2, 0, 0][By][x, y, z, t] + Derivative[1, 1, 0, 0][Bx][x, y, z, t])), -(η[x, y, z, t]*(Derivative[0, 0, 2, 0][Bz][x, y, z, t] + Derivative[0, 1, 1, 0][By][x, y, z, t] + Derivative[1, 0, 1, 0][Bx][x, y, z, t]))}

============================
STEP 3B: Divergence preservation
============================

Div(opGood)  (should be 0):
0

Div(opNaive) (generally nonzero when η varies):
Derivative[1, 0, 0, 0][η][x, y, z, t]*(Derivative[0, 0, 2, 0][Bx][x, y, z, t] + Derivative[0, 2, 0, 0][Bx][x, y, z, t] + Derivative[2, 0, 0, 0][Bx][x, y, z, t]) + Derivative[0, 1, 0, 0][η][x, y, z, t]*(Derivative[0, 0, 2, 0][By][x, y, z, t] + Derivative[0, 2, 0, 0][By][x, y, z, t] + Derivative[2, 0, 0, 0][By][x, y, z, t]) + Derivative[0, 0, 1, 0][η][x, y, z, t]*(Derivative[0, 0, 2, 0][Bz][x, y, z, t] + Derivative[0, 2, 0, 0][Bz][x, y, z, t] + Derivative[2, 0, 0, 0][Bz][x, y, z, t]) + η[x, y, z, t]*(Derivative[0, 0, 3, 0][Bz][x, y, z, t] + Derivative[0, 1, 2, 0][By][x, y, z, t] + Derivative[0, 2, 1, 0][Bz][x, y, z, t] + Derivative[0, 3, 0, 0][By][x, y, z, t] + Derivative[1, 0, 2, 0][Bx][x, y, z, t] + Derivative[1, 2, 0, 0][Bx][x, y, z, t] + Derivative[2, 0, 1, 0][Bz][x, y, z, t] + Derivative[2, 1, 0, 0][By][x, y, z, t] + Derivative[3, 0, 0, 0][Bx][x, y, z, t])

============================
STEP 3C: Magnetic energy dissipation identity
============================

Check: B·opGood - (flux + diss)  (should be 0):
0

Result (local form):  ∂t(|B|^2/2) = fluxTerm + dissTerm
Integrated (periodic/decaying boundary):  d/dt ∫|B|^2 dV = -2 ∫ η |Curl B|^2 dV <= 0

============================
STEP 3D: Vorticity analog
============================

Replace B -> ω (vorticity). Then -Curl(η Curl ω) is the diffusive/reconnective operator that breaks frozen-in topology.
"*)
