(* Invert using inlet stagnation pressure above outlet ambient. *)
ClearAll["Global`*"];
area=4 Pi a^3/3;
pull=c0+4 Pi sig a^2; speedOut=Sqrt[2 head/rho]; speedIn=m speedOut;
pin=rho (speedOut^2-speedIn^2)/2;
d=area rho speedIn;
tNet=Factor[d speedIn + area pin - d speedOut];
needed=head/.First[Solve[tNet==pull,head]];
ratio=Factor[(needed/.a->2 a)/needed];
Print["head closure and inlet-rate squared = ",ToString[Factor/@{needed,d^2/.head->needed},InputForm]];
Print["doubling ratio for head = ",ToString[ratio,InputForm]];
Print["scalar load at this balance = ",ToString[Factor[ell z tNet/.head->needed],InputForm]];
Print["curvature and tension limits = ",ToString[Factor/@{ratio/.sig->0,ratio/.c0->0},InputForm]];
If[Factor[tNet-area head (m-1)^2]=!=0 || Factor[(needed/.sig->0)area (m-1)^2-c0]=!=0,Exit[1]];
Exit[0];
