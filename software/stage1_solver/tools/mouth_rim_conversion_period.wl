(* Compact flux period of the stationary parcel impulse. Check the
   antiunitary identities for arbitrary real onsite spin fields and hops,
   and measure the period after adding a passive pin-dependent site energy. *)
ClearAll["Global`*"];
emit[label_, value_] := WriteString["stdout",label," = ",
    ToString[value,InputForm],"\n"];
id = IdentityMatrix[2]; sigma = Table[PauliMatrix[k], {k,3}];
lead = KroneckerProduct[IdentityMatrix[3], (id-sigma[[3]])/2];
up = (IdentityMatrix[6]-lead)/3; down = lead/3;
gTheta = KroneckerProduct[DiagonalMatrix[{1,-1,1}],I sigma[[2]]];
transform[m_] := gTheta.Conjugate[m].ConjugateTranspose[gTheta];

ham[p_, hops_, flux_, siteEnergy_] := Module[{t01,t12,t20,hop},
  {t01,t12,t20} = hops;
  hop = {{0,-t01,-t20 Exp[I flux]},{-t01,0,-t12},{-t20 Exp[-I flux],-t12,0}};
  lead + KroneckerProduct[hop+DiagonalMatrix[siteEnergy],id] +
    Sum[KroneckerProduct[DiagonalMatrix[UnitVector[3,j]],p[[j]].sigma], {j,3}]];
pins = {{ax,ay,az},{bx,by,bz},{cx,cy,cz}};
real = Element[Join[Flatten[pins],{u,v,w,f,b}],Reals];
h = ham[pins,{u,v,w},f,{0,0,0}];
hf = D[h,f];
zero[m_] := TrueQ[FullSimplify[And@@Thread[Flatten[m]==0],Assumptions->real]];
identities = {
  zero[(h /. f->Pi-f) + transform[h] - IdentityMatrix[6]],
  zero[(hf /. f->Pi-f) - transform[hf]],
  zero[transform[up]-down], zero[transform[lead]+lead-IdentityMatrix[6]],
  zero[{Tr[hf]}]};
emit["Hamiltonian, derivative, beam exchange, tag complement, trace derivative",identities];
site = KroneckerProduct[DiagonalMatrix[{0,0,b}],id];
shiftRemainder = Simplify[Diagonal[site+transform[site]],Assumptions->real];
emit["diagonal residual of the Hamiltonian identity with a static site shift",shiftRemainder];

(* In an eigenbasis the time-integrated density is its initial value times
   exp(-i Delta/2) sinc(Delta/2) entry by entry, for transit time one.
   Use this independent force integral at midpoint samples of the full turn. *)
impulse[p_, detuning_, f_?NumericQ] := Module[
  {m,e,rows,cols,density,derivative,delta,integrated},
  m = N[ham[-2 p,{1,1,1},f,detuning p[[All,3]]^2]];
  {e,rows} = Eigensystem[m];
  cols = Transpose[rows];
  density = ConjugateTranspose[cols].up.cols;
  derivative = ConjugateTranspose[cols].KroneckerProduct[
      {{0,0,-I Exp[I f]},{0,0,0},{I Exp[-I f],0,0}},id].cols;
  delta = Outer[Subtract,e,e]/2;
  integrated = density Exp[-I delta] Sinc[delta];
  -Re[Tr[derivative.integrated]]];
fullTurn[p_, b_, n_] := (2 Pi/n) Total[Table[
    impulse[p,b,N[(2 k+1) Pi/n]], {k,0,n-1}]];
triangle = IdentityMatrix[3];
values[n_, b_] := Table[fullTurn[RotateLeft[triangle,j],b,n],{j,0,2}];
baseline = values[48,0];
split48 = values[48,3/4]; split96 = values[96,3/4];
mirror = Table[fullTurn[RotateLeft[Reverse[triangle],j],3/4,48],{j,0,2}];
emit["zero onsite shift full turns three markings",baseline];
emit["b=3/4 full turns three markings, ordinary mean, mirror mean at 48 midpoints",
    {split48,Mean[split48],Mean[mirror]}];
emit["b=3/4 full turns three markings, ordinary mean at 96 midpoints",
    {split96,Mean[split96]}];
error = Max[Abs[Join[baseline,split48-split96,{Mean[split48]+Mean[mirror]}]]];
emit["maximum of undetuned period, refinement change and mirror-sum residuals",error];
Exit[If[And@@identities && TrueQ[error<10^-10],0,1]];
