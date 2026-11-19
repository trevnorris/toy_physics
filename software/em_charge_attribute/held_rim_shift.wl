(* Obtain the transmitted load from the x-momentum conormal across a depth
   cut and integrate its constant flux. The tensor term is kappa chi_x chi_w;
   integrate in the co-translating section to remove the odd angular term. *)
ClearAll["Global`*"];
emit[key_, v_] := Print[key," = ",ToString[v,InputForm]];
conormal=FullSimplify[2 Pi Integrate[n cr (cw-n cr slope),{n,-1,1}]];
emit["integrated x-w conormal per kappa r^2 dr",conormal];
ep=1/50; ell=1/200; len=37/20; kval=1/2; amp=1/10; delta=1/500;
rad[s_]:=1+(3s^2-2s^3)/2;
(* In order-fraction coordinates a logistic wall gives moment generator
   int u^(1+t) (1-u)^(1-t) du = Gamma[2+t] Gamma[2-t]/6 on the full line.
   Bound the extension past the physical origin separately. *)
gen=Gamma[2+t] Gamma[2-t]/6;
mom=FullSimplify[Table[D[gen,{t,j}]/.t->0,{j,0,2}]];
emit["logistic radial moments",mom];
stiff[r_]:=4 Pi ep^2/(3 ell) ({r^2,2 r ell,ell^2}.mom);
tail=4 Pi/3 ep^2 ell Exp[-2/ell]/4;
emit["upper bound for stiffness extension per kappa",N[tail,9]];
compliance=len NIntegrate[1/stiff[rad[s]],{s,0,1},WorkingPrecision->35];
kk=1/compliance;
emit["flare and cylinder stiffness per kappa",N[{kk,stiff[1]/len},16]];
emit["signed holds per kappa",N[-kk {delta,-delta,delta/2,-delta/2,0},16]];
(* With a possible local compliance, compose displacement jumps at the same
   constant transmitted load. Both endpoint laws must use that one load. *)
transmitted=force/.First[Solve[y-x==force/ck+force/lk,force]];
emit["external hold divided by far shift with x=0",Factor[-transmitted/y/.x->0]];
emit["free contact and pinned limits",FullSimplify[
    {Limit[-transmitted,ck->0],Limit[-transmitted,ck->Infinity]}/.x->0]];

(* Integrate the density deficit by decomposing it into sharp spheres with
   weights -d deficit/dr. Disk slicing gives the sharp sphere cosine moment. *)
ball=FullSimplify[Integrate[Pi(r^2-u^2) Cos[2 k u],{u,-r,r}],
   Assumptions->r>0&&k>0];
emit["cosine moment of a ball",ball];
density[s_]:=ep/(4 Cosh[s/2]^2);
form=NIntegrate[Evaluate[(ball/.{k->kval,r->1+ell s}) density[s]],
    {s,-1/ell,-10,0,10,Infinity},WorkingPrecision->34];
emit["midplane order form factor",form];
(* Fixed-action first variation (constant transverse inertia), or virtual
   work on a translated fixed shape. Additive uniform terms drop out. *)
epart=bpeak ff Cos[2 k (y-x)+2 phi]/2;
hold=D[epart,x]/.x->0;
odd=FullSimplify[((hold/.y->d)-(hold/.y->-d))/(2 d),
    Assumptions->Element[{d,phi,k},Reals]&&d!=0];
emit["signed holder slope at finite shift",odd];
emit["lab-held fringe direct mixed entry",D[epart/.x->0,x,y]];
rules={k->kval,ff->form,bpeak->amp^2 kval^2/4};
Do[emit["phase/Pi="<>ToString[ph/Pi,InputForm]<>" baseline and signed slopes",
   N[{hold/.y->0,odd/.d->delta,odd/.d->delta/2}/.phi->ph/.rules,16]],
   {ph,{0,Pi/4,-Pi/4,Pi/2}}];
If[!TrueQ[kk>0&&form>0] || !FreeQ[{kk,form},Indeterminate],Exit[1]];
