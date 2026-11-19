(* Polarization of a real standing pattern at fixed energy. *)
ClearAll["Global`*"];
$Assumptions = rho > 0 && len > 0 && om > 0 && k > 0 && en > 0 &&
  Element[{t, z, ang, x, y, px, py, g, probe, s, h1, h2}, Reals];
emit[label_, expr_, target_] := Module[{v = FullSimplify[expr]},
  Print[label, " = ", ToString[v, InputForm]];
  If[!TrueQ[FullSimplify[v == target]], Exit[1]]];

shape = Sin[Pi z/(2 len)];
norm = Integrate[rho shape^2, {z, 0, len}, Assumptions -> len > 0 && rho > 0];
profile = shape Cos[om t];
emit["mean axial wave energy flux up to a common stiffness",
  Integrate[-D[profile,t] D[profile,z],{t,0,2 Pi/om}]/(2 Pi/om),0];
emit["identical reflected integrated norm",
  Integrate[rho (shape /. z -> -z)^2,{z,-len,0}] - norm,0];

(* The scalar Hamiltonian depends on geometry and amplitude. Test the
   electric source with the matched constituent contribution removed. *)
ham = (px^2 + py^2)/(2 norm) + norm om^2 (x^2+y^2)/2;
emit["native wave source above common carrier charge",
  D[ham + probe q number,probe] - q number, 0];
poisson[a_,b_] := D[a,x] D[b,px]-D[a,px] D[b,x]+D[a,y] D[b,py]-D[a,py] D[b,y];
rot = x py - y px;
emit["free generator time derivative",poisson[rot,ham],0];
emit["source under a test coupling to the rotation",D[ham+probe g rot,probe],g rot];
emit["torque from lifting the degeneracy",poisson[rot,ham + (h1 x^2+h2 y^2)/2],(h1-h2) x y];

(* Continuous ellipse family. ang=0 is linear, +/-Pi/4 are circular.
   Hold the integrated energy fixed, allowing a finite-energy linear state. *)
amp = Sqrt[2 en/(norm om^2)];
coords = amp {Cos[ang] Cos[om t], Sin[ang] Sin[om t]};
state = Thread[{x,y,px,py} -> Join[coords,norm D[coords,t]]];
emit["energy throughout ellipse family", ham /. state,en];
emit["ellipse polarization fraction", om (rot /. state)/en,Sin[2 ang]];
emit["fractions at three preparations",Table[(om (rot /. state)/en) /. ang -> v,
  {v, {-Pi/4,0,Pi/4}}],{-1,0,1}];

(* Derive the two-radius control directly from coordinates and momenta.
   The actions only set the requested initial displacement and velocity. *)
freq[r_] := Sqrt[1 + 4/r^2];
Do[Module[{r=sample[[1]],p=sample[[2]],m=sample[[3]],iv,readEnergy,readRot,load},
  iv={x->(Sqrt[p]+Sqrt[m])/Sqrt[norm om],y->0,px->0,
      py->Sqrt[norm om] (Sqrt[p]-Sqrt[m])};
  readEnergy=FullSimplify[(ham/.iv) /. om->freq[r]];
  readRot=FullSimplify[rot/.iv];
  load=FullSimplify[-(readEnergy/freq[r]) (D[freq[rr],rr] /. rr->r)];
  Print["R, actions, energy, radial load, polarization, native charge = ",
        ToString[{r,p,m,readEnergy,load,readRot,D[ham,probe]},InputForm]]],
 {sample,{{1,1,0},{2,2,1},{2,1,2},{2,3,0}}}];
