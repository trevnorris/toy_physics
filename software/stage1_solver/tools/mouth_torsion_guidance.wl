(* Flat welded slab in unbounded lossless bulk. One field u(x,w), per face
   and polarization. Both x and w gradients load equal, independent contact
   packets. This isotropic routing is an extra material hypothesis. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[N[value, 15], InputForm]];
stiff = 2; rhoS = 4; refH = 8192; halfWidth = 32;
w0 = Sqrt[stiff/refH];
mod[om_] := stiff (1 - stiff/(stiff - refH om^2));
ks2[om_] := rhoS om^2/stiff;
kb2[rho_, om_] := rho om^2/mod[om];
omega = 2 w0;
rhoB = 4;

(* Propagate the even core solution from its center, then subtract the
   traction on an exterior decaying exponential at the upper face. *)
mismatch[k_, om_] := With[{p = Sqrt[ks2[om] - k^2],
                           beta = Sqrt[k^2 - kb2[rhoB, om]]},
  stiff p Sin[halfWidth p] - mod[om] beta Cos[halfWidth p]];
delta = ks2[omega] - kb2[rhoB, omega];
(* x = p/sqrt(delta) on the first, pole-free half lobe. *)
res[x_] := stiff x Sin[halfWidth Sqrt[delta] x] -
  mod[omega] Sqrt[1 - x^2] Cos[halfWidth Sqrt[delta] x];
lo = 0; hi = 1;
Do[mid = (lo + hi)/2;
   If[N[res[mid], 45] > 0, hi = mid, lo = mid], {100}];
p = N[Sqrt[delta] (lo + hi)/2, 35];
kk = Sqrt[ks2[omega] - p^2];
beta = Sqrt[kk^2 - kb2[rhoB, omega]];
faceAmp = Cos[halfWidth p];
emit["{L,omega/omega0,rhoS,rhoB,k_bulk,k_parallel,k_core,p,beta}", {
  2 halfWidth, omega/w0, rhoS, rhoB, Sqrt[kb2[rhoB, omega]], kk,
  Sqrt[ks2[omega]], p, beta}];
emit["{face_amplitude,tail_amplitude_decay_length,traction_residual}", {
  faceAmp, 1/beta, mismatch[kk, omega]}];

(* Cycle-mean power equals minus real torque times conjugate velocity / 2.
   Core norm is integrated over -a..a. Tail norm includes both exterior
   half-spaces. The energy of EACH free reference and its spring enters
   both x and w strain components; center displacement amplitude is one. *)
coreNorm = halfWidth + Sin[2 halfWidth p]/(2 p);
coreGradNorm = p^2 (halfWidth - Sin[2 halfWidth p]/(2 p));
tailNorm = faceAmp^2/beta;
slip = stiff/(stiff - refH omega^2);
strainEnergyCoef = stiff (1 - slip)^2 + refH omega^2 slip^2;
coreE = (rhoS omega^2 coreNorm + stiff (kk^2 coreNorm + coreGradNorm))/4;
tailE = (rhoB omega^2 + strainEnergyCoef (kk^2 + beta^2)) tailNorm/4;
coreP = omega kk stiff coreNorm/2;
tailP = omega kk mod[omega] tailNorm/2;
work[traction_, velocity_] := -Re[traction Conjugate[velocity]]/2;
outward = work[-mod[omega] beta faceAmp, -I omega faceAmp];
transport = (coreP + tailP)/(coreE + tailE);
slope = (-D[mismatch[kvar, omvar], kvar]/D[mismatch[kvar, omvar], omvar]) /.
  {kvar -> kk, omvar -> omega};
emit["{outward_face_power,core_axial_power,tail_axial_power,core_energy,tail_energy}", {
  outward, coreP, tailP, coreE, tailE}];
emit["{tail_fraction_axial_power,tail_fraction_energy,work_per_energy,dispersion_slope}", {
  tailP/(coreP + tailP), tailE/(coreE + tailE), transport, slope}];

(* A prescribed incident core plane wave supplies a positive normal wave
   number. Use traction continuity to solve the single welded face for its
   reflected and outgoing waves. Fluxes here refer to that incident wave. *)
scatter[rho_, om_, kx_] := Module[{qS, qB, zs, zb, amp, ampR, refl, enter},
  qS = Sqrt[ks2[om] - kx^2];
  qB = Sqrt[kb2[rho, om] - kx^2];
  zs = stiff qS; zb = mod[om] qB;
  amp = 2 zs/(zs + zb); ampR = amp - 1;
  refl = Abs[ampR]^2;
  enter = work[I qB mod[om] amp, -I om amp]/work[I qS stiff, -I om];
  emit["{rhoB,omega/omega0,k_parallel,Re(q_bulk),Im(q_bulk),R,entry,R+entry}", {
    rho, om/w0, kx, Re[qB], Im[qB], refl, enter, refl + enter}];
  {refl, enter}];
normal = scatter[rhoB, omega, 0];
guided = scatter[rhoB, omega, kk];
heavier = scatter[8, omega, kk];
subpole = scatter[rhoB, w0/2, 0];
emit["{heavy_k_bulk,heavy_k_core,rhoB_for_equal_cones_at_drive,heavy_cutoff_over_omega0}",
  {Sqrt[kb2[8, omega]], Sqrt[ks2[omega]],
   rhoS omega^2/(omega^2 - w0^2), Sqrt[8/(8 - rhoS)]}];
checks = {halfWidth Sqrt[delta] < Pi/2, res[0] < 0, res[1] > 0,
  Sqrt[kb2[rhoB, omega]] < kk < Sqrt[ks2[omega]],
  Abs[mismatch[kk, omega]] < 10^-25, coreE > 0, tailE > 0,
  Abs[transport - slope] < 10^-24, Abs[outward] < 10^-25,
  normal[[2]] > 0, Abs[guided[[2]]] < 10^-25, heavier[[2]] > 0,
  kb2[8, omega] > ks2[omega], Abs[subpole[[2]]] < 10^-25,
  And @@ (Abs[Total[#] - 1] < 10^-24 & /@ {normal, guided, heavier, subpole})};
Print["checks = ", ToString[checks, InputForm]];
Exit[If[And @@ checks, 0, 1]];
