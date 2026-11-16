(* The isolated conservative carrier from conversion_region_weak_flow.wl.
   Currents are per unit transverse three-area. Compare fixed cuts. *)
ClearAll["Global`*"];
hb = 1/Sqrt[6];
psi = Exp[I (w - mu t)/hb];
real = Element[{w, t, mu, b, len}, Reals] && len > 0;
dens = FullSimplify[Conjugate[psi] psi, real];
flux = FullSimplify[hb Im[Conjugate[psi] D[psi, w]], real];
muRule = First[Solve[FullSimplify[
    (I hb D[psi, t] + hb^2 D[psi, {w, 2}]/2 - 35 dens^4 psi/96)/psi, real] == 0, mu]];
Print["chemical potential and [density, flux] = ", InputForm[{muRule, dens, flux}]];
stored = Integrate[dens, {w, 0, len}];
numberBalance = (flux /. w -> 0) - (flux /. w -> len) - D[stored, t];
(* The specified passive loss law with the previous zero-drive offset.
   Integrate the source in physical w independently of the two traces. *)
fraction = (1 - Tanh[(w - b)/2])/2;
loss = Integrate[dens fraction (1 - fraction) // TrigReduce,
    {w, 0, len}, Assumptions -> real];
channels = {fraction, 1 - fraction} flux;
intake = channels /. w -> 0;
outlet = channels /. w -> len;
checks = FullSimplify[
    {numberBalance, Sequence @@ (outlet - intake + {loss, -loss})}, real];
Print["inventory per area and section fluxes = ", InputForm[{stored, flux /. w -> 0, flux /. w -> len}]];
Print["integrated net order loss = ", InputForm[loss]];
Print["constituent, ordered, disordered residuals = ", InputForm[checks]];
sample = FullSimplify[Join[intake, outlet, {loss}] /. {b -> Log[3], len -> 2 Log[3]} // TrigToExp];
Print["[ordered in, disordered in, ordered out, disordered out, converted] = ", InputForm[sample]];
If[checks =!= {0, 0, 0}, Exit[1]];
