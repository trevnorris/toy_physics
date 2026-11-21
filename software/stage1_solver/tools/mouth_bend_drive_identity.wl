(* The held neutral carrier, arbitrary interior evolution and common-frequency baths. *)
ClearAll["Global`*"];
complexQ[z_] := ComplexExpand[z];
field = aa + I bb; fluxGradient = ca + I cb;
number = complexQ[Im[Conjugate[field] fluxGradient]];
energy[rate_] := complexQ[-Re[Conjugate[-I rate field] fluxGradient]];
gauge = Simplify[energy[mu] - mu number];
(* jN into the mouth at feed; jN out at bulk. Difference is stored number. *)
boundaryPower = muF qF - muB qB;
grandPort = Simplify[(boundaryPower - mu0 (qF - qB)) /. {muF->mu0, muB->mu0}];
biasPower = Simplify[(boundaryPower /. {qF->q, qB->q}) - q (muF-muB)];
(* At a held reflector, grad_n f=0 and t·n=0. *)
wallShear = Simplify[complexQ[Re[Conjugate[tx+I ty] fluxGradient]] /. {ca->0, cb->0}];
(* Transport of linear mean contact energy under constant held slip. *)
r[z_] = s0 (1-Exp[-z/(v tau)]);
e[z_] = k (s0-r[z])^2/2;
loss[z_] = k (s0-r[z])^2/tau;
carried = Simplify[v D[e[z],z]+loss[z],Assumptions->(v>0 && tau>0)];
Print["common_frequency_energy_minus_mu_number=",gauge];
Print["equal_frequency_net_grand_energy_through_baths=",grandPort];
Print["two_head_steady_power_residual=",biasPower];
Print["held_native_wall_tangential_traction=",wallShear];
Print["loaded_contact_stream_energy_flux_plus_loss=",carried];
If[{gauge,grandPort,biasPower,wallShear,carried}==={0,0,0,0,0},Exit[0],Exit[1]];
