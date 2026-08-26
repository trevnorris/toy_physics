(* Isolated background-current form ablation.
   On a COPY of the engine defs, replace bulkVelocityZero with a formal nonzero probe.
   Faithful j0 = rho0 * v0 => current tracks probe.
   Disguised :=0 => current stays 0 under probe. *)

(* Engine rhoBulkZero form (symbolic background density — nonzero free head) *)
rhoBulkZero[coordinates_List, normal_] :=
  rhoBulkBackground @@ Append[coordinates, {normal, time}];

(* PROBE: replace bulkVelocityZero with formal nonzero 4-vector *)
bulkVelocityZeroPROBE[coordinates_List, normalCoordinate_] :=
  {vp1, vp2, vp3, vpW};

(* Engine bulkCurrentZero / currentWZero / currentXZero forms, using PROBE velocity *)
bulkCurrentZero[coordinates_List, normal_] :=
  rhoBulkZero[coordinates, normal] bulkVelocityZeroPROBE[coordinates, normal];
currentWZero[coordinates_List, normal_] :=
  Last[bulkCurrentZero[coordinates, normal]];
currentXZero[index_Integer][coordinates_List, normal_] :=
  bulkCurrentZero[coordinates, normal][[index]];

coords = {x1, x2, x3};
n = w;

Print["PROBE_bulkVelocityZero = ", bulkVelocityZeroPROBE[coords, n]];
Print["PROBE_rhoBulkZero = ", rhoBulkZero[coords, n]];
Print["PROBE_bulkCurrentZero = ", bulkCurrentZero[coords, n]];
Print["PROBE_currentWZero = ", currentWZero[coords, n]];
Print["PROBE_currentXZero[1] = ", currentXZero[1][coords, n]];
Print["PROBE_currentXZero[2] = ", currentXZero[2][coords, n]];
Print["PROBE_currentXZero[3] = ", currentXZero[3][coords, n]];

expectedW = rhoBulkZero[coords, n] * vpW;
expectedX1 = rhoBulkZero[coords, n] * vp1;
expectedX2 = rhoBulkZero[coords, n] * vp2;
expectedX3 = rhoBulkZero[coords, n] * vp3;
Print["TRACKS_probe_W = ", Simplify[currentWZero[coords, n] - expectedW] === 0];
Print["TRACKS_probe_X1 = ", Simplify[currentXZero[1][coords, n] - expectedX1] === 0];
Print["TRACKS_probe_X2 = ", Simplify[currentXZero[2][coords, n] - expectedX2] === 0];
Print["TRACKS_probe_X3 = ", Simplify[currentXZero[3][coords, n] - expectedX3] === 0];
Print["STAYS_IDENTICALLY_ZERO_UNDER_PROBE = ",
  And[
    currentWZero[coords, n] === 0,
    currentXZero[1][coords, n] === 0,
    currentXZero[2][coords, n] === 0,
    currentXZero[3][coords, n] === 0
  ]];
