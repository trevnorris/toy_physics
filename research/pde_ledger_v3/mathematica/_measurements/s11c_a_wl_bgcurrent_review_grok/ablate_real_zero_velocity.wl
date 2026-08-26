(* Reverse check: real engine bulkVelocityZero (identically zero 4-vector).
   Background current must evaluate to 0; no Part::pkspec1 on concrete OR symbolic index. *)

rhoBulkZero[coordinates_List, normal_] :=
  rhoBulkBackground @@ Append[coordinates, {normal, time}];

(* Exact engine bulkVelocityZero form from lines 117-120 *)
bulkVelocityZero[coordinates_List, normalCoordinate_] := Through[
  {zeroVelocityOne, zeroVelocityTwo, zeroVelocityThree, zeroVelocityW}[
      Sequence @@ Append[coordinates, {normalCoordinate, time}]]] /.
    {zeroVelocityOne[___] -> 0, zeroVelocityTwo[___] -> 0,
      zeroVelocityThree[___] -> 0, zeroVelocityW[___] -> 0};

bulkCurrentZero[coordinates_List, normal_] :=
  rhoBulkZero[coordinates, normal] bulkVelocityZero[coordinates, normal];
currentWZero[coordinates_List, normal_] :=
  Last[bulkCurrentZero[coordinates, normal]];
currentXZero[index_Integer][coordinates_List, normal_] :=
  bulkCurrentZero[coordinates, normal][[index]];

coords = {x1, x2, x3};
n = w;

msgs = {};
Internal`InheritedBlock[{Message},
  Unprotect[Message];
  Message[args___] := (AppendTo[msgs, HoldForm[Message[args]]];
    If[FreeQ[Hold[args], Part::pkspec1], Message[args], Null]);
  Protect[Message];

  Print["REAL_bulkVelocityZero = ", bulkVelocityZero[coords, n]];
  Print["REAL_rhoBulkZero = ", rhoBulkZero[coords, n]];
  Print["REAL_bulkCurrentZero = ", bulkCurrentZero[coords, n]];
  Print["REAL_currentWZero = ", currentWZero[coords, n]];
  Print["REAL_currentXZero[1] = ", currentXZero[1][coords, n]];
  Print["REAL_currentXZero[2] = ", currentXZero[2][coords, n]];
  Print["REAL_currentXZero[3] = ", currentXZero[3][coords, n]];

  (* symbolic-index call: pattern is index_Integer, so non-integer should stay inert *)
  inert = currentXZero[k][coords, n];
  Print["SYMBOLIC_INDEX_currentXZero[k] = ", inert];
  Print["SYMBOLIC_INDEX_stays_inert = ", MatchQ[inert, currentXZero[k][coords, n]] || Head[inert] === currentXZero];
];

pkspecMsgs = Select[msgs, !FreeQ[#, Part::pkspec1] &];
Print["Part_pkspec1_message_count = ", Length[pkspecMsgs]];
Print["Part_pkspec1_messages = ", pkspecMsgs];
Print["currentWZero_is_zero = ", currentWZero[coords, n] === 0];
Print["currentXZero_all_zero = ",
  And[
    currentXZero[1][coords, n] === 0,
    currentXZero[2][coords, n] === 0,
    currentXZero[3][coords, n] === 0
  ]];
Print["rhoBulkZero_still_nonzero_head = ", Head[rhoBulkZero[coords, n]] === rhoBulkBackground];
