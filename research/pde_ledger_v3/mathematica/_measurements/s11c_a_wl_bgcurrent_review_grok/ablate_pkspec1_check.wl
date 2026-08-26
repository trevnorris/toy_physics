rhoBulkZero[coordinates_List, normal_] :=
  rhoBulkBackground @@ Append[coordinates, {normal, time}];
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

coords = {x1, x2, x3}; n = w;
$MessageList = {};
Off[General::stop];
results = {
  currentWZero[coords, n],
  currentXZero[1][coords, n],
  currentXZero[2][coords, n],
  currentXZero[3][coords, n],
  currentXZero[k][coords, n]
};
msgs = $MessageList;
pks = Cases[msgs, HoldForm[Message[Part::pkspec1, ___]] | HoldForm[Part::pkspec1] | _?(Not@*FreeQ[Part::pkspec1])];
(* Also catch via string form *)
Print["RESULTS = ", results];
Print["MESSAGE_LIST = ", msgs];
Print["HAS_Part_pkspec1 = ", !FreeQ[msgs, Part::pkspec1]];
Print["currentWZero = ", results[[1]]];
Print["currentX concrete = ", results[[2;;4]]];
Print["symbolic index form = ", results[[5]]];
Print["symbolic stays unevaluated = ", MatchQ[results[[5]], _currentXZero[__] | currentXZero[_][__]] ||
  (Head[results[[5]]] === currentXZero)];
