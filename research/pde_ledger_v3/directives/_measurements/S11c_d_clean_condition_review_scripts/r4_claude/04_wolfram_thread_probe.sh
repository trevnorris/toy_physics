#!/usr/bin/env bash
# Leg script 04: count the tasks (threads) of a trivial Wolfram kernel, launched with the
# same single-thread environment variables the guard sets, and its peak RSS.
# Purpose: the guarded runner caps TasksMax=32 (pids.max counts threads) and MemoryMax=2 GiB;
# no Wolfram kernel has been launched under it (no wolfram command in _scratch/s11c/*/invocation.json).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
timeout 600 wolframscript -code 'x = Expand[(a + b + c)^40]; Pause[12]; Print[Length[x]]' > /tmp/s11cd_clean_review_r4_claude/04_wolfram_probe_kernel.out 2>&1 &
WS=$!
maxthreads=0; maxtasks=0
for i in $(seq 1 30); do
  sleep 1
  total=0
  desc="$WS"; frontier="$WS"
  while [ -n "$frontier" ]; do nxt=""; for q in $frontier; do c=$(pgrep -P $q | tr '\n' ' '); nxt="$nxt $c"; done; frontier=$(echo $nxt); desc="$desc $frontier"; done
  for p in $desc; do
    t=$(awk '/^Threads:/{print $2}' /proc/$p/status 2>/dev/null)
    [ -n "$t" ] && total=$((total + t))
    rss=$(awk '/^VmRSS:/{print $2}' /proc/$p/status 2>/dev/null)
    echo "SAMPLE $i PID $p CMD $(tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | cut -c1-80) THREADS $t RSS_KB $rss"
  done
  echo "SAMPLE $i TOTAL_THREADS_ALL_WOLFRAM_PROCS $total"
  [ $total -gt $maxtasks ] && maxtasks=$total
  if ! kill -0 $WS 2>/dev/null; then break; fi
done
wait $WS; echo "WOLFRAMSCRIPT_EXIT $?"
echo "MAX_TOTAL_THREADS $maxtasks"
cat /tmp/s11cd_clean_review_r4_claude/04_wolfram_probe_kernel.out
