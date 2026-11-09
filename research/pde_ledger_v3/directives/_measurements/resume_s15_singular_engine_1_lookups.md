# Lookups for resume note 1 of the S15 item-4 Singular engine build

Mechanical lookups by the orchestrator, 2026-11-09 (UTC). Each is a command and its literal output. The host's local time zone is -0700, so a local file time of 21:14:41 on 11-08 is 04:14:41Z on 11-09.

## L1. The first builder session's end marker and its time

```
$ date +%z; ls -la --time-style=+%FT%T%z _scratch/s15/singular_engine/build.DONE; cat _scratch/s15/singular_engine/build.DONE
-0700
-rw-rw-r-- 1 trevnorris trevnorris 7 2026-11-08T20:45:00-0700 _scratch/s15/singular_engine/build.DONE
exit 0
```

## L2. attempt01's command, pool and cap

```
$ python3 -c 'import json; d=json.load(open("_scratch/s11c/s15_singular/attempt01/invocation.json")); print(d["command"]); print(d["unit"], d["memoryMax"], d["memoryMax"]/2**30, "GiB")'
['/usr/bin/python3', '-B', '/var/projects/toy_physics/_scratch/s15/singular_engine/build01/attempt01_source/S15_separation_departure_item4.py', '--work-directory', '/var/projects/toy_physics/_scratch/s15/singular_engine/build01/attempt01_work']
s11c-guard-fb8466aef0ae 60129542144 56.0 GiB
```

## L3. attempt01's resource samples: first and last time, peak memory, out-of-memory events

```
$ python3 -c '
import json
rows=[json.loads(l) for l in open("_scratch/s11c/s15_singular/attempt01/resource-samples.jsonl") if l.strip()]
print(len(rows), rows[0]["utc"], rows[-1]["utc"])
print(max(int(r["memory.peak"]) for r in rows)/2**30, "GiB peak")
print([l for l in rows[-1]["memory.events"].split(chr(10)) if l.startswith("oom")])'
26708 2026-11-09T03:42:21.473968+00:00 2026-11-09T18:32:56.476047+00:00
11.229839324951172 GiB peak
['oom 0', 'oom_kill 0', 'oom_group_kill 0']
```

## L4. attempt01's output: when it last grew, and its last WORKING_ON lines

```
$ ls -la --time-style=+%FT%T%z _scratch/s11c/s15_singular/attempt01/stdout _scratch/s15/singular_engine/build01/attempt01_work/object_000.stdout; grep -c '^WORKING_ON:' _scratch/s11c/s15_singular/attempt01/stdout; grep '^WORKING_ON:' _scratch/s11c/s15_singular/attempt01/stdout | tail -3
-rw-rw-r-- 1 trevnorris trevnorris 2635203 2026-11-08T21:14:41-0700 _scratch/s11c/s15_singular/attempt01/stdout
-rw-rw-r-- 1 trevnorris trevnorris 2581582 2026-11-08T21:14:41-0700 _scratch/s15/singular_engine/build01/attempt01_work/object_000.stdout
9
WORKING_ON: differential ideal round 1
WORKING_ON: differential ideal round 2
WORKING_ON: differential ideal round 3
```

## L5. Objects in scope, and which objects attempt01 reached

```
$ python3 -c 'import json; print(len(json.load(open("_scratch/s15/singular_engine/build01/attempt01_work/inventory.json"))))'; ls _scratch/s15/singular_engine/build01/attempt01_work
144
inventory.json
object_000.sing
object_000.stdout
```

## L6. attempt01's record files (no child-outcome.json, no outcome.json)

```
$ ls _scratch/s11c/s15_singular/attempt01
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
resource-samples.jsonl
stderr
stdout
```

## L7. The unit after the orchestrator stopped it (stopped with `systemctl --user stop s11c-guard-fb8466aef0ae.service` at 2026-11-09T18:32:56Z, on the user's decision)

```
$ systemctl --user show s11c-guard-fb8466aef0ae.service -p ActiveState -p SubState -p Result -p Job; pgrep -x Singular || echo 'no Singular process'
Result=success
ActiveState=inactive
SubState=dead
Job=
no Singular process
```

## L8. The guard registry, reconciled by the guard's rule (terminal unit, no queued job, dead owner, old registry kept)

```
$ cat _scratch/s15/singular_engine/reconcile_only.log; cat /home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json; ls /home/trevnorris/.local/state/s11c-guard/resource-pool | grep stale-s11c-guard-fb8466aef0ae
PROPS {'ActiveState': 'inactive', 'Job': ''}
REGISTRY_BEFORE {"s11c-guard-fb8466aef0ae.service": {"pool": "s15-singular", "poolMemoryMax": 66571993088, "memoryMax": 60129542144, "tasksMax": 32, "cpu": 7, "owner": {"pid": 3589903, "startTicks": "146589063", "bootId": "c3b7b2c7-a5fc-4b45-adf0-6cbde3d609bb"}, "logDirectory": "/var/projects/toy_physics/_scratch/s11c/s15_singular/attempt01"}}
OWNER {'pid': 3589903, 'startTicks': '146589063', 'bootId': 'c3b7b2c7-a5fc-4b45-adf0-6cbde3d609bb'} ALIVE False
RECONCILED; old registry kept at /home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json.stale-s11c-guard-fb8466aef0ae.service-20261109T183307Z
REGISTRY_AFTER {}
{}
reservations.json.stale-s11c-guard-fb8466aef0ae.service-20261109T183307Z
```

