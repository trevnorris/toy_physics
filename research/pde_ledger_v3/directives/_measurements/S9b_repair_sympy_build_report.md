# S9b SymPy harness repair round 6 — amendment 9

**Complete; stopped for handoff; no commits.** The engine is byte-unchanged. The final harness completed baseline, K1–K10 and compact K11. Source pins, every single-site mutation and full emission coverage checked. The baseline regenerated the delta and is byte-identical to the earlier validated baseline. The completed final-source method demonstration was retained without replay.

**Final guarded outcomes:** both exit **0**, verified `--pool s9b --memory-gib 16`, zero sampled swap, no deadlines. Harness runtime **10659.089266449 s**, sampled cgroup peak **17179869184 bytes (16 GiB)**. Its `memory.events max` reached **1514**; all sampled OOM counters remained **0**. Method runtime **0.8802560490003088 s**, sampled peak **53874688 bytes**. Exact commands (recorded, not rerun):

```sh
R=/var/projects/toy_physics-s9b-r2-py
P="$R/_scratch/s9b-repair-amend9"
G="$R/_scratch/s11c"
cd "$R"
/usr/bin/python -B scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --log-directory "$G/s9b-amend9-final2-method" -- /usr/bin/python -B "$P/final2/method_probe.py"
/usr/bin/python -B scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --log-directory "$G/s9b-amend9-final3-harness" -- /usr/bin/python -B "$R/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py" --scratch "$P/final3/harness"
```

**Repair:** exact rational nonidentity certificates, admissible finite radial jets, linear-real satisfiability and set-membership witnesses establish existence of change. Value comparisons retain common definition domains; domain differences remain separate. Proved equality prints zero; unresolved comparisons print `UNDECIDED`. Failed witness slices do not establish equality.

**Inventory:** **330 shared + 49 local** baseline tags; none missing; **271** contain conditional `NOT_ESTABLISHED`. The [payload-free inventory](S9b_repair_sympy_tag_inventory.json) names every emitted tag, run/source hash and marker, retaining flags, counts, exclusions, commands and outcomes. It includes all undecided tags and exact knife sites: K1 advection, K2 metric, K3 speed, K4a/b/c derivatives, K5 return direction, K6 radar source, K7 quantifier, K8 mass density, K9a/b responses, K10 gate, K11 azimuthal flow.

**Item 13/bindings:** general radial functions; Abel every-impact reduction; coordinate mass divergence and velocity-square branches. Only `c_s0` is bound (spec line **223**); directive item **9** supplies `rho_0 > 0`. No shared profile family is selected.

**K11:** engine constructions, all 12 grades, exact arithmetic; `L=c0=1`, `delta=1/(10r²)`, `V=1/(8r)`, `xi=log(r)/6`, azimuthal coefficient `1/12`, points `(3,4,0),(4,0,3),(0,6,8)`. Three objects evaluated, all method residuals zero in value and domain; **376** tags individually listed as `NOT_EVALUATED`. Unsampled profiles/points are excluded.

**Events/open inputs:** two adapter bugs were fixed before final runs. The earlier stop was orchestrator-requested for reboot; the later interruption was an OpenAI service error after harness completion. No guard refusal/kill or item 12 event in the successful continuation. `Phi` versus `GM` and unsupplied material relations remain open. No Part D, comparator, production run or review was launched.
