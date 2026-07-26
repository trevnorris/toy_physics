# RESTART PROMPT — paste after compact / at session start

> ⭐ **Read `research/pde_ledger_v2/manifests/REWRITE_HANDOFF.md` FIRST**, then
> `notes/rewrite_reference_table.md` (per-stage data). Background = `manifests/PIVOT_TO_REWRITE.md`.
> ⛔ Do NOT read survey-era docs (`HANDOFF_NEXT_SESSION.md`, `SURVEY_DIRECTIVE.md`, `PILOT_DIRECTIVE.md`)
> — that approach is dead. Do NOT read vision docs.
>
> Branch `ledger-v2-rebuild`. **We are rewriting 30 audit scripts onto a shared dimension module,
> one stage at a time. 2 are done (`adcfbdfd`), 28 remain.**
>
> **The per-stage loop** (REWRITE_HANDOFF §2): add print-only dimension output to the stage's `.wl` →
> re-run and re-baseline its `.out` → **write the prediction down** → rewrite the `.py` onto
> `scripts/ledger_dimensions.py` → compare `.py` vs `.wl` **axis-labelled** → commit.
> Acceptance: `.py` exits 0 with an identical PASS multiset; `.wl` exits 0 with an unchanged tally and
> reproduces its `.out` byte-identically after `$NNNNN` normalisation; zero unpredicted mismatches.
>
> ⭐ **Why step (e) is mandatory:** transposing stage011's basis left the script **passing all 60
> markers, exit 0** — a script's own assertions cannot see a transposition once both ends come from the
> module. The labelled cross-engine comparison caught 7 of 10 values. It is the only thing between us
> and a silent corpus-wide transposition.
>
> **FIRST TWO ACTIONS, in order:**
> 1. **Settle the open decision in REWRITE_HANDOFF §7** — write stage038's and stage042's basis
>    declarations on paper and confirm nothing in the module structurally blocks them. Short check, not
>    a redesign. Cheap now; a re-run of groups A–C if discovered late.
> 2. **Settle whether the `.py` should emit exponent triples** (§2 note). The pilot's comparison came
>    from an **external in-process script**, not a transcript diff — so the cross-check is not currently
>    reproducible from committed artifacts. Decide before scaling.
>
> Then **group A** in order: 012, 013, 018, 016, 023, 027, 021.
> ⭐ Group A first because the 9 stages whose `.out` already renders computed values are exactly the
> `(L,M,T)` group plus 004 — no `.wl` edit and no Mathematica seat needed. 012/013 reuse stage011's
> scaffolding almost verbatim.
>
> ⚠ **Adjudication lands at group A step 2**: `K_eta` carries three different dimensions across
> 013/016/023. **Stop and bring it to the user** — it is a question about the model, not the code. There
> are **13** such conflicts, not the 9 listed in PIVOT §3.
>
> **Gotchas:** `grep -c '^PASS'` over-counts by exactly 1 (the tally line self-matches) · Codex needs
> `--sandbox danger-full-access` for Mathematica · **never >2 concurrent `math -script`** · never wait
> on a `pgrep` pattern your own waiter contains (cost 6h43m).
>
> **Landmines** (REWRITE_HANDOFF §5): stage003's `.wl` is `(M,L,T)` and says so nowhere · stage042's
> `.wl` comment mislabels its basis · stage021 emits three renderings under two conventions · **five
> stages where the `.wl` print step is impossible**: 037, 036, 035, 044, 042.
>
> **Operating model:** thin conductor, agents return ≤12-line receipts, `model: opus` on every dispatch,
> gate every Codex launch on a check that precedes it. **Do not build a corpus-wide inventory, oracle,
> or completeness proof** — three such gates were specified and rejected; the per-stage loop replaces
> them. **Do not generalise the module beyond what the stage in front of you needs.**

---

## State

| Item | State |
|---|---|
| Module | ✅ `scripts/ledger_dimensions.py`, 179 lines, axis-labelled + per-stage basis |
| Stages done | **2 / 30** — stage004 `(L,T,M)`, stage011 `(L,M,T)` (`adcfbdfd`) |
| Remaining | **28** (dim lines 3,752 → 3,675) |
| Cross-check | ✅ proven able-to-fail (7/10 caught on a deliberate transposition) |
| Open decisions | §7 four-axis/stiffness paper check · `.py` triple emission |
| Conflicts | **13**, first at group A step 2 (`K_eta`, 013/016/023) |

## Do not redo
`.out` determinism (44/44) · SymPy baseline (43/43 exit 0) · the two measurements
(`notes/measure_*.md`) · the reference table · repair round 8 (parked, survey-era).

## Parked, only if the survey is revived
`schemas/` + validator (`35dc6aa0`, 173 fixtures green), `SURVEY_DIRECTIVE.md` r18, tasks #12/#13/#16–#19.
