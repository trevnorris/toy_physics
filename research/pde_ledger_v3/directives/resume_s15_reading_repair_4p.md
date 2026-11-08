Note from the orchestrator. Review round 4 of your repair is done. Both legs found the same one defect, and I confirmed it. The reviewed version is committed at `0d9675eb`. Its two transcripts are kept as `research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r4.out` and `research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r4.out`. Do not change those two copies.

## The finding

The finding concerns the six region objects that carry both drains, `region/component_{1,2,3}/{momentum,rate}_labels_{drain1, drain2}`.

For each of them, your final Wolfram run started an item-4 kernel. That kernel computed the sufficient conditions and `Q` under each, and printed them before the third-kind label. Review found these conditions true.

The label and the run record describe something else. The label says: "The retained sufficient conditions are printed above." The run record says the six "were not attempted in these final segments; each retained its previously printed sufficient conditions and corresponding Q objects."

The cause is in my notes, which conflict:
- `resume_s15_reading_repair_4.md` says: "Every condition printed now that is true stays printed."
- `resume_s15_reading_repair_4n.md` says the six "are not attempted again in this build." It also asks for a label that says "the object was not attempted in that run."

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4o.md`, and in the notes it keeps standing, stands, with these changes for the six objects in Wolfram:
- The sufficient conditions, and `Q` under each, stay printed.
- The six objects' lines in the final Wolfram transcript, and the run record, say what the final run did for each of the six. The final run computed the printed conditions and `Q` under each. It made no exponent-sensitive attempt. It did not obtain the complete conditions.
- The label is the third kind, set by the user's decision.
- Nothing says that the object was not attempted in that run, or that its conditions come from an earlier run.
- These rules replace, for the six, the requirement in `resume_s15_reading_repair_4n.md` that the label say the object was not attempted in that run.

Two rules cover every object:
- Every other printed object, in both transcripts, stays as it is by value.
- The run record and the two run manifests say how each final transcript was produced. That covers any part rerun for this repair: which runs made it, from which source and runner (sha256), and how the parts were joined.

## What to do now, and where to stop

- Make the change, and run what it requires through the guard.
- Update the run record and the run manifests.
- Write a short final report. Stop there; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
