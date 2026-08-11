# Independent review — the SUBTRACTIVE `C17`/`C18` revision

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree).
⭐ Diff against its base:
```
git -C /var/projects/toy_physics diff 9229bd65 -- research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
```

⚠ **Shared spec both engines read.** An error makes both engines agree on the same wrong thing, or
**disagree for a reason that is not physics**.

⚠⚠ **THIS FILE HAS NOW HAD THREE REVISIONS AND EVERY ONE BRED DEFECTS** — thirteen in total, and **every
one of them lived in a mechanism invented to stop two engines describing a locus differently**: clear a
denominator, pin an elimination, define a donor scope, spell a token set, add an emission axis.

⭐⭐ **This revision is therefore SUBTRACTIVE. It REMOVES mechanism.** ⛔ Your job is to find what the
removals broke, not only what the additions got wrong.

You are one of two independent legs; the other is not visible to you.

## ⛔ DO NOT READ
Any `…_fix_round1_*` or `…_round2_review_prompt` file; anything under `_scratch/`.

## ⭐⭐ READING ORDER — ⛔ artifact last
1. `DEFECT_REGISTER.md`, `C17` and `C18` in full.
2. `directives/S11_C17_C18_spec_repair_decisions_v2.md` — decisions `T1`–`T7`.
3. The committed outputs `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` and
   `scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⭐ **Open them.**
4. **Write down** how *you* would make two blind engines produce comparable stratum information. Keep it.
5. **Only now** the artifact and the diff.

## ⭐⭐⭐ WHAT WAS REMOVED, AND WHAT REPLACED IT — ⛔ judge each

| | removed | replaced by |
|---|---|---|
| **A** | the supplied **locus selector**, the **donor scope**, and the orchestrator's **root-index translation** | ⭐ `WITNESS<w>_OWN_LOCUS_RESIDUALS` — the receiver **computes** which of its own loci the point satisfies, by evaluating **every** one of them at that point |
| **B** | the orchestrator's **pairing test** before supplying a point | ⭐ **every** exact point either engine emitted is supplied to the other, unconditionally |
| **C** | the **pinned elimination** rule | ⭐ symbolic component payloads declared **inspection-only**; **counts and statuses** are the comparison rows |
| **D** | the second count payload shape (`NOT_COMPUTED_COMPONENT_WIDE`) | ⭐ one record in all statuses |

## ⭐⭐⭐ ATTACK THE REMOVALS

1. ⛔⛔ **Is `A` computable and bounded?** The receiver must evaluate **every locus it emitted at that
   package and dimension** at the received point. ⭐ **Count them from the committed outputs** — how many
   loci per `(package, D)`, and how many points would be supplied under `B`? ⛔ Is the resulting object
   size reasonable, or has a combinatorial explosion been specified? ⭐ Give the numbers.
2. ⛔⛔ **Does `A` actually replace the selector?** ⭐ Trace `XFORM_EXTRA, D = 2`: supply the counterpart
   point, evaluate every one of the receiver's loci at it, and report what comes out. ⛔ Can an orchestrator
   tell from that object which locus the donor meant? ⚠ If several of the receiver's loci vanish at the
   point, is that ambiguity or information?
3. ⛔ **Does `B` create work that is not physics?** Every point to every engine — ⭐ is any of it
   meaningless (a point from a locus the receiver has no analogue of)? ⛔ Is meaningless different from
   harmful?
4. ⛔⛔ **Did `C` give up something that was comparable?** ⭐ The claim is that **counts and statuses** are
   invariant under which variable was eliminated, but certificates and change loci are not. ⛔ **Test the
   invariance**: take a real component from the committed outputs, restrict under two different
   eliminations in both CASes, and report whether every count agrees. ⚠ Is any physics-bearing object now
   uncompared that could have been compared?
5. ⛔ **Did `D`'s single record lose the distinction** between *could not compute* and *not defined there*?

## What else to check

### 1. ⛔⛔ Does the file state or imply an answer?
⭐ Clause by clause on added lines. **Derivability, ⛔ not literal presence.**

### 2. ⭐⭐ Is `C18`'s measured case now closed — and can silence still be reached?
⛔ Trace `XFORM_EXTRA, D = 2` end to end from the committed payloads. ⭐ Is the outcome a **computation**,
a **finding**, or **silence**? ⛔ Construct a compliant run that reaches silence, or state that you cannot.

### 3. ⭐⭐ Is the blind build intact?
⚠ The input is now **one field**. ⛔ Does anything still let an engine receive a computed result? ⭐ Is the
Wolfram engine still importing nothing of the other engine's computation?

### 4. ⚠ Stale prose above the fix — ⛔ THREE FOR THREE SO FAR
⭐ For every edited section read its introduction, and search the whole file for sentences that referred to
the removed mechanisms. ⛔ Name any that survive.

### 5. ⭐ What of `T1`–`T7` is still not discharged?

## Physics filter
Report a finding only if it catches a way the **physics could be wrong**, a **claim unsupported**, or
**wrong or incomparable artifacts**. ⛔ Not style.

## Method
- ⭐⭐ Quote both sides of every finding.
- ⭐ State what you opened or ran for every claim about a file, line or count.
- ⭐⭐ Attacks 1, 2 and 4, and check 2, are settled by **running it** on real committed data. ⛔ Prose is
  discarded there.
- ⭐ Any CAS run: `timeout 600`, ⛔ one kernel at a time (two seats), ⛔ copy to `/tmp` and run the copy.
- ⛔ Read-only.
- ⭐ End with: A–D held or not; whether any removal broke something; whether silence is still reachable;
  whether the file leaks; what `T1`–`T7` still lacks; and what you checked that could have failed and did
  not.
