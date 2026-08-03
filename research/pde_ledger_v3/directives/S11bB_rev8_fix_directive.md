# S11b-B — WRITE REVISION 8 OF THE BUILD DIRECTIVE

⛔⛔ **DOCUMENT TASK. ⛔ Do NOT write, run, or modify any `.wl` or `.py` script.**

You authored rev 7. Two independent reviewers have now assessed it. ⭐ **You are the author again for rev 8.**

## Artifacts

**Edit only:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md`, then
reassemble and verify byte-identity exactly as you did for rev 7, and report the new `sha256`.
⛔ **Do not edit the headers.**

## The review to fold

`/var/projects/toy_physics/research/pde_ledger_v3/_scratch/S11bB_dir_rev7_agent.txt`

⚠ Findings are **F1–F4**, plus the "Also-checks" paragraph. A second reviewer returned **no findings**;
⛔ that is not evidence F1–F4 are wrong — ⚠ the same reviewer cleared rev 6 entirely while five defects
were present, and **its own affinity derivation is cited inside F1 as the leak operating in practice.**

⭐ **Fix every finding that survives your own judgement, and say so if one does not.** ⛔ Do not transcribe
remedies — **derive and verify them** with SymPy/numerics.

## ⭐⭐ TWO DECISIONS ALREADY TAKEN — implement these, ⛔ do not re-litigate

### 1 · F1 — the affinity leak. ⭐ Structural fix required; ⭐ YOU choose which, and justify it

The reviewer offers two: **(a)** remove the `(v_bulk, 𝒜)` half of the power identity from the shared spec so
engines must derive the second pairing; **(b)** declare `𝒜` **SUPPLIED** — as §1b and §3b already are — and
delete the "construct it" instruction, the `S11BB_AFFINITY` derivation claim, and the affinity-power
check's "wrong derivation caught" claim, since it would then be pass-by-construction.

⚠ **Weigh this against the project's measured history:** revisions that left a load-bearing object
**under-determined** produced blockers (the branch rule, the derivation route), and both were ultimately
fixed by **supplying** a version that independent parties had verified. ⚠ `𝒜 = μ_s − δp/ρ_m` has now been
derived independently by **two** parties who agree. ⛔ But do not treat that as an instruction — **decide,
and state in one line what your choice costs.**

⚠ Whichever you choose: ⛔ **a check that can no longer fail must be DELETED, not weakened.**

### 2 · F2 — the reciprocal traction. ⭐ The project owner has decided: CARRY BOTH SYMBOLICALLY

The reviewer showed that with the traction fixed as `t_± = −δp_± n̂_±` and added face-force terms banned,
interfacial entropy production is a single flux–force pair `𝒜_±J_±` with **no force conjugate to `V_±`** —
whence admissibility forces `Λ_V⁰ = 0`. ⚠ That would void a headline result of the preceding sub-step.

⭐ **Implement this instead:** add a **reciprocal traction term with its own free coefficient**, so that

```
t_±  =  −δp_± n̂_±  −  Λ_X(ω) 𝒜_± n̂_±              (Λ_X free, Λ_X⁰ real)
```

or whatever equivalent form your own derivation shows is correct — ⭐ **derive it; do not transcribe the
line above.** ⛔ **`Λ_X⁰ = 0` must recover rev 7's specification exactly**, so nothing is asserted.

Then require the engines to **compute and report**:
- the **Onsager reciprocity relation** between `Λ_X` and `Λ_V`, if one is forced;
- whether admissibility still forces `Λ_V⁰ = 0`, **or** whether the reciprocal term rescues it;
- how this term changes the thickness equation and B5's roots.

⛔⛔ **Assert neither outcome.** ⭐ Whether the velocity coupling survives thermodynamics must be an
**output of this step**, not a premise. ⚠ `Λ_X` must also appear in the two-port power identity, the
passivity analysis, the dimension list, and as its own **FORM control** in B8.

## The remaining findings

- **F3** — §1b defines `q_out` only in the lower half-plane while §0 makes `Im ω > 0` first-class.
  ⭐ Extend the definition so it is unambiguous in the **upper** half-plane too. ⛔ Do not change the
  prescription itself — three independent parties have verified it.
- **F4** — §2b states the form of B2's answer, contradicting TASKS rule (4). ⭐ Either scope rule (4)
  explicitly to exclude supplied-route disambiguation, or remove the presupposing phrases. ⛔ Do not leave
  the rule silently violated.

## ⛔⛔ CONSTRAINTS

1. ⛔ **MINIMAL AND SURGICAL.** ⚠ Wholesale rewrites have bred every round's new defects.
2. ⛔ **NEVER SUPPLY AN EXPECTED RESULT OR ITS REASON.** ⛔ Do not state or hint what the longitudinal mode
   does, what sign any imaginary part takes, whether the transverse coupling vanishes, or whether `Λ_V⁰`
   survives.
3. ⛔⛔ **KEEP THE FALSIFICATION CHANNEL OPEN.** A **growing** root remains a first-class outcome. ⚠ Any
   diagnostic that could reject one must be **explicitly scope-bounded** to a sub-case where growth would
   necessarily be a derivation error.
4. ⛔ **EVERY CHECK MUST BE ABLE TO FAIL.** ⭐ For each check you add, keep, or modify, state in one line
   what wrong derivation it catches. ⛔ If you cannot, delete it.
5. ⛔ **DO NOT RE-OPEN** what independent legs cleared: the scope boundary; header symmetry and the
   `reduction/` bar; `B1`'s mass balance; the `A/B/C` split; `B8` controls B/C/D; §1b's branch prescription
   **as such**; §3b's virtual-displacement rule (verified correct by computation this round).
6. ⛔ **TWO ENGINES MUST NOT BE ABLE TO DIVERGE.**
7. ⛔ **No new scope** beyond what a finding or decision above requires.

## Output

The edited files, plus a report **under 60 lines**: each finding, what you changed, the `sha256`, your F1
choice with its one-line cost, and for every check added/kept/modified the one-line statement of what wrong
derivation it catches.

⛔ Then stop and exit.
