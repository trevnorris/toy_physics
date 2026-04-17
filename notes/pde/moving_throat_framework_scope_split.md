# What the current moving-throat PDE paper does **not** cover, and what should become separate papers

## Bottom line

No — we did **not** cover everything in the current framework paper, and that is the right outcome.

The compact program file is much larger than a single publishable framework paper. The paper we drafted is the **reader-facing PDE framework / usage paper**. It explains the exact parent theory, the moving-throat lift, the reduced wall--support--gauge stack, the grouped real `P2` bundle, the isotropic normalization target, the first weak-axisymmetric observable, and how a reader is supposed to use the PDE.

That is already a full paper.

The compact master still contains several additional theorem clusters that are **not** just appendix material and should not be forced into this framework paper. Some of them belong in the technical derivation supplement. Others deserve their own higher-level papers.

---

## What the current framework paper **does** cover

The current paper is the right home for the following material:

- the **claim-status map** and notation firewall,
- the **exact parent `4+1` theory**,
- the **moving-throat geometry lift**,
- the **linearized moving-throat PDE**,
- the **reduced conservative wall--support--gauge system**,
- the **full grouped real `P2` bundle** as the first nontrivial payload,
- the **outgoing quadrupole bridge and isotropic normalization**,
- the first **weak-axisymmetric departure** through `Xi_1 = P_1 / P_0`,
- and a **how-to-use** section that tells people what branch data they must actually compute.

That scope is coherent. It tells the reader what the PDE is, what quantities it controls, what is exact versus reduced, and what the main operational theorem target is.

In compact-file language, this paper is basically built from:

- Section `2` (parent theory),
- Section `3` (geometry lift),
- Section `4` (reduced wall/support engine),
- the **grouped-`P2` front end** of Section `7`,
- plus a paper-facing discussion of the realization gap.

That is the right boundary for Paper 1.

---

## Why everything should **not** be forced into this one paper

There are three reasons.

### 1. The current paper is a **framework** paper, not a realization paper

Its job is to explain the PDE and the observable extraction chain. Once we start including the later branch-selection, same-charge, rigid-mouth, twin-support, and open-system barrier machinery, the center of gravity shifts from “what the PDE is” to “what one very specific realized branch does.” Those are different papers.

### 2. The compact file mixes multiple theorem families

The compact document is not one narrative. It is a full program ledger. It contains:

- the exact parent theory,
- reduced PDE structure,
- branch-specific mouth/core closures,
- coherent invariant/orbit machinery,
- same-charge and rigid-mouth selection,
- relaxed open-system/barrier continuations,
- and even export/materials companions.

Those are related, but they are not one journal article.

### 3. Keeping the claim-status hierarchy visible is a strength

One of the best features of the compact program is that it clearly separates:

- exact statements,
- exact-within-closure statements,
- reduced results,
- effective closures,
- numerically located branches,
- and genuinely open realization gaps.

Trying to compress all of that into one paper would blur exactly the distinction that currently makes the project readable and referee-safe.

---

## What the current framework paper does **not** cover

## 1. Family-1 mouth/core closure and the co-evolving mouth/core branch

This is the material in compact Sections `5` and `6`.

At a high level, this block contains:

- the normalized **Family-1 core/mouth variables**,
- the **positive localized mouth-source theorem**,
- the explicit **mouth boundary-layer law**,
- the explicit **core-to-mouth gain map**,
- the self-matched mouth susceptibility closure,
- the canonical mouth branch,
- the **co-evolving mouth/core fixed-point map**,
- and the renormalized co-evolving branch.

Why it should not live in Paper 1:

- This is not needed for a reader to understand the PDE itself.
- It is a **branch-realization** paper, not a PDE-definition paper.
- It is already substantial enough to read like a separate narrative: “given the PDE framework, how do explicit mouth/core closures realize an admissible branch?”

### Recommended separate paper

**Paper C — Family-1 Mouth/Core Closure and the Co-Evolving Branch**

Suggested emphasis:

- mouth/core variables and compensation structure,
- explicit mouth-layer law,
- gain map,
- co-evolving fixed point,
- renormalized canonical branch,
- and what this means for realized throat data.

---

## 2. The coherent local-kernel / invariant / orbit / realization-compiler program

This is the large middle block in compact Section `8`, especially `8.1` through `8.15`.

At a high level, this block contains:

- the coherent local-kernel effective ratios,
- the coherent placement map,
- the support-compensation theorem,
- the microscopic slippage ledger,
- the triangular defect normal form,
- the direct monomial invariants,
- the similarity orbit and quotient closure,
- the branch-observable packet,
- the transfer-shape / outgoing-prefactor compiler,
- the canonical orbit projection,
- the free-quintuple target graph,
- the graph-error packet,
- and the finite local mixed-ray search / sieve.

Why it should not live in Paper 1:

- This is a **realization and admissibility** program built on top of the PDE, not part of the first explanation of the PDE.
- It introduces a new language of invariants, orbit packets, graph slices, and local search that would completely dominate the framework narrative.
- It is higher-level and important, but it is a different paper-level question: “how are realized branches classified and searched once the PDE front end is fixed?”

### Recommended separate paper

**Paper D — Coherent Local-Kernel Reduction, Invariants, Orbit Lock, and the Realization Compiler**

Suggested emphasis:

- coherent placement,
- slippage and monomial invariants,
- similarity-orbit / quotient closure,
- orbit packet,
- free-quintuple realization compiler,
- finite local search,
- and how these organize actual branch realization.

---

## 3. The same-charge / weak-axisymmetric / rigid-mouth / selected twin-support packet

This is the later part of compact Section `8`, especially `8.16` through `8.20`.

At a high level, this block contains:

- the **one-port same-charge mixed-kernel verdict**,
- the actual-branch prefactor ceiling,
- the weak-axisymmetric transported prefactor slope `Xi_1`,
- the exact weak-axisymmetric line `b_{P_0} = 3 a_{P_0}`,
- the **rigid-mouth physical normal form** in `(U,V)`,
- Cartesian orbit lock,
- the **selected twin-support slice**,
- primitive ranking,
- actual twin-support placement,
- and the support-versus-orbit split.

Why it should not live in Paper 1:

- Paper 1 already uses `Xi_1` as the **first weak-axisymmetric observable**, and that is enough.
- The full same-charge / rigid-mouth / support-selection story is a separate theorem family about **which physical branches are actually admissible**.
- That material is not “just more PDE.” It is a distinct branch-selection paper.

### Recommended separate paper

This material can either remain part of **Paper D** above, or be broken out again as:

**Paper E — Same-Charge Actual-Branch Selection, Rigid-Mouth Normal Form, and Twin-Support Placement**

That split is worth making if you want the realization/orbit paper to stay mathematical and the same-charge/selection paper to stay physical.

---

## 4. The relaxed/open-system barrier extension and the physical companions

This is compact Section `9`.

At a high level, this block contains:

- the codimension-three **relaxed branch**,
- the exact recovery slice,
- the selected-branch leakage/work compiler,
- the non-rigid mouth/dressing packet,
- the compensated multimode source compiler,
- the lowered stationary barrier front end,
- the dynamic event-chain / threshold / WKB / Goldilocks continuations,
- the microscopic cubic/quintic export kernel,
- the vacuum/lattice heat partition,
- and the physical calibration / material-screening companion.

Why it should not live in Paper 1:

- This is an **open-system / barrier / export / materials** companion.
- It is downstream of the core PDE framework.
- It mixes physics questions that are interesting in their own right, but they are not part of the clean “what is the PDE and how do I use it?” story.

### Recommended separate paper

**Paper F — Relaxed Open-System Branch, Barrier Physics, Export, and Materials Calibration**

Suggested emphasis:

- relaxed branch declaration,
- lowered front end,
- dynamic threshold/event-chain/WKB compilers,
- export kernel,
- heat partition,
- and materials-screening/calibration dictionary.

---

## 5. The actual realization theorem is still open

This is the substance of compact Section `11`.

The framework paper should absolutely say that there is an open realization gap, but it should **not** pretend to close it.

What remains genuinely open includes, among other things:

- the realized one-port branch data,
- the actual `Delta_norm` / `Xi_1` packet,
- whether one and the same realized branch satisfies the outgoing normalization, the weak-axisymmetric ceiling, and rigid-mouth orbit lock together,
- the actual selected twin-support placement state,
- and, if the relaxed branch is kept, the realized leakage/export/materials packet.

That means there is still a real difference between:

- a completed **framework paper**,
- a completed **derivation supplement**,
- and a completed **branch-realization theorem**.

They should not be written as if they are the same accomplishment.

---

## 6. The full derivation chain still belongs in its own companion paper

Even after the higher-level papers are split off, the full algebraic derivation still belongs in a separate document.

That paper should be built from the long derivation record and the symbolic notebooks/scripts. Its job is:

- to prove the reduced formulas used in the framework paper,
- to expose the Schur complements, projector algebra, and normalization bridges,
- and to give referees a reproducible technical audit trail.

That is a very different reading experience from the framework paper.

### Recommended separate paper

**Paper B — Technical Derivation and Symbolic Verification Companion**

---

## Recommended publication stack

## Minimum sensible stack

1. **Paper A — Moving-Throat PDE Framework and Usage**  
   The paper we just drafted.

2. **Paper B — Technical Derivation and Symbolic Verification Companion**  
   Full derivation walk-through and script-backed proof chain.

3. **Paper C — Family-1 Mouth/Core Closure and the Co-Evolving Branch**  
   The first realized mouth/core branch paper.

4. **Paper D — Coherent Local-Kernel Reduction, Orbit/Invariant Structure, and Actual-Branch Realization**  
   Includes the monomial invariants, similarity orbit, realization compiler, and optionally the same-charge / rigid-mouth / twin-support packet.

5. **Paper E — Relaxed Open-System Barrier, Export, and Materials Companion**  
   Keeps the relaxed/open-system physics and materials-facing story out of the framework paper.

## Lower-count alternative

A lower-count alternative is:

- keep **Paper A** and **Paper B** unchanged,
- merge **Paper C** and **Paper D** into one larger **branch-realization paper**,
- keep **Paper E** as the open-system/barrier/materials companion.

That gives a four-paper stack instead of five.

I do **not** recommend compressing below that unless the goal is a very long monograph rather than a clean paper sequence.

---

## Direct answer to the question

You were **not** wrong.

The current framework paper does **not** cover everything in the compact project file, and it should not try to.

What we have now is the right **first paper**: it teaches the PDE, the reduction chain, the grouped `P2` observables, the isotropic normalization target, and the first weak-axisymmetric deformation.

What remains outside it is real and important, but it belongs to other papers:

- the **Family-1 mouth/core realization story**,
- the **coherent invariant/orbit/realization story**,
- the **same-charge / rigid-mouth / support-selection story**,
- the **relaxed open-system / barrier / export / materials story**,
- and the **full derivation proof chain**.

That is a healthy split, not a failure of coverage.
