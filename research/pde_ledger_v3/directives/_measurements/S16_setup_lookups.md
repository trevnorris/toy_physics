# Measurements: S16 decision list, setup lookups (generated 2026-10-21 19:57)

Generator: `_scratch/s14a/gen/s16_setup_lookups.sh` (sed/grep/git/sha256sum only; this file is written by
it). Mechanical lookups only: existence, verbatim retrieval by line, literal-match counts.

## L0. The list

```
$ cd /var/projects/toy_physics && sha256sum research/pde_ledger_v3/directives/S16_decision_list.md && git log -1 --format='%h %s'
29d3225e559741b5ae71a677bd926de52903ee702a0fc24df0be9f2ef39a9980  research/pde_ledger_v3/directives/S16_decision_list.md
6ad2d816 Record the S16 decision list review pass 2 (spec content)
```

## L1. The plan's S16 section

```
$ grep -nE '^### S1[5-9]' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
652:### S15 · Two drains attract ⚠ RUN S16 FIRST for regime context — ⛔ NOT as S15's premise
670:### S16 · ⭐⭐ The worldtube reduction — the theorem the whole scope boundary rests on
747:### S17 · The PN ladder — cite-only
760:### S18 · 2.5PN Burke–Thorne matchback
766:### S19 · Gravitomagnetism ⛔ NOT a cite-and-record step — it has a known open problem
```

```
$ sed -n '739,741p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
⇒ **Rewrite S16 as a conditional response-side monopole approximation** — ⛔ not a source-generation or
interior-independence theorem — and list `M`, `a`, the multipoles, compactness, source normalization and
boundary-flux suppression as **supplied inputs** that belong in S22's debt table.
```

```
$ sed -n '700,703p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
⛔ **QUOTATION INTEGRITY — the block above is the source AS WRITTEN.** `4d_2_5pn.tex:605-614` writes the remainder **without** the mass factor,
which makes it dimensionally inhomogeneous (force on the left, acceleration in the remainder). ⛔ **Do
not present the repaired formula as a quotation.** Quote the source as written, record the missing-mass
**typo**, and cite the dimensionally correct form from `4d_1pn_full.tex:886`.
```

```
$ sed -n '705p;729,731p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
⛔⛔ **GRAVITY-SECTOR PLUG — REVISIT BEFORE COMMITTING THE SPECIES/INTERIOR-DEPENDENT RESPONSE (banked 2026-09-15).**
  medium chain — a warning to not corner ourselves here, ⛔ not a step to run now. Outcome 5 (can't compute
  without importing κ=1 ⇒ no muonium prediction) is a legitimate honest result.

```

```
$ sed -n '732,737p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
⛔ **A11 GATE — this step must not substitute one length for another.** The source's `a` is a
**Gaussian/profile support width** controlling `Q_ij = M a² δ_ij / 2`; the model's **mouth radius** is a
different object, and the bridge report says outright it *"is not an invariant reduction width"*.
⇒ **Use distinct symbols — `a_WT` (worldtube profile width) and `a_mouth`** — and list the **missing
bridge between them** as debt in S22. ⛔ Never substitute one for the other. Same dimension is not same
object.
```

## L2. Framing, "Keep S12" and R-S14a-01

```
$ sed -n '729,731p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
- **Flow comparison decision, user 2026-10-19** — Each force states the flow it needs; S21
  compares those requirements, and a mismatch reopens S9b. S12 states the source compatibility
  requirements and leaves forms open; it fits no flow and adds no optical calculation or
```

```
$ sed -n '15,16p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md
**User decision, 2026-10-20: "Keep S12".** **"S12 stays open for the forms. This record closes only
S12's requirements part, and the forms come once the brane is defined, after the force sectors.
```

```
$ cd /var/projects/toy_physics && git log -1 --format='%h %s' 79511272
79511272 Accept the S14a part 2 record and its register pass (round 2 clear)
```

```
$ grep -nE '^### R-S14a-01|^- \*\*status\*\*' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | grep -A1 'R-S14a-01'
938:### R-S14a-01 — the dynamically selected slow limit and far behaviour
944:- **status** OPEN
```

## L3. The accepted part 1 record's projected objects

```
$ cd /var/projects/toy_physics && git log -1 --format='%h %s' 4a4dfaf5 && git diff --stat 79511272 -- research/pde_ledger_v3/steps/S14a_1_projection.md | wc -l
4a4dfaf5 Accept the S14a part 1 record and its register pass (round 1 clear)
0
```

```
$ sed -n '62p;84p;87p;91p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S14a_1_projection.md
| Material face factor; accepted S12 record 2026-10-20, supplied by the 2026-10-21 list | **"It is supplied and unfalsifiable in this build. The build's face terms are reached by its own integration."** | Printed one-sided factor `a_{r,s}=u_r^w−∂t w_s−u_r^i∂i w_s`; it is not a derivation target substituted to build the face terms. |
| Projected material density and in-plane flux | `Σ=∫n₂dw`, `F^i=∫n₂u₂^i dw` | b51/b52: each compared component has exact residual `0`. |
| Material leaving each face, before / after full concentrated part | Lower: `−n₂a₂,₋` / `−n₁a₁,₋`; upper: `n₂a₂,₊` / `n₃a₃,₊` | b53/b54/b56/b57: exact `0`. Both forms are retained; no material face condition is solved/substituted. |
| Projected material balance | `∂tΣ+∂iF^i−n₁a₁,₋+n₃a₃,₊=0` | b59 lhs and rhs: exact `0` between engines. This does not say the lhs vanishes for arbitrary unconstrained fields. |
```

```
$ sed -n '45p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S12_decision_list.md
- `n` is the material number density on the 4D domain (in-plane `x^i`, bulk direction `w`), `u` its velocity,
```

## L4. The source's statements

```
$ sed -n '695,720p' /var/projects/toy_physics/research/pde_ledger_v3/../4d_1pn_full/paper/4d_1pn_full.tex
The reduction carried out below uses the following explicit assumptions.
\begin{enumerate}[leftmargin=1.5em]
  \item \textbf{Compact support / small body.}
  The defect occupies a region of characteristic size \(a\) that is small compared with
  the inter-body separation and with the scale over which the external field varies:
  \(a\ll r\).

  \item \textbf{Smooth external field across the support.}
  Over the defect worldtube, the external potential may be Taylor-expanded about the
  center-of-mass position.

  \item \textbf{Negligible worldtube boundary fluxes at the retained order.}
  Any residual leakage, entrainment, or stress flux through the defect boundary is beyond
  the Newtonian order kept here.

  \item \textbf{Defect-centered coordinates.}
  The spatial origin is chosen so that the defect dipole moment vanishes.

  \item \textbf{Quasi-static near zone.}
  The dropped terms in the exact longitudinal identity remain parametrically small, so
  the Poisson hook of \cref{eq:dictionary-poisson-hook} is valid.

  \item \textbf{Conservative sector only.}
  Radiation, dissipative exchange, and higher-response effects are deferred to later
  work and are not part of the Newtonian limit derived in this section.
\end{enumerate}
```

```
$ sed -n '761,765p' /var/projects/toy_physics/research/pde_ledger_v3/../4d_1pn_full/paper/4d_1pn_full.tex
For the isotropic Gaussian \eqref{eq:defect-gaussian-profile},
\begin{equation}
Q_{ij}=\frac{Ma^2}{2}\,\delta_{ij}.
\label{eq:defect-isotropic-second-moment}
\end{equation}
```

```
$ sed -n '885,891p' /var/projects/toy_physics/research/pde_ledger_v3/../4d_1pn_full/paper/4d_1pn_full.tex
\begin{equation}
M\ddot{\mathbf X}_{\rm cm}
=
- M\,\gradThree\PhiN(\mathbf X_{\rm cm},t)
+
O\!\left(M\,\frac{a^2}{r^2}\,\frac{\PhiN}{r}\right).
\label{eq:newtonian-center-of-mass-eom}
```

```
$ sed -n '605,614p' /var/projects/toy_physics/research/pde_ledger_v3/../4d_2_5pn/paper/4d_2_5pn.tex
For a compact defect of size $a\ll r$ in a smooth external field, the center-of-mass
reduction gives the usual point-particle force law
\[
M_A\,\ddot{\mathbf X}_A
=
-M_A\nabla\Phi_{\rm ext}(\mathbf X_A)
+O\!\left(\frac{a^2}{r^2}\,\frac{\Phi_{\rm ext}}{r}\right),
\]
so that the universal point-particle source at leading order is controlled by the
worldline/worldtube multipoles rather than by arbitrary internal details of the defect.
```

## L5. Charter, S22 debt rows, C3

```
$ sed -n '74,80p' /var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md
S16/S19 may show the premise fails (see §1.1). The quotations below are that conditional result, not a
licence:

> *"controlled by the worldline/worldtube multipoles rather than by arbitrary internal details of the
> defect"* — `research/4d_2_5pn/paper/4d_2_5pn.tex:613`
> *"we do not attempt to solve the fully dynamical moving-throat problem in the bulk"* —
> `research/4d_1pn_full/paper/4d_1pn_full.tex:110`
```

```
$ sed -n '109,119p' /var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md
the conditional response-side conclusion S16 actually permits. ⛔ It does **not** establish that the far
field is independent of the interior.

⛔⛔ **AMENDED 2026-07-31 — the boundary is weaker than stated above, and may not hold at all.**
Reading the worldtube theorem's actual statement rather than its slogan
(`research/4d_2_5pn/paper/4d_2_5pn.tex:605-614`), the independence is **leading-order only** and
**conditional on the defect being compact (`a ≪ r`)**, with a first correction `O(a_WT²/r²)` — controlled by the
**worldtube profile width**, which is a different object from the mouth radius (A11) and is itself
undetermined.

⚠ And the Spin Problem (**C10**) concludes you *cannot* get frame dragging from a compact defect —
```

```
$ sed -n '1303,1305p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
| `a_WT` (worldtube profile width) | controls the `O(a_WT^2/r^2)` correction |
| the `a_WT` / `a_mouth` BRIDGE | DOES NOT EXIST (A11). Never close it by substitution |
| `M` and the supplied multipoles | inputs of the S16 reduction, not outputs |
```

```
$ sed -n '171p' /var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md
| **C3** | **`m_defect` is a GAP** — `INFLOW_MASS_SOURCE_MISSING` / `BLOCKS_MASS_EMERGENCE`. Only a dimensional bridge `α_J ħJ/c_γ²` exists; pathA_21's verdict is `MASS_BRIDGE_FORM_NOT_DERIVED`, *"No action-level, boundary-source, Noether, or Hamiltonian equation"* | `stage004:181`; `parameter_register.md:152-153`; `pathA_21_emergent_G_mass_bridge.md:98` | **OPEN** |
```

```
$ sed -n '21,25p' /var/projects/toy_physics/notes/muonium_gravity_corrections.md
2. **Three-mass ledger — keep them separate.** Inertial `m_i`, **passive** `m_p` (response: how it falls — what
   LEMING measures), **active** `m_a` (source: the field it produces). Define each independently; report
   `κ^pass = m_p/m_i` and `κ^act = m_a/m_i` separately. Equality `m_i=m_p=m_a` must **emerge as a Ward-identity
   theorem** if it holds — ⛔ never assumed. (This is "keep every varying quantity live; a required freeze is the
   finding.")
```

## L6. O2's S16 rows, its inertia and density rows, and the S12 owner row

```
$ cd /var/projects/toy_physics && git log -1 --format='%h %s' 72866fcf
72866fcf Accept the O2 record r4 (sub-step 7 closed; O2 closed)
```

```
$ sed -n '133p;138p' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
| `ℐ_br^live`; **OPEN**, C §3 | Flowing/embedded inertial and kinetic response, entering momentum storage and transport. No identification with `ρ_br V`, a quadratic live kinetic energy, or a constant inertia is supplied. | S8, with S1.5 antecedents when available. |
| `ℛ_br` (O7); **OPEN**, C §6 | General brane-density response with live bulk, flow, embedding, thickness/projection dependence and gradients. It is input to the density in mass/momentum/energy accounting; it is not inferred from the bulk EOS or slab factorization. | S8, substrate reduction and S22 completion. |
```

```
$ sed -n '607p' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md
| **S16** | Response-side connection to the independent orbital `GM`. The inventory's designation of S16 as this matching interface is an **inference**. Its compact-body, supplied external-field/worldtube approximation retains calibrated potential, supplied mass/multipoles and discarded boundary-flux assumptions; its normalization is target-matched, not derived. A response-side match does not derive the drain's source or coupling. |
```

```
$ sed -n '164p' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md
| Orbital/test-matter response and matching to independent `GM` | **S16**, `:670–698`: the worldtube reduction is a conditional **response-side** approximation for compact matter in a supplied external field. v9:132–134 assigns `GM` to the gravity sector; treating S16 as the response-side owner of its orbital connection is this inventory's **inference**, not a plan assignment. S16 assumes **a calibrated potential `−Gm/r`** and **supplied mass and multipole moments**, with boundary flux discarded (`:696–698`), and records `κ_ρ=1` as target-matched, not derived (`:716`). A route matching `F_drive` to `GM` through S16 therefore inherits that calibration and does not derive the coupling. S16 is not a derivation of the field's source. Source and response matching remain separate obligations. |
```

```
$ sed -n '292p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md
| Interior-order law/history/strain; faces' law; source/flux forms | No separate owner named for the first two; source/flux forms remain S12's. Slab material completion remains separately S8/S22. None of the forms is delivered here (M01, M11, M18). |
```

## L7. The plan's S19 and its register line

```
$ sed -n '766p;786p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
### S19 · Gravitomagnetism ⛔ NOT a cite-and-record step — it has a known open problem
**Register:** **C10**. ⭐ **Two independent legs, both blocking** — this is new physics, not bookkeeping.
```

## L8. Literal-match counts in the list (leak probe; each hit is read by the orchestrator)

Patterns: the source's potential and force-law symbols, the calibrated potential, named prior-art results,
and the withdrawn lepton results named in the plan's guardrail.

```
$ for p in 'Φ' 'Phi' 'Gm' 'G M' '−M' '-M' 'Lagally' 'Bjerknes' 'κ' 'kappa' '1.85' '11:2:5' '57/64' 'Ward'; do printf '%s: ' "$p"; grep -cF -- "$p" /var/projects/toy_physics/research/pde_ledger_v3/directives/S16_decision_list.md; done
Φ: 0
Phi: 0
Gm: 0
G M: 0
−M: 0
-M: 0
Lagally: 0
Bjerknes: 0
κ: 0
kappa: 0
1.85: 0
11:2:5: 0
57/64: 0
Ward: 0
```

```
$ grep -nF -- 'potential' /var/projects/toy_physics/research/pde_ledger_v3/directives/S16_decision_list.md
24:  - No source normalization is used, because no potential is supplied.
42:- **Outside field — "No potential":** "The outside field enters only as the live flow and stress around the region.
54:The two later decisions narrow the last two sentences of "No potential". The gradient-form condition and the two
137:  compactness, vanishing dipole, boundary-flux suppression, potential, mass-to-count ratio, or symmetry of the
205:- no potential;
```
