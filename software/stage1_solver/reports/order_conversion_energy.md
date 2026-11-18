# Can the order equations supply the two-drain energy?

**Deliverable:** identify what the existing conversion equations determine about the required 0PN power, and the effective mouth datum still needed, without constructing a throat interior.

**Result.** They identify possible energy contributions but do not select an available power source. The specified bare order potential has equal-depth phases. Loss of order alone therefore establishes neither released energy nor its transmission into the mouth's mechanical motion. This remains a requirement on the conversion response.

**What is supplied.** The accepted balances conserve constituent number and track order:

$$
D_t\chi=\Gamma_{\rm bal}-\frac{\nabla_4\cdot J_\chi}{n}.
$$

The historical real-fraction energy adjunct supplies the partial order work

$$
P_{\rm order}=\int\mu_\chi D_t\chi\,d^4X,
\qquad \mu_\chi=\delta F/\delta\chi.
$$

Here $\mu_\chi$ has energy-per-four-volume units; no extra factor of $n$ multiplies it. This contribution is not a complete moving-region energy law. The live source functions, energy partners and reference response remain open, and S12 excludes the puncture interior. The separately labelled relaxation adjunct cannot be added to $\Gamma_{\rm bal}$ without a consistent rate convention. [Accepted source scope](../../../research/pde_ledger_v3/steps/S12_order_source.md) · [Order-work scope](../../../research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md#86-order-work-and-energy-reference-records-c-8) · [Prior work assessment](../../../docs/conversion_work_assessment_and_next_step.md)

**What the specified energy terms can say.** For the postulated $f_B=a_B\chi^2(1-\chi)^2$, both homogeneous endpoints have zero bare order energy. In a homogeneous diagnostic with prescribed $D_t\chi=-\gamma<0$, partial order power is $+3a_B\gamma/16$ at $\chi=3/4$ and $-3a_B\gamma/16$ at $\chi=1/4$. Thus the same signed loss rate can raise or lower this energy. Its integrated change over $1\to0$ is zero; the barrier is $a_B/16$. This kinematic diagnostic is not a physical conversion trajectory or a thermodynamic latent-heat calculation.

| Existing term | Energy information still needed |
| --- | --- |
| $\chi f_{\rm shear}$ | Stored strain energy could change during conversion. Its live strain/reference history and transfer to motion, heat or other fields are unselected. The shear gate does not create stored strain. |
| $\kappa_B|\nabla\chi|^2/2$ | Changing an interface can change its stored energy; the physical profile, geometry and boundary work are needed. A stationary interface is not automatically a continuing energy supply. |
| $U(n)$ and pressure | Actual feeding and receiving states must determine their energy difference. A common conserved-number reference $C_{\rm ref}n$ cancels between equal-number transfers; choosing its zero supplies no power. |
| $f_{\rm throat},f_{\rm mix}$ | These remain placeholders in the cited energy. No sign, magnitude or mechanical coupling is supplied. |

The [historical functional and provenance](../../../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md#the-postulated-action-and-balance-laws-all-postulated-labels-printed) are conditional inputs. No term was promoted to a physical power measurement.

**Minimal requirement without the interior.** Let $\ell_i$ denote net energy made available per unit mouth throughput, with storage, other outputs and any external supply identified and counted once. The [ideal pair result](two_drain_energy.md) requires, at fixed $D_i=D(a_i)$,

$$
D_i\Delta\ell_i=P_{\rm mech}\quad(i=1,2),\qquad
\sum_iD_i\Delta\ell_i=2P_{\rm mech}.
$$

A fixed energy yield at unchanged throughput supplies unchanged baseline power. The interaction increment requires a changed yield, changed storage, or redistribution from another output. A derived or independently measured mouth energy-response relation could supply this datum without resolving its interior; the existing equations do not yet give it.

**Evidence and stop.** SymPy and Mathematica agree on the homogeneous power diagnostic, live shear contribution, conserved-reference cancellation and required per-mass response. [Sources, scripts and guarded outputs](../../../docs/_measurements/order_conversion_energy.json). No new conversion law or interior calculation was introduced. Base `fd667ea9291b4430a96afe3fd932926232f5507d`; working, uncommitted pending independent review.
