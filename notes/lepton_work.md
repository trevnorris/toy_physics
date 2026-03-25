## 1. Executive summary

### 1.1 What changed this session

This session took the lepton audit from a broad “mixed-sector Berry rotor is the last live corridor” statement to a much sharper conditional closure. We started from the corrected ontology, accepted the audit’s conclusion that the mixed ((a-w)) sector was the only live same-charge corridor, and then pushed that corridor through a sequence of increasingly strict tests: finite-throat replacement of the old Gaussian/noncompact picture, exact D/N mouth-branch reconstruction, dynamic pinning of the transverse odd radius, induced (P_{22}) splitting from area-preserving deformation, repeated topological no-go tests, selective (\tau)-subbundle analysis, recirculation/plumbing holonomy, and finally standing-wave and autonomous-soliton closure tests. The overall effect was not to produce a free-standing lepton theorem from the frozen files alone, but to isolate one very specific closure under which the half-integer same-charge quantizer does go through.

### 1.2 Main result in one paragraph

The main result of the session is that the minimal charged-lepton route is now **conditionally viable** in a much sharper sense than before. The corrected stack already supplies the electric branch ((\eta_Q,q_*,q_{\rm eff})), and the lepton audit had already isolated the mixed-sector Berry rotor as the only live same-charge corridor. What this session added was a candidate closure chain: if the same-charge mixed sector lives on a selective (\tau)-subbundle, if its loop holonomy closes with the central sign (-I_2), and if the isolated particle is an autonomous self-reproducing standing-wave soliton so that the common scalar phase cancels, then the Berry Wilson loop forces
[
\alpha \in \mathbb Z+\frac12,
]
which is exactly the half-integer discretizer the audit said was missing. That is the strongest positive result reached so far, but it is still conditional on the autonomous-eigenmode closure rather than a theorem already frozen in the files.

### 1.3 Final conditional verdict

The session did not prove that the current conservative file stack already contains a charged lepton. It did something more useful for the next session: it turned the old vague stall into a precise theorem-plus-caveat structure. The theorem is conditional: **if** an isolated defect is an autonomous self-reproducing GNLS soliton, and **if** the same-charge mixed sector closes with the central sign on its selective subbundle, then the common scalar phase becomes trivial and the mixed-sector Berry flux is locked to the half-integer lattice. The caveat is equally precise: the current files still leave the full particle-scale recirculation/plumbing law, dissipative leakage/radiation reaction, and the fully solved moving-throat dynamics unresolved, so the autonomous-eigenmode closure remains the only nontrivial assumption still separating the present work from a fully file-grounded theorem.

## 2. Starting point and frozen ontology

### 2.1 Corrected variable dictionary

This session began from the corrected ontology, not from the older circulation-based charge dictionary. The frozen charge bookkeeping is
[
\eta_Q=\pm1,\qquad
q_*=\eta_Q e_*,\qquad
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int}=\int Z(w),dw,
]
with (\kappa_\rho=1) on the gravity-side Newtonian scalar ledger. In this corrected stack, (\eta_Q) is the electric charge sign, (q_*) is the microscopic branch charge, and (q_{\rm eff}) is the observable brane charge. Quantized circulation/winding is no longer the electric-charge definition; it is kept separately as the magnetic/vortical label (n_\Gamma). The proposed same-charge internal label is (\tau=\pm1), which belongs to the mixed sector rather than the electric branch.

### 2.2 Non-negotiable ontology fixes

Several corrections were treated as non-negotiable throughout the session. The old relation (q\sim a^2\Gamma) was not allowed back in as the charge dictionary. Circulation/winding was treated as belonging to the magnetic/vortical sector rather than defining electric charge. Likewise, (\kappa_\rho=1) was treated as a gravity-side mass-dressing coefficient, not as electric charge. The mixed channels
[
A_w,\qquad J^w,\qquad F_{\mu w},\qquad F_{aw}
]
were treated as suppressed only in the controlled far-field brane limit, not erased from the microscopic ontology. That last point mattered because the entire session’s live corridor depended on those mixed channels remaining real hidden-core degrees of freedom rather than forbidden variables.

### 2.3 Frozen carry-forward constraints

The frozen coefficients and geometric constraints carried into this session were those already preserved by the corrected 1PN/2PN stack:
[
\kappa_\rho=1,\qquad
n=5,\qquad
\kappa_{\rm add}=\frac12,\qquad
\kappa_{\rm PV}=\frac32,\qquad
\beta_{\rm 1PN}=3,
]
together with the parity-even wake data
[
\alpha^2=\frac34,\qquad
a_H=0,\qquad
K_{\rm vec}=\frac{2}{\pi^2},
]
and the preferred throat aspect ratio
[
L/a\approx1.85.
]
The 2PN particle was also treated as already more than a scalar monopole: it carries a dipole wake, a real (P_0\oplus P_2) mouth/support layer, and a separate geometry-closure channel. These facts mattered because the session was not trying to invent particle structure from nothing; it was trying to determine whether the structure already visible in the corrected stack could plausibly host the minimal charged-lepton package.

### 2.4 Frozen versus proposed objects

One of the most important discipline points this session kept was the separation between frozen and proposed quantities. Frozen objects included the electric branch variables ((\eta_Q,q_*,q_{\rm eff},Z_{\rm int})), the magnetic/vortical label (n_\Gamma), and the carry-forward 1PN/2PN constants. Proposed objects included the same-charge mixed-twist sign (\tau), the transverse odd collective coordinates ((b_x,b_y)), the pinned transverse radius (s_*=b_x^2+b_y^2), the mixed-sector Berry coefficient (\nu_0), and all later selective-subbundle / recirculation closure laws built on top of them. This distinction mattered because many of the successful formulas derived later in the session were **conditional** formulas built on proposed closures rather than direct frozen theorems of the current file set.

## 3. Audit baseline at session start

### 3.1 What was already alive

At the start of this session, the lepton audit had already reduced the live same-charge possibilities to a very short list. First, the corrected stack already supplied the electric branch cleanly through (\eta_Q), (q_*), and (q_{\rm eff}). Second, once the mixed ((a-w)) sector was kept alive, the model could support a same-charge nonorbital internal odd coordinate. Third, the mixed-sector route could support a nonzero Berry/gyroscopic coefficient because the parent matter action is first order in time and the first odd (w)-mode survives in derivative overlaps even when it is hidden from a direct brane source. Finally, if a transverse radius (s_*) ever pinned, the low-energy dynamics would become those of a genuine internal rotor rather than a mere relabeling of orbital circulation. That was the exact launch point for this session’s work.

### 3.2 What was already dead or stalled

The audit endpoint had also already ruled out or stalled several routes before this session began. Ordinary circulation and orbital relabelings failed the intrinsic-spin bar. The scalar-bundle rescue of the magnetic sector failed to produce the missing Dirac piece. The static side/orientation route failed because a pure (+w/-w) switch is just the opposite-charge branch (\eta_Q), not a same-charge spin label, and because the frozen files did not provide a driven (P_{22}) anisotropy source. The axisymmetric Family-1 wall route could at best produce a continuous transverse rotor or an Ising-like odd sign, not a Pauli-like same-charge doublet. The Gaussian noncompact localization scheme also provided no built-in half-integer or topological Berry-flux quantization. In other words, by session start the audit had already cleared away the easy interpretations and left only one serious corridor alive.

### 3.3 Why the mixed-sector Berry rotor was the only live corridor

The reason the mixed-sector Berry rotor survived while the others failed was structural rather than ad hoc. The parent matter theory is first order in time, so a genuine collective-coordinate Berry/Magnus term is kinematically allowed from the start. At the same time, the 4+1 Maxwell/plasma ontology explicitly retains the mixed channels
[
A_w,\qquad
F_{\mu w},\qquad
J^w,\qquad
E_w,\qquad
F_{aw},
]
outside the strict zero-mode brane Maxwell limit. Those are precisely the channels that can supply a first-order internal phase law even when the static axisymmetric wall sector cannot. This is why the audit treated the reduced mixed-sector model
[
L_{\rm tr}
==========

\frac{M_b}{2}\dot b_i\dot b_i
+
\frac{\tau\nu_0}{2}\epsilon_{ij}b_i\dot b_j
-------------------------------------------

V_{\rm stat}(b_x^2+b_y^2)
]
as the only remaining same-charge corridor: it was the first point in the whole program where a same-charge internal sign and a genuine first-order internal symplectic term coexisted in one reduced model.

### 3.4 What the audit had already isolated as missing

By the time this session started, the audit had already narrowed the missing physics to a small and explicit set. To keep the minimal charged-lepton route alive, the model needed either a genuine non-axisymmetric (m=2) source capable of splitting the same-charge transverse rotor, or a compact/topological phase law that quantized the Berry flux, or some new first-order internal structure beyond the present conservative closure that effectively did one of those jobs. Put differently, the corrected stack already supplied the charge sign and the candidate mixed-sector rotor, but it did **not** yet supply the same-charge internal discretizer. The entire session can be read as an attempt to see whether that missing discretizer could be supplied by finite-throat physics, induced anisotropy, selective mixed-sector topology, or autonomous recirculation closure.

## 4. Corrected same-charge mixed-sector rotor

### 4.1 Reduced transverse model

Once the old circulation-based and static side-bit routes had failed, the only live same-charge corridor left in the corrected ontology was the mixed ((a-w)) sector. The reason is structural. The parent matter action is already first order in time, so a genuine collective-coordinate Berry/Magnus term is kinematically allowed from the start, and the full 4+1 ontology keeps the mixed channels
[
A_w,\qquad F_{\mu w},\qquad J^w,\qquad E_w,\qquad C_a=F_{aw}
]
alive outside the strict zero-mode brane Maxwell limit. That is exactly the sector in which a same-charge first-order internal law could exist without reducing to ordinary orbital circulation.

The minimal same-charge ansatz already isolated in the audit was built from the transverse odd pair ((b_x,b_y)), the first odd longitudinal mode, and a proposed mixed-twist sign (\tau=\pm1):
[
\psi(b_x,b_y,\tau)
==================

e^{i\theta_0}\Big[
\Psi_0
+b_x,u_x(\mathbf x)\phi_0(w)
+b_y,u_y(\mathbf x)\phi_0(w)
+i\tau,\zeta\bigl(-b_yu_x+b_xu_y\bigr)\partial_w\phi_1(w)
\Big].
]
At low energy this reduces to the transverse mixed-sector model
[
L_{\rm tr}
==========

\frac{M_b}{2}\dot b_i\dot b_i
+
\frac{\tau\nu_0}{2}\epsilon_{ij}b_i\dot b_j
-------------------------------------------

V_{\rm stat}(b_x^2+b_y^2),
]
and, if the radius ever pins,
[
b_x=\sqrt{s_*}\cos\phi,\qquad
b_y=\sqrt{s_*}\sin\phi,
]
so the system becomes an internal rotor,
[
L_{\rm rot}
===========

\frac{I}{2}\dot\phi^2+\tau\kappa_0\dot\phi,
\qquad
I=M_b s_*,
\qquad
\kappa_0=\frac{\nu_0 s_*}{2},
]
with spectrum
[
E_n^{(\tau)}=\frac{(n\hbar-\tau\kappa_0)^2}{2I}.
]
This was the first place in the audit where a same-charge internal sign and a genuine first-order internal phase law coexisted in one reduced model.

### 4.2 Berry curvature in corrected variables

The decisive point was that the same mixed ansatz produces a genuine Berry curvature,
[
F_{b_xb_y}
==========

2\hbar,\Im!\int d^4X,
\partial_{b_x}\psi^*,\partial_{b_y}\psi
=======================================

2\hbar,\tau,\zeta,\mathcal I,
]
with
[
\mathcal I \sim
\int d^3x,(u_x^2+u_y^2)
\int dw,\phi_0(w),\partial_w\phi_1(w).
]
So (\tau) is not an arbitrary label; it is literally the sign of the mixed-sector Berry coefficient. This is why the same-charge internal sign became plausible only in the mixed ((a-w)) sector and nowhere else in the corrected stack.

The Gaussian/Hermite unit test already in the notes made that concrete. Although odd Hermite modes do not couple directly to a strictly brane-centered source because (H_{2m+1}(0)=0), the mixed derivative overlap does not vanish:
[
\int_{-\infty}^{\infty}
e^{-w^2/\lambda^2},\phi_0(w),\partial_w\phi_1(w),dw
===================================================

\frac{\sqrt2}{\lambda}\neq0.
]
This was the first controlled demonstration that a hidden odd mode could still generate a nonzero first-order internal coefficient through the mixed sector even when it was invisible to the leading brane source.

### 4.3 Same-charge invariance at fixed electric branch

The ontology correction mattered here because it finally made the “same-charge” test clean. In the corrected stack,
[
q_*=\eta_Q e_*,
\qquad
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
]
with (\eta_Q) the electric branch sign and (n_\Gamma) the magnetic/vortical winding label. That means a candidate same-charge internal sign (\tau) can be studied while holding (\eta_Q) and (n_\Gamma) fixed, instead of silently changing electric charge through circulation or breathing geometry.

In the symbolic checks run this session, that bookkeeping reduced the same-charge condition to a very simple form. Holding (\eta_Q) fixed gives
[
\partial_b q_*=0,
]
while
[
\partial_b q_{\rm eff}
======================

-\frac{q_*}{2Z_{\rm int}^{3/2}},
\partial_b Z_{\rm int}.
]
So same-charge invariance at fixed electric branch becomes the geometric statement that the centered odd mixed deformation should not linearly change the localization thickness:
[
\partial_b Z_{\rm int}\big|*0=0.
]
This is exactly the corrected replacement for the dead (q\sim a^2\Gamma) logic. The conceptual result is that ((\eta_Q,n*\Gamma,\tau=+1)) and ((\eta_Q,n_\Gamma,\tau=-1)) are kinematically allowed same-charge, same-static-support candidates, at least to the order where the density correction remains (\tau)-even. 

### 4.4 Magnetic-moment audit and why minimal coupling gives (g=1)-type behavior

The magnetic audit was the next gate. Evaluating the parent minimal coupling in a weak external brane magnetic field gives a Zeeman shift
[
\delta H_B=-\mu_z^{(\tau)}B_z,
]
with
[
S_z^{(\tau)}=\tau,\hbar,\zeta,\mathcal I,s,
\qquad
\mu_z^{(\tau)}=\tau,\frac{q_*,\hbar,\zeta,\mathcal I}{2m_\psi},s,
\qquad
s=b_x^2+b_y^2.
]
So the same overlap that generated the Berry curvature also generates an intrinsic magnetic moment tied to the same same-charge internal sign. This was a genuine success: the mixed-sector route does produce a rest-frame intrinsic dipole law rather than reducing everything to orbital circulation. 

What it did **not** produce automatically was the Dirac factor. At this stage the relation is
[
\mu=\frac{q_*}{2m_\psi}S,
]
which is (g=1)-type relative to the parent GNLS mass parameter (m_\psi), not automatic Dirac (g=2). The audit also ruled out the tempting loophole of counting the mixed (v^wC_a) force term as a second independent Zeeman contribution; that would double-count the same minimal interaction. So the mixed-sector construction can reach an intrinsic magnetic moment of the right qualitative form, but it does not by itself close the Dirac magnetic package.

## 5. Why the original Gaussian/noncompact route stalled

### 5.1 No built-in half-integer Berry-flux quantizer

The original mixed-sector route stalled not because the Berry coefficient vanished, but because the Gaussian/noncompact localization scheme supplied no built-in discretizer for it. The frozen EM localization program uses
[
Z(w)=e^{-w^2/\lambda^2},
\qquad
f_n(w)=H_n(w/\lambda),
\qquad
m_n^2=\frac{2n}{\lambda^2},
]
which is a noncompact line with localized Hermite modes, not a compact internal circle or double-cover bundle. That means the Berry offset
[
\alpha\equiv \frac{\kappa_0}{\hbar}=\frac{\nu_0 s_*}{2\hbar}
]
remains a continuous overlap parameter rather than a topologically protected half-integer. This was one of the session’s most important negative results: the ontology supports, at best, a same-charge Berry rotor, but not yet a protected lepton-like doublet.

### 5.2 Why passive (P_{22}) support was insufficient

The second stall was static rather than topological. The 2PN operator already contains (P_{22c}) and (P_{22s}) as support channels, but the solved source data remain purely axisymmetric:
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
]
and the Family-1 wall/source law depends only on (\mu=\cos\theta). So the (P_{22}) channels are present only as passive support capacity, not as a driven anisotropy source. Their contribution reduces to
[
|Q_2|^2
=======

# (b_x^2-b_y^2)^2+(2b_xb_y)^2

# (b_x^2+b_y^2)^2

s^2,
]
which can soften or stiffen the transverse radius but cannot select an angle in the plane. That is why the static side/orientation route failed: the current freeze contains (P_{22}) support slots, but not the (-V_2\cos2(\phi-\phi_0)) term that would be needed to split the rotor into a genuine same-charge doublet.

### 5.3 What exactly was missing

By the time those two failures were clear, the missing physics had been narrowed to a very small set. The model needed either a genuine non-axisymmetric (m=2) source that could split the same-charge transverse rotor, or a compact/topological law that could quantize the Berry flux, or some new first-order internal structure beyond the present conservative closure that effectively did one of those jobs. In the corrected ontology, charge sign was already supplied; the missing structure was no longer “spin in general” but a **same-charge internal discretizer/quantizer**. That sharper statement is the real reason this session could make progress: it stopped treating the problem as vague and started treating it as the search for one very specific missing object.

## 6. Finite-throat reinterpretation

### 6.1 Geometric tube versus EM soft-wall mismatch

The first major reinterpretation of the session was to recognize a geometric scope mismatch in the old route. On the mechanics/particle side, the throat is treated as a real finite structure with radius (a), length (L), and a preferred aspect ratio (L/a\approx1.85). The 1PN bridge also explicitly interprets inertial mass as trapped wave energy and treats the defect as supported by at least one stable standing mode. On the EM/lepton side, however, the localization sector is written with a noncompact Gaussian soft wall (Z(w)) on the whole line. That means the same particle was being described as finite when geometry mattered and effectively infinite when the internal odd-mode/Berry problem was computed.

This did not by itself prove a new quantizer, but it did justify replacing the old Gaussian/noncompact picture with a finite-throat reinterpretation before declaring the Berry route “classically continuous.” In other words, before concluding that (\nu_0) was an unavoidably continuous overlap coefficient, we had to ask whether we were even integrating it over the correct geometric domain.

### 6.2 Bounded-throat replacement of the Gaussian/Hermite overlap

The first bounded-throat replacement treated the internal corridor as a finite interval (w\in[-L/2,L/2]) and replaced the Gaussian/Hermite basis by ordinary bounded modes. In the simplest Dirichlet test basis, the mixed overlap became
[
\mathcal I_L=\int_{-L/2}^{L/2}\phi_0,\partial_w\phi_1,dw=\frac{8}{3L},
]
while other standard closures gave similar (1/L)-type results. This already changed the interpretation of (\nu_0): it was no longer “an overlap in an infinite soft wall,” but a finite-throat spectral datum tied to the actual geometric length (L).

The more exact finite-throat version then used the 2PN cylinder / Neumann-bottom mouth branch. The exact mouth operator on that branch is
[
Z_{00}(\omega)=-\frac{\omega}{c_s}\tan!\left(\frac{\omega L}{c_s}\right),
]
with a half-shifted pole ladder
[
\omega_n^{\rm pole}=\frac{\pi c_s}{L}\left(n+\frac12\right),
]
so the natural throat basis is the D/N family
[
\chi_n(z)=\sqrt{\frac{2}{L}}\sin!\left(\frac{(n+\tfrac12)\pi z}{L}\right),
\qquad z\in[0,L].
]
In that exact finite-throat basis, the derivative-overlap matrix becomes
[
I_{mn}^{\rm DN}=\int_0^L \chi_m,\partial_z\chi_n,dz,
]
with an explicit closed form and, in particular,
[
I_{nn}^{\rm DN}=\frac1L.
]
That was the session’s cleanest replacement of the old Gaussian overlap by a genuine finite-throat throat-response object. 

### 6.3 What finite support fixes

Finite-throat reinterpretation fixed several real problems. It removed the artificial mismatch between a finite geometric throat and an infinite localization integral. It discretized the longitudinal basis and made the Berry coefficient mode-indexed rather than purely continuum-indexed. It also showed that a finite open throat can support a first-order mixed coefficient without depending on the special Gaussian parity trick: in the exact D/N basis, the derivative-overlap matrix is generically nonzero. Finally, it tied the same internal direction (w) that appears in the mixed-sector Berry route back to the same finite throat length (L) that already appears in the geometry and cavity side of the model.

Put differently, finite support made the mixed-sector rotor much more physically honest. It no longer lived in a mathematically permissive infinite line when the rest of the particle ontology lived in a finite defect.

### 6.4 What finite support does not fix

What finite support did **not** do was close the quantizer by itself. Even on the exact D/N branch, the rebuilt Berry coefficient remained of the form
[
\nu_0^{\rm DN}=2\hbar,\zeta,U_2,\mathcal I_{\rm DN},
]
and the Berry offset remained continuous:
[
\alpha(\theta)=\frac{\zeta U_2 s_*}{L},(2\sin^2\theta-3)
]
in the simplest normalized two-mode family. So finite support converted the old Gaussian continuum into a discrete mode basis, but it still did not force
[
\alpha\in\mathbb Z+\frac12.
]
Likewise, the open-throat half-shifted D/N pole ladder turned out to be only a **spectral** half-shift, not yet a fermionic or spinorial one. It is quarter-wave cavity physics, not automatic (4\pi) holonomy. That is why the finite-throat reinterpretation was a real advance, but not yet the finish line: it corrected the geometric scope error without yet supplying the missing same-charge discretizer.

## 7. Boundary-condition sweep on the finite throat

### 7.1 Neumann, Dirichlet, Robin, periodic, and antiperiodic closures

Once the Gaussian/noncompact localization scheme was recognized as too permissive, the next step was to replace the internal (w)-direction by a finite throat interval and test the mixed-sector overlap under explicit boundary conditions. The point of this sweep was not to commit immediately to one compactification, but to determine which closures preserve a nonzero mixed derivative overlap and whether any of them automatically quantize the Berry offset.

For a centered finite interval (w\in[-L/2,L/2]), the simplest representative closures gave the following results.

For an open-segment **Neumann**-type basis,
[
\phi_0=\frac{1}{\sqrt L},
\qquad
\phi_1=\sqrt{\frac{2}{L}}\sin!\left(\frac{\pi w}{L}\right),
]
the mixed overlap becomes
[
I_N=\int_{-L/2}^{L/2}\phi_0,\partial_w\phi_1,dw
===============================================

\frac{2\sqrt2}{L}.
]

For a **Dirichlet**-type basis,
[
\phi_0=\sqrt{\frac{2}{L}}\cos!\left(\frac{\pi w}{L}\right),
\qquad
\phi_1=\sqrt{\frac{2}{L}}\sin!\left(\frac{2\pi w}{L}\right),
]
the overlap is
[
I_D=\frac{8}{3L}.
]

For a **periodic** compact closure in the simplest even/odd real basis, the overlap vanishes:
[
I_P=0.
]

For an **antiperiodic** compact closure,
[
\phi_0=\sqrt{\frac{2}{L}}\cos!\left(\frac{\pi w}{L}\right),
\qquad
\phi_1=\sqrt{\frac{2}{L}}\sin!\left(\frac{\pi w}{L}\right),
]
the overlap is
[
I_A=\frac{\pi}{L}.
]

Finally, for a symmetric **Robin** family on the interval, with dimensionless boundary parameter
[
\eta=\frac{hL}{2},
]
the parity-preserving roots satisfy
[
x_e\tan x_e=\eta,
\qquad
x_o\cot x_o=-\eta,
]
and the mixed overlap takes the form
[
I_{\rm Robin}=\frac{F(\eta)}{L},
]
where (F(\eta)) is a continuous dimensionless function. On the lowest branch it runs between the Neumann and Dirichlet values, so the Robin family interpolates continuously between the open-segment cases rather than producing a new topological quantizer.

### 7.2 Which closures preserve a nonzero mixed overlap

The sweep immediately split the closures into two physically different classes.

The open-segment **Neumann**, **Dirichlet**, and **Robin** closures all preserve a nonzero mixed derivative overlap, so they keep the same-charge Berry corridor alive. The **antiperiodic** compact closure also preserves a nonzero overlap, and therefore survives as a compact candidate. By contrast, the simplest **periodic** compactification is actually hostile to the live corridor because in the natural real basis it kills the overlap altogether:
[
I_P=0.
]

This was already useful conceptually. It showed that “compactification” is not a magic word by itself. Some compact closures destroy the only live same-charge corridor, while some open-segment closures preserve it perfectly well.

### 7.3 Why boundedness alone does not force (\alpha\in\mathbb Z+\tfrac12)

The most important negative result of the boundary sweep was that finite support does not by itself quantize the Berry offset. In all successful finite-throat cases, the overlap takes the form
[
I_{\rm BC}=\frac{C_{\rm BC}}{L}
]
or, in the Robin family,
[
I_{\rm Robin}=\frac{F(\eta)}{L}
]
with (F(\eta)) continuous in the boundary parameter. If the mixed-sector coefficient inherits one power of this overlap, then
[
\nu_0 \sim C_\nu I_{\rm BC},
\qquad
\alpha=\frac{\nu_0 s_*}{2\hbar}
\sim
\frac{C_\nu s_*}{2\hbar L},C_{\rm BC}.
]
So finite geometry discretizes the **mode basis**, but it leaves the Berry offset continuous in the geometric and response data unless an additional law fixes the combination (C_\nu s_*/L).

This was the clean boundary-condition version of the earlier Gaussian/noncompact no-go. The old route stalled because the internal line was noncompact; the new sweep showed that finite support corrects that geometric mismatch, but still does not by itself produce the missing half-integer theorem.

### 7.4 Why antiperiodic closure was the strongest compact candidate

Among the compact candidates, the antiperiodic closure emerged as the strongest. Unlike the periodic compact case, it preserves a nonzero mixed overlap. It also shifts the longitudinal spectrum to the half-integer lattice,
[
k_n=\frac{(2n+1)\pi}{L},
]
which is exactly the kind of spectral half-shift one would want if a spinor-like quantizer were hiding nearby.

But the sweep also made clear what antiperiodicity does **not** do. It half-shifts the longitudinal momentum labels; it does not automatically force
[
\alpha\in\mathbb Z+\frac12.
]
So antiperiodic compactification survived as the strongest compact **candidate**, but not yet as the missing quantizer itself. The lesson of the whole sweep was therefore:

[
\text{finite throat} + \text{good boundary condition}
]
can make the same-charge Berry corridor physically honest,

but only an additional topological or holonomy law can make it spinorial.

## 8. Exact D/N mouth branch and finite-throat Berry rebuild

### 8.1 Exact 2PN mouth operator and half-shifted pole ladder

The next step was to stop working with generic bounded test functions and instead rebuild the mixed-sector corridor on the **exact** finite-throat branch already singled out by the 2PN mouth analysis. On that branch, the throat coordinate (z\in[0,L]) carries a mouth datum at (z=0) and a Neumann bottom at (z=L). Solving
[
\psi''+k^2\psi=0,
\qquad
k=\frac{\omega}{c_s},
]
with (\psi(0)=\psi_m) and (\psi'(L)=0) gives
[
\psi(z,\omega)=\psi_m,\frac{\cos(k(L-z))}{\cos(kL)}.
]
So the outward mouth derivative is
[
-\partial_z\psi(0,\omega)
=========================

-k\tan(kL),\psi_m,
]
which defines the exact mouth operator
[
Z_{00}(\omega)=-\frac{\omega}{c_s}\tan!\left(\frac{\omega L}{c_s}\right).
]

Its pole ladder is
[
\omega_n^{\rm pole}
===================

\frac{\pi c_s}{L}\left(n+\frac12\right),
]
while its zeros sit at
[
\omega_n^{\rm zero}
===================

\frac{\pi c_s n}{L}.
]
This immediately showed that the exact conservative throat branch already carries a half-shifted **spectral** structure. But the session also established that this half-shift is quarter-wave open-throat cavity physics, not yet a fermionic or spinorial identity law.

### 8.2 D/N overlap matrix (I_{mn}^{\rm DN})

The natural longitudinal basis on this D/N branch is
[
\chi_n(z)
=========

\sqrt{\frac{2}{L}}
\sin!\left(\frac{(n+\tfrac12)\pi z}{L}\right).
]
For this exact basis, the mixed derivative-overlap matrix is
[
I_{mn}^{\rm DN}
\equiv
\int_0^L \chi_m(z),\partial_z\chi_n(z),dz.
]
This turned out to have a simple closed form:
[
I_{mn}^{\rm DN}
===============

\begin{cases}
\dfrac{2n+1}{L(m+n+1)}, & m+n\ \text{even},[6pt]
\dfrac{2n+1}{L(m-n)}, & m+n\ \text{odd}.
\end{cases}
]

Several important qualitative facts follow immediately.

First,
[
I_{nn}^{\rm DN}=\frac1L
]
for every mode (n), so the open-throat D/N geometry automatically supports a nonzero first-order mixed coefficient even on the diagonal. Second, the matrix is generically nonzero off-diagonal as well; for the lowest two-mode block,
[
I^{\rm DN}_{(0,1)}
==================

\frac1L
\begin{pmatrix}
1 & -3\
1 & 1
\end{pmatrix}.
]
This meant that the finite-throat corridor no longer depended on the special Gaussian “centered even times odd derivative” trick. The D/N throat geometry itself builds in a strong mixed derivative structure.

### 8.3 Rebuilt (\nu_0^{\rm DN}) and (\alpha(\theta))

With the exact D/N basis in hand, the mixed-sector Berry coefficient could be rebuilt in a genuinely finite-throat form. Writing the longitudinal mixed profiles as
[
\Phi(z)=\sum_m a_m\chi_m(z),
\qquad
X(z)=\sum_n b_n\chi_n(z),
]
the same mixed-sector algebra as before gives
[
F_{b_xb_y}
==========

2\hbar,\tau,\zeta,U_2,\mathcal I_{\rm DN},
\qquad
\mathcal I_{\rm DN}=\sum_{m,n}a_m b_n I_{mn}^{\rm DN},
]
where
[
U_2=\int d^3x,(u_x^2+u_y^2).
]
So the rebuilt Berry coefficient is
[
\nu_0^{\rm DN}=2\hbar,\zeta,U_2,\mathcal I_{\rm DN}.
]

At the level of fixed mode pairs this already gives discrete (1/L)-weighted outcomes. For example,
[
I_{01}^{\rm DN}=-\frac{3}{L}
]
implies
[
\nu_0^{(0,1)}=-\frac{6\hbar\zeta U_2}{L}.
]
But the decisive test was not a fixed pair. It was a normalized two-mode family,
[
\Phi_\theta
===========

\cos\theta,\chi_0+\sin\theta,\chi_1,
\qquad
X_\theta
========

-\sin\theta,\chi_0+\cos\theta,\chi_1.
]
For this family the mixed overlap becomes
[
\mathcal I_\theta
=================

\frac{2\sin^2\theta-3}{L},
]
so
[
\nu_0(\theta)
=============

\frac{2\hbar\zeta U_2}{L},(2\sin^2\theta-3),
]
and therefore
[
\alpha(\theta)
==============

# \frac{\nu_0(\theta)s_*}{2\hbar}

\frac{\zeta U_2 s_*}{L},(2\sin^2\theta-3).
]

### 8.4 Result: finite throat improves ontology but leaves (\alpha) continuous

This was the strongest finite-throat reconstruction result of the whole session. The D/N rebuild improved the ontology in every physically relevant sense: it put the mixed-sector coefficient on the actual finite throat, made it explicitly mode-indexed, and showed that the exact open-throat branch naturally supports the first-order mixed structure the lepton corridor requires. It also clarified that the half-shifted D/N pole ladder is a genuine feature of the exact mouth operator.

But it still did **not** quantize the Berry offset. Even after replacing the infinite Gaussian line by the exact D/N throat basis, the Berry offset remained a continuous function of the mode-mixing angle and the continuous response combination (\zeta U_2 s_*/L). So the outcome of the rebuild was:

[
\text{finite throat fixes the geometry of the corridor,}
]
but

[
\text{finite throat alone does not supply the missing same-charge quantizer.}
]

This sharpened the problem again. The old “continuous (\nu_0)” complaint could no longer be blamed on the Gaussian/noncompact mismatch. After the D/N rebuild, whatever quantizer was still missing really was missing from the current ontology.

## 9. Dynamic pinning of the transverse branch

### 9.1 Why (s=b_x^2+b_y^2) is not the same as throat radius (a)

Once the D/N rebuild made clear that (\alpha) was still continuous, the next obvious loophole was the transverse radius (s_*). The earlier mixed-sector derivation had treated
[
s=b_x^2+b_y^2
]
as a free low-energy radius. Another AI correctly pointed out that leaving (s_*) free was too kinematic. But the session also clarified a crucial bookkeeping point: (s) is **not** the same object as the geometric throat radius (a). The radius (a) is the actual particle/throat geometry variable already present in the mechanical and 1PN/2PN closures. By contrast, (s) is the squared amplitude of the transverse odd (|m|=1) branch that lives inside the same-charge mixed-sector reduction.

So “fix the defect core radius” is not yet the right statement for the lepton problem. What the Berry corridor needs is not merely the geometry equilibrium for (a), but a genuine pinning law for the **transverse odd branch** itself.

### 9.2 Generic Landau pinning law (V_{\rm eff}(s))

The most conservative file-grounded way to parameterize that pinning is the even-under-(b\to-b) Landau reduction of the odd branch:
[
V_{\rm eff}(s)
==============

\frac12\kappa_s s
+
\frac14 u_{\rm eff}s^2
+
\frac16 v_{\rm eff}s^3,
\qquad
s=b_x^2+b_y^2.
]
The stationarity condition is
[
\frac{dV_{\rm eff}}{ds}
=======================

\frac12\bigl(\kappa_s+u_{\rm eff}s+v_{\rm eff}s^2\bigr)=0,
]
so the nonzero stationary radii are
[
s_*^\pm
=======

\frac{-u_{\rm eff}\pm\sqrt{u_{\rm eff}^2-4\kappa_s v_{\rm eff}}}{2v_{\rm eff}}.
]
In the quartic pitchfork limit (v_{\rm eff}\to0), this reduces to the familiar
[
s_*=-\frac{\kappa_s}{u_{\rm eff}}.
]

This was the right replacement for the too-crude “centrifugal core balance” idea. The session did not find a frozen theorem equating (s) directly with the ordinary defect radius. What it found was the correct generic **odd-branch** pinning law that any same-charge transverse rotor must satisfy if it is to remain inside the current conservative closure.

### 9.3 Possible (n_\Gamma)-dependence of the static branch

Because the corrected ontology now keeps (n_\Gamma) as the magnetic/vortical label rather than the electric branch itself, it also becomes natural to let the static odd-branch stiffness depend on (n_\Gamma) in an even way. The safest proposed insertion tested this by writing
[
\kappa_s(n_\Gamma)=\kappa_0+c_\Gamma n_\Gamma^2.
]
This preserves the corrected bookkeeping: the winding label can affect the static energy of the branch without redefining electric charge. With that insertion, the pinned radius becomes
[
s_*(n_\Gamma)
=============

\frac{-u_{\rm eff}+\sqrt{u_{\rm eff}^2-4v_{\rm eff}(\kappa_0+c_\Gamma n_\Gamma^2)}}{2v_{\rm eff}}.
]

This was only a proposal, not a frozen theorem of the current files. But it was the least aggressive way to see how the magnetic/vortical sector might enter the static branch without violating the corrected ontology.

### 9.4 Result: pinning is necessary but not sufficient

The key outcome of the pinning test was that a dynamically locked (s_*) is indeed necessary for the Berry corridor, but still not sufficient to quantize it. Substituting the pinned radius into the D/N rebuilt Berry offset leaves
[
\alpha_*(n_\Gamma,\theta)
=========================

-\frac{U_2\zeta}{2Lv_{\rm eff}}
\Bigl(
u_{\rm eff}
-----------

\sqrt{u_{\rm eff}^2-4v_{\rm eff}(\kappa_0+c_\Gamma n_\Gamma^2)}
\Bigr)
(2\sin^2\theta-3).
]
So even after the branch is dynamically pinned, the Berry offset still depends on continuous response data and the continuous mode-mixing angle.

This meant the dynamic-pinning loophole could not close the quantizer by itself. The honest conclusion was therefore:

* if the odd branch never pins with (s_*>0), the same-charge Berry corridor dies dynamically;
* if it does pin, the corridor survives dynamically, but the missing same-charge discretizer is still missing.

That result was important because it prevented another false shortcut. It showed that “fixing the radius” is a real part of the problem, but not the final one. The quantizer still had to come from anisotropy, topology, or recirculation closure beyond the pinning law itself.

## 10. Route A revived: area-preserving (P_{22}) splitter

### 10.1 Constant-area squish and exact ellipse construction

After the dynamic-pinning analysis, the next serious question was whether the missing (P_{22}) splitter could appear **induced** rather than already frozen. The key physical correction was to impose a constant-area deformation of the throat cross-section, motivated by the idea that the particle’s mass-supporting inflow should not change under a same-charge spin reorganization. If the throat is squeezed along one transverse axis and elongated along the perpendicular one while preserving area exactly, the natural parametrization is
[
a_\parallel=a,e^{-\xi},
\qquad
a_\perp=a,e^{\xi},
\qquad
A=\pi a_\parallel a_\perp=\pi a^2.
]
So the deformation is not an arbitrary ellipse approximation; it is the exact constant-area one-parameter family.

In polar form, with principal axis (\phi_0), the boundary is
[
r(\vartheta)
============

\frac{a}{\sqrt{e^{2\xi}\cos^2(\vartheta-\phi_0)+e^{-2\xi}\sin^2(\vartheta-\phi_0)}}.
]
At small (\xi), this is already an explicit (m=2) distortion. The importance of this construction is that it finally supplies a geometric mechanism for a real in-plane quadrupole without changing the total cross-sectional area.

### 10.2 Quadrupole tensor and induced (-V_2\cos2(\phi-\phi_0))

Once the exact ellipse is written down, the second-moment tensor of the cross-section immediately produces a traceless in-plane quadrupole. In the body frame, the second moments are proportional to (a_\parallel^3a_\perp) and (a_\parallel a_\perp^3), so after rotation into the lab frame the traceless quadrupole components become
[
Q_c\equiv Q_{xx}-Q_{yy}
=======================

-\frac{\pi a^4}{2}\cos2\phi,\sinh2\xi,
]
[
Q_s\equiv 2Q_{xy}
=================

-\frac{\pi a^4}{2}\sin2\phi,\sinh2\xi.
]
So the constant-area ellipse is not merely “sort of anisotropic.” It is a genuine real (P_{22}) object with amplitude proportional to (\sinh2\xi).

Coupling this quadrupole to an external traceless transverse squish tensor (\Sigma_{ij}) with axis (\phi_0) gives the leading symmetry-allowed anisotropy energy
[
\delta E_{\rm aniso}
====================

# -\lambda_{\rm cpl},\Sigma_{ij}Q_{ij}

-\frac{\pi}{4}\lambda_{\rm cpl}\Sigma_0 a^4\sinh2\xi,
\cos2(\phi-\phi_0).
]
This is exactly the static splitter the audit had identified as missing:
[
\delta E_{\rm aniso}=-V_2\cos2(\phi-\phi_0),
\qquad
V_2=\frac{\pi}{4}\lambda_{\rm cpl}\Sigma_0 a^4\sinh2\xi.
]
So the session’s symbolic work showed that an area-preserving thermodynamic squish really can generate the missing (P_{22}) term. The important point is not the particular radial scaling of (V_2), but the fact that an honest (\cos2\phi) splitter emerges from the geometry as soon as the induced deformation is forced to preserve area.

### 10.3 Why this rescues Route A only as an induced anisotropy mechanism

This revived Route A, but only in a very specific way. The frozen 2PN notes already said there is no file-derived static (P_{22}) drive in the isolated conservative particle: the support operator contains the (P_{22c}) and (P_{22s}) channels, but the solved source vector remains axisymmetric, and the Family-1 wall/source law stays axisymmetric as well. So the current conservative one-particle freeze still does **not** contain the anisotropy source required to split the transverse rotor statically.

What the area-preserving construction did was show that a (P_{22}) source can arise **induced** by an external partner, a finite-size bound-state response, or some other environment that exerts a traceless transverse squish while preserving the total cross-sectional area. In other words, the static side-bit rescue is still dead for the isolated particle, but Route A becomes alive again as an **interaction-induced anisotropy mechanism**. That distinction matters for the lepton program: it means the session found a real route to a two-well same-charge rotor, but not yet a free-particle same-charge spinor theorem.

### 10.4 Why it creates a two-well rotor but does not quantize (\alpha)

Injecting this splitter into the mixed-sector rotor gives
[
H
=

## \frac{1}{2I}\left(-i\hbar\partial_\phi-\tau\hbar\alpha\right)^2

V_2\cos2(\phi-\phi_0).
]
For (V_2>0), the continuous rotor angle is no longer flat. The static minima sit at
[
\phi=\phi_0,
\qquad
\phi=\phi_0+\pi,
]
so the rotor is broken from a continuous (O(2)) family to a genuine two-well system. This was the first clear rescue of Route A in the whole session.

But the tunnel splitting between the two wells still depends on the Berry phase. In the deep-well limit the two semiclassical tunneling paths interfere with phases (e^{\pm i\pi\alpha}), giving
[
t_{\rm eff}=2t_0\cos(\pi\alpha).
]
So the same induced (P_{22}) term that creates the doublet does **not** force
[
\alpha\in\mathbb Z+\frac12.
]
It only makes the half-integer points physically distinguished, because those are the points where
[
t_{\rm eff}=0.
]
That was the exact result of the combined rotor-plus-(P_{22}) test: the induced anisotropy discretizes the angle (\phi), but it does not quantize the Berry flux. The problem of fixing (\alpha) therefore survived the Route A rescue and had to be attacked separately.

## 11. Spinor/topological stress tests that failed

### 11.1 (2\pi) mouth torque no-go

One of the most tempting ideas tested during the session was that a literal (2\pi) physical rotation of the mouth might induce a (-1) spinor-like sign on the internal throat mode. The clean stress test used the most general closed-throat identification for a separated mode
[
\Psi(w,\varphi)=e^{iqw}e^{im\varphi},
\qquad
\Psi(w+L,\varphi)=e^{i\delta}\Psi(w,\varphi+\beta),
]
with (\beta) the geometric twist after one traversal and (\delta) an optional bundle phase. This gives
[
q_{n,m}=\frac{2\pi n+m\beta+\delta}{L}.
]
Now apply a literal (2\pi) mouth rotation by (\beta\to\beta+2\pi). The result is
[
q_{n,m}(\beta+2\pi)=q_{n+m,m}(\beta),
]
so the effect is just a relabeling of the longitudinal level (n). Equivalently, the angular factor picks up
[
e^{im2\pi}=1
]
for ordinary integer mouth harmonics. So a classical (2\pi) mouth tumble cannot by itself produce the (-1) spinor holonomy we were after.

This failure was important because it killed a seductive but incorrect intuition: dynamic torsion of a scalar throat mode is not the missing source of (720^\circ) behavior. A (2\pi) loop can only become nontrivial if some additional (\mathbb Z_2) monodromy is already built into the internal bundle or topology.

### 11.2 Frame-matching no-go on the open throat

The next stress test asked whether the puncture orientation sign (\eta_Q) could generate the missing same-charge sign through frame matching between the 3D mouth and the 4D bulk. The strongest interval-level endpoint law of that kind is
[
\Psi(L,\varphi)=\eta_Q,e^{i\delta},\Psi(0,\varphi+\beta),
]
which, for a separated mode (\Psi=u(w)e^{im\varphi}), implies
[
u(L)=e^{i\Theta_{\rm end}}u(0),
\qquad
\Theta_{\rm end}=m\beta+\delta+\frac{\pi}{2}(1-\eta_Q).
]
But on an **open** throat interval this endpoint phase is gauge-removable. Defining
[
\widetilde u(w)=e^{-i\Theta_{\rm end}w/L}u(w)
]
gives
[
\widetilde u(0)=\widetilde u(L).
]
So any sign produced by interval-level frame matching can be removed by rephasing on the interval. It is not an intrinsic (\mathbb Z_2) monodromy.

The second no-go in the same test was even more decisive:
[
\frac{\partial \Theta_{\rm end}}{\partial \tau}=0.
]
So frame matching based on (\eta_Q) does not generate the missing same-charge sign (\tau) at all. It only reorganizes the already-frozen electric branch. This clarified something that mattered for the rest of the session: (\eta_Q) is the electric branch, not the same-charge spinor quantizer.

### 11.3 Cross-cap no-go for the actual mixed vector-odd sector

After the (2\pi) and frame-matching failures, the strongest remaining topological proposal was a cross-cap or antipodal reflection closure. The corresponding gluing law is
[
\Psi(w+L,\varphi)=p_w,e^{i\delta},\Psi(w,\varphi+\pi),
\qquad
p_w=\pm1,
]
with monodromy
[
M=p_w,e^{i(m\pi+\delta)}.
]
This is already stronger than the previous tests because it includes explicit non-orientable/antipodal structure. But it still failed for the **actual** mixed carrier.

The reason is that the relevant mixed mode is not a pure odd scalar. In the lepton ansatz, it is the first odd longitudinal factor multiplied by the transverse in-plane vector pair ((u_x,u_y)), so it behaves effectively like a (w)-odd, (m=1) object. For that sector, pure geometric cross-cap gluing gives
[
M_{\rm mixed}=(-1)\times(-1)=+1.
]
The (-1) from the odd (w)-parity is canceled by the (-1) from the antipodal (m=1) factor. So the actual mixed mode remains **periodic**, not antiperiodic, under pure cross-cap geometry.

A pure odd scalar with (m=0) would get the desired (-1), but that is not the same object as the same-charge mixed sector we actually need. Adding an extra antiperiodic bundle phase (\delta=\pi) can force the mixed sector to flip, but then the even scalar branch flips too unless the bundle is made sector-specific. So the cross-cap test showed that geometry alone does not selectively quantize the same-charge mixed mode.

It also showed one more thing: even the strongest compact cross-cap gluing only shifts the compact momentum lattice. It changes the compact spectrum from periodic to antiperiodic form, but it does not by itself force
[
\alpha\in\mathbb Z+\frac12.
]
So cross-cap topology was not the missing quantizer either.

### 11.4 What each failure taught us

These no-go tests were not wasted effort. Each one removed a different false shortcut and sharpened the search.

The (2\pi) mouth-torque test showed that classical dynamic torsion is not enough; the missing sign must be a static internal monodromy or bundle feature, not a simple rotational phase of a scalar mode.

The frame-matching test showed that interval-level endpoint signs are not intrinsic on the open throat and that (\eta_Q)-based constructions target the electric branch, not the same-charge mixed sector.

The cross-cap test showed that even genuine geometric/non-orientable gluing acts on the **total** parity (p_w(-1)^m), not on (w)-oddness alone, so the actual mixed vector-odd carrier evades the sign flip we wanted.

Taken together, these failures taught us that any viable same-charge quantizer must be **selective**: it must live in the mixed (\tau)-sector itself, not in the electric branch, not in a universal scalar bundle, and not in a generic geometrical parity law. That lesson is exactly what led to the next route.

## 12. Selective (\tau)-subbundle route

### 12.1 Why the quantizer must live in the mixed sector, not in (\eta_Q)

By the time the topological no-go tests were finished, one structural point had become unavoidable: any viable same-charge quantizer must live in the mixed sector itself. The corrected ontology already fixes
[
\eta_Q
]
as the electric branch sign, and the lepton audit had already shown that (\eta_Q)-based side/orientation switches are opposite-charge labels, not same-charge spin labels. The same-charge mixed-twist sign
[
\tau=\pm1
]
appears only in the mixed ((a-w)) sector, where it is literally the sign of the Berry coefficient. The even scalar/support bundle is (\tau)-even and cannot distinguish the two same-charge branches.

So the right bundle decomposition is not “put everything on one twisted bundle.” It is
[
\mathcal E_{\rm tot}=\mathcal L_{\rm even}\oplus \mathcal E_\tau,
]
where the even electric branch (\mathcal L_{\rm even}) remains trivial and all potentially nontrivial gluing is confined to the mixed same-charge subspace (\mathcal E_\tau). This was the cleanest way to avoid contaminating the electric charge branch with the missing spinor-like sign.

### 12.2 Diagonal mixed-subbundle gluing

On the mixed (\tau)-doublet,
[
\chi(\phi)=
\begin{pmatrix}
\chi_+(\phi)\
\chi_-(\phi)
\end{pmatrix},
]
the current finite-throat Berry rotor is already diagonal in the (\tau) basis. Its connection on the rotor circle is
[
A_\phi=\alpha,\sigma_3,
]
so its Wilson loop is
[
W(\alpha)
=========

# e^{i2\pi\alpha\sigma_3}

\begin{pmatrix}
e^{i2\pi\alpha}&0\
0&e^{-i2\pi\alpha}
\end{pmatrix}.
]

The first selective extension tested was an **independent diagonal gluing law**
[
U_\tau=\mathrm{diag}(e^{i\theta_+},e^{i\theta_-}).
]
After gauging away the connection on the interval, the compact monodromy is just
[
M_{\rm diag}=U_\tau W.
]
This yields the compact spectrum
[
E_{n,+}
=======

\frac{\hbar^2}{2I}
\left(n+\frac{\theta_+}{2\pi}+\alpha\right)^2,
\qquad
E_{n,-}
=======

\frac{\hbar^2}{2I}
\left(n+\frac{\theta_-}{2\pi}-\alpha\right)^2.
]
So diagonal selective gluing does not quantize (\alpha). It merely shifts the compact momentum labels of the two (\tau)-branches.

This was the first clear indication that a selective (\tau)-subbundle is not enough by itself. If its gluing is independent of the Berry connection, it only relabels the compact ladder.

### 12.3 Central-sign gluing (U_\tau=-I)

The strongest same-charge sign flip one can impose on the mixed subspace without introducing off-diagonal transport is the central sign
[
U_\tau=-I_2.
]
As an **independent** gluing law, this gives
[
M_{-I}=-W,
]
so the compact spectra become
[
E_{n,+}^{(-I)}
==============

\frac{\hbar^2}{2I}
\left(n+\frac12+\alpha\right)^2,
\qquad
E_{n,-}^{(-I)}
==============

\frac{\hbar^2}{2I}
\left(n+\frac12-\alpha\right)^2.
]
Again, (\alpha) remains continuous; the central sign simply produces the usual antiperiodic half-shift of the compact ladder.

But this same test also revealed the one way the current Abelian rotor **can** force a half-integer theorem. If the physical mixed-sector gluing is not extra boundary data but is literally the Berry Wilson loop itself, then requiring
[
W(\alpha)=-I_2
]
gives
[
e^{i2\pi\alpha}=-1
\quad\Longrightarrow\quad
\alpha\in\mathbb Z+\frac12.
]
So the central-sign route only succeeds if the missing (\tau)-subbundle law is **identified with the current Berry holonomy itself**, not imposed independently.

That became the crucial distinction for everything that followed.

### 12.4 Exchange gluing and why it needs a new non-Abelian structure

The next natural idea was a genuine (\mathbb Z_2) exchange law,
[
U_\tau=\sigma_x,
]
which would swap (\tau=+) and (\tau=-) after one internal loop. This is attractive because it would give a real same-charge two-state exchange structure rather than just a common sign.

However, with the current diagonal Wilson loop (W(\alpha)), the combined monodromy
[
M_{\rm swap}=\sigma_x W
]
has eigenvalues
[
\mathrm{spec}(M_{\rm swap})={+1,-1},
]
independent of (\alpha). So an exchange subbundle can produce a same-charge doublet structure, but it does **not** quantize the Berry offset; it removes (\alpha) from the compact monodromy instead.

More importantly, because the current Berry connection is Abelian and diagonal in (\tau), its Wilson loop can never itself become (\sigma_x). Any true exchange holonomy would therefore require a new off-diagonal non-Abelian transport law with explicit (\sigma_1) or (\sigma_2) structure beyond the current mixed rotor
[
L_{\rm rot}
===========

\frac{I}{2}\dot\phi^2+\tau\kappa_0\dot\phi.
]
So the exchange route was not ruled out forever, but it was ruled out **inside the present reduced closure**.

### 12.5 Key result: independent selective gluing only relabels compact ladders

The overall conclusion of the (\tau)-subbundle analysis was very sharp. A selective mixed subbundle is exactly the right place to look, because it keeps the even electric branch trivial and localizes all nontrivial holonomy inside the same-charge mixed sector. But as long as its gluing is independent of the Berry connection, it only **relabels** the compact momentum ladders. Diagonal gluing shifts the two (\tau)-branches; central sign gluing gives the usual antiperiodic shift; exchange gluing produces a (\mathbb Z_2) doublet but erases (\alpha) from the monodromy.

So the route succeeds only in one very specific case:

[
\text{the selective mixed-sector gluing must be the Berry Wilson loop itself.}
]

That is why the next stage of the session moved away from arbitrary bundle assignments and into recirculation/plumbing holonomy. The problem was no longer “find any selective gluing law,” but

[
\text{derive }U_\tau=W(\alpha)
]
as a real physical loop law of the mixed sector. That is the exact bridge to the next section.

## 13. Recirculation/plumbing holonomy

### 13.1 Current diagonal loop law (U_{\rm loop}=e^{i\phi_0}W(\alpha))

Once the selective (\tau)-subbundle route had isolated the mixed sector as the only viable home of the missing same-charge quantizer, the next question was whether the model’s unresolved recirculation/plumbing physics could generate the required gluing law there. The crucial simplification was that the current closure does not support a general non-Abelian (2\times2) transport law on the (\tau)-doublet. Instead, everything already frozen in the ontology points to a very narrow diagonal structure.

On the mixed (\tau)-subbundle, the Berry part of the loop law is already known from the current rotor:
[
A_\phi=\alpha,\sigma_3,
\qquad
W(\alpha)=e^{i2\pi\alpha\sigma_3}.
]
Any unresolved scalar recirculation/support transport can only enter as a common (\tau)-even phase. So the strongest loop operator compatible with the current closure is
[
U_{\rm loop}
============

# \exp!\big(i\Phi_0 I_2+i\Phi_3\sigma_3\big)

e^{i\phi_0}W(\alpha),
]
where (\phi_0) is the common scalar recirculation phase and (\Phi_3=2\pi\alpha) is the same-charge Berry phase. Explicitly,
[
U_{\rm loop}
============

\begin{pmatrix}
e^{i(\phi_0+2\pi\alpha)}&0\
0&e^{i(\phi_0-2\pi\alpha)}
\end{pmatrix}.
]

This was the key structural reduction of the entire endgame. Once written in this form, the problem of the lepton discretizer ceased to be “find some topology that might help” and became the much narrower question: can the plumbing/recirculation law force this diagonal loop operator to close with the central sign on the mixed subbundle?

### 13.2 What the existing ledgers can feed: common phase versus Berry phase

The reason the loop operator takes this restricted form is that the project’s existing support/leakage ledgers are all naturally (\tau)-even. The notes already contain exact scalar diagnostics of through-flow, back-pressure, leakage, and mixed work, such as
[
\Phi_w(t),\qquad
\Delta h(t),\qquad
S_{\rm leak}^{(s)},\qquad
J^wE_w.
]
Within the present conservative closure, these quantities can renormalize or constrain the **common scalar phase** (\phi_0), but they do not provide the missing (\tau)-flip matrices (\sigma_1) or (\sigma_2). So they can feed the (I_2) part of the loop law, not a new off-diagonal exchange structure.

That mattered because it killed one more tempting shortcut. It meant the current recirculation problem is not “maybe a hidden non-Abelian bundle appears if one looks hard enough.” Inside the present closure it is a strictly diagonal holonomy problem:
[
U_{\rm loop}=e^{i\phi_0}e^{i2\pi\alpha\sigma_3}.
]
So the final unknown could no longer be buried in general bundle language. It had a very specific form: a common scalar phase (\phi_0) on top of the already-known Berry Wilson loop.

### 13.3 Central-sign theorem and why it only implies (\alpha\in\tfrac12\mathbb Z) unless (\phi_0) cancels

The natural target of the recirculation route was the central mixed-sector sign
[
U_{\rm loop}=-I_2.
]
Imposing that condition on the current diagonal loop law gives
[
e^{i(\phi_0+2\pi\alpha)}=-1,
\qquad
e^{i(\phi_0-2\pi\alpha)}=-1.
]
So there must exist integers (p,q) such that
[
\phi_0+2\pi\alpha=(2p+1)\pi,
\qquad
\phi_0-2\pi\alpha=(2q+1)\pi.
]
Solving these simultaneously gives
[
\alpha=\frac{p-q}{2}\in \frac12\mathbb Z,
\qquad
\phi_0=\pi(p+q+1).
]

This is one of the most important refinements of the session. The central-sign condition by itself does **not** yet force the lepton target
[
\alpha\in\mathbb Z+\frac12.
]
It only implies the weaker statement
[
\alpha\in \frac12\mathbb Z.
]
In other words, if the common scalar phase remains free, both integer and half-integer Berry flux remain compatible with the same central sign on the mixed subbundle.

This was the precise point where the old “find a (-1) sign and the problem is solved” intuition failed. The session showed that even a successful mixed-sector (-I_2) holonomy is still too weak unless the common scalar phase is simultaneously forced into the trivial class.

### 13.4 Scalar-phase cancellation as the real bottleneck

That observation immediately identified the true bottleneck:
[
\text{the real missing theorem is not just }U_{\rm loop}=-I_2,
\text{ but }U_{\rm loop}=-I_2\text{ with }\phi_0\equiv0\pmod{2\pi}.
]
If
[
\phi_0\equiv0\pmod{2\pi},
]
then the central-sign condition collapses to
[
W(\alpha)=-I_2
\quad\Longrightarrow\quad
\alpha\in\mathbb Z+\frac12.
]
But if, for example,
[
\phi_0\equiv\pi\pmod{2\pi},
]
then the same central-sign law instead gives
[
\alpha\in\mathbb Z.
]

So by the end of the recirculation/plumbing analysis, the whole lepton corridor had been compressed to one narrow issue: whether the common scalar recirculation phase on the mixed subbundle is forced to cancel. This was the exact point where the session pivoted from general holonomy language to standing-wave and autonomous-soliton closure.

## 14. Standing-wave closure and autonomous soliton closure

### 14.1 Stable standing-wave round-trip condition (R_{\rm rt}=1)

The first attempt to kill the scalar phase (\phi_0) was to treat the common support channel as a standing wave on a full internal loop. Writing the scalar round-trip factor as
[
R_{\rm rt}=r_0,r_L,e^{2ikL},
]
with
[
r_0=\rho_0e^{i\delta_0},\qquad
r_L=\rho_Le^{i\delta_L},
]
the full-loop scalar phase is
[
\phi_0=\delta_0+\delta_L+2kL.
]
If the support wave is a leakage-free trapped standing wave, then after one full scalar round trip it must reproduce itself:
[
R_{\rm rt}=1.
]
This is the usual cavity/eigenmode condition. It implies
[
\rho_0\rho_L=1,
\qquad
\phi_0=2\pi N,
]
so
[
e^{i\phi_0}=1,
\qquad
\phi_0\equiv0\pmod{2\pi}.
]

This was a crucial positive result. It showed that the desired scalar-phase cancellation is not some arbitrary aesthetic choice. It is exactly what one gets if the common support channel is a true leakage-free standing-wave eigenmode.

### 14.2 Why generic steady state is too weak

However, the next thermodynamic test showed that one must be careful not to overstate that result. The handoff suggestion that “thermodynamic steady state” by itself forces
[
R_{\rm rt}=1
]
turned out to be too strong.

The simplest one-cycle amplitude map is
[
A_{n+1}=\Lambda A_n + D,
\qquad
\Lambda=\rho e^{i\phi_0},
]
where (\Lambda) is the intrinsic round-trip coefficient of the support loop and (D) is a replenishment/pump term over one cycle. The exact fixed point is
[
A_*=\frac{D}{1-\Lambda}.
]
So a non-decaying steady state with (A_*= \text{const}) does **not** imply
[
\Lambda=1
]
unless
[
D=0.
]

This is the key loophole the session had to close. A pump-balanced or leakage-balanced steady state can keep the amplitude constant while the intrinsic round-trip coefficient remains different from 1. So “thermodynamic steady state” is not enough by itself to prove the standing-wave closure.

The same distinction appears in the continuous envelope model
[
\dot A=\left[\frac{\gamma_{\rm in}-\gamma_{\rm leak}}{2}+i\omega\right]A.
]
Steady state only forces
[
\gamma_{\rm in}=\gamma_{\rm leak},
]
which fixes the modulus. It does not stop the phase from advancing. To make the wave reproduce itself after one loop of duration (T), one still needs the resonance condition
[
\omega T=2\pi N.
]

So the session sharpened the physical alternatives into two genuinely different closures:
[
\text{steady amplitude}
]
versus
[
\text{self-reproducing standing-wave eigenmode}.
]
Only the second one kills (\phi_0).

### 14.3 Difference between driven steady state and autonomous eigenmode

This distinction then became the core conceptual divide in the endgame.

A **driven-dissipative steady state** allows external replenishment or pump balance. In discrete form,
[
A_{n+1}=\Lambda A_n + D,
]
and in continuous form,
[
\dot A=\left[\frac{\gamma_{\rm in}-\gamma_{\rm leak}}{2}+i\omega\right]A.
]
Such a state can be perfectly stable in amplitude while carrying a nontrivial phase drift around the loop. It is a fixed point, but not necessarily an eigenmode.

An **autonomous eigenmode**, by contrast, is self-reproducing with no external replenishment:
[
D=0,
\qquad
A_{n+1}=\Lambda A_n,
\qquad
A_{n+1}=A_n\neq0.
]
This is stronger than “no pumping” and stronger than “constant modulus.” It means the internal support wave closes on itself exactly after one full loop.

This distinction turned out to be the final conceptual bottleneck between a plausible spinor-like closure and an actual conditional theorem.

### 14.4 Autonomous soliton closure (D=0) and (A_{n+1}=A_n\neq0)

The final handoff of the session proposed reading the isolated free particle as an autonomous GNLS soliton. The symbolic test showed that this works, but only in the strong eigenmode sense just described.

Injecting the soliton condition (D=0) into the cycle map gives
[
A_{n+1}=\Lambda A_n.
]
By itself this does not force (\Lambda=1). It only removes external pumping.

If one asks merely for non-decay of the modulus,
[
|A_{n+1}|=|A_n|,
]
then one gets
[
|\Lambda|=1,
]
but the phase (\phi_0) remains free.

Only the full self-reproducing eigenmode condition
[
A_{n+1}=A_n\neq0
]
forces
[
\Lambda=1.
]
Since
[
\Lambda=\rho e^{i\phi_0},
]
this gives
[
\rho=1,
\qquad
\phi_0=2\pi N.
]

That was the strongest positive result of the closure chain. It showed that the common scalar phase is indeed killed by the autonomous soliton/eigenmode closure, but not by weaker notions of “stability” alone.

## 15. Final conditional theorem

### 15.1 Derivation of (\phi_0\equiv0\ (\mathrm{mod}\ 2\pi)) under autonomous eigenmode closure

At this point the endgame becomes very short. The current mixed-subspace loop law had already been reduced to
[
U_{\rm loop}=e^{i\phi_0}W(\alpha),
\qquad
W(\alpha)=e^{i2\pi\alpha\sigma_3}.
]
The autonomous-soliton closure now supplies
[
D=0,
\qquad
A_{n+1}=A_n\neq0
\quad\Longrightarrow\quad
\Lambda=1
\quad\Longrightarrow\quad
\phi_0=2\pi N.
]
So on the mixed subbundle,
[
e^{i\phi_0}=1,
\qquad
U_{\rm loop}=W(\alpha).
]

This is the exact point where the old ambiguity disappears. The central-sign law no longer has to compete with an arbitrary scalar phase. The mixed-sector holonomy is now nothing but the Berry Wilson loop itself.

### 15.2 Final lock (U_{\rm loop}=-I_2 \Rightarrow \alpha\in\mathbb Z+\tfrac12)

Once the common scalar phase has been eliminated, the selective mixed-sector central-sign target becomes
[
U_{\rm loop}=-I_2
\quad\Longrightarrow\quad
W(\alpha)=-I_2.
]
But
[
W(\alpha)=
\begin{pmatrix}
e^{i2\pi\alpha} & 0\
0 & e^{-i2\pi\alpha}
\end{pmatrix},
]
so
[
W(\alpha)=-I_2
\quad\Longrightarrow\quad
e^{i2\pi\alpha}=-1.
]
Therefore
[
\alpha\in\mathbb Z+\frac12.
]

Equivalently,
[
W!\left(k+\frac12\right)
========================

\begin{pmatrix}
-1&0\
0&-1
\end{pmatrix},
\qquad
W(k)=
\begin{pmatrix}
1&0\
0&1
\end{pmatrix}.
]
So the half-integer same-charge spin discretizer is not merely plausible under this closure; it is mathematically unavoidable.

This is the strongest positive result reached in the entire session.

### 15.3 Exact status: theorem conditional on autonomous self-reproducing soliton closure

The final status must still be stated carefully.

What is now proved is a **conditional theorem**:

[
\boxed{
D=0
;+;
A_{n+1}=A_n\neq0
;+;
U_{\rm loop}=-I_2
;\Longrightarrow;
\phi_0\equiv0\ (\mathrm{mod}\ 2\pi),
\quad
\alpha\in\mathbb Z+\tfrac12.
}
]

So if an isolated defect is an autonomous self-reproducing GNLS soliton, and if the selective mixed-sector loop closes with the central sign, then the half-integer same-charge Berry-flux quantizer follows rigorously.

What is **not yet** upgraded to frozen theorem status is the physical step that an isolated particle really is such an autonomous eigenmode. The project-side handoff remains explicit that the current 1PN/2PN program is solved only within a declared conservative closure hierarchy and still leaves spin couplings, dissipative leakage/radiation reaction, and the fully solved moving-throat PDE unresolved.

So the clean final verdict of this session is:

[
\boxed{
\text{The minimal charged-lepton route is provisionally viable.}
}
]

More precisely, it is viable up to one remaining nontrivial closure assumption: that the isolated particle’s support channel is a self-reproducing autonomous eigenmode rather than merely a pump-free state with constant norm. If that final physical closure can be derived from the actual moving-throat/support dynamics, the half-integer same-charge spin discretizer becomes a genuine file-grounded theorem. If it cannot, then the present result remains a rigorous but conditional closure theorem rather than a completed particle derivation.

## 16. What is now established versus still assumed

### 16.1 What is established by this session’s symbolic work

This session established a much sharper internal theorem structure than existed at the start. First, the corrected ontology cleanly separates the electric branch ((\eta_Q,q_*,q_{\rm eff})) from the magnetic/vortical label (n_\Gamma), which means a same-charge internal label can now be discussed without silently redefining electric charge. Second, the mixed ((a-w)) sector remains the only live same-charge corridor: the first odd longitudinal mode survives in mixed derivative overlaps, the reduced transverse odd pair can support a first-order Berry term, and if the transverse radius pins the resulting dynamics are those of an internal rotor rather than an orbital relabeling. Third, the exact finite-throat D/N rebuild showed that the old Gaussian/noncompact treatment was too permissive geometrically, but that the same-charge Berry corridor survives on the exact finite throat as a genuine mode-indexed throat response. Fourth, an induced area-preserving ellipse supplies a real (P_{22}) anisotropy mechanism, so Route A is no longer dead in an interaction-induced setting even though it remains absent in the isolated conservative freeze.

At the symbolic level, the strongest concrete result is the following conditional chain. If the mixed same-charge subspace carries the central sign on its loop holonomy, and if the common scalar recirculation/support phase cancels, then the Berry Wilson loop itself forces
[
W(\alpha)=-I_2
\quad\Longrightarrow\quad
\alpha\in\mathbb Z+\frac12.
]
That statement was not available at session start. It is the first exact formulation of a half-integer same-charge discretizer compatible with the current mixed-sector rotor. In that sense, the session did not merely produce another candidate mechanism; it produced the clearest conditional endpoint of the entire lepton audit so far.

### 16.2 Strong conditional assumptions still needed

The route, however, closes only under a narrow set of additional conditions. The first is that the mixed-sector gluing must be the **central** sign on the selective (\tau)-subbundle, not merely an independent diagonal phase choice and not a noncentral exchange law. Independent diagonal gluing only shifts the compact momentum ladders, while exchange gluing would require a new non-Abelian transport structure that is not present in the current reduced closure. The second is that the common scalar loop phase must cancel, so that the mixed-sector loop law reduces from
[
U_{\rm loop}=e^{i\phi_0}W(\alpha)
]
to
[
U_{\rm loop}=W(\alpha).
]
The third is that the isolated particle must be interpreted in the strong autonomous-eigenmode sense, i.e.
[
D=0,\qquad A_{n+1}=A_n\neq0,
]
not merely as a pump-free state with constant modulus. Only under that stronger closure does the common scalar phase become trivial mod (2\pi). These assumptions are physically coherent inside the session’s logic, but they are still assumptions rather than frozen file-derived theorems.

### 16.3 What is still not file-grounded theorem status

The remaining gap is now very specific. The current project-side handoffs are explicit that the conservative 1PN/2PN program is solved only within a declared closure hierarchy and still leaves spin couplings, dissipative leakage / radiation reaction, and the fully solved moving-throat PDE unresolved. That means the session did **not** prove, directly from the frozen files alone, that an isolated particle is an autonomous self-reproducing eigenmode. It also did not prove a free-particle, file-derived static (P_{22}) splitter inside the conservative freeze, and it did not derive a genuine non-Abelian (\tau)-exchange transport law. So the half-integer discretizer is now best described as a rigorous **conditional theorem** rather than a finished file-grounded particle theorem.

The honest current status is therefore:

[
\boxed{
\text{minimal charged-lepton route: provisionally viable, with one remaining nontrivial closure assumption.}
}
]

That assumption is no longer vague. It is precisely the autonomous self-reproducing eigenmode closure of the isolated support channel, together with the central-sign closure on the mixed subbundle. If those two pieces can be derived rather than imposed, the session’s half-integer lock becomes a true theorem of the model. If not, the route remains an elegant but conditional possibility.

## 17. Highest-leverage next step for the next session

### 17.1 Derive autonomous eigenmode closure from the actual support/leakage dynamics

The next session should not reopen the whole lepton audit from scratch. It should start from the end-state reached here and ask one sharply framed question:

[
\boxed{
\text{Can the actual support/leakage dynamics of an isolated particle derive the autonomous eigenmode closure?}
}
]

In the language of the one-cycle amplitude map, that means testing whether the isolated particle’s support channel is genuinely a self-reproducing trapped eigenmode
[
D=0,\qquad A_{n+1}=A_n\neq0,
]
rather than merely a pump-balanced or leakage-balanced steady state. The exact ledgers already singled out by the session — through-flow, back-pressure, scalar leakage, and mixed work — are the right place to do this, because they are the only file-grounded scalar channels that can feed the common phase sector. So the highest-leverage next move is not “derive more spin structure,” but “derive the closure that kills (\phi_0).”

### 17.2 What would count as success

A successful next step would show, using the project’s own support/leakage language, that an isolated stable defect does not merely have constant norm or balanced flux, but actually satisfies the self-reproducing loop condition. In practical terms, that means deriving a file-grounded reason why the support channel reproduces itself after one full internal cycle, so that the round-trip coefficient satisfies
[
\Lambda=1
\quad\Longrightarrow\quad
e^{i\phi_0}=1.
]
Once that is in hand, the session’s existing central-sign result immediately finishes the discretizer:
[
U_{\rm loop}=-I_2
\quad\Longrightarrow\quad
\alpha\in\mathbb Z+\frac12.
]
That would upgrade the present closure from a conditional theorem to the first really file-grounded same-charge half-flux theorem in the program.

### 17.3 What failure would mean for the minimal charged-lepton program

Failure is also now cleanly interpretable. If the support/leakage dynamics show only a pump-free state with constant modulus, or a leakage-balanced steady state with residual phase drift, then the common scalar phase (\phi_0) remains free and the central-sign route collapses back to the weaker statement
[
\alpha\in\frac12\mathbb Z,
]
with no rigid half-integer lock. In that case the minimal charged-lepton route would still have an interesting mixed-sector rotor and possibly an induced (P_{22}) doublet mechanism, but it would still lack a file-grounded same-charge spinor quantizer. At that point the program would either have to find a new non-Abelian mixed transport law or accept that the present ontology cannot yet finish the charged-lepton target. This is exactly the sort of falsification-first endpoint the handoff methodology asked for from the beginning.
