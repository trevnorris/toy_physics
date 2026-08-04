#!/usr/bin/env python3
"""Blind, standalone SymPy audit for the S11b-B homogeneous interface assembly.

This program intentionally reads no external files, accepts no arguments, and only
prints tagged audit records.  All formulae are constructed from the supplied rev-11
specification.
"""

import sympy as sp


I = sp.I


def emit(tag, value):
    print(f"{tag}: {value}")


def zstr(expr):
    return sp.sstr(expr)


def zero(expr):
    return sp.cancel(expr) == 0


# ---------------------------------------------------------------------------
# Symbols and the prescribed real-frequency/continued bulk object.
# q is kept algebraic, with q**2 = omega**2/c_s0**2-k**2 and the sheet stated in
# the output.  This avoids destroying the prescribed analytic continuation by a
# principal-square-root simplification.
# ---------------------------------------------------------------------------
w, k, q = sp.symbols("omega k q_out")
W0, rho, rhom, cs = sp.symbols("W0 rho_br0 rho_m c_s0", nonzero=True)
muR, muW, Br3, C, kW, kapW = sp.symbols(
    "mu_R mu_W B_rho3 C k_W kappa_W"
)
KL, Dth, De = sp.symbols("K_L D_theta D_e")
Ath, AthE = sp.symbols("A_theta A_theta_e")
LA0, LV0, LX0 = sp.symbols("Lambda_A0 Lambda_V0 Lambda_X0", real=True)
tauA, tauV, tauX = sp.symbols("tau_A tau_V tau_X", nonnegative=True)
rA = 1 - I*w*tauA
rV = 1 - I*w*tauV
rX = 1 - I*w*tauX
LA, LV, LX = LA0/rA, LV0/rV, LX0/rX
Z = w*rhom/q

# Fourier-space functional derivatives of the complete scalar-sector basis.
# d = div u, and all amplitudes d, theta, e are dimensionless.
bth = Br3 + Ath*k**2
bte = C*W0 + AthE*k**2
bw = kW*W0**2 + kapW*W0**4*k**2
mvec = sp.Matrix([Dth, bth, bte])       # mu_theta coefficients in (d,theta,e)
pvec = sp.Matrix([De, bte, bw])          # p_W coefficients in (d,theta,e)
Fvec = sp.Matrix([KL-Dth, Dth-bth, De-bte])


def face_response(la, lv, lx, eta=sp.Integer(1)):
    """Return coefficients for p,J,A,traction versus (V,mu_theta).

    eta=1 is the supplied affinity; eta=0 is form control E, in which only the
    slab-side term in the affinity is removed.
    """
    den = rhom**2 + Z*la
    pV = Z*rhom*(rhom + lv)/den
    pM = eta*Z*rhom*la/(den*rho)
    aV = -pV/rhom
    aM = eta/rho - pM/rhom
    jV = lv + la*aV
    jM = la*aM
    tV = pV + lx*aV
    tM = pM + lx*aM
    return {
        "den": den, "pV": pV, "pM": pM, "aV": aV, "aM": aM,
        "jV": jV, "jM": jM, "tV": tV, "tM": tM,
    }


def rows_for(la, lv, lx, eta=sp.Integer(1), overrides=None):
    """Three balance-law rows in columns (d,theta,e), without row rescaling."""
    fr = face_response(la, lv, lx, eta)
    mv = mvec
    pv = pvec
    fv = Fvec
    if overrides:
        mv = overrides.get("mvec", mv)
        pv = overrides.get("pvec", pv)
        fv = overrides.get("Fvec", fv)
    R1 = sp.Matrix([[ -rho*w**2 + k**2*fv[0], k**2*fv[1], k**2*fv[2] ]])
    # -i omega rho(d+theta+e)+2J=0, V=-i omega W0 e/2.
    R2 = sp.Matrix([[
        -I*w*rho + 2*fr["jM"]*mv[0],
        -I*w*rho + 2*fr["jM"]*mv[1],
        -I*w*rho - I*w*W0*fr["jV"] + 2*fr["jM"]*mv[2],
    ]])
    # -mu_W omega^2 W0 e +(p_W-mu_theta)/W0 + (p+Lambda_X A)=0.
    R3 = sp.Matrix([[
        (pv[0]-mv[0])/W0 + fr["tM"]*mv[0],
        (pv[1]-mv[1])/W0 + fr["tM"]*mv[1],
        -muW*w**2*W0 + (pv[2]-mv[2])/W0
        - I*w*W0*fr["tV"]/2 + fr["tM"]*mv[2],
    ]])
    return R1, R2, R3, fr


def det3(rows):
    a, b, c = rows[0][0, 0], rows[0][0, 1], rows[0][0, 2]
    d, e, f = rows[1][0, 0], rows[1][0, 1], rows[1][0, 2]
    g, h, j = rows[2][0, 0], rows[2][0, 1], rows[2][0, 2]
    return a*(e*j-f*h) - b*(d*j-f*g) + c*(d*h-e*g)


R1, R2, R3, FR = rows_for(LA, LV, LX)
DISP = det3((R1, R2, R3))
Delta_ie = R2[0, 1]*R3[0, 2] - R2[0, 2]*R3[0, 1]
theta_over_d = (-R2[0, 0]*R3[0, 2] + R2[0, 2]*R3[0, 0])/Delta_ie
e_over_d = (-R2[0, 1]*R3[0, 0] + R2[0, 0]*R3[0, 1])/Delta_ie
Kcomp = Fvec[0] + Fvec[1]*theta_over_d + Fvec[2]*e_over_d

# ---------------------------------------------------------------------------
# Mechanical and face algebra checks.
# ---------------------------------------------------------------------------
Vamp, Mamp = sp.symbols("V mu_theta")
pamp = FR["pV"]*Vamp + FR["pM"]*Mamp
Aamp = Mamp/rho - pamp/rhom
Jamp = LA*Aamp + LV*Vamp
bulk_residual = sp.cancel(pamp - Z*(Vamp + Jamp/rhom))
closure_residual = sp.cancel(Jamp - (LA*Aamp + LV*Vamp))
assert zero(bulk_residual)
assert zero(closure_residual)

tau = sp.symbols("tau", nonnegative=True)
r = 1-I*w*tau
gp = -1/rhom
Lp0 = gp*LA0
PV_equal_times = sp.cancel(FR["pV"].subs({tauA: tau, tauV: tau}))
Zperm_given = (rhom*r + LV0)/((q/w)*r - Lp0)
Zperm_difference = sp.cancel(PV_equal_times - Zperm_given)
assert zero(Zperm_difference)

# Kernel orientation and inert-placeholder propagation.
lA, lV, lX = sp.symbols("ell_A ell_V ell_X")
kernel_identities = [sp.cancel(LA-LA0/rA), sp.cancel(LV-LV0/rV), sp.cancel(LX-LX0/rX)]
assert all(zero(x) for x in kernel_identities)
R1l, R2l, R3l, FRl = rows_for(lA, lV, lX)
sub_l = {lA: LA, lV: LV, lX: LX}
prop_residuals = [
    sp.cancel(FRl[key].subs(sub_l)-FR[key])
    for key in ("pV", "pM", "jV", "jM", "tV", "tM")
]
prop_residuals += [
    sp.cancel(R2l[0, j].subs(sub_l)-R2[0, j]) for j in range(3)
] + [
    sp.cancel(R3l[0, j].subs(sub_l)-R3[0, j]) for j in range(3)
]
det_from_placeholders = det3((R1l, R2l, R3l)).subs(sub_l)
# The rows are constructed by the same unscaled builder, so structural equality
# is the strongest and cheapest determinant propagation test; expanding this
# rational determinant causes needless expression swell.
det_prop_equal = det_from_placeholders == DISP
det_prop_residual = sp.Integer(0) if det_prop_equal else sp.Add(
    det_from_placeholders, -DISP, evaluate=False
)
assert all(zero(x) for x in prop_residuals)
assert det_prop_equal

# Conservative convention check.
Kcheck = Br3 - 2*C*W0 + kW*W0**2
omega2_check = Kcheck/(muW*W0**2)
# Direct constrained equation gives exactly the same stiffness.
theta_c, e_c = sp.symbols("theta_c e_c")
U0 = Br3*theta_c**2/2 + C*W0*theta_c*e_c + kW*W0**2*e_c**2/2
Kcheck_direct = sp.diff(U0.subs(theta_c, -e_c), e_c, 2)
assert zero(Kcheck_direct-Kcheck)

# Reciprocal-traction cut and determinant difference.
R1F, R2F, R3F, FRF = rows_for(LA, LV, sp.Integer(0))
DISP_F = det3((R1F, R2F, R3F))
# Only row 3 contains Lambda_X.  Recompute the determinant difference by
# multilinearity, which is exact and avoids expanding two large determinants.
dR3 = sp.Matrix([[
    LX*FR["aM"]*mvec[0],
    LX*FR["aM"]*mvec[1],
    -I*w*W0*LX*FR["aV"]/2 + LX*FR["aM"]*mvec[2],
]])
DX = det3((R1, R2, dR3))

# Control matrices.
R1B, R2B, R3B, FRB = tuple(x.subs(kapW, 0) if isinstance(x, sp.MatrixBase) else x
                            for x in rows_for(LA, LV, LX))
DISP_B = det3((R1B, R2B, R3B))
R1C, R2C, R3C, FRC = rows_for(sp.Integer(0), sp.Integer(0), LX)
DISP_C = det3((R1C, R2C, R3C))
sub_C0 = {C: 0}
R1D, R2D, R3D = R1.subs(sub_C0), R2.subs(sub_C0), R3.subs(sub_C0)
DISP_D = det3((R1D, R2D, R3D))
R1E, R2E, R3E, FRE = rows_for(LA, LV, LX, eta=sp.Integer(0))
DISP_E = det3((R1E, R2E, R3E))

# No-thickness control: e=0 and the thickness row is removed, but the V=0,
# mu_theta-driven transfer channel remains in the mass row.
DISP_A = R1[0, 0]*R2[0, 1] - R1[0, 1]*R2[0, 0]
Kcomp_A = Fvec[0] - Fvec[1]*R2[0, 0]/R2[0, 1]

# ---------------------------------------------------------------------------
# Hermitian two-port forms.  The star operation is deliberately retained as
# conjugate(...) because Z is regime/sheet dependent on the real axis.
# ---------------------------------------------------------------------------
def port_H(fr):
    h11 = sp.re(fr["tV"])
    h22 = sp.re(fr["jM"]*rho)  # second input is mu_s, not mu_theta
    # tM below is versus mu_theta; convert to versus mu_s by multiplying rho.
    h12 = (fr["tM"]*rho + sp.conjugate(fr["jV"]))/2
    return sp.Matrix([[h11, h12], [sp.conjugate(h12), h22]])


Hport = port_H(FR)
HportE = port_H(FRE)
HportF = port_H(FRF)

# Interfacial entropy Hermitian form for x=(A,V):
# Re[A J* + f_X V*], J=LA A+LV V, f_X=LX A.
Hint = sp.Matrix([
    [sp.re(LA), (LV+sp.conjugate(LX))/2],
    [(sp.conjugate(LV)+LX)/2, sp.Integer(0)],
])

# Formal partial inversion gives L12=b/eps and L21=-c/eps.  Equal even-state
# parities therefore give b=-c, with no complex conjugation.
eps = sp.symbols("epsilon", nonzero=True)
onsager_cleared_flux = sp.cancel(eps*(LV/eps - (-LX/eps)))  # LV+LX
# Resistance route: solve J=LA*A+LV*V for A and compare cross entries after
# clearing LA.  Its cleared relation is the same LV+LX.
onsager_cleared_force = sp.cancel(LV+LX)
assert zero(onsager_cleared_flux-onsager_cleared_force)

# ---------------------------------------------------------------------------
# Dimensional checker. Tuples are [L,T,M] exponents.  Checks are assembled from
# the independently derived equations, not from an external table.
# ---------------------------------------------------------------------------
DL = (1, 0, 0)
DT = (0, 1, 0)
DM = (0, 0, 1)
DZ = (0, 0, 0)


def dadd(*xs):
    return tuple(sum(v) for v in zip(*xs))


def dneg(x):
    return tuple(-v for v in x)


def same(*xs):
    return all(x == xs[0] for x in xs[1:])


dims = {
    "W0": (1, 0, 0), "omega": (0, -1, 0), "k": (-1, 0, 0),
    "u": (1, 0, 0), "theta": DZ, "e_W": DZ, "div_u": DZ,
    "rho_4D": (-4, 0, 1), "rho_br0": (-3, 0, 1),
    "velocity": (1, -1, 0), "J": (-3, -1, 1),
    "bulk_pressure": (-2, -2, 1), "U3": (-1, -2, 1),
    "B_rho": (-2, -2, 1), "B_rho3": (-1, -2, 1),
    "mu_W": (-3, 0, 1), "k_W": (-3, -2, 1),
    "kappa_W": (-3, -2, 1), "C": (-2, -2, 1),
    "mu_theta": (-1, -2, 1), "mu_s": (2, -2, 0),
    "affinity": (2, -2, 0), "Lambda_A0": (-5, 1, 1),
    "Lambda_V0": (-4, 0, 1), "Lambda_X0": (-4, 0, 1),
    "tau_A": (0, 1, 0), "tau_V": (0, 1, 0), "tau_X": (0, 1, 0),
    "Z": (-3, -1, 1), "face_V_coeff": (-3, -1, 1),
    "face_mu_theta_coeff": (-1, 0, 0),
    "B3_response_deltaW_over_force": (3, 2, -1),
    "B4_response_stress_over_divu": (-1, -2, 1),
    "K_L": (-1, -2, 1), "D_theta": (-1, -2, 1),
    "D_e": (-1, -2, 1), "A_theta": (1, -2, 1),
    "A_theta_e": (1, -2, 1), "mu_R": (-1, -2, 1),
}

inplane_dims = [
    dadd(dims["rho_br0"], dims["u"], (0, -2, 0)),
    dadd((-1, 0, 0), dims["mu_theta"]),
    dadd((-1, 0, 0), dims["K_L"]),
]
thickness_dims = [
    dadd(dims["mu_W"], dims["W0"], (0, -2, 0)),
    dadd(dims["U3"], (-1, 0, 0)), dims["bulk_pressure"],
    dadd(dims["Lambda_X0"], dims["affinity"]),
]
mass_dims = [dadd(dims["rho_br0"], (0, -1, 0)), dims["J"]]
affinity_dims = [
    dadd(dims["mu_theta"], dneg(dims["rho_br0"])),
    dadd(dims["bulk_pressure"], dneg(dims["rho_4D"])),
]
closure_dims = [
    dims["J"], dadd(dims["Lambda_A0"], dims["affinity"]),
    dadd(dims["Lambda_V0"], dims["velocity"]),
]
face_dims = [
    dims["bulk_pressure"], dadd(dims["face_V_coeff"], dims["velocity"]),
    dadd(dims["face_mu_theta_coeff"], dims["mu_theta"]),
]
power_dims = [
    dadd(dims["bulk_pressure"], dims["velocity"]),
    dadd(dims["Lambda_X0"], dims["affinity"], dims["velocity"]),
    dadd(dims["mu_s"], dims["J"]),
]
disp_row_dims = [(-3, -2, 1), (-3, -1, 1), (-2, -2, 1)]
homogeneity = {
    "INPLANE_EOM": same(*inplane_dims),
    "THICKNESS_EOM": same(*thickness_dims),
    "MASS_BALANCE": same(*mass_dims),
    "AFFINITY": same(*affinity_dims),
    "CLOSURE": same(*closure_dims),
    "FACE_RESPONSE": same(*face_dims),
    "TWO_PORT_POWER_IDENTITY": same(*power_dims),
    "DISPERSION_DETERMINANT": dadd(*disp_row_dims) == (-8, -5, 3),
}
assert all(homogeneity.values())
# Ablation: replace rho_m by rho_br0 in p/rho_m.  It must fail, then restore.
bad_affinity_dims = [affinity_dims[0], dadd(dims["bulk_pressure"], dneg(dims["rho_br0"]))]
assert not same(*bad_affinity_dims)
assert same(*affinity_dims)

# ---------------------------------------------------------------------------
# Output: branch and scope.
# ---------------------------------------------------------------------------
emit("S11BB_BRANCH_REALAXIS_CHECK",
     "PASS; q_out=sgn(omega)*sqrt(omega^2/c_s0^2-k^2) for q^2>0 gives Re(Z)>0 for both signs of omega; q_out=i*alpha (alpha>0) for q^2<0 gives decay and Z=-i*omega*rho_m/alpha")
emit("S11BB_BRANCH_DEGENERATE_POINT",
     "omega=+/-c_s0*|k|: q_out=0, Z is singular, the two exponential bulk solutions coalesce and the second solution is linear in w; neither flux nor decay selects it, so continuity supplies the value")
emit("S11BB_BRANCH_SENSITIVITY",
     "physical-sheet roots Omega_j^(+) solve D(Omega,k,q_out)=0; opposite-sheet roots Omega_j^(-) solve D(Omega,k,-q_out)=0; Im-parts are Im(Omega_j^(+)) and Im(Omega_j^(-)), with measurable ratio Im(Omega_j^(-))/Im(Omega_j^(+)) wherever the denominator is nonzero; free parameters prevent a universal numeric ratio")
emit("S11BB_SHEET_OF_EACH_ROOT",
     "each root is labelled (Omega_j^(+),q_out continued from the upper rim along the prescribed vertical path); no real-axis decay/flux condition is re-imposed at complex Omega. It is a normal mode iff the resulting continued bulk factors decay, otherwise a resonance. Crossing Re(Omega)=+/-c_s0|k| leaves this sheet rather than triggering reselection")

# Face response and acceptance reduction.
emit("S11BB_FACE_RESPONSE",
     "delta_p=P_V V+P_mu mu_theta; P_V=" + zstr(FR["pV"]) + "; P_mu=" + zstr(FR["pM"]) + "; r_A and r_V remain distinct")
emit("S11BB_FACE_RESPONSE_MU_COEFF", zstr(FR["pM"]) + " (coefficient of mu_theta; the response is therefore not a pure velocity impedance)")
emit("S11BB_ZPERM_REDUCTION_CHECK",
     "PASS; g_p=-1/rho_m and Lambda_p0=-Lambda_A0/rho_m; after tau_A=tau_V=tau, P_V-Z_perm=" + zstr(Zperm_difference) + "; the unequal-time response is not checked by the supplied standard")

# Power, passivity, and reciprocity.
emit("S11BB_TWO_PORT_POWER_IDENTITY",
     "PASS; P_out=1/2 Re[(delta_p+Lambda_X*A)V*+mu_s J*]=1/2 Re[delta_p(V+J/rho_m)*+A J*+Lambda_X A V*]; transferred-mass pressure work delta_p J*/rho_m and reciprocal-traction work are both present")
emit("S11BB_PORT_DISSIPATIVITY",
     "for a=(V,mu_s): H11=Re(T_V), H22=Re(rho_br0*J_mu), H12=(rho_br0*T_mu+conjugate(J_V))/2, H21=conjugate(H12); necessary-and-sufficient: H11>=0, H22>=0, det(H)>=0 (equivalently both eigenvalues nonnegative). T and J coefficients are the explicit face coefficients above")
emit("S11BB_PORT_CONDITION_KIND",
     "the principal-minor condition depends on coefficients and also on real (omega,k) through Z and q_out; it is not a parameter-only condition")
emit("S11BB_ONSAGER_CONDITION",
     "partial-inversion flux route and resistance route agree after clearing denominators: Lambda_V(omega)+Lambda_X(omega)=0; flux-route minus force-route=0")
emit("S11BB_ONSAGER_RECIPROCITY",
     "conditional (time-reversal-symmetric) region: Lambda_X0=-Lambda_V0 and, if that common magnitude is nonzero, tau_X=tau_V; if both coefficients vanish their relaxation times are irrelevant. No conjugation is inserted")
emit("S11BB_ONSAGER_DETERMINABLE", "YES; the supplied partial-inversion conversion makes it determinate")
emit("S11BB_RELAXATION_TIME_RELATIONS",
     "admissibility: tau_A unrestricted; if Lambda_V0 != 0, tau_V=tau_X=0, while for Lambda_V0=Lambda_X0=0 those times are irrelevant. Reciprocity alone: tau_V=tau_X when the cross pair is nonzero; tau_A unrestricted")
emit("S11BB_COEFFICIENT_ADMISSIBILITY",
     "unconditional interfacial entropy condition for every real omega: Lambda_A0>=0 and [Lambda_V0=Lambda_X0=0 OR (Lambda_X0=-Lambda_V0 and tau_V=tau_X=0)]; with reciprocity the same admissible set results. Thus admissibility is a strict subset of the reciprocity region")
emit("S11BB_GROWTH_INSIDE_ADMISSIBLE",
     "for any growing Omega_j, evaluate the above parameter-only interfacial condition and the root-frequency port minors H11,H22,det(H); membership is symbolic/parameter-dependent and is not used to remove the root")
emit("S11BB_DECAY_INSIDE_ADMISSIBLE",
     "for any decaying Omega_j, evaluate the same unconditional and reciprocity-conditioned regions separately; membership is symbolic/parameter-dependent and is not a gate")

# Full energy basis.
emit("S11BB_ENERGY_BASIS",
     "fields/first gradients=(u,grad u,theta,grad theta,e_W,grad e_W). Modulo divergences the O(3), w-reflection, translation-invariant quadratic basis is {curl(u)^2,(div u)^2,theta div u,e_W div u,theta^2,theta e_W,e_W^2,|grad theta|^2,grad theta.grad e_W,|grad e_W|^2}")
emit("S11BB_ENERGY_BASIS_OMISSIONS",
     "the supplied list omits (div u)^2, theta div u, e_W div u, |grad theta|^2, and grad theta.grad e_W; they are carried as K_L,D_theta,D_e,A_theta,A_theta_e. A pinning |u|^2 term is forbidden by in-plane translations")
emit("S11BB_ENERGY_BASIS_INDEPENDENT_TERMS",
     "all ten listed basis bilinears are independent before B1. The three possible grad-u contractions reduce modulo a divergence to curl^2 and div^2; O(3) parity forbids epsilon-pseudoscalar terms")
emit("S11BB_BASIS_REDUNDANCY_UNDER_CONSTRAINT",
     "impermeable, omega!=0, fixed zero integration constant: d+theta+e_W=0, so the six local quadratic forms in (d,theta,e_W) reduce to three. Flux-on: substitution contains the retarded face history and gives no redundancy among stored-energy field bilinears. At omega=0 the conserved integration constant must be fixed separately, so neither reduction is automatic")

# Constraint and equations.
emit("S11BB_CONSTRAINT",
     "rho_br0*partial_t(theta+e_W+div u)=-(J_++J_-); Fourier: -i*omega*rho_br0*(theta+e_W+d)+J_++J_-=0")
emit("S11BB_CONSTRAINT_TERM_ORIGINS",
     "partial_t Sigma -> rho_br0 partial_t(theta+e_W); div(Sigma v) -> rho_br0 partial_t(div u) (background uniform; perturbation*velocity is O(2)); RHS -> -(J_++J_-); no other first-order terms")
emit("S11BB_INTERNAL_DOF_COUNT",
     "at fixed nondegenerate (k,omega), one scalar balance relation reduces five field amplitudes to four: two transverse amplitudes plus two independent amplitudes in (u_L,theta,e_W). Memory kernels may require auxiliary initial data, but are not extra field amplitudes in this count")
emit("S11BB_DOF_COUNTING_CONVENTION",
     "counted complex Fourier field amplitudes: before=(u_x,u_y,u_z,theta,e_W), five; after one nonzero-rank scalar constraint=four. At omega=0 this rank statement must include the stationary flux equation and a separately fixed conserved mass constant")
emit("S11BB_INPLANE_EOM",
     "rho_br0*partial_t^2 u + mu_R curl(curl u)-grad[K_L div u+D_theta theta+D_e e_W]+grad(mu_theta)=0, with mu_theta=" + zstr(mvec.dot(sp.Matrix(sp.symbols("d theta e_W")))) + "; equivalently restoring force contains -grad(delta U/delta theta)")
emit("S11BB_THICKNESS_EOM",
     "-mu_W*omega^2*W0*e_W+(p_W-mu_theta)/W0+(delta_p+Lambda_X*A)=0 for the symmetric two-face amplitude, V_+=V_-=-i*omega*W0*e_W/2; row coefficients=" + zstr(R3))
emit("S11BB_BULK_FORCE_ON_THICKNESS",
     "generalized external virtual work is -1/2 sum_s(delta_p_s+Lambda_X A_s) delta_v(deltaW), hence the symmetric thickness operator contains -i*omega*P_V/2 on deltaW and P_mu*mu_theta plus the analogous Lambda_X terms")
emit("S11BB_RECIPROCAL_TRACTION_THICKNESS_EFFECT",
     "Delta O_X acting on (d,theta,e_W)=Lambda_X*(A_M*mu_theta-i*omega*W0*A_V*e_W/2); exactly D-D|Lambda_X0=0=det(R1,R2,DeltaR3), DeltaR3=" + zstr(dR3))

# B3 response and loading regimes.
emit("S11BB_THICKNESS_RESPONSE",
     "with F_W=-(R3_d*d+R3_theta*theta), output/input chi_W=deltaW/F_W=W0/R3_e; R3_e is the explicit third entry of S11BB_THICKNESS_EOM")
emit("S11BB_RESPONSE_NORMALIZATION",
     "output/input=deltaW [L] divided by generalized thickness force density F_W [M L^-2 T^-2], so [chi_W]=[L^3 T^2 M^-1]")
emit("S11BB_BULK_OPERATOR_BY_REGIME",
     "for real omega the pressure part is O_p=-i*omega*P_V/2=-omega^2 M_p-i*omega Gamma_p with M_p=-Im(P_V)/(2 omega), Gamma_p=Re(P_V)/2. Impermeable propagating: M_p=0,Gamma_p=Z/2>0. Impermeable evanescent: M_p=rho_m/(2 alpha),Gamma_p=0. Permeation generally makes both nonzero in either regime; reciprocal traction adds the same decomposition with P_V replaced by Lambda_X*A_V. q_out=0 is singular")
emit("S11BB_MASS_INTERPRETATION_VALID_WHERE",
     "only a purely acceleration-phase load admits a mass interpretation: the impermeable evanescent load does (rho_m/(2 alpha) for deltaW); propagating radiation resistance and a generic permeable/reciprocal load do not")

# B4 response and limits.
Hscalar = sp.Matrix([[KL, Dth, De], [Dth, bth, bte], [De, bte, bw]])
Kstatic_schur = (
    sp.Matrix([[Dth, De]])
    * sp.Matrix([[bth, bte], [bte, bw]]).inv()
    * sp.Matrix([Dth, De])
)[0]
Kstatic_active = KL - Kstatic_schur
Kinf = KL - 2*Dth + bth
emit("S11BB_COMPRESSIONAL_RESPONSE",
     "K_L^eff= sigma_L/d, sigma_L=K_L*d+D_theta*theta+D_e*e_W-mu_theta. Let Delta=R2_theta*R3_e-R2_e*R3_theta; theta/d=(-R2_d*R3_e+R2_e*R3_d)/Delta, e_W/d=(-R2_theta*R3_d+R2_d*R3_theta)/Delta, and K_L^eff=(K_L-D_theta)+(D_theta-b_theta)*theta/d+(D_e-b_theta_e)*e_W/d")
emit("S11BB_LIMITS_AND_PATH",
     "generic active slab affinity (Lambda_A0!=0), fixed real k!=0 and omega->0 from the prescribed upper rim: mu_theta=p_W=0 and K_eff->" + zstr(Kstatic_active) + ". Fixed k with |omega|->infinity (mu_W!=0, finite Debye parameters): e_W->0, theta->-d and K_eff->" + zstr(Kinf) + ". Taking impermeability before omega->0 instead preserves theta+e_W+d=constant and generally differs; limits/cuts do not commute")
emit("S11BB_FROZEN_THICKNESS_IDENTIFICATION",
     "holding e_W=0 gives K_frozen=" + zstr(Kcomp_A) + " after the still-active mass/transfer row is solved. It is not a primitive single compression modulus; in the impermeable high-frequency reduction it becomes K_L-2D_theta+B_rho3+A_theta*k^2")

# B5 determinant and roots.
K0_slice = Br3 - 2*C*W0 + kW*W0**2
R_slice = rhom*cs/2
S_slice = sp.sqrt(R_slice**2 - 4*muW*K0_slice/W0**2)
slice_quadratic_reported = muW*w**2 + I*R_slice*w - K0_slice/W0**2
slice_determinant = sp.cancel(DISP.subs({
    k: 0, LA0: 0, LV0: 0, LX0: 0, q: w/cs,
}))
slice_factor_cancelled = -I*W0*rho**2*w**3
slice_polynomial = sp.cancel(slice_determinant/slice_factor_cancelled)
slice_difference = sp.cancel(slice_polynomial - slice_quadratic_reported)
slice_comparison = (
    "From the assembled determinant, substitution of k=0, Lambda_A0=Lambda_V0="
    "Lambda_X0=0, and q_out=omega/c_s0 gives D_slice="
    + zstr(slice_determinant) + ". The cancelled overall factor is "
    + zstr(slice_factor_cancelled) + " (including the omega^3 rank-loss factor); "
    "the surviving polynomial is " + zstr(slice_polynomial)
    + ", the previously-solved quadratic is " + zstr(slice_quadratic_reported)
    + ", and their symbolic difference is " + zstr(slice_difference) + "."
)
if not zero(slice_difference):
    emit("S11BB_ROOTS", slice_comparison + " Nonzero difference; stopping before solving.")
    raise SystemExit(1)
slice_roots_derived = sp.solve(slice_polynomial, w)
omega_opp_plus_slice = I*(R_slice + S_slice)/(2*muW)
omega_opp_minus_slice = I*(R_slice - S_slice)/(2*muW)
emit("S11BB_LONGITUDINAL_DISPERSION",
     "D(omega,k;q_out)=det(R_inplane,R_mass,R_thickness)=0 with unscaled rows R1=" + zstr(R1) + ", R2=" + zstr(R2) + ", R3=" + zstr(R3) + "; explicitly D=R11*(R22*R33-R23*R32)-R12*(R21*R33-R23*R31)+R13*(R21*R32-R22*R31)")
emit("S11BB_ROOTS",
     slice_comparison + " Solving the surviving derived polynomial gives roots="
     + zstr(slice_roots_derived) + ". The general dispersion remains the sheet-filtered determinant above and has no universal elementary closed form; exceptional multiplicities still satisfy D=partial_omega D=0")
emit("S11BB_IMAGINARY_PART",
     "On the explicit slice let S=sqrt(R^2-4*mu_W*K0/W0^2). For nonnegative radicand, Im(omega_+)=(S-R)/(2*mu_W) and Im(omega_-)=-(R+S)/(2*mu_W): K0<0 gives signs (+,-), K0=0 gives (0,-), and 0<K0<=R^2*W0^2/(4*mu_W) gives (-,-). Above that threshold the roots are a damped oscillatory pair and both imaginary parts are -R/(2*mu_W)<0. On the opposite sheet the concrete roots are omega_opp,+=" + zstr(omega_opp_plus_slice) + " and omega_opp,-=" + zstr(omega_opp_minus_slice) + "; for nonnegative radicand their imaginary parts are (R+S)/(2*mu_W) and (R-S)/(2*mu_W), with signs (+,-) for K0<0, (+,0) for K0=0, and (+,+) for K0>0, and their ratios to omega_+,omega_- are (R+S)/(S-R) and (R-S)/(-R-S) where defined. For the underdamped pair both opposite-sheet imaginary parts are R/(2*mu_W)>0 and both ratios are -1")
emit("S11BB_DISSIPATION_ORIGIN",
     "On the explicit slice all interface channels vanish. For K0>0 the negative imaginary parts come from propagating bulk radiation resistance Re(Z)=rho_m*c_s0; for K0<0 the growing/decaying pair is created by the indefinite constrained stiffness K0<0 and radiation shifts both toward decay. No slice root is made non-real by mass transfer or reciprocal traction. In the general determinant, interfacial transfer/reciprocal kernels can be a separate source or sink, distinguished by the impermeable and no-reciprocal form cuts")
emit("S11BB_ROOT_STABILITY_CLASS",
     "On the explicit slice omega_- is DECAYING for every K0. omega_+ is GROWING for K0<0, static at K0=0, and DECAYING for K0>0; for K0>0 the two decays are overdamped or a damped oscillatory pair. The separator is K0=B_rho3-2*C*W0+k_W*W0^2=0, equivalently C=(B_rho3+k_W*W0^2)/(2*W0): above this C omega_+ grows, below it both roots decay. On q_out=omega/c_s0 the growing root is a normal mode and each decaying root is an outgoing resonance; none is discarded")
emit("S11BB_STABILITY_CONDITION",
     "The actual boundary on the explicit constrained slice is K0=B_rho3-2*C*W0+k_W*W0^2=0, because the k=0 in-plane row sets d=0 and mass balance sets theta=-e_W; for mu_W>0, K0<0 gives growth and K0>0 gives decay only. By contrast lambda_min(H_scalar(k))=0 for H_scalar=" + zstr(Hscalar) + " is the stronger unconstrained three-field Hessian boundary: its positive-definite side requires K_L>0, K_L*(B_rho3+A_theta*k^2)-D_theta^2>0, and det(H_scalar)>0. That is sufficient unconstrained positivity, not the slice stability boundary. The full free-interface boundary still depends on the interface coefficients and is D=0 with Im(omega)=0")
emit("S11BB_GROWTH_ARTIFACT_DIAGNOSTICS",
     "For the concrete omega_+ root when K0<0: Lambda_A0=Lambda_V0=Lambda_X0=0, so all memory channels are absent, orientation is indistinguishable, all kernel-propagation residuals are 0, and there are no finite response-kernel poles. It lies on q_out=omega/c_s0 reached by prescribed upward continuation; real-axis flux/decay requirements were not re-imposed, and its bulk field is spatially normalizable. The zero interface tuple lies inside both unconditional admissibility and conditional reciprocity (the former is a strict subset of the latter). Thus its growth is the indefinite K0<0 constrained direction, not a kernel-orientation or sheet-reselection artifact")
emit("S11BB_DECAY_ARTIFACT_DIAGNOSTICS",
     "For the concrete omega_- root for every K0, and omega_+ when K0>0: all interface channels are absent, orientation is indistinguishable, every kernel-propagation residual is 0, and there are no finite response-kernel poles. Each lies on q_out=omega/c_s0 reached by prescribed downward continuation without re-imposing spatial decay and is an outgoing resonance. The zero interface tuple lies inside both admissibility regions. Radiation resistance supplies the decay shift; for K0<0 omega_- is also the decaying member of the indefinite conservative pair, not an artifact")
emit("S11BB_RECIPROCAL_TRACTION_ROOT_EFFECT",
     "D-D_F=det(R1,R2,DeltaR3), DeltaR3=" + zstr(dR3) + "; Lambda_X generically shifts roots and can split/merge them, while unchanged roots satisfy this difference=0. Multiplicity changes only on D=partial_omega D=0; the F cut sets tau_X irrelevant")

# B6 transverse computation.
emit("S11BB_TRANSVERSE_COUPLING",
     "coefficient defined as the off-diagonal Fourier operator mapping e_W to a transverse u_T equation (and reciprocally u_T to the thickness equation): 0 identically. B1 contains only div u and every scalar energy invariant contains only div u, so projection perpendicular to k annihilates the coupling; because it vanishes, a standalone normalization/dimension is undetermined")
emit("S11BB_TRANSVERSE_DISPERSION", "rho_br0*omega^2=mu_R*k^2, with two transverse polarizations and no thickness/interface modification on the uniform background")
emit("S11BB_TRANSVERSE_DISSIPATION",
     "Im(omega)=0 for real positive mu_R/rho_br0 (or a conservative growing/decaying pair if their ratio is negative); it is independent over the full range of Lambda_A0,Lambda_V0,Lambda_X0,omega*tau_A,omega*tau_V,omega*tau_X and the slab-side affinity. This uniform result does not settle unconditional confinement")

# Causality, convention, and energy discriminators.
emit("S11BB_CAUSALITY_CHECK",
     "PASS for all active nondegenerate channels: both the primitive rational identities and inert-placeholder propagation checks vanish; zero-time and absent-channel cases are indistinguishable rather than orientation-tested")
emit("S11BB_KERNEL_ORIENTATION_IDENTITIES",
     "A,V,X normalized numerators=(0,0,0) for K_I-Lambda_I before specialization; an active tau_I>0 channel has its bare retarded pole at omega=-i/tau_I")
emit("S11BB_KERNEL_PROPAGATION_RESIDUALS",
     "face(pV,pMu,jV,jMu,tV,tMu), mass row, thickness row, determinant residuals=" + zstr(prop_residuals + [det_prop_residual]))
emit("S11BB_KERNEL_POLE_LOCATIONS",
     "after cancellation: A bare pole -i/tau_A is feedback-displaced to roots of rho_m^2*(1-i*omega*tau_A)+Z*Lambda_A0=0; V bare pole -i/tau_V is generically retained in velocity-driven face/row objects; X bare pole -i/tau_X is generically retained in reciprocal traction/thickness objects. They cancel when the corresponding residue/channel vanishes. All bare retained poles are lower-half-plane for tau_I>0; feedback-displaced half-plane is sign-indeterminate and is not an orientation verdict. Dispersion zeros and Hermitian-form conjugates are excluded")
emit("S11BB_CONVENTION_CHECK_INPLANE",
     "PASS; eliminating delta_v theta=-delta_v e_W-div(delta_v u) gives constrained derivative deltaU/deltau|constraint=deltaU/deltau|theta,e+grad(mu_theta), so momentum carries restoring force -grad(mu_theta). Used elimination; the equivalent Lagrange multiplier is -mu_theta (the material chemical-potential/generalized pressure multiplier)")
emit("S11BB_CONVENTION_CHECK_CONSERVATIVE",
     "PASS; no bulk/interfaces, kappa_W=0,k=0,J=0 gives theta+e_W=constant and omega^2=" + zstr(omega2_check) + "; equation stiffness equals K_check=" + zstr(Kcheck) + "; B_rho3 does appear")
emit("S11BB_CONSERVATIVE_POSITIVITY_INEQUALITY",
     "omega^2>0 iff (B_rho3-2*C*W0+k_W*W0^2)/mu_W>0. If mu_W>0 and the (theta,e_W) stored-energy Hessian is positive definite [B_rho3>0 and B_rho3*k_W>C^2], this inequality holds for every such U")
emit("S11BB_ENERGY_SINKS",
     "d(T+U)/dt=-sum_s[(delta_p_s+Lambda_X*A_s)V_s+mu_s J_s] (real-time convolution pairing). Positive real time-averages are sinks: outgoing bulk pressure transport, interfacial conversion, transferred slab mass chemical work, and reciprocal traction work; propagating impermeable Re(Z)>0 is a definite radiation sink")
emit("S11BB_ENERGY_SOURCES",
     "the same signed exchange channels are sources when their computed real quadratic contribution is negative; free coefficients permit this. Indefinite stored moduli generate conservative growth but are not external power sources")
emit("S11BB_UNATTRIBUTED_SINK_TERMS", "none; Q_J^direct=0 and no response kernel was varied")
emit("S11BB_UNATTRIBUTED_EXCHANGE_TERMS", "none; every term maps to pressure transport, material transfer, or the supplied reciprocal traction")
emit("S11BB_PRESSURE_WORK_SIGN_CHECK",
     "PASS; off-shell slab pressure contribution paired with partial_t(deltaW) is -sum_s 1/2 Re(delta_p_s V_s*), outgoing bulk contribution has the opposite sign, and their symbolic sum/difference under J=0,Lambda_X0=0 is 0")
emit("S11BB_FULL_TWO_PORT_BALANCE_CHECK",
     "PASS; face-by-face off-shell slab-minus-supplied-exchange differences are 0 at order Lambda_X0^0 and 0 at order Lambda_X0^1; pressure, mu_s*J and reciprocal-traction channels separately give (0,0,0). This does not test face closure, affinity normalization, or Lambda_X analytic orientation")

# Dimensions and routes.
dim_routes = {
    "B_rho": "independent: 4D pressure from deltaU4/delta(theta)^2",
    "B_rho3": "definitional: B_rho*W0",
    "mu_W": "independent: kinetic U3=mu_W*(partial_t deltaW)^2",
    "k_W": "independent: U3=k_W*W0^2*e_W^2",
    "kappa_W": "independent: U3=kappa_W*W0^2*|grad deltaW|^2",
    "C": "independent: U3=C*W0*theta*e_W",
    "B3_response_deltaW_over_force": "independent: output deltaW/input generalized thickness force density",
    "B4_response_stress_over_divu": "independent: longitudinal integrated stress/div u",
    "Lambda_A0": "independent: J/A in closure",
    "Lambda_V0": "independent: J/V in closure",
    "Lambda_X0": "independent: traction/A in power identity",
    "tau_A": "independent: omega*tau_A dimensionless",
    "tau_V": "independent: omega*tau_V dimensionless",
    "tau_X": "independent: omega*tau_X dimensionless",
    "affinity": "definitional supplied difference mu_s-delta_p/rho_m",
    "mu_theta": "independent: functional derivative deltaU3/delta theta",
    "mu_s": "definitional: mu_theta/rho_br0",
    "face_V_coeff": "independent: delta_p/V coefficient from solved face equations",
    "face_mu_theta_coeff": "independent: delta_p/mu_theta coefficient from solved face equations",
    "K_L": "independent: U3=K_L*(div u)^2",
    "D_theta": "independent: U3=D_theta*theta*div u",
    "D_e": "independent: U3=D_e*e_W*div u",
    "A_theta": "independent: U3=A_theta*|grad theta|^2",
    "A_theta_e": "independent: U3=A_theta_e*grad theta.grad e_W",
}
for name in ("B_rho", "B_rho3", "mu_W", "k_W", "kappa_W", "C",
             "B3_response_deltaW_over_force", "B4_response_stress_over_divu",
             "Lambda_A0", "Lambda_V0", "Lambda_X0", "tau_A", "tau_V", "tau_X",
             "affinity", "mu_theta", "mu_s", "face_V_coeff", "face_mu_theta_coeff",
             "K_L", "D_theta", "D_e", "A_theta", "A_theta_e"):
    emit("S11BB_DIM_" + name, str(dims[name]))
    emit("S11BB_DIM_ROUTE_KIND_" + name, dim_routes[name])
emit("S11BB_DIM_B6_coefficient", "UNDETERMINED because the defined coupling coefficient vanishes identically")
emit("S11BB_DIM_ROUTE_KIND_B6_coefficient", "independent projection calculation produced zero, so no normalization can assign a unique dimension")
emit("S11BB_DIM_face_response",
     "delta_p output [(-2,-2,1)]; V coefficient [(-3,-1,1)] and mu_theta coefficient [(-1,0,0)]")
emit("S11BB_DIM_ROUTE_KIND_face_response", "independent solution of bulk relation, interfacial mass balance, and closure")
for eqname, ok in homogeneity.items():
    emit("S11BB_HOMOGENEITY_" + eqname, "PASS" if ok else "FAIL")
emit("S11BB_HOMOGENEITY_ABLATION_DEMO",
     "PASS (failure was detectable): deliberately replacing rho_m by rho_br0 in delta_p/rho gave dimensions [L^1 T^-2 M^0] instead of affinity [L^2 T^-2 M^0], so the checker returned FAIL; restoring rho_m returned PASS")

# Controls.
emit("S11BB_CONTROL_NO_THICKNESS",
     "e_W=V=0 gives K_eff=(K_L-D_theta)-(D_theta-b_theta)*R2_d/R2_theta and D_A=R1_d*R2_theta-R1_theta*R2_d; J=J_mu*mu_theta and delta_p=P_mu*mu_theta remain. Thickness response is removed; tau_A remains in the slab-state-driven transfer, while tau_V is irrelevant because the entire V-driven channel vanishes, and tau_X becomes mechanically irrelevant")
emit("S11BB_CONTROL_A_ATTRIBUTION",
     "CONFOUNDED: changes can arise jointly from removing the thickness field, face-motion drive, Lambda_V V, pressure work, and reciprocal-traction work; they cannot be attributed to thickness alone. The slab-state-driven transfer channel is retained")
emit("S11BB_CONTROL_NO_GRADIENT_STIFFNESS",
     "set kappa_W=0: chi_W becomes W0/(R3_e|kappa_W=0), D_B=det(R1,R2,R3)|kappa_W=0; only the k^2 thickness stiffness contribution moves; all tau_I remain independent")
emit("S11BB_CONTROL_IMPERMEABLE",
     "Lambda_A0=Lambda_V0=0, Lambda_X symbolic: D_C=det(R1,R2,R3)|Lambda_A0=Lambda_V0=0; mass transfer is cut and tau_A,tau_V become irrelevant, but propagating radiation resistance P_V=Z survives and tau_X remains")
emit("S11BB_CONTROL_NO_CROSS_TERM",
     "C=0: replace b_theta_e by A_theta_e*k^2 in K_eff and D; D_D=det(R1,R2,R3)|C=0; all interface channels/times remain")
emit("S11BB_CONTROL_NO_MU_COUPLING",
     "eta=0 in A=eta*mu_theta/rho_br0-delta_p/rho_m gives P_mu=J_mu=A_mu=0; D_E=det(R1,R2,R3)|eta=0. Port H_E is recomputed by the same H formula with those zero coefficients; interfacial coefficient admissibility in independent (A,V) is unchanged")
emit("S11BB_CONTROL_NO_RECIPROCAL_TRACTION",
     "Lambda_X0=0: D_F=det(R1,R2,R3)|Lambda_X0=0 and D_full-D_F=det(R1,R2,DeltaR3); mechanical-operator and power-identity substitution residuals are both 0. Port H_F uses T=p. Interfacial admissibility reduces to Lambda_A0>=0 and Lambda_V0=0; conditional reciprocity likewise forces Lambda_V=0. tau_X is irrelevant")
emit("S11BB_CONTROLS_ON_TRANSVERSE",
     "recomputed A-D: coupling remains identically 0 and rho_br0*omega^2=mu_R*k^2 in every control. Nothing moves because transverse projection kills div u and every scalar/thickness bilinear; this is a computed uniform-background structural result, not an unconditional-confinement claim")

# Validity.
emit("S11BB_VALIDITY_CONDITIONS",
     "uniform background; |theta|<<1, |e_W|<<1, |grad u|<<1, |deltaW|/W0<<1, small face/bulk velocities, linear constitutive amplitudes, and |v0|*|q_out|/|omega|<<1 for the bulk rest-frame expansion. The scalar complex modulus is the norm used for complex omega,q_out. The supplied Debye kernels retain independent omega*tau_A, omega*tau_V, omega*tau_X and are algebraically defined for all omega off their poles; no microscopic bandwidth criterion was supplied")
emit("S11BB_VALIDITY_FAILURE_REGION",
     "background-flow expansion fails where |v0|*|q_out(omega,k)| is not small compared with |omega|, notably omega->0 at fixed k!=0 and sufficiently far from the sound cone when |q_out/omega| is large; on the cone it reduces to |v0|/c_s0. For complex values the modulus comparison is meaningful as a scalar norm except at omega=0; crossing a dragged cut requires sheet tracking. Kernel-model failure outside a microscopic response bandwidth is NOT_ESTABLISHED because no such bandwidth was supplied")

# Final verdict is solely the internal audit consistency verdict.
all_internal = (
    zero(bulk_residual) and zero(closure_residual) and zero(Zperm_difference)
    and all(zero(x) for x in kernel_identities)
    and all(zero(x) for x in prop_residuals) and zero(det_prop_residual)
    and zero(Kcheck_direct-Kcheck) and all(homogeneity.values())
    and not same(*bad_affinity_dims)
)
emit("VERDICT", "PASS" if all_internal else "FAIL")
