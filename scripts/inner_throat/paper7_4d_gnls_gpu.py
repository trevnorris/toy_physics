from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Dict, Any, Optional
import json, math, time

# -----------------------------
# Backend selection (NumPy/CuPy)
# -----------------------------
def select_backend(prefer_gpu: bool = True):
    if prefer_gpu:
        try:
            import cupy as cp
            print("[System] CuPy detected. GPU Mode Engaged.")
            return cp, cp.fft, cp.asnumpy, True
        except Exception as e:
            print(f"[System] GPU unavailable ({type(e).__name__}: {e}). Falling back to NumPy.")
    import numpy as np
    print("[System] CPU Mode (NumPy).")
    return np, np.fft, (lambda a: a), False

# -----------------------------
# Frozen functional forms (mirror definitions.wl)
# -----------------------------
def smooth_step(xp, u):      # (1+tanh)/2
    return 0.5 * (1.0 + xp.tanh(u))

def smooth_step_prime(xp, u):  # 0.5*sech^2(u)
    c = xp.cosh(u)
    return 0.5 / (c * c)

def smooth_abs(xp, u, eps):  # sqrt(u^2+eps^2)
    return xp.sqrt(u*u + eps*eps)

def smooth_abs_prime(xp, u, eps):  # u / sqrt(u^2+eps^2)
    return u / smooth_abs(xp, u, eps)

def r3(xp, x, y, z):
    return xp.sqrt(x*x + y*y + z*z)

def gate_throat(xp, R, w, a, L, dr, dw, epsW):
    GateR = 1.0 - smooth_step(xp, (R - a) / dr)
    GateW = 1.0 - smooth_step(xp, (smooth_abs(xp, w, epsW) - 0.5 * L) / dw)
    return GateR * GateW

@dataclass
class Params:
    # core
    hbar: float
    m: float
    K: float
    # Vconf(F1)
    dr: float
    dw: float
    epsW: float
    V0: float
    p: float
    OmOut: float
    OmIn: float
    # geometry energy
    Pvac: float
    Sig: float
    # misc/freeze
    rPort: float
    wProjNSig: float
    Wcut: float
    Rmeasure: float
    epsAng: float
    wDrive: float
    rDrive: float
    freezeSHA256: str = ""

    dtype_real: str = "float32"
    dtype_cplx: str = "complex64"

    @staticmethod
    def from_json(path: str) -> "Params":
        import json
        with open(path, "r") as f:
            d = json.load(f)

        # Support wrapped schema: {freezeSHA256, operational, params, ...}
        if isinstance(d, dict) and "params" in d and isinstance(d["params"], dict):
            params = dict(d["params"])
            # carry over freeze signature if present at top-level
            if "freezeSHA256" in d and "freezeSHA256" not in params:
                params["freezeSHA256"] = d["freezeSHA256"]
            d = params

        # Drop any unknown keys safely
        allowed = set(Params.__dataclass_fields__.keys())
        d = {k: v for k, v in d.items() if k in allowed}

        return Params(**d)

class Grid4D:
    def __init__(self, xp, fft, Nx, Ny, Nz, Nw, Lx, Ly, Lz, Lw, dtype_real="float32"):
        self.xp, self.fft = xp, fft
        self.Nx, self.Ny, self.Nz, self.Nw = Nx, Ny, Nz, Nw
        self.Lx, self.Ly, self.Lz, self.Lw = Lx, Ly, Lz, Lw
        self.dx, self.dy, self.dz, self.dw = Lx/Nx, Ly/Ny, Lz/Nz, Lw/Nw

        self.x = xp.linspace(-Lx/2, Lx/2 - self.dx, Nx, dtype=dtype_real)
        self.y = xp.linspace(-Ly/2, Ly/2 - self.dy, Ny, dtype=dtype_real)
        self.z = xp.linspace(-Lz/2, Lz/2 - self.dz, Nz, dtype=dtype_real)
        self.w = xp.linspace(-Lw/2, Lw/2 - self.dw, Nw, dtype=dtype_real)

        kx = 2.0*xp.pi*fft.fftfreq(Nx, d=self.dx).astype(dtype_real)
        ky = 2.0*xp.pi*fft.fftfreq(Ny, d=self.dy).astype(dtype_real)
        kz = 2.0*xp.pi*fft.fftfreq(Nz, d=self.dz).astype(dtype_real)
        kw = 2.0*xp.pi*fft.fftfreq(Nw, d=self.dw).astype(dtype_real)

        self.k2 = (kx[:,None,None,None]**2 +
                   ky[None,:,None,None]**2 +
                   kz[None,None,:,None]**2 +
                   kw[None,None,None,:]**2).astype(dtype_real)

    @property
    def dV(self) -> float:
        return float(self.dx*self.dy*self.dz*self.dw)

    def bc(self):
        # broadcasted coordinate views
        x4 = self.x[:,None,None,None]
        y4 = self.y[None,:,None,None]
        z4 = self.z[None,None,:,None]
        w4 = self.w[None,None,None,:]
        return x4,y4,z4,w4

def build_vconf_F1(grid: Grid4D, a: float, L: float, P: Params, out_dtype: Optional[str] = None):
    xp = grid.xp
    if out_dtype is None:
        out_dtype = P.dtype_real
    x = xp.asarray(grid.x, dtype=out_dtype)
    y = xp.asarray(grid.y, dtype=out_dtype)
    z = xp.asarray(grid.z, dtype=out_dtype)
    w = xp.asarray(grid.w, dtype=out_dtype)
    x4 = x[:,None,None,None]
    y4 = y[None,:,None,None]
    z4 = z[None,None,:,None]
    w4 = w[None,None,None,:]
    R = r3(xp, x4,y4,z4)
    G = gate_throat(xp, R, w4, a, L, P.dr, P.dw, P.epsW)

    p_int = int(P.p)
    OmegaSq = (P.OmOut**2) - (P.OmOut**2 - P.OmIn**2)*G
    Vbrane  = 0.5*P.m*OmegaSq*(w4**2)
    Vwall   = P.V0*(smooth_step(xp, (R - a)/P.dr)**p_int)
    Vcap    = P.V0*(smooth_step(xp, (smooth_abs(xp, w4, P.epsW) - 0.5*L)/P.dw)**p_int)
    return (Vbrane + Vwall + Vcap).astype(out_dtype)

def build_vconf_F1_and_derivs(grid: Grid4D, a: float, L: float, P: Params,
                              deriv_dtype: str = "float64"):
    xp = grid.xp
    x = xp.asarray(grid.x, dtype=deriv_dtype)
    y = xp.asarray(grid.y, dtype=deriv_dtype)
    z = xp.asarray(grid.z, dtype=deriv_dtype)
    w = xp.asarray(grid.w, dtype=deriv_dtype)
    x4 = x[:,None,None,None]
    y4 = y[None,:,None,None]
    z4 = z[None,None,:,None]
    w4 = w[None,None,None,:]

    R = r3(xp, x4, y4, z4)
    uR = (R - a) / P.dr
    uW = (smooth_abs(xp, w4, P.epsW) - 0.5 * L) / P.dw

    S_R = smooth_step(xp, uR)
    S_W = smooth_step(xp, uW)
    S_R_p = smooth_step_prime(xp, uR)
    S_W_p = smooth_step_prime(xp, uW)

    GateR = 1.0 - S_R
    GateW = 1.0 - S_W
    G = GateR * GateW

    OmegaSq = (P.OmOut**2) - (P.OmOut**2 - P.OmIn**2) * G
    Vbrane = 0.5 * P.m * OmegaSq * (w4**2)
    p_int = int(P.p)
    Vwall = P.V0 * (S_R**p_int)
    Vcap = P.V0 * (S_W**p_int)
    Vconf = (Vbrane + Vwall + Vcap).astype(P.dtype_real)

    dGateR_da = S_R_p / P.dr
    dGateW_dL = S_W_p * (0.5 / P.dw)
    dG_da = dGateR_da * GateW
    dG_dL = GateR * dGateW_dL

    p1 = p_int - 1
    dVwall_da = -P.V0 * p_int * (S_R**p1) * S_R_p * (1.0 / P.dr)
    dVcap_dL = -P.V0 * p_int * (S_W**p1) * S_W_p * (0.5 / P.dw)

    delta_Om = (P.OmOut**2 - P.OmIn**2)
    dOmegaSq_da = -delta_Om * dG_da
    dOmegaSq_dL = -delta_Om * dG_dL

    dVbrane_da = 0.5 * P.m * (w4**2) * dOmegaSq_da
    dVbrane_dL = 0.5 * P.m * (w4**2) * dOmegaSq_dL

    dVda = (dVbrane_da + dVwall_da).astype(deriv_dtype)
    dVdL = (dVbrane_dL + dVcap_dL).astype(deriv_dtype)
    return Vconf, dVda, dVdL

def _build_vconf_F1_dVdL(grid: Grid4D, a: float, L: float, P: Params,
                         deriv_dtype: str = "float64", L_factor: float = 0.5):
    xp = grid.xp
    x = xp.asarray(grid.x, dtype=deriv_dtype)
    y = xp.asarray(grid.y, dtype=deriv_dtype)
    z = xp.asarray(grid.z, dtype=deriv_dtype)
    w = xp.asarray(grid.w, dtype=deriv_dtype)
    x4 = x[:,None,None,None]
    y4 = y[None,:,None,None]
    z4 = z[None,None,:,None]
    w4 = w[None,None,None,:]

    R = r3(xp, x4, y4, z4)
    uR = (R - a) / P.dr
    uW = (smooth_abs(xp, w4, P.epsW) - 0.5 * L) / P.dw

    S_R = smooth_step(xp, uR)
    S_W = smooth_step(xp, uW)
    S_W_p = smooth_step_prime(xp, uW)

    GateR = 1.0 - S_R
    GateW = 1.0 - S_W
    G = GateR * GateW

    dGateW_dL = S_W_p * (L_factor / P.dw)
    dG_dL = GateR * dGateW_dL

    delta_Om = (P.OmOut**2 - P.OmIn**2)
    dOmegaSq_dL = -delta_Om * dG_dL
    dVbrane_dL = 0.5 * P.m * (w4**2) * dOmegaSq_dL
    p_int = int(P.p)
    p1 = p_int - 1
    dVcap_dL = -P.V0 * p_int * (S_W**p1) * S_W_p * (L_factor / P.dw)
    return (dVbrane_dL + dVcap_dL).astype(deriv_dtype)

class SplitStepGNLS4D:
    def __init__(self, grid: Grid4D, P: Params, dt: float, imaginary: bool):
        self.grid, self.P, self.dt, self.imag = grid, P, float(dt), bool(imaginary)
        xp = grid.xp
        coef = (P.hbar*self.dt)/(2.0*P.m)
        if self.imag:
            self.Kprop = xp.exp(-coef*grid.k2).astype(P.dtype_cplx)
        else:
            self.Kprop = xp.exp(-1j*coef*grid.k2).astype(P.dtype_cplx)

    def _Veff(self, Vconf, psi):
        xp = self.grid.xp
        rho = (psi.real*psi.real + psi.imag*psi.imag).astype(self.P.dtype_real)
        nonlin = (5.0*self.P.K/4.0) * (rho**4)  # |psi|^8 = (|psi|^2)^4
        return (Vconf + nonlin).astype(self.P.dtype_real)

    def step(self, psi, Vconf):
        xp = self.grid.xp
        dt, hbar = self.dt, self.P.hbar

        Veff = self._Veff(Vconf, psi)
        if self.imag:
            Pprop = xp.exp(-(0.5*dt/hbar)*Veff).astype(self.P.dtype_cplx)
        else:
            Pprop = xp.exp(-1j*(0.5*dt/hbar)*Veff).astype(self.P.dtype_cplx)
        psi = psi * Pprop

        psi_k = self.grid.fft.fftn(psi, axes=(0,1,2,3))
        psi_k = psi_k * self.Kprop
        psi = self.grid.fft.ifftn(psi_k, axes=(0,1,2,3))

        Veff = self._Veff(Vconf, psi)
        if self.imag:
            Pprop = xp.exp(-(0.5*dt/hbar)*Veff).astype(self.P.dtype_cplx)
        else:
            Pprop = xp.exp(-1j*(0.5*dt/hbar)*Veff).astype(self.P.dtype_cplx)
        psi = psi * Pprop
        return psi

def norm_L2(grid: Grid4D, psi) -> float:
    xp = grid.xp
    rho = psi.real*psi.real + psi.imag*psi.imag
    s = xp.sum(rho)
    if hasattr(s, "get"): s = s.get()
    return float(s)*grid.dV

def normalize_to(grid: Grid4D, psi, Ntarget: float):
    N = norm_L2(grid, psi)
    psi = psi * math.sqrt(Ntarget/max(1e-30, N))
    return psi

def init_gaussian(grid: Grid4D, P: Params, sigma_xyz: float, sigma_w: float, seed: int = 0):
    xp = grid.xp
    x4,y4,z4,w4 = grid.bc()
    amp = xp.exp(-(x4*x4+y4*y4+z4*z4)/(2*sigma_xyz**2)) * xp.exp(-(w4*w4)/(2*sigma_w**2))
    # small random phase noise
    if xp.__name__ == "cupy":
        rs = xp.random.RandomState(seed)
    else:
        import numpy as np
        rs = np.random.RandomState(seed)
    phase = 0.2*rs.standard_normal(size=amp.shape)
    return (amp * xp.exp(1j*phase)).astype(P.dtype_cplx)

def find_ground_state_imagtime(grid: Grid4D, P: Params, Vconf, Ntarget: float,
                               dtau: float, max_steps: int = 2000,
                               tol: float = 1e-7, report_every: int = 50,
                               psi0=None):
    xp = grid.xp
    if psi0 is None:
        psi = init_gaussian(grid, P, sigma_xyz=0.5*P.rPort, sigma_w=0.5, seed=0)
    else:
        psi = psi0.astype(P.dtype_cplx, copy=True)
    psi = normalize_to(grid, psi, Ntarget)

    stepper = SplitStepGNLS4D(grid, P, dt=dtau, imaginary=True)
    prev = psi
    t0 = time.time()

    for n in range(1, max_steps+1):
        psi = stepper.step(psi, Vconf)
        psi = normalize_to(grid, psi, Ntarget)

        if n % report_every == 0:
            diff = psi - prev
            rel = norm_L2(grid, diff) / max(1e-30, norm_L2(grid, psi))
            print(f"[ImagTime] step={n:5d} rel_change={rel:.3e} elapsed={time.time()-t0:.1f}s")
            if rel < tol:
                break
            prev = psi
    return psi

def dEgeom_da(a: float, L: float, P: Params) -> float:
    # V=(4π/3)a^3L, A=(4πa^2)L + 2*(4π/3)a^3
    return (P.Pvac*(4.0*math.pi)*a*a*L +
            P.Sig*((8.0*math.pi)*a*L + (8.0*math.pi)*a*a))

def dEgeom_dL(a: float, L: float, P: Params) -> float:
    return (P.Pvac*(4.0*math.pi/3.0)*a**3 +
            P.Sig*(4.0*math.pi)*a*a)

def force_residuals(grid: Grid4D, P: Params, psi, a: float, L: float,
                    deriv_method: str = "analytic",
                    da: Optional[float] = None, dL: Optional[float] = None,
                    deriv_dtype: str = "float64") -> Tuple[float,float,Dict[str,float]]:
    xp = grid.xp
    rho = (psi.real*psi.real + psi.imag*psi.imag).astype(P.dtype_real)

    deriv_method = str(deriv_method).lower()
    info: Dict[str, float] = {}

    if deriv_method == "analytic":
        _, dVda, dVdL = build_vconf_F1_and_derivs(grid, a, L, P, deriv_dtype=deriv_dtype)
    elif deriv_method == "fd64":
        if dL is None:
            dL = max(grid.dw, 0.5 * P.dw)
        if da is None:
            da = max(grid.dx, 0.5 * P.dr)
        Vp_a = build_vconf_F1(grid, a + da, L, P, out_dtype=deriv_dtype)
        Vm_a = build_vconf_F1(grid, a - da, L, P, out_dtype=deriv_dtype)
        Vp_L = build_vconf_F1(grid, a, L + dL, P, out_dtype=deriv_dtype)
        Vm_L = build_vconf_F1(grid, a, L - dL, P, out_dtype=deriv_dtype)
        dVda = (Vp_a - Vm_a) / (2.0 * da)
        dVdL = (Vp_L - Vm_L) / (2.0 * dL)
    elif deriv_method == "fd32":
        if dL is None:
            dL = 0.01
        if da is None:
            da = 0.01
        Vp_a = build_vconf_F1(grid, a + da, L, P, out_dtype=P.dtype_real)
        Vm_a = build_vconf_F1(grid, a - da, L, P, out_dtype=P.dtype_real)
        Vp_L = build_vconf_F1(grid, a, L + dL, P, out_dtype=P.dtype_real)
        Vm_L = build_vconf_F1(grid, a, L - dL, P, out_dtype=P.dtype_real)
        dVda = (Vp_a - Vm_a) / (2.0 * da)
        dVdL = (Vp_L - Vm_L) / (2.0 * dL)
    else:
        raise ValueError(f"Unknown deriv_method='{deriv_method}'")

    Ia = xp.sum(rho*dVda); IL = xp.sum(rho*dVdL)
    if hasattr(Ia, "get"): Ia, IL = Ia.get(), IL.get()
    Ia, IL = float(Ia)*grid.dV, float(IL)*grid.dV

    Ra = (-Ia) - dEgeom_da(a,L,P)
    RL = (-IL) - dEgeom_dL(a,L,P)

    w4 = grid.w[None, None, None, :]
    acc_dtype = xp.result_type(rho, dVdL)
    mask_plus = (w4 > 0).astype(acc_dtype)
    mask_minus = (w4 < 0).astype(acc_dtype)
    IL_plus = xp.sum(rho * dVdL * mask_plus, dtype=acc_dtype)
    IL_minus = xp.sum(rho * dVdL * mask_minus, dtype=acc_dtype)
    if hasattr(IL_plus, "get"):
        IL_plus, IL_minus = IL_plus.get(), IL_minus.get()
    IL_plus = float(IL_plus) * grid.dV
    IL_minus = float(IL_minus) * grid.dV
    IL_split_rel = abs(IL_plus - IL_minus) / max(1e-30, abs(IL))

    info.update({
        "Ia": Ia,
        "IL": IL,
        "dEda": dEgeom_da(a, L, P),
        "dEdL": dEgeom_dL(a, L, P),
        "IL_plus": IL_plus,
        "IL_minus": IL_minus,
        "IL_split_rel": IL_split_rel,
    })

    if deriv_method == "analytic":
        dVdL_half = _build_vconf_F1_dVdL(grid, a, L, P, deriv_dtype=deriv_dtype, L_factor=1.0)
        IL_half = xp.sum(rho * dVdL_half)
        if hasattr(IL_half, "get"):
            IL_half = IL_half.get()
        IL_half = float(IL_half) * grid.dV
        info["IL_Lhalf"] = IL_half
        den = max(1e-30, abs(2.0 * IL))
        info["IL_Lhalf_err"] = abs(abs(IL_half) / den - 1.0)
        info["IL_Lhalf_ratio_signed"] = IL_half / (2.0 * IL) if abs(2.0 * IL) > 1e-30 else float("nan")
        info["IL_Lhalf_ratio_abs"] = abs(IL_half) / den

    return Ra, RL, info
