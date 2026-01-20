# 4D Harris Sheet Reconnection Simulation Specification

## Overview

This document specifies the initial conditions and equations for simulating magnetic reconnection in the 4D hydrodynamic framework from Paper 7 (4d.tex). The key hypothesis is that the extra dimension $w$ provides a natural channel for reconnection that doesn't require anomalous resistivity.

**Central Question:** Does the w-dimension provide a "fast reconnection" channel that standard 3D MHD lacks?

---

## 1. Field Content

| Field | Symbol | Components | Description |
|-------|--------|------------|-------------|
| Matter | $\psi(x,y,z,w,t)$ | 1 complex | Order parameter, $\rho = |\psi|^2$ |
| Gauge | $A_i$ | 4 real | Vector potential $(A_x, A_y, A_z, A_w)$ |
| Electric | $E_i$ | 4 real | Electric field (in temporal gauge $A_0=0$) |

**Derived quantities:**
- Density: $\rho = |\psi|^2$
- Phase: $\theta = \arg(\psi)$
- Magnetic field: $B_{ij} = \partial_i A_j - \partial_j A_i$
- Current: $J_i = \rho(\partial_i\theta - A_i)$

---

## 2. Evolution Equations (Dimensionless)

### 2.1 Normalizations

| Quantity | Normalization |
|----------|---------------|
| Length | $\delta$ (sheet half-width) |
| Time | $\tau_A = \delta/v_A$ where $v_A = B_0/\sqrt{\mu_0\rho_0}$ |
| Density | $\rho_0$ |
| Magnetic field | $B_0$ |
| Velocity | $v_A$ |
| Electric field | $v_A B_0$ |

### 2.2 Dimensionless Parameters

| Parameter | Definition | Typical Value |
|-----------|------------|---------------|
| $\lambda$ | $\lambda_w/\delta$ | 2.0 |
| $\Omega$ | $\Omega_w \tau_A$ | 0.5 |
| $\beta$ | $2\mu_0 K \rho_0^5/B_0^2$ | 1.0 |
| $\sigma$ | $c/v_A$ | 10-100 |

### 2.3 Matter Field (Gauged GNLS)

$$\partial_t \psi = \frac{i}{2}\left[\nabla_4^2\psi - 2i(\mathbf{A}\cdot\nabla_4)\psi - i(\nabla_4\cdot\mathbf{A})\psi + |\mathbf{A}|^2\psi\right] - i\left[\frac{\Omega^2 w^2}{2} + \frac{5\beta}{8}\rho^4\right]\psi$$

where $\nabla_4^2 = \partial_x^2 + \partial_y^2 + \partial_z^2 + \partial_w^2$.

### 2.4 Gauge Field (Localized Maxwell)

$$\partial_t A_i = -E_i$$

$$\partial_t E_i = \sigma^2\left[Z(w)\partial_j F_{ji} + Z'(w)F_{wi}\right] - J_i$$

where:
- $F_{ij} = \partial_i A_j - \partial_j A_i$
- $Z(w) = \exp(-w^2/\lambda^2)$
- $Z'(w) = -\frac{2w}{\lambda^2}Z(w)$

### 2.5 Current

$$J_i = \rho\,\partial_i\theta - \rho\,A_i$$

Computed from $\psi$ as: $J_i = \text{Im}(\psi^*\partial_i\psi) - |\psi|^2 A_i$

---

## 3. Initial Conditions (4D Harris Sheet)

### 3.1 Gauge Field

**Harris equilibrium:**
$$A_z(y) = B_0\,\delta\,\ln(\cosh(y/\delta))$$
$$A_x = A_y = A_w = 0$$

This produces: $B_x(y) = B_0\tanh(y/\delta)$

**Tearing mode perturbation:**
$$\delta A_z = \varepsilon\cos(k_x x)\exp(-y^2/\delta^2)\exp(-w^2/\lambda^2)$$

where $k_x = 2\pi/L_x$ and $\varepsilon = 0.01$.

### 3.2 Matter Field

**Density (w-localized):**
$$\rho(w) = \rho_0\exp(-w^2/\lambda_w^2)$$

**Phase (produces Harris current, thin z-domain):**
$$\theta(y,z,w) = \theta_0\,\text{sech}^2(y/\delta)\exp(-w^2/\lambda_w^2)\cdot z$$

where $\theta_0 \sim 1$ (adjustable).

**Wave function:**
$$\psi = \sqrt{\rho}\exp(i\theta)$$

### 3.3 Electric Field

$$E_x = E_y = E_z = E_w = 0$$

(Equilibrium has no electric field)

---

## 4. Computational Setup

### 4.1 Domain

| Direction | Extent | BC | Purpose |
|-----------|--------|-----|---------|
| $x$ | $[-2\pi, 2\pi]$ | Periodic | Along sheet |
| $y$ | $[-5, 5]$ | Open/absorbing | Across sheet |
| $z$ | $[-0.05, 0.05]$ | Periodic | Thin (2.5D) |
| $w$ | $[-4, 4]$ | Soft walls | Extra dimension |

### 4.2 Resolution and Memory

| Config (Nx, Ny, Nz, Nw) | Points | Memory | Use Case |
|-------------------------|--------|--------|----------|
| (64, 128, 2, 32) | 0.5M | 0.1 GB | Quick tests |
| (128, 256, 4, 64) | 8M | 1.5 GB | **Development** |
| (256, 512, 4, 128) | 67M | 13 GB | Production |
| (512, 1024, 4, 256) | 537M | 103 GB | High-res |

### 4.3 Time Stepping

**CFL condition:**
$$\Delta t < \frac{\min(\Delta x, \Delta y, \Delta z, \Delta w)}{\sigma}$$

The light speed $\sigma = c/v_A$ dominates. For $\sigma = 100$ and $\Delta x = 0.1$, need $\Delta t < 0.001$.

**Recommended schemes:**
- RK4 for explicit time stepping
- Split-step for GNLS (spectral kinetic, real-space potential)
- Consider implicit for stiff EM waves if $\sigma$ is large

---

## 5. Implementation Template

```python
import numpy as np
from scipy.fft import fftn, ifftn

# Grid
Nx, Ny, Nz, Nw = 128, 256, 4, 64
Lx, Ly, Lz, Lw = 4*np.pi, 10.0, 0.1, 8.0

x = np.linspace(-Lx/2, Lx/2, Nx, endpoint=False)
y = np.linspace(-Ly/2, Ly/2, Ny, endpoint=False)
z = np.linspace(-Lz/2, Lz/2, Nz, endpoint=False)
w = np.linspace(-Lw/2, Lw/2, Nw, endpoint=False)
X, Y, Z, W = np.meshgrid(x, y, z, w, indexing='ij')

# Parameters
delta = 1.0
lambda_w = 2.0
B0 = 1.0
rho0 = 1.0
epsilon = 0.01
theta0 = 1.0
sigma = 10.0  # c/v_A ratio
Omega = 0.5
beta = 1.0

# Initial gauge field
Az = B0 * delta * np.log(np.cosh(Y/delta))
kx = 2*np.pi / Lx
Az += epsilon * np.cos(kx*X) * np.exp(-(Y/delta)**2) * np.exp(-W**2/lambda_w**2)
Ax = Ay = Aw = np.zeros_like(Az)

# Initial matter field  
rho = rho0 * np.exp(-W**2 / lambda_w**2)
phase = theta0 * (1/np.cosh(Y/delta))**2 * np.exp(-W**2/lambda_w**2) * Z
psi = np.sqrt(rho) * np.exp(1j * phase)

# Initial electric field
Ex = Ey = Ez = Ew = np.zeros_like(Az)

# EM localization
Z_em = np.exp(-W**2 / lambda_w**2)
Z_em_prime = -2*W/lambda_w**2 * Z_em

# Confinement
V_conf = 0.5 * Omega**2 * W**2

# Time evolution (pseudocode)
dt = 0.001
for step in range(num_steps):
    # Compute density and phase from psi
    rho = np.abs(psi)**2
    
    # Compute current J_i = Im(psi* d_i psi) - rho * A_i
    # ... (use FFT for derivatives)
    
    # Update E from Maxwell
    # dE/dt = sigma^2 * [Z * curl(curl(A)) + Z' * F_wi] - J
    
    # Update A from E
    # dA/dt = -E
    
    # Update psi from GNLS
    # dpsi/dt = (i/2)[laplacian - gauge terms] - i[V + h(rho)] psi
    
    # Diagnostics every N steps
    if step % 100 == 0:
        measure_reconnection_rate()
        measure_leakage()
        save_checkpoint()
```

---

## 6. Diagnostics

### 6.1 Reconnection Rate

**Primary metric:**
$$\Gamma = E_z(x_{\text{X-point}}, y=0, w=0)$$

**Alternative (flux-based):**
$$\Gamma = \frac{d\Psi}{dt} \quad\text{where}\quad \Psi = \int_{-\infty}^{0} B_x(y)\,dy$$

**4D version (brane-projected):**
$$\Gamma_{\text{brane}} = \langle E_z \rangle_W = \int W(w) E_z \,dw$$

### 6.2 Leakage (Unique to 4D)

$$S_{\text{leak}} = \int W'(w) j_w \,dw$$

This measures flux exchange with the w-dimension.

### 6.3 Structure Diagnostics

- X-point location: $(x_*, y_*)$ where $B_x = B_y = 0$
- Current sheet width: $\delta(t)$ from half-width of $J_z(y)$
- Outflow velocity: $v_y(y=\pm\delta)$
- Energy partition: kinetic, magnetic, internal

---

## 7. Expected Outcomes

### If 4D provides fast reconnection:
- $\Gamma \sim 0.01 - 0.1$ (Petschek-like) without explicit resistivity
- Significant $S_{\text{leak}}$ at the X-point
- Scaling independent of system size (unlike Sweet-Parker)

### If it doesn't:
- $\Gamma \sim S^{-1/2}$ (resistive MHD scaling)
- $S_{\text{leak}}$ negligible
- w-dimension is a spectator

**Either outcome is scientifically valuable.** The first validates 4D hydrodynamics as a plasma tool. The second constrains where the framework applies.

---

## 8. Comparison Points

| Quantity | Sweet-Parker | Petschek | This Model |
|----------|--------------|----------|------------|
| Rate $\Gamma$ | $S^{-1/2} \sim 10^{-6}$ | $\sim 0.01-0.1$ | **?** |
| Mechanism | Resistive diffusion | Slow shocks | w-leakage? |
| Sheet width | Thins to $\delta/\sqrt{S}$ | Localized X-point | **?** |
| Requires | Resistivity $\eta$ | Open geometry | 4D structure |

---

## References

- **4d.tex** (Paper 7): Framework definition
- **toy_model_summary.md**: Connection to GR matching  
- **1pn_derivation_val.md**: Derivation of $\beta=3$, $n=5$, $\alpha^2=3/4$
- Biskamp, D. (2000). *Magnetic Reconnection in Plasmas*
- Yamada, M. et al. (2010). Rev. Mod. Phys. 82, 603
