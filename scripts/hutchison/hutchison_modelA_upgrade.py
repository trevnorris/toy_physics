
"""
Hutchison Toy Model (Model A) -- upgraded coupled dynamics
----------------------------------------------------------
State:
  z(t), v(t)      : center-of-mass vertical motion under gravity + averaged ponderomotive force
  u(t), vq(t)     : internal collective mode quadratures (slow envelope, rotating-wave approx)

Notes:
- This is a thought-experiment model, not a blueprint for apparatus.
- Drive p'(z,t) = pA cos(k z) cos(omega t)
- COM feels averaged force from gradient of <p'^2>
- Internal mode is driven near resonance; its amplitude can optionally feed back into COM damping
"""

import math
import numpy as np

def simulate(
    L=10e-3,
    rho_m=2700.0,
    rho_star=1.2,
    c_star=578.0,
    chi=1.0,
    Q=300.0,
    pA=5000.0,
    f_drive=None,
    z0=0.0, v0=0.0,
    u0=0.0, vq0=0.0,
    gamma_z=8.0,
    t_end=0.25,
    dt=2e-5,
    k_lock=True,
    k=None,
    epscrit=1e-4,
    feedback_kappa=0.0,
    eta_lift=0.0,
    burn_frac=0.5,
    sample_stride=25
):
    g = 9.80665
    if k_lock:
        k = math.pi / L
    elif k is None:
        raise ValueError("Provide k if k_lock=False")

    f0 = chi*c_star/(2*L)
    if f_drive is None:
        f_drive = f0
    w = 2*math.pi*f_drive
    w0 = 2*math.pi*f0

    gamma_int = w0/(2*Q)
    Delta = (w0*w0 - w*w)/(2*w)

    # averaged ponderomotive acceleration coefficient
    a0_base = (pA*pA * k) / (4.0 * rho_m * rho_star * c_star*c_star)

    n = int(t_end/dt)
    burn = int(burn_frac*n)

    z=float(z0); vz=float(v0); u=float(u0); vq=float(vq0)

    # diagnostics
    count=0
    z_sum=z2_sum=0.0
    lift_sum=lift2_sum=0.0
    eps_sum=0.0
    eps_peak=0.0
    z_hist=[]

    def deriv(z, vz, u, vq):
        amp = math.sqrt(u*u + vq*vq)
        eps = (amp/1.41421356237)/L  # rms strain proxy
        gz = gamma_z*(1.0 + feedback_kappa*(eps/epscrit)**2) if (feedback_kappa!=0.0 and epscrit>0) else gamma_z

        a0 = a0_base*(1.0 + eta_lift*(eps/epscrit)**2) if (eta_lift!=0.0 and epscrit>0) else a0_base
        a_rad = a0*math.sin(2.0*k*z)

        dz=vz
        dvz=-g + a_rad - gz*vz

        p0 = pA*math.cos(k*z)             # local pressure amplitude
        f_drv = p0/(rho_m*L)              # assumes (β/α)=1 and S/V~1/L
        du = -gamma_int*u + Delta*vq + f_drv/(2.0*w)
        dvq = -gamma_int*vq - Delta*u

        return dz,dvz,du,dvq,a_rad,eps,a0

    for i in range(n):
        dz1,dvz1,du1,dvq1,_,_,_ = deriv(z,vz,u,vq)
        z2=z+0.5*dt*dz1; vz2=vz+0.5*dt*dvz1; u2=u+0.5*dt*du1; vq2=vq+0.5*dt*dvq1
        dz2,dvz2,du2,dvq2,_,_,_ = deriv(z2,vz2,u2,vq2)
        z3=z+0.5*dt*dz2; vz3=vz+0.5*dt*dvz2; u3=u+0.5*dt*du2; vq3=vq+0.5*dt*dvq2
        dz3,dvz3,du3,dvq3,_,_,_ = deriv(z3,vz3,u3,vq3)
        z4=z+dt*dz3; vz4=vz+dt*dvz3; u4=u+dt*du3; vq4=vq+dt*dvq3
        dz4,dvz4,du4,dvq4,_,eps,a0 = deriv(z4,vz4,u4,vq4)

        z  += (dt/6.0)*(dz1+2*dz2+2*dz3+dz4)
        vz += (dt/6.0)*(dvz1+2*dvz2+2*dvz3+dvz4)
        u  += (dt/6.0)*(du1+2*du2+2*du3+du4)
        vq += (dt/6.0)*(dvq1+2*dvq2+2*dvq3+dvq4)

        if sample_stride and (i % sample_stride == 0):
            z_hist.append(z)

        if i >= burn:
            count += 1
            z_sum += z; z2_sum += z*z
            lift = a0*math.sin(2.0*k*z)/g
            lift_sum += lift; lift2_sum += lift*lift
            eps_sum += eps
            if eps > eps_peak: eps_peak = eps

    z_mean=z_sum/count
    z_std=math.sqrt(max(0.0,z2_sum/count - z_mean*z_mean))

    return {
        "L_mm": L*1000,
        "f0_kHz": f0/1000,
        "f_drive_kHz": f_drive/1000,
        "Q": Q,
        "pA_Pa": pA,
        "z_mean_m": z_mean,
        "z_std_m": z_std,
        "lift_mean_g": lift_sum/count,
        "lift_rms_g": math.sqrt(lift2_sum/count),
        "eps_rms_mean": eps_sum/count,
        "eps_rms_peak": eps_peak,
        "z_hist": np.array(z_hist)
    }
