"""Same resting attachment, two physically different moving mouth contacts.

An internal angular reference s on the comoving, nonrotating support carries
inertia m. Its two positive springs are G1(y-s)^2/2 and h*s^2/2. Here y
is the angular mismatch relative to the preferred texture on the moving cut.
Optional internal drag transfers power to heat. Radius is the length unit.
"""
import numpy as np

from mouth_pair_translating_cuts import AMPLITUDE, MODULUS, SPEED, Pair, polynomial


def experiment(centers, motion, k, softness, loss, resolution=10):
    static_g = 8.
    pair = Pair(np.array(centers),resolution,static_g)
    omega = k*SPEED
    contact = MODULUS*static_g
    # h=G1=2G_static. softness=omega/omega_r, loss=D*omega_r/h.
    h = g1 = 2*contact
    reference = 1/(2-softness**2-1j*loss*softness)
    dynamic = g1*(1-reference)
    ell = np.arange(1,resolution+1)
    degrees = ell*(ell+1)
    tangential = np.stack([polynomial(l,pair.u)[1] for l in ell])
    tangential *= ((2*ell+1)/(2*degrees))[:,None]*pair.weights
    angular_norm = 4*np.pi*degrees/(2*ell+1)
    mat, rhs, geometry = [], [], []
    step=2e-5
    for z, displacement in zip(centers,motion):
        # Independently use tangential ANGLES and torques in the boundary
        # equation. The static calculator uses normal-flux coefficients.
        waves = pair.field(z,1,k)
        mat.append(tangential@(MODULUS*waves[1]+dynamic*waves[2]))
        derivative = (pair.background(z,step)-pair.background(z,-step))/(2*step)
        v = (pair.pr@derivative[0])/degrees
        t = tangential@derivative[1]
        rhs.append(-displacement*(MODULUS*t+dynamic*v))
        geometry.append((v,t))
    coeff = np.linalg.solve(np.vstack(mat),np.concatenate(rhs))
    supplied = heat = 0.
    for z,x,(v,t) in zip(centers,motion,geometry):
        y_lab = tangential@(pair.field(z,1,k)[2]@coeff)
        mismatch = y_lab+x*v
        reaction = MODULUS*t*y_lab + v*dynamic*mismatch
        supplied += angular_norm@(reaction*(1j*omega*x)).real/2
        # D*omega^2=h*loss*softness*omega, s=reference*mismatch.
        heat += h*loss*softness*omega*(angular_norm@abs(reference*mismatch)**2)/2
    far = pair.outgoing(coeff,k)
    near = sum(pair.flux(coeff,z,1.18,k) for z in centers)
    scale = max(abs(supplied),far,heat,1e-26)
    return far, supplied, heat, abs(supplied-far-heat)/scale, abs(near-far)/scale


def main():
    wr = .06*SPEED
    null_ratio = np.sqrt(18/17)
    print(f"fixed_DC_g=8 internal_h_over_C=16 internal_G1_over_C=16 omega_r={wr:.12f} "
          f"isolated_lossless_null_omega_over_omega_r={null_ratio:.12f}")
    for centers,x in (([0.],[AMPLITUDE]), ([-2.,2.],[AMPLITUDE,-AMPLITUDE])):
        for ratio in (.05,.5,.95,1.,null_ratio,1.15,2.):
            k = .06*ratio
            # A direct conservative spring with the same equilibrium modulus.
            base = Pair(np.array(centers),10,8.)
            reference_work = base.outgoing(base.emit(k,x),k)
            for loss in (0.,.02):
                result = experiment(centers,x,k,ratio,loss)
                power, drive, heat, closure, flux = result
                coarse = experiment(centers,x,k,ratio,loss,6)[0]
                scaled_error = abs(coarse-power)/reference_work
                print(f"mouths={len(x)} ratio={ratio:.9f} loss={loss:.2f} "
                      f"wave_over_static_spring={power/reference_work:.10g} "
                      f"wave_fraction_of_drive={power/drive:.9g} heat_rate={heat:.9g} "
                      f"drive={drive:.9g} work_residual={closure:.2g} shell_residual={flux:.2g} "
                      f"truncation_over_static_spring={scaled_error:.3g}")
                # At the exact isolated antiresonance only differencing error
                # survives. Compare it to the nonzero direct-contact reference.
                if len(x)==1 and loss==0 and ratio==null_ratio:
                    assert max(power,abs(drive))/reference_work < 1e-12
                else:
                    assert closure<1e-5 and flux <1e-6 and drive>0
                assert scaled_error<1e-4 and heat>=0


if __name__ == "__main__":
    main()
