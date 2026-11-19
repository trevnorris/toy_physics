"""Two contacts with identical held scattering and different translation loads.

Each co-moving reference has kinetic energy m|sdot-alpha (Xdot cross D*)|²/2.
The signed radial reference D* is known on each sphere. Alpha is dimensionless.
This extra inertial coupling is passive. It is absent from the held experiment.
"""
import numpy as np

from mouth_pair_translating_cuts import AMPLITUDE, FLUX, MODULUS, SPEED, Pair, hankel, polynomial
from mouth_held_contact_prediction import held_measurement


def components(pair,k):
    ell = np.arange(1,pair.lmax+1)
    degrees = ell*(ell+1)
    norm = 4*np.pi*degrees/(2*ell+1)
    proj = np.stack([polynomial(l,pair.u)[1] for l in ell])
    proj *= ((2*ell+1)/(2*degrees))[:,None]*pair.weights
    items=[]
    for center in pair.centers:
        field = pair.field(center,1,k)
        step=.0001
        der = ((pair.background(center,-2*step)-pair.background(center,2*step))
               +8*(pair.background(center,step)-pair.background(center,-step)))/(12*step)
        items.append((proj@field[1],proj@field[2],proj@der[1],pair.pr@der[0]/degrees))
    return norm,items


def full(pair,k,ratio,loss,alpha,motion):
    # Solve the exterior and every local reference simultaneously, without
    # eliminating s. Independent of the held and translating inverse fits.
    norm,items = components(pair,k)
    n = pair.lmax*len(items)
    matrix = np.zeros((2*n,2*n),complex)
    rhs = np.zeros(2*n,complex)
    g1=h=16*MODULUS
    recoil=[]
    for i,((theta,y,t,v),x,q) in enumerate(zip(items,motion,pair.charges)):
        rows = np.arange(i*pair.lmax,(i+1)*pair.lmax)
        b = np.zeros(pair.lmax)
        b[0] = alpha*q/(4*np.pi)
        recoil.append(b)
        matrix[rows,:n] = MODULUS*theta+g1*y
        matrix[rows,rows+n] = -g1
        rhs[rows] = -x*(MODULUS*t+g1*v)
        matrix[rows+n,:n] = -g1*y
        matrix[rows+n,rows+n] = g1+h-h*(ratio**2+1j*ratio*loss)
        rhs[rows+n] = x*(g1*v-h*ratio**2*b)
    out,s = np.split(np.linalg.solve(matrix,rhs),2)
    heat = drive = 0.
    for (theta,y,t,v),x,b,slip in zip(items,motion,recoil,s.reshape(-1,pair.lmax)):
        lab = y@out
        torque = g1*(lab+x*v-slip)
        load = MODULUS*t*lab+v*torque+h*ratio**2*b*(slip-b*x)
        drive += norm@(load*(1j*SPEED*k*x)).real/2
        heat += h*loss*ratio*k*SPEED*(norm@abs(slip)**2)/2
    power = pair.outgoing(out,k)
    assert abs(power+heat-drive)<1e-7*max(drive,1e-22)
    return out,power,heat,drive


def predict(pair,k,g,extra,motion):
    _,items = components(pair,k)
    A,rhs=[],[]
    for (theta,y,t,v),x,q in zip(items,motion,pair.charges):
        A.append(MODULUS*(theta+g*y))
        source = -x*MODULUS*(t+g*v)
        source[0] -= x*extra*q/(4*np.pi)
        rhs.append(source)
    out = np.linalg.solve(np.vstack(A),np.concatenate(rhs))
    return pair.outgoing(out,k)


def main():
    single = Pair(np.array([0.]),8,8.)
    pair = Pair(np.array([-2.,2.]),8,8.)
    motion = [AMPLITUDE,-AMPLITUDE]
    loss=.02
    for ratio in (.5,1.,1.15):
        k=.06*ratio
        _,g = held_measurement(k,.06*SPEED,loss,1)
        # Calibration: prescribe translation of ONE positive cut and recover
        # its outgoing l=1 amplitude from an angle measured on a fixed sphere.
        # f(r)=f(R)h1(kr)/h1(kR), with R=1 and no incident angular wave.
        for alpha in (0.,.5,1.):
            out,_,_,_ = full(single,k,ratio,loss,alpha,[AMPLITUDE])
            # An independent quadrature reads the l=1 signal on a stationary
            # sphere of radius 1.55; higher harmonics are discarded.
            angular = (single.field(0.,1.55,k)[2]@out)
            read = .75*sum(single.weights*single.sin*angular)
            amp = read*hankel(1,k)/hankel(1,1.55*k)
            t_lab = -(1+k*hankel(1,k,True)/hankel(1,k))*amp
            b0 = FLUX/(4*np.pi)
            # For the isolated radial field: (d_z D_theta)_1=-b0 and v_1=-b0.
            extra = -MODULUS*((t_lab-b0*AMPLITUDE)
                              +g*(amp-b0*AMPLITUDE))/(b0*AMPLITUDE)
            _,target,heat,work = full(pair,k,ratio,loss,alpha,motion)
            held_only = predict(pair,k,g,0.,motion)
            both = predict(pair,k,g,extra,motion)
            print(f"ratio={ratio:g} alpha={alpha:g} held_g=({g.real:.8g},{g.imag:.8g}) "
                  f"motion_torque_over_Dstar_X=({extra.real:.8g},{extra.imag:.8g}) "
                  f"actual_over_held_only={target/held_only:.9g} wave_power={target:.9g} "
                  f"both_relative_error={abs(both/target-1):.3g} wave_over_work={target/work:.9g}")
            assert abs(both/target-1)<1e-6
            if alpha==0:
                assert abs(target/held_only-1)<1e-7


if __name__ == "__main__":
    main()
