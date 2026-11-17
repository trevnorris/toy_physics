"""Axial integration of primitive energy/momentum fluxes and exterior storage."""
from math import acos, cos, fsum, pi, sin


AREA = {2: 2*pi, 3: 4*pi, 4: 2*pi**2}
D = 2.7
SPEED = 0.37
CHANNELS = [(3, 0.87, -0.74, -1.12), (4, 1.41, 0.74, 1.12)]
CONTROLS = [(1, 1, 1), (-1, 1, 0), (-1, 0, 1), (1, 0, 0)]


def integrate(f, lo, hi, count):
    """Composite Simpson quadrature for an even number of panels."""
    step = (hi-lo)/count
    inner = fsum((4 if k % 2 else 2)*f(lo+k*step) for k in range(1, count))
    return step*(f(lo)+inner+f(hi))/3


def field(q, dx, y, n, rho):
    c = q/(rho*AREA[n]*(dx*dx+y*y)**(n/2))
    return c*dx, c*y


def surface(channel, center, radius, speed, observable, panels=160):
    n, rho, q1, q2 = channel
    w = -speed if center == D else 0.0

    def integrand(t):
        nx, ny = cos(t), sin(t)
        x, y = center+radius*nx, radius*ny
        v1x, v1y = field(q1, x, y, n, rho)
        v2x, v2y = field(q2, x-D, y, n, rho)
        readings = []
        for sign, first, second in CONTROLS:
            ux, uy = first*v1x+second*v2x, first*v1y+second*v2y
            un = ux*nx+uy*ny
            kinetic = rho*(ux*ux+uy*uy)/2
            dt_phi = second*speed*v2x
            p = 0.31-kinetic-rho*dt_phi
            if observable == "force":
                value = -(rho*ux*un+p*nx)
            else:
                if observable == "omit_dt_phi":
                    p += rho*dt_phi
                value = (kinetic+p)*un-kinetic*w*nx
            readings.append(sign*value)
        return fsum(readings)*ny**(n-2)

    return AREA[n-1]*radius**(n-1)*integrate(integrand, 0, pi, panels)


def storage(channel, d, r1, r2):
    """Volume integral of rho v1.v2 on shells centered on mouth 1."""
    n, rho, q1, q2 = channel
    assert r1+r2 < d
    # Intact shells r1<r<d-r2 have zero companion flux by Gauss.
    # Where shell intersects ball 2, integrate the surviving angular segment.
    def partial_shell(t):
        r = d-r2*cos(t)
        edge = (r*r+d*d-r2*r2)/(2*r*d)
        theta0 = acos(max(-1, min(1, edge)))

        def angular(theta):
            nx, ny = cos(theta), sin(theta)
            vx, vy = field(q2, r*nx-d, r*ny, n, rho)
            return (vx*nx+vy*ny)*ny**(n-2)

        return r2*sin(t)*integrate(angular, theta0, pi, 320)

    near = q1*AREA[n-1]/AREA[n]*integrate(partial_shell, 0, pi, 160)
    # Intact shells r>d+r2 enclose q2. Integrate the resulting radial tail.
    tail = q1*q2/(rho*AREA[n]*(n-2)*(d+r2)**(n-2))
    return near+tail


def main():
    for r1, r2 in [(0.43, 0.58), (0.93, 0.82)]:
        for ch in CHANNELS:
            n = ch[0]
            force = surface(ch, 0, r1, 0, "force")
            other_force = surface(ch, D, r2, 0, "force")
            assert force > 0 and abs(force+other_force) < 1e-8
            print(f"n={n} r1={r1} r2={r2} F1={force:+.10g} F2={other_force:+.10g}")
            for s in (0, SPEED, -SPEED, SPEED/2, -SPEED/2):
                j1 = surface(ch, 0, r1, s, "energy")
                j2 = surface(ch, D, r2, s, "energy")
                work = -s*other_force
                print(f"  s={s:+.3f}  J1={j1:+.10g}  J2={j2:+.3e}  centerWork={work:+.10g}")
                assert abs(j1-force*s) < 1e-8 and abs(j2) < 1e-10
            broken = surface(ch, 0, r1, SPEED, "omit_dt_phi")
            h = 0.0003
            near = storage(ch, D-h, r1, r2)
            far = storage(ch, D+h, r1, r2)
            kdot = SPEED*(near-far)/(2*h)
            print(f"  omit_dt_phi J1={broken:+.3e}  storage Knear={near:.12g} Kfar={far:.12g}  Kdot={kdot:+.10g}")
            assert abs(broken) < 1e-10 and abs(broken-force*SPEED) > 0.0005
            assert near > far > 0 and abs(kdot-force*SPEED) < 1e-8
    print("PASS")


if __name__ == "__main__":
    main()
