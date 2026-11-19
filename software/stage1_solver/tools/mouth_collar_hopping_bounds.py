"""Free-link energy witness for ordinary spin-blind hopping across prepared sites.

Common, real emission-record overlap u ranges over [0,1]. The ground and
excited amplitudes at each point are the projections of the same upstream
spinor. A trial with links equal to the actual complex hopping phases has
zero attachment penalty for any nonnegative individual bond stiffness.
"""
from itertools import product

from mpmath import iv

from mouth_collar_link_bounds import (
    R, ONE, ZERO, PI, TAU, dot, mesh, rotated, rounded_bound, spring
)


def conjugate(a):
    return (a[0], -a[1])


def multiply(a,b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def spinors(n,h):
    x,y,z = n
    # First spin component is real for the global source |up>.
    excited = ((1-z)/2, -x/2,-y/2)
    vacuum = (1-(1-h)*excited[0], (1-h)*x/2,(1-h)*y/2)
    return vacuum, excited


def inner(a,b):
    return (dot(a,b), a[1]*b[2]-a[2]*b[1])


def edge(n_i,n_j,s_i,s_j,t):
    vi,bi = s_i
    vj,bj = s_j
    vac = inner(vi,vj)
    factor = (1-t)*(1+dot(n_i,n_j))/2
    emitted = inner(bi,bj)
    common = tuple(a+factor*b for a,b in zip(vac,emitted))
    assert min(vac[0].a, common[0].a) > 0
    return vac, common


def loop_cost(links):
    # All u share the same nonnegative Bernstein weights. Every product
    # coefficient with m common-record factors is summed before measuring
    # its argument. Its binomial normalization is positive and immaterial.
    d = len(links)
    coefficients = [[ZERO,ZERO] for _ in range(d+1)]
    for bits in product((0,1), repeat=d):
        z = (ONE,ZERO)
        for pair,bit in zip(links,bits):
            z = multiply(z,pair[bit])
        total = coefficients[sum(bits)]
        total[0] += z[0]
        total[1] += z[1]
    assert all(x>0 for x,y in coefficients)
    max_angle = max(abs(iv.atan2(y,x)).b for x,y in coefficients)
    return spring(max_angle,4).b


def calculate(level,rotation,eta):
    verts,faces = mesh(level)
    axes = [tuple(eta*x for x in n) for n in rotated(verts,rotation)]
    edges = sorted({tuple(sorted((i,j))) for f in faces for i,j in zip(f,f[1:]+f[:1])})
    r = iv.cos(iv.sqrt(R(2))*R(".33"))**2
    surface, side = ZERO,ZERO
    previous = {}
    for k in range(7):
        t = r**k
        rays = [spinors(n,iv.sqrt(t)) for n in axes]
        links = {(i,j): edge(axes[i],axes[j],rays[i],rays[j],t) for i,j in edges}
        for (i,j),pair in list(links.items()):
            links[j,i] = tuple(conjugate(z) for z in pair)
        surface += sum(loop_cost([links[i,j] for i,j in zip(f,f[1:]+f[:1])]) for f in faces)
        if k:
            # Same axis and two different encounter counts give real positive
            # vacuum overlap and a nonnegative matched-emission contribution.
            # The two radial links consequently drop from the quadrilateral.
            side += sum(loop_cost([links[i,j],previous[j,i]]) for i,j in edges)
        previous = links
    total = surface+side
    tie = spring(PI,4)
    cell = 5*spring(TAU/5,4)
    tail = PI/32
    all_cuts = 7*len(faces)*spring(TAU/len(faces),4)
    floor = min(tie.a,cell.a,(tail+all_cuts).a)
    print(f"spinblind level={level} rotation={','.join(rotation)} eta={eta:+} u=[0,1] "
          f"tangent_upper={rounded_bound(surface,8,True):.8f} "
          f"side_upper={rounded_bound(side,8,True):.8f} "
          f"neutral_trial_upper={rounded_bound(total,8,True):.8f} "
          f"monopole_tail_lower={rounded_bound(tail,8):.8f} "
          f"any_charged_cut_or_tie_lower={rounded_bound(floor,8):.8f}",flush=True)
    assert total < floor


if __name__ == "__main__":
    for level in (2,3):
        for rotation in ((".26","-.19",".17"),("-.13",".24",".31")):
            for eta in (1,-1):
                calculate(level,rotation,eta)
