"""Independent variational certificate. All omitted Gaussian probability has gate >= 0."""
from math import erf, exp, pi, sqrt, sin, atan
from scipy.special import expit

def shape(s):
    if s>1.2:
        u=(s-1.2)/1.85
        return 1+.5*u*u*(3-2*u)
    t=(s-.5)/.7
    neck=.5+.5*sin(atan(.25/(4*pi*.4)))
    return neck+(1-neck)*t**3*(10-15*t+6*t*t)+.49*t*t*(1-t)**3

def radius_cdf(r,a):
    return erf(r/a)-2*r*exp(-(r/a)**2)/(a*sqrt(pi))
def certificate(R,H):
    # Widths of the wave amplitude, not its density.
    ar,aw = R, 1
    w0=2*H
    cell_count=16
    if R not in (1,2,4): raise ValueError(R)
    edges=[.5+2.55*k/cell_count for k in range(cell_count+1)]
    rad=[1.50*R*k/cell_count for k in range(cell_count+1)]
    v=0.
    for s0,s1 in zip(edges,edges[1:]):
        wmass=(erf((H*s1-w0)/aw)-erf((H*s0-w0)/aw))/2
        axial=expit((H*s0-.5*H)/.025)*expit((3.05*H-H*s1)/.025)
        for r0,r1 in zip(rad,rad[1:]):
            rmass=radius_cdf(r1,ar)-radius_cdf(r0,ar)
            v+=rmass*wmass*axial*expit((R*shape(s0)-r1)/.025)
    kin=1.5/ar**2+.5/aw**2
    shifted=kin-6.21*v
    print('R',R,'H',H,'cellcount',cell_count,'gate_lower',v,'gradient',kin,'rayleigh_shift_upper',shifted,'frequency_upper',sqrt(6.25+shifted))
    assert shifted < 0
for radius,depth in [(1,1),(2,0.9947506011480172),(4,0.950896767876907)]:
    certificate(radius,depth)
