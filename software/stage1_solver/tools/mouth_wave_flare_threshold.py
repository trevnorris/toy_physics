"""Estimate the material contrast required for the first bound wave.

At the bulk continuum edge the 4D radial s-wave, u=r phi, has the free
Dirichlet half-plane Green function. Its image potential fixes the additive
constant of the free-space logarithm. There is no outer reflecting wall.
"""
import argparse
import numpy as np
import yaml
from scipy import fft
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.special import expit

from mouth_wave_flare_work import ARCHIVE, EXIT, ROUNDING, emit, radial_profile, requirements


def critical(r, h, step, margin, p):
    nr = int(np.ceil((1.5*r+margin)/step))
    nw = int(np.ceil(((EXIT-ROUNDING)*h+2*margin)/step))
    radius = (np.arange(nr)[:,None]+.5)*step
    w = ROUNDING*h-margin+(np.arange(nw)[None,:]+.5)*step
    gate = (expit((r*radial_profile(w/h)[0]-radius)/(p['delta_r']/2)) *
            expit((w-ROUNDING*h)/(p['delta_parallel']/2))*
            expit((EXIT*h-w)/(p['delta_parallel']/2)))
    root = np.sqrt(gate)
    # Toeplitz product as zero-padded FFT on the odd radial extension. Padding
    # ensures none of the retained target points see a periodic image.
    shape = (fft.next_fast_len(4*nr), fft.next_fast_len(2*nw))
    i, j = (np.minimum(np.arange(n), n-np.arange(n)) for n in shape)
    distance = step*np.hypot(i[:,None],j[None,:])
    distance[0,0] = 1.
    log_r = np.log(distance)
    log_r[0,0] = np.log(step/2)+np.log(2)/2+np.pi/4-1.5
    kernel = fft.rfftn(-step**2*log_r/(2*np.pi*p['c_w']**2))
    def apply(vector):
        a = root*vector.reshape(nr,nw)
        odd = np.concatenate((-a[::-1],a),axis=0)
        potential = fft.irfftn(kernel*fft.rfftn(odd,s=shape),s=shape)
        return (root*potential[nr:2*nr,:nw]).ravel()
    op = LinearOperator((nr*nw,nr*nw),matvec=apply,dtype=float)
    val, vec = eigsh(op,k=1,which='LA',v0=root.ravel(),tol=1e-10)
    residual = np.linalg.norm(apply(vec[:,0])-val[0]*vec[:,0])
    if residual>1e-8:
        raise ValueError('Threshold eigenvector residual')
    return 1/val[0]


def main(step, margin):
    p=yaml.safe_load(ARCHIVE.read_text())['parameters']
    strength=p['mu_out']**2-p['mu_in']**2
    data=requirements()
    for r in (.5,1.):
        h=data[r][0]
        value=critical(r,h,step,margin,p)
        emit(test='threshold',step=step,margin=margin,radius=r,depth=h,
             contrast_threshold=value,archive_contrast=strength,
             equivalent_mu_out_if_mu_in_fixed=np.sqrt(p['mu_in']**2+value))
    # Hold physical depth fixed for this separate radius sweep.
    root=brentq(lambda r: critical(r,1.,step,margin,p)-strength,.5,1.,xtol=2e-6)
    emit(test='radius_threshold',step=step,margin=margin,depth=1.,
         contrast=strength,critical_radius=root)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step',type=float,required=True)
    parser.add_argument('--margin',type=float,required=True)
    args=parser.parse_args()
    main(args.step,args.margin)
