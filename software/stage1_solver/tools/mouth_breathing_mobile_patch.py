"""Drive the half-order surface of a flowing planar Route A layer with both acoustic baths open.

This is a material response per unit area, in mouth units. It imposes the measured mean
normal current on the flat sample, not a pressure or flux source on the curved mouth.
"""
import argparse
import hashlib
import json

import numpy as np
from scipy import sparse
from scipy.integrate import solve_bvp
from scipy.linalg import null_space
from scipy.sparse.linalg import spsolve
from scipy.special import expit

from mouth_live_waist_fields import ORDER_INERTIA, WIDTH
from mouth_live_waist_order import WELL
from mouth_poly_shallow_release import ROOT, emit

KAPPA = 2 * WELL * WIDTH**2


def preparation(radius):
    record = ROOT/'research/pde_ledger_v3/_measurements/mouth_gallery_action.json'
    run = json.loads(record.read_bytes())['runs']['mouth_action_breathing01/run']
    raw = run['literal_stdout'].encode()
    if hashlib.sha256(raw).hexdigest() != run['stdout_sha256'] or run['child_exit'] or run['guard_exit']:
        raise ValueError('Unverified reference')
    rows = [r for r in map(json.loads, raw.splitlines()) if r['n'] == 72 and r['radius'] == radius]
    shape, = [r for r in rows if r['kind'] == 'setup']
    mode, = [r for r in rows if r['kind'] == 'three_shape_trial' and r['index'] == 0]
    return shape, mode['freq'], run['stdout_sha256']


def stationary_front(j, end):
    # Density one at infinity, moving fluid and a stationary centered front.
    mu = j**2/2
    def rhs(z, y):
        f, df, c, g = y
        u = j/f**2
        pi = ORDER_INERTIA*u*g
        ddf = 2*f*((f**8-1)/4 + u*u/2 + u*pi*g/f**2 - mu)
        ddc = (2*WELL*c*(1-c)*(1-2*c) - 4*ORDER_INERTIA*u*u*g*df/f)/(KAPPA-ORDER_INERTIA*u*u)
        return np.array([df, ddf, g, ddc])
    def bc(left, right):
        return np.array([left[1], left[2]-.5, right[0]-1, right[2]])
    z = np.unique(np.r_[np.linspace(0, 24*WIDTH, 241), np.linspace(24*WIDTH, end, 201)])
    c = expit(-z/WIDTH)
    y = np.array([np.ones_like(z), np.zeros_like(z), c, -c*(1-c)/WIDTH])
    sol = solve_bvp(rhs, bc, z, y, tol=2e-8, max_nodes=40000)
    if not sol.success or abs(sol.y[1,-1]) > 3e-7 or abs(sol.y[3,-1]) > 1e-7:
        raise ArithmeticError('Unable to establish uniform incoming and receiving material')
    def fields(z):
        f, df, c, g = sol.sol(abs(z))
        c = np.where(z < 0, 1-c, c)
        df = np.sign(z)*df
        ddf = rhs(z, [f,df,c,g])[1]
        return f, ddf, c, g
    return fields


def evolution(z, fields, omega, j):
    # Perturbation variables: density amplitude delta(sqrt n), its normal derivative,
    # phase, absolute mass flux, shear order and delta(kappa grad chi - pi u).
    f, ddf, c, g = fields(z)
    u = j/f**2
    pi = ORDER_INERTIA*u*g
    Wcc = WELL*(12*(c-.5)**2-1)
    eye = np.eye(6)
    u, pi, g, f, ddf, Wcc = [np.atleast_1d(v)[:, None] for v in (u,pi,g,f,ddf,Wcc)]
    a, ap, phase, flux, order, traction = eye
    dn = 2*f*a
    dv = (flux-u*dn)/f**2
    dc = (traction-1j*omega*ORDER_INERTIA*u*order+(ORDER_INERTIA*u*g+pi)*dv)/(KAPPA-ORDER_INERTIA*u*u)
    dp = ORDER_INERTIA*(-1j*omega*order+g*dv+u*dc)
    extra = (u*g*dp + pi*g*dv + pi*u*dc)/f**2 - u*pi*g*dn/f**4
    return np.stack(np.broadcast_arrays(ap, ddf*a/f + 2*f*(u*dv+extra+f**6*dn-1j*omega*phase),
                        dv+(g*dp+pi*dc)/f**2-pi*g*dn/f**4,
                        1j*omega*dn, dc, Wcc*order-1j*omega*dp), axis=1)


def reservoir(omega, j, side):
    # Acoustic Bogoliubov propagating + evanescent branch; one decaying order branch.
    roots = np.roots([.25, 0., 1-j*j, 2*omega*j, -omega**2])
    modes = []
    for k in roots:
        if abs(k.imag) < 1e-9:
            k = k.real
            group = j + k*(1+k*k/2)/(omega-j*k)
            if group*side > 0:
                modes.append((True,k,group))
        elif side*k.imag > 0:
            modes.append((False,k,0))
    modes.sort(reverse=True, key=lambda x:x[0])
    if len(modes)!=2 or not modes[0][0]:
        raise ArithmeticError('Unresolved outgoing acoustic branches')
    columns = []
    for _, k, _ in modes:
        columns.append([1,1j*k, 2*(omega/k-j)/(1j*k),2*omega/k,0,0])
    order_roots = np.roots([KAPPA-ORDER_INERTIA*j*j,2*ORDER_INERTIA*j*omega,2*WELL-ORDER_INERTIA*omega**2])
    ko, = [k for k in order_roots if side*k.imag > 0]
    columns.append([0,0,0,0,1,1j*((KAPPA-ORDER_INERTIA*j*j)*ko+ORDER_INERTIA*j*omega)])
    basis = np.array(columns).T
    _, kp, vg = modes[0]
    # Outgoing grand-canonical Hamiltonian energy in uniform flowing density one.
    power_per_a2 = side*vg*2*(omega/kp-j)*(omega/kp)
    return basis, power_per_a2


def response(j, omega, end, refinement):
    fields = stationary_front(j,end)
    # Resolve separately the order layer, the density healing tail and both acoustic wavelengths.
    half = np.r_[np.linspace(0,24*WIDTH,120*refinement+1),
                 np.linspace(24*WIDTH,end,round((end-24*WIDTH)*18*refinement)+1)[1:]]
    z = np.r_[-half[::-1],half]
    mid = len(half)-1
    starts = np.delete(np.arange(z.size-1),mid)
    step = np.diff(z)[starts,None,None]
    mat = evolution(z,fields,omega,j)
    center = evolution((z[starts]+z[starts+1])/2,fields,omega,j)
    eye = np.eye(6)
    # Fourth-order three-point Lobatto collocation, using cubic Hermite midpoints.
    left = -eye-step/6*(mat[starts]+4*center@(.5*eye+step/8*mat[starts]))
    right = eye-step/6*(mat[starts+1]+4*center@(.5*eye-step/8*mat[starts+1]))
    ids = 6*np.arange(z.size)[:,None]+np.arange(6)
    rows = []; cols = []; data = []
    def block(rd, cd, matrix):
        rows.append(np.broadcast_to(rd[..., :, None],matrix.shape).ravel())
        cols.append(np.broadcast_to(cd[..., None, :],matrix.shape).ravel())
        data.append(matrix.ravel())
    block(ids[:len(starts)],ids[starts],left)
    block(ids[:len(starts)],ids[starts+1],right)
    row = 6*len(starts)
    ends = []
    for side, node in [(-1,0),(1,-1)]:
        basis, scale = reservoir(omega,j,side)
        block(np.arange(row,row+3), ids[node], null_space(basis.conj().T).conj().T)
        ends.append((basis,scale))
        row += 3
    block(np.arange(row,row+5),ids[mid],eye[:5])
    block(np.arange(row,row+5),ids[mid+1],-eye[:5])
    row += 5
    block(np.array([row]),ids[mid],eye[4:5])
    source = np.zeros(z.size*6)
    g0 = fields(0.)[-1].item()
    source[row] = -g0   # Unit displacement of the half-order surface.
    system = sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),
                               shape=(len(source),len(source))).tocsc()
    states = spsolve(system,source).reshape(-1,6)
    external = g0*(states[mid+1,5]-states[mid,5])
    work = .5*(external*(1j*omega)).real
    power = []
    for node,(basis,scale) in zip([0,-1],ends):
        coeff = np.linalg.lstsq(basis,states[node],rcond=None)[0]
        power.append(float(scale*abs(coeff[0])**2))
    impedance = external/(-1j*omega)
    error = abs(work-sum(power))/max(work,sum(power),1e-15)
    n0 = fields(0.)[0].item()**2
    # Relative linear mass crossing at the moving midplane includes the background flux gradient
    # (zero in this stationary patch) and minus n0 times the prescribed interface velocity.
    relative = states[mid,3]+1j*omega*n0
    return dict(frequency=omega, normal_current=j, half_domain=end, mesh_refinement=refinement,
                center_density=n0, center_order_slope=g0, force_re=external.real,force_im=external.imag,
                normal_impedance_re=impedance.real,normal_impedance_im=impedance.imag,
                power_per_unit_displacement_squared=work, toward_feed=power[0],toward_receiver=power[1],
                acoustic_work_relative_error=error, relative_current_re=relative.real,relative_current_im=relative.imag,
                peak_linear_order_per_unit_displacement=float(np.max(abs(states[:,4]))))


def run(radius):
    background,omega,sha=preparation(radius)
    j=background['current']/background['cap_three_area']
    emit(kind='patch_input',radius=radius,source_sha=sha,shape_hash=background['endpoint_sha256'],
         occupation='original degree-four background trial frequency; uniform translating planar layer', omega=omega,
         mean_flux=j,width=WIDTH,I=ORDER_INERTIA,kappa=KAPPA,well_scale=WELL)
    cases={}
    for current in (0.,j,-j):
        for end,refinement in [(14.,1),(21.,2)]:
            row=response(current,omega,end,refinement)
            emit(kind='two_sided_patch',**row)
            if row['acoustic_work_relative_error']>2e-5 or min(row['toward_feed'],row['toward_receiver'])<=0:
                raise ArithmeticError('Unresolved two-sided power balance')
            cases[current,end]=row
        a,b=cases[current,14.],cases[current,21.]
        for key in ('normal_impedance_re','normal_impedance_im','toward_receiver','toward_feed'):
            if abs(a[key]-b[key])>3e-4*max(abs(b[key]),1e-12):
                raise ArithmeticError('Window or spatial resolution affected response')
    pos, neg=cases[j,21.],cases[-j,21.]
    if abs(pos['toward_feed']/neg['toward_receiver']-1)>2e-5 or abs(pos['toward_receiver']/neg['toward_feed']-1)>2e-5:
        raise ArithmeticError('Reversal of the current failed to interchange outgoing ports')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,choices=(1,2),required=True)
    run(parser.parse_args().radius)
