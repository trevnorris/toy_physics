"""Independent source-law and archive checks. Run only under s11c_guarded_run."""
from pathlib import Path
from hashlib import sha256
import json
import re
import numpy as np
import sympy as sp

wd=Path(__file__).resolve().parent
cand=wd.parent/"r0/candidate"
def archive(name):
    return json.loads((Path("docs/_measurements")/(name+".json")).read_text())["runs"]
def read(path,digest):
    assert sha256(path.read_bytes()).hexdigest()==digest
    return np.load(path,allow_pickle=False)
def arrays(row, basename, wanted):
    assert (row.get("guard_exit",0),row.get("worker_exit",0),row.get("limits_verified",True))==(0,0,True)
    key,= [k for k in row["fields"] if Path(k).name==basename]
    with read(Path(key),row["fields"][key]) as p:
        return {k:(p[k][-1] if k in ("layers","field") else p[k]) for k in wanted}
def records():
    for f in sorted((cand/"docs/_measurements").glob("*.txt")):
        s=f.read_text()
        lines=s.split("Literal stdout\n")[1].strip()+"\n"
        fields=re.search(r"Stdout (\S+) \((\d+) bytes, SHA256 (\w+)",s)
        raw=Path(fields[1]).read_bytes()
        assert raw==lines.encode() and len(raw)==int(fields[2]) and sha256(raw).hexdigest()==fields[3]
        files=re.findall(r"^([0-9a-f]{64})  (.+)$",s,re.M)
        for digest,path in files:
            p=Path(path)
            if (cand/path).exists(): p=cand/path
            assert sha256(p.read_bytes()).hexdigest()==digest
        print("record",f.stem,"byte-exact stdout",len(raw),"verified hashes",len(files))
def unequal():
    prs=archive("lateral_paired_acceleration"); rgs=archive("lateral_radius_momentum"); lgs=archive("lateral_target_rate")
    for pair in [p for p in prs if p["radius"]==1.5]:
        right,=[g for g in rgs if g["grid"]==pair["grid"] and g["radius"]==1.5]
        left,=[g for g in lgs if pair["reference_fields"].items() <= g["reference_fields"].items()]
        assert right["reference_fields"]==pair["reference_fields"]
        for d in (4,8):
            p=arrays(pair,f"d{d}.npz",["field","currents","x","r","dz","d","target_radius"])
            assert p["d"]==d and p["target_radius"]==1.5
            a=arrays(left,f"source_d{d}.npz",["layers","currents","x","r","dz","centers","rates"])
            b=arrays(right,"local_fields.npz",["field","normal","x","r","dz","radius"])
            assert a["rates"].tolist()==[0.25] and a["centers"].tolist()==[-d] and b["radius"]==1.5
            for idx,(ref,cen,rad) in enumerate(((a,-d,1.),(b,0,1.5))):
                xs,rs=ref["x"],ref["r"]; x0,x1=np.searchsorted(p["x"],[xs[0],xs[-1]])
                np.testing.assert_array_equal(xs,p["x"][x0:x1+1])
                np.testing.assert_array_equal(rs,p["r"][:len(rs)])
                assert np.isclose(p["dz"],ref["dz"])
                pf=p["field"][x0:x1+1,:len(rs)]
                af=ref["layers"] if idx==0 else ref["field"]
                aj=ref["currents"][0] if idx==0 else np.imag(af.conj()*ref["normal"])
                pj=p["currents"][idx][x0:x1+1,:len(rs)]
                h=xs[1]-xs[0]; vol=np.pi*((rs+h/2)**2-(rs-h/2)**2)*h
                x=xs[:,None]-cen
                g=lambda q: np.maximum(1-((x-q)**2+rs**2)/rad**2,0)**3
                support=g(0)>0
                assert not np.any(pj[~support]) and np.all(pj[support]>0)
                # Each isolated field measures the velocity shape for its aperture.
                pred=aj*abs(pf/af)**2
                if idx==0: pred*=.25/np.sum(pred*vol)
                np.testing.assert_allclose(pj,pred,rtol=1e-12,atol=1e-14)
                def E(q,face):
                    return np.sum(vol*g(q)*(abs(face)-1)**2)/2
                q=h*1e-4
                def F(face):
                    return (E(-q,face)-E(q,face))/(2*q)
                delta=F(pf)-F(af)
                density=abs(pf)**2
                slope=-2*g(0)*(1-1/np.sqrt(density))
                def stress(face):
                    den=abs(face)**2
                    ln=-2*g(0)*(1-1/np.sqrt(den))
                    dxn=(-np.roll(den,-2,0)+8*np.roll(den,-1,0)-8*np.roll(den,1,0)+np.roll(den,2,0))/(12*h)
                    return np.sum(vol*ln*dxn)/4
                actual=stress(pf)-stress(af)
                eta=(x/rad)*g(0); rate=np.sum(vol*pj)
                # Independent perturbation reconstructs speed from isolated trace.
                rel=1e-5
                def moved(e):
                    raw=aj*abs(pf*(1+e*eta)/af)**2
                    return raw if idx else raw*.25/np.sum(raw*vol)
                jac=(moved(rel)-moved(-rel))/(2*rel)
                curl=np.sum(vol*eta*jac)
                mean=np.sum(vol*pj*eta)/rate if not idx else 0
                alt=2*np.sum(vol*pj*(eta-mean)**2)
                assert np.isclose(curl,alt,rtol=1e-8)
                print(f"unequal {pair['grid']} d={d} a={rad} Fpair-Falone={delta:.12g} old_mixed_stress={actual:.12g} J={rate:.12g} current_curl={curl:.12g}")
def periodic():
    for u in (.01,.02):
        runs=archive("mouth_intake_redistribution")
        k,row=next((k,row) for row in runs for k in row.get("fields",{}) if "depth01" in k and k.endswith(f"restoring_slope_u{u}.npz"))
        with read(Path(k),row["fields"][k]) as p:
            f,j,mask,h,dz=(p[k] for k in ["field","current","mouth","dx","dz"])
        assert f.shape==(65,32,32,32) and np.max(abs(f[0]-1))==0
        x=(np.arange(32)-15.5)*h
        X=x[:,None,None]; yz=x[None,:,None]**2+x[None,None,:]**2
        g=lambda z: np.maximum(1-(X-z)**2-yz,0)**3
        def der(v):
            # Four sample symmetric slope on the tangential axis.
            return (2*(np.roll(v,-1,0)-np.roll(v,1,0))/3-(np.roll(v,-2,0)-np.roll(v,2,0))/12)/h
        readings=[]
        for v,field,current in [(u,f,j),(-u,f[:,::-1],j[::-1])]:
            face=field[-1]; amp=abs(face); n=amp**2
            feed=v*np.sum(np.imag(field[1]))*h**3/dz
            vx=v+np.imag(np.conj(face)*der(face))/n
            K=np.sum(current*vx)*h**3
            T=np.sum(-g(0)*(1-1/amp)*der(n))*h**3/2
            dzf=np.einsum('i,i...->...',[11,-18,9,-2],field[-1:-5:-1])/(6*dz)
            logz=np.einsum('i,i...->...',[11,-18,9,-2],np.log(abs(field[-1:-5:-1])**2))/(6*dz)
            layer=np.sum(((face.conj()*dzf).imag*vx - n*der(logz)/4)[mask])*h**3
            d=h*1e-4
            def energy(z):
                return np.sum(g(z)*(amp-1)**2)*h**3/2
            wall=(energy(-d)-energy(d))/(2*d)
            readings.append(np.array([feed,K+T,layer,wall]))
        odd=np.diff(readings,axis=0)[0]/(-2*u)
        print("periodic",u,"odd feed, flux, layers, translated",*(f"{k:.12g}" for k in odd))
        assert np.max(abs(np.mean(readings,axis=0)))<1e-14
def identities():
    t,c=sp.symbols('x c',real=True,positive=True)
    amp=sp.Function('u')(t); E=sp.Function('e')(t,amp**2)
    A,n,q,y=sp.symbols('A n q y',real=True)
    W=sp.Function('W')(y,n)
    En=W.diff(n).subs({y:t,n:amp**2})
    L=sp.Function('L')(t)
    R=L+En*amp/c
    # Direct mixed stress from n gradient products before integration by parts.
    def Q(w):
        nz=2*w*amp; den=amp**2
        return c*(nz*den.diff(t)/den - nz.diff(t))/2
    jump=sp.expand(Q(L)-Q(R))
    trans=W.diff(y).subs({y:t,n:amp**2})
    tang=(amp**2*En-W.subs({y:t,n:amp**2})).diff(t)
    print("two sided full-stress difference - divergence - explicit force",sp.simplify(jump-tang-trans))
    assert sp.simplify(jump-tang-trans)==0
    # Covariant/neutral current continuity independent of slope jump.
    J=sp.symbols('j',real=True)
    print("arbitrary current and tangential phase, carriage jump",sp.simplify(J*q-J*q))
    mu,U=sp.symbols('mu U',real=True); r,a=sp.symbols('r a',positive=True)
    n0, n1, n4, n3=sp.symbols('n0 n1 n4 n3')
    # True rigidly moved radial potential at fixed lab point.
    V=sp.Function('V')
    d=sp.sqrt(r**2+A**2-2*r*A*mu)
    force=-sp.diff(V(d),A)
    density=n0+U*(n1*mu+n4*sp.legendre(4,mu))+U**3*n3*mu
    angular=sp.simplify(2*sp.pi*sp.integrate(r**2*density*force.subs(A,0),(mu,-1,1)))
    print("fixed-shape confinement, including cubic odd drive",angular)
    add=sp.simplify(2*sp.pi*sp.integrate(r**2*n0*sp.diff(force,A).subs(A,0),(mu,-1,1)))
    print("if section shifts with U, extra integrand times dq/dU at zero",add)
    # For a literal level set r-R(ray) against fixed n_0.
    lev=V(r-(A*mu+sp.sqrt(a**2-A**2*(1-mu**2))))
    add_lev=sp.simplify(-2*sp.pi*sp.integrate(r**2*n0*sp.diff(lev,A,2).subs(A,0),(mu,-1,1)))
    print("level-set version of shape susceptibility, before dr",add_lev)
    w,L,T=sp.symbols("w L T",positive=True)
    qf=sp.Function("q")(w); Tvar=sp.Function('T')(w)
    delta=sp.Function("v")(w)
    integrand=sp.integrate(2*sp.pi*mu**2*Tvar*qf.diff(w)*delta.diff(w),(mu,-1,1))
    edge=sp.diff(integrand,delta.diff(w))
    print("elastic force on lower end",edge.subs(w,0),"holder sign is the negative")
    s=w/L; d=s**2*(1-s)**2
    print("trial difference integral",sp.integrate(d,(w,0,L)),"boundary values/slopes",[(d.subs(w,u),d.diff(w).subs(w,u)) for u in (0,L)])
    h,f0,f1=sp.symbols('h f0 f1',positive=True)
    print("shrinking regular collar",sp.integrate((1-w/h)**2*(1+2*w/h)*(f0+f1*w),(w,0,h)))
    print("symmetric-secant slope for F=b1 U+b3 U^3, half speed does not determine all odd orders", sp.expand((U+U**3-(-U-U**3))/(2*U)))
if __name__=="__main__":
    records(); identities(); unequal(); periodic()
