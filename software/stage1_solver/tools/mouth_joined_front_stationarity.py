"""Resolve terminal field stationarity when total-energy subtraction loses significant digits."""
import argparse
from pathlib import Path
import numpy as np
from mouth_joined_front_audit import run
from mouth_supercore_joined_front import Joined,check_guard,emit

def solve(self,rest):
    base=self.seed.copy()
    if rest:base[1]=0.
    if self.warm_path:
        check_guard(self.warm_path)
        with np.load(self.warm_path) as data:
            old=data['rest' if rest else 'fields']
        if old.shape!=base.shape:raise ValueError('Input warm field topology mismatch')
        base[self.masks]=old[self.masks]
    x=self.vars(base)
    E,all_g=self.evaluate(base)
    g=self.vars(all_g);history=[]
    small=0;restarts=0
    for iteration in range(700):
        z=self.norm(g);err=max(abs(z))
        if err<2e-7:
            return self.unpack(x,base),dict(iterations=iteration,energy=E,stationarity=err,
                                             gradient_integral_steps=small,restarts=restarts)
        q=g.copy();coeffs=[]
        for dx,dg,scale in reversed(history):
            c=scale*(dx@q);q-=c*dg;coeffs.append(c)
        d=self.norm(q)
        for (dx,dg,scale),c in zip(history,reversed(coeffs)):
            d+=dx*(c-scale*(dg@d))
        success=False
        for which,search in enumerate((d,z)):
            slope=g@search
            if slope<=0:continue
            alpha=min(1,.15/max(abs(search)))
            while alpha>1.e-9:
                step=-alpha*search
                new=self.unpack(x+step,base)
                En,gn=self.evaluate(new)
                row=self.vars(gn)
                predicted=alpha*slope
                if predicted<1.e-9*max(1.,abs(E)):
                    # Integrate the directional gradient along this actual step
                    # to avoid subtracting large, nearly identical field energies.
                    # Simpson and trapezoid agreement also checks the local expansion.
                    mid=self.unpack(x+.5*step,base)
                    _,gm=self.evaluate(mid)
                    edge=(g+row)@step/2
                    delta=(edge+2*(self.vars(gm)@step))/3
                    accept=(delta < -.0001*predicted and abs(edge-delta)<.02*predicted
                            and En<=E+1.e-10*max(1.,abs(E)))
                else:
                    accept=En<=E-.0001*predicted
                if accept:
                    success=True;break
                alpha*=.5
            if success:
                if which:history=[];restarts+=1
                break
        if not success:raise ArithmeticError(f'No energy-descending direction at {iteration}, {err}')
        dg=row-g
        curvature=step@dg
        if curvature>max(0,1.e-11*np.linalg.norm(step)*np.linalg.norm(dg)):
            history.append((step,dg,1/curvature));history=history[-12:]
        if predicted<1.e-9*max(1.,abs(E)):small+=1
        E,g,x=En,row,x+step
    raise ArithmeticError(f'Unresolved stationarity {err}')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('--grid',type=int,nargs=3,required=True)
    p.add_argument('--winding',type=int,choices=(0,1),required=True)
    p.add_argument('--work',action='store_true');p.add_argument('--warm',type=Path);p.add_argument('--save',type=Path)
    a=p.parse_args()
    Joined.solve=solve;Joined.warm_path=a.warm
    emit(kind='method',name='directional gradient work for small energy differences',warm=str(a.warm))
    run(a)
