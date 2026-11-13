"""References from S15 definitions; uses neither engine parser nor comparator code."""
import sympy as s

def references():
    # Translating held sphere, oriented n = sigma*e. This integral uses dOmega.
    theta, phi = s.symbols('theta phi', real=True)
    e = s.Matrix([s.sin(theta)*s.cos(phi), s.sin(theta)*s.sin(phi), s.cos(theta)])
    moment = (e*e.T*s.sin(theta)).applyfunc(
        lambda value: s.integrate(value, (phi, 0, 2*s.pi), (theta, 0, s.pi)))
    print('angular second moment', moment)
    xdot = s.Matrix([2, -3, 5])
    eps = s.symbols('eps')
    for surf, b, sigma in [('inner', 2, -1), ('outer', 3, 1)]:
        # Independent instantaneous field choices at t=1, both profiles and changes live.
        rho = 5 + b**2 + eps*(7 - b)
        u = 3*b + eps*(b + 4)
        q = sigma*b**2*rho*u*(moment*xdot)[0]
        q0 = q.subs(eps, 0)
        response = s.diff(q, eps)
        # The selected drain2-only object has no X_1 or Xdot_1 entry.
        print(surf, 'transport_centre component_1 drain2 Q, geometric, response, v1part, v2part:',
              q0, 2*q0, response.subs(eps, 0), 0, q0)
        A,B,p,h = s.symbols('A B a b', real=True)
        radius = s.symbols('radius', positive=True)
        print(surf,'closed Q, E:',
              sigma*radius**(2+p+h)*A*B*(moment*xdot)[0],
              2*sigma*radius**(2+p+h)*A*B*(moment*xdot)[0])
        # Pi_o^(1j) affine at x = X2+b*e. Row j has arbitrary derivatives
        # with respect to all three independent coordinates; integral contracts
        # row j with column j. An off-diagonal 12->13 changes this contraction.
        M = s.Matrix([[2,3,5], [7,11,13], [17,19,23]])
        constant = s.Matrix([1,4,-2])
        other = -sigma*b**3*(M*moment).trace()
        changed = M.copy()
        changed[1,:] = M[2,:]
        bad = -sigma*b**3*(changed*moment).trace()
        print(surf, 'other last row independent; Q, sign flip, component 12<-13:', other,-other,bad)
    # Full shell rate dimension rho * velocity^2 * area.
    print('crossing dimension (density,length,time)', (1,4,-2))
    # Surface/radial force Q family has r1^k moment, smooth as d>R.
    # r1^2 = d^2+b^2-2*d*b*cos(theta), dOmega=2*pi*d(cos theta).
    b,d,p=s.symbols('r d p',positive=True)
    print('angular branch k=-2 from r1^(-2)', s.pi/(d*b)*s.log((d+b)/(d-b)))
    print('angular branch k=-2 d=5 r=2', s.pi*s.log(s.Rational(7,3))/10)
    print('angular k=0',4*s.pi)
if __name__=='__main__':
    references()
