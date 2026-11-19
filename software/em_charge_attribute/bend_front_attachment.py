"""One planar bend, with separate conversion and reflecting-wave assumptions."""

import sympy as sp


def show(label, value):
    print(f"{label} = {value}", flush=True)


def check(label, value, comparison):
    error = sp.simplify(value - comparison)
    show(label, (value, comparison, error))
    if error != 0:
        raise ValueError(label)


# Positive q points from the ordered side into bulk. Keep the flux relative
# to the rim fixed while its lab displacement z changes.
q, t, follow = sp.symbols("q t lambda", real=True)
n, V, L, width = sp.symbols("n V L delta", positive=True)
z = sp.Function("z")(t)
front = L + follow * z
argument = q - front
h = sp.Function("h")
chi = h(argument)
speed = V + sp.diff(z, t)
order_flux = sp.Function("J")(q, t)
source = sp.diff(n*chi, t) + sp.diff(n*speed*chi + order_flux, q)
gradient = sp.diff(chi, q)
required_speed = sp.factor((source - sp.diff(order_flux, q)) / (n*gradient))
show("material flux relative to rim", sp.simplify(n*(speed - z.diff(t))))
show("net local order source / (n chi_q)", required_speed)
for fraction in (0, 1):
    show(f"lambda={fraction} source coefficient", required_speed.subs(follow, fraction))
local_speed = sp.solve(sp.diff(chi, t) + sp.Symbol("u_front")*gradient,
                       sp.Symbol("u_front"))[0]
check("front speed from order balance", sp.simplify(speed-required_speed), local_speed)
# Unit drop for a smooth layer with finite width. No endpoint flux is
# absorbed into the material intake when integrating the order balance.
s = sp.symbols("s", real=True)
profile = 1/(1 + sp.exp(s/width))
drop = sp.limit(profile, s, -sp.oo) - sp.limit(profile, s, sp.oo)
effective_loss = sp.factor(n * drop * required_speed)
check("endpoint-corrected loss minus rim intake",
      sp.factor(effective_loss - n*V), n*(z.diff(t)-front.diff(t)))

# Independent end coordinates for shape derivatives at fixed mode action.
rim, tip = sp.symbols("rim tip", real=True)
length, c, action = sp.symbols("ell c I", positive=True)
x = sp.symbols("x", real=True)
k = sp.pi/(2*length)  # fundamental of the node / free-end interval
mode = sp.sin(k*x)
check("wave equation", sp.diff(mode,x,2), -k**2*mode)
show("node and free-end derivative", (mode.subs(x,0), sp.diff(mode,x).subs(x,length)))
strain_weight = sp.diff(mode,x)**2
centroid = sp.simplify(sp.integrate(x*strain_weight, (x,0,length))
                       / sp.integrate(strain_weight, (x,0,length)))
fraction = sp.factor(centroid/length)
show("strain energy centroid fraction from rim", fraction)
omega = c*k.subs(length, tip-rim)
energy = action*omega
fr = sp.factor(-sp.diff(energy,rim)/energy)
fy = sp.factor(-sp.diff(energy,tip)/energy)
show("rim and front forces / energy", (fr,fy))
check("end force sum / energy", fr+fy, 0)
bent = {tip: L+follow*z, rim: z}
ell = (tip-rim).subs(bent)
omega_bent = omega.subs(bent)
show("strain centroid change per rim displacement", sp.diff(z+fraction*ell,z))
check("wave work per energy per rim displacement",
      sp.factor((fr+follow*fy).subs(bent)), -sp.diff(omega_bent,z)/omega_bent)
check("frequency drift from corrected conversion rate",
      sp.factor(sp.diff(omega_bent,t)/omega_bent), (effective_loss/n-V)/ell)
