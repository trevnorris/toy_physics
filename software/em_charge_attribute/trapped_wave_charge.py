"""A trapped real mode, a degenerate transverse pair, and their electric readout."""

import sympy as sp


def check(label, expression, expected):
    value = sp.factor(sp.trigsimp(expression))
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise AssertionError(label)


t, a0, g = sp.symbols("t A0 g", real=True)
length, z, omega, inertia = sp.symbols("L z omega M", positive=True)
x, y, vx, vy = sp.symbols("x y vx vy", real=True)
shape = sp.sin(sp.pi*z/(2*length))
check("Dirichlet value", shape.subs(z, 0), 0)
check("Neumann slope", shape.diff(z).subs(z, length), 0)
check("mode squared norm along path", sp.integrate(shape**2, (z, 0, length)), length/2)
check("standing mode longitudinal flux factor", sp.im(shape*shape.diff(z)), 0)

# M includes the spatial norm and brane inertia. Both components are in a
# fixed tangent-plane basis. There is no coupling to A0 in the wave energy.
lag = inertia*(vx**2 + vy**2 - omega**2*(x**2 + y**2))/2
generator = x*lag.diff(vy) - y*lag.diff(vx)
check("generator of tangent polarization rotations", generator, inertia*(x*vy-y*vx))
check("native direct electric source", -lag.diff(a0), 0)
fx, fy = sp.symbols("fx fy", real=True)
rate = (generator.diff(x)*vx + generator.diff(y)*vy
        + generator.diff(vx)*(fx/inertia - omega**2*x)
        + generator.diff(vy)*(fy/inertia - omega**2*y))
check("change of generator with applied forces", rate, x*fy-y*fx)
trial = lag.subs({vx: vx+g*a0*y, vy: vy-g*a0*x}, simultaneous=True)
check("electric source if polarization rotation is gauged",
      -trial.diff(a0).subs(a0, 0), g*generator)

# Both circular solutions have the same standing profile. Their sum gives
# an ellipse, with opposite senses labelled in the same laboratory frame.
plus, minus = sp.symbols("P B", real=True)
field = sp.Matrix([(plus+minus)*sp.cos(omega*t), (plus-minus)*sp.sin(omega*t)])
state = dict(zip((x, y, vx, vy), (*field, *field.diff(t))))
energy = vx*lag.diff(vx) + vy*lag.diff(vy) - lag
ipos = inertia*omega*plus**2
ineg = inertia*omega*minus**2
check("total action from energy", energy.subs(state)/omega, ipos+ineg)
check("polarization from the standing pair", generator.subs(state), ipos-ineg)
check("single linear mode", generator.subs(state).subs(minus, plus), 0)
check("unchanged under normal reflection of tangent components",
      generator.subs(state)*shape.subs(z, -z)**2 - generator.subs(state)*shape**2, 0)

radius, vr, vz, gap = sp.symbols("R v_r v_z gap", positive=True)
freq = sp.sqrt(gap**2 + (vr/radius)**2 + (vz/length)**2)
right, left = sp.symbols("I_plus I_minus", nonnegative=True)
h = freq*(right+left)
support = -h.diff(radius)
check("radial load from adding to one action and removing from the other",
      support.diff(right)-support.diff(left), 0)
check("equal additions to both actions change directed action by",
      (right-left).diff(right)+(right-left).diff(left), 0)
for rad, i_plus, i_minus in [(1, 1, 0), (2, 2, 1), (2, 1, 2), (2, 3, 0)]:
    sample = {radius: rad, length: 1, vr: 2, vz: 1, gap: 0, right: i_plus, left: i_minus}
    readings = [sp.simplify(v.subs(sample)) for v in (h, support, right-left, -lag.diff(a0))]
    print(f"R={rad} actions=({i_plus},{i_minus}) [energy, radial load, polarization, native charge] = {readings}")
