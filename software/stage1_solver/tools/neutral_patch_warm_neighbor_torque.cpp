// Mobile-site Langevin paths with an applied torque conjugate to the
// difference of two unwrapped Procrustes angles. Python supplies Wiener noise.
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

extern "C" int warm_turn(int count, int layered, int steps, int stride, double dt,
                         double softening, double gap, double pitch, double torque,
                         const double *within, const double *between, const int *group,
                         const double *initial, const double *reference, const double *noise,
                         double *angles, double *last_before, double *final, double *checks) {
    const int n = 2*count;
    std::vector<double> xy(initial, initial+2*n), drift(2*n);
    std::array<double,2> previous{}, unwrapped{}, initial_moment{};
    std::array<int,2> members{};
    for (int i = 0; i < n; ++i) {
        ++members[group[i]];
        initial_moment[group[i]] += reference[2*i]*reference[2*i]+reference[2*i+1]*reference[2*i+1];
    }
    checks[0] = 1.; // smallest normalized Procrustes overlap
    checks[1] = checks[2] = checks[3] = 0.; // max angle step, force and torque defect
    for (int k = 0; k <= steps; ++k) {
        std::array<double,2> a{}, b{}, cx{}, cy{}, moment{}, group_torque{}, fx{}, fy{};
        for (int i = 0; i < n; ++i) {
            const int g = group[i];
            cx[g] += xy[2*i]/members[g];
            cy[g] += xy[2*i+1]/members[g];
        }
        for (int i = 0; i < n; ++i) {
            const int g = group[i];
            const double x = xy[2*i]-cx[g], y = xy[2*i+1]-cy[g];
            const double rx = reference[2*i], ry = reference[2*i+1];
            a[g] += rx*x+ry*y;
            b[g] += rx*y-ry*x;
            moment[g] += x*x+y*y;
        }
        for (int g = 0; g < 2; ++g) {
            const double overlap = std::hypot(a[g], b[g])/std::sqrt(moment[g]*initial_moment[g]);
            if (!std::isfinite(overlap) || overlap < 1.e-8) return 1;
            checks[0] = std::min(checks[0], overlap);
            const double phase = std::atan2(b[g],a[g]);
            const double delta = std::remainder(phase-previous[g],2.*std::acos(-1.));
            unwrapped[g] += delta;
            previous[g] = phase;
            checks[1] = std::max(checks[1],std::abs(delta));
            if (k % stride == 0) angles[2*(k/stride)+g] = unwrapped[g];
        }
        if (k == steps) break;
        if (k == steps-1) std::copy(xy.begin(),xy.end(),last_before);
        // Gradient of the fitted angle at the *current* configuration. The
        // generalized load gives +torque to the interior, -torque to its rim.
        for (int i = 0; i < n; ++i) {
            const int g = group[i];
            const double rx = reference[2*i], ry = reference[2*i+1];
            const double scale = (g == 0 ? 1. : -1.)/(a[g]*a[g]+b[g]*b[g]);
            const double gx = scale*(-a[g]*ry-b[g]*rx), gy = scale*(a[g]*rx-b[g]*ry);
            drift[2*i] = torque*gx;
            drift[2*i+1] = torque*gy;
            fx[g] += gx;
            fy[g] += gy;
            group_torque[g] += (xy[2*i]-cx[g])*gy-(xy[2*i+1]-cy[g])*gx;
        }
        for (int g = 0; g < 2; ++g) {
            checks[2] = std::max(checks[2],std::hypot(fx[g],fy[g]));
            checks[3] = std::max(checks[3],std::abs(group_torque[g]-(g == 0 ? 1. : -1.)));
        }
        // All native Morse pairs, common mobility one. Neither plane is held.
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                const bool cross = i < count && j >= count;
                const double *well = cross ? between : within;
                const double x = xy[2*j]-xy[2*i], y = xy[2*j+1]-xy[2*i+1];
                const double s = x*x+y*y;
                const double r = std::sqrt(cross ? gap*gap+(layered ? s*s/(pitch*pitch) : s) : s);
                if (!(r > 0.)) return 2;
                const double weight = std::exp(-well[1]*(r-well[2]));
                const double projection = cross ? softening*(layered ? 2.*s/(pitch*pitch) : 1.) : 1.;
                const double c = 2.*well[0]*well[1]*projection*weight*(1.-weight)/r;
                drift[2*i] += c*x;
                drift[2*i+1] += c*y;
                drift[2*j] -= c*x;
                drift[2*j+1] -= c*y;
            }
        }
        for (int i = 0; i < 2*n; ++i) xy[i] += dt*drift[i]+noise[k*2*n+i];
    }
    std::copy(xy.begin(),xy.end(),final);
    return 0;
}
