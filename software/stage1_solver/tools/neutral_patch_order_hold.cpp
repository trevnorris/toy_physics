// Overdamped time stepping for neutral_patch_order_hold.py. Python supplies
// each Wiener path, the archived pair-potential parameters and the drive.
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

extern "C" int assay(int count, int layered, int steps, int stride, double dt,
                     double softening, double gap, double pitch,
                     const double *within, const double *between,
                     const double *initial, const double *noise, const double *drive,
                     double *force, double *final, double *clamp_error) {
    const int n = 2 * count;
    std::vector<double> xy(initial, initial + 2*n), gradient(2*n);
    std::array<double, 2> base{};
    for (int i = 0; i < n; ++i) {
        const int plane = i / count;
        base[plane] += initial[2*i] / count;
        xy[2*i] += (plane - .5) * drive[0];
    }
    *clamp_error = 0;
    for (int step = 0; step <= steps; ++step) {
        std::fill(gradient.begin(), gradient.end(), 0.);
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                const bool cross = i < count && j >= count;
                const double *well = cross ? between : within;
                const double x = xy[2*j] - xy[2*i];
                const double y = xy[2*j+1] - xy[2*i+1];
                const double s = x*x + y*y;
                double r2 = s, factor = 1.;
                if (cross) {
                    r2 = gap*gap + (layered ? s*s/(pitch*pitch) : s);
                    factor = softening * (layered ? 2.*s/(pitch*pitch) : 1.);
                }
                const double r = std::sqrt(r2);
                if (!(r > 0.)) return 1;
                const double weight = std::exp(-well[1]*(r-well[2]));
                const double c = 2.*well[0]*well[1]*factor*weight*(1.-weight)/r;
                gradient[2*i] -= c*x;
                gradient[2*j] += c*x;
                gradient[2*i+1] -= c*y;
                gradient[2*j+1] += c*y;
            }
        }
        if (step % stride == 0) {
            double reaction = 0.;
            std::array<double, 2> mean{};
            for (int i = 0; i < n; ++i) {
                const int plane = i / count;
                reaction += (plane-.5)*gradient[2*i];
                mean[plane] += xy[2*i]/count;
            }
            if (!std::isfinite(reaction)) return 2;
            force[step/stride] = reaction;
            for (int plane = 0; plane < 2; ++plane) {
                const double expected = base[plane] + (plane-.5)*drive[step];
                *clamp_error = std::max(*clamp_error, std::abs(mean[plane]-expected));
            }
        }
        if (step == steps) break;
        const double *increment = noise + step*2*n;
        std::array<double, 2> mean{};
        for (int i = 0; i < n; ++i) {
            xy[2*i] += -dt*gradient[2*i] + increment[2*i];
            xy[2*i+1] += -dt*gradient[2*i+1] + increment[2*i+1];
            mean[i/count] += xy[2*i]/count;
        }
        // Orthogonal constraint on both plane means in x. The individual x
        // residuals and every y degree of freedom rearrange at unit mobility.
        for (int i = 0; i < n; ++i) {
            const int plane = i / count;
            xy[2*i] += base[plane] + (plane-.5)*drive[step+1] - mean[plane];
        }
    }
    std::copy(xy.begin(), xy.end(), final);
    return 0;
}
