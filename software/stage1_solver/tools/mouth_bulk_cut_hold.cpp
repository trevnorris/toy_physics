// Use the archived thermal material's all-image kernel and unwrapped coordinates.
#include "mouth_bulk_thermal_cell.cpp"

// Hold a single coordinate Q = sum_i D_i t.x_i / W, where
// D_i = label_i - mean(label), W = sum_i D_i^2, and t is a unit tangent.
// The mobility of every Cartesian site coordinate is one.  The constraint
// moves site i by D_i t (Qtarget-Qtrial). Its total force is zero and its
// generalized force is W (Qtarget-Qtrial)/dt.
extern "C" int cut_advance(void* ptr, int steps, double dt, const double* noise,
                          const double* d, const double* tangent, double w,
                          double target, double* x, double* stats) {
    try {
        Cell& cell = *static_cast<Cell*>(ptr);
        const int count = 4*cell.n;
        std::vector<double> gradient(count), delta(count);
        std::fill_n(stats, 8, 0.);
        for (int step=0; step<steps; ++step) {
            if (!std::isfinite(cell.evaluate<false>(x, gradient.data(), nullptr))) return step+1;
            double trial_q = 0., old_q = 0., noise_moment = 0., gradient_moment = 0.;
            for (int i=0; i<cell.n; ++i) {
                for (int k=0; k<4; ++k) {
                    const int j = 4*i+k;
                    const double weight = d[i]*tangent[k], eta = noise[step*count+j];
                    delta[j] = eta-dt*gradient[j];
                    old_q += weight*x[j];
                    trial_q += weight*(x[j]+delta[j]);
                    noise_moment += weight*eta;
                    gradient_moment += weight*gradient[j];
                }
            }
            const double correction = target-trial_q/w;
            const double applied = w*correction/dt;
            const double thermal = -noise_moment/dt;
            const double conditional = applied-thermal;
            for (int i=0; i<cell.n; ++i) {
                for (int k=0; k<4; ++k) {
                    double& increment = delta[4*i+k];
                    increment += correction*d[i]*tangent[k];
                    if (!std::isfinite(increment) || std::abs(increment)>1.) return step+1;
                }
            }
            for (int j=0; j<count; ++j) x[j] += delta[j];
            stats[0] += applied;
            stats[1] += thermal;
            stats[2] += conditional;
            stats[3] += gradient_moment;
            stats[4] = std::max(stats[4], std::abs(conditional-gradient_moment));
            stats[5] = std::max(stats[5], std::abs(old_q/w-target));
            if (step == 0) stats[6] = conditional;
            stats[7] += thermal*thermal;
        }
        return 0;
    } catch (...) {return -1;}
}
