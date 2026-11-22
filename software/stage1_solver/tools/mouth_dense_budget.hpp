// Discrete momentum balance for mobile, unit-mass sites in a periodic 4D coupon.
// The plus volume is [0,L/2) in w. Each call receives integrator impulses
// at the actual kick position or the unwrapped endpoints of a constant-v drift.
#pragma once

#include <array>
#include <cmath>
#include <iomanip>
#include <ostream>

namespace dense_budget {

using Vec = std::array<double, 4>;

inline int crossings(double from, double to, double cut, double length) {
    return static_cast<int>(std::floor((to - cut) / length)
                          - std::floor((from - cut) / length));
}

inline bool inside(double w, double length) {
    return std::floor(w / length) != std::floor((w - length / 2) / length);
}

inline void add(Vec& out, const Vec& in, double scale = 1.) {
    for (int k = 0; k < 4; ++k) out[k] += scale * in[k];
}

inline void write_vector(std::ostream& stream, const Vec& v) {
    stream << '[' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << ']';
}

struct Snapshot {
    long long count = 0;
    Vec p{};

    void record(const Vec& velocity) {
        ++count;
        add(p, velocity);
    }
};

struct Cut {
    long long count_pos = 0, count_neg = 0;
    Vec p_pos{}, p_neg{}, pair_into_plus{};

    void record(int turns, const Vec& v) {
        if (turns > 0) {
            count_pos += turns;
            add(p_pos, v, turns);
        } else if (turns < 0) {
            count_neg -= turns;
            add(p_neg, v, -turns);
        }
    }

    void json(std::ostream& os) const {
        os << "{\"count_pos\":" << count_pos << ",\"count_neg\":" << count_neg
           << ",\"p_pos\":";
        write_vector(os, p_pos);
        os << ",\"p_neg\":";
        write_vector(os, p_neg);
        os << ",\"mobile_pair_impulse_into_plus\":";
        write_vector(os, pair_into_plus);
        os << '}';
    }
};

struct Region {
    Snapshot initial{}, final{};
    Vec body{}, bath{}, wall{};

    void json(std::ostream& os, const Vec& advected_in, const Vec& pairs_in,
              long long advected_count) const {
        Vec residual = final.p;
        add(residual, initial.p, -1);
        add(residual, body, -1);
        add(residual, bath, -1);
        add(residual, wall, -1);
        add(residual, advected_in, -1);
        add(residual, pairs_in, -1);
        os << "{\"initial_n\":" << initial.count << ",\"final_n\":" << final.count
           << ",\"count_residual\":" << final.count - initial.count - advected_count
           << ",\"initial_p\":";
        write_vector(os, initial.p);
        os << ",\"final_p\":";
        write_vector(os, final.p);
        os << ",\"body\":";
        write_vector(os, body);
        os << ",\"bath\":";
        write_vector(os, bath);
        os << ",\"wall_on_mobile\":";
        write_vector(os, wall);
        os << ",\"advected_in\":";
        write_vector(os, advected_in);
        os << ",\"mobile_pair_in\":";
        write_vector(os, pairs_in);
        os << ",\"residual\":";
        write_vector(os, residual);
        os << '}';
    }
};

struct Budget {
    Region global{}, plus{};
    Cut sheet{}, middle{};

    void start(int, double w, const Vec& v, double length) {
        global.initial.record(v);
        if (inside(w, length)) plus.initial.record(v);
    }

    void finish(int, double w, const Vec& v, double length) {
        global.final.record(v);
        if (inside(w, length)) plus.final.record(v);
    }

    void drift(double old_w, double new_w, const Vec& v, double length) {
        sheet.record(crossings(old_w, new_w, 0., length), v);
        middle.record(crossings(old_w, new_w, length / 2, length), v);
    }

    void bath(double w, const Vec& impulse, double length) {
        add(global.bath, impulse);
        if (inside(w, length)) add(plus.bath, impulse);
    }

    void kick_mobile(double w, const Vec& body_impulse, const Vec& wall_impulse,
                     double length) {
        add(global.body, body_impulse);
        add(global.wall, wall_impulse);
        if (inside(w, length)) {
            add(plus.body, body_impulse);
            add(plus.wall, wall_impulse);
        }
    }

    // F is the impulse on i from j, n is oriented along the image bond i->j.
    // Crossing the sheet in +w puts j on the inside, so its impulse is -F.
    void pair_kick(int n0, int nh, const Vec& impulse_on_i) {
        add(sheet.pair_into_plus, impulse_on_i, -n0);
        add(middle.pair_into_plus, impulse_on_i, nh);
    }

    void json(std::ostream& os) const {
        Vec advection = sheet.p_pos;
        add(advection, sheet.p_neg, -1);
        add(advection, middle.p_pos, -1);
        add(advection, middle.p_neg);
        Vec pair = sheet.pair_into_plus;
        add(pair, middle.pair_into_plus);
        Vec fluid_on_wall{};
        add(fluid_on_wall, global.wall, -1);
        os << std::setprecision(17) << "{\"sheet\":";
        sheet.json(os);
        os << ",\"middle\":";
        middle.json(os);
        os << ",\"global\":";
        global.json(os, Vec{}, Vec{}, 0);
        os << ",\"plus\":";
        plus.json(os, advection, pair, sheet.count_pos - sheet.count_neg
                                        - middle.count_pos + middle.count_neg);
        os << ",\"mobile_on_sheet\":";
        write_vector(os, fluid_on_wall);
        os << '}';
    }
};

}  // namespace dense_budget
