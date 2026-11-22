// Four-dimensional periodic WITHIN material in a fixed row-vector lattice.
// Coordinates and Brownian increments are unwrapped Cartesian vectors.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace {
using Vec = std::array<double, 4>;
using Tile = std::array<int, 4>;
using Matrix = std::array<Vec, 4>;

struct Radial {
    long double energy, first, second;  // derivatives with respect to r^2
};
struct Polynomial {
    // Quintic Hermite energy in t = 256*(r^2 - left), with its first five powers.
    double e0, c1, c2, c3, c4, c5;
};
struct Potential {
    double cutoff, depth, alpha, rest;
    static constexpr int scale = 256;
    std::vector<Polynomial> table;
    int switch_bin, end_bin;

    Radial exact(long double x) const {
        if (x >= (long double)cutoff*cutoff) return {};
        const long double r = std::sqrt(x), w = std::exp(-alpha*(r-rest));
        long double u = depth*w*(w-2), du = 2*alpha*depth*w*(1-w);
        long double ddu = 2*alpha*alpha*depth*w*(2*w-1);
        if (r > cutoff-.25) {
            const long double t = 4*(cutoff-r);
            const long double s = t*t*t*(10+t*(-15+6*t));
            const long double ds = -120*t*t*(1-t)*(1-t);
            const long double dds = 960*t*(1-t)*(1-2*t);
            ddu = ddu*s+2*du*ds+u*dds;
            du = du*s+u*ds;
            u *= s;
        }
        return {u, du/(2*r), (ddu-du/r)/(4*x)};
    }

    Potential(double rc, const double* p) :
        cutoff(rc), depth(p[0]), alpha(p[1]), rest(p[2]) {
        const double start = (rc-.25)*(rc-.25)*scale, end = rc*rc*scale;
        // Evaluate a bin straddling a switch endpoint analytically. Both
        // endpoints are exact grid points at the specified cutoff of 2.5.
        switch_bin = start == std::floor(start) ? -1 : int(start);
        end_bin = end == std::floor(end) ? -1 : int(end);
        table.resize(std::ceil(end));
        const long double h = 1.L/scale;
        for (size_t k=scale/4; k<table.size(); ++k) {
            const Radial a = exact(k*h), b = exact((k+1)*h);
            const long double c1 = h*a.first, c2 = h*h*a.second/2;
            const long double e = b.energy-a.energy-c1-c2;
            const long double g = h*b.first-c1-2*c2;
            const long double h2 = h*h*b.second-2*c2;
            table[k] = {double(a.energy), double(c1), double(c2), double(10*e-4*g+h2/2),
                        double(-15*e+7*g-h2), double(6*e-3*g+h2/2)};
        }
    }

    template<bool read>
    double force_prefactor(double x, double& energy) const {
        const double bin = x*scale;
        const int k = int(bin);
        if (k < scale/4 || k == switch_bin || k == end_bin) {
            const Radial z = exact(x);
            if constexpr (read) energy = z.energy;
            return 2*z.first;
        }
        const double t = bin-k;
        const auto& c = table[k];
        if constexpr (read)
            energy = c.e0+t*(c.c1+t*(c.c2+t*(c.c3+t*(c.c4+t*c.c5))));
        return (2*scale)*(c.c1+t*(2*c.c2+t*(3*c.c3+t*(4*c.c4+t*5*c.c5))));
    }
};

// Each shift is shared by all j in a group for a single i.
struct Group {
    Vec shift;
    int begin, end;
};
struct Cell {
    int n;
    Matrix lattice{}, inverse{};
    Potential potential;
    static constexpr double skin = .22;
    std::vector<Vec> reference;
    std::vector<int> neighbors, first;
    std::vector<Group> groups;
    double self_energy = 0.;
    Matrix self_virial{};

    Cell(int count, const double* a, const double* inv, double rc, const double* p) :
        n(count), potential(rc,p), first(n+1) {
        for (int k=0; k<4; ++k) {
            std::copy_n(a+4*k,4,lattice[k].begin());
            std::copy_n(inv+4*k,4,inverse[k].begin());
        }
        // Opposite periodic images of each site contribute a single bond.
        each_image(Vec{}, rc, [&](const Tile& tile, const Vec& shift, double r2) {
            if (tile <= Tile{}) return;
            double e = 0.;
            const double g = potential.force_prefactor<true>(r2,e);
            self_energy += n*e;
            for (int k=0;k<4;++k) for (int l=0;l<4;++l)
                self_virial[k][l] += n*g*shift[k]*shift[l];
        });
    }

    template<class F> void each_image(const Vec& raw, double reach, F visit) const {
        Tile lo, hi;
        for (int k=0;k<4;++k) {
            double center = 0., norm2 = 0.;
            for (int l=0;l<4;++l) {
                center -= raw[l]*inverse[l][k];
                norm2 += inverse[l][k]*inverse[l][k];
            }
            // If |raw+tile*A|<reach, every fractional component lies in
            // +/-reach times the corresponding reciprocal column norm.
            const double width = reach*std::sqrt(norm2)+1e-12;
            lo[k] = int(std::ceil(center-width));
            hi[k] = int(std::floor(center+width));
        }
        for (int a=lo[0];a<=hi[0];++a) for (int b=lo[1];b<=hi[1];++b)
        for (int c=lo[2];c<=hi[2];++c) for (int d=lo[3];d<=hi[3];++d) {
            Tile tile{a,b,c,d};
            Vec shift;
            double r2 = 0.;
            for (int k=0;k<4;++k) {
                shift[k] = a*lattice[0][k]+b*lattice[1][k]+c*lattice[2][k]+d*lattice[3][k];
                const double v = raw[k]+shift[k];
                r2 += v*v;
            }
            if (r2<reach*reach) visit(tile,shift,r2);
        }
    }

    void rebuild_if_needed(const double* x) {
        if (!reference.empty()) {
            bool rebuild = false;
            for (int i=0;i<n&&!rebuild;++i) {
                double travel2 = 0.;
                for (int k=0;k<4;++k) {
                    const double dx = x[4*i+k]-reference[i][k];
                    travel2 += dx*dx;
                }
                rebuild = travel2 >= .249*skin*skin;
            }
            if (!rebuild) return;
        }
        reference.resize(n);
        neighbors.clear();
        groups.clear();
        for (int i=0;i<n;++i) {
            std::copy_n(x+4*i,4,reference[i].begin());
            first[i] = groups.size();
            std::map<Tile, std::pair<Vec,std::vector<int>>> images;
            for (int j=i+1;j<n;++j) {
                Vec raw;
                for (int k=0;k<4;++k) raw[k]=x[4*j+k]-x[4*i+k];
                each_image(raw,potential.cutoff+skin,[&](const Tile& tile, const Vec& shift, double) {
                    auto& group = images[tile];
                    group.first = shift;
                    group.second.push_back(j);
                });
            }
            for (const auto& entry:images) {
                const auto& [shift,sites] = entry.second;
                const int begin = neighbors.size();
                neighbors.insert(neighbors.end(),sites.begin(),sites.end());
                groups.push_back({shift,begin,int(neighbors.size())});
            }
        }
        first[n] = groups.size();
    }

    template<bool read>
    double evaluate(const double* x, double* grad, double* virial) {
        rebuild_if_needed(x);
        double energy = self_energy;
        const double rc2 = potential.cutoff*potential.cutoff;
        std::fill_n(grad,4*n,0.);
        if constexpr (read) {
            for (int k=0;k<4;++k) std::copy_n(self_virial[k].begin(),4,virial+4*k);
        }
        for (int i=0;i<n;++i) {
            Vec gi{};
            for (int h=first[i];h<first[i+1];++h) {
                const auto& group = groups[h];
                Vec origin;
                for (int k=0;k<4;++k) origin[k] = group.shift[k]-x[4*i+k];
                for (int b=group.begin;b<group.end;++b) {
                    const int j = neighbors[b];
                    Vec v;
                    double r2 = 0.;
                    for (int k=0;k<4;++k) {v[k] = x[4*j+k]+origin[k]; r2 += v[k]*v[k];}
                    if (r2 >= rc2) continue;
                    if (r2 == 0.) return std::numeric_limits<double>::quiet_NaN();
                    double e;
                    const double g = potential.force_prefactor<read>(r2,e);
                    if constexpr (read) energy += e;
                    for (int k=0;k<4;++k) {
                        const double f = g*v[k];
                        grad[4*j+k] += f;
                        gi[k] -= f;
                        if constexpr (read) for (int l=0;l<4;++l) virial[4*k+l] += f*v[l];
                    }
                }
            }
            for (int k=0;k<4;++k) grad[4*i+k] += gi[k];
        }
        return energy;
    }
};
}  // namespace

extern "C" {
void* thermal_make(int n,const double* a,const double* inv,double rc,const double* p) {
    try {return new Cell(n,a,inv,rc,p);} catch (...) {return nullptr;}
}
void thermal_free(void* ptr) {delete static_cast<Cell*>(ptr);}
double thermal_read(void* ptr,const double* x,double* grad,double* virial) {
    try {return static_cast<Cell*>(ptr)->evaluate<true>(x,grad,virial);}
    catch (...) {return std::numeric_limits<double>::quiet_NaN();}
}
// Negative means allocation failure. A positive result is the failed step,
// one-based. Positions contain only the preceding complete steps.
int thermal_advance(void* ptr,int steps,double dt,const double* noise,double* x) {
    try {
        Cell& cell = *static_cast<Cell*>(ptr);
        const size_t count = 4*cell.n;
        std::vector<double> grad(count);
        for (int step=0;step<steps;++step) {
            if (!std::isfinite(cell.evaluate<false>(x,grad.data(),nullptr))) return step+1;
            for (size_t k=0;k<count;++k) {
                grad[k] = noise[count*step+k]-dt*grad[k];
                if (!std::isfinite(grad[k]) || std::abs(grad[k])>1.) return step+1;
            }
            for (size_t k=0;k<count;++k) x[k] += grad[k];
        }
        return 0;
    } catch (...) {return -1;}
}
}
