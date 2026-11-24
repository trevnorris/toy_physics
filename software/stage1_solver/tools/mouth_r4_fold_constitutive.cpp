// Point readback on a specified saved material prism. The included file is
// selected from the source pin of the producing solver.
#include <mouth_free_shape_solid.cpp>

extern "C" void fold_point(const Data *data, int i, int j, int k, int upper,
                           double x, double t, double a, double *out) {
    const Data &d = *data;
    const double lo = d.x.v[i], dx = d.x.v[i+1]-lo;
    const double v = (t-d.t.v[j])/(d.t.v[j+1]-d.t.v[j]);
    // On a shared diagonal use the requested prism's polynomial extension.
    double anchor = x;
    if ((((x-lo)/dx+v)>1.) != bool(upper))
        anchor = lo+dx*(1.-v+(upper ? .125 : -.125));
    const TriSlice s(d,i,j,k,t,a,anchor);
    Vec X = s.X;
    Mat M = s.M;
    for (int c=0; c<3; ++c) {
        X[c] += s.M[c][0]*(x-anchor);
        M[c][2] += s.Max[c]*(x-anchor);
    }
    const auto ref = reference(d,x,t,profile(d,t,true),i<(d.x.n-1)*2/3);
    const auto wall = radius(d,X);
    const Mat C = {wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}};
    const Mat F = product(C,product(M,ref.inv));
    // Analytic polar limit sin(theta')/sin(theta) when the vertex lies
    // on the symmetry axis. The material map fixes that axis.
    const double angle_ratio = a==0. ? M[2][2] : std::sin(X[2])/std::sin(a);
    const double hoop = wall.R*angle_ratio/ref.S;
    const double scale = wall.grad[0]*d.length*wall.R*hoop/(ref.det2*ref.S);
    double J = scale*determinant(M);
    double W = J>0. ? .5*(norm2(F)+hoop*hoop-4.)-std::log(J)+.25*(J-1.)*(J-1.) : NAN;
    int off=0;
    for (double z: X) out[off++]=z;
    for (auto matrix: {M,C,ref.inv,F}) for(auto row: matrix) for(double z: row) out[off++]=z;
    for (double z: {hoop,scale,determinant(M),J,W,ref.S}) out[off++]=z;
    if (x>0. && x<1. && t>0. && t<1. && a>0. && a<pi/2.) {
        std::array<double,30> stats{};
        Results results{nullptr,nullptr,stats.data()};
        double weighted=0., dx_density=0.;
        int err = sample(d,results,profile(d,t,true),s,t,a,x,1.,i<(d.x.n-1)*2/3,weighted,dx_density);
        out[off++]=err;
        out[off++]=err ? NAN : weighted/(4*pi*ref.det2*ref.S*ref.S*std::sin(a));
    } else {
        out[off++]=-1.;
        out[off++]=NAN;
    }
}
