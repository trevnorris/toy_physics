// Included in the promoted signed-energy translation unit, sharing its radial
// cuts and quadrature. Compare actual branch IDs, without comparing coordinates
// of cuts that legitimately move inside the same interval.
namespace endpoint_depth {
using Sample = std::array<int, 3>;

void branches(const Data &d, const Q1 &q, Grid rule, int it, int ia,
              std::vector<Sample> &result) {
    const auto left = q.at(0,rule.v[it],rule.v[ia]).xi;
    const auto right = q.at(1,rule.v[it],rule.v[ia]).xi;
    const Real weld = d.x.v[(d.x.n-1)*2/3];
    euler::Cuts cuts;
    euler::radial_cuts(cuts,d,left,right,weld);
    result.clear();
    for (size_t i=1; i<cuts.size(); ++i) {
        for (int k=0; k<rule.n; ++k) {
            const Real u = cuts[i-1].u + rule.v[k]*(cuts[i].u-cuts[i-1].u);
            Vec p;
            for (int a=0; a<3; ++a) p[a]=left[a]+u*(right[a]-left[a]);
            Sample id={-1,-1,-1};
            if (p[0]>weld) {
                const Real eta=-std::log(p[0])/d.pressure;
                const int c=cell(d.reg_e,eta), z=cell(d.reg_t,p[1]);
                const Real mu=(eta-d.reg_e.v[c])/(d.reg_e.v[c+1]-d.reg_e.v[c]);
                const Real nu=(p[1]-d.reg_t.v[z])/(d.reg_t.v[z+1]-d.reg_t.v[z]);
                id={c,z,int(mu+nu>1)};
            }
            // Adjacent integration cuts can sample the very same reference
            // triangle. Only changing the stored reference metric matters here.
            if (result.empty() || result.back()!=id) result.push_back(id);
        }
    }
}
}

extern "C" int endpoint_changed_rays(const Data *data, const Real *labels, Grid rule,
                                     const int *active, const Real *dt, int *first) {
    const auto &d=*data;
    const int nx=d.field_x.n, nt=d.field_t.n, na=d.field_a.n;
    int changed=0;
    std::vector<endpoint_depth::Sample> at_zero, at_shift;
    for (int i=0; i<nx-1; ++i) for (int j=0; j<nt-1; ++j) for (int k=0; k<na-1; ++k) {
        if (!active[(i*(nt-1)+j)*(na-1)+k]) continue;
        Q1 base{}, moved{};
        for (int x=0; x<2; ++x) for (int t=0; t<2; ++t) for (int a=0; a<2; ++a) {
            const int n=((i+x)*nt+j+t)*na+k+a;
            for (int c=0; c<3; ++c) {
                base.labels[x][t][a][c]=labels[3*n+c];
                moved.labels[x][t][a][c]=labels[3*n+c]+(c==1 ? dt[j+t] : 0);
            }
        }
        for (int it=0; it<rule.n; ++it) for (int ia=0; ia<rule.n; ++ia) {
            endpoint_depth::branches(d,base,rule,it,ia,at_zero);
            endpoint_depth::branches(d,moved,rule,it,ia,at_shift);
            if (at_zero!=at_shift) {
                if (!changed) {
                    const int site[]={i,j,k,it,ia};
                    std::copy(site,site+5,first);
                }
                ++changed;
            }
        }
    }
    return changed;
}
