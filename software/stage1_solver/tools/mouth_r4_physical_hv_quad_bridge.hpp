// Audit bridge. Inputs are exact binary64 values widened in C++. Both displaced
// gradients and their subtraction are computed in binary128 before any narrowing.
struct WireGrid { int n; const double *v; };
struct WireData {
    WireGrid x,t,a,depth_edges,angle_edges,base_x,base_t,field_x,field_t,field_a,cut_x,cut_t,cut_a,reg_e,reg_t,gauss;
    const double *weights,*labels,*base,*correction,*offset,*spline,*shape;
    int spline_n,diagnose;
    double length,ref_length,pressure,core_ratio,epsilon;
    int exact_melt_support;
    const double *collar_coeff;
};

extern "C" int r4_quad_derivative(const WireData *source, const double *labels,
        const double *order, const double *direction, double h, WireGrid gp,
        const double *weights, const int *selected, double *product) {
    std::vector<Real> tables[12];
    Data d{};
    int bank=0;
    auto widen=[&](WireGrid input) {
        auto &stored=tables[bank++];
        stored.assign(input.v,input.v+input.n);
        return Grid{input.n,stored.data()};
    };
    d.x=widen(source->x); d.t=widen(source->t); d.a=widen(source->a);
    d.field_x=widen(source->field_x); d.field_t=widen(source->field_t);
    d.field_a=widen(source->field_a); d.reg_e=widen(source->reg_e); d.reg_t=widen(source->reg_t);
    d.labels=widen({2*d.reg_e.n*d.reg_t.n,source->labels}).v;
    d.spline_n=source->spline_n; d.spline=widen({4*d.spline_n,source->spline}).v;
    d.shape=widen({8,source->shape}).v;
    if (!source->collar_coeff || !finiteq(h) || h<=0) return 10;
    const int collar_cells=(d.reg_e.n-1)/2*(d.reg_t.n-1);
    d.collar_coeff=widen({32*collar_cells,source->collar_coeff}).v;
    d.length=source->length; d.ref_length=source->ref_length;
    d.pressure=source->pressure; d.core_ratio=source->core_ratio; d.epsilon=source->epsilon;

    const int n=d.field_x.n*d.field_a.n*d.field_t.n;
    std::vector<Real> ys(3*n), chi(order,order+n), gradients[2], load(n), intact(4*d.field_t.n);
    std::vector<Real> points(gp.v,gp.v+gp.n), wt(weights,weights+gp.n);
    Real diagnostics[23], energy=0;
    for (int arm=0;arm<2;++arm) {
        Real shift=(arm?1:-1)*Real(h);
        for (int i=0;i<3*n;++i) ys[i]=Real(labels[i])+shift*Real(direction[i]);
        gradients[arm].resize(3*n);
        std::fill(load.begin(),load.end(),Real(0));
        std::fill(intact.begin(),intact.end(),Real(-1));
        int rc=r4_eulerian_material(&d,ys.data(),chi.data(),{gp.n,points.data()},wt.data(),
                      &energy,gradients[arm].data(),load.data(),diagnostics,intact.data(),selected);
        if (rc) return rc;
    }
    for (int i=0;i<3*n;++i)
        product[i]=static_cast<double>((gradients[1][i]-gradients[0][i])/(2*Real(h)));
    return 0;
}
