// Widen every frozen grid entry exactly from the production binary64/80 ABI.
// Subtract the reference cell energy in binary128, before conversion for reporting.
struct InputGrid {int n;const long double *v;};
struct InputData {
    InputGrid x,t,a,depth_edges,angle_edges,base_x,base_t,field_x,field_t,field_a,cut_x,cut_t,cut_a,reg_e,reg_t,gauss;
    const long double *weights,*labels,*base,*correction,*offset,*spline,*shape;
    int spline_n,diagnose;long double length,ref_length,pressure,core_ratio,epsilon;int exact_melt_support;const long double *collar_coeff;
};
extern "C" int high_precision_banks(const InputData *h,const long double *ys,const long double *chi,InputGrid rule,
    const long double *wt,const int *enabled,long double *out,long double *lower){
    static std::vector<Real> store[12],baseline;
    static Data d{};
    if(!d.x.v){
        // Only inputs referenced by this held geometry's physical-cell energy.
        auto copy=[&](int i,InputGrid v){store[i].assign(v.v,v.v+v.n);return Grid{v.n,store[i].data()};};
        d.x=copy(0,h->x);d.t=copy(1,h->t);d.a=copy(2,h->a);
        d.field_x=copy(3,h->field_x);d.field_t=copy(4,h->field_t);d.field_a=copy(5,h->field_a);
        d.reg_e=copy(6,h->reg_e);d.reg_t=copy(7,h->reg_t);
        d.labels=copy(8,{2*d.reg_e.n*d.reg_t.n,h->labels}).v;
        d.spline_n=h->spline_n;d.spline=copy(9,{4*h->spline_n,h->spline}).v;
        d.shape=copy(10,{8,h->shape}).v;
        if(h->collar_coeff){
            const int cells=(d.reg_e.n-1)/2*(d.reg_t.n-1);
            d.collar_coeff=copy(11,{32*cells,h->collar_coeff}).v;
        }
        d.length=h->length;d.ref_length=h->ref_length;d.epsilon=h->epsilon;
        d.pressure=h->pressure;d.core_ratio=h->core_ratio;
    }
    const int n=d.field_x.n*d.field_t.n*d.field_a.n;
    const int m=(d.field_x.n-1)*(d.field_t.n-1)*(d.field_a.n-1);
    std::vector<Real> labels(ys,ys+3*n),order(chi,chi+n),qp(rule.v,rule.v+rule.n),weights(wt,wt+rule.n),values(m);
    Real min_det;
    int rc=r4_cell_work(&d,labels.data(),order.data(),{rule.n,qp.data()},weights.data(),enabled,values.data(),&min_det);
    if(rc)return rc;
    *lower=static_cast<long double>(min_det);
    if(baseline.empty())baseline=values;
    for(int i=0;i<m;i++)out[i]=static_cast<long double>(values[i]-baseline[i]);
    return 0;
}
