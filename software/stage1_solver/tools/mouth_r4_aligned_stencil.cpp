// Q1 chart determinant quadrature at 27 Lobatto nodes and 12 moving edge nodes.
// All four evaluation modes instantiate the same cell energy and its gradient.
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

using Real = long double;
using Complex = std::complex<Real>;

// Differential<false> carries a first variation. Differential<true> carries the
// finite increment exactly through the algebra before rounding the two states.
template<bool finite> struct Differential {
    Real base, change;
    Differential(Real value=0, Real delta=0) : base(value), change(delta) {}
};
template<bool f> using D = Differential<f>;
template<bool f> D<f> operator+(D<f> x,D<f> y) { return {x.base+y.base,x.change+y.change}; }
template<bool f> D<f> operator-(D<f> x,D<f> y) { return {x.base-y.base,x.change-y.change}; }
template<bool f> D<f> operator*(D<f> x,D<f> y) {
    return {x.base*y.base,x.base*y.change+x.change*y.base+(f?x.change*y.change:0)};
}
template<bool f> D<f> operator/(D<f> x,D<f> y) {
    const Real q=x.base/y.base;
    return {q,(x.change-q*y.change)/(y.base+(f?y.change:0))};
}
using std::sqrt;
using std::log;
template<bool f> D<f> sqrt(D<f> x) {
    const Real r=sqrt(x.base);
    return {r,x.change/(r+(f?sqrt(x.base+x.change):r))};
}
template<bool f> D<f> log(D<f> x) {
    const Real q=x.change/x.base;
    return {log(x.base),f?std::log1p(q):q};
}
bool positive(Real x) { return x>0 && std::isfinite(x); }
bool positive(Complex x) { return positive(x.real()) && std::isfinite(x.imag()); }
template<bool f> bool positive(D<f> x) {
    return positive(x.base) && std::isfinite(x.change) && (!f || positive(x.base+x.change));
}

constexpr int index(int i,int j,int k) { return 9*i+3*j+k; }
struct Rule {
    Real basis[27][8][3]{};
    Real weights[27]{};
    std::array<std::array<int,3>,12> edges{};
    Rule() {
        int e=0;
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) for(int k=0;k<3;++k) {
            const int q=index(i,j,k);
            const std::array<int,3> ijk={i,j,k};
            int mids=(i==1)+(j==1)+(k==1);
            // Split each edge midpoint's mass equally with its moving sample.
            weights[q]=(mids==0?1:mids==1?2:mids==2?16:64)/216.L;
            if(mids==1) {
                int stride=i==1?9:j==1?3:1;
                edges[e++]={q-stride,q,q+stride};
            }
            for(int v=0;v<8;++v) {
                int bits[3]={(v>>2)&1,(v>>1)&1,v&1};
                for(int c=0;c<3;++c) {
                    Real b=bits[c]?1:-1;
                    for(int h=0;h<3;++h) if(h!=c)
                        b*=bits[h] ? .5L*ijk[h] : 1-.5L*ijk[h];
                    basis[q][v][c]=b;
                }
            }
        }
    }
};
const Rule rule;

template<class T> using Vec = std::array<T,3>;
template<class T> Vec<T> cross(const Vec<T>& a,const Vec<T>& b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

template<class T,bool gradient>
int cell(const T* local,const T* width,T& energy,T* nodal) {
    T det[27],adjoint[27]{};
    std::array<Vec<T>,3> cof[27];
    for(int q=0;q<27;++q) {
        std::array<Vec<T>,3> columns{};
        for(int v=1;v<8;++v) for(int d=0;d<3;++d) {
            const Real b=rule.basis[q][v][d];
            if(b==0) continue;
            for(int c=0;c<3;++c)
                columns[d][c]=columns[d][c]+T(b)*(local[3*v+c]-local[c]);
        }
        for(int d=0;d<3;++d) for(int c=0;c<3;++c) columns[d][c]=columns[d][c]/width[d];
        cof[q][0]=cross(columns[1],columns[2]);
        if constexpr(gradient) {
            cof[q][1]=cross(columns[2],columns[0]);
            cof[q][2]=cross(columns[0],columns[1]);
        }
        T j{};
        for(int c=0;c<3;++c) j=j+columns[0][c]*cof[q][0][c];
        det[q]=j;
        if(!positive(j)) return q+1;
        energy=energy+T(rule.weights[q])*(j-T(1)-log(j));
        if constexpr(gradient) adjoint[q]=T(rule.weights[q])*(T(1)-T(1)/j);
    }
    for(int e=0;e<12;++e) {
        const auto [l,m,r]=rule.edges[e];
        const T a=det[l], b=det[m], c=det[r];
        const T sa=sqrt(a), sc=sqrt(c);
        const T t=sa/(sa+sc), u=sc/(sa+sc);
        const T n0=u*(T(1)-T(2)*t), nm=T(4)*t*u, n1=t*(T(2)*t-T(1));
        const T j=n0*a+nm*b+n1*c;
        if(!positive(j)) return 28+e;
        const T weight=T(rule.weights[m]);
        energy=energy+weight*(j-T(1)-log(j));
        if constexpr(gradient) {
            const T dJdt=(T(4)*t-T(3))*a+(T(4)-T(8)*t)*b+(T(4)*t-T(1))*c;
            const T shift=dJdt*t*u/T(2);
            const T multiplier=weight*(T(1)-T(1)/j);
            adjoint[l]=adjoint[l]+multiplier*(n0+shift/a);
            adjoint[m]=adjoint[m]+multiplier*nm;
            adjoint[r]=adjoint[r]+multiplier*(n1-shift/c);
        }
    }
    if constexpr(gradient) {
        for(int q=0;q<27;++q) for(int d=0;d<3;++d) {
            const T scale=adjoint[q]/width[d];
            for(int v=0;v<8;++v) {
                const Real b=rule.basis[q][v][d];
                if(b==0) continue;
                const T coefficient=T(b)*scale;
                for(int c=0;c<3;++c)
                    nodal[3*v+c]=nodal[3*v+c]+coefficient*cof[q][d][c];
            }
        }
    }
    return 0;
}

template<class T>
int values(int count,const T* local,const T* width,T* energy,T* grad,int* bad) {
    for(int c=0;c<count;++c) {
        const int rc=grad?cell<T,true>(local+24*c,width+3*c,energy[c],grad+24*c):
                          cell<T,false>(local+24*c,width+3*c,energy[c],nullptr);
        if(rc) { *bad=c; return rc; }
    }
    return 0;
}
extern "C" int real_values(int n,const Real* x,const Real* w,Real* e,Real* g,int* bad) {
    return values(n,x,w,e,g,bad);
}
extern "C" int complex_values(int n,const Complex* x,const Complex* w,Complex* e,Complex* g,int* bad) {
    return values(n,x,w,e,g,bad);
}

template<bool finite>
int increments(int count,const Real* from,const Real* change,const Real* width,
               Real* result,int* bad) {
    for(int c=0;c<count;++c) {
        D<finite> local[24],w[3],energy{},nodal[24]{};
        for(int k=0;k<24;++k) local[k]={from[c*24+k],change[c*24+k]};
        for(int d=0;d<3;++d) w[d]=width[c*3+d];
        const int rc=cell<D<finite>,!finite>(local,w,energy,nodal);
        if(rc) { *bad=c; return rc; }
        if constexpr(finite) result[c]=energy.change;
        else for(int k=0;k<24;++k) result[24*c+k]=nodal[k].change;
    }
    return 0;
}
extern "C" int signed_work(int n,const Real* x,const Real* dx,const Real* w,Real* delta,int* bad) {
    return increments<true>(n,x,dx,w,delta,bad);
}
extern "C" int hessian_action(int n,const Real* x,const Real* v,const Real* w,Real* hv,int* bad) {
    return increments<false>(n,x,v,w,hv,bad);
}
