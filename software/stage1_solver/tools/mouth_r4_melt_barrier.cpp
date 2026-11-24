// Positive cell-volume corner quadrature for the dimensionless chart-volume barrier.
// Compute the chart determinant from contiguous vertex differences to avoid
// cancellation of the untranslated material coordinates at the limiting cell.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

using Edge = std::array<long double,3>;

static Edge cross3(const Edge& a,const Edge& b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
extern "C" int melt_barrier(int count,const int* ids,const double* width,
        const double* weight,const double* labels,double* gradient,double* determinants,
        long double* total) {
    *total=0;
    for(int c=0;c<count;++c){
        for(int v=0;v<8;++v){
            int node=ids[8*c+v];
            std::array<Edge,3> column,cofactor;
            for(int d=0;d<3;++d){
                int bit=4>>d; int other=ids[8*c+(v^bit)];
                long double sign=(v&bit)?-1.:1.;
                for(int a=0;a<3;++a){
                    column[d][a]=sign*(static_cast<long double>(labels[3*other+a])-labels[3*node+a])/width[3*c+d];
                }
            }
            for(int d=0;d<3;++d) cofactor[d]=cross3(column[(d+1)%3],column[(d+2)%3]);
            long double j=0;
            for(int a=0;a<3;++a) j+=column[0][a]*cofactor[0][a];
            if(!(j>0) || !std::isfinite(j)) return 1;
            determinants[c*8+v]=static_cast<double>(j);
            const long double z=j-1;
            // Use log(j) when the determinant is small. Subtracting one and
            // restoring it inside log1p would discard its significant digits.
            *total+=weight[c]*(z-(j<.5L ? std::log(j):std::log1p(z)));
            for(int d=0;d<3;++d){
                int bit=4>>d; int other=ids[8*c+(v^bit)];
                long double a=weight[c]*(1.-1./j)*((v&bit)?-1.:1.)/width[3*c+d];
                for(int k=0;k<3;++k){
                    double f=static_cast<double>(a*cofactor[d][k]);
                    gradient[3*other+k]+=f; gradient[3*node+k]-=f;
                }
            }
        }
    }
    return 0;
}
