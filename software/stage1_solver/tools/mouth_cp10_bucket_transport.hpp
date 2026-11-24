// A common physical-chart cell receives the pair of signed energies. The same
// checkpoint-10 cell and spatial mask apply to the native and projected maps.
struct CellWork {
    enum Partition { common, original, split } partition=common;
    const Data &chart;
    using Pair=std::array<long double,2>;
    std::vector<Pair> archived,exact_cut,spans_weld;
    Grid refined_x,refined_t,refined_a;
    const unsigned char *weld;

    CellWork(const Data &d,Grid x,Grid t,Grid a,const unsigned char *mask):
        chart(d),refined_x(x),refined_t(t),refined_a(a),weld(mask) {
        size_t n=(d.field_x.n-1)*(d.field_t.n-1)*(d.field_a.n-1);
        archived.resize(n);exact_cut.resize(n);spans_weld.resize(n);
    }
    size_t locate(const Vec &p)const {
        return (cell(chart.field_x,p[0])*(chart.field_t.n-1)+cell(chart.field_t,p[1]))
               *(chart.field_a.n-1)+cell(chart.field_a,p[2]);
    }
    bool straddles(const Vec &p) const {
        return weld[(cell(refined_x,p[0])*(refined_t.n-1)+cell(refined_t,p[1]))
                    *(refined_a.n-1)+cell(refined_a,p[2])];
    }
    void deposit(size_t id,bool straddles_weld,double wt,const double (&energy)[2][2],double h) {
        for(int map=0;map<2;map++){
            double df=wt*((energy[0][map]-energy[1][map])/(2*h));
            if(partition!=split){
                archived[id][map]+=df;
                if(straddles_weld)spans_weld[id][map]+=df;
            }
            if(partition!=original)exact_cut[id][map]+=df;
        }
    }
    void export_cells(double *out) const {
        for(size_t i=0;i<archived.size();i++)for(int j=0;j<2;j++){
            out[6*i+j]=archived[i][j];out[6*i+2+j]=exact_cut[i][j];
            out[6*i+4+j]=spans_weld[i][j];
        }
    }
};
