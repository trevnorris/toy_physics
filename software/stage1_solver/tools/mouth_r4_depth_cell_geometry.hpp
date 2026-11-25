// Read-only, promoted-precision trace of an existing physical cell-work quadrature.
// Included after the unmodified promoted crystal kernel. One row per production sample.
extern "C" int depth_cell_trace(const Data *data,const Real *labels,const Real *order,Grid rule,
        const Real *weights,const int *index,Real *rows,int capacity,Real *center){
    const Data &d=*data;
    const Grid axes[]={d.field_x,d.field_t,d.field_a};
    Geometry g{d,{},{}};
    for(int i=0;i<3;++i){
        if(index[i]<0 || index[i]+1>=axes[i].n)return -2;
        g.lo[i]=axes[i].v[index[i]];
        g.width[i]=axes[i].v[index[i]+1]-g.lo[i];
    }
    Vec pos{};
    for(int i=0;i<3;++i)pos[i]=g.lo[i]+g.width[i]/2;
    center[0]=radius(d,pos).R;center[1]=d.length*pos[1];center[2]=pos[2];
    Q1 q{};
    for(int x=0;x<2;++x)for(int t=0;t<2;++t)for(int a=0;a<2;++a){
        int id=((index[0]+x)*axes[1].n+index[1]+t)*axes[2].n+index[2]+a;
        q.order[x][t][a]=order[id];
        for(int z=0;z<3;++z)q.labels[x][t][a][z]=labels[3*id+z];
    }
    const Real weld=d.x.v[(d.x.n-1)*2/3];
    int count=0;
    for(int t=0;t<rule.n;++t)for(int a=0;a<rule.n;++a){
        const Real v=rule.v[t],theta=rule.v[a];
        const auto left=q.at(0,v,theta),right=q.at(1,v,theta);
        euler::Cuts cuts;
        euler::radial_cuts(cuts,d,left.xi,right.xi,weld);
        for(size_t s=1;s<cuts.size();++s){
            Real width=cuts[s].u-cuts[s-1].u;
            for(int ix=0;ix<rule.n;++ix){
                if(count==capacity)return -1;
                Real *row=rows+56*count++;
                const Real u=cuts[s-1].u+width*rule.v[ix];
                Vec local{u,v,theta};
                auto p=left;
                p.chi+=u*(right.chi-left.chi);
                for(int c=0;c<3;++c){
                    p.xi[c]+=u*(right.xi[c]-left.xi[c]);
                    for(int k=0;k<3;++k)p.D[c][k]+=u*(right.D[c][k]-left.D[c][k]);
                    pos[c]=g.lo[c]+local[c]*g.width[c];
                }
                const auto state=euler::material(g,local,p,false);
                const auto wall=radius(d,pos);
                const bool core=p.xi[0]<=weld;
                const Ref ref=reference(d,p.xi[0],p.xi[1],profile(d,p.xi[1],true),core);
                const Mat invF=product(invert(ref.inv),product(p.D,invert(physical_basis(wall,d.length,g.width))));
                row[0]=t;row[1]=a;row[2]=s-1;row[3]=ix;
                row[4]=u;row[5]=v;row[6]=theta;
                row[7]=cuts[s-1].u;row[8]=cuts[s].u;
                row[9]=wall.R;row[10]=d.length*pos[1];row[11]=pos[2];
                row[12]=p.chi;row[13]=state.Jinv;row[14]=state.W;
                row[15]=weights[t]*weights[a]*weights[ix]*width*p.chi*state.unit;
                row[16]=core; row[17]=std::clamp(int(pos[1]*d.spline_n),0,d.spline_n-1);
                row[18]=core?std::clamp(int(p.xi[1]*d.spline_n),0,d.spline_n-1):-1;
                row[19]=pos[1]>.06 && pos[1]<.62;
                row[20]=row[21]=row[22]=-1;
                if(!core){
                    const Real eta=-std::log(p.xi[0])/d.pressure;
                    const int i=cell(d.reg_e,eta),j=cell(d.reg_t,p.xi[1]);
                    const Real uu=(eta-d.reg_e.v[i])/(d.reg_e.v[i+1]-d.reg_e.v[i]);
                    const Real vv=(p.xi[1]-d.reg_t.v[j])/(d.reg_t.v[j+1]-d.reg_t.v[j]);
                    row[20]=i;row[21]=j;row[22]=uu+vv>1.;
                }
                for(int k=0;k<3;++k){
                    row[23+k]=p.xi[k];row[26+k]=left.xi[k];row[29+k]=right.xi[k];
                    row[32+k]=cuts[s-1].normal[k];row[35+k]=cuts[s].normal[k];
                    for(int j=0;j<3;++j)row[38+3*k+j]=invF[k][j];
                }
                row[47]=ref.S*std::sin(p.xi[2])/(wall.R*std::sin(pos[2]));
                row[48]=state.unit*p.chi;row[49]=p.D[0][0];
                row[50]=pos[0];row[51]=pos[1];row[52]=p.xi[0]-weld;
                row[53]=cell(d.x,p.xi[0]);row[54]=cell(d.t,p.xi[1]);row[55]=cell(d.a,p.xi[2]);
            }
        }
    }
    return count;
}
