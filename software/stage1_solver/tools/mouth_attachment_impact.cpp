// Dilute tag scattering on the cold fixed-area triangular Surface. Each input
// row starts with a new undisturbed sheet. Time, energy and mass are native.
// Build with g++ -O3 -std=c++17. Input: id u v z vx vy vz chi [ax ay gate].
// u,v multiply the primitive vectors; all velocities are Cartesian.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double sin60 = .86602540378443864676;
constexpr double pi = 3.14159265358979323846;
constexpr double switch_start = 4.5, cutoff = 5.;
constexpr int host_shell = 5;  // Reference-lattice units, exactly as in Bilayer.
struct V2 { double x=0, y=0; };
struct V3 { double x=0, y=0, z=0; };
double norm2(V3 v) { return v.x*v.x + v.y*v.y + v.z*v.z; }
int wrap(long long i, int n) { int r = i % n; return r < 0 ? r+n : r; }

struct Pair { double energy, radial; };  // radial=(dV/dr)/r
struct Cubic { double a, b, c, d; };

// Tabulate each *potential* as a cubic Hermite polynomial in q=r^2. Use
// its own exact polynomial derivative for the force. --direct 1 evaluates
// both from Morse exponentials instead. Extreme separations always use Morse.
class Morse {
    static constexpr int bins = 2048;
    static constexpr double low = .25, high = 32.;
    double depth, alpha, rest;
    std::vector<Cubic> table;
    Pair exact(double q) const {
        if (!(q > 0.) || !std::isfinite(q)) throw std::runtime_error("invalid pair separation");
        double r = std::sqrt(q), w = std::exp(-alpha*(r-rest));
        return {depth*w*(w-2.), 2.*depth*alpha*w*(1.-w)/r};
    }
public:
    Morse(double depth_, double alpha_, double rest_, bool direct)
        : depth(depth_), alpha(alpha_), rest(rest_) {
        if (direct) return;
        int n = static_cast<int>((high-low)*bins);
        table.reserve(n);
        Pair left = exact(low);
        for (int i=0; i<n; ++i) {
            Pair right = exact(low+double(i+1)/bins);
            double b = left.radial/(2*bins), e = right.radial/(2*bins);
            double delta = right.energy-left.energy;
            table.push_back({left.energy, b, 3*delta-2*b-e, b+e-2*delta});
            left=right;
        }
    }
    Pair operator()(double q) const {
        if (table.empty() || !(q >= low && q < high)) return exact(q);
        double index = (q-low)*bins;
        int i = static_cast<int>(index);
        double t = index-i;
        const Cubic &p = table[i];
        return {p.a+t*(p.b+t*(p.c+t*p.d)), 2*bins*(p.b+t*(2*p.c+3*t*p.d))};
    }
};

struct Params {
    int side=0, steps=0;
    double pitch=0, dt=0, gate_slope=0;
    bool mobile=true, mixed=true, direct=false;
    std::string input="-";
};
bool flag(const std::string &s) {
    if (s=="1" || s=="0") return s=="1";
    throw std::runtime_error("boolean must be 0 or 1");
}
Params parse(int argc, char **argv) {
    Params p;
    for (int i=1; i<argc; i+=2) {
        if (i+1==argc) throw std::runtime_error("missing option value");
        std::string k=argv[i], v=argv[i+1];
        if (k=="--side") p.side=std::stoi(v);
        else if (k=="--dt") p.dt=std::stod(v);
        else if (k=="--steps") p.steps=std::stoi(v);
        else if (k=="--pitch") p.pitch=std::stod(v);
        else if (k=="--gate-slope") p.gate_slope=std::stod(v);
        else if (k=="--host" && (v=="mobile" || v=="rigid")) p.mobile=(v=="mobile");
        else if (k=="--mixed") p.mixed=flag(v);
        else if (k=="--direct") p.direct=flag(v);
        else if (k=="--input") p.input=v;
        else throw std::runtime_error("invalid option/value: "+k+" "+v);
    }
    if (p.side<1 || p.side>1000 || !(p.pitch>0.) || !(p.dt>0.) || p.steps<=0 ||
        !std::isfinite(p.pitch) || !std::isfinite(p.dt*p.steps) || !std::isfinite(p.gate_slope))
        throw std::runtime_error("supply valid --side, --pitch, --dt and --steps");
    return p;
}

struct Input { long long id; double u, v, z; V3 vel; double chi, ax, ay, gate; };
struct Shift {
    V2 ref;
    double prepared_energy;
    std::vector<int> neighbor;
};
class Surface {
    Params p;
    Morse within, between;
    std::vector<Shift> bonds, images;
    double image_allowance=-1;
    const int n;
    Shift shift(int i, int j, double prepared) const {
        Shift s{{p.pitch*(i+.5*j), p.pitch*sin60*j}, prepared, {}};
        s.neighbor.reserve(n);
        for (int a=0;a<p.side;++a)
            for (int b=0;b<p.side;++b)
                s.neighbor.push_back(wrap(a+i,p.side)*p.side+wrap(b+j,p.side));
        return s;
    }
    void cover_images(double displacement) {
        if (displacement<=image_allowance) return;
        // Every undeformed origin within this disk is represented, even if
        // several images share one mobile site. The tag is within sin60*pitch
        // of its folded primitive-cell centre. Rebuild after large host motion.
        image_allowance=displacement+.5*p.pitch;
        double span=cutoff+image_allowance+sin60*p.pitch;
        int limit=static_cast<int>(std::ceil(span/(p.pitch*sin60)));
        images.clear();
        for (int a=-limit;a<=limit;++a) for (int b=-limit;b<=limit;++b) {
            double x=p.pitch*(a+.5*b), y=p.pitch*sin60*b;
            if (x*x+y*y<=span*span) images.push_back(shift(a,b,0.));
        }
    }
public:
    explicit Surface(const Params &params)
        : p(params), within(1.,7.,1.,p.direct), between(.6,10.,.94,p.direct), n(p.side*p.side) {
        int limit=static_cast<int>(std::ceil(host_shell/sin60));
        for (int a=0;a<=limit;++a) for (int b=-limit;b<=limit;++b) {
            if ((a==0 && b<=0) || a*a+a*b+b*b > host_shell*host_shell) continue;
            double x=p.pitch*(a+.5*b), y=p.pitch*sin60*b;
            bonds.push_back(shift(a,b,within(x*x+y*y).energy));
        }
        cover_images(0.);
    }
    int sites() const { return n; }
    size_t directions() const { return bonds.size(); }

    struct Forces {
        explicit Forces(int n) : host(n), contact(n) {}
        std::vector<V2> host, contact;
        V3 tag;
        double strain, interaction;
    };
    void evaluate(V3 pos, const std::vector<V2> &u, double max_disp, double chi, Forces &f) {
        f.tag={}; f.strain=0.; f.interaction=0.;
        std::fill(f.host.begin(),f.host.end(),V2{});
        std::fill(f.contact.begin(),f.contact.end(),V2{});
        if (p.mobile) for (const Shift &s : bonds) {
            double strain=0.;
            for (int i=0;i<n;++i) {
                int j=s.neighbor[i];
                double x=s.ref.x+(u[j].x-u[i].x), y=s.ref.y+(u[j].y-u[i].y);
                Pair pair=within(x*x+y*y);
                strain+=pair.energy-s.prepared_energy;
                double fx=pair.radial*x, fy=pair.radial*y;
                f.host[i].x+=fx; f.host[i].y+=fy;
                f.host[j].x-=fx; f.host[j].y-=fy;
            }
            f.strain+=strain;
        }
        if (chi>0. && std::abs(pos.z)<cutoff) {
            cover_images(max_disp);
            double b=pos.y/(p.pitch*sin60), a=pos.x/p.pitch-.5*b;
            long long ia=std::llround(a), ib=std::llround(b);
            double rx=p.pitch*((a-ia)+.5*(b-ib)), ry=p.pitch*sin60*(b-ib);
            int anchor=wrap(ia,p.side)*p.side+wrap(ib,p.side);
            for (const Shift &s : images) {
                int j=s.neighbor[anchor];
                double x=rx-s.ref.x-u[j].x, y=ry-s.ref.y-u[j].y;
                double q=x*x+y*y+pos.z*pos.z;
                if (q>=cutoff*cutoff) continue;
                Pair pair=between(q);
                if (q>switch_start*switch_start) {
                    // C2 smootherstep sends the already tiny far tail to zero.
                    double r=std::sqrt(q), t=(r-switch_start)/(cutoff-switch_start);
                    double a=1.-t, sw=a*a*a*(1.+3.*t+6.*t*t);
                    double ds=-30.*t*t*a*a/(cutoff-switch_start);
                    pair.radial=sw*pair.radial+(ds/r)*pair.energy;
                    pair.energy*=sw;
                }
                double k=chi*pair.radial, fx=k*x, fy=k*y;
                f.tag.x-=fx; f.tag.y-=fy; f.tag.z-=k*pos.z;
                f.contact[j].x+=fx; f.contact[j].y+=fy;
                f.interaction+=chi*pair.energy;
            }
        }
        if (p.mixed && chi<1.) {
            constexpr double d=.6, a=10., rest=.94;
            double z=std::abs(pos.z), w=std::exp(-a*(z-rest));
            double fac=(1.-chi)*2.*pi*d/(sin60*p.pitch*p.pitch);
            f.interaction+=fac*(w*w*(z/(2*a)+1./(4*a*a))-2.*w*(z/a+1./(a*a)));
            f.tag.z+=fac*pos.z*w*(w-2.);
        }
        if (p.mobile) for (int i=0;i<n;++i) {
            f.host[i].x+=f.contact[i].x; f.host[i].y+=f.contact[i].y;
        }
    }
};

// Residence in the exposed window [0,hi), using a linear crossing within
// each drift step. Track the latest completed or ongoing residence separately
// from the continuous dwell ending at the final observation.
struct Window {
    double hi, entry=-1., exit=-1., last_exit=-1., active=-1., last=0., total=0.;
    bool inside(double z) const { return z>=0. && z<hi; }
    Window(double cut, double z) : hi(cut) {
        if (inside(z)) entry=active=0.;
    }
    void update(double t, double dt, double from, double to) {
        if (from==to) { if (inside(from)) total+=dt; return; }
        double below=(0.-from)/(to-from), above=(hi-from)/(to-from);
        double lo=std::max(0.,std::min(below,above)), up=std::min(1.,std::max(below,above));
        if (up>lo) {
            if (active<0.) { active=t+dt*lo; if (entry<0.) entry=active; }
            total+=dt*(up-lo);
        }
        if (inside(to) && active<0.) { active=t+dt; if (entry<0.) entry=active; }
        if (!inside(to) && active>=0.) {
            last_exit=t+dt*std::clamp(up,0.,1.);
            if (exit<0.) exit=last_exit;
            last=last_exit-active;
            active=-1.;
        }
    }
    void emit(const char *suffix, double t, double z) const {
        double tail=active<0. ? 0. : t-active;
        std::cout<<",\"entry"<<suffix<<"\":"<<entry<<",\"exit"<<suffix<<"\":"<<exit
                 <<",\"last_exit"<<suffix<<"\":"<<last_exit<<",\"residence"<<suffix<<"\":"<<total
                 <<",\"last_dwell"<<suffix<<"\":"<<(active<0. ? last : tail)
                 <<",\"tail_dwell"<<suffix<<"\":"<<tail<<",\"inside"<<suffix<<"\":"<<inside(z);
    }
};

struct Interval {
    double lo=0., hi=1.;
    // Intersect with the affine half-space a+(b-a)*fraction<0. An entire
    // interval sitting on the lower surface z=0 belongs to the exposed face.
    void below(double a, double b, bool closed=false) {
        if (a==b) {
            if (closed ? a>0. : a>=0.) hi=lo;
        } else if (a<b) hi=std::min(hi,-a/(b-a));
        else lo=std::max(lo,-a/(b-a));
    }
};

struct BindingReadout {
    double joint_first=-1., residence=0., gated=0.;
    double forward=-1., backward=-1., first_exit=-1.;
    int side=0;
    bool half=false;
    void update(double t, double h, double z0, double z1, double e0, double e1,
                double x0, double x1, double gate, double midpoint) {
        if (gate>0.) {
            if (forward<0. && x0<gate && x1>=gate) forward=t+h*(gate-x0)/(x1-x0);
            if (backward<0. && x0>=0. && x1<0.) backward=t+h*x0/(x0-x1);
            if (first_exit<0. && (forward>=0. || backward>=0.)) {
                side=forward>=0. ? 1 : -1;
                first_exit=side==1 ? forward : backward;
            }
        }
        Interval joint;
        joint.below(-z0,-z1,true);
        joint.below(z0-1.5,z1-1.5);
        joint.below(e0,e1);
        if (joint.hi>joint.lo) {
            double start=t+h*joint.lo, end=t+h*joint.hi;
            if (joint_first<0.) joint_first=start;
            residence+=end-start;
            gated+=std::max(0.,(first_exit<0. ? end : std::min(first_exit,end))-start);
        }
        if (t<=midpoint && midpoint<=t+h) {
            double f=(midpoint-t)/h, z=z0+f*(z1-z0), e=e0+f*(e1-e0);
            half=z>=0. && z<1.5 && e<0.;
        }
    }
};

void trajectory(const Params &p, Surface &s, Input in) {
    const int n=s.sites();
    V3 pos{p.pitch*(in.u+.5*in.v),p.pitch*sin60*in.v,in.z}, vel=in.vel;
    V3 start=pos;
    std::vector<V2> u(n), v(n);
    Surface::Forces cur(n), next(n);
    s.evaluate(pos,u,0.,in.chi,cur);
    double u0=cur.interaction, k0=.5*norm2(vel), e0=k0+u0, e=e0, emin=e0, k=k0;
    double hk=0., work=0., max_error=0., max_work_error=0., max_p_error=0.;
    double zmin=in.z, first_negative=(e0<0. ? 0. : -1.), escape=-1., through=-1.;
    double host_px=0., host_py=0.;
    V2 support{};
    Window w15(1.5,in.z), w20(2.,in.z);
    BindingReadout binding;
    double axial=0.;
    const double h=p.dt, half=.5*h;
    for (int step=0;step<p.steps;++step) {
        double t=step*h, old_z=pos.z, old_e=e;
        vel.x+=half*cur.tag.x; vel.y+=half*cur.tag.y; vel.z+=half*cur.tag.z;
        pos.x+=h*vel.x; pos.y+=h*vel.y; pos.z+=h*vel.z;
        double max_disp=0.;
        if (p.mobile) for (int i=0;i<n;++i) {
            v[i].x+=half*cur.host[i].x; v[i].y+=half*cur.host[i].y;
            u[i].x+=h*v[i].x; u[i].y+=h*v[i].y;
            max_disp=std::max(max_disp,std::abs(u[i].x)+std::abs(u[i].y));
        }
        s.evaluate(pos,u,max_disp,in.chi,next);
        vel.x+=half*next.tag.x; vel.y+=half*next.tag.y; vel.z+=half*next.tag.z;
        hk=0.; host_px=0.; host_py=0.;
        if (p.mobile) for (int i=0;i<n;++i) {
            // Trapezoidal contact force integrated along actual site motion.
            work+=half*((cur.contact[i].x+next.contact[i].x)*v[i].x+
                        (cur.contact[i].y+next.contact[i].y)*v[i].y);
            v[i].x+=half*next.host[i].x; v[i].y+=half*next.host[i].y;
            host_px+=v[i].x; host_py+=v[i].y;
            hk+=.5*(v[i].x*v[i].x+v[i].y*v[i].y);
        }
        else {
            // Impulse taken up by the rigid support (tag exerts -F_tag).
            support.x-=half*(cur.tag.x+next.tag.x);
            support.y-=half*(cur.tag.y+next.tag.y);
        }
        std::swap(cur,next);
        k=.5*norm2(vel); e=k+cur.interaction;
        if (!std::isfinite(e+hk+cur.strain+work)) throw std::runtime_error("nonfinite MD energy");
        max_error=std::max(max_error,std::abs(e+hk+cur.strain-e0));
        max_work_error=std::max(max_work_error,std::abs(hk+cur.strain-work));
        max_p_error=std::max(max_p_error,std::hypot(vel.x+host_px+support.x-in.vel.x,
                                                   vel.y+host_py+support.y-in.vel.y));
        emin=std::min(emin,e); zmin=std::min(zmin,pos.z);
        w15.update(t,h,old_z,pos.z); w20.update(t,h,old_z,pos.z);
        double old_axial=axial;
        // A fixed mouth-depth slice in the local tangent and normal chart.
        // Its two boundary planes satisfy delta_a + R' delta_z = 0, gate.
        axial=(pos.x-start.x)*in.ax+(pos.y-start.y)*in.ay+p.gate_slope*(pos.z-start.z);
        binding.update(t,h,old_z,pos.z,old_e,e,old_axial,axial,in.gate,.5*(p.steps*h));
        if (first_negative<0. && e<0.) first_negative=t+h*old_e/(old_e-e);
        if (escape<0. && w20.entry>=0. && old_z<in.z && pos.z>=in.z)
            escape=t+h*(in.z-old_z)/(pos.z-old_z);
        if (through<0. && old_z>=0. && pos.z<0.) through=t+h*old_z/(old_z-pos.z);
    }
    double time=p.steps*h;
    bool joint_end=w15.inside(pos.z) && e<0.;
    std::cout<<std::setprecision(17)<<std::boolalpha
        <<"{\"id\":"<<in.id<<",\"side\":"<<p.side<<",\"host_sites\":"<<n
        <<",\"host_bond_directions\":"<<s.directions()<<",\"mobile\":"<<p.mobile
        <<",\"mixed\":"<<p.mixed<<",\"direct\":"<<p.direct
        <<",\"pitch\":"<<p.pitch<<",\"dt\":"<<h<<",\"steps\":"<<p.steps<<",\"t\":"<<time
        <<",\"u\":"<<in.u<<",\"v\":"<<in.v<<",\"chi\":"<<in.chi<<",\"z0\":"<<in.z
        <<",\"vx0\":"<<in.vel.x<<",\"vy0\":"<<in.vel.y<<",\"vz0\":"<<in.vel.z
        <<",\"x_end\":"<<pos.x<<",\"y_end\":"<<pos.y<<",\"z_end\":"<<pos.z
        <<",\"vx_end\":"<<vel.x<<",\"vy_end\":"<<vel.y<<",\"vz_end\":"<<vel.z
        <<",\"interaction_initial\":"<<u0<<",\"kinetic_initial\":"<<k0<<",\"energy_tag_initial\":"<<e0
        <<",\"interaction_final\":"<<cur.interaction<<",\"kinetic_final\":"<<k<<",\"energy_tag_final\":"<<e
        <<",\"min_energy_tag\":"<<emin<<",\"first_negative_energy_time\":"<<first_negative
        <<",\"bound_at_end\":"<<(e<0.)<<",\"min_z\":"<<zmin
        <<",\"host_kinetic\":"<<hk<<",\"host_strain\":"<<cur.strain<<",\"host_px\":"<<host_px
        <<",\"host_py\":"<<host_py<<",\"support_impulse_x\":"<<support.x<<",\"support_impulse_y\":"<<support.y
        <<",\"work_to_host\":"<<work<<",\"energy_error\":"<<(e+hk+cur.strain-e0)
        <<",\"max_energy_error\":"<<max_error<<",\"max_work_error\":"<<max_work_error
        <<",\"max_momentum_error\":"<<max_p_error<<",\"escape_time\":"<<escape
        <<",\"plane_crossing_time\":"<<through<<",\"dx_end\":"<<(pos.x-start.x)<<",\"dy_end\":"<<(pos.y-start.y)
        <<",\"ax\":"<<in.ax<<",\"ay\":"<<in.ay<<",\"gate\":"<<in.gate<<",\"gate_slope\":"<<p.gate_slope
        <<",\"axial_displacement\":"<<axial<<",\"gate_forward_time\":"<<binding.forward
        <<",\"gate_backward_time\":"<<binding.backward<<",\"gate_exit_time\":"<<binding.first_exit
        <<",\"gate_exit_side\":"<<binding.side<<",\"gate_before_binding\":"
        <<(binding.first_exit>=0. && (binding.joint_first<0. || binding.first_exit<binding.joint_first))
        <<",\"first_joint_time\":"<<binding.joint_first<<",\"joint_residence\":"<<binding.residence
        <<",\"gated_joint_residence\":"<<binding.gated<<",\"joint_half\":"<<binding.half
        <<",\"joint_end\":"<<joint_end<<",\"joint_both\":"<<(binding.half && joint_end);
    w15.emit("15",time,pos.z); w20.emit("20",time,pos.z);
    std::cout<<"}\n";
}
} // namespace

int main(int argc, char **argv) {
    try {
        Params p=parse(argc,argv);
        std::ifstream file;
        if (p.input!="-") {
            file.open(p.input);
            if (!file) throw std::runtime_error("could not open "+p.input);
        }
        std::istream &stream=(p.input=="-") ? std::cin : file;
        Surface s(p);
        std::string line; int count=0;
        while (std::getline(stream,line)) {
            line=line.substr(0,line.find('#'));
            std::istringstream row(line); row>>std::ws;
            if (row.eof()) continue;
            Input in{}; std::string extra;
            if (!(row>>in.id>>in.u>>in.v>>in.z>>in.vel.x>>in.vel.y>>in.vel.z>>in.chi))
                throw std::runtime_error("invalid tag row: "+line);
            row>>std::ws;
            if (!row.eof() && (!(row>>in.ax>>in.ay>>in.gate) || row>>extra))
                throw std::runtime_error("invalid optional directional gate: "+line);
            if (in.id<0 || !std::isfinite(in.u+in.v+in.z+norm2(in.vel)+in.chi+in.ax+in.ay+in.gate) ||
                in.z<0. || in.chi<0. || in.chi>1. || in.gate<0. ||
                (in.gate>0. && (std::abs(std::hypot(in.ax,in.ay)-1.)>1e-9 ||
                                in.ax*in.vel.x+in.ay*in.vel.y+p.gate_slope*in.vel.z<0.)))
                throw std::runtime_error("invalid tag row: "+line);
            trajectory(p,s,in); ++count;
        }
        if (!count || !stream.eof()) throw std::runtime_error("empty or unreadable tag input");
    } catch (const std::exception &e) {
        std::cerr<<e.what()<<"\n"; return EXIT_FAILURE;
    }
}
