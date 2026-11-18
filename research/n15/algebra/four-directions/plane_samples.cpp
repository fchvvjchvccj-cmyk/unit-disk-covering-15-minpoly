// 用完整整数系数表计算模截面；每条截面均按次数上界插值。
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using U = uint32_t;
constexpr U p = 998244353, grid = 65536, block = 128, cosets = grid/block;
constexpr size_t required = 27739;
U add(U a,U b) { U c=a+b; return c>=p ? c-p : c; }
U sub(U a,U b) { return a>=b ? a-b : a+p-b; }
U mul(U a,U b) { return uint64_t(a)*b%p; }
U power(U a,U n) {
    U b=1;
    for (;n;n>>=1) { if (n&1) b=mul(b,a); a=mul(a,a); }
    return b;
}

struct Term { U basis, r, a, d, coefficient; };
struct Grid {
    std::array<U,grid> roots;
    std::array<U,block> reverse;
    std::array<std::array<U,67>,cosets> scale;
    Grid() {
        U omega=power(3,(p-1)/grid);
        roots[0]=1;
        for (U i=1;i<grid;++i) roots[i]=mul(roots[i-1],omega);
        for (U i=0;i<block;++i) {
            U n=i, reversed=0;
            for (U j=1;j<block;j*=2) { reversed=2*reversed+(n&1); n>>=1; }
            reverse[i]=reversed;
        }
        for (U j=0;j<cosets;++j)
            for (U i=0;i<67;++i) scale[j][i]=roots[i*j%grid];
    }
    void evaluate(const std::array<U,67>& poly, U* values) const {
        for (U j=0;j<cosets;++j) {
            std::array<U,block> buffer{};
            for (U i=0;i<67;++i) buffer[reverse[i]]=mul(poly[i],scale[j][i]);
            for (U length=2;length<=block;length*=2) {
                for (U i=0;i<block;i+=length) for (U k=0;k<length/2;++k) {
                    U u=buffer[i+k], v=mul(buffer[i+k+length/2], roots[k*(grid/length)]);
                    buffer[i+k]=add(u,v); buffer[i+k+length/2]=sub(u,v);
                }
            }
            std::copy(buffer.begin(),buffer.end(),values+j*block);
        }
    }
};

struct Polynomial {
    nmod_poly_t value;
    Polynomial() { nmod_poly_init2(value,p,105); }
    ~Polynomial() { nmod_poly_clear(value); }
    Polynomial(const Polynomial&)=delete;
};

U resultant_at(const std::vector<U>& values, size_t node, U w,
               Polynomial& f, Polynomial& df) {
    for (size_t i=0;i<105;++i)
        f.value->coeffs[i]=add(values[i*grid+node],mul(w,values[(105+i)*grid+node]));
    if (!f.value->coeffs[104]) return 0;
    _nmod_poly_set_length(f.value,105);
    nmod_poly_derivative(df.value,f.value);
    return nmod_poly_resultant(f.value,df.value);
}

std::vector<U> compute(U r, const std::vector<Term>& terms, const Grid& nodes) {
    std::array<U,84> rp; rp[0]=1;
    for (size_t i=1;i<rp.size();++i) rp[i]=mul(rp[i-1],r);
    std::array<std::array<U,67>,210> coeffs{};
    for (auto t: terms) {
        auto& c=coeffs[105*t.basis+t.d][t.a]; c=add(c,mul(rp[t.r],t.coefficient));
    }
    std::vector<U> values(210*grid);
    for (U i=0;i<210;++i) nodes.evaluate(coeffs[i],values.data()+i*grid);
    std::vector<ulong> ax, z; ax.reserve(required+4); z.reserve(required+4);
    Polynomial f,df;
    for (U node=0;node<grid && ax.size()<required+4;++node) {
        U k=node/block+(node%block)*cosets;
        U a=nodes.roots[k], s=mul(r,sub(a,nodes.roots[(grid-k)%grid]));
        U root=n_sqrtmod(add(mul(s,s),4),p);
        if (!root) continue;
        U half_a=mul(a,(p+1)/2);
        U v1=resultant_at(values,node,mul(half_a,add(s,root)),f,df);
        U v2=resultant_at(values,node,mul(half_a,sub(s,root)),f,df);
        ax.push_back(a); z.push_back(mul(v1,v2));
    }
    if (ax.size()<required+4) return {};
    Polynomial plane,derivative,repeated,radical,common,layer;
    nmod_poly_interpolate_nmod_vec_fast(plane.value,ax.data(),z.data(),required);
    for (size_t i=required;i<ax.size();++i)
        if (nmod_poly_evaluate_nmod(plane.value,ax[i])!=z[i]) throw std::runtime_error("截面校验失败");
    if (nmod_poly_is_zero(plane.value)) return {};
    nmod_poly_derivative(derivative.value,plane.value);
    nmod_poly_gcd(repeated.value,plane.value,derivative.value);
    nmod_poly_div(radical.value,plane.value,repeated.value);
    nmod_poly_gcd(common.value,radical.value,repeated.value);
    nmod_poly_div(layer.value,radical.value,common.value);
    nmod_poly_make_monic(layer.value,layer.value);
    std::vector<U> out(nmod_poly_length(layer.value));
    for (size_t i=0;i<out.size();++i) out[i]=layer.value->coeffs[i];
    return out;
}

int main(int argc,char** argv) {
    if (argc!=6) return 1;
    std::ifstream data(argv[1],std::ios::binary);
    std::vector<Term> terms; Term t;
    while (data.read(reinterpret_cast<char*>(&t),sizeof(t))) {
        if (t.basis>1 || t.d>104 || t.r>83 || t.a>66 || t.coefficient>=p) return 2;
        terms.push_back(t);
    }
    if (terms.size()!=644413) return 3;
    U first=std::stoul(argv[3]); int count=std::stoi(argv[4]), threads=std::stoi(argv[5]);
    if (!first || count<1 || uint64_t(first)+count>=p || threads<1 || threads>6) return 4;
    std::ofstream out(argv[2],std::ios::binary|std::ios::app);
    const Grid nodes;
    int successes=0;
    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (int i=0;i<count;++i) {
        U r=first+i;
        auto row=compute(r,terms,nodes);
        U length=row.size();
        #pragma omp critical
        {
            out.write(reinterpret_cast<char*>(&r),sizeof(r));
            out.write(reinterpret_cast<char*>(&length),sizeof(length));
            out.write(reinterpret_cast<char*>(row.data()),length*sizeof(U));
            out.flush();
            if (length==1321) ++successes;
            if ((i+1)%64==0) std::cout << "slice " << r << " degree " << int(length)-1 << std::endl;
        }
    }
    if (!out) return 5;
    std::cout << "complete " << count << " regular " << successes << std::endl;
}
