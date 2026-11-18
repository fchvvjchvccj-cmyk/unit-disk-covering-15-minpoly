// 将单个半径上的首一方向多项式统一到 A=3 的截面，再逐列插值。
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

using U=uint32_t;
constexpr U mod=998244353;

U eval(const std::vector<U>& f,U x) {
    U out=0;
    for(auto it=f.rbegin();it!=f.rend();++it) out=(uint64_t(out)*x+*it)%mod;
    return out;
}

int main(int argc,char** argv) {
    if(argc<7) return 1;
    const U width=std::stoul(argv[4]), anchor=std::stoul(argv[5]);
    std::ifstream cross(argv[1],std::ios::binary);
    U c;
    std::vector<U> fixed;
    while(cross.read(reinterpret_cast<char*>(&c),4)) fixed.push_back(c);
    if(fixed.empty()) return 2;
    std::map<U,std::vector<U>> rows;
    for(int file=6;file<argc;++file) {
        std::ifstream data(argv[file],std::ios::binary);
        if(!data) return 3;
        U r,n;
        while(data.read(reinterpret_cast<char*>(&r),4)) {
            if(!data.read(reinterpret_cast<char*>(&n),4) || n>34777) return 4;
            std::vector<U> row(n);
            if(!data.read(reinterpret_cast<char*>(row.data()),4*n)) return 4;
            if(n!=width) continue;
            U at3=eval(row,anchor), wanted=eval(fixed,r);
            if(!at3 || !wanted) continue;
            U scale=uint64_t(wanted)*n_invmod(at3,mod)%mod;
            for(auto& v:row) v=uint64_t(v)*scale%mod;
            if(!rows.emplace(r,std::move(row)).second) return 5;
        }
    }
    size_t training=std::stoul(argv[3]);
    if(rows.size()<training+4) return 6;
    std::vector<ulong> xs,ys(rows.size());
    for(auto& entry:rows) xs.push_back(entry.first);
    std::ofstream output(argv[2],std::ios::binary);
    nmod_poly_t poly;
    nmod_poly_init(poly,mod);
    slong max_degree=-1;
    for(U i=0;i<width;++i) {
        size_t j=0;
        for(auto& entry:rows) ys[j++]=entry.second[i];
        nmod_poly_interpolate_nmod_vec_fast(poly,xs.data(),ys.data(),training);
        for(j=training;j<rows.size();++j) {
            if(nmod_poly_evaluate_nmod(poly,xs[j])!=ys[j]) {
                std::cerr << "mismatch A exponent " << i << " r " << xs[j] << std::endl;
                nmod_poly_clear(poly);
                return 7;
            }
        }
        max_degree=std::max(max_degree,nmod_poly_degree(poly));
        for(j=0;j<training;++j) {
            U coefficient=nmod_poly_get_coeff_ui(poly,j);
            output.write(reinterpret_cast<char*>(&coefficient),4);
        }
    }
    nmod_poly_clear(poly);
    std::cout << rows.size() << " " << max_degree << std::endl;
    return output ? 0 : 8;
}
