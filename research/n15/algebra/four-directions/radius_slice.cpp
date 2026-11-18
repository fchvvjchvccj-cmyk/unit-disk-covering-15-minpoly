// 在一个明确的素域和半径截面上，对两个角方向的消元关系取范数。
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr uint64_t prime = 998244353;

struct Field {
    uint64_t v;
    Field(int64_t x = 0) : v((x % int64_t(prime) + prime) % prime) {}
    static Field residue(uint64_t x) { Field f; f.v = x; return f; }
    Field operator+(Field b) const {
        uint64_t sum = v + b.v;
        return residue(sum >= prime ? sum - prime : sum);
    }
    Field operator-() const { return residue(v ? prime - v : 0); }
    Field operator-(Field b) const { return *this + (-b); }
    Field operator*(Field b) const { return residue(v * b.v % prime); }
    Field power(uint64_t n) const {
        Field base = *this, out = 1;
        for (; n; n >>= 1) {
            if (n & 1) out = out * base;
            base = base * base;
        }
        return out;
    }
    Field inv() const {
        if (!v) throw std::domain_error("分母为零");
        return power(prime - 2);
    }
};

struct Pair { Field a, b; };
Pair operator+(Pair a, Pair b) { return {a.a+b.a, a.b+b.b}; }
Pair operator-(Pair a, Pair b) { return {a.a-b.a, a.b-b.b}; }
Pair operator*(Pair a, Field k) { return {a.a*k, a.b*k}; }

struct Algebra {
    Field constant, linear;
    Pair mul(Pair a, Pair b) const {
        Field product = a.b*b.b;
        return {a.a*b.a+product*constant, a.a*b.b+a.b*b.a+product*linear};
    }
    Field norm(Pair a) const {
        return a.a*(a.a+linear*a.b)-constant*a.b*a.b;
    }
};

struct DiskSample {
    Field d, shift, right, prefactor;
    DiskSample(Field value, Field r) : d(value) {
        shift = d+d.inv()-1;
        right = (Field(1)+r*d)*(Field(1)+r*d.inv()).inv();
        prefactor = d.power(36)*(d+r).power(16);
    }
};

struct Polynomial {
    nmod_poly_t data;
    Polynomial() { nmod_poly_init(data, prime); }
    Polynomial(const Polynomial&) = delete;
    ~Polynomial() { nmod_poly_clear(data); }
};

bool evaluate_branch(Field e, Field a, Field r, const std::vector<DiskSample>& points,
                     const std::vector<ulong>& abscissas, Field& out) {
    const Field ar = a.inv(), unit = Field(1)-r*r, tbase = (e-1)*r.inv()-a;
    std::array<ulong,105> values;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        Field inner = Field(1)+r*(p.shift+a), bar = Field(1)+r*(p.shift+ar);
        Field k = inner-p.right, h = p.right*bar-1, t = tbase-p.shift;
        Field h0 = -inner*k*h-r*r*k*k, h1 = (inner*bar+1-r*r*2)*k*h;
        Field h2 = -bar*k*h-r*r*h*h;
        if (!h2.v) return false;
        Field inverse = h2.inv();
        Algebra alg{-h0*inverse, -h1*inverse};
        Pair x{0, r*r*t-r*unit*e}, y{r*unit,r*r*t};
        Pair z{-r*r*t*e,-r*r*t}, w{-unit*e+r*t,unit-r*t*e+r*r*t*t};
        auto sq = [&](Pair a) { return alg.mul(a,a); };
        Pair xy = alg.mul(x,y), wz = alg.mul(w,z), xz = alg.mul(x+y,z);
        Pair c0 = sq(z), c3 = sq(x-y);
        Pair c1 = wz*t+sq(w)-xy*4-xz*2+c0;
        Pair c2 = xy*(t*t)-alg.mul(w,x+y)*t+sq(w)+c3-xy*4-xz*2;
        Pair c12 = alg.mul(c1,c2), c03 = alg.mul(c0,c3);
        Pair disc = sq(c12)-alg.mul(alg.mul(c3,c1),sq(c1))*4
                    -alg.mul(alg.mul(c0,c2),sq(c2))*4-sq(c03)*27+alg.mul(c03,c12)*18;
        values[i] = (p.prefactor*h2.power(8)*alg.norm(disc)).v;
    }
    Polynomial f, derivative;
    nmod_poly_interpolate_nmod_vec_fast(f.data, abscissas.data(), values.data(), values.size());
    if (nmod_poly_degree(f.data) != 104) {
        out = 0;
    } else {
        nmod_poly_derivative(derivative.data, f.data);
        out = Field::residue(nmod_poly_resultant(f.data, derivative.data));
    }
    return true;
}

std::vector<DiskSample> disk_samples(Field r) {
    std::vector<DiskSample> points;
    for (ulong d=1; d<=105; ++d) points.emplace_back(Field(d),r);
    return points;
}

int main(int argc, char** argv) {
    if (argc != 3) return 4;
    const Field fixed(std::stoll(argv[1]));
    const bool vary_radius = std::string(argv[2]) == "r";
    const size_t samples = vary_radius ? 57961 : 34777;
    std::vector<DiskSample> points = vary_radius ? std::vector<DiskSample>() : disk_samples(fixed);
    std::vector<ulong> abscissas, inputs, outputs;
    for (ulong d=1; d<=105; ++d) {
        abscissas.push_back(d);
    }
    // 选择分裂且分母有效的插值点；退化截面在有限预算内明确失败。
    for (ulong i=2; inputs.size()<samples+8 && i<8*samples; ++i) {
        Field a = vary_radius ? fixed : Field(i);
        Field r = vary_radius ? Field(i) : fixed;
        Field sum = r*(a-a.inv()), radical = sum*sum+4;
        ulong root = n_sqrtmod(radical.v, prime);
        if (!root) continue;
        if (vary_radius) points = disk_samples(r);
        Field e = (sum+Field(root))*Field(2).inv(), u, v;
        if (!evaluate_branch(e,a,r,points,abscissas,u) ||
            !evaluate_branch(-e.inv(),a,r,points,abscissas,v)) continue;
        inputs.push_back(i);
        Field normalization = a.power(14904);
        if (vary_radius) normalization = normalization*r.power(9936);
        outputs.push_back((normalization*u*v).v);
    }
    if (inputs.size()!=samples+8) return 1;
    Polynomial projection;
    nmod_poly_interpolate_nmod_vec_fast(projection.data,inputs.data(),outputs.data(),samples);
    for (size_t i=samples; i<inputs.size(); ++i) {
        if (nmod_poly_evaluate_nmod(projection.data,inputs[i])!=outputs[i]) return 2;
    }
    long length = nmod_poly_length(projection.data);
    if (!length) return 3;
    std::cout << "[";
    for (long i=0; i<length; ++i) {
        if (i) std::cout << ",";
        std::cout << nmod_poly_get_coeff_ui(projection.data,i);
    }
    std::cout << "]\n";
}
